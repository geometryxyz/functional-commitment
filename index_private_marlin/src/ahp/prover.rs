use ark_ff::{Field, PrimeField, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations as EvaluationsOnDomain,
    GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

use ac_compiler::R1CSfIndex;
use ark_std::{cfg_iter_mut, end_timer, rand::RngCore, start_timer};

use crate::ahp::UnnormalizedBivariateLagrangePoly;

use super::{
    constraint_systems::LabeledPolynomial,
    indexer::{Index, Matrix},
    verifier::{VerifierFirstMsg, VerifierSecondMsg},
    AHPForR1CS, Error,
};

/// State for the AHP prover.
pub struct ProverState<'a, F: PrimeField> {
    assignment: Vec<F>,
    public_inputs: Vec<F>,
    output: Vec<F>,
    // witness_assignment: Vec<F>,
    /// Az
    z_a: Option<Vec<F>>,
    /// Bz
    z_b: Option<Vec<F>>,
    /// query bound b
    zk_bound: usize,

    z_poly: Option<LabeledPolynomial<F>>,
    mz_polys: Option<(LabeledPolynomial<F>, LabeledPolynomial<F>)>,

    index: &'a Index<F>,

    /// the random values sent by the verifier in the first round
    verifier_first_msg: Option<VerifierFirstMsg<F>>,

    /// the blinding polynomial for the first round
    mask_poly: Option<LabeledPolynomial<F>>,

    /// the blinding polynomial for the inner sumcheck
    inner_mask_poly: Option<LabeledPolynomial<F>>,

    // domain X, sized for the public input
    // domain_x: GeneralEvaluationDomain<F>,
    /// domain H, sized for constraints
    pub domain_h: GeneralEvaluationDomain<F>,

    /// domain K, sized for matrix nonzero elements
    pub domain_k: GeneralEvaluationDomain<F>, // TODO: remove pub, it's just for easier testing
}

impl<'a, F: PrimeField> ProverState<'a, F> {
    /// Get the public input.
    pub fn public_input(&self) -> Vec<F> {
        self.public_inputs.clone()
    }

    pub fn get_assignment(&self) -> Vec<F> {
        self.assignment.clone()
    }

    pub fn output(&self) -> Vec<F> {
        self.output.clone()
    }
}

/// Each prover message that is not a list of oracles is a list of field elements.
#[derive(Clone)]
pub enum ProverMsg<F: Field> {
    /// Some rounds, the prover sends only oracles. (This is actually the case for all
    /// rounds in Marlin.)
    EmptyMessage,
    /// Otherwise, it's one or more field elements.
    FieldElements(Vec<F>),
}
impl<F: Field> ark_ff::ToBytes for ProverMsg<F> {
    fn write<W: Write>(&self, w: W) -> ark_std::io::Result<()> {
        match self {
            ProverMsg::EmptyMessage => Ok(()),
            ProverMsg::FieldElements(field_elems) => field_elems.write(w),
        }
    }
}

impl<F: Field> CanonicalSerialize for ProverMsg<F> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let res: Option<Vec<F>> = match self {
            ProverMsg::EmptyMessage => None,
            ProverMsg::FieldElements(v) => Some(v.clone()),
        };
        res.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        let res: Option<Vec<F>> = match self {
            ProverMsg::EmptyMessage => None,
            ProverMsg::FieldElements(v) => Some(v.clone()),
        };
        res.serialized_size()
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let res: Option<Vec<F>> = match self {
            ProverMsg::EmptyMessage => None,
            ProverMsg::FieldElements(v) => Some(v.clone()),
        };
        res.serialize_unchecked(&mut writer)
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let res: Option<Vec<F>> = match self {
            ProverMsg::EmptyMessage => None,
            ProverMsg::FieldElements(v) => Some(v.clone()),
        };
        res.serialize_uncompressed(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        let res: Option<Vec<F>> = match self {
            ProverMsg::EmptyMessage => None,
            ProverMsg::FieldElements(v) => Some(v.clone()),
        };
        res.uncompressed_size()
    }
}

impl<F: Field> CanonicalDeserialize for ProverMsg<F> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let res = Option::<Vec<F>>::deserialize(&mut reader)?;

        if let Some(res) = res {
            Ok(ProverMsg::FieldElements(res))
        } else {
            Ok(ProverMsg::EmptyMessage)
        }
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let res = Option::<Vec<F>>::deserialize_unchecked(&mut reader)?;

        if let Some(res) = res {
            Ok(ProverMsg::FieldElements(res))
        } else {
            Ok(ProverMsg::EmptyMessage)
        }
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let res = Option::<Vec<F>>::deserialize_uncompressed(&mut reader)?;

        if let Some(res) = res {
            Ok(ProverMsg::FieldElements(res))
        } else {
            Ok(ProverMsg::EmptyMessage)
        }
    }
}

/// The first set of prover oracles.
pub struct ProverFirstOracles<F: Field> {
    /// The LDE of `z`.
    pub z: LabeledPolynomial<F>,
    /// The LDE of `Az`.
    pub z_a: LabeledPolynomial<F>,
    /// The LDE of `Bz`.
    pub z_b: LabeledPolynomial<F>,
    /// The sum-check hiding polynomial.
    pub mask_poly: LabeledPolynomial<F>,

    /// The inner sum-check hiding polynomial.
    pub inner_mask_poly: LabeledPolynomial<F>,
}

impl<F: Field> ProverFirstOracles<F> {
    /// Iterate over the polynomials output by the prover in the index private first round.
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        vec![
            &self.z,
            &self.z_a,
            &self.z_b,
            &self.mask_poly,
            &self.inner_mask_poly,
        ]
        .into_iter()
    }

    pub fn get_z(&self) -> LabeledPolynomial<F> {
        self.z.clone()
    }
}

/// The second set of prover oracles.
pub struct ProverSecondOracles<F: Field> {
    /// The polynomial `t` that is produced in the first round.
    pub t: LabeledPolynomial<F>,
    /// The polynomial `g` resulting from the first sumcheck.
    pub g_1: LabeledPolynomial<F>,
    /// The polynomial `h` resulting from the first sumcheck.
    pub h_1: LabeledPolynomial<F>,
}

impl<F: Field> ProverSecondOracles<F> {
    /// Iterate over the polynomials output by the prover in the second round.
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        vec![&self.t, &self.g_1, &self.h_1].into_iter()
    }
}

/// The third set of prover oracles.
pub struct ProverThirdOracles<F: Field> {
    /// The polynomial `f` resulting from the a / b
    pub f: LabeledPolynomial<F>,
    /// The polynomial `g` resulting from the second sumcheck.
    pub g_2: LabeledPolynomial<F>,
    // The polynomial `h` resulting from the second sumcheck.
    // pub h_2: LabeledPolynomial<F>, //TODO NOT H2
}

impl<F: Field> ProverThirdOracles<F> {
    /// Iterate over the polynomials output by the prover in the third round.
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        // &self.h_2
        vec![&self.f, &self.g_2].into_iter()
    }
}

impl<F: PrimeField> AHPForR1CS<F> {
    pub fn prover_init<'a>(
        index: &'a Index<F>,
        assignment: Vec<F>,
    ) -> Result<ProverState<'a, F>, Error> {
        // Perform matrix multiplications
        let inner_prod_fn = |row: &[(F, usize)]| {
            let mut acc = F::zero();
            for &(_, i) in row {
                acc += assignment[i];
            }
            acc
        };

        let z_a = index.a.iter().map(|row| inner_prod_fn(row)).collect();
        let z_b = index.b.iter().map(|row| inner_prod_fn(row)).collect();

        let zk_bound = 1; // One query is sufficient for our desired soundness

        if !index.index_info.check_domains_sizes::<F>() {
            return Err(Error::DomainHLargerThanDomainK);
        }

        let domain_k =
            GeneralEvaluationDomain::<F>::new(index.index_info.number_of_non_zero_entries)
                .ok_or(Error::DomainTooLarge)?;
        let domain_h = GeneralEvaluationDomain::<F>::new(index.index_info.number_of_constraints)
            .ok_or(Error::DomainTooLarge)?;

        // let domain_h = GeneralEvaluationDomain::new(r1csf_index.number_of_constraints)
        //     .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        // let domain_k = GeneralEvaluationDomain::new(r1csf_index.number_of_non_zero_entries)
        //     .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        // let domain_x = GeneralEvaluationDomain::new(r1csf_index.number_of_input_rows)
        //     .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        let public_inputs = assignment[..index.index_info.number_of_input_rows].to_vec();
        let output = assignment
            [index.index_info.number_of_constraints - index.index_info.number_of_outputs..]
            .to_vec(); // num constraints is always
                       // let witness_assignment = assignment[r1csf_index.number_of_input_rows..].to_vec();

        Ok(ProverState {
            assignment,
            public_inputs,
            output,
            z_a: Some(z_a),
            z_b: Some(z_b),
            z_poly: None,
            mz_polys: None,
            zk_bound,
            index,
            verifier_first_msg: None,
            mask_poly: None,
            inner_mask_poly: None,
            domain_h,
            domain_k,
        })
    }

    pub fn prover_first_round<'a, R: RngCore>(
        mut state: ProverState<'a, F>,
        rng: &mut R,
    ) -> Result<(ProverMsg<F>, ProverFirstOracles<F>, ProverState<'a, F>), Error> {
        let round_time = start_timer!(|| "AHP::Prover::FirstRound");
        let domain_h = state.domain_h;
        let domain_k = state.domain_k;
        let zk_bound = state.zk_bound;

        let v_h = domain_h.vanishing_polynomial().into();

        let mut z_poly =
            DensePolynomial::from_coefficients_slice(&domain_h.ifft(&state.get_assignment()));
        z_poly += &(&DensePolynomial::from_coefficients_slice(&[F::rand(rng)]) * &v_h);

        let z_a_poly_time = start_timer!(|| "Computing z_A polynomial");
        let z_a = state.z_a.clone().unwrap();
        let z_a_poly = &EvaluationsOnDomain::from_vec_and_domain(z_a, domain_h).interpolate()
            + &(&DensePolynomial::from_coefficients_slice(&[F::rand(rng)]) * &v_h);
        end_timer!(z_a_poly_time);

        let z_b_poly_time = start_timer!(|| "Computing z_B polynomial");
        let z_b = state.z_b.clone().unwrap();
        let z_b_poly = &EvaluationsOnDomain::from_vec_and_domain(z_b, domain_h).interpolate()
            + &(&DensePolynomial::from_coefficients_slice(&[F::rand(rng)]) * &v_h);
        end_timer!(z_b_poly_time);

        let mask_poly_time = start_timer!(|| "Computing mask polynomial");
        let mask_poly_degree = 3 * domain_h.size() + 2 * zk_bound - 1; // TODO: I calculated -1, paper says -2 and marlin uses -3
        let mut mask_poly = DensePolynomial::rand(mask_poly_degree, rng);

        let nh = domain_h.size();
        let upper_bound = mask_poly_degree / nh;
        let mut r_0 = F::zero();
        for i in 0..upper_bound + 1 {
            r_0 += mask_poly[nh * i];
        }

        mask_poly[0] -= &r_0;
        end_timer!(mask_poly_time);

        let inner_mask_poly_time =
            start_timer!(|| "Computing mask polynomial for inner zk sumcheck");
        let inner_mask_poly_degree = domain_k.size() - 1;
        let mut inner_mask_poly = DensePolynomial::rand(inner_mask_poly_degree, rng);
        inner_mask_poly[0] = F::zero();
        end_timer!(inner_mask_poly_time);

        let msg = ProverMsg::EmptyMessage;

        // assert!(w_poly.degree() < domain_h.size() - domain_x.size() + zk_bound);
        assert!(z_a_poly.degree() < domain_h.size() + zk_bound);
        assert!(z_b_poly.degree() < domain_h.size() + zk_bound);
        // assert!(mask_poly.degree() <= 3 * domain_h.size() + 2 * zk_bound - 3);

        let z = LabeledPolynomial::new("z".to_string(), z_poly, None, Some(1));
        let z_a = LabeledPolynomial::new("z_a".to_string(), z_a_poly, None, Some(1));
        let z_b = LabeledPolynomial::new("z_b".to_string(), z_b_poly, None, Some(1));
        let mask_poly =
            LabeledPolynomial::new("mask_poly".to_string(), mask_poly.clone(), None, None);
        let inner_mask_poly =
            LabeledPolynomial::new("inner_mask_poly".to_string(), inner_mask_poly, None, None);

        let oracles = ProverFirstOracles {
            z: z.clone(),
            z_a: z_a.clone(),
            z_b: z_b.clone(),
            mask_poly: mask_poly.clone(),
            inner_mask_poly: inner_mask_poly.clone(),
        };

        state.z_poly = Some(z);
        state.mz_polys = Some((z_a, z_b));
        state.mask_poly = Some(mask_poly);
        state.inner_mask_poly = Some(inner_mask_poly);
        end_timer!(round_time);

        Ok((msg, oracles, state))
    }

    /// Output the number of oracles sent by the prover in the first round.
    pub fn prover_num_first_round_oracles() -> usize {
        5
    }

    pub fn prover_first_round_degree_bounds(
        _info: &R1CSfIndex,
    ) -> impl Iterator<Item = Option<usize>> {
        vec![None; 5].into_iter()
    }

    fn calculate_t<'a>(
        matrices: impl Iterator<Item = &'a Matrix<F>>,
        matrix_randomizers: &[F],
        domain_h: GeneralEvaluationDomain<F>,
        r_alpha_x_on_h: Vec<F>,
    ) -> DensePolynomial<F> {
        let mut t_evals_on_h = vec![F::zero(); domain_h.size()];
        for (matrix, eta) in matrices.zip(matrix_randomizers) {
            for (r, row) in matrix.iter().enumerate() {
                for (coeff, c) in row.iter() {
                    t_evals_on_h[*c] += *eta * coeff * r_alpha_x_on_h[r];
                }
            }
        }
        EvaluationsOnDomain::from_vec_and_domain(t_evals_on_h, domain_h).interpolate()
    }

    /// Output the second round message and the next state for second round.
    pub fn prover_second_round<'a, R: RngCore>(
        ver_message: &VerifierFirstMsg<F>,
        mut state: ProverState<'a, F>,
        _r: &mut R,
    ) -> (ProverMsg<F>, ProverSecondOracles<F>, ProverState<'a, F>) {
        let round_time = start_timer!(|| "AHP::Prover::SecondRound");

        let domain_h = state.domain_h;
        let zk_bound = state.zk_bound;

        let mask_poly = state
            .mask_poly
            .as_ref()
            .expect("ProverState should include mask_poly when prover_second_round is called");

        let VerifierFirstMsg {
            alpha,
            eta_a,
            eta_b,
            eta_c,
        } = *ver_message;

        let summed_z_m_poly_time = start_timer!(|| "Compute z_m poly");
        let (z_a_poly, z_b_poly) = state.mz_polys.as_ref().unwrap();
        let z_c_poly = z_a_poly.polynomial() * z_b_poly.polynomial();

        let mut summed_z_m_coeffs = z_c_poly.coeffs;
        // Note: Can't combine these two loops, because z_c_poly has 2x the degree
        // of z_a_poly and z_b_poly, so the second loop gets truncated due to
        // the `zip`s.
        cfg_iter_mut!(summed_z_m_coeffs).for_each(|c| *c *= &eta_c);
        cfg_iter_mut!(summed_z_m_coeffs)
            .zip(&z_a_poly.polynomial().coeffs)
            .zip(&z_b_poly.polynomial().coeffs)
            .for_each(|((c, a), b)| *c += &(eta_a * a + &(eta_b * b)));

        let summed_z_m = DensePolynomial::from_coefficients_vec(summed_z_m_coeffs);
        end_timer!(summed_z_m_poly_time);

        let r_alpha_x_evals_time = start_timer!(|| "Compute r_alpha_x evals");
        let r_alpha_x_evals =
            domain_h.batch_eval_unnormalized_bivariate_lagrange_poly_with_diff_inputs(alpha);
        end_timer!(r_alpha_x_evals_time);

        let r_alpha_poly_time = start_timer!(|| "Compute r_alpha_x poly");
        let r_alpha_poly = DensePolynomial::from_coefficients_vec(domain_h.ifft(&r_alpha_x_evals));
        end_timer!(r_alpha_poly_time);

        let t_poly_time = start_timer!(|| "Compute t poly");
        let t_poly = Self::calculate_t(
            vec![&state.index.a, &state.index.b, &state.index.c].into_iter(),
            &[eta_a, eta_b, eta_c],
            state.domain_h,
            r_alpha_x_evals,
        );
        end_timer!(t_poly_time);

        let z_poly = state.z_poly.expect("Z must be calculated");
        let q_1_time = start_timer!(|| "Compute q_1 poly");

        let mul_domain_size = *[
            mask_poly.len(),
            r_alpha_poly.coeffs.len() + summed_z_m.coeffs.len(),
            t_poly.coeffs.len() + z_poly.len(),
        ]
        .iter()
        .max()
        .unwrap();
        let mul_domain = GeneralEvaluationDomain::new(mul_domain_size)
            .expect("field is not smooth enough to construct domain");
        let mut r_alpha_evals = r_alpha_poly.evaluate_over_domain_by_ref(mul_domain);
        let summed_z_m_evals = summed_z_m.evaluate_over_domain_by_ref(mul_domain);
        let z_poly_evals = z_poly.evaluate_over_domain_by_ref(mul_domain);
        let t_poly_m_evals = t_poly.evaluate_over_domain_by_ref(mul_domain);

        cfg_iter_mut!(r_alpha_evals.evals)
            .zip(&summed_z_m_evals.evals)
            .zip(&z_poly_evals.evals)
            .zip(&t_poly_m_evals.evals)
            .for_each(|(((a, b), &c), d)| {
                *a *= b;
                *a -= c * d;
            });
        let rhs = r_alpha_evals.interpolate();
        let q_1 = mask_poly.polynomial() + &rhs;
        end_timer!(q_1_time);

        let sumcheck_time = start_timer!(|| "Compute sumcheck h and g polys");
        let (h_1, x_g_1) = q_1.divide_by_vanishing_poly(domain_h).unwrap();
        let g_1 = DensePolynomial::from_coefficients_slice(&x_g_1.coeffs[1..]);
        end_timer!(sumcheck_time);

        let msg = ProverMsg::EmptyMessage;

        assert!(g_1.degree() <= domain_h.size() - 2);
        assert!(h_1.degree() <= 2 * domain_h.size() + 2 * zk_bound - 1);

        let oracles = ProverSecondOracles {
            t: LabeledPolynomial::new("t".into(), t_poly, None, None),
            g_1: LabeledPolynomial::new("g_1".into(), g_1, Some(domain_h.size() - 2), Some(1)),
            h_1: LabeledPolynomial::new("h_1".into(), h_1, None, None),
        };

        state.z_poly = Some(z_poly);
        state.verifier_first_msg = Some(*ver_message);
        end_timer!(round_time);

        (msg, oracles, state)
    }

    /// Output the number of oracles sent by the prover in the second round.
    pub fn prover_num_second_round_oracles() -> usize {
        3
    }

    /// Output the degree bounds of oracles in the second round.
    pub fn prover_second_round_degree_bounds(
        info: &R1CSfIndex,
    ) -> impl Iterator<Item = Option<usize>> {
        let h_domain_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_constraints)
                .unwrap();

        vec![None, Some(h_domain_size - 2), None].into_iter()
    }

    /// Output the second round message and the next state for third round.
    pub fn prover_third_round<'a>(
        // PC: HomomorphicPolynomialCommitment<F>, FS: FiatShamirRng>
        ver_message: &VerifierSecondMsg<F>,
        mut prover_state: ProverState<'a, F>,
    ) -> Result<(ProverMsg<F>, ProverThirdOracles<F>, ProverState<'a, F>), Error> {
        let ProverState {
            index,
            verifier_first_msg,
            domain_h,
            domain_k,
            inner_mask_poly,
            ..
        } = prover_state;

        let inner_mask_poly = inner_mask_poly.expect("Inner masking should be in the state");

        let VerifierFirstMsg {
            eta_a,
            eta_b,
            eta_c,
            alpha,
        } = verifier_first_msg.expect(
            "ProverState should include verifier_first_msg when prover_third_round is called",
        );

        let beta = ver_message.beta;

        let v_h_at_alpha = domain_h.evaluate_vanishing_polynomial(alpha);
        let v_h_at_beta = domain_h.evaluate_vanishing_polynomial(beta);

        let v_h_alpha_v_h_beta = v_h_at_alpha * v_h_at_beta;
        let eta_a_times_v_h_alpha_v_h_beta = eta_a * v_h_alpha_v_h_beta;
        let eta_b_times_v_h_alpha_v_h_beta = eta_b * v_h_alpha_v_h_beta;
        let eta_c_times_v_h_alpha_v_h_beta = eta_c * v_h_alpha_v_h_beta;

        let a_row = index.a_arith.row.polynomial();
        let a_col = index.a_arith.col.polynomial();
        let a_val = index.a_arith.val.polynomial();

        let b_row = index.b_arith.row.polynomial();
        let b_col = index.b_arith.col.polynomial();
        let b_val = index.b_arith.val.polynomial();

        let c_row = index.c_arith.row.polynomial();
        let c_col = index.c_arith.col.polynomial();
        let c_val = index.c_arith.val.polynomial();

        let a_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[beta])) - a_row)
            * &(&(DensePolynomial::from_coefficients_slice(&[alpha])) - a_col);

        let b_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[beta])) - b_row)
            * &(&(DensePolynomial::from_coefficients_slice(&[alpha])) - b_col);

        let c_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[beta])) - c_row)
            * &(&(DensePolynomial::from_coefficients_slice(&[alpha])) - c_col);

        let b_poly = &(&a_part_denom * &b_part_denom) * &c_part_denom;

        let a_part_nom =
            &(DensePolynomial::from_coefficients_slice(&[eta_a_times_v_h_alpha_v_h_beta])) * a_val;
        let b_part_nom =
            &(DensePolynomial::from_coefficients_slice(&[eta_b_times_v_h_alpha_v_h_beta])) * b_val;
        let c_part_nom =
            &(DensePolynomial::from_coefficients_slice(&[eta_c_times_v_h_alpha_v_h_beta])) * c_val;

        let a_poly = {
            let summand_0 = &a_part_nom * &(&b_part_denom * &c_part_denom);
            let summand_1 = &b_part_nom * &(&a_part_denom * &c_part_denom);
            let summand_2 = &c_part_nom * &(&a_part_denom * &b_part_denom);

            summand_0 + summand_1 + summand_2
        };

        let f_evals = domain_k
            .elements()
            .map(|elem| a_poly.evaluate(&elem) * b_poly.evaluate(&elem).inverse().unwrap())
            .collect::<Vec<F>>();

        let f_poly = EvaluationsOnDomain::from_vec_and_domain(f_evals, domain_k).interpolate();

        let (h_2, reminder) = (&a_poly - &(&b_poly * &f_poly))
            .divide_by_vanishing_poly(domain_k)
            .unwrap();

        // sanity check
        assert_eq!(reminder, DensePolynomial::zero());

        let zk_f_poly = f_poly.clone() + inner_mask_poly.polynomial().clone();
        let g_2 = DensePolynomial::from_coefficients_slice(&zk_f_poly.coeffs[1..]);

        // let g_2 = DensePolynomial::from_coefficients_slice(&f_poly.coeffs[1..]);

        /* BEGIN SANITY CHECKS */
        let r_alpha_x_evals =
            domain_h.batch_eval_unnormalized_bivariate_lagrange_poly_with_diff_inputs(alpha);
        let t_poly = Self::calculate_t(
            vec![
                &prover_state.index.a,
                &prover_state.index.b,
                &prover_state.index.c,
            ]
            .into_iter(),
            &[eta_a, eta_b, eta_c],
            prover_state.domain_h,
            r_alpha_x_evals.clone(),
        );

        let beta2 = F::from(123131u64);

        let x_poly = DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);

        let sumcheck_eval =
            t_poly.evaluate(&beta) * domain_k.size_as_field_element().inverse().unwrap();

        assert_eq!(zk_f_poly[0], sumcheck_eval);
        assert_eq!(
            zk_f_poly,
            (&x_poly * &g_2) + DensePolynomial::from_coefficients_slice(&[sumcheck_eval])
        );

        assert_eq!(
            zk_f_poly.evaluate(&beta2),
            beta2 * g_2.evaluate(&beta2) + sumcheck_eval
        );

        /* END SANITY CHECKS */

        let msg = ProverMsg::EmptyMessage;

        assert!(h_2.degree() < 6 * domain_k.size() - 6);
        assert!(g_2.degree() <= domain_k.size() - 2);

        let oracles = ProverThirdOracles {
            f: LabeledPolynomial::new("f".to_string(), f_poly.clone(), None, None),
            g_2: LabeledPolynomial::new(
                "g_2".to_string(),
                g_2.clone(),
                Some(domain_k.size() - 2),
                None,
            ),
        };

        prover_state.inner_mask_poly = None;
        Ok((msg, oracles, prover_state))
    }

    /// Output the number of oracles sent by the prover in the third round.
    pub fn prover_num_third_round_oracles() -> usize {
        2
    }

    /// Output the degree bounds of oracles in the third round.
    pub fn prover_third_round_degree_bounds(
        info: &R1CSfIndex,
    ) -> impl Iterator<Item = Option<usize>> {
        let k_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_non_zero_entries)
                .unwrap();

        vec![None, Some(k_size - 2)].into_iter()
    }
}
