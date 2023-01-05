use super::constraint_systems::LabeledPolynomial;
use super::Error;
use crate::well_formation::{construct_lagrange_basis, construct_vanishing};
use ac_compiler::R1CSfIndex;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::GeneralEvaluationDomain;
use ark_poly_commit::QuerySet;
use ark_std::rand::RngCore;

use super::AHPForR1CS;

/// State of the AHP verifier
pub struct VerifierState<F: PrimeField> {
    pub(crate) domain_h: GeneralEvaluationDomain<F>,
    pub(crate) domain_k: GeneralEvaluationDomain<F>,

    pub(crate) first_round_msg: Option<VerifierFirstMsg<F>>,
    pub(crate) second_round_msg: Option<VerifierSecondMsg<F>>,

    pub(crate) gamma: Option<F>,
}

/// First message of the verifier.
#[derive(Copy, Clone)]
pub struct VerifierFirstMsg<F> {
    /// Query for the random polynomial.
    pub alpha: F,
    /// Randomizer for the lincheck for `A`.
    pub eta_a: F,
    /// Randomizer for the lincheck for `B`.
    pub eta_b: F,
    /// Randomizer for the lincheck for `C`.
    pub eta_c: F,
}

/// Second verifier message.
#[derive(Copy, Clone)]
pub struct VerifierSecondMsg<F> {
    /// Query for the second round of polynomials.
    pub beta: F,
}

/// The first set of prover oracles.
pub struct VerifierWellFormationOracles<F: PrimeField> {
    /// The LDE of `pi`.
    pub x: LabeledPolynomial<F>,
    /// The LDE of `output`.
    pub y: LabeledPolynomial<F>,

    pub vh_gt_x: LabeledPolynomial<F>,
    pub vh_lt_y: LabeledPolynomial<F>,
}

impl<F: PrimeField> VerifierWellFormationOracles<F> {
    /// Iterate over the polynomials output by the prover in the index private first round.
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        vec![&self.x, &self.vh_gt_x, &self.y, &self.vh_lt_y].into_iter()
    }
}

impl<F: PrimeField> AHPForR1CS<F> {
    /// Output the first message and next round state.
    pub fn verifier_first_round<R: RngCore>(
        index_info: &R1CSfIndex,
        rng: &mut R,
    ) -> Result<(VerifierFirstMsg<F>, VerifierState<F>), Error> {
        if !index_info.check_domains_sizes::<F>() {
            return Err(Error::DomainHLargerThanDomainK);
        }
        let domain_k = GeneralEvaluationDomain::<F>::new(index_info.number_of_non_zero_entries)
            .ok_or(Error::DomainTooLarge)?;
        let domain_h = GeneralEvaluationDomain::<F>::new(index_info.number_of_constraints)
            .ok_or(Error::DomainTooLarge)?;

        let alpha = domain_h.sample_element_outside_domain(rng);
        let eta_a = F::rand(rng);
        let eta_b = F::rand(rng);
        let eta_c = F::rand(rng);

        let msg = VerifierFirstMsg {
            alpha,
            eta_a,
            eta_b,
            eta_c,
        };

        let new_state = VerifierState {
            domain_h,
            domain_k,
            first_round_msg: Some(msg),
            second_round_msg: None,
            gamma: None,
        };

        Ok((msg, new_state))
    }

    /// Output the second message and next round state.
    pub fn verifier_second_round<R: RngCore>(
        mut state: VerifierState<F>,
        rng: &mut R,
    ) -> (VerifierSecondMsg<F>, VerifierState<F>) {
        let beta = state.domain_h.sample_element_outside_domain(rng);
        let msg = VerifierSecondMsg { beta };
        state.second_round_msg = Some(msg);

        (msg, state)
    }

    /// Output the third message and next round state.
    pub fn verifier_third_round<R: RngCore>(
        mut state: VerifierState<F>,
        rng: &mut R,
    ) -> VerifierState<F> {
        state.gamma = Some(F::rand(rng));
        state
    }

    /// Output the query state and next round state.
    pub fn verifier_query_set<'a, R: RngCore>(
        state: VerifierState<F>,
        _: &'a mut R,
    ) -> (QuerySet<F>, VerifierState<F>) {
        let beta = state.second_round_msg.unwrap().beta;

        let gamma = state.gamma.unwrap();

        let mut query_set = QuerySet::new();

        query_set.insert(("g_1".into(), ("beta".into(), beta)));
        query_set.insert(("z_b".into(), ("beta".into(), beta)));
        query_set.insert(("t".into(), ("beta".into(), beta)));
        query_set.insert(("outer_sumcheck".into(), ("beta".into(), beta)));

        query_set.insert(("g_2".into(), ("gamma".into(), gamma)));
        query_set.insert(("f_sumcheck".into(), ("gamma".into(), gamma)));

        (query_set, state)
    }

    pub fn verifier_well_formation_oracles<'a>(
        info: &R1CSfIndex,
        public_input: &Vec<F>,
        output: &Vec<F>,
        state: &VerifierState<F>,
    ) -> VerifierWellFormationOracles<F> {
        let elems: Vec<F> = state.domain_h.elements().collect();

        let pi_roots_of_unity = &elems[..info.number_of_input_rows];
        let output_roots_of_unity = &elems[info.number_of_constraints - info.number_of_outputs..];

        let pi_bases = construct_lagrange_basis(pi_roots_of_unity);
        // let pi_evals = state.public_input;

        let output_bases = construct_lagrange_basis(output_roots_of_unity);
        // let output_evals = state.output();

        let vh_gt_x = construct_vanishing(&elems[info.number_of_input_rows..]);
        let vh_lt_y = construct_vanishing(&elems[..elems.len() - info.number_of_outputs]);

        let mut x_poly = DensePolynomial::<F>::zero();
        for (l_i, x_i) in pi_bases.iter().zip(public_input.iter()) {
            x_poly += &(l_i * *x_i);
        }

        let mut y_poly = DensePolynomial::<F>::zero();
        for (l_i, y_i) in output_bases.iter().zip(output.iter()) {
            y_poly += &(l_i * *y_i);
        }

        // let z_poly = prover_state.z_poly.expect("Z must be calculated");

        // Sanity checks
        // for elem in &elems {
        //     let z_i = z_poly.evaluate(elem);
        //     let x_i = x_poly.evaluate(elem);
        //     let vhx_i = vh_gt_x.evaluate(elem);
        //     assert_eq!((z_i - x_i) * vhx_i, F::zero());

        //     let y_i = y_poly.evaluate(elem);
        //     let vhy_i = vh_lt_y.evaluate(elem);
        //     assert_eq!((z_i - y_i) * vhy_i, F::zero());
        // }

        let x = LabeledPolynomial::new("pi_lde".to_string(), x_poly, None, None);
        let y = LabeledPolynomial::new("output_lde".to_string(), y_poly, None, None);
        let vh_gt_x = LabeledPolynomial::new("vh_gt_x".to_string(), vh_gt_x, None, None);
        let vh_lt_y = LabeledPolynomial::new("vh_lt_y".to_string(), vh_lt_y, None, None);

        VerifierWellFormationOracles {
            x,
            y,
            vh_gt_x,
            vh_lt_y,
        }
    }
}
