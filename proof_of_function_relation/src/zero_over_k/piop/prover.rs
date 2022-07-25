use super::PIOPforZeroOverK;
use crate::error::Error;
use crate::util::*;
use crate::virtual_oracle::VirtualOracle;
use crate::zero_over_k::piop::{verifier::VerifierFirstMsg, LabeledPolynomial};
use ark_ff::{PrimeField, Zero};
use ark_marlin::ahp::prover::ProverMsg;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial, univariate::DenseOrSparsePolynomial
};
use ark_std::rand::Rng;
use std::iter;

// TODO: change to use the new VirtualOracle implementation
pub struct ProverState<'a, F: PrimeField, VO: VirtualOracle<F>> {
    all_concrete_oracles: &'a [LabeledPolynomial<F>],

    maximum_oracle_degree_bound: Option<usize>,

    alphas: &'a Vec<F>,

    pub virtual_oracle: &'a VO, // TODO: made public for debugging. Private later

    /// domain K over which a virtual oracle should be equal to 0
    domain_k: &'a GeneralEvaluationDomain<F>,

    masking_polynomials: Option<Vec<LabeledPolynomial<F>>>,

    random_polynomials: Option<Vec<LabeledPolynomial<F>>>,

    // this variable is made public to avoid recomputing the masked oracles at the PIOP to SNARK compiler stage
    pub masked_oracles: Option<Vec<LabeledPolynomial<F>>>,

    q_1: Option<LabeledPolynomial<F>>,

    q_2: Option<LabeledPolynomial<F>>,

    verifier_message: Option<VerifierFirstMsg<F>>,
}

/// The first set of prover oracles
pub struct ProverFirstOracles<F: PrimeField> {
    pub masking_polynomials: Vec<LabeledPolynomial<F>>,

    pub random_polynomials: Vec<LabeledPolynomial<F>>,

    pub q_1: LabeledPolynomial<F>,
}

#[allow(dead_code)]
impl<F: PrimeField> ProverFirstOracles<F> {
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        self.masking_polynomials
            .iter()
            .chain(self.random_polynomials.iter())
            .chain(iter::once(&self.q_1))
    }
}

/// The second set of prover oracles
pub struct ProverSecondOracles<F: PrimeField> {
    pub q_2: LabeledPolynomial<F>,
}

#[allow(dead_code)]
impl<F: PrimeField> ProverSecondOracles<F> {
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        iter::once(&self.q_2)
    }
}

impl<F: PrimeField, VO: VirtualOracle<F>> PIOPforZeroOverK<F, VO> {
    /// Return the initial prover state
    pub fn prover_init<'a>(
        domain: &'a GeneralEvaluationDomain<F>,
        all_concrete_oracles: &'a [LabeledPolynomial<F>],
        maximum_oracle_degree_bound: Option<usize>,
        virtual_oracle: &'a VO,
        alphas: &'a Vec<F>,
    ) -> Result<ProverState<'a, F, VO>, Error> {
        Ok(ProverState {
            all_concrete_oracles,
            maximum_oracle_degree_bound,
            alphas,
            virtual_oracle,
            domain_k: domain,
            masking_polynomials: None,
            random_polynomials: None,
            masked_oracles: None,
            q_1: None,
            q_2: None,
            verifier_message: None,
        })
    }

    pub fn prover_first_round<'a, R: Rng>(
        mut state: ProverState<'a, F, VO>,
        rng: &mut R,
    ) -> Result<(ProverMsg<F>, ProverFirstOracles<F>, ProverState<'a, F, VO>), Error> {
        let domain = state.domain_k;
        let alphas = state.alphas;

        // generate the masking polynomials and keep a record of the random polynomials that were used
        let (random_polynomials, masking_polynomials) = compute_maskings(
            state.virtual_oracle,
            &domain,
            &alphas,
            state.maximum_oracle_degree_bound,
            rng,
        );

        let mapping_vector = state.virtual_oracle.mapping_vector();
        let mut hs = Vec::with_capacity(state.virtual_oracle.num_of_oracles());
        for concrete_oracle_index in mapping_vector {
            hs.push(state.all_concrete_oracles[concrete_oracle_index].clone())
        }

        // compute the masked oracles
        let h_primes = hs
            .iter()
            .zip(masking_polynomials.iter())
            .enumerate()
            .map(|(i, (oracle, masking_poly))| {
                LabeledPolynomial::new(
                    format!("h_prime_{}", i),
                    oracle.polynomial() + masking_poly.polynomial(),
                    state.maximum_oracle_degree_bound,
                    Some(1),
                )
            })
            .collect::<Vec<_>>();

        ///////////////////////////////////////////////////////////
        // HOW TO COMPUTE Q IN EVALS FORM
        // let k = state.virtual_oracle.compute_scaling_factor(&domain);
        // let domain_kn = GeneralEvaluationDomain::<F>::new(k * domain.size()).unwrap();
        // let vh = compute_vanishing_poly_over_coset(&domain_kn, domain.size() as u64);

        // let f_prime_evals =
        //     state
        //         .virtual_oracle
        //         .instantiate_in_evals_form(h_primes.as_slice(), &alphas, domain)?;

        // let quotient_evals = f_prime_evals
        //     .iter()
        //     .zip(vh.evals.iter())
        //     .map(|(&nominator_eval, &denominator_eval)| {
        //         nominator_eval * denominator_eval.inverse().unwrap()
        //     })
        //     .collect::<Vec<_>>();

        // let quotient =
        //     DensePolynomial::from_coefficients_slice(&domain_kn.coset_ifft(&quotient_evals));
        ///////////////////////////////////////////////////////////
        

        ///////////////////////////////////////////////////////////
        // HOW TO COMPUTE Q IN COEFFS FORM
        let f_prime = state
            .virtual_oracle
            .instantiate_in_coeffs_form(h_primes.as_slice(), &alphas)?;

        let (quotient, _r) = DenseOrSparsePolynomial::from(&f_prime)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(
                &domain.vanishing_polynomial(),
            ))
            .unwrap();
        
        // sanity check
        // assert_eq!(_r, DensePolynomial::<F>::zero());
        ///////////////////////////////////////////////////////////

        let msg = ProverMsg::EmptyMessage;

        let q_1 = LabeledPolynomial::new(String::from("q_1"), quotient, None, None);

        let oracles = ProverFirstOracles {
            masking_polynomials: masking_polynomials.clone(),
            random_polynomials: random_polynomials.clone(),
            q_1: q_1.clone(),
        };

        state.q_1 = Some(q_1);
        state.masking_polynomials = Some(masking_polynomials);
        state.random_polynomials = Some(random_polynomials);
        state.masked_oracles = Some(h_primes);

        Ok((msg, oracles, state))
    }

    pub fn prover_second_round<'a, R: Rng>(
        ver_message: &VerifierFirstMsg<F>,
        mut state: ProverState<'a, F, VO>,
        _r: &mut R,
    ) -> (ProverMsg<F>, ProverSecondOracles<F>, ProverState<'a, F, VO>) {
        let random_polynomials = state.random_polynomials.clone().expect("ProverState should include the random polynomials that were used to create the maskings in round 1");

        // q_2 is defined as r1 + c*r2 + c^2r3 + ...
        let q_2 = random_polynomials
            .iter()
            .zip(powers_of(ver_message.c))
            .fold(DensePolynomial::<F>::zero(), |acc_poly, (r, c_power)| {
                acc_poly + (r.polynomial() * c_power)
            });

        let q_2 = LabeledPolynomial::new(String::from("q_2"), q_2, None, None);

        let msg = ProverMsg::EmptyMessage;

        let oracles = ProverSecondOracles { q_2: q_2.clone() };

        state.q_2 = Some(q_2);
        state.verifier_message = Some(*ver_message);

        (msg, oracles, state)
    }
}

/// computes array of m_i = ri(alpha_i^-1) * zk(alpha_i^-1)
fn compute_maskings<R: Rng, F: PrimeField, VO: VirtualOracle<F>>(
    virtual_oracle: &VO,
    domain: &GeneralEvaluationDomain<F>,
    alphas: &[F],
    masking_bound: Option<usize>,
    rng: &mut R,
) -> (Vec<LabeledPolynomial<F>>, Vec<LabeledPolynomial<F>>) {
    let num_of_concrete_oracles = virtual_oracle.num_of_oracles();
    //r is defined as polynomial degree < 2
    let degree = 1;

    let mut random_polynomials = Vec::with_capacity(num_of_concrete_oracles);
    let mut masking_polynomials = Vec::with_capacity(num_of_concrete_oracles);

    for i in 0..num_of_concrete_oracles {
        let r = DensePolynomial::<F>::rand(degree, rng);
        let shifting_factor = alphas[i].inverse().unwrap();
        let r_shifted = shift_dense_poly(&r, &shifting_factor);
        let vanishing_shifted =
            shift_dense_poly(&domain.vanishing_polynomial().into(), &shifting_factor);

        random_polynomials.push(LabeledPolynomial::new(
            format!("r_{}", i),
            r,
            Some(2),
            Some(1),
        ));
        // random_polynomials.push(LabeledPolynomial::new(format!("r_{}", i), r, None, None));
        masking_polynomials.push(LabeledPolynomial::new(
            format!("m_{}", i),
            &r_shifted * &vanishing_shifted,
            masking_bound,
            Some(1),
        ));
    }

    (random_polynomials, masking_polynomials)
}
