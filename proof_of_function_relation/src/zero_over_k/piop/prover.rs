use super::PIOPforZeroOverK;
use crate::error::Error;
use crate::util::*;
use crate::zero_over_k::piop::{verifier::VerifierFirstMsg, LabeledPolynomial};
use crate::zero_over_k::VirtualOracle;
use ark_ff::{PrimeField, Zero};
use ark_marlin::ahp::prover::ProverMsg;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_std::rand::Rng;
use std::iter;

// TODO: change to use the new VirtualOracle implementation
pub struct ProverState<'a, F: PrimeField, VO: VirtualOracle<F>> {
    all_concrete_oracles: &'a [LabeledPolynomial<F>],

    virtual_oracle: &'a VO,

    /// domain K over which a virtual oracle should be equal to 0
    domain_k: GeneralEvaluationDomain<F>,

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

impl<F: PrimeField> ProverSecondOracles<F> {
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F>> {
        iter::once(&self.q_2)
    }
}

impl<F: PrimeField> PIOPforZeroOverK<F> {
    pub fn prover_init<'a, VO: VirtualOracle<F>>(
        domain: GeneralEvaluationDomain<F>,
        all_concrete_oracles: &'a [LabeledPolynomial<F>],
        virtual_oracle: &'a VO,
    ) -> Result<ProverState<'a, F, VO>, Error> {
        Ok(ProverState {
            all_concrete_oracles,
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

    pub fn prover_first_round<'a, VO: VirtualOracle<F>, R: Rng>(
        mut state: ProverState<'a, F, VO>,
        rng: &mut R,
    ) -> Result<(ProverMsg<F>, ProverFirstOracles<F>, ProverState<'a, F, VO>), Error> {
        let domain = state.domain_k;
        let alphas = state.virtual_oracle.alphas();

        // generate the masking polynomials and keep a record of the random polynomials that were used
        let (random_polynomials, masking_polynomials) = compute_maskings(&domain, &alphas, rng);

        // compute the masked oracles
        let h_primes = state
            .all_concrete_oracles
            .iter()
            .zip(masking_polynomials.iter())
            .enumerate()
            .map(|(i, (oracle, masking_poly))| {
                LabeledPolynomial::new(
                    format!("h_prime_{}", i),
                    oracle.polynomial() + masking_poly.polynomial(),
                    None,
                    None,
                )
            })
            .collect::<Vec<_>>();

        // compute the blinded virtual oracle F'[X]
        let f_prime = VO::instantiate(&h_primes, &alphas);
        let (q_1, _) = DenseOrSparsePolynomial::from(&f_prime)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(
                &domain.vanishing_polynomial(),
            ))
            .unwrap();

        // // sanity check
        // assert_eq!(r, DensePolynomial::<F>::zero());

        let msg = ProverMsg::EmptyMessage;

        let q_1 = LabeledPolynomial::new(String::from("q_1"), q_1, None, None);

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

    pub fn prover_second_round<'a, R: Rng, VO: VirtualOracle<F>>(
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
fn compute_maskings<R: Rng, F: PrimeField>(
    domain: &GeneralEvaluationDomain<F>,
    alphas: &[F],
    rng: &mut R,
) -> (Vec<LabeledPolynomial<F>>, Vec<LabeledPolynomial<F>>) {
    let num_of_concrete_oracles = alphas.len();
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

        random_polynomials.push(LabeledPolynomial::new(format!("r_{}", i), r, None, None));
        masking_polynomials.push(LabeledPolynomial::new(
            format!("m_{}", i),
            &r_shifted * &vanishing_shifted,
            None,
            None,
        ));
    }

    (random_polynomials, masking_polynomials)
}
