use super::PIOPforZeroOverK;
use crate::error::Error;
use crate::zero_over_k::{piop::LabeledPolynomial, VirtualOracle};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use rand::Rng;

pub struct VerifierState<'a, F: PrimeField, VO: VirtualOracle<F>> {
    virtual_oracle: &'a VO,

    /// domain K over which a virtual oracle should be equal to 0
    domain_k: GeneralEvaluationDomain<F>,

    verifier_first_message: Option<VerifierFirstMsg<F>>,

    beta_1: Option<F>,

    beta_2: Option<F>,
}

/// Verifier message
#[derive(Copy, Clone)]
pub struct VerifierFirstMsg<F: PrimeField> {
    /// Linearisation challenge for the maskings
    pub c: F,
}

pub struct VerifierQuerySet<F: PrimeField> {
    pub beta_1: F,
    pub beta_2: F,
}

impl<F: PrimeField> PIOPforZeroOverK<F> {
    pub fn verifier_init<'a, VO: VirtualOracle<F>>(
        domain_k: GeneralEvaluationDomain<F>,
        virtual_oracle: &'a VO,
    ) -> Result<VerifierState<'a, F, VO>, Error> {
        Ok(VerifierState {
            virtual_oracle,
            domain_k,
            verifier_first_message: None,
            beta_1: None,
            beta_2: None,
        })
    }

    pub fn verifier_first_round<'a, R: Rng, VO: VirtualOracle<F>>(
        mut state: VerifierState<'a, F, VO>,
        rng: &mut R,
    ) -> Result<(VerifierFirstMsg<F>, VerifierState<'a, F, VO>), Error> {
        let beta_1 = state.domain_k.sample_element_outside_domain(rng);
        let beta_2 = state.domain_k.sample_element_outside_domain(rng);
        let c = state.domain_k.sample_element_outside_domain(rng);

        let msg = VerifierFirstMsg { c };

        state.verifier_first_message = Some(msg);
        state.beta_1 = Some(beta_1);
        state.beta_2 = Some(beta_2);

        Ok((msg, state))
    }

    pub fn verifier_query_set<VO: VirtualOracle<F>>(
        state: VerifierState<F, VO>,
    ) -> Result<VerifierQuerySet<F>, Error> {
        let beta_1 = state
            .beta_1
            .expect("Verifier should have computed beta 1 at this stage");
        let beta_2 = state
            .beta_1
            .expect("Verifier should have computed beta 2 at this stage");

        Ok(VerifierQuerySet { beta_1, beta_2 })
    }
}
