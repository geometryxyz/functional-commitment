use super::PIOPforZeroOverK;
use crate::error::Error;
use crate::zero_over_k::{piop::LabeledPolynomial, VirtualOracle};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::QuerySet;
use rand::Rng;
use std::collections::HashMap;

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
        ///map for alpha_value => point_label
        let mut evaluation_points = HashMap::new();

        evaluation_points.insert(F::one(), "My favorite book.".to_string());

        let beta_1 = state
            .beta_1
            .expect("Verifier should have computed beta 1 at this stage");
        let beta_2 = state
            .beta_2
            .expect("Verifier should have computed beta 2 at this stage");

        Ok(VerifierQuerySet { beta_1, beta_2 })
    }

    pub fn verifier_query_set_new<VO: VirtualOracle<F>>(
        state: &VerifierState<F, VO>,
    ) -> Result<QuerySet<F>, Error> {
        let alphas = state.virtual_oracle.alphas();
        let beta_1 = state
            .beta_1
            .expect("Verifier should have computed beta 1 at this stage");
        let beta_2 = state
            .beta_2
            .expect("Verifier should have computed beta 2 at this stage");

        let mut query_set = QuerySet::<F>::new();
        let mut point_evaluations: HashMap<F, String> = HashMap::new();

        point_evaluations.insert(beta_1, String::from("beta_1"));
        point_evaluations.insert(beta_2, String::from("beta_2"));

        for (i, alpha) in alphas.iter().enumerate() {
            let test_point = *alpha * beta_1;
            let label = match point_evaluations.get(&test_point) {
                Some(label) => label.clone(),
                None => {
                    let label = format!("alpha_{}_beta_1", i);
                    point_evaluations.insert(test_point, label.clone());

                    label
                }
            };
            query_set.insert((format!("h_prime_{}", i), (label, *alpha * beta_1)));

            let test_point = *alpha * beta_2;
            let label = match point_evaluations.get(&test_point) {
                Some(label) => label.clone(),
                None => {
                    let label = format!("alpha_{}_beta_2", i);
                    point_evaluations.insert(test_point, label.clone());

                    label
                }
            };

            query_set.insert((format!("m_{}", i), (label, *alpha * beta_2)));
        }

        query_set.insert((String::from("q_1"), (String::from("beta_1"), beta_1)));
        query_set.insert((String::from("q_2"), (String::from("beta_2"), beta_2)));

        for (poly_label, (point_label, point)) in &query_set {
            println!("evaluate poly: {}, in point: {} with value {}", poly_label, point_label, point);
        }
        Ok(query_set)
    }

    // evaluate poly: h_prime_0, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: h_prime_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: m_0, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
    // evaluate poly: m_1, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
    // evaluate poly: q_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: q_2, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
}
