use super::PIOPforZeroOverK;
use crate::error::Error;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::QuerySet;
use rand::Rng;
use std::collections::HashMap;

#[derive(Copy, Clone)]
pub struct VerifierState<'a, F: PrimeField, VO: VirtualOracle<F>> {
    virtual_oracle: &'a VO,

    maximum_oracle_degree_bound: Option<usize>,

    /// domain K over which a virtual oracle should be equal to 0
    domain_k: &'a GeneralEvaluationDomain<F>,

    verifier_first_message: Option<VerifierFirstMsg<F>>,

    pub beta_1: Option<F>,

    pub beta_2: Option<F>,
}

/// Verifier message
#[derive(Copy, Clone)]
pub struct VerifierFirstMsg<F: PrimeField> {
    /// Linearisation challenge for the maskings
    pub c: F,
}

impl<F: PrimeField, VO: VirtualOracle<F>> PIOPforZeroOverK<F, VO> {
    /// Return the initial verifier state
    pub fn verifier_init<'a>(
        virtual_oracle: &'a VO,
        maximum_oracle_degree_bound: Option<usize>,
        domain_k: &'a GeneralEvaluationDomain<F>,
    ) -> Result<VerifierState<'a, F, VO>, Error> {
        Ok(VerifierState {
            virtual_oracle,
            maximum_oracle_degree_bound,
            domain_k,
            verifier_first_message: None,
            beta_1: None,
            beta_2: None,
        })
    }

    pub fn verifier_first_round<'a, R: Rng>(
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

    pub fn verifier_query_set(
        state: &VerifierState<F, VO>,
        alphas: &[F],
    ) -> Result<QuerySet<F>, Error> {
        let beta_1 = state
            .beta_1
            .expect("Verifier should have computed beta 1 at this stage");
        let beta_2 = state
            .beta_2
            .expect("Verifier should have computed beta 2 at this stage");

        // let num_of_oracles = state.virtual_oracle.num_of_oracles();

        // let mut h_prime_labels = Vec::with_capacity(num_of_oracles);
        // let mut m_labels = Vec::with_capacity(num_of_oracles);

        // for i in 0..num_of_oracles {
        //     h_prime_labels.push(format!("h_prime_{}", i));
        //     m_labels.push(format!("m_{}", i));
        // }

        // let alphas_labeled = alphas
        //     .iter()
        //     .enumerate()
        //     .map(|(i, alpha)| (format!("alpha_{}", i), *alpha))
        //     .collect::<Vec<_>>();

        // let mut query_set = QuerySet::new();
        // let h_primes_set = state.virtual_oracle.get_query_set(
        //     &h_prime_labels,
        //     &alphas_labeled.to_vec(),
        //     &(String::from("beta_1"), beta_1),
        // );
        // let m_query_set = state.virtual_oracle.get_query_set(
        //     &m_labels,
        //     &alphas_labeled.to_vec(),
        //     &(String::from("beta_2"), beta_2),
        // );

        // query_set.extend(h_primes_set);
        // query_set.extend(m_query_set);

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

        /*
         * What do we get with geo seq virtual oracle
         * we have h_primes_0 and m_primes_0
         * h_prime_0 in point alpha_0_beta_1
         * h_prime_0 in point beta_1
         *
         * m_0 in point alpha_0_beta_2
         * m_0 in point beta_2
         *
         * q_1, in point: beta_1
         * q_2, in point: beta_2
         */

        Ok(query_set)
    }

    // evaluate poly: h_prime_0, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: h_prime_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: m_0, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
    // evaluate poly: m_1, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
    // evaluate poly: q_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
    // evaluate poly: q_2, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
}
