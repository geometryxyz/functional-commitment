use crate::{
    util::powers_of,
    virtual_oracle::{get_term_labels, VirtualOracle},
};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{LinearCombination, PolynomialLabel};
use ark_std::marker::PhantomData;

mod prover;
mod verifier;

/// A labeled DensePolynomial with coefficients over `F`
pub type LabeledPolynomial<F> = ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>;

#[allow(dead_code)]
pub struct PIOPforZeroOverK<F: PrimeField, VO: VirtualOracle<F>> {
    _field: PhantomData<F>,
    _oracle: PhantomData<VO>,
}

impl<F: PrimeField, VO: VirtualOracle<F>> PIOPforZeroOverK<F, VO> {
    pub fn get_h_prime_labels(virtual_oracle: &VO) -> impl Iterator<Item = PolynomialLabel> {
        (0..virtual_oracle.num_of_variable_terms())
            .enumerate()
            .map(|(i, _)| format!("h_prime_{}", i))
    }

    pub fn get_r_labels(virtual_oracle: &VO) -> impl Iterator<Item = PolynomialLabel> {
        (0..virtual_oracle.num_of_variable_terms())
            .enumerate()
            .map(|(i, _)| format!("r_{}", i))
    }

    pub fn generate_h_prime_linear_combinations(
        virtual_oracle: &VO,
        concrete_oracle_labels: &[PolynomialLabel],
    ) -> Vec<LinearCombination<F>> {
        let h_labels = get_term_labels(virtual_oracle, concrete_oracle_labels);
        let h_prime_labels = Self::get_h_prime_labels(virtual_oracle);
        let mut linear_combinations: Vec<LinearCombination<F>> = Vec::new();

        // Generate all the h_primes
        h_prime_labels.enumerate().for_each(|(i, label)| {
            let lc = LinearCombination::new(
                label,
                vec![
                    (F::one(), h_labels[i].clone()),
                    (F::one(), format!("m_{}", i)),
                ],
            );
            linear_combinations.push(lc)
        });

        linear_combinations
    }

    pub fn generate_q2_linear_combination(
        virtual_oracle: &VO,
        q2_separation_challenge: F,
    ) -> LinearCombination<F> {
        // Generate q2
        let r_labels = Self::get_r_labels(virtual_oracle);
        let terms: Vec<_> = r_labels
            .zip(powers_of(q2_separation_challenge))
            .map(|(r_label, challenge_power)| (challenge_power, r_label))
            .collect();
        LinearCombination::new("q_2", terms)
    }
}
