use crate::error::Error;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{LabeledPolynomial, PolynomialLabel};

pub mod generic_shifting_vo;

pub trait VirtualOracle<F: PrimeField> {
    /// Returns the list of concrete oracle labels ordered according to the mapping vector
    fn get_term_labels(&self, concrete_oracle_labels: &[PolynomialLabel]) -> Vec<PolynomialLabel> {
        let mut h_labels = Vec::with_capacity(self.num_of_oracles());
        for concrete_oracle_index in self.mapping_vector() {
            h_labels.push(concrete_oracle_labels[concrete_oracle_index].clone());
        }

        h_labels
    }
    /// returns instantiation of virtual oracle from concrete oracles and shifting factors
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error>;

    /// verifier always have the access to the oracle X^n through point
    fn query(&self, evals: &[F], point: F) -> Result<F, Error>;

    /// each new (f, alpha) pair should be mapped to new h (new concrete oracle)
    /// this function provides mapping from concrete oracle indices to h indices
    fn mapping_vector(&self) -> Vec<usize>;

    fn num_of_oracles(&self) -> usize;
}
