use crate::error::Error;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{LabeledPolynomial, PolynomialLabel, QuerySet, Evaluations};

use self::generic_shifting_vo::vo_term::VOTerm;

pub mod generic_shifting_vo;

pub trait VirtualOracle<F: PrimeField> {
    /// Returns the list of concrete oracle labels ordered according to the mapping vector
    fn get_term_labels(&self, concrete_oracle_labels: &[PolynomialLabel]) -> Vec<PolynomialLabel> {
        let mut h_labels = Vec::with_capacity(self.num_of_variable_terms());
        for concrete_oracle_index in self.mapping_vector() {
            h_labels.push(concrete_oracle_labels[concrete_oracle_index].clone());
        }

        h_labels
    }

    /// Returns the polynomial that results from the combination of the given concrete oracles
    fn compute_polynomial(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
    ) -> Result<DensePolynomial<F>, Error>;

    /// each new (f, alpha) pair should be mapped to new h (new concrete oracle)
    /// this function provides mapping from concrete oracle indices to h indices
    fn mapping_vector(&self) -> Vec<usize>;

    /// returns the shifting coefficients (denoted alpha in the paper)
    fn shifting_coefficients(&self) -> Vec<F>;

    fn apply_evaluation_function(&self, terms: &[VOTerm<F>]) -> VOTerm<F>;

    /// Gives a count of all the terms expected by the VO function excluding the X term
    fn num_of_variable_terms(&self) -> usize;

    /// Generate a query set that will allow to evaluate the VO at a requested (labeled) point
    fn generate_query_set(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        query_point: &(String, F),
    ) -> Result<QuerySet<F>, Error>;

    /// Given evalutations of each of the concrete oracles, produce the corresponding evaluation for the virtual oracle
    fn evaluate_from_concrete_evals(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        eval_point: &F,
        evaluations: &Evaluations<F, F>,
    ) -> Result<F, Error>;
}
