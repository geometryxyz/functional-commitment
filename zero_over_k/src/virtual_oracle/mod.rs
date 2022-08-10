use crate::error::Error;
use ark_ff::Field;
use ark_poly_commit::{Evaluations, PolynomialLabel, QuerySet};

use self::generic_shifting_vo::vo_term::VOTerm;

pub mod generic_shifting_vo;

pub trait VirtualOracle<F: Field> {
    /// maps input concrete oracles to internal terms, e.g.:
    /// mapping_vector = [0, 0, 2] means h_0 = concrete_0, h_1 = concrete_0, h_2 = concrete_2
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

/// Returns the list of concrete oracle labels ordered according to the mapping vector
pub fn get_term_labels<F: Field, VO: VirtualOracle<F>>(
    virtual_oracle: &VO,
    concrete_oracle_labels: &[PolynomialLabel],
) -> Vec<PolynomialLabel> {
    let mut h_labels = Vec::with_capacity(virtual_oracle.num_of_variable_terms());
    for concrete_oracle_index in virtual_oracle.mapping_vector() {
        h_labels.push(concrete_oracle_labels[concrete_oracle_index].clone());
    }

    h_labels
}
