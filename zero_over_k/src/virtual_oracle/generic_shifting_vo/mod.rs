use crate::error::Error;
use crate::util::shift_dense_poly;
use crate::virtual_oracle::generic_shifting_vo::vo_term::VOTerm;
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::{Evaluations, PolynomialLabel, QuerySet};

use super::VirtualOracle;

pub mod presets;
mod tests;
pub mod vo_term;

/// A virtual oracle which shifts the concrete oracles and combines them according to a user-defined function.
/// The function can be specified as a function or a closure if one needs to capture environment variables. Presets are
/// available in this crate. See equation (2) in https://eprint.iacr.org/2021/1342 for an explicit definition.
pub struct GenericShiftingVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    mapping_vector: Vec<usize>,
    shifting_coefficients: Vec<F>,
    combine_function: T,
    minimum_oracle_length: usize,
}

impl<F, T> GenericShiftingVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    /// Constructor for an input-shifting virtual oracle
    pub fn new(
        mapping_vector: &[usize],
        shifting_coefficients: &[F],
        combine_function: T,
    ) -> Result<Self, Error> {
        let number_of_user_oracles = mapping_vector.len();

        if shifting_coefficients.len() != number_of_user_oracles {
            return Err(Error::InputLengthError(String::from(
                "mapping vector and shifting coefficients do not match",
            )));
        }

        let max_index = mapping_vector
            .iter()
            .max()
            .expect("mapping vector is empty")
            .clone();

        let minimum_oracle_length = max_index;

        Ok(Self {
            mapping_vector: mapping_vector.to_vec(),
            shifting_coefficients: shifting_coefficients.to_vec(),
            combine_function,
            minimum_oracle_length,
        })
    }

    /// Returns the polynomial that results from the combination of the given concrete oracles
    pub fn compute_polynomial(
        &self,
        concrete_oracles: &[ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>],
    ) -> Result<DensePolynomial<F>, Error> {
        self.check_conrete_oracle_length(concrete_oracles.len())?;

        let x_poly = DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);

        // Put the indeterminate variable X as terms[0]
        let mut terms: Vec<VOTerm<F>> = vec![VOTerm::Polynomial(x_poly)];

        // For each item in the mapping vector, we select the corresponding concrete oracle, apply the desired
        // shift and push the resulting polynomial as a term.
        self.mapping_vector
            .iter()
            .enumerate()
            .for_each(|(term_index, &mapped_index)| {
                let shifted = shift_dense_poly(
                    &concrete_oracles[mapped_index],
                    &self.shifting_coefficients[term_index],
                );
                terms.push(VOTerm::Polynomial(shifted))
            });

        let combined = (self.combine_function)(&terms);
        match combined {
            VOTerm::Evaluation(_) => Err(Error::VOFailedToInstantiate),
            VOTerm::Polynomial(poly) => Ok(poly),
        }
    }

    /// Check that enough oracles were provided.
    fn check_conrete_oracle_length(&self, input_length: usize) -> Result<(), Error> {
        if input_length < self.minimum_oracle_length {
            return Err(Error::InputLengthError(format!(
                "Mapping vector requires {} oracles/evaluations but only {} were provided",
                self.minimum_oracle_length, input_length
            )));
        }
        Ok(())
    }
}

impl<F, T> VirtualOracle<F> for GenericShiftingVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    fn generate_query_set(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        query_point: &(String, F),
    ) -> Result<QuerySet<F>, Error> {
        self.check_conrete_oracle_length(concrete_oracle_labels.len())?;

        let mut query_set = QuerySet::new();

        // for each term, extract the concrete oracle's label, choose the desired eval point (based on the query point and alpha),
        // and create a label for this point. Insert the resulting tuple into the query set.
        self.mapping_vector
            .iter()
            .enumerate()
            .for_each(|(term_index, &mapped_index)| {
                let poly_label = concrete_oracle_labels[mapped_index].clone();
                let eval_point = self.shifting_coefficients[term_index] * query_point.1;
                let point_label = format!("{}_times_alpha{}", query_point.0, term_index);

                query_set.insert((poly_label, (point_label, eval_point)));
            });

        Ok(query_set)
    }

    fn evaluate_from_concrete_evals(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        eval_point: &F,
        evaluations: &Evaluations<F, F>,
    ) -> Result<F, Error> {
        let mut terms: Vec<VOTerm<_>> = self
            .mapping_vector
            .iter()
            .enumerate()
            .map(|(term_index, &mapped_index)| {
                let poly_label = concrete_oracle_labels[mapped_index].clone();
                let shifted_eval_point = self.shifting_coefficients[term_index] * eval_point;
                let key = (poly_label, shifted_eval_point);

                VOTerm::Evaluation(
                    evaluations
                        .get(&key)
                        .expect("Missing a concrete oracle evaluation for VO computation")
                        .clone(),
                )
            })
            .collect();

        terms.insert(0, VOTerm::Evaluation(eval_point.clone()));

        let combined = (self.combine_function)(&terms);
        match combined {
            VOTerm::Evaluation(eval) => Ok(eval),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
        }
    }

    fn mapping_vector(&self) -> Vec<usize> {
        self.mapping_vector.clone()
    }

    fn shifting_coefficients(&self) -> Vec<F> {
        self.shifting_coefficients.clone()
    }

    fn apply_evaluation_function(&self, terms: &[VOTerm<F>]) -> VOTerm<F> {
        (self.combine_function)(terms)
    }

    fn num_of_variable_terms(&self) -> usize {
        self.mapping_vector.len()
    }
}
