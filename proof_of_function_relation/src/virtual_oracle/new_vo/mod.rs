use crate::error::Error;
use crate::util::shift_dense_poly;
use crate::virtual_oracle::new_vo::vo_term::VOTerm;
use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::{Evaluations, PolynomialLabel, QuerySet};

use super::VirtualOracle;

pub mod presets;
mod tests;
pub mod vo_term;

pub struct NewVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    pub mapping_vector: Vec<usize>,
    pub shifting_coefficients: Vec<F>,
    pub combine_function: T,
    minimum_oracle_length: usize,
}

impl<F, T> NewVO<F, T>
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

    pub fn number_of_internal_terms(&self) -> usize {
        self.mapping_vector.len()
    }

    pub fn get_term_labels(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
    ) -> Vec<PolynomialLabel> {
        self.mapping_vector
            .iter()
            .map(|&mapped_index| concrete_oracle_labels[mapped_index].clone())
            .collect()
    }

    /// Returns the polynomial that results from the combination of the given concrete oracles
    pub fn compute_polynomial(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
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

    pub fn generate_query_set(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        labeled_point: &(String, F),
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
                let eval_point = self.shifting_coefficients[term_index] * labeled_point.1;
                let point_label = format!("{}_times_alpha{}", labeled_point.0, term_index);

                query_set.insert((poly_label, (point_label, eval_point)));
            });

        Ok(query_set)
    }

    /// Given evalutations of each of the concrete oracles, produce the corresponding evaluation for the virtual oracle
    pub fn evaluate_from_concrete_evals(
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

        // match terms[0] {
        //     VOTerm::Evaluation(eval) => {assert_eq!(&eval, eval_point); println!("all good")},
        //     _ => println!("bad branch"),
        // }

        let combined = (self.combine_function)(&terms);
        match combined {
            VOTerm::Evaluation(eval) => Ok(eval),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
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

impl<F, T> VirtualOracle<F> for NewVO<F, T>
where
    F: PrimeField,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>],
        _alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        let oracle_polys: Vec<_> = concrete_oracles
            .iter()
            .map(|p| p.polynomial().clone())
            .collect();
        self.compute_polynomial(&oracle_polys)
    }

    fn get_term_labels(&self, concrete_oracle_labels: &[PolynomialLabel]) -> Vec<PolynomialLabel> {
        self.get_term_labels(concrete_oracle_labels)
    }

    fn mapping_vector(&self) -> Vec<usize> {
        self.mapping_vector.clone()
    }

    fn num_of_oracles(&self) -> usize {
        self.mapping_vector.len()
    }

    fn query(&self, evals: &[F], point: F) -> Result<F, Error> {
        let terms: Vec<_> = vec![point]
            .iter()
            .chain(evals.iter())
            .map(|e| VOTerm::Evaluation(e.clone()))
            .collect();

        match (self.combine_function)(&terms) {
            VOTerm::Evaluation(res) => Ok(res),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
        }
    }
}
