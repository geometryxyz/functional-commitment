use std::{ops::Add, ops::Mul};

use crate::error::Error;
use crate::util::shift_dense_poly;
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;

pub struct NewVO<F: Field> {
    mapping_vector: Vec<usize>,
    shifting_coefficients: Vec<F>,
    combine_function: fn(&[VOTerm<F>]) -> Result<VOTerm<F>, Error>,
    minimum_oracle_length: usize,
}

impl<F: Field> NewVO<F> {
    /// Constructor for an input-shifting virtual oracle
    pub fn new(
        mapping_vector: Vec<usize>,
        shifting_coefficients: Vec<F>,
        combine_function: fn(&[VOTerm<F>]) -> Result<VOTerm<F>, Error>,
    ) -> Result<Self, Error> {
        let number_of_terms = mapping_vector.len();

        if shifting_coefficients.len() != number_of_terms {
            return Err(Error::InputLengthError(String::from(
                "mapping vetcor and shifting coefficients do not match",
            )));
        }

        let max_index = mapping_vector
            .iter()
            .max()
            .expect("mapping vector is empty")
            .clone();

        let minimum_oracle_length = max_index + 1;

        Ok(Self {
            mapping_vector,
            shifting_coefficients,
            combine_function,
            minimum_oracle_length,
        })
    }

    /// Returns the polynomial that results from the combination of the given concrete oracles
    pub fn instantiate(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() < self.minimum_oracle_length {
            return Err(Error::InputLengthError(format!(
                "Mapping vector requires {} oracles/evaluations but only {} were provided",
                self.minimum_oracle_length, concrete_oracles.len()
            )));
        }

        let mut terms: Vec<VOTerm<F>> = Vec::new();

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

        let combined = (self.combine_function)(&terms)?;
        match combined {
            VOTerm::Evaluation(_) => Err(Error::VOFailedToInstantiate),
            VOTerm::Polynomial(poly) => Ok(poly),
        }
    }

    /// Given evalutations of each of the concrete oracles, produce the corresponding evaluation for the virtual oracle
    pub fn evaluate_from_concrete_evals(&self, ordered_evaluated_terms: &[F]) -> Result<F, Error> {
        let terms: Vec<VOTerm<_>> = ordered_evaluated_terms
            .iter()
            .map(|eval| VOTerm::Evaluation(eval.clone()))
            .collect();

        let combined = (self.combine_function)(&terms)?;
        match combined {
            VOTerm::Evaluation(eval) => Ok(eval),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
        }
    }

    // /// Check that enough oracles or evaluations were provided.
    // fn check_input_length(&self, input_length: usize) -> Result<(), Error> {
    //     if input_length < self.minimum_input_length {
    //         return Err(Error::InputLengthError(format!(
    //             "Mapping vector requires {} oracles/evaluations but only {} were provided",
    //             self.minimum_input_length, input_length
    //         )));
    //     }
    //     Ok(())
    // }
}

/// A term to be manipulated by the virtual oracle's function, this is either a (shifted) polynomial or an evaluation thereof.
#[derive(Clone)]
pub enum VOTerm<F: Field> {
    Evaluation(F),
    Polynomial(DensePolynomial<F>),
}

impl<F: Field> Add for VOTerm<F> {
    type Output = Result<VOTerm<F>, Error>;

    fn add(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Ok(Self::Evaluation(eval + rhs_eval)),
                Self::Polynomial(_) => Err(Error::VOMismatchedVariants),
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(_) => Err(Error::VOMismatchedVariants),
                Self::Polynomial(rhs_poly) => Ok(Self::Polynomial(poly + rhs_poly)),
            },
        }
    }
}

impl<F: FftField> Mul for VOTerm<F> {
    type Output = Result<VOTerm<F>, Error>;

    fn mul(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Ok(Self::Evaluation(eval * rhs_eval)),
                Self::Polynomial(_) => Err(Error::VOMismatchedVariants),
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(_) => Err(Error::VOMismatchedVariants),
                Self::Polynomial(rhs_poly) => Ok(Self::Polynomial(&poly * &rhs_poly)),
            },
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{error::Error, util::shift_dense_poly};
    use crate::util::sample_vector;
    use ark_bn254::Fr;
    use ark_ff::{One, UniformRand};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use rand::thread_rng;

    use super::{NewVO, VOTerm};

    type F = Fr;

    pub fn simple_addition(concrete_terms: &[VOTerm<F>]) -> Result<VOTerm<F>, Error> {
        concrete_terms[0].clone() + concrete_terms[1].clone()
    }

    pub fn simple_mul(concrete_terms: &[VOTerm<F>]) -> Result<VOTerm<F>, Error> {
        concrete_terms[0].clone() * concrete_terms[1].clone()
    }

    #[test]
    fn test_add_oracle() {
        let rng = &mut thread_rng();
        let a_poly = DensePolynomial::<F>::rand(4, rng);
        let b_poly = DensePolynomial::<F>::rand(4, rng);
        let c_poly = DensePolynomial::<F>::rand(4, rng);
        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let expected = &shifted_c + &shifted_a;

        let add_oracle = NewVO::new(mapping_vector, shifting_coefficients, simple_addition).unwrap();

        // Check that we get the right polynomial
        let sum = add_oracle.instantiate(concrete_oracles).unwrap();
        assert_eq!(expected, sum);

        // Check that we combine evaluations correctly
        let eval_point = F::rand(rng);
        let sum_from_vo_evals = add_oracle
            .evaluate_from_concrete_evals(&[
                shifted_c.evaluate(&eval_point),
                shifted_a.evaluate(&eval_point),
            ])
            .unwrap();
        assert_eq!(expected.evaluate(&eval_point), sum_from_vo_evals)
    }

    #[test]
    fn test_mul_oracle() {
        let rng = &mut thread_rng();
        let a_poly = DensePolynomial::<F>::rand(4, rng);
        let b_poly = DensePolynomial::<F>::rand(4, rng);
        let c_poly = DensePolynomial::<F>::rand(4, rng);
        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let expected = &shifted_c * &shifted_a;

        let mul_oracle = NewVO::new(mapping_vector, shifting_coefficients, simple_mul).unwrap();

        // Check that we get the right polynomial
        let prod = mul_oracle.instantiate(concrete_oracles).unwrap();
        assert_eq!(expected, prod);

        // Check that we combine evaluations correctly
        let eval_point = F::rand(rng);
        let prod_from_vo_evals = mul_oracle
            .evaluate_from_concrete_evals(&[
                shifted_c.evaluate(&eval_point),
                shifted_a.evaluate(&eval_point),
            ])
            .unwrap();
        assert_eq!(expected.evaluate(&eval_point), prod_from_vo_evals)
    }

    #[test]
    fn test_short_input_vec() {
        // mapping vector expects there to be a concrete oracle with index 1; effectively expected at last 2 concrete oracles
        let mapping_vector = vec![1];
        let shift_coefficients = vec![F::one()];
        let add_oracle = NewVO::new(mapping_vector, shift_coefficients, simple_addition).unwrap();

        // We only provide one concrete oracle
        let err_poly = add_oracle.instantiate(&vec![DensePolynomial::<F>::default()]);
        assert!(err_poly.is_err());
        assert_eq!(
            err_poly.unwrap_err(),
            Error::InputLengthError(String::from(
                "Mapping vector requires 2 oracles/evaluations but only 1 were provided"
            ))
        );

        // // We only provide one concrete evaluation
        // let err_eval = add_oracle.evaluate_from_concrete_evals(&vec![F::default()]);
        // assert!(err_eval.is_err());
        // assert_eq!(
        //     err_eval.unwrap_err(),
        //     Error::InputLengthError(String::from(
        //         "Mapping vector requires 2 oracles/evaluations but only 1 were provided"
        //     ))
        // )
    }
}
