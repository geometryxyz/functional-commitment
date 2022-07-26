use std::{ops::Add, ops::Mul};

use crate::error::Error;
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;

pub struct NewVO<F: Field> {
    // mapping_vector: Vec<usize>,
    // shifting_coefficients: Vec<F>,
    combine_function: fn(&[VOTerm<F>]) -> Result<VOTerm<F>, Error>,
}

impl<F: Field> NewVO<F> {
    pub fn new(
        // mapping_vector: Vec<usize>,
        // shifting_coefficients: Vec<F>,
        combine_function: fn(&[VOTerm<F>]) -> Result<VOTerm<F>, Error>,
    ) -> Self {
        Self {
            // mapping_vector,
            // shifting_coefficients,
            combine_function,
        }
    }

    pub fn instantiate(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
    ) -> Result<DensePolynomial<F>, Error> {
        let terms: Vec<VOTerm<_>> = concrete_oracles
            .iter()
            .map(|poly| VOTerm::Polynomial(poly.clone()))
            .collect();
        let combined = (self.combine_function)(&terms)?;
        match combined {
            VOTerm::Evaluation(_) => Err(Error::VOFailedToInstantiate),
            VOTerm::Polynomial(poly) => Ok(poly),
        }
    }

    pub fn evaluate_from_concrete_evals(&self, concrete_evals: &[F]) -> Result<F, Error> {
        let terms: Vec<VOTerm<_>> = concrete_evals
            .iter()
            .map(|eval| VOTerm::Evaluation(eval.clone()))
            .collect();
        let combined = (self.combine_function)(&terms)?;
        match combined {
            VOTerm::Evaluation(eval) => Ok(eval),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
        }
    }
}

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
    use crate::error::Error;
    use ark_bn254::Fr;
    use ark_ff::UniformRand;
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
        let expected = &a_poly + &b_poly;

        let concrete_oracles = &[a_poly.clone(), b_poly.clone()];
        let add_oracle = NewVO::new(simple_addition);

        // Check that we get the right polynomial
        let sum = add_oracle.instantiate(concrete_oracles).unwrap();
        assert_eq!(expected, sum);

        // CHeck that we combine evaluations correctly
        let eval_point = F::rand(rng);
        let sum_from_vo_evals = add_oracle
            .evaluate_from_concrete_evals(&[
                a_poly.evaluate(&eval_point),
                b_poly.evaluate(&eval_point),
            ])
            .unwrap();
        assert_eq!(expected.evaluate(&eval_point), sum_from_vo_evals)
    }

    #[test]
    fn test_mul_oracle() {
        let rng = &mut thread_rng();
        let a_poly = DensePolynomial::<F>::rand(4, rng);
        let b_poly = DensePolynomial::<F>::rand(4, rng);
        let expected = &a_poly * &b_poly;

        let concrete_oracles = &[a_poly.clone(), b_poly.clone()];
        let mul_oracle = NewVO::new(simple_mul);

        // Check that we get the right polynomial
        let prod = mul_oracle.instantiate(concrete_oracles).unwrap();
        assert_eq!(expected, prod);

        // CHeck that we combine evaluations correctly
        let eval_point = F::rand(rng);
        let prod_from_vo_evals = mul_oracle
            .evaluate_from_concrete_evals(&[
                a_poly.evaluate(&eval_point),
                b_poly.evaluate(&eval_point),
            ])
            .unwrap();
        assert_eq!(expected.evaluate(&eval_point), prod_from_vo_evals)
    }
}
