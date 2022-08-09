use std::ops::{Add, Div, Mul, Sub};

use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};

/// A term to be manipulated by the virtual oracle's function, this is either a (shifted) polynomial or an evaluation thereof.
#[derive(Clone, Debug)]
pub enum VOTerm<F: Field> {
    Evaluation(F),
    Polynomial(DensePolynomial<F>),
}

impl<F: Field> Add for VOTerm<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval + rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(DensePolynomial::from_coefficients_slice(&[eval]) + rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(poly + DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(poly + rhs_poly),
            },
        }
    }
}

impl<F: Field> Sub for VOTerm<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval - rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(&DensePolynomial::from_coefficients_slice(&[eval]) - &rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(&poly - &DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly - &rhs_poly),
            },
        }
    }
}

impl<F: FftField> Mul for VOTerm<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval * rhs_eval),
                Self::Polynomial(rhs_poly) => Self::Polynomial(&rhs_poly * eval),
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Polynomial(&poly * rhs_eval),
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly * &rhs_poly),
            },
        }
    }
}

impl<F: FftField> Div for VOTerm<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval / rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(&DensePolynomial::from_coefficients_slice(&[eval]) / &rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(&poly / &DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly / &rhs_poly),
            },
        }
    }
}

#[macro_export]
macro_rules! vo_constant {
    ($field_element:expr) => {
        VOTerm::Evaluation($field_element)
    };
}
