use crate::error::Error;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial};
use ark_poly_commit::{LabeledPolynomial, QuerySet};

pub mod add_vo;
pub mod eq_vo;
pub mod geometric_sequence_vo;
pub mod inverse_check_oracle;
pub mod normalized_vo;
pub mod prod_vo;
pub mod product_check_oracle;
pub mod square_check_oracle;

pub trait VirtualOracle<F: PrimeField> {
    /// abstract function which will return instantiation of virtual oracle from concrete oracles and shifting factors
    fn instantiate(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error>;

    fn num_of_oracles(&self) -> usize;

    fn query(&self, evals: &[F], point: F) -> Result<F, Error>;

    /// each new (f, alpha) pair should be mapped to new h (new concrete oracle)
    /// this function provides mapping from concrete oracle inidices to h indices
    fn mapping_vector(&self) -> Vec<usize>;
}

pub trait EvaluationsProvider<F: PrimeField> {
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error>;
}

impl<F: PrimeField> EvaluationsProvider<F> for Vec<LabeledPolynomial<F, DensePolynomial<F>>> {
    /// Instantiate and evaluate the virtual oracle. Returns an error if the length of
    /// concrete_oracles differs from the number of concrete oracles in the Description.
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error> {
        if self.len() != virtual_oracle.num_of_oracles() {
            return Err(Error::EvaluationError);
        }

        let poly = virtual_oracle.instantiate(&self, alpha_coeffs).unwrap();
        return Ok(poly.evaluate(&point));
    }
}

//1. add query funct to oracle vec<f>.evaluate => vo.query(vec<f>)
//2. return function from oracle and then return func(self) of vec<f>
impl<F: PrimeField> EvaluationsProvider<F> for Vec<F> {
    /// Return the evaluation of the virtual oracle given a list of evaluations of each concrete
    /// oracle. The length of the input vector of evaluations must equal to the number of concrete
    /// oracles in the Description.
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        _: &Vec<F>,
    ) -> Result<F, Error> {
        virtual_oracle.query(&self, point)
    }
}
