use crate::error::Error;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::LabeledPolynomial;

pub struct ProdVO {}

impl ProdVO {
    pub fn new() -> Self {
        Self {}
    }
}

impl<F: PrimeField> VirtualOracle<F> for ProdVO {
    fn instantiate(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        _alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() != 2 {
            return Err(Error::InstantiationError);
        }

        Ok(concrete_oracles[0].polynomial() * concrete_oracles[1].polynomial())
    }

    fn num_of_oracles(&self) -> usize {
        return 2;
    }

    fn query(&self, evals: &[F], _point: F) -> Result<F, Error> {
        if evals.len() != 2 {
            return Err(Error::EvaluationError);
        }

        Ok(evals[0] * evals[1])
    }

    /// this map encodes at which concrete oracle should h_i point
    fn mapping_vector(&self) -> Vec<usize> {
        // h0 = f0, h1 = f1
        Vec::from([0, 1])
    }
}
