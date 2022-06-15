use crate::{
    error::Error,
    virtual_oracle::VirtualOracle,
};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain};
use ark_poly_commit::LabeledPolynomial;

pub struct EqVO {}

impl EqVO {
    pub fn new() -> Self {
        Self {}
    }
}

impl<F: PrimeField> VirtualOracle<F> for EqVO {
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        _alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() != 2 {
            return Err(Error::InstantiationError);
        }

        Ok(concrete_oracles[0].polynomial() - concrete_oracles[1].polynomial())
    }

    fn instantiate_in_evals_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        _alphas: &[F],
        domain: &GeneralEvaluationDomain<F>,
    ) -> Result<Vec<F>, Error> {
        if concrete_oracles.len() != 2 {
            return Err(Error::InstantiationError);
        }

        let n = domain.size();
        let k = self.compute_scaling_factor(domain);

        let domain_kn = GeneralEvaluationDomain::<F>::new(k * n).unwrap();

        let f_evals = domain_kn.coset_fft(concrete_oracles[0].polynomial());
        let g_evals = domain_kn.coset_fft(concrete_oracles[1].polynomial());

        let vo_evals = (0..domain_kn.size())
            .map(|i| f_evals[i] - g_evals[i])
            .collect::<Vec<_>>();

        Ok(vo_evals)
    }

    fn num_of_oracles(&self) -> usize {
        return 2;
    }

    fn query(&self, evals: &[F], _point: F) -> Result<F, Error> {
        if evals.len() != 2 {
            return Err(Error::EvaluationError);
        }

        Ok(evals[0] - evals[1])
    }

    /// this map encodes at which concrete oracle should h_i point
    fn mapping_vector(&self) -> Vec<usize> {
        // h0 = f0, h1 = f1
        Vec::from([0, 1])
    }

    fn degree_bound(&self, domai_size: usize) -> usize {
        domai_size + 1
    }

    fn name(&self) -> String {
        String::from("eq")
    }
}