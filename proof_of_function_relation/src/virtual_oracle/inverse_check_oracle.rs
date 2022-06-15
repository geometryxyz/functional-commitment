use crate::error::Error;
use crate::to_poly;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, UVPolynomial, GeneralEvaluationDomain, EvaluationDomain};
use ark_poly_commit::LabeledPolynomial;

pub struct InverseCheckOracle {}

impl<F: PrimeField> VirtualOracle<F> for InverseCheckOracle {
    /// abstract function which will return instantiation of virtual oracle from concrete oracles and shifting factors
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        _alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() != 2 {
            return Err(Error::InstantiationError);
        }

        Ok(
            concrete_oracles[0].polynomial() * concrete_oracles[1].polynomial()
                + to_poly!(-F::one()),
        )
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
            .map(|i| f_evals[i] * g_evals[i] - F::one())
            .collect::<Vec<_>>();

        Ok(vo_evals)
    }

    fn num_of_oracles(&self) -> usize {
        2
    }

    fn query(&self, evals: &[F], _point: F) -> Result<F, Error> {
        if evals.len() != 2 {
            return Err(Error::EvaluationError);
        }
        Ok(evals[0] * evals[1] - F::one())
    }

    /// each new (f, alpha) pair should be mapped to new h (new concrete oracle)
    /// this function provides mapping from concrete oracle inidices to h indices
    fn mapping_vector(&self) -> Vec<usize> {
        Vec::from([0, 1])
    }

    fn degree_bound(&self, domain_size: usize) -> usize {
        2 * domain_size + 2
    }
}

#[cfg(test)]
mod test {
    use super::InverseCheckOracle;
    use crate::{label_polynomial, util::sample_vector, virtual_oracle::VirtualOracle};
    use ark_bn254::Fr;
    use ark_ff::{Field, One, Zero};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_std::rand::thread_rng;

    type F = Fr;

    #[test]
    fn inverse_vo_instantiation() {
        let mut rng = thread_rng();
        let n = 8;

        let f_evals: Vec<F> = sample_vector(&mut rng, n);

        let g_evals = f_evals
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&f_evals));
        let g = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&g_evals));

        let inverse_check_vo = InverseCheckOracle {};

        let instantiated_poly = inverse_check_vo
            .instantiate_in_coeffs_form(
                &[label_polynomial!(f), label_polynomial!(g)],
                &[F::one(), F::one()],
            )
            .unwrap();
        for root_of_unity in domain.elements() {
            let eval = instantiated_poly.evaluate(&root_of_unity);
            assert_eq!(eval, F::zero());
        }
    }
}
