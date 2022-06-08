use crate::error::Error;
use crate::to_poly;
use crate::util::shift_dense_poly;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::LabeledPolynomial;
use std::iter;

pub struct GeoSequenceVO<F: PrimeField> {
    pi_s: Vec<usize>,
    ci_s: Vec<usize>,
    gamma: F,
    r: F,
}

impl<F: PrimeField> GeoSequenceVO<F> {
    pub fn new(ci_s: &Vec<usize>, gamma: F, r: F) -> Self {
        let pi_s = iter::once(0)
            .chain(ci_s.iter().scan(0, |st, elem| {
                *st += elem;
                Some(*st)
            }))
            .collect::<Vec<_>>();

        Self {
            pi_s: pi_s.clone(),
            ci_s: ci_s.clone(),
            gamma,
            r,
        }
    }

    pub fn get_pi_s(&self) -> Vec<usize> {
        self.pi_s.clone()
    }
}

impl<F: PrimeField> VirtualOracle<F> for GeoSequenceVO<F> {
    fn instantiate(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() != 2 || alphas.len() != 2 {
            return Err(Error::InstantiationError);
        }

        // construct (f(gamma * x) - r * f(x))
        let mut instantiation_poly = shift_dense_poly(&concrete_oracles[1], &alphas[1])
            + (&shift_dense_poly(&concrete_oracles[0], &alphas[0]) * -self.r);

        let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
            // construct x - y^(pi + ci - 1)
            let stitch_i = &x_poly + &to_poly!(-self.gamma.pow([(pi + ci - 1) as u64]));
            instantiation_poly = &instantiation_poly * &stitch_i;
        }

        Ok(instantiation_poly)
    }

    fn num_of_oracles(&self) -> usize {
        return 2;
    }

    fn query(&self, evals: &[F], point: F) -> Result<F, Error> {
        if evals.len() != 2 {
            return Err(Error::EvaluationError);
        }

        // evals = [shifted_f, f]
        let mut eval = evals[1] - self.r * evals[0];
        for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
            // construct x - y^(pi + ci - 1)
            let stitch_i = point - self.gamma.pow([(pi + ci - 1) as u64]);
            eval *= stitch_i;
        }

        Ok(eval)
    }

    /// this map encodes at which concrete oracle should h_i point
    fn mapping_vector(&self) -> Vec<usize> {
        // in geo sequence test we have just one concrete oracle f0, hence
        // h0 = f0, h1 = f0
        Vec::from([0, 0])
    }
}

#[cfg(test)]
mod test {

    use super::{GeoSequenceVO, VirtualOracle};
    use crate::{label_polynomial, util::generate_sequence};
    use ark_bn254::Fr;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };

    /// @param a_s The initial value of each sequence

    #[test]
    fn test_geo_seq_instantiation() {
        let r = Fr::from(2u64);
        let mut a_s = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut c_s = vec![5, 3, 10, 30];

        let m = c_s.iter().sum();

        // domain is a subset of the scalar field F.
        // it's significant because domain.element(1) = gamma, domain.element(2) = gamma^2, etc
        // m is also called the order of some subgroup (TBD)
        // m is the length that we want (it must be a power of 2 if you want to ffts and iffts; otherwise
        // you don't need it to be so).
        let domain = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let to_pad = domain.size() - m;
        if to_pad > 0 {
            a_s.push(Fr::from(0u64));
            c_s.push(to_pad);
        }

        let seq = generate_sequence::<Fr>(r, &a_s.as_slice(), &c_s.as_slice());
        let f = DensePolynomial::<Fr>::from_coefficients_slice(&domain.ifft(&seq));

        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);
        let instantiatied_geo_seq_vo = geo_seq_vo
            .instantiate(
                &[label_polynomial!(f), label_polynomial!(f)],
                &[Fr::from(1u64), domain.element(1)],
            )
            .unwrap();

        for root_of_unity in domain.elements() {
            let eval = instantiatied_geo_seq_vo.evaluate(&root_of_unity);
            assert_eq!(eval, Fr::from(0u64));
        }
    }
}
