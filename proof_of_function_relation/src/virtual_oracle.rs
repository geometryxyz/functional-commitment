use crate::to_poly;
use crate::util::shift_dense_poly;
use ark_bn254::Fr;
use ark_ff::ToBytes;
use ark_ff::{Field, One};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{test_rng, UniformRand};

/// Encode a virtual oracle. `alphas` allow to keep track of the shift applied to the input point for each concrete oracle.
/// The `evaluate` function encodes the function `G` defined in the paper.
pub trait VirtualOracle<F: Field> {
    fn alphas(&self) -> Vec<F>;
    fn evaluate(concrete_oracles: &[DensePolynomial<F>], alphas: &[F], input: &F) -> F;
    fn query(concrete_oracle_evals: &[F]) -> F;

    fn instantiate(concrete_oracles: &[DensePolynomial<F>], alphas: &[F]) -> DensePolynomial<F>;
}

//this virtual oracle will encode f(alpha_1*x) + 2*g(alpha_2*x) - 2
#[derive(CanonicalSerialize)]
pub struct TestVirtualOracle<F: Field> {
    pub oracles: Vec<DensePolynomial<F>>,
    pub alphas: Vec<F>,
}

impl<F: Field> VirtualOracle<F> for TestVirtualOracle<F> {
    fn alphas(&self) -> Vec<F> {
        self.alphas.to_vec()
    }

    fn evaluate(concrete_oracles: &[DensePolynomial<F>], alphas: &[F], input: &F) -> F {
        concrete_oracles[0].evaluate(&(alphas[0] * input))
            + F::from(2 as u64) * concrete_oracles[1].evaluate(&(alphas[1] * input))
            - F::from(2 as u64)
    }

    fn query(concrete_oracle_evals: &[F]) -> F {
        concrete_oracle_evals[0] + F::from(2 as u64) * concrete_oracle_evals[1] - F::from(2 as u64)
    }

    fn instantiate(concrete_oracles: &[DensePolynomial<F>], alphas: &[F]) -> DensePolynomial<F> {
        shift_dense_poly(&concrete_oracles[0], &alphas[0])
            + &shift_dense_poly(&concrete_oracles[1], &alphas[1]) * F::from(2 as u64)
            + to_poly!(-F::from(2 as u64))
    }
}

#[cfg(test)]
mod tests {

    use super::{TestVirtualOracle, VirtualOracle};
    use crate::to_poly;
    use crate::util::shift_dense_poly;
    use ark_bn254::Fr;
    use ark_ff::{Field, One};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn create_virtual_oracle() {
        let rng = &mut test_rng();

        let f_coeffs = vec![
            Fr::from(-3i64),
            Fr::from(1u64),
            Fr::from(-2i64),
            Fr::from(1u64),
        ];
        let f_poly = DensePolynomial::from_coefficients_slice(&f_coeffs);

        let g_coeffs = vec![
            Fr::from(-2i64),
            Fr::from(0u64),
            Fr::from(1u64),
            Fr::from(0u64),
        ];
        let g_poly = DensePolynomial::from_coefficients_slice(&g_coeffs);

        let test_point = Fr::rand(rng);

        let f_part = f_poly.evaluate(&(Fr::from(3 as u64) * test_point));
        let g_part = g_poly.evaluate(&test_point);
        let virtual_oracle_eval_by_hand = f_part + Fr::from(2 as u64) * g_part - Fr::from(2 as u64);

        let alphas = &[Fr::from(3 as u64), Fr::one()];
        let oracles = &[f_poly, g_poly];
        let virtual_oracle_eval = TestVirtualOracle::<Fr>::evaluate(oracles, alphas, &test_point);

        let virtual_oracle_query_eval = TestVirtualOracle::<Fr>::query(&[f_part, g_part]);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_eval);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_query_eval);
    }
}
