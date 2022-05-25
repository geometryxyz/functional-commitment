#[cfg(test)]
mod test {
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        error::Error,
        virtual_oracle::{TestVirtualOracle, VirtualOracle},
        zero_over_k::{proof::Proof, ZeroOverK},
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::PrimeField;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use rand_core::OsRng;

    const NUM_OF_CONCRETE_ORACLES: usize = 2;

    fn zero_over_k_template<
        F: PrimeField,
        PC: HomomorphicPolynomialCommitment<F>,
        VO: VirtualOracle<F, NUM_OF_CONCRETE_ORACLES>,
    >() -> Result<(), Error> {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let alphas = &[domain.element(1), F::one()];

        //a_evals (2, 4, 6, 8) -> f(w * x) = (4, 6, 8, 2)
        //b_evals (-1, -2, -3, 0) -> 2 * g(x) = (-2, -4, -6, 0)

        //test oracle to be zero at roots of unity
        let a_evals = vec![
            F::from(2u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let a_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&a_evals));

        let b_evals = vec![
            -F::from(1u64),
            -F::from(2u64),
            -F::from(3u64),
            -F::from(0u64),
        ];
        let b_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&b_evals));

        let check_virtual_oracle = VO::instantiate(&[a_poly.clone(), b_poly.clone()], &alphas);

        let eval = check_virtual_oracle.evaluate(&domain.element(1));
        assert_eq!(eval, F::zero());

        let maximum_degree: usize = 16;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, _) = PC::trim(&pp, maximum_degree, 0, None).unwrap();

        let proof = ZeroOverK::<F, PC>::prove::<VO, OsRng, NUM_OF_CONCRETE_ORACLES>(
            &[a_poly, b_poly],
            alphas,
            &domain,
            &ck,
            &mut OsRng,
        )?;

        ZeroOverK::<F, PC>::verify::<VO, NUM_OF_CONCRETE_ORACLES>(proof, &domain, &alphas)
    }

    #[test]
    fn test_zero_over_k() {
        assert_eq!(zero_over_k_template::<Fr, KZG10<Bn254>, TestVirtualOracle<Fr, NUM_OF_CONCRETE_ORACLES>>().is_ok(), true);
    }
}
