#[cfg(test)]
mod test {
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        label_polynomial,
        t_ft::TFT,
    };

    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, FftField, Field, One, Zero};
    use ark_marlin::rng::FiatShamirRng;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use ark_std::UniformRand;
    use blake2::Blake2s;
    use rand::Rng;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    #[test]
    fn test_one() {
        // Copied from t_diag/tests.rs
        let mut rng = thread_rng();
        let m = 8;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let rowM_evals = vec![
            omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
        ];
        let colM_evals = vec![
            omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
        ];
        let valM_evals = vec![
            F::from(1u64),
            F::from(1u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
        ];

        let t = 2;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&rowM_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&colM_evals));
        let val_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&valM_evals));

        let row_poly = label_polynomial!(row_poly);
        let col_poly = label_polynomial!(col_poly);
        let val_poly = label_polynomial!(val_poly);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(&mut rng),
        )
        .unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &col_poly,
            &val_poly,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &rands[0],
            &rands[1],
            &rands[2],
            &mut fs_rng,
            &mut rng,
        ).unwrap();

        let is_valid = TFT::<F, PC, D>::verify(
            &vk,
            &ck,
            t,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_ok());

    }
}
