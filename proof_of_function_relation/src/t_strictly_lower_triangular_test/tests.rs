#[cfg(test)]
mod test {
    use crate::{
        commitment::KZG10,
        label_polynomial,
        t_strictly_lower_triangular_test::TStrictlyLowerTriangular,
        error::Error,
    };

    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes};
    use ark_marlin::rng::FiatShamirRng;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
        UVPolynomial,
    };
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    #[test]
    fn test_valid_matrix() {
        // M indices
        /*
            00, 01, 02, 03
            10, 11, 12, 13
            20, 21, 22, 23
            30, 31, 32, 33
        */

        // M values with t = 2
        /*
            0, 0, 0, 0
            0, 0, 0, 0
            1, 2, 0, 0
            0, 3, 5, 0
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // row_m_evals   =  w^2    w^2       w^3       w^3      w^3      w^3      w^3     w^3
        // col_m_evals   =  w^0    w^1       w^1       w^2      w^2      w^2      w^2     w^2
        //
        // or should it be:
        //
        // col_m_evals   =  w^0    w^1       w^1       w^2      w^2      w^2      w^2     w^2

        let mut rng = thread_rng();
        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let _gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_m_evals = vec![
            omega_2, omega_2, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3,
        ];
        let col_m_evals = vec![
            omega_0, omega_1, omega_2, omega_2, omega_2, omega_2, omega_2, omega_2,
            // or should it be:
            //omega_0, omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let t = 2;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_m_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_m_evals));

        let row_poly = label_polynomial!(row_poly);
        let col_poly = label_polynomial!(col_poly);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, _) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(&mut rng)).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &col_poly,
            &commitments[0].clone(),
            &commitments[1].clone(),
            &mut fs_rng,
            &mut rng,
        )
        .unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            TStrictlyLowerTriangular::<F, PC, D>::verify(
                &vk,
                &ck,
                t,
                &domain_k,
                &domain_h,
                &commitments[0].clone(),
                &commitments[1].clone(),
                proof,
                &mut fs_rng,
            )
            .is_ok(),
            true
        );
    }

    #[test]
    fn test_outside_of_lower_triangle() {
        // M indices
        /*
            00, 01, 02, 03
            10, 11, 12, 13
            20, 21, 22, 23
            30, 31, 32, 33
        */

        // M values with t = 2
        /*
            0, 0, 0, 0
            0, 0, 0, 0
            1, 2, 5, 0
            0, 3, 0, 0
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // row_m_evals   =  w^2    w^2       w^3       w^3      w^2      w^3      w^3     w^3
        // col_m_evals   =  w^0    w^1       w^1       w^2      w^2      w^2      w^2     w^2

        let mut rng = thread_rng();

        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let _gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_m_evals = vec![
            omega_2, omega_2, omega_3, omega_3, omega_2, omega_3, omega_3, omega_3,
        ];
        let col_m_evals = vec![
            omega_0, omega_1, omega_2, omega_2, omega_2, omega_2, omega_2, omega_2,
        ];

        let t = 2;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_m_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_m_evals));

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let row_poly = label_polynomial!(row_poly);
        let col_poly = label_polynomial!(col_poly);

        let (commitments, _) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(&mut rng)).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &col_poly,
            &commitments[0].clone(),
            &commitments[1].clone(),
            &mut fs_rng,
            &mut rng,
        );

        // Test for a specific error
        assert_eq!(
            proof.err().unwrap(),
            Error::FEvalIsZero
        );
    }

    #[test]
    #[should_panic]
    fn test_not_t() {
        // M indices
        /*
            00, 01, 02, 03
            10, 11, 12, 13
            20, 21, 22, 23
            30, 31, 32, 33
        */

        // M values with t = 2
        /*
            0, 0, 0, 0
            1, 0, 0, 0
            1, 2, 0, 0
            0, 3, 0, 0
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // row_m_evals   =  w^1    w^2       w^2       w^3      w^3      w^3      w^3     w^3
        // col_m_evals   =  w^0    w^0       w^1       w^1      w^1      w^1      w^1     w^1

        let mut rng = thread_rng();

        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let _gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_m_evals = vec![
            omega_1, omega_2, omega_2, omega_3, omega_2, omega_3, omega_3, omega_3,
        ];
        let col_m_evals = vec![
            omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let t = 1;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_m_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_m_evals));

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let row_poly = label_polynomial!(row_poly);
        let col_poly = label_polynomial!(col_poly);

        let (commitments, _) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(&mut rng)).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &col_poly,
            &commitments[0].clone(),
            &commitments[1].clone(),
            &mut fs_rng,
            &mut rng,
        )
        .unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            TStrictlyLowerTriangular::<F, PC, D>::verify(
                &vk,
                &ck,
                t,
                &domain_k,
                &domain_h,
                &commitments[0].clone(),
                &commitments[1].clone(),
                proof,
                &mut fs_rng,
            )
            .is_ok(),
            true
        );
    }
}
