#[cfg(test)]
mod test {
    use crate::{error::Error, t_strictly_lower_triangular_test::TStrictlyLowerTriangular};

    use ark_bn254::{Bn254, Fr};
    use ark_ff::to_bytes;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    type F = Fr;
    type PC = KZG10<Bn254>;

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
        // i.e. the position of the non-zero elements are:
        // (2, 0), (2, 1), (3, 1), (3, 2)

        let rng = &mut thread_rng();
        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let _gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_m_evals = vec![
            omega_2, omega_2, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3,
        ];
        let col_m_evals = vec![
            omega_0, omega_1, omega_2, omega_2, omega_2, omega_2, omega_2,
            omega_2,
            // or should it be:
            //omega_0, omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let t = 2;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_m_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_m_evals));

        let row_poly = LabeledPolynomial::new(
            String::from("row_poly"),
            row_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );
        let col_poly = LabeledPolynomial::new(
            String::from("col_poly"),
            col_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &commitments[0].clone(),
            &rands[0].clone(),
            &col_poly,
            &commitments[1].clone(),
            &rands[1].clone(),
            Some(enforced_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            TStrictlyLowerTriangular::<F, PC, FS>::verify(
                &vk,
                &ck,
                t,
                &domain_k,
                &domain_h,
                &commitments[0].clone(),
                &commitments[1].clone(),
                Some(enforced_degree_bound),
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

        let rng = &mut thread_rng();

        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

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
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let row_poly = LabeledPolynomial::new(
            String::from("row_poly"),
            row_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );
        let col_poly = LabeledPolynomial::new(
            String::from("col_poly"),
            col_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let (commitments, rands) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &commitments[0].clone(),
            &rands[0].clone(),
            &col_poly,
            &commitments[1].clone(),
            &rands[1].clone(),
            Some(enforced_degree_bound),
            &mut fs_rng,
            rng,
        );

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    #[test]
    fn test_reject_wrong_degree() {
        let rng = &mut thread_rng();
        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() - 1;
        let other_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let _gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_m_evals = vec![
            omega_2, omega_2, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3,
        ];
        let col_m_evals = vec![
            omega_0, omega_1, omega_2, omega_2, omega_2, omega_2, omega_2,
            omega_2,
            // or should it be:
            //omega_0, omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let t = 2;
        let row_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_m_evals));
        let col_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_m_evals));

        let row_poly = LabeledPolynomial::new(
            String::from("row_poly"),
            row_poly,
            Some(other_degree_bound),
            Some(enforced_hiding_bound),
        );
        let col_poly = LabeledPolynomial::new(
            String::from("col_poly"),
            col_poly,
            Some(other_degree_bound),
            Some(enforced_hiding_bound),
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound, other_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &commitments[0].clone(),
            &rands[0].clone(),
            &col_poly,
            &commitments[1].clone(),
            &rands[1].clone(),
            Some(other_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let res = TStrictlyLowerTriangular::<F, PC, FS>::verify(
            &vk,
            &ck,
            t,
            &domain_k,
            &domain_h,
            &commitments[0].clone(),
            &commitments[1].clone(),
            Some(enforced_degree_bound),
            proof,
            &mut fs_rng,
        );

        assert!(res.is_err());

        assert_eq!(res.err().unwrap(), Error::BatchCheckError)
    }

    #[test]
    #[should_panic]
    #[ignore] // TODO: include subset test and remove this ignore
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

        let rng = &mut thread_rng();

        let m = 6;
        let n = 4;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

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
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let row_poly = LabeledPolynomial::new(
            String::from("row_poly"),
            row_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );
        let col_poly = LabeledPolynomial::new(
            String::from("col_poly"),
            col_poly,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let (commitments, rands) =
            PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            &row_poly,
            &commitments[0].clone(),
            &rands[0].clone(),
            &col_poly,
            &commitments[1].clone(),
            &rands[1].clone(),
            Some(enforced_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            TStrictlyLowerTriangular::<F, PC, FS>::verify(
                &vk,
                &ck,
                t,
                &domain_k,
                &domain_h,
                &commitments[0].clone(),
                &commitments[1].clone(),
                Some(enforced_degree_bound),
                proof,
                &mut fs_rng,
            )
            .is_ok(),
            true
        );
    }
}
