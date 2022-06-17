#[cfg(test)]
mod test {
    use crate::{
        commitment::{KZG10},
        label_polynomial,
        t_functional_triple::TFT,
        error::{Error},
    };

    use ark_poly_commit::kzg10::Randomness;
    use ark_poly_commit::LabeledPolynomial;
    use ark_bn254::{Bn254, Fr};

    use ark_ec::bn::Bn;
    use ark_bn254::{FrParameters, Parameters};
    use ark_ff::{to_bytes, Fp256};
    use ark_marlin::rng::FiatShamirRng;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
        UVPolynomial,
    };
    use ark_poly_commit::{PolynomialCommitment, LabeledCommitment};
    use ark_poly_commit::kzg10::Commitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    fn gen_polys(
        domain_k: GeneralEvaluationDomain<F>,
        domain_h: GeneralEvaluationDomain<F>,
    ) -> Vec::<LabeledPolynomial<F, DensePolynomial<F>>> {
        // M indices
        /*
            00, 01, 02, 03
            10, 11, 12, 13
            20, 21, 22, 23
            30, 31, 32, 33
        */

        // A (slt), t = 2
        /*
            0, 0, 0, 0
            1, 0, 0, 0
            1, 2, 0, 0
            0, 3, 0, 0
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // rowA_evals   =  w^1    w^2       w^2       w^3      w^3      w^3      w^3     w^3
        // colA_evals   =  w^0    w^0       w^1       w^1      w^1      w^1      w^1     w^1

        // B (stl), t = 1
        /*
            0, 0, 0, 0
            7, 0, 0, 0
            0, 0, 0, 0
            2, 2, 0, 0
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // rowB_evals   =  w^1    w^3       w^3       w^3      w^3      w^3      w^3     w^3
        // colB_evals   =  w^0    w^0       w^1       w^1      w^1      w^1      w^1     w^1

        // C (diag)
        /*
            0, 0, 0, 0
            0, 0, 0, 0
            0, 0, 2, 0
            0, 0, 0, 2
        */

        // rowM and colM are vectors that encode position of each non-zero element

        // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
        // rowC_evals   =  w^2    w^3       w^0       w^0      w^0      w^0      w^0     w^0
        // colC_evals   =  w^2    w^3       w^0       w^0      w^0      w^0      w^0     w^0
        // vaLC_evals   =   2      2         0         0        0        0        0       0

        let gamma = domain_k.element(1);

        let omega_0 = domain_h.element(0);
        let omega_1 = domain_h.element(1);
        let omega_2 = domain_h.element(2);
        let omega_3 = domain_h.element(3);

        let row_a_evals = vec![
            omega_1, omega_2, omega_2, omega_3, omega_3, omega_3, omega_3, omega_3,
        ];
        let col_a_evals = vec![
            omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let row_b_evals = vec![
            omega_1, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3,
        ];
        let col_b_evals = vec![
            omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1, omega_1,
        ];

        let row_c_evals = vec![
            omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
        ];
        let col_c_evals = vec![
            omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
        ];

        let val_c_evals = vec![
            F::from(2u64),
            F::from(2u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
            F::from(0u64),
        ];

        let row_a_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_a_evals));
        let col_a_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_a_evals));
        let row_a_poly = label_polynomial!(row_a_poly);
        let col_a_poly = label_polynomial!(col_a_poly);

        let row_b_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_b_evals));
        let col_b_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_b_evals));
        let row_b_poly = label_polynomial!(row_b_poly);
        let col_b_poly = label_polynomial!(col_b_poly);

        let row_c_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_c_evals));
        let col_c_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_c_evals));
        let val_c_poly =
            DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&val_c_evals));
        let row_c_poly = label_polynomial!(row_c_poly);
        let col_c_poly = label_polynomial!(col_c_poly);
        let val_c_poly = label_polynomial!(val_c_poly);

        vec![
            row_a_poly,
            col_a_poly,
            row_b_poly,
            col_b_poly,
            row_c_poly,
            col_c_poly,
            val_c_poly,
        ]
    }

    fn gen_commitments_and_rands(
        ck: &ark_poly_commit::sonic_pc::CommitterKey<ark_ec::bn::Bn<ark_bn254::Parameters>>,
        rng: &mut rand::prelude::ThreadRng,
        polys: Vec::<LabeledPolynomial<F, DensePolynomial<F>>>,
    ) -> Vec::<(
            Vec<LabeledCommitment<Commitment<Bn<Parameters>>>>,
            Vec<Randomness<Fp256<FrParameters>, DensePolynomial<F>>>
    )> {
        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();
        vec![
            PC::commit(
                &ck,
                &[row_a_poly.clone(), col_a_poly.clone()],
                Some(rng),
            )
            .unwrap(),
            PC::commit(
                &ck,
                &[row_b_poly.clone(), col_b_poly.clone()],
                Some(rng),
            )
            .unwrap(),
            PC::commit(
                &ck,
                &[row_c_poly.clone(), col_c_poly.clone(), val_c_poly.clone()],
                Some(rng),
            )
            .unwrap()
        ]
    }

    #[test]
    fn test_tft() {
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_polys(domain_k, domain_h);

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let mut rng = thread_rng();
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let commitments_and_rands = gen_commitments_and_rands(
            &ck,
            &mut rng,
            polys
        );

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands       = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands       = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands       = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            //a
            &row_a_poly,
            &col_a_poly,
            &a_commitments[0],
            &a_commitments[1],
            &a_rands[0],
            &a_rands[1],
            //b
            &row_b_poly,
            &col_b_poly,
            &b_commitments[0],
            &b_commitments[1],
            &b_rands[0],
            &b_rands[1],
            //c
            &row_c_poly,
            &col_c_poly,
            &val_c_poly,
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            &c_rands[0],
            &c_rands[1],
            &c_rands[2],
            &mut fs_rng,
            &mut rng,
        )
        .unwrap();

        let is_valid = TFT::<F, PC, D>::verify(
            &vk,
            &ck,
            t,
            &a_commitments[0],
            &a_commitments[1],
            &b_commitments[0],
            &b_commitments[1],
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_tft_tslt_a_proof_error() {
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_polys(domain_k, domain_h);

        let row_a_poly = polys[1].clone(); // this will make the first proof step return an error
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let mut rng = thread_rng();
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let commitments_and_rands = gen_commitments_and_rands(
            &ck,
            &mut rng,
            polys
        );

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands       = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands       = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands       = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            //a
            &row_a_poly,
            &col_a_poly,
            &a_commitments[0],
            &a_commitments[1],
            &a_rands[0],
            &a_rands[1],
            //b
            &row_b_poly,
            &col_b_poly,
            &b_commitments[0],
            &b_commitments[1],
            &b_rands[0],
            &b_rands[1],
            //c
            &row_c_poly,
            &col_c_poly,
            &val_c_poly,
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            &c_rands[0],
            &c_rands[1],
            &c_rands[2],
            &mut fs_rng,
            &mut rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(
            proof.err().unwrap(),
            Error::FEvalIsZero
        );
    }

    #[test]
    fn test_tft_tslt_b_proof_error() {
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_polys(domain_k, domain_h);

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[3].clone(); // this will make the second proof step return an error
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let mut rng = thread_rng();
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let commitments_and_rands = gen_commitments_and_rands(
            &ck,
            &mut rng,
            polys
        );

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands       = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands       = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands       = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            //a
            &row_a_poly,
            &col_a_poly,
            &a_commitments[0],
            &a_commitments[1],
            &a_rands[0],
            &a_rands[1],
            //b
            &row_b_poly,
            &col_b_poly,
            &b_commitments[0],
            &b_commitments[1],
            &b_rands[0],
            &b_rands[1],
            //c
            &row_c_poly,
            &col_c_poly,
            &val_c_poly,
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            &c_rands[0],
            &c_rands[1],
            &c_rands[2],
            &mut fs_rng,
            &mut rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(
            proof.err().unwrap(),
            Error::FEvalIsZero
        );
    }

    #[test]
    fn test_tft_tslt_c_proof_error() {
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_polys(domain_k, domain_h);

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[0].clone(); // this will make the third proof step return an error
        let col_c_poly = polys[0].clone();
        let val_c_poly = polys[0].clone();

        let max_degree = 20;
        let mut rng = thread_rng();
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let commitments_and_rands = gen_commitments_and_rands(
            &ck,
            &mut rng,
            polys
        );

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands       = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands       = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands       = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            //a
            &row_a_poly,
            &col_a_poly,
            &a_commitments[0],
            &a_commitments[1],
            &a_rands[0],
            &a_rands[1],
            //b
            &row_b_poly,
            &col_b_poly,
            &b_commitments[0],
            &b_commitments[1],
            &b_rands[0],
            &b_rands[1],
            //c
            &row_c_poly,
            &col_c_poly,
            &val_c_poly,
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            &c_rands[0],
            &c_rands[1],
            &c_rands[2],
            &mut fs_rng,
            &mut rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(
            proof.err().unwrap(),
            Error::FEvalIsZero
        );
    }
}
