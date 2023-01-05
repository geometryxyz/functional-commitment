#[cfg(test)]
mod test {
    use crate::{error::Error, t_functional_triple::TFT, util::gen_t_diag_test_polys};

    use ark_bn254::{Bn254, Fr};
    use ark_poly_commit::kzg10::Randomness;
    use ark_poly_commit::LabeledPolynomial;

    use ark_bn254::{FrParameters, Parameters};
    use ark_ec::bn::Bn;
    use ark_ff::{to_bytes, Fp256};
    use ark_marlin::rng::FiatShamirRng;
    use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly_commit::kzg10::Commitment;
    use ark_poly_commit::{LabeledCommitment, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use homomorphic_poly_commit::marlin_kzg::KZG10;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    fn gen_commitments_and_rands(
        ck: &PC::CommitterKey,
        rng: &mut rand::prelude::ThreadRng,
        polys: Vec<LabeledPolynomial<F, DensePolynomial<F>>>,
    ) -> Vec<(
        Vec<LabeledCommitment<Commitment<Bn<Parameters>>>>,
        Vec<Randomness<Fp256<FrParameters>, DensePolynomial<F>>>,
    )> {
        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();
        vec![
            PC::commit(&ck, &[row_a_poly.clone(), col_a_poly.clone()], Some(rng)).unwrap(),
            PC::commit(&ck, &[row_b_poly.clone(), col_b_poly.clone()], Some(rng)).unwrap(),
            PC::commit(
                &ck,
                &[row_c_poly.clone(), col_c_poly.clone(), val_c_poly.clone()],
                Some(rng),
            )
            .unwrap(),
        ]
    }

    #[test]
    fn test_tft() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
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
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_tft_tslt_a_proof_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = LabeledPolynomial::new(
            polys[0].label().clone(),
            polys[1].polynomial().clone(),
            polys[1].degree_bound(),
            polys[1].hiding_bound(),
        ); // this will make the first proof step return an error
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    #[test]
    fn test_tft_tslt_b_proof_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = LabeledPolynomial::new(
            polys[2].label().clone(),
            polys[3].polynomial().clone(),
            polys[3].degree_bound(),
            polys[3].hiding_bound(),
        ); // this will make the second proof step return an error
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    #[test]
    fn test_tft_tslt_c_proof_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = LabeledPolynomial::new(
            polys[4].label().clone(),
            polys[0].polynomial().clone(),
            polys[0].degree_bound(),
            polys[0].hiding_bound(),
        ); // this will make the third proof step return an error
        let col_c_poly = LabeledPolynomial::new(
            polys[5].label().clone(),
            polys[0].polynomial().clone(),
            polys[0].degree_bound(),
            polys[0].hiding_bound(),
        );
        let val_c_poly = LabeledPolynomial::new(
            polys[6].label().clone(),
            polys[0].polynomial().clone(),
            polys[0].degree_bound(),
            polys[0].hiding_bound(),
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    #[test]
    fn test_tft_a_verify_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
        )
        .unwrap();

        let is_valid = TFT::<F, PC, D>::verify(
            &vk,
            &ck,
            t,
            &a_commitments[1], // Delibrately incorrect
            &a_commitments[0],
            &b_commitments[0],
            &b_commitments[1],
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_err());

        // Test for a specific error. TODO: figure out how to bubble up errors...
        assert_eq!(
            is_valid.err().unwrap(),
            Error::ZeroOverKError(String::from("BatchCheckError"))
        );
    }

    #[test]
    fn test_tft_b_verify_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
        )
        .unwrap();

        let is_valid = TFT::<F, PC, D>::verify(
            &vk,
            &ck,
            t,
            &a_commitments[0],
            &a_commitments[1],
            &b_commitments[1], // Delibrately incorrect
            &b_commitments[0],
            &c_commitments[0],
            &c_commitments[1],
            &c_commitments[2],
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_err());

        // Test for a specific error. TODO: figure out how to bubble up errors...
        assert_eq!(
            is_valid.err().unwrap(),
            Error::ZeroOverKError(String::from("BatchCheckError"))
        );
    }

    #[test]
    fn test_tft_c_verify_error() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(enforced_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_a_poly = polys[0].clone();
        let col_a_poly = polys[1].clone();
        let row_b_poly = polys[2].clone();
        let col_b_poly = polys[3].clone();
        let row_c_poly = polys[4].clone();
        let col_c_poly = polys[5].clone();
        let val_c_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let commitments_and_rands = gen_commitments_and_rands(&ck, rng, polys);

        let a_commitments = &commitments_and_rands[0].0;
        let a_rands = &commitments_and_rands[0].1;
        let b_commitments = &commitments_and_rands[1].0;
        let b_rands = &commitments_and_rands[1].1;
        let c_commitments = &commitments_and_rands[2].0;
        let c_rands = &commitments_and_rands[2].1;

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = TFT::<F, PC, D>::prove(
            &ck,
            t,
            &domain_k,
            &domain_h,
            Some(enforced_degree_bound),
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
            rng,
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
            &c_commitments[2], // Delibrately incorrect
            &c_commitments[1],
            &c_commitments[0],
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_err());

        // Test for a specific error. TODO: figure out how to bubble up errors...
        assert_eq!(
            is_valid.err().unwrap(),
            Error::ZeroOverKError(String::from("BatchCheckError"))
        );
    }
}
