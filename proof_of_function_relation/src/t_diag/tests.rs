#[cfg(test)]
mod test {
    use crate::{error::Error, t_diag::TDiag, util::gen_t_diag_test_polys};

    use ark_bn254::{Bn254, Fr};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::SimpleHashFiatShamirRng;
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_diag_matrix() {
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

        let row_poly = polys[4].clone();
        let col_poly = polys[5].clone();
        let val_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 1, Some(&[2, enforced_degree_bound])).unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(rng),
        )
        .unwrap();

        let proof = TDiag::<F, PC, FS>::prove(
            &ck,
            t,
            &row_poly,
            &col_poly,
            &val_poly,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &rands[0],
            &rands[1],
            &rands[2],
            Some(enforced_degree_bound),
            &domain_k,
            &domain_h,
            domain_h.size(),
            rng,
        )
        .unwrap();

        let is_valid = TDiag::<F, PC, FS>::verify(
            &vk,
            t,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            domain_h.size(),
            proof,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_diag_matrix_error() {
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

        let row_poly = polys[0].clone(); // This will cause an error
        let col_poly = polys[0].clone();
        let val_poly = polys[0].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(rng),
        )
        .unwrap();

        let proof = TDiag::<F, PC, FS>::prove(
            &ck,
            t,
            &row_poly,
            &col_poly,
            &val_poly,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &rands[0],
            &rands[1],
            &rands[2],
            Some(enforced_degree_bound),
            &domain_k,
            &domain_h,
            domain_h.size(),
            rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    #[test]
    fn test_reject_wrong_degree() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let other_degree_bound = domain_k.size() + 5;
        let enforced_hiding_bound = 1;

        let polys = gen_t_diag_test_polys(
            domain_k,
            domain_h,
            Some(other_degree_bound),
            Some(enforced_hiding_bound),
        );

        let row_poly = polys[4].clone();
        let col_poly = polys[5].clone();
        let val_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            1,
            Some(&[2, enforced_degree_bound, other_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(rng),
        )
        .unwrap();

        let proof = TDiag::<F, PC, FS>::prove(
            &ck,
            t,
            &row_poly,
            &col_poly,
            &val_poly,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &rands[0],
            &rands[1],
            &rands[2],
            Some(other_degree_bound),
            &domain_k,
            &domain_h,
            domain_h.size(),
            rng,
        )
        .unwrap();

        let is_valid = TDiag::<F, PC, FS>::verify(
            &vk,
            t,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            Some(enforced_degree_bound),
            &domain_h,
            &domain_k,
            domain_h.size(),
            proof,
        );

        assert!(is_valid.is_err());

        assert_eq!(is_valid.err().unwrap(), Error::BatchCheckError)
    }
}
