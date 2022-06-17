#[cfg(test)]
mod test {
    use crate::{
        commitment::{KZG10},
        t_diag::TDiag,
        util::gen_t_diag_test_polys,
        error::Error,
    };

    use ark_bn254::{Bn254, Fr};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    #[test]
    fn test_diag_matrix() {
        let mut rng = thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_t_diag_test_polys(domain_k, domain_h);

        let row_poly = polys[4].clone();
        let col_poly = polys[5].clone();
        let val_poly = polys[6].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(&mut rng),
        )
        .unwrap();

        let proof = TDiag::<F, PC, D>::prove(
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
            &domain_k,
            &domain_h,
            &mut rng,
        )
        .unwrap();

        let is_valid = TDiag::<F, PC, D>::verify(
            &vk,
            t,
            &commitments[0],
            &commitments[1],
            &commitments[2],
            &domain_h,
            &domain_k,
            proof,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_diag_matrix_error() {
        let mut rng = thread_rng();
        let m = 8;
        let n = 4;
        let t = 2;

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let polys = gen_t_diag_test_polys(domain_k, domain_h);

        let row_poly = polys[0].clone(); // This will cause an error
        let col_poly = polys[0].clone();
        let val_poly = polys[0].clone();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, rands) = PC::commit(
            &ck,
            &[row_poly.clone(), col_poly.clone(), val_poly.clone()],
            Some(&mut rng),
        )
        .unwrap();

        let proof = TDiag::<F, PC, D>::prove(
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
            &domain_k,
            &domain_h,
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
