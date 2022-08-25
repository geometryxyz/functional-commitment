#[cfg(test)]
mod tests {
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, Field, SquareRootField};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use homomorphic_poly_commit::marlin_kzg::KZG10;

    use crate::{discrete_log_comparison::DLComparison, error::Error};
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_square_root_friendly() {
        let m = 6;
        let n = 4;

        let domain = GeneralEvaluationDomain::new(m).unwrap();
        let domain2 = GeneralEvaluationDomain::new(n).unwrap();

        let generator: F = domain.element(1);
        let g2: F = domain2.element(1);

        assert_eq!(g2.pow(&[n as u64]), F::from(1u64));

        // assert_eq!(generator, g2);

        let sq_root = generator.sqrt().unwrap();
        let sq_root2 = g2.sqrt().unwrap();

        assert_eq!(sq_root2.pow(&[2 * n as u64]), F::from(1u64));

        // println!("{}", sq_root);

        assert_eq!(generator, sq_root * sq_root);
        assert_eq!(g2, sq_root2 * sq_root2);
    }

    #[test]
    fn test_discrete_log_proof() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;

        let enforced_degree_bound = m + 1;
        let enforced_hiding_bound = Some(1);

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        // For the test to pass, each value in f_evals must be less than its corresponding value in
        // g_evals, mod n (since we use domain_h)

        let f_evals = vec![
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
        ];
        let g_evals = vec![
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
        ];

        let f_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&f_evals));
        let f_poly = LabeledPolynomial::new(
            String::from("f_poly"),
            f_poly,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );
        let g_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&g_evals));
        let g_poly = LabeledPolynomial::new(
            String::from("g_poly"),
            g_poly,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 1, Some(&[2, enforced_degree_bound])).unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[f_poly.clone(), g_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = DLComparison::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &domain_h,
            &f_poly,
            &commitments[0],
            &rands[0],
            &g_poly,
            &commitments[1],
            &rands[1],
            Some(enforced_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let res = DLComparison::verify(
            &vk,
            &ck,
            &domain_k,
            &domain_h,
            &commitments[0],
            &commitments[1],
            Some(enforced_degree_bound),
            proof,
            &mut fs_rng,
        )
        .unwrap();

        assert_eq!((), res)
    }

    #[test]
    fn test_malicious_discrete_log_proof() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;

        let enforced_degree_bound = m + 1;
        let enforced_hiding_bound = Some(1);

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let f_evals = vec![
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
        ];
        let g_evals = vec![
            domain_h.element(3), // here log(g) > log(f)
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
        ];

        let f_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&f_evals));
        let f_poly = LabeledPolynomial::new(
            String::from("f_poly"),
            f_poly,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );
        let g_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&g_evals));
        let g_poly = LabeledPolynomial::new(
            String::from("g_poly"),
            g_poly,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 1, Some(&[2, enforced_degree_bound])).unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[f_poly.clone(), g_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = DLComparison::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &domain_h,
            &f_poly,
            &commitments[0],
            &rands[0],
            &g_poly,
            &commitments[1],
            &rands[1],
            Some(enforced_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let res = DLComparison::verify(
            &vk,
            &ck,
            &domain_k,
            &domain_h,
            &commitments[0],
            &commitments[1],
            Some(enforced_degree_bound),
            proof,
            &mut fs_rng,
        );
        assert!(res.is_err());

        // Test for a specific error
        assert_eq!(
            res.err().unwrap(),
            Error::ZeroOverKError(String::from("Check2Failed"))
        );
    }

    #[test]
    fn test_reject_large_degree() {
        let rng = &mut thread_rng();
        let m = 8;
        let n = 4;

        let enforced_degree_bound = m + 1;
        let other_degree_bound = m + 5;
        let enforced_hiding_bound = Some(1);

        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(n).unwrap();

        // For the test to pass, each value in f_evals must be less than its corresponding value in
        // g_evals, mod n (since we use domain_h)

        let f_evals = vec![
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(3),
            domain_h.element(3),
        ];
        let g_evals = vec![
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
        ];

        let f_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&f_evals));
        let f_poly = LabeledPolynomial::new(
            String::from("f_poly"),
            f_poly,
            Some(other_degree_bound),
            enforced_hiding_bound,
        );
        let g_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&g_evals));
        let g_poly = LabeledPolynomial::new(
            String::from("g_poly"),
            g_poly,
            Some(other_degree_bound),
            enforced_hiding_bound,
        );

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            1,
            Some(&[2, enforced_degree_bound, other_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[f_poly.clone(), g_poly.clone()], Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let proof = DLComparison::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &domain_h,
            &f_poly,
            &commitments[0],
            &rands[0],
            &g_poly,
            &commitments[1],
            &rands[1],
            Some(other_degree_bound),
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let res = DLComparison::verify(
            &vk,
            &ck,
            &domain_k,
            &domain_h,
            &commitments[0],
            &commitments[1],
            Some(enforced_degree_bound),
            proof,
            &mut fs_rng,
        );

        assert!(res.is_err());

        // Test for a specific error
        assert_eq!(
            res.err().unwrap(),
            Error::ZeroOverKError(String::from("BatchCheckError"))
        );
    }
}
