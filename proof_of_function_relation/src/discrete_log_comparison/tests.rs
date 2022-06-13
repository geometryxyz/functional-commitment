#[cfg(test)]
mod tests {
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, FftField, Field, SquareRootField};
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

    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        discrete_log_comparison::DLComparison,
        label_polynomial,
    };

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

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

        println!("{}", sq_root);

        assert_eq!(generator, sq_root * sq_root);
        assert_eq!(g2, sq_root2 * sq_root2);
    }

    #[test]
    fn test_discrete_log_proof() {
        let mut rng = thread_rng();
        let m = 8;
        let n = 4;

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
        let f_poly = label_polynomial!(f_poly);
        let g_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&g_evals));
        let g_poly = label_polynomial!(g_poly);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[f_poly.clone(), g_poly.clone()], Some(&mut rng)).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = DLComparison::<F, PC, D>::prove(
            &ck,
            &domain_k,
            &domain_h,
            &f_poly,
            &g_poly,
            &commitments[0],
            &commitments[1],
            &mut fs_rng,
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            DLComparison::verify(
                &vk,
                &ck,
                &domain_k,
                &domain_h,
                &commitments[0],
                &commitments[1],
                proof,
                &mut fs_rng,
            ).is_ok(), 
            true 
        );
    }

    #[test]
    #[should_panic]
    fn test_malicious_discrete_log_proof() {
        let mut rng = thread_rng();
        let m = 8;
        let n = 4;

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
            domain_h.element(3),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
            domain_h.element(0),
            domain_h.element(1),
            domain_h.element(1),
            domain_h.element(2),
        ];

        let f_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&f_evals));
        let f_poly = label_polynomial!(f_poly);
        let g_poly = DensePolynomial::<F>::from_coefficients_vec(domain_k.ifft(&g_evals));
        let g_poly = label_polynomial!(g_poly);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (commitments, rands) =
            PC::commit(&ck, &[f_poly.clone(), g_poly.clone()], Some(&mut rng)).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

        let proof = DLComparison::<F, PC, D>::prove(
            &ck,
            &domain_k,
            &domain_h,
            &f_poly,
            &g_poly,
            &commitments[0],
            &commitments[1],
            &mut fs_rng,
            &mut rng,
        )
        .unwrap();

        DLComparison::verify(
            &vk,
            &ck,
            &domain_k,
            &domain_h,
            &commitments[0],
            &commitments[1],
            proof,
            &mut fs_rng,
        )
        .unwrap();
    }
}
