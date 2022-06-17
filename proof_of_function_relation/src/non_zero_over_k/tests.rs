#[cfg(test)]
mod test {
    use crate::{
        commitment::KZG10,
        label_polynomial,
        non_zero_over_k::NonZeroOverK,
        util::sample_vector,
        error::Error,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::Zero;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
        UVPolynomial,
    };
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_zero_over_k_for_non_zero_over_k() {
        let mut rng = thread_rng();
        let n = 8;

        let f_evals: Vec<F> = sample_vector(&mut rng, n);

        for &eval in f_evals.iter() {
            assert_ne!(eval, F::zero())
        }

        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&f_evals));
        let f = label_polynomial!(f);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (f_commit, _) = PC::commit(&ck, &[f.clone()], Some(&mut rng)).unwrap();

        let proof = NonZeroOverK::<F, PC, Blake2s>::prove(&ck, &domain, f, &mut rng).unwrap();

        assert_eq!(
            true,
            NonZeroOverK::<F, PC, Blake2s>::verify(&vk, &domain, f_commit[0].clone(), proof)
                .is_ok()
        );
    }

    #[test]
    fn rand_failing_test() {
        let mut rng = thread_rng();
        let n = 8;

        let mut f_evals: Vec<F> = sample_vector(&mut rng, n);
        f_evals[4] = F::zero();

        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&f_evals));
        let f = label_polynomial!(f);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, _) = PC::trim(&pp, max_degree, 0, None).unwrap();

        // panic here because there is no inverse of 0
        let proof = NonZeroOverK::<F, PC, Blake2s>::prove(&ck, &domain, f, &mut rng);

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(
            proof.err().unwrap(),
            Error::FEvalIsZero
        );

    }
}
