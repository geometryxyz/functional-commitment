#[cfg(test)]
mod test {
    use crate::virtual_oracle::VirtualOracle;
    use crate::{
        commitment::KZG10, error::Error, label_polynomial, non_zero_over_k::NonZeroOverK,
        util::sample_vector, virtual_oracle::inverse_check_oracle::InverseCheckOracle,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::Field;
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

    type F = Fr;
    type PC = KZG10<Bn254>;

    // This test should pass because it the randomly generated virtual oracle should
    // not evalute to 0 - if it does, just rerun the test
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

        let proof = NonZeroOverK::<F, PC, Blake2s>::prove(&ck, &domain, &f, &mut rng).unwrap();

        assert_eq!(
            true,
            NonZeroOverK::<F, PC, Blake2s>::verify(&vk, &domain, f_commit[0].clone(), proof)
                .is_ok()
        );
    }

    // This test will fail because one of the f_evals is 0, which has no inverse.
    #[test]
    fn test_f_eval_is_zero() {
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

        // This proof cannot be generated because there is no inverse of 0
        let proof = NonZeroOverK::<F, PC, Blake2s>::prove(&ck, &domain, &f, &mut rng);

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    // This test uses a virtual oracle that evalutes to zero over k.
    // Therefore, non_zero_over_k should fail.
    #[test]
    fn test_using_zero_over_k_vo() {
        let mut rng = thread_rng();
        let n = 8;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let a_evals: Vec<F> = vec![
            F::from(1u64),
            F::from(2u64),
            F::from(3u64),
            F::from(4u64),
            F::from(5u64),
            F::from(6u64),
            F::from(7u64),
            F::from(8u64),
        ];

        let a = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&a_evals));
        let a = label_polynomial!(a);

        let a_evals = domain.fft(a.coeffs());

        let b_evals = a_evals
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        let b = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&b_evals));
        let b = LabeledPolynomial::new(String::from("b"), b.clone(), None, None);

        let concrete_oracles = [a, b];
        let alphas = vec![F::one(), F::one()];

        let zero_over_k_vo = InverseCheckOracle::new();

        let f = zero_over_k_vo
            .instantiate_in_coeffs_form(&concrete_oracles, alphas.as_slice())
            .unwrap();
        let f = LabeledPolynomial::new(String::from("f"), f.clone(), None, None);

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

        let (f_commit, _) = PC::commit(&ck, &[f.clone()], Some(&mut rng)).unwrap();

        let proof = NonZeroOverK::<F, PC, Blake2s>::prove(&ck, &domain, &f, &mut rng).unwrap();
        let is_valid =
            NonZeroOverK::<F, PC, Blake2s>::verify(&vk, &domain, f_commit[0].clone(), proof);

        assert!(is_valid.is_err());
        // Test for a specific error
        assert_eq!(is_valid.err().unwrap(), Error::Check2Failed);
    }
}
