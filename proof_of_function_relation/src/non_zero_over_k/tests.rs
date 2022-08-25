#[cfg(test)]
mod test {
    use crate::{error::Error, non_zero_over_k::NonZeroOverK};
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
    use fiat_shamir_rng::SimpleHashFiatShamirRng;
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;
    use zero_over_k::{
        util::sample_vector,
        virtual_oracle::generic_shifting_vo::{presets, GenericShiftingVO},
    };

    type F = Fr;
    type PC = KZG10<Bn254>;
    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;

    // This test should pass because it the randomly generated virtual oracle should
    // not evalute to 0 - if it does, just rerun the test
    #[test]
    fn test_non_zero_over_k() {
        let m = 8;
        let rng = &mut thread_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        let max_degree = 20;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = 14;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let f_unlabeled: DensePolynomial<F> = DensePolynomial::rand(7, rng);
        let f = LabeledPolynomial::new(
            String::from("f"),
            f_unlabeled,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let (f_commit, f_rand) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let proof = NonZeroOverK::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &f,
            &f_commit[0].clone(),
            &f_rand[0].clone(),
            rng,
        )
        .unwrap();

        let res = NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            f_commit[0].commitment().clone(),
            Some(enforced_degree_bound),
            proof,
        )
        .unwrap();

        assert_eq!(res, ());
    }

    // This test will fail because one of the f_evals is 0. The prover will fail to run.
    #[test]
    fn test_f_eval_is_zero() {
        let m = 8;
        let rng = &mut thread_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        let max_degree = 20;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = 14;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, _vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        // choose a random function f and set one of its evaluations in K to 0
        let mut f_evals: Vec<F> = sample_vector(rng, m);
        f_evals[4] = F::zero();
        let f = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&f_evals));
        let f = LabeledPolynomial::new(
            String::from("f"),
            f,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let (f_commit, f_rand) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let proof = NonZeroOverK::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &f,
            &f_commit[0].clone(),
            &f_rand[0].clone(),
            rng,
        );

        assert!(proof.is_err());

        // Test for a specific error
        assert_eq!(proof.err().unwrap(), Error::FEvalIsZero);
    }

    // This test uses a virtual oracle that evalutes to zero over k.
    // Therefore, non_zero_over_k should fail.
    #[test]
    fn test_using_zero_over_k_vo() {
        let m = 8;
        let rng = &mut thread_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        let max_degree = 20;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = 14;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

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

        let a = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&a_evals));
        let a = LabeledPolynomial::new(
            String::from("a"),
            a,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let a_evals = domain_k.fft(a.coeffs());

        let b_evals = a_evals
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        let b = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&b_evals));
        let b = LabeledPolynomial::new(
            String::from("b"),
            b.clone(),
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let concrete_oracles = [a, b];
        let alphas = vec![F::one(), F::one()];

        let zero_over_k_vo =
            GenericShiftingVO::new(&vec![0, 1], &alphas, presets::inverse_check).unwrap();

        let f = zero_over_k_vo
            .compute_polynomial(&concrete_oracles)
            .unwrap();
        let f = LabeledPolynomial::new(
            String::from("f"),
            f.clone(),
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let (f_commit, f_rand) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let proof = NonZeroOverK::<F, PC, FS>::prove(
            &ck,
            &domain_k,
            &f,
            &f_commit[0].clone(),
            &f_rand[0].clone(),
            rng,
        )
        .unwrap();
        let res = NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            f_commit[0].commitment().clone(),
            Some(enforced_degree_bound),
            proof,
        );

        assert!(res.is_err());
        // Test for a specific error
        assert_eq!(
            res.err().unwrap(),
            Error::ZeroOverKError(String::from("Check2Failed"))
        );
    }
}
