#[cfg(test)]
mod test {
    use crate::{
        error::{to_pc_error, Error},
        virtual_oracle::generic_shifting_vo::{presets, GenericShiftingVO},
        zero_over_k::ZeroOverK,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::Field;
    use ark_ff::One;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
        UVPolynomial,
    };
    use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, PolynomialCommitment};
    use ark_std::{rand::thread_rng, test_rng};
    use blake2::Blake2s;
    use fiat_shamir_rng::SimpleHashFiatShamirRng;
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;
    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_zero_over_k_inverse_check_oracle() {
        let m = 8;
        let rng = &mut test_rng();
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

        // Step 1: choose a random polynomial
        let f_unlabeled: DensePolynomial<F> = DensePolynomial::rand(7, rng);
        let f = LabeledPolynomial::new(
            String::from("f"),
            f_unlabeled,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        // Step 2: evaluate it
        let f_evals = f.evaluate_over_domain_by_ref(domain_k);

        // Step 3: find the inverse at each of these points
        let desired_g_evals = f_evals
            .evals
            .iter()
            .map(|&x| x.inverse().unwrap())
            .collect::<Vec<_>>();
        let desired_g_evals = Evaluations::from_vec_and_domain(desired_g_evals, domain_k);

        // Step 4: interpolate a polynomial from the inverses
        let g = desired_g_evals.clone().interpolate();
        let g = LabeledPolynomial::new(
            String::from("g"),
            g.clone(),
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        // Step 5: commit to the concrete oracles
        let concrete_oracles = [f, g];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng))
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        // Step 6: Derive the desired virtual oracle
        let alphas = vec![F::one(), F::one()];
        let inverse_check_oracle =
            GenericShiftingVO::new(&vec![0, 1], &alphas, presets::inverse_check).unwrap();

        // Step 7: prove
        let zero_over_k_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            Some(enforced_degree_bound),
            &inverse_check_oracle,
            &domain_k,
            &ck,
            rng,
        );

        // Step 8: verify
        let is_valid = ZeroOverK::<F, PC, FS>::verify(
            zero_over_k_proof.unwrap(),
            &commitments,
            Some(enforced_degree_bound),
            &inverse_check_oracle,
            &domain_k,
            &vk,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_error_on_invalid_proof() {
        let m = 8;
        let rng = &mut thread_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        let max_degree = 20;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = 14;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, max_hiding, Some(&[2, 14])).unwrap();

        let f_evals: Vec<F> = vec![
            F::from(1u64),
            F::from(2u64),
            F::from(3u64),
            F::from(4u64),
            F::from(5u64),
            F::from(6u64),
            F::from(7u64),
            F::from(8u64),
        ];

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&f_evals));
        let f = LabeledPolynomial::new(
            String::from("f"),
            f,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let g_evals: Vec<F> = vec![
            F::from(8u64),
            F::from(7u64),
            F::from(6u64),
            F::from(5u64),
            F::from(4u64),
            F::from(3u64),
            F::from(2u64),
            F::from(1u64),
        ];

        let g = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&g_evals));
        let g = LabeledPolynomial::new(
            String::from("g"),
            g.clone(),
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let concrete_oracles = [f, g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng))
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let zero_over_k_vo =
            GenericShiftingVO::new(&[0, 1], &alphas, presets::equality_check).unwrap();

        let zero_over_k_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            Some(enforced_degree_bound),
            &zero_over_k_vo,
            &domain_k,
            &ck,
            rng,
        )
        .unwrap();

        let res = ZeroOverK::<F, PC, FS>::verify(
            zero_over_k_proof,
            &commitments,
            Some(enforced_degree_bound),
            &zero_over_k_vo,
            &domain_k,
            &vk,
        );

        assert!(res.is_err());

        // Test for a specific error
        assert_eq!(res.unwrap_err(), Error::Check2Failed);
    }

    #[test]
    fn test_degree_bound_not_respected() {
        let m = 8;
        let rng = &mut test_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        let max_degree = 20;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let prover_degree_bound = 15;
        let verifier_degree_bound = 14;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, verifier_degree_bound, prover_degree_bound]),
        )
        .unwrap();

        // Step 1: choose a random polynomial
        let f_unlabeled: DensePolynomial<F> = DensePolynomial::rand(7, rng);
        let f = LabeledPolynomial::new(
            String::from("f"),
            f_unlabeled,
            Some(prover_degree_bound),
            enforced_hiding_bound,
        );

        // Step 2: evaluate it
        let f_evals = f.evaluate_over_domain_by_ref(domain_k);

        // Step 3: find the inverse at each of these points
        let desired_g_evals = f_evals
            .evals
            .iter()
            .map(|&x| x.inverse().unwrap())
            .collect::<Vec<_>>();
        let desired_g_evals = Evaluations::from_vec_and_domain(desired_g_evals, domain_k);

        // Step 4: interpolate a polynomial from the inverses
        let g = desired_g_evals.clone().interpolate();
        let g = LabeledPolynomial::new(
            String::from("g"),
            g.clone(),
            Some(prover_degree_bound),
            enforced_hiding_bound,
        );

        // Step 5: commit to the concrete oracles
        let concrete_oracles = [f, g];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng))
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        // Step 6: Derive the desired virtual oracle
        let alphas = vec![F::one(), F::one()];
        let inverse_check_oracle =
            GenericShiftingVO::new(&vec![0, 1], &alphas, presets::inverse_check).unwrap();

        // Step 7: prove
        let zero_over_k_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            Some(prover_degree_bound),
            &inverse_check_oracle,
            &domain_k,
            &ck,
            rng,
        )
        .unwrap();

        // Step 8: verify
        let verifier_commitments: Vec<LabeledCommitment<_>> = commitments
            .iter()
            .map(|c| {
                let label = c.label().clone();
                let comm = c.commitment().clone();
                LabeledCommitment::new(label, comm, Some(verifier_degree_bound))
            })
            .collect();
        let res = ZeroOverK::<F, PC, FS>::verify(
            zero_over_k_proof,
            &verifier_commitments,
            Some(verifier_degree_bound),
            &inverse_check_oracle,
            &domain_k,
            &vk,
        );

        assert!(res.is_err());

        assert_eq!(res.unwrap_err(), Error::BatchCheckError);
    }
}
