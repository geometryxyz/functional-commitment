#[cfg(test)]
mod test {
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        error::Error,
        virtual_oracle::{Description, EvaluationsProvider, Term, VirtualOracle},
        zero_over_k::ZeroOverK,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{One, Zero};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_poly_commit::{
        LabeledCommitment, LabeledPolynomial, PCRandomness, PolynomialCommitment,
    };
    use ark_std::UniformRand;
    use blake2::Blake2s;
    use rand::thread_rng;
    use rand_core::OsRng;
    use std::iter;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_zero_over_k() {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let alphas = [domain.element(1), F::one()];
        // let alphas = [F::one(), F::one()];

        //a_evals (2, 4, 6, 8) -> f(w * x) = (4, 6, 8, 2)
        //b_evals (-1, -2, -3, 0) -> 2 * g(x) = (-2, -4, -6, 0)

        //test oracle to be zero at roots of unity
        let a_evals = vec![
            F::from(2u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let a_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&a_evals));
        let a_poly = LabeledPolynomial::new(String::from("a"), a_poly, None, None);

        let b_evals = vec![
            -F::from(1u64),
            -F::from(2u64),
            -F::from(3u64),
            -F::from(0u64),
        ];
        let b_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&b_evals));
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, None, None);

        let concrete_oracles = [a_poly.clone(), b_poly.clone()];

        // concrete_oracles = [f, g]
        // alpha_coeffs = [alpha_1, alpha_2]
        //1*f(alpha_1X) + 2*g(alpha_2X) - 2
        let term0 = Term {
            concrete_oracle_indices: vec![0],
            alpha_coeff_indices: vec![0],
            constant: Fr::from(1 as u64),
        };
        let term1 = Term {
            concrete_oracle_indices: vec![1],
            alpha_coeff_indices: vec![1],
            constant: Fr::from(2 as u64),
        };

        let term2 = Term {
            concrete_oracle_indices: vec![],
            alpha_coeff_indices: vec![],
            constant: -F::from(2u64),
        };

        let description = Description::<Fr> {
            terms: vec![term0, term1, term2],
        };
        let vo = VirtualOracle::new(description).unwrap();

        let concrete_oracles2 = [a_poly.clone(), b_poly.clone()].to_vec();

        let eval2 = concrete_oracles2
            .evaluate(&vo, domain.element(1), &alphas.to_vec())
            .unwrap();
        assert_eq!(eval2, F::zero());

        // The proof of zero over K
        let maximum_degree: usize = 30;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, Some(&[2, 5])).unwrap();

        let (concrete_oracles_commitments, concrete_oracle_rands) =
            PC::commit(&ck, &concrete_oracles, None).unwrap();

        let proof = ZeroOverK::<F, KZG10<Bn254>, Blake2s>::prove(
            &concrete_oracles,
            &concrete_oracles_commitments,
            &concrete_oracle_rands,
            &vo,
            &alphas.to_vec(),
            domain,
            &ck,
            &mut OsRng,
        )
        .unwrap();

        assert_eq!(
            true,
            ZeroOverK::<F, KZG10<Bn254>, Blake2s>::verify(
                proof,
                &concrete_oracles_commitments,
                &vo,
                domain,
                &alphas,
                &vk,
            )
            .is_ok()
        );
    }

    #[test]
    fn test_failure_on_malicious_virtual_oracle() {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let alphas = [F::one(), F::one()];
        // let alphas = [F::one(), F::one()];

        //a_evals (2, 4, 6, 8) -> f(w * x) = (4, 6, 8, 2)
        //b_evals (-1, -2, -3, 0) -> 2 * g(x) = (-2, -4, -6, 0)

        //test oracle to be zero at roots of unity
        let a_evals = vec![
            F::from(2u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let a_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&a_evals));
        let a_poly = LabeledPolynomial::new(String::from("b"), a_poly, None, None);

        let b_evals = vec![
            -F::from(1u64),
            -F::from(2u64),
            -F::from(3u64),
            -F::from(0u64),
        ];
        let b_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&b_evals));
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, None, None);

        let concrete_oracles = [a_poly.clone(), b_poly.clone()];

        // let instantiated_virtual_oracle =
        //     TestVirtualOracle::instantiate(&concrete_oracles, &alphas);

        // let _ = instantiated_virtual_oracle.evaluate(&domain.element(1));
        //
        // concrete_oracles = [f, g]
        // alpha_coeffs = [alpha_1, alpha_2]
        let term0 = Term {
            concrete_oracle_indices: vec![0],
            alpha_coeff_indices: vec![0],
            constant: Fr::from(1 as u64),
        };
        let term1 = Term {
            concrete_oracle_indices: vec![1],
            alpha_coeff_indices: vec![1],
            constant: Fr::from(2 as u64),
        };

        let term2 = Term {
            concrete_oracle_indices: vec![],
            alpha_coeff_indices: vec![],
            constant: -F::from(2u64),
        };

        let description = Description::<Fr> {
            terms: vec![term0, term1, term2],
        };
        let vo = VirtualOracle::new(description).unwrap();

        let maximum_degree: usize = 16;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, None).unwrap();

        // println!("vk: {}", vk);

        let (concrete_oracles_commitments, concrete_oracle_rands) =
            PC::commit(&ck, &concrete_oracles, None).unwrap();

        let proof = ZeroOverK::<F, KZG10<Bn254>, Blake2s>::prove(
            &concrete_oracles,
            &concrete_oracles_commitments,
            &concrete_oracle_rands,
            &vo,
            &alphas.to_vec(),
            domain,
            &ck,
            &mut OsRng,
        )
        .unwrap();

        let verification_result = ZeroOverK::<F, KZG10<Bn254>, Blake2s>::verify(
            proof,
            &concrete_oracles_commitments,
            &vo,
            domain,
            &alphas,
            &vk,
        );

        assert_eq!(Err(Error::Check2Failed), verification_result);
    }

    #[test]
    fn test_commit_with_bounds() {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let rng = &mut thread_rng();

        //test oracle to be zero at roots of unity
        let a_evals = vec![
            F::from(2u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let a_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&a_evals));
        let a_poly = LabeledPolynomial::new(String::from("a"), a_poly, Some(4), None);

        let b_evals = vec![
            F::from(2u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let b_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&b_evals));
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, Some(4), None);

        let a_plus_b_poly = a_poly.polynomial() + b_poly.polynomial();
        let a_plus_b_poly =
            LabeledPolynomial::new(String::from("a_plus_b"), a_plus_b_poly, Some(4), None);

        let point = Fr::from(10u64);
        let eval = a_plus_b_poly.evaluate(&point);

        let maximum_degree: usize = 16;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, Some(&[4])).unwrap();

        let (commitments, rands) = PC::commit(&ck, &[a_poly, b_poly], None).unwrap();
        let a_commit = commitments[0].clone();
        let a_rand = rands[0].clone();

        let b_commit = commitments[1].clone();
        let b_rand = rands[1].clone();

        let one = Fr::one();
        let a_plus_b_commit = PC::multi_scalar_mul(&[a_commit, b_commit], &[one, one]);
        let a_plus_b_commit =
            LabeledCommitment::new(String::from("a_plus_b"), a_plus_b_commit, Some(4));

        let challenge = Fr::from(1u64);

        let homomorphic_randomness =
            ark_poly_commit::kzg10::Randomness::<Fr, DensePolynomial<Fr>>::empty();

        let opening = PC::open(
            &ck,
            &[a_plus_b_poly.clone()],
            &[a_plus_b_commit.clone()],
            &point,
            challenge,
            &[homomorphic_randomness],
            None,
        )
        .unwrap();

        let res = PC::check(
            &vk,
            &[a_plus_b_commit],
            &point,
            iter::once(eval),
            &opening,
            challenge,
            None,
        );

        println!("{:?}", res);

        // let (commitments, randoms) = PC::commit(&ck, &[a_poly_no_bound, a_poly_with_bound], None).unwrap();

        // for c in commitments {
        //     println!("{:?}", c.commitment());
        // }
    }
}
