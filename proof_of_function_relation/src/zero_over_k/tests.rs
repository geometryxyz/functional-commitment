#[cfg(test)]
mod test {
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        error::{to_pc_error, Error},
        label_polynomial,
        virtual_oracle::inverse_check_oracle::InverseCheckOracle,
        virtual_oracle::prod_vo::ProdVO,
        zero_over_k::ZeroOverK,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::One;
    use ark_ff::{Field, UniformRand};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
        UVPolynomial,
    };
    use ark_poly_commit::{
        LabeledCommitment, LabeledPolynomial, PCRandomness, PolynomialCommitment,
    };
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use rand_core::OsRng;
    use std::iter;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type D = Blake2s;

    // TODO
    // Write more tests. e.g. non_zero_over_k uses zero_over_k with an InverseCheckOracle.
    // Try other virtual oracles like add_vo.

    #[test]
    fn test_zero_over_k_inverse_check_oracle() {
        let n = 8;
        let mut rng = thread_rng();
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

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

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&f_evals));
        let f = label_polynomial!(f);

        let f_evals = domain.fft(f.coeffs());

        let g_evals = f_evals
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        let g = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&g_evals));
        let g = LabeledPolynomial::new(String::from("g"), g.clone(), None, None);

        let concrete_oracles = [f, g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, None)
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let zero_over_k_vo = InverseCheckOracle {};

        let zero_over_k_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            &zero_over_k_vo,
            &alphas,
            &domain,
            &ck,
            &mut rng,
        );

        assert!(zero_over_k_proof.is_ok());

        let is_valid = ZeroOverK::<F, PC, D>::verify(
            zero_over_k_proof.unwrap(),
            &commitments,
            &zero_over_k_vo,
            &domain,
            &alphas,
            &vk,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_zero_over_k_prod_vo() {
        let n = 8;
        let mut rng = thread_rng();
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, &mut rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

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

        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&f_evals));
        let f = label_polynomial!(f);

        let g_evals: Vec<F> = vec![
            F::from(1u64),
            F::from(2u64),
            F::from(3u64),
            F::from(4u64),
            F::from(5u64),
            F::from(6u64),
            F::from(7u64),
            F::from(8u64),
        ];

        let g = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&g_evals));
        let g = LabeledPolynomial::new(String::from("g"), g.clone(), None, None);

        let concrete_oracles = [f, g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, None)
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let zero_over_k_vo = ProdVO {};

        let zero_over_k_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            &zero_over_k_vo,
            &alphas,
            &domain,
            &ck,
            &mut rng,
        );

        assert!(zero_over_k_proof.is_ok());

        let is_valid = ZeroOverK::<F, PC, D>::verify(
            zero_over_k_proof.unwrap(),
            &commitments,
            &zero_over_k_vo,
            &domain,
            &alphas,
            &vk,
        );

        assert!(is_valid.is_err());

        // Test for a specific error
        assert_eq!(is_valid.err().unwrap(), Error::Check2Failed);
    }

    #[test]
    fn test_commit_with_bounds() {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let rng = &mut thread_rng();

        //test oracle to be zero at roots of unity
        let a_evals = vec![F::from(2u64), F::from(4u64)];
        let a_poly = Evaluations::from_vec_and_domain(a_evals, domain).interpolate();
        let a_poly = LabeledPolynomial::new(String::from("a"), a_poly, Some(4), None);

        let b_evals = vec![
            F::from(5u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let b_poly = Evaluations::from_vec_and_domain(b_evals, domain).interpolate();
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, Some(4), None);

        let a_plus_b_poly = a_poly.polynomial() + b_poly.polynomial();
        let a_plus_b_poly =
            LabeledPolynomial::new(String::from("a_plus_b"), a_plus_b_poly, Some(4), None);

        let point = Fr::rand(rng);
        let eval = a_plus_b_poly.evaluate(&point);

        let maximum_degree: usize = 16;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, Some(&[4])).unwrap();

        let (commitments, _rands) = PC::commit(&ck, &[a_poly, b_poly], None).unwrap();
        let a_commit = commitments[0].clone();
        //let a_rand = rands[0].clone();

        let b_commit = commitments[1].clone();
        // let b_rand = rands[1].clone();

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
