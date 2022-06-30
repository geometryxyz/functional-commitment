#[cfg(test)]
mod test {
    use crate::{
        commitment::KZG10,
        error::{to_pc_error, Error},
        label_polynomial,
        virtual_oracle::inverse_check_oracle::InverseCheckOracle,
        virtual_oracle::prod_vo::ProdVO,
        zero_over_k::ZeroOverK,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::Field;
    use ark_ff::One;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;

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

        let zero_over_k_vo = InverseCheckOracle::new();

        // let concrete_oracle_labels: Vec<_> =
        //     concrete_oracles.iter().map(|f| f.label().clone()).collect();
        // zero_over_k_vo
        //     .get_h_labels(&concrete_oracle_labels)
        //     .iter()
        //     .for_each(|label| println!("{label}"));

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

        // assert!(zero_over_k_proof.is_ok());

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
}
