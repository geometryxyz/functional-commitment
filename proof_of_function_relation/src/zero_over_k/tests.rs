#[cfg(test)]
mod test {
    use crate::{
        commitment::KZG10,
        error::{to_pc_error, Error},
        label_polynomial,
        util::random_deg_n_polynomial,
        virtual_oracle::inverse_check_oracle::InverseCheckOracle,
        virtual_oracle::prod_vo::ProdVO,
        virtual_oracle::VirtualOracle,
        zero_over_k::ZeroOverK,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::One;
    use ark_ff::{Field, UniformRand};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
        Polynomial, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::{rand::thread_rng, test_rng};
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
        let rng = &mut test_rng();
        let domain_k = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let hiding_bound = Some(1);
        let degree_bound = 14;

        let max_degree = 20;
        let max_hiding = 1;
        let pp = PC::setup(max_degree, None, rng).unwrap();

        let (ck, vk) = PC::trim(&pp, max_degree, max_hiding, Some(&[2, degree_bound])).unwrap();

        // Step 1: choose a random polynomial
        let f_unlabeled: DensePolynomial<F> = random_deg_n_polynomial(7, rng);
        // let f_unlabeled =
        // DensePolynomial::<F>::from_coefficients_vec(vec![F::from(5), F::from(2), F::from(3)]);
        let f = LabeledPolynomial::new(
            String::from("f"),
            f_unlabeled,
            Some(degree_bound),
            hiding_bound,
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
            Some(degree_bound),
            hiding_bound,
        );

        // // PRINT G AND ITS DEGREE
        // println!("{:?}", g);
        // println!("degree of g: {}", g.degree());

        // // CHECK THAT G INDEED EVALUATES TO THE DESIRED VALUES
        // let g_evals = g.evaluate_over_domain_by_ref(domain_k);
        // g_evals
        //     .evals
        //     .iter()
        //     .zip(desired_g_evals.evals.iter())
        //     .for_each(|(&have, &want)| assert_eq!(have, want));

        // // CHECK THAT EVALUATE_OVER_DOMAIN DOES THE SAME AS EVALUATING AT EACH POINT
        // let g_manual_evals: Vec<_> = (0..n)
        //     .map(|index| g.evaluate(&domain_k.element(index)))
        //     .collect();
        // println!("mega dummy check");
        // g_manual_evals
        //     .iter()
        //     .zip(g_evals.evals.iter())
        //     .for_each(|(&manual, &g)| {
        //         assert_eq!(manual, g);
        //         println!("{}, {}", manual, g)
        //     });

        // // CHECK THAT THE PRODUCT OF EVALUATIONS OF F WITH EVALUATIONS OF G IS STILL ONE
        // let product_of_the_evaluations: Vec<_> = f
        //     .evaluate_over_domain_by_ref(domain_k)
        //     .evals
        //     .iter()
        //     .zip(g.evaluate_over_domain_by_ref(domain_k).evals.iter())
        //     .map(|(&f_eval, g_eval)| f_eval * g_eval)
        //     .collect();
        // println!("PRODUCT OF THE EVALUATIONS");
        // product_of_the_evaluations
        //     .iter()
        //     .for_each(|&p| println!("{p}"));

        // Step 5: Manually compute the product polynomial of f and g
        let f_times_g = LabeledPolynomial::new(
            String::from("f_times_g"),
            f.polynomial() * g.polynomial(),
            None,
            hiding_bound,
        );

        // // CHECK THAT THE EVALUTION OF THE PRODUCT MATCHES THE PRODUCT OF THE EVALUATIONS
        // let evaluations_of_product_poly = f_times_g.evaluate_over_domain_by_ref(domain_k);
        // let products_of_poly_evaluations: Vec<_> = f
        //     .evaluate_over_domain_by_ref(domain_k)
        //     .evals
        //     .iter()
        //     .zip(g.evaluate_over_domain_by_ref(domain_k).evals.iter())
        //     .map(|(&f_eval, &g_eval)| f_eval * g_eval)
        //     .collect();

        // evaluations_of_product_poly
        //     .evals
        //     .iter()
        //     .zip(products_of_poly_evaluations.iter())
        //     .for_each(|(eval_of_prod, prod_of_evals)| println!("{eval_of_prod}, {prod_of_evals}"));

        let concrete_oracles = [f, g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng))
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let zero_over_k_vo = InverseCheckOracle::new();

        let sanity_check = zero_over_k_vo
            .instantiate_in_coeffs_form(&concrete_oracles, &alphas)
            .unwrap();
        let check_evals = sanity_check.evaluate_over_domain(domain_k);
        check_evals.evals.iter().for_each(|eval| println!("{eval}"));

        let zero_over_k_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            Some(degree_bound),
            &zero_over_k_vo,
            &alphas,
            &domain_k,
            &ck,
            rng,
        );

        let is_valid = ZeroOverK::<F, PC, D>::verify(
            zero_over_k_proof.unwrap(),
            &commitments,
            Some(degree_bound),
            &zero_over_k_vo,
            &domain_k,
            &alphas,
            &vk,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_zero_over_k_prod_vo() {
        let n = 8;
        let rng = &mut thread_rng();
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let max_degree = 20;
        let degree_bound = 14;
        let hiding_bound = Some(1);
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(&pp, max_degree, 1, Some(&[2, 14])).unwrap();

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
        let f = LabeledPolynomial::new(String::from("g"), f, Some(degree_bound), hiding_bound);

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
        let g = LabeledPolynomial::new(
            String::from("g"),
            g.clone(),
            Some(degree_bound),
            hiding_bound,
        );

        let concrete_oracles = [f, g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng))
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let zero_over_k_vo = ProdVO {};

        let zero_over_k_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            Some(degree_bound),
            &zero_over_k_vo,
            &alphas,
            &domain,
            &ck,
            rng,
        )
        .unwrap();

        let is_valid = ZeroOverK::<F, PC, D>::verify(
            zero_over_k_proof,
            &commitments,
            Some(degree_bound),
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
