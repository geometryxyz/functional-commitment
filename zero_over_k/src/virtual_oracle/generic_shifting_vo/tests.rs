#[cfg(test)]
mod test {
    use crate::util::sample_vector;
    use crate::virtual_oracle::generic_shifting_vo::presets;
    use crate::virtual_oracle::VirtualOracle;
    use crate::{error::to_pc_error, zero_over_k::ZeroOverK};
    use crate::{error::Error, util::shift_dense_poly, vo_constant};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{Field, One, UniformRand};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
        Polynomial, UVPolynomial,
    };
    use ark_poly_commit::{evaluate_query_set, LabeledPolynomial, PolynomialCommitment};
    use ark_std::{rand::thread_rng, test_rng};
    use blake2::Blake2s;
    use fiat_shamir_rng::SimpleHashFiatShamirRng;
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;

    use super::super::{GenericShiftingVO, VOTerm};

    type F = Fr;
    type PC = KZG10<Bn254>;
    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;

    pub fn simple_addition(terms: &[VOTerm<F>]) -> VOTerm<F> {
        terms[1].clone() + vo_constant!(F::from(2u64)) * terms[2].clone()
    }

    pub fn simple_mul(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
        concrete_terms[1].clone() * concrete_terms[2].clone()
    }

    pub fn harder_addition(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
        concrete_terms[0].clone()
            + vo_constant!(F::from(2u64)) * concrete_terms[2].clone()
            + vo_constant!(F::from(3u64)) * concrete_terms[3].clone()
    }

    #[test]
    fn test_add_oracle() {
        let rng = &mut thread_rng();
        let a_poly = LabeledPolynomial::new(
            String::from("a"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let b_poly = LabeledPolynomial::new(
            String::from("b"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let c_poly = LabeledPolynomial::new(
            String::from("c"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let two_constant_poly = DensePolynomial::from_coefficients_vec(vec![F::from(2u64)]);
        let expected = &shifted_c + &(&two_constant_poly * &shifted_a);

        let add_oracle =
            GenericShiftingVO::new(&mapping_vector, &shifting_coefficients, simple_addition)
                .unwrap();

        // Check that we get the right polynomial
        let sum = add_oracle.compute_polynomial(concrete_oracles).unwrap();
        assert_eq!(expected, sum);
    }

    #[test]
    fn test_mul_oracle() {
        let rng = &mut thread_rng();
        let a_poly = LabeledPolynomial::new(
            String::from("a"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let b_poly = LabeledPolynomial::new(
            String::from("b"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let c_poly = LabeledPolynomial::new(
            String::from("c"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );

        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let expected = &shifted_c * &shifted_a;

        let mul_oracle =
            GenericShiftingVO::new(&mapping_vector, &shifting_coefficients, simple_mul).unwrap();

        // Check that we get the right polynomial
        let prod = mul_oracle.compute_polynomial(concrete_oracles).unwrap();
        assert_eq!(expected, prod);
    }

    #[test]
    fn test_short_input_vec() {
        // mapping vector expects there to be a concrete oracle with index 1; effectively expected at last 2 concrete oracles
        let mapping_vector = vec![2];
        let shift_coefficients = vec![F::one()];
        let add_oracle =
            GenericShiftingVO::new(&mapping_vector, &shift_coefficients, simple_addition).unwrap();

        // We only provide one concrete oracle
        let default_labeled = LabeledPolynomial::new(
            String::from(""),
            DensePolynomial::<F>::default(),
            None,
            None,
        );
        let err_poly = add_oracle.compute_polynomial(&vec![default_labeled]);
        assert!(err_poly.is_err());
        assert_eq!(
            err_poly.unwrap_err(),
            Error::InputLengthError(String::from(
                "Mapping vector requires 2 oracles/evaluations but only 1 were provided"
            ))
        );
    }

    #[test]
    fn test_query_set() {
        let rng = &mut thread_rng();
        let a_poly = LabeledPolynomial::new(
            String::from("a"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let b_poly = LabeledPolynomial::new(
            String::from("b"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let c_poly = LabeledPolynomial::new(
            String::from("c"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );

        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];
        let oracle_labels: Vec<_> = concrete_oracles
            .iter()
            .map(|oracle| oracle.label().clone())
            .collect();

        let shifting_coefficients: Vec<F> = sample_vector(rng, 3);
        let mapping_vector = vec![2, 2, 0];

        let add_oracle =
            GenericShiftingVO::new(&mapping_vector, &shifting_coefficients, harder_addition)
                .unwrap();

        let eval_point = (String::from("beta"), F::rand(rng));
        let query_set = add_oracle
            .generate_query_set(&oracle_labels, &eval_point)
            .unwrap();

        let evals = evaluate_query_set(concrete_oracles.iter(), &query_set);

        let evaluated = add_oracle
            .evaluate_from_concrete_evals(&oracle_labels, &eval_point.1, &evals)
            .unwrap();
        let eval_from_poly = add_oracle
            .compute_polynomial(concrete_oracles)
            .unwrap()
            .evaluate(&eval_point.1);

        assert_eq!(evaluated, eval_from_poly)
    }

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
            GenericShiftingVO::new(&vec![0, 1], &alphas.clone(), presets::inverse_check).unwrap();

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
}
