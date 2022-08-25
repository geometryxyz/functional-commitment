#[cfg(test)]
mod tests {
    use crate::{error::Error, geo_seq::GeoSeqTest, util::generate_sequence};
    use ark_bn254::{Bn254, Fr};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::SimpleHashFiatShamirRng;
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    type F = Fr;
    type PC = KZG10<Bn254>;

    /// Test that geometric_sequence() works correctly
    #[test]
    fn test_generate_sequence_0() {
        let common_ratio = F::from(2u64);
        let sequence_initial_values = &[F::from(1u64), F::from(2u64)];
        let sequence_lengths = &[3, 3];
        let seq = generate_sequence(common_ratio, sequence_initial_values, sequence_lengths);

        let expected = [1, 2, 4, 2, 4, 8]
            .iter()
            .map(|x| F::from(*x as u64))
            .collect::<Vec<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }

        assert!(GeoSeqTest::<F, KZG10<Bn254>, FS>::naive_verify(
            &seq,
            common_ratio,
            sequence_initial_values,
            sequence_lengths
        ));
    }

    #[test]
    fn test_generate_sequence_1() {
        let common_ratio = F::from(1u64);
        let sequence_initial_values = &[F::from(1u64), F::from(1u64)];
        let sequence_lengths = &[1, 1];

        let seq = generate_sequence(common_ratio, sequence_initial_values, sequence_lengths);
        let expected = [1, 1]
            .iter()
            .map(|x| F::from(*x as u64))
            .collect::<Vec<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }
        assert!(GeoSeqTest::<F, KZG10<Bn254>, FS>::naive_verify(
            &seq,
            common_ratio,
            sequence_initial_values,
            sequence_lengths
        ));
    }

    #[test]
    fn test_geo_seq_proof() {
        let rng = &mut thread_rng();

        // define a sequence
        let common_ratio = Fr::from(9u64);
        let mut sequence_initial_values = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut sequence_lengths = vec![5, 3, 10, 30];

        // choose an appropriate domain for our sequence
        let m = sequence_lengths.iter().sum();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        // pad sequence to fit the domain
        let to_pad = domain_k.size() - m;
        if to_pad > 0 {
            sequence_initial_values.push(Fr::from(0u64));
            sequence_lengths.push(to_pad);
        }

        // generate the sequence
        let seq = generate_sequence::<F>(
            common_ratio,
            &sequence_initial_values.as_slice(),
            &sequence_lengths.as_slice(),
        );

        // Setup our polynomial commitment scheme
        let max_degree = 80;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = domain_k.size() + 1; // masking polynomials in zero over k have degree |K|+1 by definition

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        // Generate a polynomial from the sequence defined above
        let f = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let f = LabeledPolynomial::new(
            String::from("f"),
            f,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let (commitment, rands) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let proof = GeoSeqTest::<F, KZG10<Bn254>, FS>::prove(
            &ck,
            common_ratio,
            &f,
            &commitment[0].clone(),
            &rands[0].clone(),
            &mut sequence_initial_values,
            &mut sequence_lengths,
            &domain_k,
            rng,
        )
        .unwrap();

        let res = GeoSeqTest::<F, KZG10<Bn254>, FS>::verify(
            common_ratio,
            &sequence_initial_values,
            &sequence_lengths,
            &domain_k,
            &commitment[0],
            Some(enforced_degree_bound),
            proof,
            &vk,
        )
        .unwrap();

        assert_eq!((), res)
    }

    #[test]
    fn test_geo_seq_invalid() {
        let rng = &mut thread_rng();

        // define a sequence
        let common_ratio = Fr::from(9u64);
        let mut sequence_initial_values = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut sequence_lengths = vec![5, 3, 10, 30];

        // choose an appropriate domain for our sequence
        let m = sequence_lengths.iter().sum();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        // pad sequence to fit the domain
        let to_pad = domain_k.size() - m;
        if to_pad > 0 {
            sequence_initial_values.push(Fr::from(0u64));
            sequence_lengths.push(to_pad);
        }

        // generate the sequence
        let seq = generate_sequence::<F>(
            common_ratio,
            &sequence_initial_values.as_slice(),
            &sequence_lengths.as_slice(),
        );

        // Setup our polynomial commitment scheme
        let max_degree = 80;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = domain_k.size() + 1; // masking polynomials in zero over k have degree |K|+1 by definition

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        // Generate a polynomial from the sequence defined above
        let f = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let f = LabeledPolynomial::new(
            String::from("f"),
            f,
            Some(enforced_degree_bound),
            enforced_hiding_bound,
        );

        let (commitment, rands) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let wrong_common_ratio = Fr::from(2u64);

        let proof = GeoSeqTest::<F, KZG10<Bn254>, FS>::prove(
            &ck,
            wrong_common_ratio,
            &f,
            &commitment[0].clone(),
            &rands[0].clone(),
            &mut sequence_initial_values,
            &mut sequence_lengths,
            &domain_k,
            rng,
        )
        .unwrap();

        let res = GeoSeqTest::<F, KZG10<Bn254>, FS>::verify(
            common_ratio,
            &sequence_initial_values,
            &sequence_lengths,
            &domain_k,
            &commitment[0],
            Some(enforced_degree_bound),
            proof,
            &vk,
        );

        assert!(common_ratio != wrong_common_ratio);

        assert!(res.is_err());

        // Test for a specific error
        assert_eq!(res.err().unwrap(), Error::BatchCheckError);
    }

    #[test]
    fn test_reject_wrong_degree() {
        let rng = &mut thread_rng();

        // define a sequence
        let common_ratio = Fr::from(9u64);
        let mut sequence_initial_values = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut sequence_lengths = vec![5, 3, 10, 30];

        // choose an appropriate domain for our sequence
        let m = sequence_lengths.iter().sum();
        let domain_k = GeneralEvaluationDomain::<F>::new(m).unwrap();

        // pad sequence to fit the domain
        let to_pad = domain_k.size() - m;
        if to_pad > 0 {
            sequence_initial_values.push(Fr::from(0u64));
            sequence_lengths.push(to_pad);
        }

        // generate the sequence
        let seq = generate_sequence::<F>(
            common_ratio,
            &sequence_initial_values.as_slice(),
            &sequence_lengths.as_slice(),
        );

        // Setup our polynomial commitment scheme
        let max_degree = 80;
        let max_hiding = 1;

        let enforced_hiding_bound = Some(1);
        let enforced_degree_bound = domain_k.size() + 1; // masking polynomials in zero over k have degree |K|+1 by definition
        let other_degree = max_degree;

        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            max_hiding,
            Some(&[2, enforced_degree_bound, other_degree]),
        )
        .unwrap();

        // Generate a polynomial from the sequence defined above
        let f = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let f = LabeledPolynomial::new(
            String::from("f"),
            f,
            Some(other_degree),
            enforced_hiding_bound,
        );

        let (commitment, rands) = PC::commit(&ck, &[f.clone()], Some(rng)).unwrap();

        let proof = GeoSeqTest::<F, KZG10<Bn254>, FS>::prove(
            &ck,
            common_ratio,
            &f,
            &commitment[0].clone(),
            &rands[0].clone(),
            &mut sequence_initial_values,
            &mut sequence_lengths,
            &domain_k,
            rng,
        )
        .unwrap();

        let res = GeoSeqTest::<F, KZG10<Bn254>, FS>::verify(
            common_ratio,
            &sequence_initial_values,
            &sequence_lengths,
            &domain_k,
            &commitment[0],
            Some(enforced_degree_bound),
            proof,
            &vk,
        );

        assert!(res.is_err());

        assert_eq!(res.err().unwrap(), Error::BatchCheckError);
    }
}
