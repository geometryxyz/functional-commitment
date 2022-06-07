#[cfg(test)]
mod tests {
    use crate::virtual_oracle::{geometric_sequence_vo::GeoSequenceVO, VirtualOracle};
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        geo_seq::GeoSeqTest,
    };
    use crate::{label_polynomial, util::generate_sequence, zero_over_k::ZeroOverK};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::PrimeField;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_poly_commit::{
        LabeledCommitment, LabeledPolynomial, PCRandomness, PolynomialCommitment,
    };
    use blake2::Blake2s;
    use rand_core::OsRng;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_geo_seq() {
        let r = F::from(2u64);
        let a_s = &[F::from(1u64), F::from(2u64)];
        let c_s = &[3, 3];

        let seq = generate_sequence(r, a_s, c_s);

        let proof = GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::prove(seq, r, a_s, c_s);
    }

    /// Test that geometric_sequence() works correctly
    #[test]
    fn test_generate_sequence_0() {
        let r = F::from(2u64);
        let a_s = &[F::from(1u64), F::from(2u64)];
        let c_s = &[3, 3];
        let seq = generate_sequence(r, a_s, c_s);

        let expected = [1, 2, 4, 2, 4, 8]
            .iter()
            .map(|x| F::from(*x as u64))
            .collect::<Vec<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }

        assert!(GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::naive_verify(
            seq, r, a_s, c_s
        ));
    }

    #[test]
    fn test_generate_sequence_1() {
        let r = F::from(1u64);
        let a_s = &[F::from(1u64), F::from(1u64)];
        let c_s = &[1, 1];

        let seq = generate_sequence(r, a_s, c_s);
        let expected = [1, 1]
            .iter()
            .map(|x| F::from(*x as u64))
            .collect::<Vec<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }
        assert!(GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::naive_verify(
            seq, r, a_s, c_s
        ));
    }

    #[test]
    fn test_zero_over_k_for_geo_seq() {
        let r = Fr::from(2u64);
        let mut a_s = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut c_s = vec![5, 3, 10, 30];
        let m = c_s.iter().sum();

        let domain = GeneralEvaluationDomain::<Fr>::new(m).unwrap();
        let to_pad = domain.size() - m;
        if to_pad > 0 {
            a_s.push(Fr::from(0u64));
            c_s.push(to_pad);
        }

        let seq = generate_sequence::<Fr>(r, &a_s.as_slice(), &c_s.as_slice());
        let f = DensePolynomial::<Fr>::from_coefficients_slice(&domain.ifft(&seq));

        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);

        let concrete_oracles = [label_polynomial!(f)];
        let alphas = [Fr::from(1u64), domain.element(1)];

        let maximum_degree: usize = 80;

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, Some(&[2, 5])).unwrap();

        let (concrete_oracles_commitments, concrete_oracle_rands) =
            PC::commit(&ck, &concrete_oracles, None).unwrap();

        let proof = ZeroOverK::<F, KZG10<Bn254>, Blake2s>::prove(
            &concrete_oracles,
            &concrete_oracles_commitments,
            &concrete_oracle_rands,
            &geo_seq_vo,
            &alphas.to_vec(),
            &domain,
            &ck,
            &mut OsRng,
        )
        .unwrap();

        assert_eq!(
            true,
            ZeroOverK::<F, KZG10<Bn254>, Blake2s>::verify(
                proof,
                &concrete_oracles_commitments,
                &geo_seq_vo,
                &domain,
                &alphas,
                &vk,
            )
            .is_ok()
        );
    }
}
