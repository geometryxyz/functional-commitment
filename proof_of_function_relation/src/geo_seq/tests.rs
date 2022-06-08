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
        
        // TODO: pull from dev!
        // TODO: implement the following. need to refactor things
        // Let's describe how to do a proper back and forth between the prover and verifier

        // 1. We have a sequence defined by r, a_s, and c_s
        // 2. Both the verifier and prover receive r, a_s, and c_s
        // 3. Verifier has oracle access to the function f (i.e. verifier already hold a commitment)
        //   - Test generates seq, interpolate the polynomial out of it to get f
        // 4. Prover generates the VO. Next, the prover runs zero over k for
        //    this VO. The prover sends the following to the verifier:
        //    - zero over k proof
        // 5. Verifier generates VO (note that it already knows f). Verifier
        //    verifies the zero over k proof. It also checks that:
        //    - for all i in n, check that f(gamma^p_i) = a_i

        // TODO: implement this after the above! It's necessary to make the proof succinct.
        // Verifier emits a query set (for the f_gamma check)
        // Prover will evaluate f at those points and return opening proof
        // Prover runs zero over k prove
        // Verifier runs zero over k verify
        // TODO: qn: can this be made noninteractive? I want the prover to be bound to the query
        // set that the verifier will emit so that the verifier doesn't need to send anything
        // before the prover. It should be as simple as "prover sends proof, verifier verifies
        // proof" instead of "verifier emits query set, prover sends proof, verifier verifies
        // proof"

        //
        // fn prove(r, a_s, c_s) {
        //     generate virtual oracle
        //     run zero over k
        //     output proof
        // }
        //
        // fn verify(proof, r, a_s, c_s) {
        //     generate seq from r, a_s, c_s
        //     derive f
        //     for all i in n, check that f(gamma^p_i) = a_i
        //     (todo: emit a query set)
        //     generate virtual oracle
        //     verify the zero over k proof
        // }


        let proof = GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::prove(
            &seq,
            r,
            a_s,
            c_s,
        );

        let is_valid = GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::verify(
            proof.unwrap(),
            &seq,
            r,
            a_s,
            c_s,
        );
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

        assert!(GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::naive_verify(&seq, r, a_s, c_s));
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
        assert!(GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::naive_verify(&seq, r, a_s, c_s));
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
