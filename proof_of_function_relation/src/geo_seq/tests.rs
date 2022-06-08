#[cfg(test)]
mod tests {
    use crate::{commitment::KZG10, geo_seq::GeoSeqTest};
    use ark_bn254::{Bn254, Fr};
    use blake2::Blake2s;

    type F = Fr;

    /// Returns a concatenation of geometric sequences.
    /// @param r The multiplicative factor
    /// @param c_s The size of each sequence
    /// @param a_s The initial value of each sequence
    fn generate_sequence(r: F, a_s: &[F], c_s: &[usize]) -> Vec<F> {
        assert_eq!(c_s.len(), a_s.len());
        let mut concatenation = Vec::<F>::default();

        for (i, a) in a_s.iter().enumerate() {
            for j in 0..c_s[i] {
                // r ** j (is there a pow() library function in Fp256?)
                let mut r_pow = F::from(1 as u64);
                for _ in 0..j {
                    r_pow = r_pow * r;
                }
                let val = *a * r_pow;
                concatenation.push(val);
            }
        }

        concatenation
    }

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
}
