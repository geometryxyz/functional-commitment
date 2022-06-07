#[cfg(test)]
mod tests {
    use crate::{
        geo_seq::GeoSeqTest,
        commitment::{KZG10},
    };
    use ark_bn254::{Bn254, Fr};
    use blake2::Blake2s;

    type F = Fr;

    /// Returns a concatenation of geometric sequences.
    /// @param r The multiplicative factor
    /// @param c_s The size of each sequence
    /// @param a_s The initial value of each sequence
    fn generate_sequence(
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> Vec::<F> {
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
        assert!(is_valid);
    }

    /// Test that geometric_sequence() works correctly
    #[test]
    fn test_generate_sequence_0() {
        let r = F::from(2u64);
        let a_s = &[F::from(1u64), F::from(2u64)];
        let c_s = &[3, 3];
        let seq = generate_sequence(r, a_s, c_s);

        let expected = [1, 2, 4, 2, 4, 8].iter().map(|x| F::from(*x as u64)).collect::<Vec::<F>>();
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
        let expected = [1, 1].iter().map(|x| F::from(*x as u64)).collect::<Vec::<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }
        assert!(GeoSeqTest::<F, KZG10<Bn254>, Blake2s>::naive_verify(&seq, r, a_s, c_s));
    }
}
