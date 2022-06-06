#[cfg(test)]
mod tests {
    use crate::{
        geo_seq::GeoSeqTest
    };
    use ark_bn254::{Bn254, Fr};

    type F = Fr;

    /// Returns a concatenation of geometric sequences.
    /// @param r The multiplicative factor
    /// @param c_s The size of each sequence
    /// @param a_s The initial value of each sequence
    fn generate_sequence(r: F, c_s: &[usize], a_s: &[F]) -> Vec::<F> {
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
        let c_s = &[3, 3];
        let a_s = &[F::from(1u64), F::from(2u64)];

        let seq = generate_sequence(r, c_s, a_s);
    }

    /// Test that geometric_sequence() works correctly
    #[test]
    fn test_generate_sequence() {
        let mut seq = generate_sequence(
            F::from(2u64),
            &[3, 3],
            &[F::from(1u64), F::from(2u64)],
        );

        let mut expected = [1, 2, 4, 2, 4, 8].iter().map(|x| F::from(*x as u64)).collect::<Vec::<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }

        seq = generate_sequence(
            F::from(1u64),
            &[1, 1],
            &[F::from(1u64), F::from(1u64)],
        );
        expected = [1, 1].iter().map(|x| F::from(*x as u64)).collect::<Vec::<F>>();
        assert_eq!(expected.len(), seq.len());
        for (i, s) in seq.iter().enumerate() {
            assert_eq!(&expected[i], s);
        }
    }
}
