use crate::zero_over_k::ZeroOverK;
use crate::geo_seq::proof::Proof;
use crate::error::{Error};
use crate::commitment::HomomorphicPolynomialCommitment;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use ark_ff::PrimeField;
mod tests;
mod proof;

pub struct GeoSeqTest <F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _zero_over_k: ZeroOverK<F, PC, D>
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> GeoSeqTest<F, PC, D> {
    pub fn prove(
        seq: &Vec::<F>,
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> Result<Proof<F>, Error>{
        // Check that seq is valid; otherwise, return an error
        if !GeoSeqTest::<F, PC, D>::naive_verify(&seq, r, a_s, c_s) {
            return Err(Error::InvalidGeoSeq);
        }

        // TODO: generate f() such that f(w^n) = a_i*r^n 
        //
        // Generate the GeoSequenceVO virtual oracle using instantiate()

        // dummy value for now
        let proof = Proof::<F> {
            blah: F::from(1u64)
        };
        Ok(proof)
    }

    pub fn verify(
        proof: Proof<F>,
        seq: &Vec::<F>,
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> bool {
        // TODO
        false
    }

    /// Inefficiently verify that the sequence is valid
    pub fn naive_verify(
        seq: &Vec::<F>,
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> bool {
        if a_s.len() != c_s.len() {
            return false;
        }

        if c_s.iter().fold(0, |x, y| x + y) != seq.len() {
            return false;
        }

        let mut i = 0;
        for (a_i, a) in a_s.iter().enumerate() {
            let mut r_pow = F::from(1 as u64);
            for _ in 0..c_s[a_i] {
                let expected = *a * r_pow;

                if expected != seq[i] {
                    return false;
                }

                r_pow = r * r_pow;
                i += 1;
            }
        }

        return true;
    }
}
