use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::Error;
use crate::label_polynomial;
use crate::geo_seq::proof::Proof;
use crate::zero_over_k::ZeroOverK;
use crate::zero_over_k::proof::{Proof as Z_Proof};
use crate::virtual_oracle::{
    Description, EvaluationsProvider, GeoSequenceVO, NormalizedVirtualOracle, Term,
    VirtualOracle,
};
use ark_ff::PrimeField;
use ark_bn254::Fr;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};
use ark_poly_commit::{LabeledPolynomial, LabeledCommitment};
use rand_core::OsRng;

mod proof;
mod tests;

pub struct GeoSeqTest<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _zero_over_k: ZeroOverK<F, PC, D>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> GeoSeqTest<F, PC, D> {
    pub fn prove(
        seq: &Vec<F>,
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> Result<Proof<F, PC>, Error>{
        // Check that seq is valid; otherwise, return an error
        if !GeoSeqTest::<F, PC, D>::naive_verify(&seq, r, a_s, c_s) {
            return Err(Error::InvalidGeoSeq);
        }

        // Generate f() such that f(w^n) = a_i*r^n 
        let m: usize = c_s.iter().sum();
        let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let to_pad: usize = domain.size() - m;

        // Generate the GeoSequenceVO virtual oracle
        let coeffs = domain.ifft(seq);
        let f = DensePolynomial::<F>::from_coefficients_vec(coeffs);
        let mut a_vec: Vec<F> = a_s.to_vec();
        let mut c_vec = c_s.to_vec();

        if to_pad > 0 {
            a_vec.push(F::from(0u64));
            c_vec.push(to_pad);
        }

        let alphas = [F::from(1u64)].to_vec();
        let concrete_oracles = [label_polynomial!(f)];

        let geo_seq_vo = GeoSequenceVO::new(&c_vec, domain.element(1), r);

        let instantiatied_geo_seq_vo = geo_seq_vo
            .instantiate(&concrete_oracles, &[domain.element(1)])
            .unwrap();

        // perform zero over k
        let maximum_degree: usize = instantiatied_geo_seq_vo.coeffs.len();
        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();

        // TODO: why [2, 5]? this is the degree bounds... needs to be implented correctly in ZOK
        // upstream.
        let (ck, vk) = PC::trim(&pp, maximum_degree, 0, Some(&[2, 5])).unwrap();

        let (concrete_oracle_commitments, concrete_oracle_rands) =
            PC::commit(&ck, &concrete_oracles, None).unwrap();

        assert_eq!(concrete_oracle_commitments.len(), 1);

        let z_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &concrete_oracle_commitments,
            &concrete_oracle_rands,
            &geo_seq_vo,
            &alphas,
            domain,
            &ck,
            &mut OsRng,
        );

        let proof = Proof::<F, PC> {
            z_proof: z_proof.unwrap(),
            concrete_oracle_commitment: concrete_oracle_commitments[0].clone(),
            seq: seq.clone(),
            r: r,
            a_vec: a_vec,
            c_vec: c_vec,
            vk: vk,
        };
        Ok(proof)
    }

    pub fn verify(
        proof: Proof<F, PC>,
        seq: &Vec<F>,
        r: F,
        a_s: &[F],
        c_s: &[usize],
    ) -> Result<(), Error> {
        // Generate f() such that f(w^n) = a_i*r^n 
        let m: usize = c_s.iter().sum();
        let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let to_pad: usize = domain.size() - m;

        // Generate the GeoSequenceVO virtual oracle
        let coeffs = domain.ifft(seq);
        let f = DensePolynomial::<F>::from_coefficients_vec(coeffs);
        let mut a_vec: Vec<F> = a_s.to_vec();
        let mut c_vec = c_s.to_vec();

        if to_pad > 0 {
            a_vec.push(F::from(0u64));
            c_vec.push(to_pad);
        }

        let alphas = [F::from(1u64)].to_vec();
        let concrete_oracles = [label_polynomial!(f)];

        let geo_seq_vo = GeoSequenceVO::new(&c_vec, domain.element(1), r);

        let is_valid = ZeroOverK::<F, PC, D>::verify(
            proof.z_proof,
            &[proof.concrete_oracle_commitment].to_vec(),
            &geo_seq_vo,
            domain,
            &alphas,
            &proof.vk,
        );

        match is_valid {
            Ok(r) => r,
            Err(e) => panic!("{:?}", e),
        };

        Ok(())
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
