use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use crate::geo_seq::proof::Proof;
use crate::virtual_oracle::geometric_sequence_vo::GeoSequenceVO;
use crate::zero_over_k::ZeroOverK;
use ark_ff::{to_bytes, PrimeField};
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, QuerySet};
use ark_std::marker::PhantomData;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use rand_core::OsRng;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use std::io::{BufReader, BufWriter};

pub mod proof;
mod tests;

pub struct GeoSeqTest<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> GeoSeqTest<F, PC, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Geometric Sequence Test";
    // TODO: for both prove() and verify:
    // TODO: degree bounds!
    // TODO: have an assertion that domain is large enough given m
    // TODO: move the padding outside and the check that the length is correct
    // TODO: verifier should check that the size of the sequence is correct given the domain
    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        r: F,
        f: &LabeledPolynomial<F, DensePolynomial<F>>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        f_rand: &PC::Randomness,
        a_s: &Vec<F>,
        c_s: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        rng: &mut R,
    //) -> Result<Proof<F, PC>, Error> {
    ) -> Result<Vec<u8>, Error> {
        // Generate the GeoSequenceVO virtual oracle
        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);

        let alphas = [F::from(1u64), domain.element(1)];

        let fs_bytes = &to_bytes![
            &Self::PROTOCOL_NAME,
            a_s,
            c_s.iter().map(|&x| x as u64).collect::<Vec<_>>(),
            r,
            &[f_commit.clone()].to_vec(),
            &alphas.to_vec()
        ]
        .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FiatShamirRng::<D>::from_seed(fs_bytes);

        let mut query_set = QuerySet::new();
        let pi_s = geo_seq_vo.get_pi_s();
        for (i, &pi) in pi_s.iter().enumerate() {
            query_set.insert((
                f.label().clone(),
                (
                    format!("gamma_pi_{}", i),
                    domain.element(1).pow([pi as u64]),
                ),
            ));
        }

        let separation_challenge = F::rand(&mut fs_rng);
        let opening_proof = PC::batch_open(
            ck,
            &[f.clone()],
            &[f_commit.clone()],
            &query_set,
            separation_challenge,
            &[f_rand.clone()],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let z_proof = ZeroOverK::<F, PC, D>::prove(
            &[f.clone()],
            &[f_commit.clone()],
            &[f_rand.clone()],
            &geo_seq_vo,
            &alphas.to_vec(),
            &domain,
            &ck,
            rng,
        )?;

        let proof = Proof::<F, PC> {
            z_proof,
            opening_proof,
        };

        let proof_bytes = Vec::new();
        let writer = BufWriter::new(proof_bytes.clone());
        let _ = proof.serialize(writer).map_err(|_| Error::ProofSerializationError)?;

        Ok(proof_bytes)
        //Ok(proof)
    }

    pub fn verify(
        r: F,
        a_s: &Vec<F>,
        c_s: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        //proof: Proof<F, PC>,
        proof_bytes: Vec<u8>,
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        let reader = BufReader::new(proof_bytes.as_slice());
        let proof: Proof::<F, PC> = Proof::<F, PC>::deserialize(reader)
            .map_err(|_| Error::ProofDeserializationError)?;

        let alphas = [F::from(1u64), domain.element(1)];
        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);

        // Test that for all i in n, check that f(gamma^p_i) = a_i
        let pi_s = geo_seq_vo.get_pi_s();
        let points = pi_s
            .iter()
            .map(|&pi| domain.element(1).pow([pi as u64]))
            .collect::<Vec<_>>();

        let fs_bytes = &to_bytes![
            &Self::PROTOCOL_NAME,
            a_s,
            c_s.iter().map(|&x| x as u64).collect::<Vec<_>>(),
            r,
            &[&f_commit].to_vec(),
            &alphas.to_vec()
        ]
        .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FiatShamirRng::<D>::from_seed(fs_bytes);

        let mut query_set = QuerySet::new();
        for (i, &point_i) in points.iter().enumerate() {
            query_set.insert((
                f_commit.label().clone(),
                (format!("gamma_pi_{}", i), point_i),
            ));
        }
        let mut evaluations = ark_poly_commit::Evaluations::new();
        for (&point_i, &a_i) in points.iter().zip(a_s.iter()) {
            evaluations.insert((f_commit.label().clone(), point_i), a_i);
        }

        let separation_challenge = F::rand(&mut fs_rng);
        match PC::batch_check(
            vk,
            &[f_commit.clone()],
            &query_set,
            &evaluations,
            &proof.opening_proof,
            separation_challenge,
            &mut OsRng,
        ) {
            Ok(true) => Ok(()),
            Ok(false) => Err(Error::BatchCheckError),
            Err(e) => panic!("{:?}", e),
        }?;

        // TODO: is this check done or does the function return before it should if the above batch
        // check passes?

        // let seq = generate_sequence::<F>(r, &a_s.as_slice(), &c_s.as_slice());
        // let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&seq));

        // TODO: raise a different error?
        ZeroOverK::<F, PC, D>::verify(
            proof.z_proof,
            &[f_commit.clone()],
            &geo_seq_vo,
            &domain,
            &alphas,
            vk,
        )?;

        Ok(())
    }

    /// Inefficiently verify that the sequence is valid
    pub fn naive_verify(seq: &Vec<F>, r: F, a_s: &[F], c_s: &[usize]) -> bool {
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
