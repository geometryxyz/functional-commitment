use crate::error::{to_pc_error, Error};
use crate::geo_seq::proof::Proof;
use ark_ff::{to_bytes, PrimeField};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, QuerySet};
use ark_std::marker::PhantomData;
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use rand_core::OsRng;
use std::iter;
use zero_over_k::{
    virtual_oracle::generic_shifting_vo::{vo_term::VOTerm, GenericShiftingVO},
    zero_over_k::ZeroOverK,
    {geometric_seq_check, vo_constant},
};

pub mod proof;
mod tests;

pub struct GeoSeqTest<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> GeoSeqTest<F, PC, FS> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Geometric Sequence Test";
    // TODO: for both prove() and verify:
    // TODO: have an assertion that domain is large enough given m
    // TODO: move the padding outside and the check that the length is correct
    // TODO: verifier should check that the size of the sequence is correct given the domain
    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        common_ratio: F,
        f: &LabeledPolynomial<F, DensePolynomial<F>>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        f_rand: &PC::Randomness,
        sequence_initial_values: &Vec<F>,
        sequence_lengths: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        // Generate the GeoSequenceVO virtual oracle
        let alphas = [F::one(), domain.element(1)];
        let geo_seq_vo = GenericShiftingVO::new(
            &[0, 0],
            &alphas,
            geometric_seq_check!(common_ratio, sequence_lengths, domain),
        )?;

        let fs_bytes = &to_bytes![
            &Self::PROTOCOL_NAME,
            sequence_initial_values,
            sequence_lengths
                .iter()
                .map(|&x| x as u64)
                .collect::<Vec<_>>(),
            common_ratio,
            &[f_commit.clone()].to_vec(),
            &alphas.to_vec()
        ]
        .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FS::initialize(fs_bytes);

        let mut query_set = QuerySet::new();
        let sequence_starting_indices = iter::once(0)
            .chain(sequence_lengths.iter().scan(0, |st, elem| {
                *st += elem;
                Some(*st)
            }))
            .collect::<Vec<_>>();
        for (i, &pi) in sequence_starting_indices.iter().enumerate() {
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
            Some(rng),
        )
        .map_err(to_pc_error::<F, PC>)?;

        let z_proof = ZeroOverK::<F, PC, FS>::prove(
            &[f.clone()],
            &[f_commit.clone()],
            &[f_rand.clone()],
            f.degree_bound(),
            &geo_seq_vo,
            &domain,
            &ck,
            rng,
        )?;

        let proof = Proof::<F, PC> {
            z_proof,
            opening_proof,
        };
        Ok(proof)
    }

    pub fn verify(
        common_ratio: F,
        sequence_initial_values: &Vec<F>,
        sequence_lengths: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        enforced_degree_bound: Option<usize>,
        proof: Proof<F, PC>,
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        let bounded_f_commit = LabeledCommitment::new(
            f_commit.label().clone(),
            f_commit.commitment().clone(),
            enforced_degree_bound,
        );

        let alphas = [F::one(), domain.element(1)];
        let geo_seq_vo = GenericShiftingVO::new(
            &[0, 0],
            &alphas,
            geometric_seq_check!(common_ratio, sequence_lengths, domain),
        )?;

        // Test that for all i in n, check that f(gamma^p_i) = a_i
        let sequence_starting_indices = iter::once(0)
            .chain(sequence_lengths.iter().scan(0, |st, elem| {
                *st += elem;
                Some(*st)
            }))
            .collect::<Vec<_>>();
        let points = sequence_starting_indices
            .iter()
            .map(|&pi| domain.element(1).pow([pi as u64]))
            .collect::<Vec<_>>();

        let fs_bytes = &to_bytes![
            &Self::PROTOCOL_NAME,
            sequence_initial_values,
            sequence_lengths
                .iter()
                .map(|&x| x as u64)
                .collect::<Vec<_>>(),
            common_ratio,
            &[f_commit.clone()].to_vec(),
            &alphas.to_vec()
        ]
        .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FS::initialize(fs_bytes);

        let mut query_set = QuerySet::new();
        for (i, &point_i) in points.iter().enumerate() {
            query_set.insert((
                bounded_f_commit.label().clone(),
                (format!("gamma_pi_{}", i), point_i),
            ));
        }
        let mut evaluations = ark_poly_commit::Evaluations::new();
        for (&point_i, &a_i) in points.iter().zip(sequence_initial_values.iter()) {
            evaluations.insert((bounded_f_commit.label().clone(), point_i), a_i);
        }

        let separation_challenge = F::rand(&mut fs_rng);
        match PC::batch_check(
            vk,
            &[bounded_f_commit.clone()],
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
        ZeroOverK::<F, PC, FS>::verify(
            proof.z_proof,
            &[bounded_f_commit.clone()],
            enforced_degree_bound,
            &geo_seq_vo,
            &domain,
            vk,
        )?;

        Ok(())
    }

    #[allow(dead_code)]
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
