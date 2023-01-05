use crate::{
    discrete_log_comparison::DLComparison,
    error::{to_pc_error, Error},
    geo_seq::GeoSeqTest,
    subset_over_k::SubsetOverK,
    t_strictly_lower_triangular_test::proof::Proof,
    util::generate_sequence,
};
use ark_ff::{to_bytes, PrimeField, SquareRootField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use std::marker::PhantomData;

pub mod proof;
mod tests;

pub struct TStrictlyLowerTriangular<
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F, PC, FS> TStrictlyLowerTriangular<F, PC, FS>
where
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"t-Strictly Lower Triangular Test";

    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        row_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_commit: &LabeledCommitment<PC::Commitment>,
        row_random: &PC::Randomness,
        col_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_commit: &LabeledCommitment<PC::Commitment>,
        col_random: &PC::Randomness,
        enforced_degree_bound: Option<usize>,
        fs_rng: &mut FS,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        let fs_bytes = &to_bytes![&Self::PROTOCOL_NAME].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let r = domain_h.element(1);

        if t > domain_h.size() {
            return Err(Error::T2Large);
        }

        // Step 1: interpolate h
        let mut a_s = vec![domain_h.element(t)];
        let mut c_s = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }

        let seq = generate_sequence::<F>(r, &a_s.as_slice(), &c_s.as_slice());
        let h = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let h = LabeledPolynomial::new(String::from("h"), h, enforced_degree_bound, Some(1));

        let (commitment, rands) =
            PC::commit(&ck, &[h.clone()], Some(rng)).map_err(to_pc_error::<F, PC>)?;

        let h_commit = commitment[0].clone();

        // Step 2: Geometric sequence test on h
        let geo_seq_proof = GeoSeqTest::<F, PC, FS>::prove(
            ck, r, &h, &h_commit, &rands[0], &a_s, &c_s, domain_k, rng,
        )?;

        // Step 3: Subset over K between row_M and h
        let subset_proof = SubsetOverK::<F, PC, FS>::prove();

        // Step 4: Discrete Log Comparison between row_M and col_M
        let dl_proof = DLComparison::<F, PC, FS>::prove(
            ck,
            domain_k,
            domain_h,
            row_poly,
            row_commit,
            row_random,
            col_poly,
            col_commit,
            col_random,
            enforced_degree_bound,
            fs_rng,
            rng,
        )?;

        let proof = Proof {
            h_commit: h_commit.commitment().clone(),
            dl_proof,
            geo_seq_proof,
            subset_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        row_commit: &LabeledCommitment<PC::Commitment>,
        col_commit: &LabeledCommitment<PC::Commitment>,
        enforced_degree_bound: Option<usize>,
        proof: Proof<F, PC>,
        fs_rng: &mut FS,
    ) -> Result<(), Error> {
        // re-label the oracle commitments with the enforced degree bound
        let row_commit = LabeledCommitment::new(
            row_commit.label().clone(),
            row_commit.commitment().clone(),
            enforced_degree_bound,
        );
        let col_commit = LabeledCommitment::new(
            col_commit.label().clone(),
            col_commit.commitment().clone(),
            enforced_degree_bound,
        );

        // Step 2: Geometric sequence test on h
        let mut a_s = vec![domain_h.element(t)];
        let mut c_s = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }

        let h_commit =
            LabeledCommitment::new(String::from("h"), proof.h_commit, enforced_degree_bound);

        GeoSeqTest::<F, PC, FS>::verify(
            domain_h.element(1),
            &a_s,
            &c_s,
            domain_k,
            &h_commit,
            enforced_degree_bound,
            proof.geo_seq_proof,
            vk,
        )?;

        // Step 3: Subset over K between row_M and h
        SubsetOverK::<F, PC, FS>::verify(proof.subset_proof)?;

        // Step 4: Discrete Log Comparison between row_M and col_M
        DLComparison::<F, PC, FS>::verify(
            vk,
            ck,
            domain_k,
            domain_h,
            &row_commit,
            &col_commit,
            enforced_degree_bound,
            proof.dl_proof,
            fs_rng,
        )?;

        Ok(())
    }
}
