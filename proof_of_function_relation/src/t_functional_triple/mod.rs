use crate::{
    commitment::HomomorphicPolynomialCommitment,
    error::Error,
    t_diag::TDiag,
    t_functional_triple::proof::Proof,
    t_strictly_lower_triangular_test::TStrictlyLowerTriangular,
};
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{
    univariate::DensePolynomial, GeneralEvaluationDomain,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use rand::Rng;
use std::marker::PhantomData;

use ark_ff::{PrimeField, SquareRootField};
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`

pub mod proof;
mod tests;

struct TFT<F: PrimeField + SquareRootField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> TFT<F, PC, D>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"t-FT Test";

    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        // a
        row_a_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_a_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_a_commit: &LabeledCommitment<PC::Commitment>,
        col_a_commit: &LabeledCommitment<PC::Commitment>,
        _row_a_random: &PC::Randomness,
        _col_a_random: &PC::Randomness,
        // b
        row_b_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_b_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_b_commit: &LabeledCommitment<PC::Commitment>,
        col_b_commit: &LabeledCommitment<PC::Commitment>,
        _row_b_random: &PC::Randomness,
        _col_b_random: &PC::Randomness,
        // c
        row_c_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_c_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        val_c_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_c_commit: &LabeledCommitment<PC::Commitment>,
        col_c_commit: &LabeledCommitment<PC::Commitment>,
        val_c_commit: &LabeledCommitment<PC::Commitment>,
        row_c_random: &PC::Randomness,
        col_c_random: &PC::Randomness,
        val_c_random: &PC::Randomness,
        //rands
        fs_rng: &mut FiatShamirRng<D>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        // ) -> Result<Proof, Error> {
        // 1. t-SLT test on A
        let a_slt_proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            ck,
            t,
            domain_k,
            domain_h,
            row_a_poly,
            col_a_poly,
            row_a_commit,
            col_a_commit,
            fs_rng,
            rng,
        )?;

        // 2. t-SLT test on B
        let b_slt_proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            ck,
            t,
            domain_k,
            domain_h,
            row_b_poly,
            col_b_poly,
            row_b_commit,
            col_b_commit,
            fs_rng,
            rng,
        )?;

        // 3. t-Diag test on C
        let c_diag_proof = TDiag::<F, PC, D>::prove(
            ck,
            t,
            row_c_poly,
            col_c_poly,
            val_c_poly,
            row_c_commit,
            col_c_commit,
            val_c_commit,
            row_c_random,
            col_c_random,
            val_c_random,
            domain_k,
            domain_h,
            rng,
        )?;

        let proof = Proof {
            a_slt_proof,
            b_slt_proof,
            c_diag_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        ck: &PC::CommitterKey,
        t: usize,
        row_a_commitment: &LabeledCommitment<PC::Commitment>,
        col_a_commitment: &LabeledCommitment<PC::Commitment>,
        row_b_commitment: &LabeledCommitment<PC::Commitment>,
        col_b_commitment: &LabeledCommitment<PC::Commitment>,
        row_c_commitment: &LabeledCommitment<PC::Commitment>,
        col_c_commitment: &LabeledCommitment<PC::Commitment>,
        val_c_commitment: &LabeledCommitment<PC::Commitment>,
        domain_h: &GeneralEvaluationDomain<F>,
        domain_k: &GeneralEvaluationDomain<F>,
        proof: Proof<F, PC>,
        fs_rng: &mut FiatShamirRng<D>,
    ) -> Result<(), Error> {
        TStrictlyLowerTriangular::<F, PC, D>::verify(
            vk,
            ck,
            t,
            domain_k,
            domain_h,
            row_a_commitment,
            col_a_commitment,
            proof.a_slt_proof,
            fs_rng,
        )?;

        TStrictlyLowerTriangular::<F, PC, D>::verify(
            vk,
            ck,
            t,
            domain_k,
            domain_h,
            row_b_commitment,
            col_b_commitment,
            proof.b_slt_proof,
            fs_rng,
        )?;

        TDiag::<F, PC, D>::verify(
            vk,
            t,
            row_c_commitment,
            col_c_commitment,
            val_c_commitment,
            domain_h,
            domain_k,
            proof.c_diag_proof,
        )?;

        Ok(())
    }
}
