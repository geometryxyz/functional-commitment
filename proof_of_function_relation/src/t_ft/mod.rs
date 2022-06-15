use crate::{
    commitment::HomomorphicPolynomialCommitment,
    t_ft::proof::Proof,
    error::{to_pc_error, Error},
    t_strictly_lower_triangular_test::{ TStrictlyLowerTriangular },
    t_strictly_lower_triangular_test::proof::{ Proof as TSLTProof },
    t_diag::{ TDiag },
    t_diag::proof::{ Proof as TDiagProof },
};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use std::marker::PhantomData;
use rand::Rng;
use ark_marlin::rng::FiatShamirRng;

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
        row_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        val_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_commit: &LabeledCommitment<PC::Commitment>,
        col_commit: &LabeledCommitment<PC::Commitment>,
        val_commit: &LabeledCommitment<PC::Commitment>,
        _row_m_random: &PC::Randomness,
        _col_m_random: &PC::Randomness,
        _val_m_random: &PC::Randomness,
        fs_rng: &mut FiatShamirRng<D>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        // 1. t-SLT test on row, col
        // TODO: need to do the test on val but t-SLT for row, col, val isn't implemented?
        let t_slt_proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
            ck,
            t,
            domain_k,
            domain_h,
            row_poly,
            col_poly,
            row_commit,
            col_commit,
            fs_rng,
            rng,
        ).unwrap();

        // 2. t-Diag test on row, col, val
        let t_diag_proof = TDiag::<F, PC, D>::prove(
            ck,
            t,
            row_poly,
            col_poly,
            val_poly,
            row_commit,
            col_commit,
            val_commit,
            _row_m_random,
            _col_m_random,
            _val_m_random,
            domain_k,
            domain_h,
            rng,
        ).unwrap();

        let proof = Proof {
            t_slt_proof,
            t_diag_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        ck: &PC::CommitterKey,
        t: usize,
        row_m_commitment: &LabeledCommitment<PC::Commitment>,
        col_m_commitment: &LabeledCommitment<PC::Commitment>,
        val_m_commitment: &LabeledCommitment<PC::Commitment>,
        domain_h: &GeneralEvaluationDomain<F>,
        domain_k: &GeneralEvaluationDomain<F>,
        proof: Proof<F, PC>,
        fs_rng: &mut FiatShamirRng<D>,
    ) -> Result<(), Error> {
        let t_slt_is_valid = TStrictlyLowerTriangular::<F, PC, D>::verify(
            vk,
            ck,
            t,
            domain_k,
            domain_h,
            row_m_commitment,
            col_m_commitment,
            proof.t_slt_proof,
            fs_rng,
        );

        if t_slt_is_valid.is_err() {
            return Err(Error::ProofVerificationError);
        }

        let t_diag_is_valid = TDiag::<F, PC, D>::verify(
            vk,
            t,
            row_m_commitment,
            col_m_commitment,
            val_m_commitment,
            domain_h,
            domain_k,
            proof.t_diag_proof,
        );

        if t_diag_is_valid.is_err() {
            return Err(Error::ProofVerificationError);
        }

        Ok(())
    }
}
