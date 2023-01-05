use crate::{
    error::Error, t_diag::TDiag, t_functional_triple::proof::Proof,
    t_strictly_lower_triangular_test::TStrictlyLowerTriangular,
};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use rand::Rng;
use std::marker::PhantomData;

use ark_ff::{PrimeField, SquareRootField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use std::io::BufReader;

pub mod proof;
// mod tests;

pub struct TFT<F: PrimeField + SquareRootField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng>
{
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<FS>,
}

impl<F, PC, FS> TFT<F, PC, FS>
where
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"t-FT Test";

    // TODO: change to use ark-marlin Index. (wait for a new release?)
    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        enforced_degree_bound: Option<usize>,
        // a
        row_a_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_a_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_a_commit: &LabeledCommitment<PC::Commitment>,
        col_a_commit: &LabeledCommitment<PC::Commitment>,
        row_a_random: &PC::Randomness,
        col_a_random: &PC::Randomness,
        // b
        row_b_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_b_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_b_commit: &LabeledCommitment<PC::Commitment>,
        col_b_commit: &LabeledCommitment<PC::Commitment>,
        row_b_random: &PC::Randomness,
        col_b_random: &PC::Randomness,
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
        fs_rng: &mut FS,
        rng: &mut R,
    ) -> Result<Vec<u8>, Error> {
        // 1. t-SLT test on A
        let a_slt_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            ck,
            t,
            domain_k,
            domain_h,
            row_a_poly,
            row_a_commit,
            row_a_random,
            col_a_poly,
            col_a_commit,
            col_a_random,
            enforced_degree_bound,
            fs_rng,
            rng,
        )?;

        // 2. t-SLT test on B
        let b_slt_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            ck,
            t,
            domain_k,
            domain_h,
            row_b_poly,
            row_b_commit,
            row_b_random,
            col_b_poly,
            col_b_commit,
            col_b_random,
            enforced_degree_bound,
            fs_rng,
            rng,
        )?;

        // 3. t-Diag test on C
        let c_diag_proof = TDiag::<F, PC, FS>::prove(
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
            enforced_degree_bound,
            domain_k,
            domain_h,
            domain_h.size(),
            rng,
        )?;

        let proof = Proof {
            a_slt_proof,
            b_slt_proof,
            c_diag_proof,
        };

        let mut writer = Vec::<u8>::new();
        let _ = proof
            .serialize(&mut writer)
            .map_err(|_| Error::ProofSerializationError)
            .unwrap();

        Ok(Vec::from(writer.as_slice()))
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
        enforced_degree_bound: Option<usize>,
        domain_h: &GeneralEvaluationDomain<F>,
        domain_k: &GeneralEvaluationDomain<F>,
        proof_bytes: Vec<u8>,
        fs_rng: &mut FS,
    ) -> Result<(), Error> {
        let reader = BufReader::new(proof_bytes.as_slice());
        let proof: Proof<F, PC> = Proof::<F, PC>::deserialize(reader).unwrap();

        TStrictlyLowerTriangular::<F, PC, FS>::verify(
            vk,
            ck,
            t,
            domain_k,
            domain_h,
            row_a_commitment,
            col_a_commitment,
            enforced_degree_bound,
            proof.a_slt_proof,
            fs_rng,
        )?;

        TStrictlyLowerTriangular::<F, PC, FS>::verify(
            vk,
            ck,
            t,
            domain_k,
            domain_h,
            row_b_commitment,
            col_b_commitment,
            enforced_degree_bound,
            proof.b_slt_proof,
            fs_rng,
        )?;

        TDiag::<F, PC, FS>::verify(
            vk,
            t,
            row_c_commitment,
            col_c_commitment,
            val_c_commitment,
            enforced_degree_bound,
            domain_h,
            domain_k,
            domain_h.size(),
            proof.c_diag_proof,
        )?;

        Ok(())
    }
}
