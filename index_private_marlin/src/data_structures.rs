use crate::ahp::indexer::*;
use crate::ahp::prover::ProverMsg;
use crate::Vec;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{BatchLCProof, PolynomialCommitment};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

use ::zero_over_k::zero_over_k::proof::Proof as ZeroOverKProof;
use ac_compiler::R1CSfIndex;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// The universal public parameters for the argument system.
pub type UniversalSRS<F, PC> = <PC as PolynomialCommitment<F, DensePolynomial<F>>>::UniversalParams;

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// Proving key for a specific index (i.e., R1CS matrices).
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ProverKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// The index verifier key.
    pub vk: VerifierKey<F, PC>,
    /// The randomness used for hiding matrix ldes
    pub rands: Vec<PC::Randomness>,
    /// The index itself.
    pub index: Index<F>,
    /// The committer key for this index, trimmed from the universal SRS.
    pub committer_key: PC::CommitterKey,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for ProverKey<F, PC>
where
    PC::Commitment: Clone,
{
    fn clone(&self) -> Self {
        Self {
            vk: self.vk.clone(),
            rands: self.rands.clone(),
            index: self.index.clone(),
            committer_key: self.committer_key.clone(),
        }
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> ProverKey<F, PC> {
    pub fn get_rands(&self) -> Vec<PC::Randomness> {
        self.rands.clone()
    }
}

/// Verifier key that will be used in index private version
/// Prover will commit to matrix arithmetizations and this data will be used for
/// slt and diag testing
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    // /// matrix a row commitment
    // pub a_row_commit: PC::Commitment,

    // /// matrix a col commitment
    // pub a_col_commit: PC::Commitment,

    // /// matrix a val commitment
    // pub a_val_commit: PC::Commitment,

    // /// matrix b row commitment
    // pub b_row_commit: PC::Commitment,

    // /// matrix b col commitment
    // pub b_col_commit: PC::Commitment,

    // /// matrix b val commitment
    // pub b_val_commit: PC::Commitment,

    // /// matrix c row commitment
    // pub c_row_commit: PC::Commitment,

    // /// matrix c col commitment
    // pub c_col_commit: PC::Commitment,

    // /// matrix c val commitment
    // pub c_val_commit: PC::Commitment,
    /// a(row, col, val), b(row, col, val), c(row, col, val),

    /// commitments of (row, col, val) for each matrix
    pub commits: Vec<PC::Commitment>,

    /// verifier key
    pub verifier_key: PC::VerifierKey,

    /// Stores information about the size of the index, as well as its field of
    /// definition.
    pub index_info: R1CSfIndex,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for VerifierKey<F, PC> {
    fn clone(&self) -> Self {
        Self {
            commits: self.commits.clone(),
            index_info: self.index_info.clone(),
            verifier_key: self.verifier_key.clone(),
        }
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> ark_ff::ToBytes for VerifierKey<F, PC> {
    fn write<W: Write>(&self, mut w: W) -> ark_std::io::Result<()> {
        self.index_info.write(&mut w)?;
        self.commits.write(&mut w)
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> VerifierKey<F, PC> {
    /// Iterate over the commitments to indexed polynomials in `self`.
    pub fn iter(&self) -> impl Iterator<Item = &PC::Commitment> {
        self.commits.iter()
    }
}

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// A zkSNARK index private proof.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// Commitments to the polynomials produced by the AHP prover.
    pub commitments: Vec<Vec<PC::Commitment>>,
    /// Evaluations of these polynomials.
    pub evaluations: Vec<F>,
    /// The field elements sent by the prover.
    pub prover_messages: Vec<ProverMsg<F>>,
    /// An evaluation proof from the polynomial commitment.
    pub pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,

    pub rational_sumcheck_zero_over_k_proof: ZeroOverKProof<F, PC>,
    pub well_formation_proof: ZeroOverKProof<F, PC>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Proof<F, PC> {
    /// Construct a new proof.
    pub fn new(
        commitments: Vec<Vec<PC::Commitment>>,
        evaluations: Vec<F>,
        prover_messages: Vec<ProverMsg<F>>,
        pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,
        rational_sumcheck_zero_over_k_proof: ZeroOverKProof<F, PC>,
        well_formation_proof: ZeroOverKProof<F, PC>,
    ) -> Self {
        Self {
            commitments,
            evaluations,
            prover_messages,
            pc_proof,
            rational_sumcheck_zero_over_k_proof,
            well_formation_proof,
        }
    }
}
