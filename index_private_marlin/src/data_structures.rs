use crate::ahp::indexer::*;
use crate::ahp::prover::ProverMsg;
use crate::Vec;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{BatchLCProof, PolynomialCommitment};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::{
    format,
    io::{Read, Write},
};

use ::zero_over_k::zero_over_k::proof::Proof as ZeroOverKProof;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// The universal public parameters for the argument system.
pub type UniversalSRS<F, PC> = <PC as PolynomialCommitment<F, DensePolynomial<F>>>::UniversalParams;

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// Verification key for a specific index (i.e., R1CS matrices).
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IndexVerifierKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// Stores information about the size of the index, as well as its field of
    /// definition.
    pub index_info: IndexInfo<F>,
    /// Commitments to the indexed polynomials.
    pub index_comms: Vec<PC::Commitment>,
    /// The verifier key for this index, trimmed from the universal SRS.
    pub verifier_key: PC::VerifierKey,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> ark_ff::ToBytes for IndexVerifierKey<F, PC> {
    fn write<W: Write>(&self, mut w: W) -> ark_std::io::Result<()> {
        self.index_info.write(&mut w)?;
        self.index_comms.write(&mut w)
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for IndexVerifierKey<F, PC> {
    fn clone(&self) -> Self {
        Self {
            index_comms: self.index_comms.clone(),
            index_info: self.index_info.clone(),
            verifier_key: self.verifier_key.clone(),
        }
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> IndexVerifierKey<F, PC> {
    /// Iterate over the commitments to indexed polynomials in `self`.
    pub fn iter(&self) -> impl Iterator<Item = &PC::Commitment> {
        self.index_comms.iter()
    }
}

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// Proving key for a specific index (i.e., R1CS matrices).
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IndexProverKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// The index verifier key.
    pub index_vk: IndexVerifierKey<F, PC>,
    /// The randomness for the index polynomial commitments.
    pub index_comm_rands: Vec<PC::Randomness>,
    /// The index itself.
    pub index: Index<F>,
    /// The committer key for this index, trimmed from the universal SRS.
    pub committer_key: PC::CommitterKey,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for IndexProverKey<F, PC>
where
    PC::Commitment: Clone,
{
    fn clone(&self) -> Self {
        Self {
            index_vk: self.index_vk.clone(),
            index_comm_rands: self.index_comm_rands.clone(),
            index: self.index.clone(),
            committer_key: self.committer_key.clone(),
        }
    }
}

/// Proving key for a index private specific index (i.e., R1CS matrices).
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IndexPrivateProverKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// The index verifier key.
    pub index_private_vk: IndexPrivateVerifierKey<F, PC>,
    /// The index itself.
    pub index: IndexPrivateIndex<F>,
    /// TODO: this should be changed to just PrivateIndex
    /// The committer key for this index, trimmed from the universal SRS.
    pub committer_key: PC::CommitterKey,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for IndexPrivateProverKey<F, PC>
where
    PC::Commitment: Clone,
{
    fn clone(&self) -> Self {
        Self {
            index_private_vk: self.index_private_vk.clone(),
            index: self.index.clone(),
            committer_key: self.committer_key.clone(),
        }
    }
}

/// Verifier key that will be used in index private version
/// Prover will commit to matrix arithmetizations and this data will be used for
/// slt and diag testing
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IndexPrivateVerifierKey<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
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
    pub polys: Vec<PC::Commitment>,

    /// verifier key
    pub verifier_key: PC::VerifierKey,

    /// Stores information about the size of the index, as well as its field of
    /// definition.
    pub index_info: IndexInfo<F>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Clone for IndexPrivateVerifierKey<F, PC> {
    fn clone(&self) -> Self {
        Self {
            polys: self.polys.clone(),
            index_info: self.index_info.clone(),
            verifier_key: self.verifier_key.clone(),
        }
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> ark_ff::ToBytes
    for IndexPrivateVerifierKey<F, PC>
{
    fn write<W: Write>(&self, mut w: W) -> ark_std::io::Result<()> {
        self.index_info.write(&mut w)?;
        self.polys.write(&mut w)
    }
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> IndexPrivateVerifierKey<F, PC> {
    /// Iterate over the commitments to indexed polynomials in `self`.
    pub fn iter(&self) -> impl Iterator<Item = &PC::Commitment> {
        self.polys.iter()
    }
}

/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */

/// A zkSNARK proof.
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// Commitments to the polynomials produced by the AHP prover.
    pub commitments: Vec<Vec<PC::Commitment>>,
    /// Evaluations of these polynomials.
    pub evaluations: Vec<F>,
    /// The field elements sent by the prover.
    pub prover_messages: Vec<ProverMsg<F>>,
    /// An evaluation proof from the polynomial commitment.
    pub pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> Proof<F, PC> {
    /// Construct a new proof.
    pub fn new(
        commitments: Vec<Vec<PC::Commitment>>,
        evaluations: Vec<F>,
        prover_messages: Vec<ProverMsg<F>>,
        pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,
    ) -> Self {
        Self {
            commitments,
            evaluations,
            prover_messages,
            pc_proof,
        }
    }

    /// Prints information about the size of the proof.
    pub fn print_size_info(&self) {
        use ark_poly_commit::PCCommitment;

        let mut num_comms_without_degree_bounds = 0;
        let mut num_comms_with_degree_bounds = 0;
        let mut size_bytes_comms_without_degree_bounds = 0;
        let mut size_bytes_comms_with_degree_bounds = 0;
        for c in self.commitments.iter().flat_map(|c| c) {
            if !c.has_degree_bound() {
                num_comms_without_degree_bounds += 1;
                size_bytes_comms_without_degree_bounds += c.serialized_size();
            } else {
                num_comms_with_degree_bounds += 1;
                size_bytes_comms_with_degree_bounds += c.serialized_size();
            }
        }

        let proofs: Vec<PC::Proof> = self.pc_proof.proof.clone().into();
        let num_proofs = proofs.len();
        let size_bytes_proofs = self.pc_proof.proof.serialized_size();

        let num_evals = self.evaluations.len();
        let evals_size_in_bytes = self.evaluations.serialized_size();
        let num_prover_messages: usize = self
            .prover_messages
            .iter()
            .map(|v| match v {
                ProverMsg::EmptyMessage => 0,
                ProverMsg::FieldElements(elems) => elems.len(),
            })
            .sum();
        let prover_msg_size_in_bytes = self.prover_messages.serialized_size();
        let arg_size = self.serialized_size();
        let stats = format!(
            "Argument size in bytes: {}\n\n\
             Number of commitments without degree bounds: {}\n\
             Size (in bytes) of commitments without degree bounds: {}\n\
             Number of commitments with degree bounds: {}\n\
             Size (in bytes) of commitments with degree bounds: {}\n\n\
             Number of evaluation proofs: {}\n\
             Size (in bytes) of evaluation proofs: {}\n\n\
             Number of evaluations: {}\n\
             Size (in bytes) of evaluations: {}\n\n\
             Number of field elements in prover messages: {}\n\
             Size (in bytes) of prover message: {}\n",
            arg_size,
            num_comms_without_degree_bounds,
            size_bytes_comms_without_degree_bounds,
            num_comms_with_degree_bounds,
            size_bytes_comms_with_degree_bounds,
            num_proofs,
            size_bytes_proofs,
            num_evals,
            evals_size_in_bytes,
            num_prover_messages,
            prover_msg_size_in_bytes,
        );
        add_to_trace!(|| "Statistics about proof", || stats);
    }
}

/// A zkSNARK index private proof.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct IndexPrivateProof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    /// Commitments to the polynomials produced by the AHP prover.
    pub commitments: Vec<Vec<PC::Commitment>>,
    /// Evaluations of these polynomials.
    pub evaluations: Vec<F>,
    /// The field elements sent by the prover.
    pub prover_messages: Vec<ProverMsg<F>>,
    /// An evaluation proof from the polynomial commitment.
    pub pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,

    pub rational_sumcheck_zero_over_k_proof: ZeroOverKProof<F, PC>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> IndexPrivateProof<F, PC> {
    /// Construct a new proof.
    pub fn new(
        commitments: Vec<Vec<PC::Commitment>>,
        evaluations: Vec<F>,
        prover_messages: Vec<ProverMsg<F>>,
        pc_proof: BatchLCProof<F, DensePolynomial<F>, PC>,
        rational_sumcheck_zero_over_k_proof: ZeroOverKProof<F, PC>,
    ) -> Self {
        Self {
            commitments,
            evaluations,
            prover_messages,
            pc_proof,
            rational_sumcheck_zero_over_k_proof,
        }
    }
}
