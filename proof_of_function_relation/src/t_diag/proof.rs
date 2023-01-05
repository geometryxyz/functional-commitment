use crate::{geo_seq::proof::Proof as GeoSeqProof, non_zero_over_k::proof::Proof as NonZeroProof};
use ark_ff::PrimeField;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use zero_over_k::zero_over_k::proof::Proof as ZeroProof;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    pub h1_commit: PC::Commitment,
    pub h2_commit: PC::Commitment,
    pub h1_seq_proof: GeoSeqProof<F, PC>,
    pub h2_seq_proof: GeoSeqProof<F, PC>,
    pub h_eq_row_m: ZeroProof<F, PC>,
    pub row_m_eq_col_m: ZeroProof<F, PC>,
    pub val_m_times_h2_proof: ZeroProof<F, PC>,
    pub val_plus_h2_proof: NonZeroProof<F, PC>,
}
