use crate::{
    commitment::HomomorphicPolynomialCommitment,
    geo_seq::proof::Proof as GeoSeqProof, non_zero_over_k::proof::Proof as NonZeroProof,
    zero_over_k::proof::Proof as ZeroProof,
};
use ark_ff::{PrimeField};

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    pub h1_commit: PC::Commitment,
    pub h2_commit: PC::Commitment,
    pub h1_seq_proof: GeoSeqProof<F, PC>,
    pub h2_seq_proof: GeoSeqProof<F, PC>,
    pub h_eq_row_m: ZeroProof<F, PC>,
    pub row_m_eq_col_m: ZeroProof<F, PC>,
    pub val_m_times_h2_proof: ZeroProof<F, PC>,
    pub val_plus_h2_proof: NonZeroProof<F, PC>,
}
