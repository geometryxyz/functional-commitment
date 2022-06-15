use crate::{
    commitment::HomomorphicPolynomialCommitment, discrete_log_comparison::proof::Proof as DLProof,
    geo_seq::proof::Proof as GeoSeqProof, subset_over_k::proof::Proof as SubsetProof,
};
use ark_ff::{PrimeField, SquareRootField};

pub struct Proof<F, PC>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
{
    pub h_commit: PC::Commitment,
    pub dl_proof: DLProof<F, PC>,
    pub geo_seq_proof: GeoSeqProof<F, PC>,
    pub subset_proof: SubsetProof,
}
