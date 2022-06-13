use crate::{
    commitment::HomomorphicPolynomialCommitment,
    discrete_log_comparison::{proof::Proof as DLProof},
    subset_over_k::proof::Proof as SubsetProof,
    geo_seq::{proof::Proof as GeoSeqProof},
};
use ark_ff::{PrimeField, SquareRootField};
use ark_poly_commit::LabeledCommitment;

pub struct Proof<F, PC>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
{
    pub dl_proof: DLProof<F, PC>,
    pub geo_seq_proof: GeoSeqProof<F, PC>,
    pub subset_proof: SubsetProof
}
