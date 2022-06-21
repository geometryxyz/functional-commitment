use crate::{
    commitment::HomomorphicPolynomialCommitment, discrete_log_comparison::proof::Proof as DLProof,
    geo_seq::proof::Proof as GeoSeqProof, subset_over_k::proof::Proof as SubsetProof,
};
use ark_ff::{PrimeField, SquareRootField};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F, PC>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
{
    pub h_commit: PC::Commitment,
    pub dl_proof: Vec<u8>, // DLProof<F, PC>,
    pub geo_seq_proof: Vec<u8>, // GeoSeqProof<F, PC>,
    pub subset_proof: SubsetProof,
}
