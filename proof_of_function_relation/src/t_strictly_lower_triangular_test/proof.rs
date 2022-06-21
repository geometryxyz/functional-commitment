use crate::{
    commitment::HomomorphicPolynomialCommitment, 
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
    pub dl_proof: Vec<u8>,
    pub geo_seq_proof: Vec<u8>,
    pub subset_proof: Vec<u8>,
}
