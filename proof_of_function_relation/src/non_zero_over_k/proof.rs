use crate::{commitment::HomomorphicPolynomialCommitment, zero_over_k};
use ark_ff::PrimeField;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    pub g_commit: PC::Commitment,
    pub zero_over_k_proof: zero_over_k::proof::Proof<F, PC>,
}
