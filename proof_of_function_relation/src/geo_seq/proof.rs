use ark_ff::PrimeField;
use zero_over_k::{commitment::AdditivelyHomomorphicPCS, zero_over_k::proof::Proof as ZProof};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    pub z_proof: ZProof<F, PC>,
    pub opening_proof: PC::BatchProof,
}
