use ark_ff::PrimeField;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    // commitments
    pub m_commitments: Vec<PC::Commitment>,
    pub r_commitments: Vec<PC::Commitment>,
    pub q1_commit: PC::Commitment,

    // evaluations
    pub q1_eval: F,
    pub q2_eval: F,
    pub h_prime_evals: Vec<F>,
    pub m_evals: Vec<F>,

    // opening proof
    pub opening_proof: PC::BatchProof,
}
