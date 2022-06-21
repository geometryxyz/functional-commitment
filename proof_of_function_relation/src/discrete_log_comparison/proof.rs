use crate::commitment::HomomorphicPolynomialCommitment;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    // Commitments
    pub s_commit: PC::Commitment,
    pub f_prime_commit: PC::Commitment,
    pub g_prime_commit: PC::Commitment,
    pub s_prime_commit: PC::Commitment,
    pub h_commit: PC::Commitment,

    // Proofs
    pub f_prime_square_proof: Vec<u8>,
    pub g_prime_square_proof: Vec<u8>,
    pub s_prime_square_proof: Vec<u8>,
    pub f_prime_product_proof: Vec<u8>,
    pub f_prime_subset_proof: Vec<u8>,
    pub g_prime_subset_proof: Vec<u8>,
    pub s_prime_subset_proof: Vec<u8>,
    pub h_proof: Vec<u8>,
    pub nzk_f_prime_proof: Vec<u8>,
    pub nzk_g_prime_proof: Vec<u8>,
    pub nzk_s_prime_proof: Vec<u8>,
    pub nzk_s_minus_one_proof: Vec<u8>,
}
