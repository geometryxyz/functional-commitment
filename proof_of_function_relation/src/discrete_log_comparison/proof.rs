use crate::commitment::HomomorphicPolynomialCommitment;
use crate::geo_seq::proof::Proof as GeoProof;
use crate::non_zero_over_k::proof::Proof as NonZeroProof;
use crate::subset_over_k::proof::Proof as SubsetProof;
use crate::zero_over_k::proof::Proof as ZeroProof;
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
    pub f_prime_square_proof: Vec<u8>,//ZeroProof<F, PC>,
    pub g_prime_square_proof: Vec<u8>,//ZeroProof<F, PC>,
    pub s_prime_square_proof: Vec<u8>,//ZeroProof<F, PC>,
    pub f_prime_product_proof: Vec<u8>, //ZeroProof<F, PC>,
    pub f_prime_subset_proof: SubsetProof,
    pub g_prime_subset_proof: SubsetProof,
    pub s_prime_subset_proof: SubsetProof,
    pub h_proof: Vec<u8>, //GeoProof<F, PC>,
    pub nzk_f_prime_proof: Vec<u8>, // NonZeroProof<F, PC>,
    pub nzk_g_prime_proof: Vec<u8>, // NonZeroProof<F, PC>,
    pub nzk_s_prime_proof: Vec<u8>, // NonZeroProof<F, PC>,
    pub nzk_s_minus_one_proof: Vec<u8>, // NonZeroProof<F, PC>,
}
