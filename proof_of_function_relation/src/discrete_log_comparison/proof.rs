use crate::geo_seq::proof::Proof as GeoProof;
use crate::non_zero_over_k::proof::Proof as NonZeroProof;
use crate::subset_over_k::proof::Proof as SubsetProof;
use ark_ff::PrimeField;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use zero_over_k::zero_over_k::proof::Proof as ZeroProof;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>> {
    // Commitments
    pub s_commit: PC::Commitment,
    pub f_prime_commit: PC::Commitment,
    pub g_prime_commit: PC::Commitment,
    pub s_prime_commit: PC::Commitment,
    pub h_commit: PC::Commitment,

    // Proofs
    pub f_prime_square_proof: ZeroProof<F, PC>,
    pub g_prime_square_proof: ZeroProof<F, PC>,
    pub s_prime_square_proof: ZeroProof<F, PC>,
    pub f_prime_product_proof: ZeroProof<F, PC>,
    pub f_prime_subset_proof: SubsetProof,
    pub g_prime_subset_proof: SubsetProof,
    pub s_prime_subset_proof: SubsetProof,
    pub h_proof: GeoProof<F, PC>,
    pub nzk_f_prime_proof: NonZeroProof<F, PC>,
    pub nzk_g_prime_proof: NonZeroProof<F, PC>,
    pub nzk_s_prime_proof: NonZeroProof<F, PC>,
    pub nzk_s_minus_one_proof: NonZeroProof<F, PC>,
}
