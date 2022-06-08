use crate::commitment::HomomorphicPolynomialCommitment;
use crate::geo_seq::proof::Proof as Geo_Proof;
use crate::non_zero_over_k::proof::Proof as NZ_Proof;
use crate::zero_over_k::proof::Proof as Z_Proof;
use ark_ff::PrimeField;

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    // Commitments
    pub s_commit: PC::Commitment,
    pub f_prime_commit: PC::Commitment,
    pub g_prime_commit: PC::Commitment,
    pub s_prime_commit: PC::Commitment,
    pub h_commit: PC::Commitment,

    // Proofs
    pub f_prime_square_proof: Z_Proof<F, PC>,
    pub g_prime_square_proof: Z_Proof<F, PC>,
    pub s_prime_square_proof: Z_Proof<F, PC>,
    pub f_prime_product_proof: Z_Proof<F, PC>,
    pub h_proof: Geo_Proof<F, PC>,
    pub nzk_f_prime_proof: NZ_Proof<F, PC>,
    pub nzk_g_prime_proof: NZ_Proof<F, PC>,
    pub nzk_s_prime_proof: NZ_Proof<F, PC>,
    pub nzk_s_minus_one_proof: NZ_Proof<F, PC>,
}
