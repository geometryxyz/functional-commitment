use crate::commitment::HomomorphicPolynomialCommitment;
use crate::zero_over_k::proof::Proof as Z_Proof;
use ark_ff::PrimeField;

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    pub z_proof: Z_Proof<F, PC>,
    pub f_commit: PC::Commitment,
    pub opening_proof: PC::BatchProof,
}
