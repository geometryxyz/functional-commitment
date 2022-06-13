use crate::commitment::HomomorphicPolynomialCommitment;
use crate::zero_over_k::proof::Proof as ZProof;
use ark_ff::PrimeField;

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    pub z_proof: ZProof<F, PC>,
    pub opening_proof: PC::BatchProof,
}
