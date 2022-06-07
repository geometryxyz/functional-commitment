use ark_ff::PrimeField;
use crate::commitment::HomomorphicPolynomialCommitment;
use crate::zero_over_k::proof::{Proof as Z_Proof};

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    pub z_proof: Z_Proof::<F, PC>,
    pub seq: Vec<F>,
    pub r: F,
    pub a_vec: Vec<F>,
    pub c_vec: Vec<usize>,
}
