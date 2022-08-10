use crate::{
    t_diag::proof::Proof as TDiagProof, t_strictly_lower_triangular_test::proof::Proof as TSLTProof,
};
use ark_ff::{PrimeField, SquareRootField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: PrimeField + SquareRootField, PC: AdditivelyHomomorphicPCS<F>> {
    pub a_slt_proof: TSLTProof<F, PC>,
    pub b_slt_proof: TSLTProof<F, PC>,
    pub c_diag_proof: TDiagProof<F, PC>,
}
