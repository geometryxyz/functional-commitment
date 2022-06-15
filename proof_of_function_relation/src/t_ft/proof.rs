use crate::{
    commitment::HomomorphicPolynomialCommitment,
    t_strictly_lower_triangular_test::proof::{ Proof as TSLTProof },
    t_diag::proof::{ Proof as TDiagProof },
};
use ark_ff::{PrimeField, SquareRootField};

pub struct Proof<
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>
> {
    pub t_slt_proof: TSLTProof<F, PC>,
    pub t_diag_proof: TDiagProof<F, PC>,
}
