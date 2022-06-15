use crate::{
    commitment::HomomorphicPolynomialCommitment, t_diag::proof::Proof as TDiagProof,
    t_strictly_lower_triangular_test::proof::Proof as TSLTProof,
};
use ark_ff::{PrimeField, SquareRootField};

pub struct Proof<F: PrimeField + SquareRootField, PC: HomomorphicPolynomialCommitment<F>> {
    pub a_slt_proof: TSLTProof<F, PC>,
    pub b_slt_proof: TSLTProof<F, PC>,
    pub c_diag_proof: TDiagProof<F, PC>,
}
