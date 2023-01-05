use ac_compiler::R1CSfIndex;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use derivative::Derivative;

use ark_std::io::{Read, Write};

use super::constraint_systems::MatrixArithmetization;

/// Represents a matrix.
pub type Matrix<F> = Vec<Vec<(F, usize)>>;

#[derive(Derivative)]
#[derivative(Clone(bound = "F: PrimeField"))]
/// The indexed version of the constraint system.
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Index<F: PrimeField> {
    /// Information about the index.
    pub index_info: R1CSfIndex,

    /// The A matrix arithmetization
    pub a_arith: MatrixArithmetization<F>,

    /// The B matrix arithmetization
    pub b_arith: MatrixArithmetization<F>,

    /// The C matrix arithmetization
    pub c_arith: MatrixArithmetization<F>,

    /// tmp store matrices
    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}
