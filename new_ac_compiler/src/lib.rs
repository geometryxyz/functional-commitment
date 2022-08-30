use ark_ff::Field;
use ark_relations::r1cs::Matrix;
use std::iter;

pub mod circuit;
pub mod circuit_compiler;
pub mod constraint_builder;
pub mod error;
pub mod example_circuits;
pub mod gate;
pub mod tests;
pub mod variable;

/// A structure containing the output-final R1CS encoding of an arithmetic circuit. There are `t` input rows,
/// the first is always reserved for the constant 1. All other input rows are for public data, regardless of
/// whether this is a public variable or public constant.
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct R1CSfIndex<F: Field> {
    /// Number of constrains (this is also the length of the matrices)
    pub number_of_constraints: usize,

    /// Number of rows that are reserved for the inputs. This is the `t` value in a t-functional triple.
    /// Note that this is the **number of public inputs to the circuit plus 1** (this is by construction)
    pub number_of_input_rows: usize,

    /// Number of outputs
    pub number_of_outputs: usize,

    /// The maximum number of non-zero entries
    pub number_of_non_zero_entries: usize,

    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}

impl<F: Field> R1CSfIndex<F> {
    /// Iterate through the matrices of the index: A, B, C
    pub fn iter_matrices(&self) -> impl Iterator<Item = &Matrix<F>> {
        iter::once(&self.a)
            .chain(iter::once(&self.b))
            .chain(iter::once(&self.c))
    }
}

fn empty_matrix<F: Field>(length: usize) -> Matrix<F> {
    let mut matrix = vec![];
    for _ in 0..length {
        matrix.push(vec![]);
    }
    matrix
}

#[macro_export]
/// Print a Matrix
macro_rules! printmatrix {
    ($matrix:expr) => {
        for (row_index, row) in $matrix.iter().enumerate() {
            for (value, col_index) in row {
                println!("row {}, col {}: {}", row_index, col_index, value)
            }
        }
    };
}

#[macro_export]
/// Print a Matrix
macro_rules! slt_test {
    ($matrix:expr, $num_of_pub_inputs_plus_one:expr) => {
        for (row_index, row) in $matrix.iter().enumerate() {
            for (_, col_index) in row {
                assert!(row_index >= $num_of_pub_inputs_plus_one);
                assert!(row_index > *col_index);
            }
        }
    };
}

#[macro_export]
/// Print a Matrix
macro_rules! diag_test {
    ($matrix:expr) => {
        for (row_index, row) in $matrix.iter().enumerate() {
            for (_, col_index) in row {
                assert_eq!(row_index, *col_index);
            }
        }
    };
}
