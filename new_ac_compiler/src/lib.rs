use ark_ff::Field;
use ark_relations::r1cs::Matrix;
pub mod example_circuits;
mod tests;

mod circuit;
mod circuit_compiler;
mod constraint_builder;
mod constraint_syntesizer;
mod error;
mod gate;
mod variable;

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
