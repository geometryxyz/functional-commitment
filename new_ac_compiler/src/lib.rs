use ark_ff::Field;
use ark_relations::r1cs::Matrix;
use circuit::Circuit;
use constraint_builder::ConstraintBuilder;
use error::Error;

pub mod example_circuits;
mod tests;

mod circuit;
mod circuit_compiler;
mod constraint_builder;
mod constraint_syntesizer;
mod error;
mod gate;
mod variable;

pub fn synthesize<Func, F>(f: Func, number_of_outputs: usize) -> Result<(), Error>
where
    F: Field,
    Func: FnOnce(&mut ConstraintBuilder<F>) -> Result<(), Error>,
{
    let mut cb = ConstraintBuilder::new(number_of_outputs);
    f(&mut cb)?;

    let circuit = Circuit::from_constraint_builder(&cb);
    Ok(())
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
