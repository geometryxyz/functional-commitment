use ark_ff::Field;
use ark_relations::r1cs::Matrix;
use std::iter;

pub mod example_circuits;
mod tests;

mod constraint_builder;
mod constraint_syntesizer;
mod error;
mod variable;

/// A structure containing the output-final R1CS encoding of an arithmetic circuit. There are `t` input rows,
/// the first is always reserved for the constant 1. All other input rows are for public data, regardless of
/// whether this is a public variable or public constant.
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

/// Type of an arithmetic gate
#[derive(Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub enum GateType {
    Add,
    Mul,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Gate {
    pub left_index: usize,
    pub right_index: usize,
    pub symbol: GateType,
}

impl Gate {
    pub fn new(left_index: usize, right_index: usize, symbol: GateType) -> Self {
        Self {
            left_index,
            right_index,
            symbol,
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Circuit {
    gates: Vec<Gate>,
    number_of_inputs: usize,
    number_of_outputs: usize,
}

impl Circuit {
    pub fn new(gates: Vec<Gate>, number_of_inputs: usize, number_of_outputs: usize) -> Self {
        Self {
            gates,
            number_of_inputs,
            number_of_outputs,
        }
    }
}

/// A compiler from arithmetic circuit to t-functional triple as described in Construction 2 of the functional
/// commitments paper
pub fn vanilla_ac2tft<F: Field>(circuit: Circuit) -> R1CSfIndex<F> {
    let number_of_constraints = circuit.gates.len() + circuit.number_of_inputs + 1;
    let number_of_input_rows = circuit.number_of_inputs + 1; // this is the `t` value in a t-functional triple
    let number_of_outputs = circuit.number_of_outputs;
    let mut number_of_non_zero_entries = 0;

    let mut a_matrix = empty_matrix(number_of_constraints);
    let mut b_matrix = empty_matrix(number_of_constraints);
    let mut c_matrix = empty_matrix(number_of_constraints);

    for (i, gate) in circuit.gates.iter().enumerate() {
        c_matrix[1 + circuit.number_of_inputs + i]
            .push((F::one(), 1 + circuit.number_of_inputs + i));
        match gate.symbol {
            GateType::Add => {
                a_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1));
                b_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.left_index));
                b_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.right_index));
                number_of_non_zero_entries += 2;
            }
            GateType::Mul => {
                a_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.left_index));
                b_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.right_index));
                number_of_non_zero_entries += 1;
            }
        }
    }

    R1CSfIndex {
        number_of_constraints,
        number_of_input_rows,
        number_of_outputs,
        number_of_non_zero_entries,
        a: a_matrix,
        b: b_matrix,
        c: c_matrix,
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
