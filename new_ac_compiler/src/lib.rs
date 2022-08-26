use ark_ff::Field;
use ark_relations::r1cs::Matrix;

mod tests;

/// A structure containing the output-final R1CS encoding of an arithmetic circuit. There are `t` input rows,
/// the first is always reserved for the constant 1. All other input rows are for public data, regardless of
/// whether this is a public variable or public constant.
pub struct R1CSfIndex<F> {
    pub number_of_constraints: usize,
    pub number_of_input_rows: usize,
    pub number_of_outputs: usize,
    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}

/// Type of each gate
#[derive(Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub enum GateType {
    Add,
    Mul,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Gate {
    left_index: usize,
    right_index: usize,
    symbol: GateType,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Circuit {
    gates: Vec<Gate>,
    number_of_inputs: usize,
    number_of_outputs: usize,
}

/// A compiler from arithmetic circuit to t-functional triple as described in Construction 2 of the functional
/// commitments paper
pub fn vanilla_ac2tft<F: Field>(circuit: Circuit) -> R1CSfIndex<F> {
    let number_of_constraints = circuit.gates.len() + circuit.number_of_inputs + 1;
    let number_of_input_rows = circuit.number_of_inputs + 1; // this is the `t` value in a t-functional triple
    let number_of_outputs = circuit.number_of_outputs;

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
            }
            GateType::Mul => {
                a_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.left_index));
                b_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 1 + gate.right_index));
            }
        }
    }

    R1CSfIndex {
        number_of_constraints,
        number_of_input_rows,
        number_of_outputs,
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
