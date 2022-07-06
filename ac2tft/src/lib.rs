pub mod example_circuits;

mod tests;
mod error;

use std::fmt;
use std::collections::BTreeMap;
use ark_ff::PrimeField;
use ark_relations::r1cs::{ Matrix, SynthesisError };
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_marlin_fork::ahp::indexer::sum_matrices;
use ark_marlin_fork::ahp::constraint_systems::{
    num_non_zero,
    arithmetize_matrix,
    MatrixArithmetization,
};

/// A compiler from arithmetic gates to t-FT sparse matrices, which can be converted into
/// polynomicals using Marlin's `arithmetize_matrix()` function.

/// Type of each gate
#[derive(Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub enum GateType {
    Add,
    Mul,
}

impl fmt::Display for GateType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GateType::Add => write!(f, "+"),
            GateType::Mul => write!(f, "*"),
        }
    }
}

pub type LabeledInput = String;

/// The left or right input to a gate
#[derive(Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub enum GateInput<F: PrimeField> {
    Constant(F),
    Input(LabeledInput),
    Gate(Box<Gate<F>>),
}

impl<F: PrimeField> GateInput<F> {
    pub fn into_gate(&self) -> Result<Gate<F>, error::Error> {
        match self {
            GateInput::Constant(_) => Err(error::Error::GateInputNotGate),
            GateInput::Input(_) => Err(error::Error::GateInputNotGate),
            GateInput::Gate::<F>(x) => Ok(*x.clone())
        }
    }
}

impl<F: PrimeField> fmt::Display for GateInput<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GateInput::Constant(x) => write!(f, "{}", x),
            GateInput::Input(x) => write!(f, "{}", x),
            GateInput::Gate::<F>(x) => write!(f, "{}", x),
        }
    }
}

/// An arithmetic gate
#[derive(Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub struct Gate<F: PrimeField> {
    pub left: GateInput<F>,
    pub right: GateInput<F>,
    pub symbol: GateType,
    pub label: String,
}

impl<F: PrimeField> fmt::Display for Gate<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // e.g. "g: (l, r, +)"
        //write!(f, "{}: ({}, {}, {})", self.label, self.left, self.right, self.symbol)
        write!(f, "{}: (", self.label)?;
        match &self.left {
            GateInput::Gate::<F>(_) => write!(f, "{}, ", self.left.into_gate().unwrap().label)?,
            GateInput::Constant(x) => write!(f, "{}, ", x)?,
            GateInput::Input(x) => write!(f, "{}, ", x)?,
        };
        match &self.right {
            GateInput::Gate::<F>(_) => write!(f, "{}, ", self.right.into_gate().unwrap().label)?,
            GateInput::Constant(x) => write!(f, "{}, ", x)?,
            GateInput::Input(x) => write!(f, "{} ", x)?,
        };
        write!(f, "{})", self.symbol)
    }
}

pub fn empty_matrix<F: PrimeField>(
    length: usize,
) -> Matrix<F> {
    let mut matrix = vec![];
    for _ in 0..length {
        matrix.push(vec![]);
    }
    matrix
}

pub fn gates_to_sparse_matrices<F: PrimeField>(
    gates: Vec<Gate<F>>,
) -> (
    Matrix<F>,
    Matrix<F>,
    Matrix<F>,
) {
    let mut left_input_map = BTreeMap::<GateInput<F>, usize>::new();
    let mut right_input_map = BTreeMap::<GateInput<F>, usize>::new();

    //println!("Left inputs");
    let mut j = 0;
    for gate in gates.iter() {
        if !left_input_map.contains_key(&gate.left) {
            left_input_map.insert(gate.left.clone(), j);
            //println!("{}, {}", j, gate.left);
            j += 1;
        }
    }
    //println!("");

    j = 0;
    //println!("Right inputs");
    for gate in gates.iter() {
        if !right_input_map.contains_key(&gate.right) {
            right_input_map.insert(gate.right.clone(), j);
            //println!("{}, {}", j, gate.right);
            j += 1;
        }
    }

    let n_i = left_input_map.len() + right_input_map.len();
    let n_o = gates.len();
    let matrix_width = n_i + n_o;

    let mut matrix_a = empty_matrix::<F>(matrix_width);
    let mut matrix_b = empty_matrix::<F>(matrix_width);
    let mut matrix_c = empty_matrix::<F>(matrix_width);
 
    // TODO: check for off-by-one errors here! the paper uses matrices which start from 1, but
    // in Rust we start from 0.
    for (i, gate) in gates.iter().enumerate() {
        let row = n_i + i;
        let l_i = *left_input_map.get(&gate.left).unwrap() + 1;
        let r_i = *right_input_map.get(&gate.right).unwrap() + 1;

        let col_a = match gate.symbol {
            GateType::Add => 0,
            GateType::Mul => l_i,
        };

        // Matrix A
        matrix_a[row].push((F::from(1u64), col_a));

        // Matrix B
        if gate.symbol == GateType::Add {
            if l_i != r_i {
                matrix_b[row].push((F::from(1u64), l_i));
            }
        }
        matrix_b[row].push((F::from(1u64), r_i));

        // Matrix C
        matrix_c[row].push((F::from(1u64), row));
    }

    (matrix_a, matrix_b, matrix_c)
}

pub fn sparse_matrices_to_polys<F: PrimeField>(
    a: Matrix<F>,
    b: Matrix<F>,
    c: Matrix<F>,
    num_constraints: usize,
    num_formatted_input_variables: usize,
) -> MatrixArithmetization<F> {
    let mut a = a.clone();
    let mut b = b.clone();
    let mut c = c.clone();

    let joint_matrix = sum_matrices(&a, &b, &c);
    let num_non_zero_vals = num_non_zero(&joint_matrix);

    let domain_h = GeneralEvaluationDomain::<F>::new(num_constraints)
        .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
    let domain_k = GeneralEvaluationDomain::<F>::new(num_non_zero_vals)
        .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
    let x_domain = GeneralEvaluationDomain::<F>::new(num_formatted_input_variables)
        .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();

    let joint_arith = arithmetize_matrix(
        &joint_matrix,
        &mut a,
        &mut b,
        &mut c,
        domain_k,
        domain_h,
        x_domain,
    );

    joint_arith
}

