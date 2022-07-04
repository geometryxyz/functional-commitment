mod tests;
mod error;
use std::fmt;
use ark_ff::PrimeField;
use ark_marlin::ahp::indexer::Matrix;

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
