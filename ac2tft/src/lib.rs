mod tests;
use std::fmt;
use ark_ff::PrimeField;

/// A compiler from arithmetic gates to t-FT sparse matrices, which can be converted into
/// polynomicals using Marlin's `arithmetize_matrix()` function.

/// Type of each gate
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
pub enum GateInput<F: PrimeField> {
    Constant(F),
    Input(LabeledInput),
    Gate(Box<Gate<F>>),
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

/// A Gate
pub struct Gate<F: PrimeField> {
    pub left: GateInput<F>,
    pub right: GateInput<F>,
    pub symbol: GateType,
    pub label: String,
}

impl<F: PrimeField> fmt::Display for Gate<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: ({}, {}, {})", self.label, self.left, self.right, self.symbol)
    }
}
