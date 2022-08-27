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
