use crate::gate::Gate;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Circuit {
    pub gates: Vec<Gate>,
    pub number_of_inputs: usize,
    pub number_of_outputs: usize,
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
