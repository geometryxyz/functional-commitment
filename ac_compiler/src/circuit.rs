use ark_ff::Field;

use crate::{constraint_builder::ConstraintBuilder, error::Error, gate::Gate};

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

    pub fn from_constraint_builder<F: Field>(cb: &ConstraintBuilder<F>) -> Self {
        Self {
            gates: cb.gates.clone(),
            number_of_inputs: cb.number_of_inputs,
            number_of_outputs: cb.number_of_outputs,
        }
    }

    pub fn synthesize<Func, F>(f: Func, cb: &mut ConstraintBuilder<F>) -> Result<Self, Error>
    where
        F: Field,
        Func: FnOnce(&mut ConstraintBuilder<F>) -> Result<(), Error>,
    {
        f(cb)?;
        cb.finalize();
        Ok(Self::from_constraint_builder(&cb))
    }
}
