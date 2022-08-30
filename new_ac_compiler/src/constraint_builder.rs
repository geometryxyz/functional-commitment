use std::{collections::BTreeMap, marker::PhantomData};

use crate::{
    error::Error,
    gate::{Gate, GateType},
    variable::{Variable, VariableType},
};
use ark_ff::Field;

pub struct ConstraintBuilder<F: Field> {
    pub(crate) gates: Vec<Gate>,
    label_to_var_index: BTreeMap<String, usize>,
    curr_index: usize,
    pub(crate) number_of_inputs: usize,
    pub(crate) number_of_ouputs: usize,
    pub assignment: Vec<F>,
    _f: PhantomData<F>,
}

impl<F: Field> ConstraintBuilder<F> {
    pub fn new() -> Self {
        Self {
            gates: Vec::new(),
            label_to_var_index: BTreeMap::new(),
            curr_index: 1, // we reserve first input to be dummy selector constraint for addition
            number_of_inputs: 0,
            number_of_ouputs: 0,
            assignment: vec![F::one()], // we add 1 to be the value of dummy selector
            _f: PhantomData,
        }
    }

    pub fn new_input_variable(&mut self, label: &str, value: F) -> Result<Variable<F>, Error> {
        self.register_new_var(label, value, VariableType::Input)
    }

    pub fn enforce_constraint(
        &mut self,
        lhs: &Variable<F>,
        rhs: &Variable<F>,
        constraint_type: GateType,
        variable_type: VariableType,
    ) -> Result<Variable<F>, Error> {
        let lhs_index = match self.label_to_var_index.get(&lhs.label) {
            Some(index) => Ok(*index),
            None => Err(Error::VarMissing(format!(
                "Var with label {} doesn't exists",
                lhs.label
            ))),
        }?;

        let rhs_index = match self.label_to_var_index.get(&rhs.label) {
            Some(index) => Ok(*index),
            None => Err(Error::VarMissing(format!(
                "Var with label {} doesn't exists",
                lhs.label
            ))),
        }?;

        let new_value = match constraint_type {
            GateType::Add => {
                self.gates.push(Gate {
                    left_index: lhs_index,
                    right_index: rhs_index,
                    symbol: GateType::Add,
                });
                lhs.value + rhs.value
            }
            GateType::Mul => {
                self.gates.push(Gate {
                    left_index: lhs_index,
                    right_index: rhs_index,
                    symbol: GateType::Mul,
                });
                lhs.value * rhs.value
            }
        };

        // for now automatically assign wtns_prefix to intermidiate valies
        let new_label = format!("w_{}", self.curr_index);
        self.register_new_var(&new_label, new_value, variable_type)
    }

    fn register_new_var(
        &mut self,
        label: &str,
        value: F,
        variable_type: VariableType,
    ) -> Result<Variable<F>, Error> {
        // we don't allow vars with same labels
        if self.label_to_var_index.contains_key(label.into()) {
            return Err(Error::VarAlreadyExists(format!(
                "Var with label {} already exists",
                label
            )));
        }

        let var = Variable {
            label: label.into(),
            value,
            variable_type,
        };

        self.label_to_var_index
            .insert(label.into(), self.curr_index);

        // whenever new var is added, global index is incremented by one
        // this follows from the fact that each gate introduces another intermidiate var
        self.curr_index += 1;

        match variable_type {
            VariableType::Input => {
                self.number_of_inputs += 1;
            }
            VariableType::Witness => {}
            VariableType::Output => {
                self.number_of_ouputs += 1;
            }
        };

        self.assignment.push(value);

        Ok(var)
    }
}
