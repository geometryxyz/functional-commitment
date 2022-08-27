use crate::R1CSfIndex;
use ark_ff::Field;
use std::marker::PhantomData;

use crate::{circuit::Circuit, empty_matrix, gate::GateType};

/// Given: an arithmetic circuit with ng gates, ni inputs, and no <= ng outputs, where gates are triples of (left_input_index, right_input_index, (add/mul))
/// Produces: An index for R_R1CS-f(ng + ni + 1, ni + 1, no)
pub trait CircuitCompiler<F: Field> {
    fn ac2tft(circuit: &Circuit) -> R1CSfIndex<F>;
}

pub struct VanillaCompiler<F: Field> {
    _f: PhantomData<F>,
}

impl<F: Field> CircuitCompiler<F> for VanillaCompiler<F> {
    fn ac2tft(circuit: &Circuit) -> R1CSfIndex<F> {
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
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.right_index));
                    number_of_non_zero_entries += 2;
                }
                GateType::Mul => {
                    a_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.right_index));
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
}
