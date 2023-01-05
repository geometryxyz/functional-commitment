use crate::{Matrix, R1CSfIndex};
use ark_ff::PrimeField;
// use ark_marlin::ahp::indexer::Matrix;
use std::{cmp::max, marker::PhantomData};

use crate::{circuit::Circuit, empty_matrix, gate::GateType};

/// Given: an arithmetic circuit with ng gates, ni inputs, and no <= ng outputs, where gates are triples of (left_input_index, right_input_index, (add/mul))
/// Produces: An index for R_R1CS-f(ng + ni + 1, ni + 1, no)
pub trait CircuitCompiler<F: PrimeField> {
    fn ac2tft(circuit: &Circuit) -> (R1CSfIndex, Matrix<F>, Matrix<F>, Matrix<F>);
}

pub struct VanillaCompiler<F: PrimeField> {
    _f: PhantomData<F>,
}

impl<F: PrimeField> CircuitCompiler<F> for VanillaCompiler<F> {
    fn ac2tft(circuit: &Circuit) -> (R1CSfIndex, Matrix<F>, Matrix<F>, Matrix<F>) {
        let number_of_constraints = circuit.gates.len() + circuit.number_of_inputs + 1;
        let number_of_input_rows = circuit.number_of_inputs + 1; // this is the `t` value in a t-functional triple
        let number_of_outputs = circuit.number_of_outputs;

        // num of non zero entires should be max of nom of non zero in a, b, c
        let mut number_of_non_zero_entries_a = 0;
        let mut number_of_non_zero_entries_b = 0;
        let mut number_of_non_zero_entries_c = 0;

        let mut a_matrix = empty_matrix(number_of_constraints);
        let mut b_matrix = empty_matrix(number_of_constraints);
        let mut c_matrix = empty_matrix(number_of_constraints);

        // 0 + var_index is intentionaly left to notate that indices from the compiler are already shifted by one (because of the dummy selector)
        for (i, gate) in circuit.gates.iter().enumerate() {
            c_matrix[1 + circuit.number_of_inputs + i]
                .push((F::one(), 1 + circuit.number_of_inputs + i));
            number_of_non_zero_entries_c += 1;
            match gate.symbol {
                GateType::Add => {
                    a_matrix[1 + circuit.number_of_inputs + i].push((F::one(), 0));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 0 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 0 + gate.right_index));
                    number_of_non_zero_entries_a += 1;
                    number_of_non_zero_entries_b += 2;
                }
                GateType::Mul => {
                    a_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 0 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 0 + gate.right_index));
                    number_of_non_zero_entries_a += 1;
                    number_of_non_zero_entries_b += 1;
                }
            }
        }

        let number_of_non_zero_entries = max(
            number_of_non_zero_entries_a,
            max(number_of_non_zero_entries_b, number_of_non_zero_entries_c),
        );

        let index_info = R1CSfIndex {
            number_of_constraints,
            number_of_input_rows,
            number_of_outputs,
            number_of_non_zero_entries,
        };

        (index_info, a_matrix, b_matrix, c_matrix)
    }
}
