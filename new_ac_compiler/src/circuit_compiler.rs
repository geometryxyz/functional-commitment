use ark_ff::Field;
use ark_relations::r1cs::Matrix;
use std::{iter, marker::PhantomData};

use crate::{circuit::Circuit, empty_matrix, gate::GateType};

/// A structure containing the output-final R1CS encoding of an arithmetic circuit. There are `t` input rows,
/// the first is always reserved for the constant 1. All other input rows are for public data, regardless of
/// whether this is a public variable or public constant.
pub struct R1CSfIndex<F: Field> {
    /// Number of constrains (this is also the length of the matrices)
    pub number_of_constraints: usize,

    /// Number of rows that are reserved for the inputs. This is the `t` value in a t-functional triple.
    /// Note that this is the **number of public inputs to the circuit plus 1** (this is by construction)
    pub number_of_input_rows: usize,

    /// Number of outputs
    pub number_of_outputs: usize,

    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}

impl<F: Field> R1CSfIndex<F> {
    /// Iterate through the matrices of the index: A, B, C
    pub fn iter_matrices(&self) -> impl Iterator<Item = &Matrix<F>> {
        iter::once(&self.a)
            .chain(iter::once(&self.b))
            .chain(iter::once(&self.c))
    }
}

/// Given: an arithmetic circuit with ng gates, ni inputs, and no <= ng outputs, where gates are triples of (left_input_index, right_input_index, (add/mul))
/// Produces: An index for R_R1CS-f(ng + ni + 1, ni + 1, no)
pub trait CircuitCompiler<F: Field> {
    fn ac2tft(circuit: Circuit) -> R1CSfIndex<F>;
}

pub struct VanillaCompiler<F: Field> {
    _f: PhantomData<F>,
}

impl<F: Field> CircuitCompiler<F> for VanillaCompiler<F> {
    fn ac2tft(circuit: Circuit) -> R1CSfIndex<F> {
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
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.right_index));
                }
                GateType::Mul => {
                    a_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.left_index));
                    b_matrix[1 + circuit.number_of_inputs + i]
                        .push((F::one(), 1 + gate.right_index));
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
}
