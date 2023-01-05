#[cfg(test)]
mod tests {
    use crate::variable::VariableType;
    use crate::{
        circuit::Circuit,
        circuit_compiler::{CircuitCompiler, VanillaCompiler},
        constraint_builder::ConstraintBuilder,
        diag_test,
        error::Error,
        gate::GateType,
        slt_test,
    };
    use ark_bn254::Fr;
    use ark_ff::Zero;

    type F = Fr;

    fn circuit_test_template<Func>(constraints: Func)
    where
        Func: FnOnce(&mut ConstraintBuilder<F>) -> Result<(), Error>,
    {
        let mut cb = ConstraintBuilder::<F>::new();

        let synthesized_circuit = Circuit::synthesize(constraints, &mut cb).unwrap();
        let (_r1csf_index_from_synthesized, a, b, c) =
            VanillaCompiler::<F>::ac2tft(&synthesized_circuit);

        slt_test!(a, r1csf_index_from_synthesized.number_of_input_rows);
        slt_test!(b, r1csf_index_from_synthesized.number_of_input_rows);
        diag_test!(c);

        // Perform matrix multiplications
        let inner_prod_fn = |row: &[(F, usize)]| {
            let mut acc = F::zero();
            for &(_, i) in row {
                acc += cb.assignment[i];
            }
            acc
        };

        let z_a: Vec<F> = a.iter().map(|row| inner_prod_fn(row)).collect();
        let z_b: Vec<F> = b.iter().map(|row| inner_prod_fn(row)).collect();
        let z_c: Vec<F> = c.iter().map(|row| inner_prod_fn(row)).collect();

        assert_eq!(z_a.len(), z_b.len());
        assert_eq!(z_b.len(), z_c.len());
        for ((&za_i, &zb_i), &zc_i) in z_a.iter().zip(z_b.iter()).zip(z_c.iter()) {
            assert_eq!(za_i * zb_i, zc_i);
        }
    }

    #[test]
    fn test_simple_circuit() {
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", F::from(7u64))?;

            let x_square = cb.enforce_constraint(&x, &x, GateType::Mul, VariableType::Witness)?;
            let x_cube =
                cb.enforce_constraint(&x_square, &x, GateType::Mul, VariableType::Witness)?;

            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul, VariableType::Witness)?;
            let x_qubed_plus_2x =
                cb.enforce_constraint(&x_cube, &two_x, GateType::Add, VariableType::Witness)?;

            let _ = cb.enforce_constraint(
                &x_qubed_plus_2x,
                &five,
                GateType::Add,
                VariableType::Output,
            )?;

            Ok(())
        };

        circuit_test_template(constraints)
    }
}
