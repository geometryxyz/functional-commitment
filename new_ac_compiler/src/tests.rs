#[cfg(test)]
mod tests {
    use crate::{circuit::Circuit, constraint_builder::ConstraintBuilder, gate::GateType, error::Error, circuit_compiler::{VanillaCompiler, CircuitCompiler}, printmatrix, example_circuits::sample_circuit_2};
    use ark_bn254::Fr;

    type F = Fr;

    #[test]
    fn test_simple_circuit_1() {

        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", F::from(7u64))?;
    
            let x_sqare = cb.enforce_constraint(&x, &x, GateType::Mul)?;
            let x_qube = cb.enforce_constraint(&x_sqare, &x, GateType::Mul)?;
    
            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul)?;
            let x_qubed_plus_2x = cb.enforce_constraint(&x_qube, &two_x, GateType::Add)?;
    
            let _ = cb.enforce_constraint(&x_qubed_plus_2x, &five, GateType::Add)?;
    
            Ok(())
        };

        let mut cb = ConstraintBuilder::<F>::new(1);

        let circuit = Circuit::synthesize(constraints, &mut cb).unwrap();

        let r1csf_index = VanillaCompiler::<F>::ac2tft(&circuit);
        printmatrix!(r1csf_index.a);
        println!("====================================");
        printmatrix!(r1csf_index.b);
        println!("====================================");
        printmatrix!(r1csf_index.c);

        let circuit = sample_circuit_2();
        let r1csf_index = VanillaCompiler::<F>::ac2tft(&circuit);


        println!("====================================");
        println!("====================================");
        println!("====================================");

        printmatrix!(r1csf_index.a);
        println!("====================================");
        printmatrix!(r1csf_index.b);
        println!("====================================");
        printmatrix!(r1csf_index.c);

    }
}
