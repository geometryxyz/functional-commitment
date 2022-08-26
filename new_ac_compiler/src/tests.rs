#[cfg(test)]
mod tests {
    use crate::{circuit::Circuit, constraint_builder::ConstraintBuilder, gate::GateType};
    use ark_bn254::Fr;

    type F = Fr;

    // #[test]
    // fn test_simple_circuit_1() {
    //     let mut cb = ConstraintBuilder::<F>::new(1);
    //     let circuit = Circuit::synthesize(|&mut constraint_builder| {
    //         let x = cb.new_input_variable("x", F::from(7u64))?;
    //         let two = cb.new_input_variable("two", F::from(2u64))?;
    //         let five = cb.new_input_variable("five", F::from(5u64))?;
    
    //         let x_sqare = cb.enforce_constraint(&x, &x, GateType::Mul)?;
    //         let x_qube = cb.enforce_constraint(&x, &x_sqare, GateType::Mul)?;
    
    //         let two_x = cb.enforce_constraint(&two, &x, GateType::Mul)?;
    //         let x_qubed_plus_2x = cb.enforce_constraint(&x_qube, &two_x, GateType::Add)?;
    
    //         let _ = cb.enforce_constraint(&x_qubed_plus_2x, &five, GateType::Add)?;
    
    //         Ok(())
    //     }, &mut cb).unwrap();
    // }
}
