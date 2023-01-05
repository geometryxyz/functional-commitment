use ac_compiler::constraint_builder::ConstraintBuilder;
use ac_compiler::error::Error;
use ac_compiler::gate::GateType;
use ac_compiler::variable::VariableType;
use ark_ff::fields::{Field, PrimeField};
use ac_compiler::variable::Variable;

pub fn build_mux1_circuit<F: Field>(
    cb: &mut ConstraintBuilder<F>,
    a_val: F,
    b_val: F,
    c_val: F,
) -> Result<(), Error> {
    let a = cb.new_input_variable("a", a_val)?;
    let b = cb.new_input_variable("b", b_val)?;
    let c = cb.new_input_variable("c", c_val)?;
    let minus_one = cb.new_input_variable("minus_one", F::zero() - F::one())?;

    // b * c
    let bc = cb.enforce_constraint(&b, &c, GateType::Mul, VariableType::Witness)?;

    // a * c
    let ac = cb.enforce_constraint(&a, &c, GateType::Mul, VariableType::Witness)?;

    // - (a * c)
    let neg_ac = cb.enforce_constraint(&ac, &minus_one, GateType::Mul, VariableType::Witness)?;

    // (b * c) - (a * c)
    let bcac = cb.enforce_constraint(&bc, &neg_ac, GateType::Add, VariableType::Witness)?;

    // (b * c) - (a * c) + a
    let _ = cb.enforce_constraint(&bcac, &a, GateType::Add, VariableType::Output)?;

    Ok(())
}

pub fn enforce_square<F: Field>(
    cb: &mut ConstraintBuilder<F>,
    x: &Variable<F>,
) -> Result<Variable<F>, Error> {

    cb.enforce_constraint(&x, &x, GateType::Mul, VariableType::Witness)
}

pub fn enforce_x_7<F: Field>(
    cb: &mut ConstraintBuilder<F>,
    t: &Variable<F>,
) -> Result<Variable<F>, Error> {
    let t2 = enforce_square(cb, &t)?;
    let t4 = enforce_square(cb, &t2)?;
    let t6 = cb.enforce_constraint(&t2, &t4, GateType::Mul, VariableType::Witness)?;
    let t7 = cb.enforce_constraint(&t, &t6, GateType::Mul, VariableType::Witness)?;
    Ok(t7)
}

// A circuit that composes enforce_square() such that the output = x^4
pub fn build_x4_circuit<F: Field>(
    cb: &mut ConstraintBuilder<F>,
    x_val: F,
) -> Result<(), Error> {
    let one = cb.new_input_variable("one", F::one())?;
    let x = cb.new_input_variable("x", x_val)?;

    let x2 = enforce_square(cb, &x)?;
    let x4 = enforce_square(cb, &x2)?;

    let _ = cb.enforce_constraint(&x4, &one, GateType::Mul, VariableType::Output)?;

    Ok(())
}


/// The mimc7 circuit where nRounds = 2 and k = 2
pub fn build_mimc7_circuit<F: PrimeField>(
    cb: &mut ConstraintBuilder<F>,
    x_val: F,
    c: Vec<F>,
) -> Result<(), Error> {
    let n_rounds = 2;

    let x = cb.new_input_variable("x", x_val)?;
    let k = cb.new_input_variable("k", F::from(2u64))?;

    let mut c_inputs: Vec::<Variable<F>> = vec![];
    for (i, val) in c.iter().enumerate() {
        let input = cb.new_input_variable(format!("c{}", i).as_str(), *val)?;
        c_inputs.push(input);
    }

    let mut r;
    let t = cb.enforce_constraint(&x, &k, GateType::Add, VariableType::Witness)?;
    r = enforce_x_7(cb, &t)?;
    for i in 1..n_rounds {
        let r_plus_k = cb.enforce_constraint(&r, &k, GateType::Add, VariableType::Witness)?;
        let t = cb.enforce_constraint(&r_plus_k, &c_inputs[i], GateType::Add, VariableType::Witness)?;
        r = enforce_x_7(cb, &t)?;
    }

    let _ = cb.enforce_constraint(&r, &k, GateType::Add, VariableType::Output)?;

    Ok(())
}