// use ark_ff::Field;

// use crate::{
//     error::Error, constraint_builder::ConstraintBuilder, circuit::Circuit, gate::GateType
// };

// pub trait ConstraintSynthesizer<F: Field> {
//     fn synthesize(cb: &mut ConstraintBuilder<F>) -> Result<Circuit, Error>;
// }

// pub struct SimpleCirc {}

// impl<F: Field> ConstraintSynthesizer<F> for SimpleCirc {
//     fn synthesize(cb: &mut ConstraintBuilder<F>) -> Result<Circuit, Error> {
//         let x = cb.new_input_variable("x", F::from(7u64))?;
//         let two = cb.new_input_variable("two", F::from(2u64))?;
//         let five = cb.new_input_variable("five", F::from(5u64))?;

//         let x_sqare = cb.enforce_constraint(&x, &x, GateType::Mul)?;
//         let x_qube = cb.enforce_constraint(&x, &x_sqare, GateType::Mul)?;

//         let two_x = cb.enforce_constraint(&two, &x, GateType::Mul)?;
//         let x_qubed_plus_2x = cb.enforce_constraint(&x_qube, &two_x, GateType::Add)?;

//         let _ = cb.enforce_constraint(&x_qubed_plus_2x, &five, GateType::Add)?;

//         Ok(Circuit::from_constraint_builder(cb))
//     }
// }