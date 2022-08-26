use crate::{Circuit, Gate, GateType};

/// Encode the circuit x^2 + 5
pub fn sample_circuit_1() -> Circuit {
    // input 0: x
    // input 1: 5
    // inpt 2 (computed from others): x^2
    let gate_0 = Gate::new(0, 0, GateType::Mul);
    let gate_1 = Gate::new(2, 1, GateType::Add);

    Circuit::new(vec![gate_0, gate_1], 2, 1)
}
