use crate::{
    circuit::Circuit,
    gate::{Gate, GateType},
};

/// Encode the circuit x^2 + 5
pub fn sample_circuit_1() -> Circuit {
    // index 0: x [public]
    // index 1: 5 [public]
    // index 2: x^2 [witness]
    let gate_0 = Gate::new(0, 0, GateType::Mul);
    let gate_1 = Gate::new(2, 1, GateType::Add);

    Circuit::new(vec![gate_0, gate_1], 2, 1)
}

/// Encode the circuit x^3 + 2x + 5
pub fn sample_circuit_2() -> Circuit {
    // index 0: 2 [public]
    // index 1: 5 [public]
    // index 2: x [public]
    // index 3: x^2 [witness]
    // index 4: x^3 [witness]
    // index 5: 2x [witness]
    // index 6: x^3 + 2x [witness]

    let gate_0 = Gate::new(2, 2, GateType::Mul); // produce x^2
    let gate_1 = Gate::new(3, 2, GateType::Mul); // produce x^3
    let gate_2 = Gate::new(0, 2, GateType::Mul); // produce 2x
    let gate_3 = Gate::new(4, 5, GateType::Add); // produce x^3 + 2x
    let gate_4 = Gate::new(6, 1, GateType::Add); // produce the output

    Circuit::new(vec![gate_0, gate_1, gate_2, gate_3, gate_4], 3, 1)
}
