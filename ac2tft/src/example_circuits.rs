use crate::{gates_to_sparse_matrices, Gate, GateInput, GateType, SparseMatrices};
use ark_ff::PrimeField;

pub fn sample_gates_0<F: PrimeField>() -> Vec<Gate<F>> {
    // Encodes x^3 + 2x + 5
    // g0: (x, x, *)
    // g1: (g0, x, *)
    // g2: (x, 2, *)
    // g3: (g1, g2, +)
    // g4: (g3, 5, +)

    // s: [*, *, *, +, +]
    // l: [x, g0, g1, g3]
    // r: [x, 2, g2, 5]

    let g0 = Gate::<F> {
        left: GateInput::Input(String::from("x")),
        right: GateInput::Input(String::from("x")),
        symbol: GateType::Mul,
        label: String::from("g0"),
    };

    let g1 = Gate::<F> {
        left: GateInput::Gate(Box::new(g0.clone())),
        right: GateInput::Input(String::from("x")),
        symbol: GateType::Mul,
        label: String::from("g1"),
    };

    let g2 = Gate::<F> {
        left: GateInput::Input(String::from("x")),
        right: GateInput::Constant(F::from(2u64)),
        symbol: GateType::Mul,
        label: String::from("g2"),
    };

    let g3 = Gate::<F> {
        left: GateInput::Gate(Box::new(g1.clone())),
        right: GateInput::Gate(Box::new(g2.clone())),
        symbol: GateType::Add,
        label: String::from("g3"),
    };

    let g4 = Gate::<F> {
        left: GateInput::Gate(Box::new(g3.clone())),
        right: GateInput::Constant(F::from(5u64)),
        symbol: GateType::Add,
        label: String::from("g4"),
    };

    return vec![g0, g1, g2, g3, g4];
}

pub fn sample_matrices<F: PrimeField>() -> SparseMatrices<F> {
    let gates = sample_gates_0();
    gates_to_sparse_matrices(gates)
}
