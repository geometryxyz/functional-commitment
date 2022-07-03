#[cfg(test)]
mod tests {

    use ark_bn254::Fr;
    use crate::{GateType, GateInput, Gate};

    type F = Fr;

    #[test]
    fn test_gate_formatting() {
        // g0: (1, 1, +)
        let gate = Gate {
            left: GateInput::Constant(F::from(1u64)),
            right: GateInput::Constant(F::from(1u64)),
            symbol: GateType::Add,
            label: String::from("g0"),
        };

        assert_eq!(
            format!("{}", gate),
            "g0: (Fp256 \"(0000000000000000000000000000000000000000000000000000000000000001)\", \
            Fp256 \"(0000000000000000000000000000000000000000000000000000000000000001)\", +)"
        );
        // 
    }
}
