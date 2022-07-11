#[cfg(test)]
pub mod tests {
    use ark_relations::r1cs::Matrix;
    use ark_poly::univariate::DensePolynomial;
    use ark_bn254::Fr;
    use crate::{GateType, GateInput, Gate};
    use crate::error::Error;

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
    }

    pub fn sample_gates_0() -> Vec<Gate<F>> {
        // Encodes x^3 + 2x + 5
        // g0: (x, x, *)
        // g1: (g0, x, *)
        // g2: (x, 2, *)
        // g3: (g1, g2, +)
        // g4: (g3, 5, +)

        // s: [*, *, *, +, +]
        // l: [x, g0, g1, g3]
        // r: [x, 2, g2, 5]

        let g0 = Gate::<Fr> {
            left: GateInput::Input(String::from("x")),
            right: GateInput::Input(String::from("x")),
            symbol: GateType::Mul,
            label: String::from("g0"),
        };

        let g1 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g0.clone())),
            right: GateInput::Input(String::from("x")),
            symbol: GateType::Mul,
            label: String::from("g1"),
        };

        let g2 = Gate::<Fr> {
            left: GateInput::Input(String::from("x")),
            right: GateInput::Constant(F::from(2u64)),
            symbol: GateType::Mul,
            label: String::from("g2"),
        };

        let g3 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g1.clone())),
            right: GateInput::Gate(Box::new(g2.clone())),
            symbol: GateType::Add,
            label: String::from("g3"),
        };

        let g4 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g3.clone())),
            right: GateInput::Constant(F::from(5u64)),
            symbol: GateType::Add,
            label: String::from("g4"),
        };

        return vec![g0, g1, g2, g3, g4];
    }

    #[test]
    fn test_gate_input_eq() {
        let gates = sample_gates_0();
        let g0 = &gates[0];
        let g1 = &gates[1];
        let g2 = &gates[2];
        assert_eq!(g0.left, g0.right);
        assert_eq!(g0.left, g1.right);
        assert_eq!(g2.right, GateInput::Constant(F::from(2u64)));
    }

    #[test]
    fn test_into_gate() {
        let gates = sample_gates_0();
        let g0 = &gates[0];
        let g1 = &gates[1];
        assert_eq!(g1.left.into_gate().unwrap(), *g0);
        assert_eq!(g1.right.into_gate(), Err(Error::GateInputNotGate));
    }

    fn print_matrix(
        label: &str,
        matrix: Matrix<F>,
    ) {
        println!("{}:", label);
        for (i, i_vec) in matrix.iter().enumerate() {
            print!("{}: ", i);
            for j_tuple in i_vec.iter() {
                print!("({}, {}), ", j_tuple.0, j_tuple.1);
            }
            println!("");
        }
        println!("");
    }

    fn print_joint_matrix(
        label: &str,
        matrix: Vec<Vec<usize>>,
    ) {
        println!("{}:", label);
        println!("{:?}:", matrix);
        println!("");
    }

    fn print_polynomial(
        label: &str,
        poly: &DensePolynomial<Fr>,
    ) {
        println!("{} coeffs: ", label,);
        for c in poly.coeffs.iter() {
            println!("{}", c);
        }
        println!("");
    }
}
