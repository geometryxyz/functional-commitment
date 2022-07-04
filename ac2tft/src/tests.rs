#[cfg(test)]
mod tests {
    use ark_relations::r1cs::{
        Matrix,
        SynthesisError,
    };
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly::univariate::DensePolynomial;
    use ark_marlin::ahp::indexer::{
        sum_matrices
    };

    use ark_marlin::ahp::constraint_systems::{
        num_non_zero,
        arithmetize_matrix,
        MatrixArithmetization,
    };
    use std::collections::BTreeMap;
    use ark_bn254::Fr;
    use crate::{GateType, GateInput, Gate, empty_matrix};
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

    fn sample_gates_0() -> Vec<Gate<F>> {
        // Encodes x^3 + 2x + 5
        /*
            TODO: rename these to g0 .. g4
            g1: (x, x, *)
            g2: (g1, x, *)
            g3: (x, 2, *)
            g4: (g2, g3, +)
            g5: (g4, 5, +)

            s: [*, *, *, +, +]
            l: [x, g1, g2, g4]
            r: [x, 2, g3, 5]
        */

        let g1 = Gate::<Fr> {
            left: GateInput::Input(String::from("x")),
            right: GateInput::Input(String::from("x")),
            symbol: GateType::Mul,
            label: String::from("g1"),
        };

        let g2 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g1.clone())),
            right: GateInput::Input(String::from("x")),
            symbol: GateType::Mul,
            label: String::from("g2"),
        };

        let g3 = Gate::<Fr> {
            left: GateInput::Input(String::from("x")),
            right: GateInput::Constant(F::from(2u64)),
            symbol: GateType::Mul,
            label: String::from("g3"),
        };

        let g4 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g2.clone())),
            right: GateInput::Gate(Box::new(g3.clone())),
            symbol: GateType::Add,
            label: String::from("g4"),
        };

        let g5 = Gate::<Fr> {
            left: GateInput::Gate(Box::new(g4.clone())),
            right: GateInput::Constant(F::from(5u64)),
            symbol: GateType::Add,
            label: String::from("g5"),
        };

        return vec![g1, g2, g3, g4, g5];
    }

    #[test]
    fn test_gate_input_eq() {
        let gates = sample_gates_0();
        let g1 = &gates[0];
        let g2 = &gates[1];
        let g3 = &gates[2];
        assert_eq!(g1.left, g1.right);
        assert_eq!(g1.left, g2.right);
        assert_eq!(g3.right, GateInput::Constant(F::from(2u64)));
    }

    #[test]
    fn test_into_gate() {
        let gates = sample_gates_0();
        let g1 = &gates[0];
        let g2 = &gates[1];
        assert_eq!(g2.left.into_gate().unwrap(), *g1);
        assert_eq!(g2.right.into_gate(), Err(Error::GateInputNotGate));
    }

    #[test]
    fn test_gates_to_matrices() {
        let gates = sample_gates_0();

        println!("Gates from sample_gates_0");
        for gate in gates.iter() {
            println!("{}", gate);
        }
        println!("");


        let mut left_input_map = BTreeMap::<GateInput<F>, usize>::new();
        let mut right_input_map = BTreeMap::<GateInput<F>, usize>::new();

        println!("Left inputs");
        let mut j = 0;
        for gate in gates.iter() {
            //left_input_map.entry(gate.left.clone()).or_insert(i);
            //right_input_map.entry(gate.right.clone()).or_insert(i);

            if !left_input_map.contains_key(&gate.left) {
                left_input_map.insert(gate.left.clone(), j);
                println!("{}, {}", j, gate.left);
                j += 1;
            }
        }
        println!("");

        j = 0;
        println!("Right inputs");
        for gate in gates.iter() {
            if !right_input_map.contains_key(&gate.right) {
                right_input_map.insert(gate.right.clone(), j);
                println!("{}, {}", j, gate.right);
                j += 1;
            }
        }
        let n_g = gates.len();
        let n_i = left_input_map.len() + right_input_map.len();
        let n_o = gates.len();

        assert_eq!(n_g, 5);
        assert_eq!(n_i, 8);
        assert_eq!(n_o, 5);

        let matrix_width = n_i + n_o;
        assert_eq!(matrix_width, 13);

        let mut matrix_a = empty_matrix::<F>(matrix_width);
        let mut matrix_b = empty_matrix::<F>(matrix_width);
        let mut matrix_c = empty_matrix::<F>(matrix_width);
 
        // TODO: check for off-by-one errors here! the paper uses matrices which start from 1, but
        // in Rust we start from 0.
        for (i, gate) in gates.iter().enumerate() {
            let row = n_i + i;
            let l_i = *left_input_map.get(&gate.left).unwrap() + 1;
            let r_i = *right_input_map.get(&gate.right).unwrap() + 1;

            let col_a = match gate.symbol {
                GateType::Add => 0,
                GateType::Mul => l_i,
            };

            // Matrix A
            matrix_a[row].push((F::from(1u64), col_a));
            
            // Matrix B
            if gate.symbol == GateType::Add {
                if l_i != r_i {
                    matrix_b[row].push((F::from(1u64), l_i));
                }
            }
            matrix_b[row].push((F::from(1u64), r_i));

            // Matrix C
            matrix_c[row].push((F::from(1u64), row));
        }

        println!("");
        print_matrix("matrix a", matrix_a.clone());
        print_matrix("matrix b", matrix_b.clone());
        print_matrix("matrix c", matrix_c.clone());

        let joint_arith = sparse_matrices_to_polys(
            matrix_a.clone(),
            matrix_b.clone(),
            matrix_c.clone(),
            100,
            0,
        );

        let joint_matrix = sum_matrices(
            &matrix_a,
            &matrix_b,
            &matrix_c,
        );
        print_joint_matrix("joint matrix", joint_matrix);
        print_polynomial("row", joint_arith.row.polynomial());
        print_polynomial("col", joint_arith.col.polynomial());
        print_polynomial("row_col", joint_arith.row_col.polynomial());
        print_polynomial("val_a", joint_arith.val_a.polynomial());
        print_polynomial("val_b", joint_arith.val_b.polynomial());
        print_polynomial("val_c", joint_arith.val_c.polynomial());
    }

    fn sparse_matrices_to_polys(
        a: Matrix<Fr>,
        b: Matrix<Fr>,
        c: Matrix<Fr>,
        num_constraints: usize,
        num_formatted_input_variables: usize,
    ) -> MatrixArithmetization<F> {
        let mut a = a.clone();
        let mut b = b.clone();
        let mut c = c.clone();

        let joint_matrix = sum_matrices(&a, &b, &c);
        let num_non_zero_val = num_non_zero(&joint_matrix);

        let domain_h = GeneralEvaluationDomain::<Fr>::new(num_constraints)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
        let domain_k = GeneralEvaluationDomain::<Fr>::new(num_non_zero_val)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
        let x_domain = GeneralEvaluationDomain::<Fr>::new(num_formatted_input_variables)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();

        let joint_arith = arithmetize_matrix(
            &joint_matrix,
            &mut a,
            &mut b,
            &mut c,
            domain_k,
            domain_h,
            x_domain,
        );

        joint_arith
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
