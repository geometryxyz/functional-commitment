#[cfg(test)]
mod tests {
    use ac_compiler::constraint_builder::ConstraintBuilder;
    use ac_compiler::error::Error;
    use ac_compiler::gate::GateType;
    use ac_compiler::variable::VariableType;
    use ac_compiler::{circuit::Circuit, variable::Variable};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::bytes::ToBytes;
    use ark_ff::PrimeField;
    use ark_ff::{to_bytes, Field, One};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly_commit::LabeledCommitment;
    use ark_std::test_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::marlin_kzg::KZG10;
    use index_private_marlin::Marlin;
    use proof_of_function_relation::t_functional_triple::TFT;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    use ac_compiler::circuit_compiler::{CircuitCompiler, VanillaCompiler};

    use crate::{diag_test, slt_test};

    type F = Fr;
    type PC = KZG10<Bn254>;

    type MarlinInst = Marlin<Fr, PC, FS>;

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
        let neg_ac =
            cb.enforce_constraint(&ac, &minus_one, GateType::Mul, VariableType::Witness)?;

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

        let mut c_inputs: Vec<Variable<F>> = vec![];
        for (i, val) in c.iter().enumerate() {
            let input = cb.new_input_variable(format!("c{}", i).as_str(), *val)?;
            c_inputs.push(input);
        }

        let mut r;
        let t = cb.enforce_constraint(&x, &k, GateType::Add, VariableType::Witness)?;
        r = enforce_x_7(cb, &t)?;
        for i in 1..n_rounds {
            let r_plus_k = cb.enforce_constraint(&r, &k, GateType::Add, VariableType::Witness)?;
            let t = cb.enforce_constraint(
                &r_plus_k,
                &c_inputs[i],
                GateType::Add,
                VariableType::Witness,
            )?;
            r = enforce_x_7(cb, &t)?;
        }

        let _ = cb.enforce_constraint(&r, &k, GateType::Add, VariableType::Output)?;

        Ok(())
    }

    fn circuit_test_template<Func>(constraints: Func, inputs: &Vec<F>, outputs: &Vec<F>)
    where
        Func: FnOnce(&mut ConstraintBuilder<F>) -> Result<(), Error>,
    {
        let mut cb = ConstraintBuilder::<F>::new();

        let synthesized_circuit = Circuit::synthesize(constraints, &mut cb).unwrap();
        let (index_info, a, b, c) = VanillaCompiler::<F>::ac2tft(&synthesized_circuit);

        assert_eq!(true, index_info.check_domains_sizes::<F>());

        let domain_k =
            GeneralEvaluationDomain::<F>::new(index_info.number_of_non_zero_entries).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(index_info.number_of_constraints).unwrap();

        slt_test!(a, index_info.number_of_input_rows);
        slt_test!(b, index_info.number_of_input_rows);
        diag_test!(c);

        let rng = &mut test_rng();

        let universal_srs = MarlinInst::universal_setup(&index_info, rng).unwrap();

        let (pk, vk) = MarlinInst::index(&universal_srs, &index_info, a, b, c, rng).unwrap();

        // TEST MARLIN
        let proof = MarlinInst::prove(&pk, cb.assignment, rng).unwrap();

        assert!(MarlinInst::verify(&vk, inputs, outputs, proof, rng, &pk.committer_key).unwrap());

        // TEST PROOF OF FUNCTION
        let labels = vec![
            "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
        ];
        let commits: Vec<LabeledCommitment<_>> = vk
            .commits
            .iter()
            .zip(labels.iter())
            .map(|(cm, &label)| {
                LabeledCommitment::new(label.into(), cm.clone(), Some(domain_k.size() + 1))
            })
            .collect();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        //JUST REVERSE ROW A AND COL A TO GET STRICTLY UPPER TRIANGULAR
        let tft_proof = TFT::<F, PC, FS>::prove(
            &pk.committer_key,
            index_info.number_of_input_rows,
            &domain_k,
            &domain_h,
            Some(domain_k.size() + 1), //enforced_degree_bound
            &pk.index.a_arith.col,     // row_a_poly,
            &pk.index.a_arith.row,     // col_a_poly,
            &commits[1],               // row_a_commit,
            &commits[0],               // col_a_commit,
            &pk.rands[1],              // row_a_random,
            &pk.rands[0],              // col_a_random,
            &pk.index.b_arith.col,     // row_b_poly,
            &pk.index.b_arith.row,     // col_b_poly,
            &commits[4],               // row_b_commit,
            &commits[3],               // col_b_commit,
            &pk.rands[4],              // row_b_random,
            &pk.rands[3],              // col_b_random,
            &pk.index.c_arith.row,     // row_c_poly,
            &pk.index.c_arith.col,     // col_c_poly,
            &pk.index.c_arith.val,     // val_c_poly,
            &commits[6],               // row_c_commit,
            &commits[7],               // col_c_commit,
            &commits[8],               // val_c_commit,
            &pk.rands[6],              // row_c_random,
            &pk.rands[7],              // col_c_random,
            &pk.rands[8],              // val_c_random,
            &mut fs_rng,               // fs_rng,
            rng,                       // rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let is_valid = TFT::<F, PC, FS>::verify(
            &vk.verifier_key,
            &pk.committer_key,
            index_info.number_of_input_rows,
            &commits[1],
            &commits[0],
            &commits[4],
            &commits[3],
            &commits[6],
            &commits[7],
            &commits[8],
            Some(domain_k.size() + 1),
            &domain_h,
            &domain_k,
            tft_proof,
            &mut fs_rng,
        );

        assert!(is_valid.is_ok());
    }

    #[test]
    fn test_simple_circuit() {
        // Tests a circuit which encodes the equation x^3 + 2x + 5

        let x_val = F::from(7u64);

        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", x_val)?;

            let x_square = cb.enforce_constraint(&x, &x, GateType::Mul, VariableType::Witness)?;
            let x_cube =
                cb.enforce_constraint(&x_square, &x, GateType::Mul, VariableType::Witness)?;

            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul, VariableType::Witness)?;
            let x_qubed_plus_2x =
                cb.enforce_constraint(&x_cube, &two_x, GateType::Add, VariableType::Witness)?;

            // output = dec: 362, hex: 16A
            let _ = cb.enforce_constraint(
                &x_qubed_plus_2x,
                &five,
                GateType::Add,
                VariableType::Output,
            )?;

            Ok(())
        };

        let inputs = vec![F::one(), F::from(2u64), F::from(5u64), x_val];
        let outputs = vec![F::from(362u64)];
        circuit_test_template(constraints, &inputs, &outputs);
    }

    #[test]
    fn test_mux1() {
        // Inputs: a, b, and c
        // Output (b - a) * c + a
        // If c = 0, outputs a
        // If c = 1, outputs b

        let a_val = F::from(123);
        let b_val = F::from(456);
        let c_val = F::from(1);
        let expected_output = b_val;

        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            build_mux1_circuit::<Fr>(cb, a_val, b_val, c_val)?;
            Ok(())
        };

        let inputs = vec![F::one(), a_val, b_val, c_val, F::from(-1)];
        let outputs = vec![expected_output];
        circuit_test_template(constraints, &inputs, &outputs);

        let c_val = F::from(0);
        let expected_output = a_val;
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            build_mux1_circuit::<Fr>(cb, a_val, b_val, c_val)?;
            Ok(())
        };

        let inputs = vec![F::one(), a_val, b_val, c_val, F::from(-1)];
        let outputs = vec![expected_output];
        circuit_test_template(constraints, &inputs, &outputs);
    }

    #[test]
    fn test_composed_circuit() {
        let x_val = F::from(2u64);
        let expected_output = F::from(16u64);
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            build_x4_circuit::<Fr>(cb, x_val)?;
            Ok(())
        };
        let inputs = vec![F::one(), F::one(), x_val];
        let outputs = vec![expected_output];
        circuit_test_template(constraints, &inputs, &outputs);
    }

    fn registers_to_f(registers: &[u64; 4]) -> F {
        let mut writer = vec![];
        let _ = registers.write(&mut writer);
        F::from_le_bytes_mod_order(&writer)
    }

    fn build_cts() -> Vec<F> {
        let cts_bigints = vec![
            [0, 0, 0, 0],
            [
                3366560258570492133,
                10070564347787222493,
                15604622621500992217,
                3327803545572961123,
            ],
            [
                2961805725237629769,
                17895266365305129804,
                3244298544786782209,
                2431874893373010386,
            ],
        ];
        let mut cts = vec![];
        for c in cts_bigints.iter() {
            cts.push(registers_to_f(c));
        }
        cts
    }

    fn pow7(x: F) -> F {
        x * x * x * x * x * x * x
    }

    #[test]
    fn test_mimc7() {
        // The mimc7 circuit where nRounds = 2 and k = 2
        let x_val = F::from(2u64);
        let k = F::from(2u64);
        let cts = build_cts();
        let mut inputs = vec![F::one(), x_val, k];
        for c in cts.iter() {
            inputs.push(*c);
        }
        let n_rounds = 2;

        let t = x_val + k;

        let mut r;
        r = pow7(t);

        for i in 1..n_rounds {
            let t = r + k + cts[i];
            r = pow7(t);
        }

        let result = r + k;
        let outputs = vec![result];

        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            build_mimc7_circuit(cb, x_val, cts)?;
            Ok(())
        };
        circuit_test_template(constraints, &inputs, &outputs);
    }
}
