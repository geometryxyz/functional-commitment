#[cfg(test)]
mod tests {
    use std::cmp::max;

    use ac2tft::{printmatrix, sample_matrices, SparseMatrices};
    use ark_bls12_381::Bls12_381;
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, One, PrimeField, Zero};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
    use ark_poly_commit::PolynomialCommitment;
    use ark_poly_commit::{marlin_pc::MarlinKZG10, LabeledCommitment};
    use ark_std::rand::thread_rng;
    use ark_std::test_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::{marlin_kzg::KZG10, AdditivelyHomomorphicPCS};
    // use index_private_marlin::{fc_arith, AHPForR1CS, Marlin, fc_arith_withoud_reindexing};
    use new_ac_compiler::circuit::Circuit;
    use new_ac_compiler::constraint_builder::ConstraintBuilder;
    use new_ac_compiler::error::Error;
    use new_ac_compiler::gate::GateType;
    use new_ac_compiler::variable::VariableType;
    use new_marlin::Marlin;
    use proof_of_function_relation::t_diag::TDiag;
    use proof_of_function_relation::t_functional_triple::TFT;
    use proof_of_function_relation::t_strictly_lower_triangular_test::TStrictlyLowerTriangular;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    use new_ac_compiler::circuit_compiler::{CircuitCompiler, VanillaCompiler};

    use crate::{diag_test, slt_test};

    type F = Fr;
    type PC = KZG10<Bn254>;

    type MultiPC = MarlinKZG10<Bn254, DensePolynomial<Fr>>;
    type MarlinInst = Marlin<Fr, PC, FS>;

    fn circuit_test_template<Func>(constraints: Func)
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

        let (pk, vk, rands) = MarlinInst::index(&universal_srs, &index_info, a, b, c, rng).unwrap();

        // TEST MARLIN
        let proof = MarlinInst::prove(&pk, cb.assignment, rng).unwrap();

        assert!(MarlinInst::verify(
            &vk,
            &vec![F::one(), F::from(2u64), F::from(5u64), F::from(7u64)],
            proof,
            rng
        )
        .unwrap());

        // TEST PROOF OF FUNCTION
        let labels = vec![
            "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
        ];
        let commits: Vec<LabeledCommitment<_>> = vk
            .commits
            .iter()
            .zip(labels.iter())
            .map(|(cm, &label)| LabeledCommitment::new(label.into(), cm.clone(), Some(domain_k.size() + 1)))
            .collect();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());
        
        //JUST REVERSE ROW A AND COL A TO GET STRICTLY UPPER TRIANGULAR
        let tft_proof = TFT::<F, PC, FS>::prove(
            &pk.committer_key,
            index_info.number_of_input_rows,
            &domain_k,
            &domain_h,
            Some(domain_k.size() + 1), //enforced_degree_bound
            &pk.index.a_arith.col, // row_a_poly,
            &pk.index.a_arith.row, // col_a_poly,
            &commits[1], // row_a_commit,
            &commits[0],// col_a_commit,
            &rands[1],// row_a_random,
            &rands[0],// col_a_random,
            &pk.index.b_arith.col, // row_b_poly,
            &pk.index.b_arith.row, // col_b_poly,
            &commits[4], // row_b_commit,
            &commits[3], // col_b_commit,
            &rands[4], // row_b_random,
            &rands[3], // col_b_random,
            &pk.index.c_arith.row, // row_c_poly,
            &pk.index.c_arith.col, // col_c_poly,
            &pk.index.c_arith.val, // val_c_poly,
            &commits[6], // row_c_commit,
            &commits[7], // col_c_commit,
            &commits[8], // val_c_commit,
            &rands[6],// row_c_random,
            &rands[7],// col_c_random,
            &rands[8],// val_c_random,
            &mut fs_rng, // fs_rng,
            rng // rng,
        ).unwrap();

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
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", F::from(7u64))?;

            let x_square = cb.enforce_constraint(&x, &x, GateType::Mul, VariableType::Witness)?;
            let x_cube =
                cb.enforce_constraint(&x_square, &x, GateType::Mul, VariableType::Witness)?;

            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul, VariableType::Witness)?;
            let x_qubed_plus_2x =
                cb.enforce_constraint(&x_cube, &two_x, GateType::Add, VariableType::Witness)?;

            let _ = cb.enforce_constraint(
                &x_qubed_plus_2x,
                &five,
                GateType::Add,
                VariableType::Output,
            )?;

            Ok(())
        };

        circuit_test_template(constraints);
    }
}

#[macro_export]
/// Print a Matrix
macro_rules! slt_test {
    ($matrix:expr, $num_of_pub_inputs_plus_one:expr) => {
        for (row_index, row) in $matrix.iter().enumerate() {
            for (_, col_index) in row {
                assert!(row_index >= $num_of_pub_inputs_plus_one);
                assert!(row_index > *col_index);
            }
        }
    };
}

#[macro_export]
/// Print a Matrix
macro_rules! diag_test {
    ($matrix:expr) => {
        for (row_index, row) in $matrix.iter().enumerate() {
            for (_, col_index) in row {
                assert_eq!(row_index, *col_index);
            }
        }
    };
}
