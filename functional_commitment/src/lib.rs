#[cfg(test)]
mod tests {
    use std::cmp::max;

    use ac2tft::{printmatrix, sample_matrices, SparseMatrices};
    use ark_bls12_381::Bls12_381;
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, One, PrimeField, Zero};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
    use ark_poly_commit::{marlin_pc::MarlinKZG10, LabeledCommitment};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use ark_std::test_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::{marlin_kzg::KZG10, AdditivelyHomomorphicPCS};
    // use index_private_marlin::{fc_arith, AHPForR1CS, Marlin, fc_arith_withoud_reindexing};
    use new_marlin::{Marlin};
    use new_ac_compiler::circuit::Circuit;
    use new_ac_compiler::constraint_builder::ConstraintBuilder;
    use new_ac_compiler::error::Error;
    use new_ac_compiler::gate::GateType;
    use new_ac_compiler::variable::VariableType;
    use proof_of_function_relation::t_diag::TDiag;
    use proof_of_function_relation::t_strictly_lower_triangular_test::TStrictlyLowerTriangular;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    use new_ac_compiler::circuit_compiler::{CircuitCompiler, VanillaCompiler};

    use crate::{slt_test, diag_test};

    type F = Fr;
    type PC = KZG10<Bn254>;

    type MultiPC = MarlinKZG10<Bn254, DensePolynomial<Fr>>;
    type MarlinInst = Marlin<Fr, PC, FS>;

    #[test]
    fn test_airth() {
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", F::from(7u64))?;

            let x_square = cb.enforce_constraint(&x, &x, GateType::Mul, VariableType::Witness)?;
            let x_cube = cb.enforce_constraint(&x_square, &x, GateType::Mul, VariableType::Witness)?;

            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul, VariableType::Witness)?;
            let x_qubed_plus_2x = cb.enforce_constraint(&x_cube, &two_x, GateType::Add, VariableType::Witness)?;

            let _ = cb.enforce_constraint(&x_qubed_plus_2x, &five, GateType::Add, VariableType::Output)?;

            Ok(())
        };

        let mut cb = ConstraintBuilder::<F>::new();

        let synthesized_circuit = Circuit::synthesize(constraints, &mut cb).unwrap();
        let (index_info, a, b, c) = VanillaCompiler::<F>::ac2tft(&synthesized_circuit);

        assert_eq!(true, index_info.check_domains_sizes::<F>());

        let domain_k = GeneralEvaluationDomain::<F>::new(index_info.number_of_non_zero_entries).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(index_info.number_of_constraints).unwrap();

        // let enforced_degree_bound = domain_k.size() + 1;
        // let enforced_hiding_bound = 1;

        slt_test!(a, index_info.number_of_input_rows);
        slt_test!(b, index_info.number_of_input_rows);
        diag_test!(c);

        let rng = &mut test_rng();

        let universal_srs = MarlinInst::universal_setup(&index_info, rng).unwrap();

        let (pk, vk, rands) =
            MarlinInst::index(&universal_srs, &index_info, a, b, c).unwrap();

    /*

    BEGIN A TEST

    */
        let labels = vec![
            "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
        ];
        let commits: Vec<LabeledCommitment<_>> = vk.commits.iter().zip(labels.iter()).map(|(cm, &label)| LabeledCommitment::new(label.into(), cm.clone(), None)).collect();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let a_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &pk.committer_key,
            index_info.number_of_input_rows,
            &domain_k,
            &domain_h,
            &pk.index.a_arith.col,
            &commits[1],
            &rands[1],
            &pk.index.a_arith.row,
            &commits[0],
            &rands[0],
            None,
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            true,
            TStrictlyLowerTriangular::verify(
                &vk.verifier_key,
                &pk.committer_key,
                index_info.number_of_input_rows,
                &domain_k,
                &domain_h,
                &commits[1],
                &commits[0],
                None,
                a_proof,
                &mut fs_rng
            )
            .is_ok()
        );

    /*

    END A TEST

    */
    /*

    BEGIN B TEST

    */
        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let b_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &pk.committer_key,
            index_info.number_of_input_rows,
            &domain_k,
            &domain_h,
            &pk.index.b_arith.col,
            &commits[4],
            &rands[4],
            &pk.index.b_arith.row,
            &commits[3],
            &rands[3],
            None,
            &mut fs_rng,
            rng,
        )
        .unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            true,
            TStrictlyLowerTriangular::verify(
                &vk.verifier_key,
                &pk.committer_key,
                index_info.number_of_input_rows,
                &domain_k,
                &domain_h,
                &commits[4],
                &commits[3],
                None,
                b_proof,
                &mut fs_rng
            )
            .is_ok()
        );

    /*

    END B TEST

    */
    /*

    BEGIN C TEST

    */
        let c_proof = TDiag::<F, PC, FS>::prove(
            &pk.committer_key,
            index_info.number_of_input_rows,
            &pk.index.c_arith.row,
            &pk.index.c_arith.col,
            &pk.index.c_arith.val,
            &commits[6],
            &commits[7],
            &commits[8],
            &rands[6],
            &rands[7],
            &rands[8],
            None,
            &domain_k,
            &domain_h,
            index_info.number_of_constraints,
            rng
        ).unwrap();

        let is_valid = TDiag::<F, PC, FS>::verify(
            &vk.verifier_key,
            index_info.number_of_input_rows,
            &commits[6],
            &commits[7],
            &commits[8],
            None,
            &domain_h,
            &domain_k,
            index_info.number_of_constraints,
            c_proof,
        );

        assert!(is_valid.is_ok());

    /*

    END C TEST

    */
    }

    #[test]
    fn test_index_private_marin() {
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

        let mut cb = ConstraintBuilder::<F>::new();

        let synthesized_circuit = Circuit::synthesize(constraints, &mut cb).unwrap();
        let (index_info, a, b, c) = VanillaCompiler::<F>::ac2tft(&synthesized_circuit);

        let rng = &mut test_rng();

        let universal_srs = MarlinInst::universal_setup(&index_info, rng).unwrap();

        let (pk, vk, _) =
            MarlinInst::index(&universal_srs, &index_info, a, b, c).unwrap();

        let proof =
            MarlinInst::prove(&pk, cb.assignment, rng).unwrap();

        assert!(MarlinInst::verify(&vk, &vec![F::one(), F::from(2u64), F::from(5u64), F::from(7u64)], proof, rng).unwrap());
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
