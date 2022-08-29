#[cfg(test)]
mod tests {
    use std::cmp::max;

    use ac2tft::{printmatrix, sample_matrices, SparseMatrices};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{to_bytes, PrimeField, Zero};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial, Polynomial};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use ark_std::test_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::{marlin_kzg::KZG10, AdditivelyHomomorphicPCS};
    use index_private_marlin::{AHPForR1CS, fc_arith};
    use new_ac_compiler::circuit::Circuit;
    use new_ac_compiler::constraint_builder::ConstraintBuilder;
    use new_ac_compiler::gate::GateType;
    use new_ac_compiler::variable::VariableType;
    use new_ac_compiler::error::Error;
    use proof_of_function_relation::t_diag::TDiag;
    use proof_of_function_relation::t_strictly_lower_triangular_test::TStrictlyLowerTriangular;
    use rand_chacha::ChaChaRng;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;
    use new_ac_compiler::circuit_compiler::{VanillaCompiler, CircuitCompiler};

    type F = Fr;
    type PC = KZG10<Bn254>;

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
        let r1csf_index_from_synthesized = VanillaCompiler::<F>::ac2tft(&synthesized_circuit);

        // interpolation domain must be greater or equal to output domain (dl comparison constraint)
        let interpolation_domain_size = max(r1csf_index_from_synthesized.number_of_non_zero_entries, r1csf_index_from_synthesized.number_of_constraints);

        let domain_k = GeneralEvaluationDomain::<F>::new(interpolation_domain_size).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(r1csf_index_from_synthesized.number_of_constraints).unwrap();

        println!("num non zero: {}, num of constraints: {}", r1csf_index_from_synthesized.number_of_non_zero_entries, r1csf_index_from_synthesized.number_of_constraints);
        println!("interpolation: {}, output size: {}", domain_k.size(), domain_h.size());

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let rng = &mut thread_rng();

        let max_degree = 200;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        /*
        
        BEGIN A TEST

        */
        let a_arith = fc_arith(&r1csf_index_from_synthesized.a, domain_k, domain_h, "a", false);

        let polys = [
            a_arith.row.clone(), 
            a_arith.col.clone(), 
        ];

        let (commits, rands) = PC::commit(&ck, &polys, None).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let a_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            r1csf_index_from_synthesized.number_of_input_rows,
            &domain_k,
            &domain_h,
            &a_arith.col,
            &commits[1],
            &rands[1],
            &a_arith.row,
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
                &vk,
                &ck,
                r1csf_index_from_synthesized.number_of_input_rows,
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

        let b_arith = fc_arith(&r1csf_index_from_synthesized.b, domain_k, domain_h, "b", false);

        let polys = [
            b_arith.row.clone(), 
            b_arith.col.clone(), 
        ];

        let (commits, rands) = PC::commit(&ck, &polys, None).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let b_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            r1csf_index_from_synthesized.number_of_input_rows,
            &domain_k,
            &domain_h,
            &b_arith.col,
            &commits[1],
            &rands[1],
            &b_arith.row,
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
                &vk,
                &ck,
                r1csf_index_from_synthesized.number_of_input_rows,
                &domain_k,
                &domain_h,
                &commits[1],
                &commits[0],
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

        let c_arith = fc_arith(&r1csf_index_from_synthesized.c, domain_k, domain_h, "c", true);

        // for val in c_arith.evals_on_K.val.evals {
        //     println!("{}", val);
        // }

        let polys = [
            c_arith.row.clone(), 
            c_arith.col.clone(), 
            c_arith.val.clone(), 
        ];

        let (commits, rands) = PC::commit(&ck, &polys, None).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let c_proof = TDiag::<F, PC, FS>::prove(
            &ck,
            r1csf_index_from_synthesized.number_of_input_rows,
            &c_arith.row,
            &c_arith.col,
            &c_arith.val,
            &commits[0],
            &commits[1],
            &commits[2],
            &rands[0],
            &rands[1],
            &rands[2],
            None,
            &domain_k,
            &domain_h,
            r1csf_index_from_synthesized.number_of_constraints,
            rng
        ).unwrap();


        let is_valid = TDiag::<F, PC, FS>::verify(
            &vk,
            r1csf_index_from_synthesized.number_of_input_rows,
            &commits[0],
            &commits[1],
            &commits[2],
            None,
            &domain_h,
            &domain_k,
            r1csf_index_from_synthesized.number_of_constraints,
            c_proof,
        );

        assert!(is_valid.is_ok());

        /*
        
        END C TEST
        
        */
    }
}
