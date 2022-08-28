#[cfg(test)]
mod tests {
    use crate::tests::naive_matrix_encoding;
    use crate::{
        circuit::Circuit,
        circuit_compiler::{CircuitCompiler, VanillaCompiler},
        constraint_builder::ConstraintBuilder,
        error::Error,
        example_circuits::sample_circuit_2,
        gate::GateType,
        printmatrix, R1CSfIndex,
    };
    use ark_bn254::Bn254;
    use ark_bn254::Fr;
    use ark_ff::to_bytes;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use blake2::Blake2s;
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use homomorphic_poly_commit::{marlin_kzg::KZG10, AdditivelyHomomorphicPCS};
    use proof_of_function_relation::t_functional_triple::TFT;
    use rand_chacha::ChaChaRng;

    type F = Fr;
    type PC = KZG10<Bn254>;
    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;

    #[test]
    fn test_sample_circuit_2() {
        let constraints = |cb: &mut ConstraintBuilder<F>| -> Result<(), Error> {
            let two = cb.new_input_variable("two", F::from(2u64))?;
            let five = cb.new_input_variable("five", F::from(5u64))?;
            let x = cb.new_input_variable("x", F::from(7u64))?;

            let x_square = cb.enforce_constraint(&x, &x, GateType::Mul)?;
            let x_cube = cb.enforce_constraint(&x_square, &x, GateType::Mul)?;

            let two_x = cb.enforce_constraint(&two, &x, GateType::Mul)?;
            let x_qubed_plus_2x = cb.enforce_constraint(&x_cube, &two_x, GateType::Add)?;

            let _ = cb.enforce_constraint(&x_qubed_plus_2x, &five, GateType::Add)?;

            Ok(())
        };

        let mut cb = ConstraintBuilder::<F>::new(1);

        let synthesized_circuit = Circuit::synthesize(constraints, &mut cb).unwrap();
        let r1csf_index_from_synthesized = VanillaCompiler::<F>::ac2tft(&synthesized_circuit);
        // printmatrix!(r1csf_index_from_synthesized.a);
        // println!("====================================");
        // printmatrix!(r1csf_index_from_synthesized.b);
        // println!("====================================");
        // printmatrix!(r1csf_index_from_synthesized.c);

        // println!("====================================");
        // println!("====================================");
        // println!("====================================");

        let manual_circuit = sample_circuit_2();
        let r1csf_index_from_manual = VanillaCompiler::<F>::ac2tft(&manual_circuit);
        // printmatrix!(r1csf_index_from_manual.a);
        // println!("====================================");
        // printmatrix!(r1csf_index_from_manual.b);
        // println!("====================================");
        // printmatrix!(r1csf_index_from_manual.c);

        assert_eq!(r1csf_index_from_manual, r1csf_index_from_synthesized)
    }

    #[test]
    fn test_arithmetization() {
        let rng = &mut thread_rng();

        let circuit = sample_circuit_2();
        let index: R1CSfIndex<F> = VanillaCompiler::ac2tft(&circuit);

        let domain_k = GeneralEvaluationDomain::new(index.number_of_non_zero_entries).unwrap();
        let domain_h = GeneralEvaluationDomain::new(index.number_of_constraints).unwrap();

        let concrete_oracles = naive_matrix_encoding(&index, domain_k, domain_h);

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let (commitments, rands) = PC::commit(&ck, &concrete_oracles, Some(rng)).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        println!(
            "domain_k {} - non-zero entries {},\ndomain_h {} - matrix length {},\nt {}",
            domain_k.size(),
            index.number_of_non_zero_entries,
            domain_h.size(),
            index.number_of_constraints,
            index.number_of_input_rows
        );

        let proof = TFT::<F, PC, FS>::prove(
            &ck,
            index.number_of_input_rows,
            &domain_k,
            &domain_h,
            None,
            &concrete_oracles[0], // a_row_poly
            &concrete_oracles[1], // a_col_poly
            &commitments[0],      // a_row_commit
            &commitments[1],      // a_col_commit
            &rands[0],            // a_row_rand
            &rands[1],            // a_col_rand
            &concrete_oracles[3], // b_row_poly
            &concrete_oracles[4], // b_col_poly
            &commitments[3],      // b_row_commit
            &commitments[4],      // b_col_commit
            &rands[3],            // b_row_rand
            &rands[4],            // b_col_rand
            &concrete_oracles[6], // c_row_poly
            &concrete_oracles[7], // c_row_poly
            &concrete_oracles[8], // c_row_poly
            &commitments[6],      // c_row_commit
            &commitments[7],      // c_row_commit
            &commitments[8],      // c_row_commit
            &rands[6],            // c_row_rand
            &rands[7],            // c_row_rand
            &rands[8],            // c_row_rand
            &mut fs_rng,
            rng,
        )
        .unwrap();

        // let proof = TFT::<F, PC, D>::prove(
        //     &ck,
        //     t,
        //     &domain_k,
        //     &domain_h,
        //     Some(enforced_degree_bound),
        //     //a
        //     &row_a_poly,
        //     &col_a_poly,
        //     &a_commitments[0],
        //     &a_commitments[1],
        //     &a_rands[0],
        //     &a_rands[1],
        //     //b
        //     &row_b_poly,
        //     &col_b_poly,
        //     &b_commitments[0],
        //     &b_commitments[1],
        //     &b_rands[0],
        //     &b_rands[1],
        //     //c
        //     &row_c_poly,
        //     &col_c_poly,
        //     &val_c_poly,
        //     &c_commitments[0],
        //     &c_commitments[1],
        //     &c_commitments[2],
        //     &c_rands[0],
        //     &c_rands[1],
        //     &c_rands[2],
        //     &mut fs_rng,
        //     rng,
        // )
        // .unwrap();

        // let is_valid = TFT::<F, PC, D>::verify(
        //     &vk,
        //     &ck,
        //     t,
        //     &a_commitments[0],
        //     &a_commitments[1],
        //     &b_commitments[0],
        //     &b_commitments[1],
        //     &c_commitments[0],
        //     &c_commitments[1],
        //     &c_commitments[2],
        //     Some(enforced_degree_bound),
        //     &domain_h,
        //     &domain_k,
        //     proof,
        //     &mut fs_rng,
        // );

        // assert!(is_valid.is_ok());

        // TODO: commit to these polynomials and try the proof of function relation!
    }
}

use crate::R1CSfIndex;
use ark_ff::FftField;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, UVPolynomial};
use ark_poly_commit::LabeledPolynomial;
use ark_relations::r1cs::Matrix;

/// A naive encoding of the R1CSfIndex into polynomial. For testing with proof of function relation
pub fn naive_matrix_encoding<F: FftField, D: EvaluationDomain<F>>(
    index: &R1CSfIndex<F>,
    domain_k: D,
    domain_h: D,
) -> Vec<LabeledPolynomial<F, DensePolynomial<F>>> {
    let encode_and_pad_matrix = |matrix: Matrix<F>| {
        let omega = domain_h.element(1);
        let zero = F::zero();

        let mut row_evals = vec![F::zero(); domain_k.size()];
        let mut col_evals = vec![F::zero(); domain_k.size()];
        let mut val_evals = vec![F::zero(); domain_k.size()];

        let mut non_zero_entries = Vec::new();
        for (row_index, row) in matrix.iter().enumerate() {
            for (value, column_index) in row {
                non_zero_entries.push((row_index, column_index, value));
            }
        }

        if non_zero_entries.len() < domain_k.size() {
            let padding_length = domain_k.size() - non_zero_entries.len();
            let last_row = non_zero_entries.last().unwrap().0;
            let last_col = non_zero_entries.last().unwrap().1;
            for _i in 0..padding_length {
                non_zero_entries.push((last_row, last_col, &zero));
            }
        }

        for (i, (r_i, c_i, v_i)) in non_zero_entries.iter().enumerate() {
            row_evals[i] = omega.pow(&[*r_i as u64]);
            col_evals[i] = omega.pow(&[**c_i as u64]);
            val_evals[i] = **v_i;
        }

        let row_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&row_evals));
        let col_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&col_evals));
        let val_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&val_evals));

        vec![row_poly, col_poly, val_poly]
    };

    let encode_and_pad_c = |matrix: Matrix<F>| {
        let omega = domain_h.element(1);
        let zero = F::zero();

        let mut row_evals = vec![F::zero(); domain_k.size()];
        let mut col_evals = vec![F::zero(); domain_k.size()];
        let mut val_evals = vec![F::zero(); domain_k.size()];

        let mut non_zero_entries = Vec::new();
        for (row_index, row) in matrix.iter().enumerate() {
            for (value, column_index) in row {
                non_zero_entries.push((row_index, column_index, value));
            }
        }

        if non_zero_entries.len() < domain_k.size() {
            let padding_length = domain_k.size() - non_zero_entries.len();
            for _i in 0..padding_length {
                non_zero_entries.push((1, &1, &zero));
            }
        }

        for (i, (r_i, c_i, v_i)) in non_zero_entries.iter().enumerate() {
            row_evals[i] = omega.pow(&[*r_i as u64]);
            col_evals[i] = omega.pow(&[**c_i as u64]);
            val_evals[i] = **v_i;
        }

        let row_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&row_evals));
        let col_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&col_evals));
        let val_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&val_evals));

        vec![row_poly, col_poly, val_poly]
    };

    let a_polys = encode_and_pad_matrix(index.a.clone());
    let b_polys = encode_and_pad_matrix(index.b.clone());
    let c_polys = encode_and_pad_c(index.c.clone());

    let labels = vec![
        "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
    ];
    let labels = labels.iter().map(|&str| String::from(str));

    let all_labeled_polys: Vec<LabeledPolynomial<F, DensePolynomial<F>>> = labels
        .zip(a_polys.iter().chain(b_polys.iter()).chain(c_polys.iter()))
        .map(|(label, poly)| LabeledPolynomial::new(label, poly.clone(), None, None))
        .collect();

    all_labeled_polys
}
