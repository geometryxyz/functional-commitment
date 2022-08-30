#[cfg(test)]
mod tests {
    use crate::tests::naive_matrix_encoding;
    use crate::variable::VariableType;
    use crate::{
        circuit::Circuit,
        circuit_compiler::{CircuitCompiler, VanillaCompiler},
        constraint_builder::ConstraintBuilder,
        error::Error,
        example_circuits::sample_circuit_2,
        gate::GateType,
        printmatrix, R1CSfIndex,
        slt_test, diag_test
    };
    use ark_bn254::Bn254;
    use ark_bn254::Fr;
    use ark_ff::{to_bytes, Zero};
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
    fn test_matrix_correctness() {
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

        // Perform matrix multiplications
        let inner_prod_fn = |row: &[(F, usize)]| {
            let mut acc = F::zero();
            for &(_, i) in row {
                acc += cb.assignment[i];
            }
            acc
        };

        let z_a: Vec<F> = r1csf_index_from_synthesized.a.iter().map(|row| inner_prod_fn(row)).collect();
        let z_b: Vec<F> = r1csf_index_from_synthesized.b.iter().map(|row| inner_prod_fn(row)).collect();
        let z_c: Vec<F> = r1csf_index_from_synthesized.c.iter().map(|row| inner_prod_fn(row)).collect();

        assert_eq!(z_a.len(), z_b.len());
        assert_eq!(z_b.len(), z_c.len());
        for ((&za_i, &zb_i), &zc_i) in z_a.iter().zip(z_b.iter()).zip(z_c.iter()) {
            assert_eq!(za_i * zb_i, zc_i);
        }

        let formatted_input_assignment = cb.assignment[..r1csf_index_from_synthesized.number_of_input_rows].to_vec();
        let witness_assignment = cb.assignment[r1csf_index_from_synthesized.number_of_input_rows..].to_vec();
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
