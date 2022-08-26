#[cfg(test)]
mod tests {
    use crate::{example_circuits::sample_circuit_2, printmatrix, vanilla_ac2tft, R1CSfIndex};
    use ark_bn254::Fr;
    use ark_ff::FftField;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
    };
    use ark_poly_commit::LabeledPolynomial;
    use ark_relations::r1cs::Matrix;
    use proof_of_function_relation::t_functional_triple::TFT;

    type F = Fr;

    fn naive_matrix_encoding<F: FftField, D: EvaluationDomain<F>>(
        index: R1CSfIndex<F>,
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

            (row_poly, col_poly, val_poly)
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

            (row_poly, col_poly, val_poly)
        };

        let a_polys = encode_and_pad_matrix(index.a);
        let b_polys = encode_and_pad_matrix(index.b);
        let c_polys = encode_and_pad_c(index.c);

        let mut out = Vec::new();

        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            a_polys.0,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_col"),
            a_polys.1,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            a_polys.2,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            b_polys.0,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            b_polys.1,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            b_polys.2,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            c_polys.0,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            c_polys.1,
            None,
            None,
        ));
        out.push(LabeledPolynomial::new(
            String::from("a_row"),
            c_polys.2,
            None,
            None,
        ));

        out
    }

    #[test]
    fn manual_check_simple_circuit() {
        let circuit = sample_circuit_2();

        let index: R1CSfIndex<F> = vanilla_ac2tft(circuit);

        index.iter_matrices().for_each(|matrix| {
            printmatrix!(matrix);
            println!("")
        })
    }

    #[test]
    fn test_arithmetization() {
        let circuit = sample_circuit_2();
        let index: R1CSfIndex<F> = vanilla_ac2tft(circuit);

        let domain_k = GeneralEvaluationDomain::new(index.number_of_non_zero_entries).unwrap();
        let domain_h = GeneralEvaluationDomain::new(index.number_of_constraints).unwrap();

        let concrete_oracles = naive_matrix_encoding(index, domain_k, domain_h);

        
    }
}
