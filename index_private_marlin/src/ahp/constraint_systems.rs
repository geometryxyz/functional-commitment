use std::collections::BTreeMap;

use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations as EvaluationsOnDomain,
    GeneralEvaluationDomain,
};
use ark_relations::r1cs::Matrix;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::{
    cfg_iter_mut,
    io::{Read, Write},
};
use derivative::Derivative;

use super::UnnormalizedBivariateLagrangePoly;

pub type LabeledPolynomial<F> = ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>;

/// Evaluations of various polynomials related to the constraint matrices,
/// over the same domain.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Debug(bound = "F: PrimeField"), Clone(bound = "F: PrimeField"))]
pub struct MatrixEvals<F: PrimeField> {
    /// Evaluations of the `row` polynomial.
    pub row: EvaluationsOnDomain<F>,
    /// Evaluations of the `col` polynomial.
    pub col: EvaluationsOnDomain<F>,
    /// Evaluations of the `val_a` polynomial.
    pub val: EvaluationsOnDomain<F>,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Debug(bound = "F: PrimeField"), Clone(bound = "F: PrimeField"))]
pub struct MatrixArithmetization<F: PrimeField> {
    /// LDE of the row indices of M^*.
    pub row: LabeledPolynomial<F>,
    /// LDE of the column indices of M^*.
    pub col: LabeledPolynomial<F>,
    /// LDE of the non-zero entries of A^*.
    pub val: LabeledPolynomial<F>,

    /// Evaluation of `self.row`, `self.col`, and `self.val` on the domain `K`.
    pub evals_on_k: MatrixEvals<F>,
}

pub fn arithmetize_matrix<F: PrimeField>(
    m: &Matrix<F>,
    interpolation_domain: GeneralEvaluationDomain<F>,
    output_domain: GeneralEvaluationDomain<F>,
    // input_domain: GeneralEvaluationDomain<F>,
    m_prefix: &str,
    is_diag_matrix: bool,
    degree_bound: Option<usize>,
    hiding_bound: Option<usize>,
) -> MatrixArithmetization<F> {
    let elems: Vec<_> = output_domain.elements().collect();

    let eq_poly_vals: BTreeMap<F, F> = output_domain
        .elements()
        .zip(output_domain.batch_eval_unnormalized_bivariate_lagrange_poly_with_same_inputs())
        .collect();

    let mut row_vec = Vec::with_capacity(interpolation_domain.size());
    let mut col_vec = Vec::with_capacity(interpolation_domain.size());
    let mut val_vec = Vec::with_capacity(interpolation_domain.size());
    let mut inverses = Vec::with_capacity(interpolation_domain.size());
    let mut count = 0;

    for (row_index, row) in m.into_iter().enumerate() {
        for (val, col_index) in row {
            let row_val = elems[row_index];
            let col_val = elems[*col_index];

            // We are dealing with the transpose of M
            row_vec.push(col_val);
            col_vec.push(row_val);

            val_vec.push(*val);
            inverses.push(eq_poly_vals[&col_val]);

            count += 1;
        }
    }
    ark_ff::batch_inversion::<F>(&mut inverses);
    drop(eq_poly_vals);

    cfg_iter_mut!(val_vec).zip(inverses).for_each(|(v, inv)| {
        *v *= &inv;
    });

    let last_column = *col_vec.last().unwrap();
    let last_row = *row_vec.last().unwrap();

    let val_to_append_col = if is_diag_matrix {
        F::one()
    } else {
        last_column
    };

    let val_to_append_row = if is_diag_matrix { F::one() } else { last_row };

    for _ in count..interpolation_domain.size() {
        col_vec.push(val_to_append_col);
        row_vec.push(val_to_append_row);
        val_vec.push(F::zero());
    }

    let row_evals_on_k = EvaluationsOnDomain::from_vec_and_domain(row_vec, interpolation_domain);
    let col_evals_on_k = EvaluationsOnDomain::from_vec_and_domain(col_vec, interpolation_domain);
    let val_evals_on_k = EvaluationsOnDomain::from_vec_and_domain(val_vec, interpolation_domain);

    let row = row_evals_on_k.clone().interpolate();
    let col = col_evals_on_k.clone().interpolate();
    let val = val_evals_on_k.clone().interpolate();

    let evals_on_k = MatrixEvals {
        row: row_evals_on_k,
        col: col_evals_on_k,
        val: val_evals_on_k,
    };

    MatrixArithmetization {
        row: LabeledPolynomial::new(
            format!("{}_row", m_prefix).into(),
            row,
            degree_bound,
            hiding_bound,
        ),
        col: LabeledPolynomial::new(
            format!("{}_col", m_prefix).into(),
            col,
            degree_bound,
            hiding_bound,
        ),
        val: LabeledPolynomial::new(
            format!("{}_val", m_prefix).into(),
            val,
            degree_bound,
            hiding_bound,
        ),
        evals_on_k,
    }
}
