use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::LabeledPolynomial;

#[inline]
pub fn powers_of<F>(scalar: F) -> impl Iterator<Item = F>
where
    F: Field,
{
    core::iter::successors(Some(F::one()), move |p| Some(*p * scalar))
}

/// Generates the concatenation of geometric sequences that all share a common ratio
pub fn generate_sequence<F: PrimeField>(
    common_ratio: F,
    initial_terms: &[F],
    subsequence_lengths: &[usize],
) -> Vec<F> {
    assert_ne!(initial_terms.len(), 0);
    assert_eq!(subsequence_lengths.len(), initial_terms.len());

    let mut concatenation = Vec::<F>::default();

    for (i, term) in initial_terms.iter().enumerate() {
        for ratio_power in powers_of(common_ratio).take(subsequence_lengths[i]) {
            let val = *term * ratio_power;
            concatenation.push(val);
        }
    }

    concatenation
}

pub fn gen_t_diag_test_polys<F: FftField>(
    domain_k: GeneralEvaluationDomain<F>,
    domain_h: GeneralEvaluationDomain<F>,
    degree_bound: Option<usize>,
    hiding_bound: Option<usize>,
) -> Vec<LabeledPolynomial<F, DensePolynomial<F>>> {
    // M indices
    /*
    00, 01, 02, 03
    10, 11, 12, 13
    20, 21, 22, 23
    30, 31, 32, 33
    */

    // A (slt), t = 2
    /*
    0, 0, 0, 0
    1, 0, 0, 0
    1, 2, 0, 0
    0, 3, 0, 0
    */

    // rowM and colM are vectors that encode position of each non-zero element

    // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
    // rowA_evals   =  w^1    w^2       w^2       w^3      w^3      w^3      w^3     w^3
    // colA_evals   =  w^0    w^0       w^1       w^1      w^1      w^1      w^1     w^1

    // B (stl), t = 1
    /*
    0, 0, 0, 0
    7, 0, 0, 0
    0, 0, 0, 0
    2, 2, 0, 0
    */

    // rowM and colM are vectors that encode position of each non-zero element

    // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
    // rowB_evals   =  w^1    w^3       w^3       w^3      w^3      w^3      w^3     w^3
    // colB_evals   =  w^0    w^0       w^1       w^1      w^1      w^1      w^1     w^1

    // C (diag)
    /*
    0, 0, 0, 0
    0, 0, 0, 0
    0, 0, 2, 0
    0, 0, 0, 2
    */

    // rowM and colM are vectors that encode position of each non-zero element

    // domain       =  1,    gamma,   gamma^2   gamma^3  gamma^4  gamma^5  gamma^6  gamma^7
    // rowC_evals   =  w^2    w^3       w^0       w^0      w^0      w^0      w^0     w^0
    // colC_evals   =  w^2    w^3       w^0       w^0      w^0      w^0      w^0     w^0
    // vaLC_evals   =   2      2         0         0        0        0        0       0

    let omega_0 = domain_h.element(0);
    let omega_1 = domain_h.element(1);
    let omega_2 = domain_h.element(2);
    let omega_3 = domain_h.element(3);

    let row_a_evals = vec![
        omega_1, omega_2, omega_2, omega_3, omega_3, omega_3, omega_3, omega_3,
    ];
    let col_a_evals = vec![
        omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1, omega_1,
    ];

    let row_b_evals = vec![
        omega_1, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3, omega_3,
    ];
    let col_b_evals = vec![
        omega_0, omega_0, omega_1, omega_1, omega_1, omega_1, omega_1, omega_1,
    ];

    let row_c_evals = vec![
        omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
    ];
    let col_c_evals = vec![
        omega_2, omega_3, omega_0, omega_0, omega_0, omega_0, omega_0, omega_0,
    ];

    let val_c_evals = vec![
        F::from(2u64),
        F::from(2u64),
        F::from(0u64),
        F::from(0u64),
        F::from(0u64),
        F::from(0u64),
        F::from(0u64),
        F::from(0u64),
    ];

    let row_a_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_a_evals));
    let col_a_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_a_evals));
    let row_a_poly = LabeledPolynomial::new(
        String::from("row_a_poly"),
        row_a_poly,
        degree_bound,
        hiding_bound,
    );
    let col_a_poly = LabeledPolynomial::new(
        String::from("col_a_poly"),
        col_a_poly,
        degree_bound,
        hiding_bound,
    );

    let row_b_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_b_evals));
    let col_b_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_b_evals));
    let row_b_poly = LabeledPolynomial::new(
        String::from("row_b_poly"),
        row_b_poly,
        degree_bound,
        hiding_bound,
    );
    let col_b_poly = LabeledPolynomial::new(
        String::from("col_b_poly"),
        col_b_poly,
        degree_bound,
        hiding_bound,
    );

    let row_c_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&row_c_evals));
    let col_c_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&col_c_evals));
    let val_c_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&val_c_evals));
    let row_c_poly = LabeledPolynomial::new(
        String::from("row_c_poly"),
        row_c_poly,
        degree_bound,
        hiding_bound,
    );
    let col_c_poly = LabeledPolynomial::new(
        String::from("col_c_poly"),
        col_c_poly,
        degree_bound,
        hiding_bound,
    );
    let val_c_poly = LabeledPolynomial::new(
        String::from("val_c_poly"),
        val_c_poly,
        degree_bound,
        hiding_bound,
    );

    vec![
        row_a_poly, col_a_poly, row_b_poly, col_b_poly, row_c_poly, col_c_poly, val_c_poly,
    ]
}
