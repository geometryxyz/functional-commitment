use ark_bn254::Fr;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
    UVPolynomial,
};
use ark_poly_commit::LabeledPolynomial;
use ark_std::UniformRand;
use rand::Rng;

#[inline]
pub fn powers_of<F>(scalar: F) -> impl Iterator<Item = F>
where
    F: Field,
{
    core::iter::successors(Some(F::one()), move |p| Some(*p * scalar))
}

#[macro_export]
macro_rules! to_poly {
    ($value:expr) => {
        DensePolynomial::from_coefficients_slice(&[$value])
    };
}

/// Lazy labelling of polynomials for testing.
#[macro_export]
macro_rules! label_polynomial {
    ($poly:expr) => {
        ark_poly_commit::LabeledPolynomial::new(
            stringify!($poly).to_owned(),
            $poly.clone(),
            None,
            None,
        )
    };
}

/// Extract the labels from a list of labeled elements
#[macro_export]
macro_rules! get_labels {
    ($vector:expr) => {
        $vector
            .iter()
            .map(|f| f.label().clone())
            .collect::<Vec<_>>()
    };
}

pub fn shift_dense_poly<F: Field>(
    p: &DensePolynomial<F>,
    shifting_factor: &F,
) -> DensePolynomial<F> {
    if *shifting_factor == F::one() {
        return p.clone();
    }

    let mut coeffs = p.coeffs().to_vec();
    let mut acc = F::one();
    for i in 0..coeffs.len() {
        coeffs[i] = coeffs[i] * acc;
        acc *= shifting_factor;
    }

    DensePolynomial::from_coefficients_vec(coeffs)
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

#[allow(dead_code)]
/// Sample a vector of random elements of type T
pub fn sample_vector<T: UniformRand, R: Rng>(seed: &mut R, length: usize) -> Vec<T> {
    (0..length)
        .collect::<Vec<usize>>()
        .iter()
        .map(|_| T::rand(seed))
        .collect::<Vec<_>>()
}

/// Sample a random polynomial of a defined degree
pub fn random_deg_n_polynomial<F: Field, R: Rng>(
    degree: usize,
    seed: &mut R,
) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(sample_vector(seed, degree))
}

pub fn compute_vanishing_poly_over_coset<F: FftField>(
    domain: &GeneralEvaluationDomain<F>, // domain to evaluate over
    poly_degree: u64,                    // degree of the vanishing polynomial
) -> Evaluations<F> {
    assert!(
        (domain.size() as u64) > poly_degree,
        "domain_size = {}, poly_degree = {}",
        domain.size() as u64,
        poly_degree
    );

    let group_gen = domain.element(1);
    let coset_gen = F::multiplicative_generator().pow(&[poly_degree, 0, 0, 0]);
    let v_h: Vec<_> = (0..domain.size())
        .map(|i| (coset_gen * group_gen.pow(&[poly_degree * i as u64, 0, 0, 0])) - F::one())
        .collect();
    Evaluations::from_vec_and_domain(v_h, *domain)
}

#[allow(dead_code)]
pub fn gen_t_diag_test_polys(
    domain_k: GeneralEvaluationDomain<Fr>,
    domain_h: GeneralEvaluationDomain<Fr>,
) -> Vec<LabeledPolynomial<Fr, DensePolynomial<Fr>>> {
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
        Fr::from(2u64),
        Fr::from(2u64),
        Fr::from(0u64),
        Fr::from(0u64),
        Fr::from(0u64),
        Fr::from(0u64),
        Fr::from(0u64),
        Fr::from(0u64),
    ];

    let row_a_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&row_a_evals));
    let col_a_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&col_a_evals));
    let row_a_poly = label_polynomial!(row_a_poly);
    let col_a_poly = label_polynomial!(col_a_poly);

    let row_b_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&row_b_evals));
    let col_b_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&col_b_evals));
    let row_b_poly = label_polynomial!(row_b_poly);
    let col_b_poly = label_polynomial!(col_b_poly);

    let row_c_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&row_c_evals));
    let col_c_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&col_c_evals));
    let val_c_poly = DensePolynomial::<Fr>::from_coefficients_slice(&domain_k.ifft(&val_c_evals));
    let row_c_poly = label_polynomial!(row_c_poly);
    let col_c_poly = label_polynomial!(col_c_poly);
    let val_c_poly = label_polynomial!(val_c_poly);

    vec![
        row_a_poly, col_a_poly, row_b_poly, col_b_poly, row_c_poly, col_c_poly, val_c_poly,
    ]
}
