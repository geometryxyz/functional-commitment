use ark_ff::{Field, PrimeField, FftField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial, GeneralEvaluationDomain, EvaluationDomain, Evaluations};
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

pub fn generate_sequence<F: PrimeField>(r: F, a_s: &[F], c_s: &[usize]) -> Vec<F> {
    assert_eq!(c_s.len(), a_s.len());
    let mut concatenation = Vec::<F>::default();

    for (i, a) in a_s.iter().enumerate() {
        for j in 0..c_s[i] {
            // r ** j (is there a pow() library function in Fp256?)
            let mut r_pow = F::from(1 as u64);
            for _ in 0..j {
                r_pow = r_pow * r;
            }
            let val = *a * r_pow;
            concatenation.push(val);
        }
    }

    concatenation
}

/// Sample a vector of random elements of type T
pub fn sample_vector<T: UniformRand, R: Rng>(seed: &mut R, length: usize) -> Vec<T> {
    (0..length)
        .collect::<Vec<usize>>()
        .iter()
        .map(|_| T::rand(seed))
        .collect::<Vec<_>>()
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
