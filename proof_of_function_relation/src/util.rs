use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};

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

pub fn commit_polynomial<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>>(
    ck: &PC::CommitterKey,
    poly: &DensePolynomial<F>,
) -> Result<(PC::Commitment, PC::Randomness), Error> {
    let labeled_poly = label_polynomial!(poly);

    let (poly_commit, randomness) =
        PC::commit(ck, &[labeled_poly], None).map_err(to_pc_error::<F, PC>)?;

    Ok((poly_commit[0].commitment().clone(), randomness[0].clone()))
}

pub fn shift_dense_poly<F: Field>(
    p: &DensePolynomial<F>,
    shifting_factor: &F,
) -> DensePolynomial<F> {
    let mut coeffs = p.coeffs().to_vec();
    let mut acc = F::one();
    for i in 0..coeffs.len() {
        coeffs[i] = coeffs[i] * acc;
        acc *= shifting_factor;
    }

    DensePolynomial::from_coefficients_vec(coeffs)
}

#[cfg(test)]
mod test {
    use super::shift_dense_poly;
    use ark_bn254::Fr;
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use ark_std::UniformRand;
    use rand::thread_rng;

    #[test]
    fn test_dense_shifting() {
        let rng = &mut thread_rng();
        let degree = 10;

        let r_poly = DensePolynomial::<Fr>::rand(degree, rng);

        let shift = Fr::rand(rng);
        let r_shifted = shift_dense_poly(&r_poly, &shift);

        let evaluation_point = Fr::rand(rng);

        assert_eq!(
            r_poly.evaluate(&(evaluation_point * shift)),
            r_shifted.evaluate(&evaluation_point)
        )
    }

    #[test]
    fn test_v_h_masking() {}
}
