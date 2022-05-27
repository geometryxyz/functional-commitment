use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use rand::RngCore;

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

#[macro_export]
macro_rules! label_polynomial_with_bound {
    ($poly:expr, $hiding_bound:expr) => {
        ark_poly_commit::LabeledPolynomial::new(
            stringify!($poly).to_owned(),
            $poly.clone(),
            None,
            $hiding_bound,
        )
    };
}

#[macro_export]
macro_rules! label_commitment {
    ($comm:expr) => {
        ark_poly_commit::LabeledCommitment::new(stringify!($comm).to_owned(), $comm.clone(), None)
    };
}

pub fn commit_polynomial<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>>(
    ck: &PC::CommitterKey,
    poly: &DensePolynomial<F>,
    rng: Option<&mut dyn RngCore>,
) -> Result<(PC::Commitment, PC::Randomness), Error> {
    let labeled_poly = label_polynomial!(poly);

    let (poly_commit, randomness) =
        PC::commit(ck, &[labeled_poly], rng).map_err(to_pc_error::<F, PC>)?;

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
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        label_commitment, label_polynomial,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_ff::Field;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_poly_commit::{kzg10, PCRandomness, PolynomialCommitment};
    use ark_std::UniformRand;
    use rand::thread_rng;
    use rand_core::OsRng;

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
    fn test_shift_and_unshift() {
        let rng = &mut thread_rng();
        let degree = 10;
        // let domain = GeneralEvaluationDomain::<Fr>::new(degree).unwrap();

        let r_poly = DensePolynomial::<Fr>::rand(degree, rng);
        // let r_poly = domain.vanishing_polynomial().into();

        let shift = Fr::rand(rng);
        let r_shifted = shift_dense_poly(&r_poly, &shift);

        let evaluation_point = Fr::rand(rng);
        let inverse_shift = shift.inverse().unwrap();

        assert_eq!(
            r_poly.evaluate(&evaluation_point),
            r_shifted.evaluate(&(evaluation_point * inverse_shift))
        )
    }

    #[test]
    fn check_empty_randomness() {
        let degree: usize = 16;

        let pp = KZG10::<Bn254>::setup(degree, None, &mut OsRng).unwrap();
        let (ck, vk) = KZG10::<Bn254>::trim(&pp, degree, 1, None).unwrap();

        let a_coeffs = vec![
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
            Fr::from(4 as u64),
            Fr::from(5 as u64),
        ];
        let a_poly = DensePolynomial::<Fr>::from_coefficients_vec(a_coeffs);

        let b_coeffs = vec![
            Fr::from(4 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
            Fr::from(4 as u64),
            Fr::from(5 as u64),
        ];
        let b_poly = DensePolynomial::<Fr>::from_coefficients_vec(b_coeffs);

        let point = Fr::from(15 as u64);
        let a_plus_b_poly = a_poly.clone() + b_poly.clone();
        let a_b_eval = a_plus_b_poly.evaluate(&point);

        let (a_commit, a_rand) = KZG10::<Bn254>::commit(
            &ck,
            &[label_polynomial_with_bound!(a_poly, Some(1))],
            Some(&mut OsRng),
        )
        .unwrap();

        let (b_commit, b_rand) = KZG10::<Bn254>::commit(
            &ck,
            &[label_polynomial_with_bound!(b_poly, Some(1))],
            Some(&mut OsRng),
        )
        .unwrap();

        let one = Fr::from(1 as u64);
        let homomorphic_commitment = KZG10::<Bn254>::multi_scalar_mul(
            &[
                a_commit[0].commitment().clone(),
                b_commit[0].commitment().clone(),
            ],
            &[one, one],
        );

        // let homomorphic_randomness = kzg10::Randomness::<Fr, DensePolynomial<Fr>>::empty();
        let homomorphic_randomness = a_rand[0].clone() + &b_rand[0];

        let proof = KZG10::<Bn254>::open(
            &ck,
            &[label_polynomial!(a_plus_b_poly)],
            &[label_commitment!(homomorphic_commitment)],
            &point,
            Fr::from(1 as u64),
            &[homomorphic_randomness],
            None,
        )
        .unwrap();

        let res = KZG10::<Bn254>::check(
            &vk,
            &[label_commitment!(homomorphic_commitment)],
            &point,
            [a_b_eval],
            &proof,
            Fr::from(1 as u64),
            None,
        )
        .unwrap();

        assert_eq!(res, true);
    }
}
