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

#[macro_export]
macro_rules! label_commitment {
    ($comm:expr) => {
        ark_poly_commit::LabeledCommitment::new(stringify!($comm).to_owned(), $comm.clone(), None)
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
    use crate::{
        commitment::{HomomorphicPolynomialCommitment, KZG10},
        label_commitment, label_polynomial,
        util::commit_polynomial,
    };
    use ark_bn254::{Bn254, Fr};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use ark_poly_commit::{
        kzg10, sonic_pc::SonicKZG10, LabeledCommitment, PCRandomness, PolynomialCommitment,
    };
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
    fn check_empty_randomness() {
        let degree: usize = 16;

        let pp = KZG10::<Bn254>::setup(degree, None, &mut OsRng).unwrap();
        let (ck, vk) = KZG10::<Bn254>::trim(&pp, degree, 0, None).unwrap();

        let a_coeffs = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
        ];
        let a_poly = DensePolynomial::<Fr>::from_coefficients_vec(a_coeffs);

        let b_coeffs = vec![
            Fr::from(4),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
        ];
        let b_poly = DensePolynomial::<Fr>::from_coefficients_vec(b_coeffs);

        let point = Fr::from(15 as u64);
        let a_plus_b_poly = a_poly.clone() + b_poly.clone();
        let a_b_eval = a_plus_b_poly.evaluate(&point);

        let (a_commit, _) =
            KZG10::<Bn254>::commit(&ck, &[label_polynomial!(a_poly)], None).unwrap();

        let (b_commit, _) =
            KZG10::<Bn254>::commit(&ck, &[label_polynomial!(b_poly)], None).unwrap();

        let one = Fr::from(1 as u64);
        let homomorphic_commitment = KZG10::<Bn254>::multi_scalar_mul(
            &[
                a_commit[0].commitment().clone(),
                b_commit[0].commitment().clone(),
            ],
            &[one, one],
        );

        let homomorphic_randomness = kzg10::Randomness::<Fr, DensePolynomial<Fr>>::empty();

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
