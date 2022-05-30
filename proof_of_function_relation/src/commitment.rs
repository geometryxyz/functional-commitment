//! Useful commitment stuff
use crate::to_poly;
use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine};
use ark_ff::{Field, PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::PCRandomness;
use ark_poly_commit::{sonic_pc::SonicKZG10, PolynomialCommitment};
use std::ops::Add;

/// A homomorphic polynomial commitment
pub trait HomomorphicPolynomialCommitment<F>: PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
    Self::VerifierKey: core::fmt::Debug,
{
    // type HomomorphicRandomness: Add;

    /// Combine a linear combination of homomorphic commitments
    fn multi_scalar_mul(commitments: &[Self::Commitment], scalars: &[F]) -> Self::Commitment;
    fn aggregate_randomness(rands: &[Self::Randomness]) -> Self::Randomness;
}

/// The Default KZG-style commitment scheme
pub type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;

/// A single KZG10 commitment
pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Commitment;

pub type KZGRandomness<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Randomness;

impl<E> HomomorphicPolynomialCommitment<E::Fr> for KZG10<E>
where
    E: PairingEngine,
{
    // type HomomorphicRandomness = KZGRandomness<E>;

    fn multi_scalar_mul(
        commitments: &[KZG10Commitment<E>],
        scalars: &[E::Fr],
    ) -> KZG10Commitment<E> {
        let scalars_repr = scalars
            .iter()
            .map(<E::Fr as PrimeField>::into_repr)
            .collect::<Vec<_>>();

        let points_repr = commitments.iter().map(|c| c.0).collect::<Vec<_>>();

        ark_poly_commit::kzg10::Commitment::<E>(
            VariableBaseMSM::multi_scalar_mul(&points_repr, &scalars_repr).into(),
        )
    }

    fn aggregate_randomness(rands: &[KZGRandomness<E>]) -> KZGRandomness<E> {
        if rands.len() == 0 {
            return KZGRandomness::<E>::empty();
        }
        let mut acc = rands[0].clone();
        for rand in rands.iter().skip(1) {
            acc += rand;
        }

        acc
    }
}

/// Aggregate polynomials with separation challenge: p1(x) + c*p2(x) + ... + c^n*p(x)^n
pub fn aggregate_polynomials<F: Field>(
    polynomials: &[DensePolynomial<F>],
    evals: &[F],
    challenge: F,
) -> DensePolynomial<F> {
    assert_eq!(evals.len(), polynomials.len());
    let powers = crate::util::powers_of(challenge)
        .take(evals.len())
        .collect::<Vec<_>>();

    let mut aggregated_poly = DensePolynomial::<F>::zero();
    for ((poly, &eval), &power) in polynomials.iter().zip(evals.iter()).zip(powers.iter()) {
        let ith_poly = &(poly - &to_poly!(eval)) * power;
        aggregated_poly = aggregated_poly + ith_poly;
    }

    aggregated_poly
}
