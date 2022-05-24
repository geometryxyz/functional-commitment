//! Useful commitment stuff
use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine};
use ark_ff::{Field, PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::{sonic_pc::SonicKZG10, PolynomialCommitment};

/// A homomorphic polynomial commitment
pub trait HomomorphicPolynomialCommitment<F>: PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
    Self::VerifierKey: core::fmt::Debug,
{
    /// Combine a linear combination of homomorphic commitments
    fn multi_scalar_mul(commitments: &[Self::Commitment], scalars: &[F]) -> Self::Commitment;
}

/// The Default KZG-style commitment scheme
pub type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;
/// A single KZG10 commitment
pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Commitment;

impl<E> HomomorphicPolynomialCommitment<E::Fr> for KZG10<E>
where
    E: PairingEngine,
{
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
}
