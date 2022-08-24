use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{LabeledCommitment, LinearCombination, PolynomialCommitment};

use crate::error::Error;

pub mod error;
pub mod marlin_kzg;
pub mod sonic_kzg;

/// An additively homomorphic polynomial commitment scheme
pub trait AdditivelyHomomorphicPCS<F>: PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
    Self::VerifierKey: core::fmt::Debug,
{
    /// Compute the linear combination of the provided commitments
    fn get_commitments_lc(
        commitments: &[LabeledCommitment<Self::Commitment>],
        lc: &LinearCombination<F>,
    ) -> Result<LabeledCommitment<Self::Commitment>, Error>;

    /// Compute the commitment and randomness that result of the linear combination of the provided commtiments and randomness values.
    /// We assume that commitments and randomness match 1-to-1 and are in the same order.
    fn get_commitments_lc_with_rands(
        commitments: &[LabeledCommitment<Self::Commitment>],
        hiding_rands: &[Self::Randomness],
        lc: &LinearCombination<F>,
    ) -> Result<(LabeledCommitment<Self::Commitment>, Self::Randomness), Error>;

    /// Aggregate labeled commitments according to the linear combination. If hiding bounds are enforced, the committer is expected to provide
    /// a vector of hiding randomness values, otherwise use `None`. A verifier can always aggregate with `None` for randomness.
    fn aggregate_commitments(
        commitments: &[LabeledCommitment<Self::Commitment>],
        randomness: Option<Vec<Self::Randomness>>,
        lc: &LinearCombination<F>,
    ) -> Result<(LabeledCommitment<Self::Commitment>, Self::Randomness), Error>;
}
