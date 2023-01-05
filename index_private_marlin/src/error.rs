use crate::ahp::Error as AHPError;
use ::zero_over_k::error::Error as ZeroOverKError;

/// A `enum` specifying the possible failure modes of the `SNARK`.
#[derive(Debug)]
pub enum Error<E> {
    /// The index is too large for the universal public parameters.
    IndexTooLarge,
    /// There was an error in the underlying holographic IOP.
    AHPError(AHPError),
    /// There was an error in the underlying polynomial commitment.
    PolynomialCommitmentError(E),

    /// ZeroOverKError
    ZeroOverKError(ZeroOverKError),

    /// Number of constraints is larger than number of non zero elements (we don't allow this because Discrete-log Comparison in proof of function fails)
    DomainHLargerThanDomainK,
    DomainTooLarge,
}

impl<E> From<AHPError> for Error<E> {
    fn from(err: AHPError) -> Self {
        Error::AHPError(err)
    }
}

impl<E> From<ZeroOverKError> for Error<E> {
    fn from(err: ZeroOverKError) -> Self {
        Error::ZeroOverKError(err)
    }
}

impl<E> Error<E> {
    /// Convert an error in the underlying polynomial commitment scheme
    /// to a `Error`.
    pub fn from_pc_err(err: E) -> Self {
        Error::PolynomialCommitmentError(err)
    }
}
