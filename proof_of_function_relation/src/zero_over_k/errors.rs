use crate::error::Error;
use crate::non_zero_over_k;

#[derive(Debug, PartialEq)]
pub enum ZeroOverKError {
    ProverInitError,
    ProverFirstRoundError,
    VerifierInitError,
    VerifierFirstRoundError,
    VerifierQuerySetError,
    PCErrorR,
    PCErrorM,
    PCErrorQ1,
    PCErrorBatchOpen,
}
