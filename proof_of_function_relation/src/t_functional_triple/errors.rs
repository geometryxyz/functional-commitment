use crate::error::Error;
use crate::non_zero_over_k;

#[derive(Debug, PartialEq)]
pub enum TftError {
    //TsltAProofError(non_zero_over_k::errors::NonZeroOverKError),
    TsltAProofError(Error),
    TsltBProofError(Error),
    TdiagProofError(Error),

    TsltAVerifyError(Error),
    TsltBVerifyError(Error),
    TdiagVerifyError(Error),
}
