#[derive(Debug, PartialEq)]
pub enum NonZeroOverKError {
    FEvalIsZero,
    PCError,
    ZeroOverKProofError,
}
