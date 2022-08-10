#[derive(Debug, PartialEq)]
#[allow(dead_code)]
pub enum Error {
    MaxDegreeExceeded,

    // In zero_over_k and geo_seq
    BatchCheckError,

    // In zero_over_k
    Check1Failed,
    Check2Failed,

    // In non_zero_over_k
    FEvalIsZero,
    FPrimeEvalError,

    // In zero_over_k, discrete_log_comparison,
    ToBytesError,

    // In discrete_log_comparison
    OmegaSqrtError,

    PCError {
        error: String,
    },

    // In various protocols
    FsRngAbsorbError,

    InvalidDescriptionError,
    EvaluationError,
    InstantiationError,
    InvalidGeoSeq,
    T2Large,

    ProofSerializationError,
    ProofDeserializationError,

    /// A commitment could not be found when evaluating a linear combination
    MissingCommitment(String),
    InputLengthError(String),
    MismatchedDegreeBounds(String),

    UnsupportedDegree(String),

    VOMismatchedVariants,
    VOFailedToInstantiate,
    VOFailedToCompute,

    ZeroOverKError(String),
}

/// Convert an ark_poly_commit error
pub fn to_pc_error<F, PC>(error: PC::Error) -> Error
where
    F: ark_ff::Field,
    PC: ark_poly_commit::PolynomialCommitment<F, ark_poly::univariate::DensePolynomial<F>>,
{
    println!("Polynomial Commitment Error: {:?}", error);
    Error::PCError {
        error: format!("Polynomial Commitment Error: {:?}", error),
    }
}

impl From<zero_over_k::error::Error> for Error {
    fn from(err: zero_over_k::error::Error) -> Self {
        Self::ZeroOverKError(format!("{:?}", err))
    }
}

impl From<homomorphic_poly_commit::error::Error> for Error {
    fn from(err: homomorphic_poly_commit::error::Error) -> Self {
        Self::PCError {
            error: format!("{:?}", err),
        }
    }
}
