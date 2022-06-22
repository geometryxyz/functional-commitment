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

    PCError { error: String },

    // In various protocols
    FsRngAbsorbError,

    InvalidDescriptionError,
    EvaluationError,
    InstantiationError,
    InvalidGeoSeq,
    T2Large,

    ProofSerializationError,
    ProofDeserializationError,
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
