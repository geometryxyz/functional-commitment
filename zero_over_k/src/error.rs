#[derive(Debug, PartialEq)]
pub enum Error {
    // In zero_over_k
    BatchCheckError,

    // In zero_over_k
    Check1Failed,
    Check2Failed,

    // In non_zero_over_k
    FEvalIsZero,
    FPrimeEvalError,

    // In zero_over_k,
    ToBytesError,

    PCError {
        error: String,
    },

    /// A commitment could not be found when evaluating a linear combination
    MissingCommitment(String),
    InputLengthError(String),
    MismatchedDegreeBounds(String),

    UnsupportedDegree(String),

    VOFailedToInstantiate,
    VOFailedToCompute,
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

impl From<homomorphic_poly_commit::error::Error> for Error {
    fn from(err: homomorphic_poly_commit::error::Error) -> Self {
        Self::PCError {
            error: format!("{:?}", err),
        }
    }
}
