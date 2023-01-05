#[derive(Debug, PartialEq)]
pub enum Error {
    PCError {
        error: String,
    },

    /// A commitment could not be found when evaluating a linear combination
    MissingCommitment(String),
    InputLengthError(String),
    MismatchedDegreeBounds(String),
    ConstantTermInAggregation,
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
