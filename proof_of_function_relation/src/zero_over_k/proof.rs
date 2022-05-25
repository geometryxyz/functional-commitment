use crate::{
    commitment::HomomorphicPolynomialCommitment,
    error::{to_pc_error, Error},
    label_commitment, label_polynomial,
    transcript::TranscriptProtocol,
    util::{commit_polynomial, powers_of, shift_dense_poly},
    virtual_oracle::VirtualOracle,
};
use ark_ff::{PrimeField, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_poly_commit::{LabeledPolynomial, PCRandomness};
use ark_std::marker::PhantomData;
use merlin::Transcript;
use rand::Rng;

struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    // commitments
    concrete_oracle_commitments: Vec<PC::Commitment>,
    m_commitments: Vec<PC::Commitment>,
    r_commitments: Vec<PC::Commitment>,
    q1_commit: PC::Commitment,
    
    // evaluations
    q1_eval: F,
    q2_eval: F,
    h_prime_evals: Vec<F>,
    m_evals: Vec<F>,
    
    // opening
    q1_opening: PC::Proof,
    m_openings: Vec<PC::Proof>,
    h_prime_openings: Vec<PC::Proof>,
}