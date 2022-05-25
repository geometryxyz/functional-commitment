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

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    // commitments
    pub concrete_oracles_commitments: Vec<PC::Commitment>,
    pub m_commitments: Vec<PC::Commitment>,
    pub r_commitments: Vec<PC::Commitment>,
    pub q1_commit: PC::Commitment,

    // evaluations
    pub q1_eval: F,
    pub q2_eval: F,
    pub h_prime_evals: Vec<F>,
    pub m_evals: Vec<F>,

    // opening
    pub q1_opening: PC::Proof,
    pub q2_opening: PC::Proof,
    pub m_openings: Vec<PC::Proof>,
    pub h_prime_openings: Vec<PC::Proof>,
}
