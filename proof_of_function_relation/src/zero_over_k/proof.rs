use crate::commitment::HomomorphicPolynomialCommitment;
use ark_ff::PrimeField;

pub struct Proof<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> {
    // commitments
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
