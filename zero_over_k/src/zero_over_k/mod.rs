#![allow(dead_code)]

use crate::error::{to_pc_error, Error};
use crate::get_labels;
use crate::util::powers_of;
use crate::virtual_oracle::{generic_shifting_vo::vo_term::VOTerm, VirtualOracle};
use crate::zero_over_k::piop::PIOPforZeroOverK;
use crate::zero_over_k::proof::Proof;
use ark_ff::to_bytes;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::Evaluations;
use ark_poly_commit::{
    data_structures::{PCCommitterKey, PCVerifierKey},
    LabeledCommitment, LabeledPolynomial,
};
use ark_std::marker::PhantomData;
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use rand_core::OsRng;
use std::iter;

mod piop;
pub mod proof;
mod tests;

/// zk-SNARK to prove that a virtual oracle evaluates to 0 over a given domain
pub struct ZeroOverK<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> {
    _field: PhantomData<F>,
    _polynomial_commitment_scheme: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> ZeroOverK<F, PC, FS> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Zero Over K";

    pub fn prove<R: Rng, VO: VirtualOracle<F>>(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        concrete_oracle_commitments: &[LabeledCommitment<PC::Commitment>],
        concrete_oracle_commit_rands: &[PC::Randomness],
        maximum_oracle_degree_bound: Option<usize>,
        virtual_oracle: &VO,
        domain: &GeneralEvaluationDomain<F>,
        ck: &PC::CommitterKey,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        if let Some(degree) = maximum_oracle_degree_bound {
            if degree > ck.supported_degree() {
                return Err(Error::UnsupportedDegree(format!(
                    "Cannot produce proof for max degree {} with the provided committer key",
                    degree
                )));
            }
        };

        let alphas = virtual_oracle.shifting_coefficients();

        let prover_initial_state = PIOPforZeroOverK::prover_init(
            domain,
            concrete_oracles,
            maximum_oracle_degree_bound,
            virtual_oracle,
            &alphas,
        )?;
        let verifier_initial_state = PIOPforZeroOverK::<F, VO>::verifier_init(
            virtual_oracle,
            maximum_oracle_degree_bound,
            domain,
        )?;

        let fs_bytes = &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments, alphas]
            .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FS::initialize(fs_bytes);

        //------------------------------------------------------------------
        // First Round
        let (_, prover_first_oracles, prover_state) =
            PIOPforZeroOverK::prover_first_round(prover_initial_state, rng)?;

        let random_polynomials = prover_first_oracles.random_polynomials.clone();
        let masking_polynomials = prover_first_oracles.masking_polynomials.clone();
        let q_1 = prover_first_oracles.q_1.clone();

        // commit to the random polynomials
        let (r_commitments, r_rands) =
            PC::commit(ck, random_polynomials.iter(), Some(rng)).map_err(to_pc_error::<F, PC>)?;

        // commit to the masking polynomials
        let (m_commitments, m_rands) =
            PC::commit(ck, masking_polynomials.iter(), Some(rng)).map_err(to_pc_error::<F, PC>)?;

        // commit to q_1
        let (q1_commit, q1_rand) =
            PC::commit(ck, &[q_1], Some(rng)).map_err(to_pc_error::<F, PC>)?;

        let fs_bytes =
            &to_bytes![r_commitments, m_commitments, q1_commit].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let (verifier_first_msg, verifier_state) =
            PIOPforZeroOverK::<F, VO>::verifier_first_round(verifier_initial_state, &mut fs_rng)?;
        //------------------------------------------------------------------

        //------------------------------------------------------------------
        // Second Round

        let (_prover_second_msg, prover_second_oracles, prover_state) =
            PIOPforZeroOverK::prover_second_round(&verifier_first_msg, prover_state, rng);

        let query_set = PIOPforZeroOverK::<F, VO>::verifier_query_set(&verifier_state, &alphas)?;
        //------------------------------------------------------------------

        let h_primes = prover_state
            .masked_oracles
            .expect("Prover should have computed masked oracles");

        let polynomials = prover_first_oracles
            .masking_polynomials
            .iter()
            .chain(h_primes.iter())
            .chain(iter::once(&prover_first_oracles.q_1))
            .chain(iter::once(&prover_second_oracles.q_2));

        // it gives us ((poly_label, point), evaluation)
        let evaluated_query_set =
            ark_poly_commit::evaluate_query_set(polynomials.clone(), &query_set);

        let mut h_prime_evals = Vec::new();
        let mut m_evals = Vec::new();
        let mut q1_eval: Option<F> = None;
        let mut q2_eval: Option<F> = None;

        // Extract evaluations from the `evaluated_query_set`
        for ((poly_label, _), &evaluation) in &evaluated_query_set {
            if poly_label.contains("h_prime_") {
                h_prime_evals.push((poly_label, evaluation));
            } else if poly_label.contains("m_") {
                m_evals.push((poly_label, evaluation));
            } else if poly_label == "q_1" {
                q1_eval = Some(evaluation);
            } else if poly_label == "q_2" {
                q2_eval = Some(evaluation);
            }
        }

        // sort the evaluation by poly label name => h_prime_is, m_is, q_1, q_2
        h_prime_evals.sort_by(|a, b| a.0.cmp(b.0));
        let h_prime_evals = h_prime_evals.into_iter().map(|x| x.1).collect::<Vec<F>>();
        m_evals.sort_by(|a, b| a.0.cmp(b.0));
        let m_evals = m_evals.into_iter().map(|x| x.1).collect::<Vec<F>>();

        // sanity checks
        let q1_eval = q1_eval.expect("q_1 was not evaluated");
        let q2_eval = q2_eval.expect("q_2 was not evaluated");
        assert_eq!(h_prime_evals.len(), m_evals.len());

        let fs_bytes = &to_bytes![h_prime_evals, m_evals, q1_eval, q2_eval]
            .map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        // Use commitments to all the random polynomials to compute a commitment to q2
        let q2_linear_combination =
            PIOPforZeroOverK::generate_q2_linear_combination(virtual_oracle, verifier_first_msg.c);

        let (q2_commit, q2_rand) = PC::aggregate_commitments(
            &r_commitments,
            Some(r_rands.to_vec()),
            &q2_linear_combination,
        )?;

        // Use commitments to each h and each m to homomorphically derive commitments to each h_prime
        let h_prime_lcs = PIOPforZeroOverK::generate_h_prime_linear_combinations(
            virtual_oracle,
            &get_labels!(concrete_oracles),
        );

        let h_and_m_commitments: Vec<_> = concrete_oracle_commitments
            .iter()
            .cloned()
            .chain(m_commitments.iter().cloned())
            .collect();
        let h_and_m_rands: Vec<_> = concrete_oracle_commit_rands
            .iter()
            .cloned()
            .chain(m_rands.iter().cloned())
            .collect();

        let mut h_prime_commitments = Vec::new();
        let mut h_prime_rands = Vec::new();
        for lc in h_prime_lcs {
            let (comm, rand) =
                PC::aggregate_commitments(&h_and_m_commitments, Some(h_and_m_rands.to_vec()), &lc)?;
            h_prime_commitments.push(comm);
            h_prime_rands.push(rand);
        }

        let fs_bytes =
            &to_bytes![h_prime_commitments, q2_commit].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let commitments = m_commitments
            .iter()
            .chain(h_prime_commitments.iter())
            .chain(q1_commit.iter())
            .chain(iter::once(&q2_commit));

        let rands = m_rands
            .iter()
            .chain(h_prime_rands.iter())
            .chain(q1_rand.iter())
            .chain(iter::once(&q2_rand));

        let separation_challenge = F::rand(&mut fs_rng);

        let batch_opening = PC::batch_open(
            ck,
            polynomials.clone(),
            commitments,
            &query_set,
            separation_challenge,
            rands,
            Some(rng),
        )
        .map_err(to_pc_error::<F, PC>)?;

        let proof = Proof {
            // commitments
            m_commitments: m_commitments
                .iter()
                .map(|p| p.commitment().clone())
                .collect(),
            r_commitments: r_commitments
                .iter()
                .map(|p| p.commitment().clone())
                .collect(),
            q1_commit: q1_commit[0].commitment().clone(),

            //evaluations
            q1_eval,
            q2_eval,
            h_prime_evals,
            m_evals,

            opening_proof: batch_opening,
        };

        Ok(proof)
    }

    pub fn verify<VO: VirtualOracle<F>>(
        proof: Proof<F, PC>,
        concrete_oracle_commitments: &[LabeledCommitment<PC::Commitment>],
        maximum_oracle_degree_bound: Option<usize>,
        virtual_oracle: &VO,
        domain: &GeneralEvaluationDomain<F>,
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        if let Some(degree) = maximum_oracle_degree_bound {
            if degree > vk.supported_degree() {
                return Err(Error::UnsupportedDegree(format!(
                    "Cannot verify proof for max degree {} with the provided verifier key",
                    degree
                )));
            }
        };

        let alphas = virtual_oracle.shifting_coefficients();

        let verifier_initial_state = PIOPforZeroOverK::<F, VO>::verifier_init(
            virtual_oracle,
            maximum_oracle_degree_bound,
            domain,
        )?;

        let fs_bytes = &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments, alphas]
            .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FS::initialize(fs_bytes);

        //------------------------------------------------------------------
        // First Round
        let fs_bytes = &to_bytes![proof.r_commitments, proof.m_commitments, proof.q1_commit]
            .map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let (verifier_first_msg, verifier_state) =
            PIOPforZeroOverK::<F, VO>::verifier_first_round(verifier_initial_state, &mut fs_rng)?;

        //------------------------------------------------------------------
        // Second Round
        let fs_bytes = &to_bytes![
            proof.h_prime_evals,
            proof.m_evals,
            proof.q1_eval,
            proof.q2_eval
        ]
        .map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let query_set = PIOPforZeroOverK::<F, VO>::verifier_query_set(&verifier_state, &alphas)?;

        let beta_1 = verifier_state
            .beta_1
            .expect("Verifier did not produce beta_1");
        let beta_2 = verifier_state
            .beta_2
            .expect("Verifier did not produce beta_2");

        let r_commitments = proof
            .r_commitments
            .iter()
            .enumerate()
            .map(|(i, c)| LabeledCommitment::new(format!("r_{}", i), c.clone(), Some(2)))
            .collect::<Vec<_>>();

        let m_commitments = proof
            .m_commitments
            .iter()
            .enumerate()
            .map(|(i, c)| {
                LabeledCommitment::new(format!("m_{}", i), c.clone(), maximum_oracle_degree_bound)
            })
            .collect::<Vec<_>>();

        // derive commitment to h_prime through additive homomorphism
        let h_and_m_commitments: Vec<_> = concrete_oracle_commitments
            .iter()
            .cloned()
            .chain(m_commitments.iter().cloned())
            .collect();

        let h_prime_lcs = PIOPforZeroOverK::generate_h_prime_linear_combinations(
            virtual_oracle,
            &get_labels!(concrete_oracle_commitments),
        );
        let (h_prime_commitments, _): (Vec<_>, Vec<_>) = h_prime_lcs
            .iter()
            .map(|lc| PC::aggregate_commitments(&h_and_m_commitments, None, lc))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .cloned()
            .unzip();

        let q1_commit = LabeledCommitment::new(String::from("q_1"), proof.q1_commit, None);

        // derive commitment to q2 through additive homomorphism
        let q2_linear_combination =
            PIOPforZeroOverK::generate_q2_linear_combination(virtual_oracle, verifier_first_msg.c);

        let (q2_commit, _) =
            PC::aggregate_commitments(&r_commitments, None, &q2_linear_combination)?;

        let fs_bytes =
            &to_bytes![h_prime_commitments, q2_commit].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        // concatenate all the evaluations in alphabetical order
        let evals: Vec<_> = proof
            .h_prime_evals
            .iter()
            .chain(proof.m_evals.iter())
            .chain(iter::once(&proof.q1_eval))
            .chain(iter::once(&proof.q2_eval))
            .collect();

        let mut evaluations = Evaluations::new();
        let mut evaluation_labels = Vec::new();
        for (poly_label, (_, point)) in query_set.iter().cloned() {
            evaluation_labels.push((poly_label, point));
        }

        // take advantage of the alphabetical order to re-label each evaluation
        evaluation_labels.sort_by(|a, b| a.0.cmp(&b.0));
        for (q, eval) in evaluation_labels.into_iter().zip(evals) {
            evaluations.insert(q, *eval);
        }

        let commitments = m_commitments
            .iter()
            .chain(h_prime_commitments.iter())
            .chain(iter::once(&q1_commit))
            .chain(iter::once(&q2_commit));

        let separation_challenge = F::rand(&mut fs_rng);

        match PC::batch_check(
            vk,
            commitments,
            &query_set,
            &evaluations,
            &proof.opening_proof,
            separation_challenge,
            &mut OsRng,
        ) {
            Ok(true) => Ok(()),
            Ok(false) => Err(Error::BatchCheckError),
            Err(e) => panic!("{:?}", e),
        }?;

        // compute M(beta_2)
        let big_m_at_beta_2 = &proof
            .m_evals
            .iter()
            .zip(powers_of(verifier_first_msg.c))
            .fold(F::zero(), |acc, (&m_eval, c_power)| {
                acc + (m_eval * c_power)
            });

        // evaluate z_k(beta_1), z_k(beta_2)
        let z_k_at_beta_1 = domain.evaluate_vanishing_polynomial(beta_1);
        let z_k_at_beta_2 = domain.evaluate_vanishing_polynomial(beta_2);

        // compute F_prime(beta_1)
        let f_prime_eval = compute_f_prime_eval(virtual_oracle, &proof.h_prime_evals, &beta_1)?;

        // check that M(beta_2) - q2(beta_2)*zK(beta_2) = 0
        let check_1 = *big_m_at_beta_2 - proof.q2_eval * z_k_at_beta_2;
        if check_1 != F::zero() {
            return Err(Error::Check1Failed);
        }

        // check that F_prime(beta_1) - q1(beta_1)*zK(beta_1) = 0
        let check_2 = f_prime_eval - proof.q1_eval * z_k_at_beta_1;
        if check_2 != F::zero() {
            return Err(Error::Check2Failed);
        }

        Ok(())
    }
}

fn compute_f_prime_eval<F: PrimeField, VO: VirtualOracle<F>>(
    virtual_oracle: &VO,
    evals: &[F],
    point: &F,
) -> Result<F, Error> {
    let terms: Vec<_> = vec![point.clone()]
        .iter()
        .chain(evals.iter())
        .map(|e| VOTerm::Evaluation(e.clone()))
        .collect();

    match virtual_oracle.apply_evaluation_function(&terms) {
        VOTerm::Evaluation(res) => Ok(res),
        VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
    }
}
