#![allow(dead_code)]

use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use crate::util::powers_of;
use crate::virtual_oracle::{EvaluationsProvider, VirtualOracle};
use crate::zero_over_k::piop::PIOPforZeroOverK;
use crate::zero_over_k::proof::Proof;
use ark_ff::to_bytes;
use ark_ff::PrimeField;
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::Evaluations;
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, PCRandomness};
use ark_std::marker::PhantomData;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use rand_core::OsRng;
use std::iter;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use std::io::{BufReader, BufWriter};

mod piop;
pub mod proof;
mod tests;

pub struct ZeroOverK<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _field: PhantomData<F>,
    _polynomial_commitment_scheme: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> ZeroOverK<F, PC, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Zero Over K";

    pub fn prove<R: Rng, VO: VirtualOracle<F>>(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        concrete_oracle_commitments: &[LabeledCommitment<PC::Commitment>],
        _concrete_oracle_commit_rands: &[PC::Randomness],
        virtual_oracle: &VO,
        alphas: &Vec<F>,
        domain: &GeneralEvaluationDomain<F>,
        ck: &PC::CommitterKey,
        rng: &mut R,
    ) -> Result<Vec<u8>, Error> {
    //) -> Result<Proof<F, PC>, Error> {
        let prover_initial_state =
            PIOPforZeroOverK::prover_init(domain, concrete_oracles, virtual_oracle, &alphas)?;
        let verifier_initial_state =
            PIOPforZeroOverK::<F, VO>::verifier_init(virtual_oracle, domain)?;

        let fs_bytes = &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments, alphas]
            .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FiatShamirRng::<D>::from_seed(fs_bytes);

        //------------------------------------------------------------------
        // First Round
        let (_, prover_first_oracles, prover_state) =
            PIOPforZeroOverK::prover_first_round(prover_initial_state, rng)?;

        let random_polynomials = prover_first_oracles.random_polynomials.clone();
        let masking_polynomials = prover_first_oracles.masking_polynomials.clone();
        let q_1 = prover_first_oracles.q_1.clone();

        // commit to the random polynomials
        let (r_commitments, _) =
            PC::commit(ck, random_polynomials.iter(), None).map_err(to_pc_error::<F, PC>)?;

        // commit to the masking polynomials
        let (m_commitments, _m_rands) =
            PC::commit(ck, masking_polynomials.iter(), None).map_err(to_pc_error::<F, PC>)?;

        // commit to q_1
        let (q1_commit, _q1_rand) =
            PC::commit(ck, &[q_1.clone()], None).map_err(to_pc_error::<F, PC>)?;

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

        let query_set = PIOPforZeroOverK::<F, VO>::verifier_query_set(&verifier_state, alphas)?;
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
        h_prime_evals.sort_by(|a, b| a.0.cmp(&b.0));
        let h_prime_evals = h_prime_evals.into_iter().map(|x| x.1).collect::<Vec<F>>();
        m_evals.sort_by(|a, b| a.0.cmp(&b.0));
        let m_evals = m_evals.into_iter().map(|x| x.1).collect::<Vec<F>>();

        // sanity checks
        let q1_eval = q1_eval.expect("q_1 was not evaluated");
        let q2_eval = q2_eval.expect("q_2 was not evaluated");
        assert_eq!(h_prime_evals.len(), m_evals.len());

        let fs_bytes = &to_bytes![h_prime_evals, m_evals, q1_eval, q2_eval]
            .map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let q2_commit = PC::multi_scalar_mul(
            &r_commitments,
            powers_of(verifier_first_msg.c)
                .take(random_polynomials.len())
                .collect::<Vec<_>>()
                .as_slice(),
        );
        let q2_commit = LabeledCommitment::new(String::from("q_2"), q2_commit, None);

        let fs2hs = virtual_oracle.mapping_vector();
        let mut h_commitments = Vec::with_capacity(virtual_oracle.num_of_oracles());
        for concrete_oracle_index in fs2hs {
            h_commitments.push(concrete_oracle_commitments[concrete_oracle_index].clone())
        }

        let h_prime_commitments = h_commitments
            .iter()
            .zip(m_commitments.iter())
            .enumerate()
            .map(|(i, (h_commitment, m))| {
                let c =
                    PC::multi_scalar_mul(&[h_commitment.clone(), m.clone()], &[F::one(), F::one()]);
                // In order to work with batched version of PC, commitment labels must be same as poly labels
                LabeledCommitment::new(format!("h_prime_{}", i), c, None)
            })
            .collect::<Vec<_>>();

        let fs_bytes =
            &to_bytes![h_prime_commitments, q2_commit].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let commitments = m_commitments
            .iter()
            .chain(h_prime_commitments.iter())
            .chain(q1_commit.iter())
            .chain(iter::once(&q2_commit));

        let rands =
            vec![PC::Randomness::empty(); m_commitments.len() + h_prime_commitments.len() + 2];

        let separation_challenge = F::rand(&mut fs_rng);

        let batch_opening = PC::batch_open(
            ck,
            polynomials.clone(),
            commitments,
            &query_set,
            separation_challenge,
            &rands,
            Some(&mut fs_rng),
        )
        .map_err(to_pc_error::<F, PC>)?;

        // let f_prime = virtual_oracle.instantiate(&h_primes, alphas).unwrap();
        // println!("F PRIME EVAL ON PROVER SIDE: {}", f_prime.evaluate(&verifier_state.beta_1.expect("")));

        let proof = Proof::<F, PC> {
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

        let mut writer = BufWriter::new(Vec::new());
        let _ = proof.serialize(&mut writer).map_err(|_| Error::ProofSerializationError)?;

        Ok(Vec::from(writer.buffer()))
    }

    pub fn verify<VO: VirtualOracle<F>>(
        //proof: Proof<F, PC>,
        proof_bytes: Vec<u8>,
        concrete_oracle_commitments: &[LabeledCommitment<PC::Commitment>],
        virtual_oracle: &VO,
        domain: &GeneralEvaluationDomain<F>,
        alphas: &[F],
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        let reader = BufReader::new(proof_bytes.as_slice());
        let proof: Proof::<F, PC> = Proof::<F, PC>::deserialize(reader)
            .map_err(|_| Error::ProofDeserializationError)?;

        let verifier_initial_state =
            PIOPforZeroOverK::<F, VO>::verifier_init(virtual_oracle, &domain)?;

        let fs_bytes = &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments, alphas]
            .map_err(|_| Error::ToBytesError)?;
        let mut fs_rng = FiatShamirRng::<D>::from_seed(fs_bytes);

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

        let query_set = PIOPforZeroOverK::<F, VO>::verifier_query_set(&verifier_state, alphas)?;

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
            .map(|(i, c)| LabeledCommitment::new(format!("r_{}", i), c.clone(), None))
            .collect::<Vec<_>>();
        let m_commitments = proof
            .m_commitments
            .iter()
            .enumerate()
            .map(|(i, c)| LabeledCommitment::new(format!("m_{}", i), c.clone(), None))
            .collect::<Vec<_>>();

        let fs2hs = virtual_oracle.mapping_vector();
        let mut h_commitments = Vec::with_capacity(virtual_oracle.num_of_oracles());
        for concrete_oracle_index in fs2hs {
            h_commitments.push(concrete_oracle_commitments[concrete_oracle_index].clone())
        }

        // derive commitment to h_prime through additive homomorphism
        let h_prime_commitments = h_commitments
            .iter()
            .zip(m_commitments.iter())
            .enumerate()
            .map(|(i, (h_commitment, m_commitment))| {
                let c = PC::multi_scalar_mul(
                    &[h_commitment.clone(), m_commitment.clone()],
                    &[F::one(), F::one()],
                );
                LabeledCommitment::new(format!("h_prime_{}", i), c, None)
            })
            .collect::<Vec<_>>();

        // derive commitment to q2 through additive homomorphism
        let q2_commit = PC::multi_scalar_mul(
            &r_commitments,
            powers_of(verifier_first_msg.c)
                .take(h_prime_commitments.len())
                .collect::<Vec<_>>()
                .as_slice(),
        );

        let q1_commit = LabeledCommitment::new(String::from("q_1"), proof.q1_commit, None);
        let q2_commit = LabeledCommitment::new(String::from("q_2"), q2_commit, None);

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
        let f_prime_eval_r =
            proof
                .h_prime_evals
                .evaluate(virtual_oracle, beta_1, &Vec::<F>::default());

        if f_prime_eval_r.is_err() {
            return Err(Error::FPrimeEvalError);
        }
        let f_prime_eval = f_prime_eval_r.unwrap();

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
