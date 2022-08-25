#![cfg_attr(not(feature = "std"), no_std)]
//! A crate for the Marlin preprocessing zkSNARK for R1CS.
//!
//! # Note
//!
//! Currently, Marlin only supports R1CS instances where the number of inputs
//! is the same as the number of constraints (i.e., where the constraint
//! matrices are square). Furthermore, Marlin only supports instances where the
//! public inputs are of size one less than a power of 2 (i.e., 2^n - 1).
// #![deny(unused_import_braces, unused_qualifications, trivial_casts)]
// #![deny(trivial_numeric_casts, private_in_public)]
// #![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
// #![deny(unused_attributes, unused_imports, unused_mut, missing_docs)]
// #![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
// #![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
// #![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;

use core::iter;

use ark_ff::{to_bytes, PrimeField, UniformRand, One};
use ark_poly::UVPolynomial;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{Evaluations, LabeledPolynomial, PCRandomness};
use ark_poly_commit::{LabeledCommitment, PCUniversalParams, PolynomialCommitment};
use ark_relations::r1cs::ConstraintSynthesizer;
use ark_std::rand::RngCore;
// use crate::virtual_oracle::rational_sumcheck_vo::RationalSumcheckVO;
// use crate::virtual_oracle::{AddVO, rational_sumcheck_vo};
use ::zero_over_k::{zero_over_k::ZeroOverK, virtual_oracle::{generic_shifting_vo::{GenericShiftingVO, vo_term::VOTerm}}, vo_constant};
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use fiat_shamir_rng::{SimpleHashFiatShamirRng, FiatShamirRng};

use ark_std::{
    collections::BTreeMap,
    format,
    marker::PhantomData,
    string::{String, ToString},
    vec,
    vec::Vec,
};

#[cfg(not(feature = "std"))]
macro_rules! eprintln {
    () => {};
    ($($arg: tt)*) => {};
}

// /// Implements a Fiat-Shamir based Rng that allows one to incrementally update
// /// the seed based on new messages in the proof transcript.
// use rng::FiatShamirRng;
// pub use rng::SimpleHashFiatShamirRng;

mod error;
pub use error::*;

mod data_structures;
pub use data_structures::*;

/// Implements an Algebraic Holographic Proof (AHP) for the R1CS indexed relation.
pub mod ahp;
pub use ahp::AHPForR1CS;
use ahp::EvaluationsProvider;

// mod rational_sumcheck_vo;

#[cfg(test)]
mod test;

/// The compiled argument system.
pub struct Marlin<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng>(
    #[doc(hidden)] PhantomData<F>,
    #[doc(hidden)] PhantomData<PC>,
    #[doc(hidden)] PhantomData<FS>,
);

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng>
    Marlin<F, PC, FS>
{
    /// The personalization string for this protocol. Used to personalize the
    /// Fiat-Shamir rng.
    pub const PROTOCOL_NAME: &'static [u8] = b"MARLIN-2019";

    /// Generate the universal prover and verifier keys for the
    /// argument system.
    pub fn universal_setup<R: RngCore>(
        num_constraints: usize,
        num_variables: usize,
        num_non_zero: usize,
        rng: &mut R,
    ) -> Result<UniversalSRS<F, PC>, Error<PC::Error>> {
        let max_degree = AHPForR1CS::<F>::max_degree(num_constraints, num_variables, num_non_zero)?;
        let setup_time = start_timer!(|| {
            format!(
            "Marlin::UniversalSetup with max_degree {}, computed for a maximum of {} constraints, {} vars, {} non_zero",
            max_degree, num_constraints, num_variables, num_non_zero,
        )
        });

        let srs = PC::setup(max_degree, None, rng).map_err(Error::from_pc_err);
        end_timer!(setup_time);
        srs
    }

    /// Generate the index-specific (i.e., circuit-specific) prover and verifier
    /// keys. This is a deterministic algorithm that anyone can rerun.
    pub fn index<C: ConstraintSynthesizer<F>>(
        srs: &UniversalSRS<F, PC>,
        c: C,
    ) -> Result<(IndexProverKey<F, PC>, IndexVerifierKey<F, PC>), Error<PC::Error>> {
        let index_time = start_timer!(|| "Marlin::Index");

        // TODO: Add check that c is in the correct mode.
        let index = AHPForR1CS::index(c)?;
        if srs.max_degree() < index.max_degree() {
            Err(Error::IndexTooLarge)?;
        }

        let coeff_support = AHPForR1CS::get_degree_bounds(&index.index_info);
        // Marlin only needs degree 2 random polynomials
        let supported_hiding_bound = 1;
        let (committer_key, verifier_key) = PC::trim(
            &srs,
            index.max_degree(),
            supported_hiding_bound,
            Some(&coeff_support),
        )
        .map_err(Error::from_pc_err)?;

        let commit_time = start_timer!(|| "Commit to index polynomials");
        let (index_comms, index_comm_rands): (_, _) =
            PC::commit(&committer_key, index.iter(), None).map_err(Error::from_pc_err)?;
        end_timer!(commit_time);

        let index_comms = index_comms
            .into_iter()
            .map(|c| c.commitment().clone())
            .collect();
        let index_vk = IndexVerifierKey {
            index_info: index.index_info,
            index_comms,
            verifier_key,
        };

        let index_pk = IndexProverKey {
            index: index.clone(),
            index_comm_rands,
            index_vk: index_vk.clone(),
            committer_key: committer_key.clone(),
        };

        end_timer!(index_time);

        Ok((index_pk, index_vk))
    }

    /// Generate the index-specific (i.e., circuit-specific) prover and verifier
    /// keys. This is not a deterministic algorithm since circuit is known only by the prover.
    pub fn index_for_index_private<C: ConstraintSynthesizer<F>>(
        srs: &UniversalSRS<F, PC>,
        c: C,
    ) -> Result<(IndexPrivateProverKey<F, PC>, IndexPrivateVerifierKey<F, PC>), Error<PC::Error>> {
        let index_time = start_timer!(|| "Marlin::Index");

        // TODO: Add check that c is in the correct mode.
        // TODO: change this to work only with individual matrices
        let index = AHPForR1CS::index(c)?;
        if srs.max_degree() < index.max_degree() {
            Err(Error::IndexTooLarge)?;
        }

        let coeff_support = AHPForR1CS::get_degree_bounds(&index.index_info);
        // Marlin only needs degree 2 random polynomials
        let supported_hiding_bound = 1;
        let (committer_key, verifier_key) = PC::trim(
            &srs,
            index.max_degree(),
            supported_hiding_bound,
            Some(&coeff_support),
        )
        .map_err(Error::from_pc_err)?;

        let (individual_matrix_poly_commits, _): (_, _) =
        PC::commit(&committer_key, index.iter_individual_matrices(), None)
            .map_err(Error::from_pc_err)?;

        let individual_matrix_poly_commits = individual_matrix_poly_commits
            .iter()
            .map(|c| c.commitment().clone())
            .collect::<Vec<_>>();

        let index_private_vk: IndexPrivateVerifierKey<F, PC> = IndexPrivateVerifierKey {
            polys: individual_matrix_poly_commits,
            verifier_key,
            index_info: index.index_info,
        };

        let index_private_pk = IndexPrivateProverKey {
            index: index.clone(),
            index_private_vk: index_private_vk.clone(),
            committer_key: committer_key.clone(),
        };

        end_timer!(index_time);

        Ok((index_private_pk, index_private_vk))
    }

    /// Create a zkSNARK asserting that the constraint system is satisfied.
    pub fn prove<C: ConstraintSynthesizer<F>, R: RngCore>(
        index_pk: &IndexProverKey<F, PC>,
        c: C,
        zk_rng: &mut R,
    ) -> Result<Proof<F, PC>, Error<PC::Error>> {
        let prover_time = start_timer!(|| "Marlin::Prover");
        // Add check that c is in the correct mode.

        let prover_init_state = AHPForR1CS::prover_init(&index_pk.index, c)?;
        let public_input = prover_init_state.public_input();
        let mut fs_rng = FS::initialize(
            &to_bytes![&Self::PROTOCOL_NAME, &index_pk.index_vk, &public_input].unwrap(),
        );

        // --------------------------------------------------------------------
        // First round

        let (prover_first_msg, prover_first_oracles, prover_state) =
            AHPForR1CS::prover_first_round(prover_init_state, zk_rng)?;

        let first_round_comm_time = start_timer!(|| "Committing to first round polys");
        let (first_comms, first_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_first_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(first_round_comm_time);

        fs_rng.absorb(&to_bytes![first_comms, prover_first_msg].unwrap());

        let (verifier_first_msg, verifier_state) =
            AHPForR1CS::verifier_first_round(index_pk.index_vk.index_info, &mut fs_rng)?;
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Second round

        let (prover_second_msg, prover_second_oracles, prover_state) =
            AHPForR1CS::prover_second_round(&verifier_first_msg, prover_state, zk_rng);

        let second_round_comm_time = start_timer!(|| "Committing to second round polys");
        let (second_comms, second_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_second_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(second_round_comm_time);

        fs_rng.absorb(&to_bytes![second_comms, prover_second_msg].unwrap());

        let (verifier_second_msg, verifier_state) =
            AHPForR1CS::verifier_second_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Third round
        let (prover_third_msg, prover_third_oracles) =
            AHPForR1CS::prover_third_round(&verifier_second_msg, prover_state, zk_rng)?;

        let third_round_comm_time = start_timer!(|| "Committing to third round polys");
        let (third_comms, third_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_third_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(third_round_comm_time);

        fs_rng.absorb(&to_bytes![third_comms, prover_third_msg].unwrap());

        let verifier_state = AHPForR1CS::verifier_third_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // Gather prover polynomials in one vector.
        let polynomials: Vec<_> = index_pk
            .index
            .iter()
            .chain(prover_first_oracles.iter())
            .chain(prover_second_oracles.iter())
            .chain(prover_third_oracles.iter())
            .collect();

        // Gather commitments in one vector.
        #[rustfmt::skip]
        let commitments = vec![
            first_comms.iter().map(|p| p.commitment().clone()).collect(),
            second_comms.iter().map(|p| p.commitment().clone()).collect(),
            third_comms.iter().map(|p| p.commitment().clone()).collect(),
        ];
        let labeled_comms: Vec<_> = index_pk
            .index_vk
            .iter()
            .cloned()
            .zip(&AHPForR1CS::<F>::INDEXER_POLYNOMIALS)
            .map(|(c, l)| LabeledCommitment::new(l.to_string(), c, None))
            .chain(first_comms.iter().cloned())
            .chain(second_comms.iter().cloned())
            .chain(third_comms.iter().cloned())
            .collect();

        // Gather commitment randomness together.
        let comm_rands: Vec<PC::Randomness> = index_pk
            .index_comm_rands
            .clone()
            .into_iter()
            .chain(first_comm_rands)
            .chain(second_comm_rands)
            .chain(third_comm_rands)
            .collect();

        // Compute the AHP verifier's query set.
        let (query_set, verifier_state) =
            AHPForR1CS::verifier_query_set(verifier_state, &mut fs_rng);
        let lc_s = AHPForR1CS::construct_linear_combinations(
            &public_input,
            &polynomials,
            &verifier_state,
        )?;

        let eval_time = start_timer!(|| "Evaluating linear combinations over query set");
        let mut evaluations = Vec::new();
        for (label, (_, point)) in &query_set {
            let lc = lc_s
                .iter()
                .find(|lc| &lc.label == label)
                .ok_or(ahp::Error::MissingEval(label.to_string()))?;
            let eval = polynomials.get_lc_eval(&lc, *point)?;
            if !AHPForR1CS::<F>::LC_WITH_ZERO_EVAL.contains(&lc.label.as_ref()) {
                evaluations.push((label.to_string(), eval));
            }
        }

        evaluations.sort_by(|a, b| a.0.cmp(&b.0));
        let evaluations = evaluations.into_iter().map(|x| x.1).collect::<Vec<F>>();
        end_timer!(eval_time);

        fs_rng.absorb(&evaluations);
        let opening_challenge: F = u128::rand(&mut fs_rng).into();

        let pc_proof = PC::open_combinations(
            &index_pk.committer_key,
            &lc_s,
            polynomials,
            &labeled_comms,
            &query_set,
            opening_challenge,
            &comm_rands,
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;

        // Gather prover messages together.
        let prover_messages = vec![prover_first_msg, prover_second_msg, prover_third_msg];

        let proof = Proof::new(commitments, evaluations, prover_messages, pc_proof);
        proof.print_size_info();
        end_timer!(prover_time);
        Ok(proof)
    }

    /// Create a zkSNARK asserting that the constraint system is satisfied where circuit is private.
    pub fn index_private_prove<C: ConstraintSynthesizer<F>, R: RngCore>(
        index_pk: &IndexPrivateProverKey<F, PC>,
        c: C,
        zk_rng: &mut R,
    ) -> Result<IndexPrivateProof<F, PC>, Error<PC::Error>> {
        let prover_time = start_timer!(|| "Marlin::Prover");
        // Add check that c is in the correct mode.

        let prover_init_state = AHPForR1CS::prover_init(&index_pk.index, c)?;
        let public_input = prover_init_state.public_input();
        let mut fs_rng = FS::initialize(
            &to_bytes![&Self::PROTOCOL_NAME, &index_pk.index_private_vk, &public_input].unwrap(),
        );

        // --------------------------------------------------------------------
        // First round

        let (prover_first_msg, prover_first_oracles, prover_state) =
            AHPForR1CS::prover_index_private_first_round(prover_init_state, zk_rng)?;

        let first_round_comm_time = start_timer!(|| "Committing to first round polys");
        let (first_comms, first_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_first_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(first_round_comm_time);

        fs_rng.absorb(&to_bytes![first_comms, prover_first_msg].unwrap());

        let (verifier_first_msg, verifier_state) =
            AHPForR1CS::verifier_first_round(index_pk.index_private_vk.index_info, &mut fs_rng)?;
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Second round

        let (prover_second_msg, prover_second_oracles, prover_state) =
            AHPForR1CS::prover_second_round(&verifier_first_msg, prover_state, zk_rng);

        let domain_k = prover_state.domain_k.clone();

        let second_round_comm_time = start_timer!(|| "Committing to second round polys");
        let (second_comms, second_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_second_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(second_round_comm_time);

        fs_rng.absorb(&to_bytes![second_comms, prover_second_msg].unwrap());

        let (verifier_second_msg, verifier_state) =
            AHPForR1CS::verifier_second_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Third round
        let (prover_third_msg, prover_third_oracles) =
            AHPForR1CS::prover_index_private_third_round(
                &verifier_second_msg,
                prover_state,
                zk_rng,
            )?;

        let third_round_comm_time = start_timer!(|| "Committing to third round polys");
        let (third_comms, third_comm_rands) = PC::commit(
            &index_pk.committer_key,
            prover_third_oracles.iter(),
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;
        end_timer!(third_round_comm_time);

        fs_rng.absorb(&to_bytes![third_comms, prover_third_msg].unwrap());

        let verifier_state = AHPForR1CS::verifier_third_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // Gather prover polynomials in one vector.
        let polynomials: Vec<_> = 
        // index_pk
        //     .index
        //     .iter()
            // .chain(prover_first_oracles.iter())
            prover_first_oracles.iter()
            .chain(prover_second_oracles.iter())
            .chain(prover_third_oracles.iter())
            .collect();

        // Gather commitments in one vector.
        #[rustfmt::skip]
            let commitments = vec![
                first_comms.iter().map(|p| p.commitment().clone()).collect(),
                second_comms.iter().map(|p| p.commitment().clone()).collect(),
                third_comms.iter().map(|p| p.commitment().clone()).collect(),
            ];

        let labeled_comms: Vec<_> = 
            // index_pk
            // .index_vk
            // .iter()
            // .cloned()
            // .zip(&AHPForR1CS::<F>::INDEXER_POLYNOMIALS)
            // .map(|(c, l)| LabeledCommitment::new(l.to_string(), c, None))
            // .chain(first_comms.iter().cloned())
            first_comms.iter().cloned()
            .chain(second_comms.iter().cloned())
            .chain(third_comms.iter().cloned())
            .collect();

        // Gather commitment randomness together.
        let comm_rands: Vec<PC::Randomness> = 
            // index_pk
            // .index_comm_rands
            // .clone()
            // .into_iter()
            // .chain(first_comm_rands)
            first_comm_rands.into_iter()
            .chain(second_comm_rands)
            .chain(third_comm_rands)
            .collect();

        // Compute the AHP verifier's query set.
        let (query_set, verifier_state) =
            AHPForR1CS::index_private_verifier_query_set(verifier_state, &mut fs_rng);
        let lc_s = AHPForR1CS::construct_linear_combinations_for_index_private(
            &public_input,
            &polynomials,
            &verifier_state,
        )?;

        let eval_time = start_timer!(|| "Evaluating linear combinations over query set");
        let mut evaluations = Vec::new();
        for (label, (_, point)) in &query_set {
            let lc = lc_s
                .iter()
                .find(|lc| &lc.label == label)
                .ok_or(ahp::Error::MissingEval(label.to_string()))?;
            let eval = polynomials.get_lc_eval(&lc, *point)?;
            if !AHPForR1CS::<F>::INDEX_PRIVATE_LC_WITH_ZERO_EVAL.contains(&lc.label.as_ref()) {
                evaluations.push((label.to_string(), eval));
            }
        }

        evaluations.sort_by(|a, b| a.0.cmp(&b.0));
        println!("{:?}", evaluations);
        let evaluations = evaluations.into_iter().map(|x| x.1).collect::<Vec<F>>();
        end_timer!(eval_time);

        fs_rng.absorb(&evaluations);
        let opening_challenge: F = u128::rand(&mut fs_rng).into();

        let concrete_oracles = [
            index_pk.index.a_arith.row.clone(),
            index_pk.index.a_arith.col.clone(),
            index_pk.index.a_arith.val.clone(),

            index_pk.index.b_arith.row.clone(),
            index_pk.index.b_arith.col.clone(),
            index_pk.index.b_arith.val.clone(),

            index_pk.index.c_arith.row.clone(),
            index_pk.index.c_arith.col.clone(),
            index_pk.index.c_arith.val.clone(),
            prover_third_oracles.f.clone() // f
        ];

        let rational_sumcheck_vo = GenericShiftingVO::new(
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &vec![F::one(); 10],
            rational_sumcheck_oracle!(verifier_first_msg, verifier_second_msg, domain_k)
        );

        let labels = vec!["a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val"];
        let mut rational_sumcheck_commitments = index_pk.index_private_vk.polys.iter().zip(labels.iter())
            .map(|(commitment, &label)| 
                LabeledCommitment::new(label.into(), commitment.clone(), None)
        ).collect::<Vec<_>>();

        rational_sumcheck_commitments.push(LabeledCommitment::new("f".into(), third_comms[0].commitment().clone(), None));

        let empty_rands = vec![PC::Randomness::empty(); concrete_oracles.len()];
            
        let rational_sumcheck_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &rational_sumcheck_commitments,
            empty_rands.as_slice(),
            &rational_sumcheck_vo,
            &vec![F::one(); concrete_oracles.len()],
            &domain_k,
            &index_pk.committer_key,
            zk_rng,
        )?;

        let pc_proof = PC::open_combinations(
            &index_pk.committer_key,
            &lc_s,
            polynomials.clone(),
            &labeled_comms,
            &query_set,
            opening_challenge,
            &comm_rands,
            Some(zk_rng),
        )
        .map_err(Error::from_pc_err)?;

        // Gather prover messages together.
        let prover_messages = vec![prover_first_msg, prover_second_msg, prover_third_msg];

        let proof = IndexPrivateProof::new(commitments, evaluations, prover_messages, pc_proof, rational_sumcheck_proof);
        // proof.print_size_info();
        end_timer!(prover_time);
        Ok(proof)
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied.
    pub fn verify<R: RngCore>(
        index_vk: &IndexVerifierKey<F, PC>,
        public_input: &[F],
        proof: &Proof<F, PC>,
        rng: &mut R,
    ) -> Result<bool, Error<PC::Error>> {
        let verifier_time = start_timer!(|| "Marlin::Verify");

        let public_input = {
            let domain_x = GeneralEvaluationDomain::<F>::new(public_input.len() + 1).unwrap();

            let mut unpadded_input = public_input.to_vec();
            unpadded_input.resize(
                core::cmp::max(public_input.len(), domain_x.size() - 1),
                F::zero(),
            );

            unpadded_input
        };

        let mut fs_rng =
            FS::initialize(&to_bytes![&Self::PROTOCOL_NAME, &index_vk, &public_input].unwrap());

        // --------------------------------------------------------------------
        // First round

        let first_comms = &proof.commitments[0];
        fs_rng.absorb(&to_bytes![first_comms, proof.prover_messages[0]].unwrap());

        let (_, verifier_state) =
            AHPForR1CS::verifier_first_round(index_vk.index_info, &mut fs_rng)?;
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Second round
        let second_comms = &proof.commitments[1];
        fs_rng.absorb(&to_bytes![second_comms, proof.prover_messages[1]].unwrap());

        let (_, verifier_state) = AHPForR1CS::verifier_second_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Third round
        let third_comms = &proof.commitments[2];
        fs_rng.absorb(&to_bytes![third_comms, proof.prover_messages[2]].unwrap());

        let verifier_state = AHPForR1CS::verifier_third_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // Collect degree bounds for commitments. Indexed polynomials have *no*
        // degree bounds because we know the committed index polynomial has the
        // correct degree.
        let index_info = index_vk.index_info;
        let degree_bounds = vec![None; index_vk.index_comms.len()]
            .into_iter()
            .chain(AHPForR1CS::prover_first_round_degree_bounds(&index_info))
            .chain(AHPForR1CS::prover_second_round_degree_bounds(&index_info))
            .chain(AHPForR1CS::prover_third_round_degree_bounds(&index_info))
            .collect::<Vec<_>>();

        // Gather commitments in one vector.
        let commitments: Vec<_> = index_vk
            .iter()
            .chain(first_comms)
            .chain(second_comms)
            .chain(third_comms)
            .cloned()
            .zip(AHPForR1CS::<F>::polynomial_labels())
            .zip(degree_bounds)
            .map(|((c, l), d)| LabeledCommitment::new(l, c, d))
            .collect();

        let (query_set, verifier_state) =
            AHPForR1CS::verifier_query_set(verifier_state, &mut fs_rng);

        fs_rng.absorb(&proof.evaluations);
        let opening_challenge: F = u128::rand(&mut fs_rng).into();

        let mut evaluations = Evaluations::new();
        let mut evaluation_labels = Vec::new();
        for (poly_label, (_, point)) in query_set.iter().cloned() {
            if AHPForR1CS::<F>::LC_WITH_ZERO_EVAL.contains(&poly_label.as_ref()) {
                evaluations.insert((poly_label, point), F::zero());
            } else {
                evaluation_labels.push((poly_label, point));
            }
        }
        evaluation_labels.sort_by(|a, b| a.0.cmp(&b.0));
        for (q, eval) in evaluation_labels.into_iter().zip(&proof.evaluations) {
            evaluations.insert(q, *eval);
        }

        let lc_s = AHPForR1CS::construct_linear_combinations(
            &public_input,
            &evaluations,
            &verifier_state,
        )?;

        let evaluations_are_correct = PC::check_combinations(
            &index_vk.verifier_key,
            &lc_s,
            &commitments,
            &query_set,
            &evaluations,
            &proof.pc_proof,
            opening_challenge,
            rng,
        )
        .map_err(Error::from_pc_err)?;

        if !evaluations_are_correct {
            eprintln!("PC::Check failed");
        }
        end_timer!(verifier_time, || format!(
            " PC::Check for AHP Verifier linear equations: {}",
            evaluations_are_correct
        ));
        Ok(evaluations_are_correct)
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied where circuit is private.
    pub fn verify_index_private<R: RngCore>(
        index_vk: &IndexPrivateVerifierKey<F, PC>,
        public_input: &[F],
        proof: &IndexPrivateProof<F, PC>,
        rng: &mut R,
    ) -> Result<bool, Error<PC::Error>> {
        let verifier_time = start_timer!(|| "Marlin::Verify");

        let public_input = {
            let domain_x = GeneralEvaluationDomain::<F>::new(public_input.len() + 1).unwrap();

            let mut unpadded_input = public_input.to_vec();
            unpadded_input.resize(
                core::cmp::max(public_input.len(), domain_x.size() - 1),
                F::zero(),
            );

            unpadded_input
        };

        let mut fs_rng =
            FS::initialize(&to_bytes![&Self::PROTOCOL_NAME, &index_vk, &public_input].unwrap());

        // --------------------------------------------------------------------
        // First round

        let first_comms = &proof.commitments[0];
        fs_rng.absorb(&to_bytes![first_comms, proof.prover_messages[0]].unwrap());

        let (verifier_first_msg, verifier_state) =
            AHPForR1CS::verifier_first_round(index_vk.index_info, &mut fs_rng)?;
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Second round
        let second_comms = &proof.commitments[1];
        fs_rng.absorb(&to_bytes![second_comms, proof.prover_messages[1]].unwrap());

        let (verifier_second_msg, verifier_state) = AHPForR1CS::verifier_second_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Third round
        let third_comms = &proof.commitments[2];
        fs_rng.absorb(&to_bytes![third_comms, proof.prover_messages[2]].unwrap());

        let verifier_state = AHPForR1CS::verifier_third_round(verifier_state, &mut fs_rng);
        // --------------------------------------------------------------------

        // Collect degree bounds for commitments. Indexed polynomials have *no*
        // degree bounds because we know the committed index polynomial has the
        // correct degree.
        let index_info = index_vk.index_info;
        let degree_bounds = 
            // vec![None; index_vk.index_comms.len()]
            // .into_iter()
            // .chain(AHPForR1CS::prover_index_private_first_round_degree_bounds(&index_info))
            AHPForR1CS::prover_index_private_first_round_degree_bounds(&index_info)
            .chain(AHPForR1CS::prover_second_round_degree_bounds(&index_info))
            .chain(AHPForR1CS::index_private_prover_third_round_degree_bounds(
                &index_info,
            ))
            .collect::<Vec<_>>();

        // Gather commitments in one vector.
        let commitments: Vec<_> = 
        // index_vk
        //     .iter()
        //     .chain(first_comms)
            first_comms.into_iter()
            .chain(second_comms)
            .chain(third_comms)
            .cloned()
            .zip(AHPForR1CS::<F>::index_private_polynomial_labels())
            .zip(degree_bounds)
            .map(|((c, l), d)| LabeledCommitment::new(l, c, d))
            .collect();

        let (query_set, verifier_state) =
            AHPForR1CS::index_private_verifier_query_set(verifier_state, &mut fs_rng);

        fs_rng.absorb(&proof.evaluations);
        let opening_challenge: F = u128::rand(&mut fs_rng).into();

        let mut evaluations = Evaluations::new();
        let mut evaluation_labels = Vec::new();
        for (poly_label, (_, point)) in query_set.iter().cloned() {
            if AHPForR1CS::<F>::INDEX_PRIVATE_LC_WITH_ZERO_EVAL.contains(&poly_label.as_ref()) {
                evaluations.insert((poly_label, point), F::zero());
            } else {
                evaluation_labels.push((poly_label, point));
            }
        }
        evaluation_labels.sort_by(|a, b| a.0.cmp(&b.0));
        for (q, eval) in evaluation_labels.into_iter().zip(&proof.evaluations) {
            evaluations.insert(q, *eval);
        }

        let lc_s = AHPForR1CS::construct_linear_combinations_for_index_private(
            &public_input,
            &evaluations,
            &verifier_state,
        )?;

        let evaluations_are_correct = PC::check_combinations(
            &index_vk.verifier_key,
            &lc_s,
            &commitments,
            &query_set,
            &evaluations,
            &proof.pc_proof,
            opening_challenge,
            rng,
        )
        .map_err(Error::from_pc_err)?;

        let domain_k = verifier_state.domain_k.clone();

        let rational_sumcheck_vo = GenericShiftingVO::new(
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &vec![F::one(); 10],
            rational_sumcheck_oracle!(verifier_first_msg, verifier_second_msg, domain_k)
        );

        let labels = vec!["a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val"];
        let mut rational_sumcheck_commitments = index_vk.polys.iter()
            .zip(labels.iter())
            .map(|(commitment, &label)| 
                LabeledCommitment::new(label.into(), commitment.clone(), None)
        ).collect::<Vec<_>>();

        rational_sumcheck_commitments.push(LabeledCommitment::new("f".into(), third_comms[0].clone(), None)); // TODO f should also be bounded with |K| -1

        let is_valid = ZeroOverK::<F, PC, FS>::verify(
            &proof.rational_sumcheck_zero_over_k_proof,
            &rational_sumcheck_commitments,
            &rational_sumcheck_vo,
            &domain_k,
            &vec![F::one(); 10],
            &index_vk.verifier_key,
            rng
        );

        println!("{:?}", is_valid);

        if !evaluations_are_correct {
            eprintln!("PC::Check failed");
        }
        end_timer!(verifier_time, || format!(
            " PC::Check for AHP Verifier linear equations: {}",
            evaluations_are_correct
        ));

        Ok(evaluations_are_correct)
    }
}

/// A closure to be used for the Marlin rational sumcheck virtual oracle. The expected terms are:
/// terms[0]: ignore, this by default becomes X
/// terms[1]: a_row(X)
/// terms[2]: a_col(X)
/// terms[3]: a_val(X)
/// terms[4]: b_row(X)
/// terms[5]: b_col(X)
/// terms[6]: b_val(X)
/// terms[7]: c_row(X)
/// terms[8]: c_col(X)
/// terms[9]: c_val(X)
/// terms[10]: f(X)
#[macro_export]
macro_rules! rational_sumcheck_oracle {
    ($verifier_first_msg:expr, $verifier_second_msg:expr, $domain_k:expr) => {
        |terms: &[VOTerm<F>]| {
            // define consts
            let alpha = vo_constant!($verifier_first_msg.alpha);
            let beta = vo_constant!($verifier_second_msg.beta);

            let vh_alpha = $domain_k.evaluate_vanishing_polynomial($verifier_first_msg.alpha);
            let vh_beta = $domain_k.evaluate_vanishing_polynomial($verifier_second_msg.beta);
            let v_H_alpha_v_H_beta = vo_constant!(vh_alpha * vh_beta);

            let eta_a_times_v_H_alpha_v_H_beta = vo_constant!($verifier_first_msg.eta_a) * v_H_alpha_v_H_beta;
            let eta_b_times_v_H_alpha_v_H_beta = vo_constant!($verifier_first_msg.eta_b) * v_H_alpha_v_H_beta;
            let eta_c_times_v_H_alpha_v_H_beta = vo_constant!($verifier_first_msg.eta_c) * v_H_alpha_v_H_beta;

            let alpha_beta = alpha * beta;

            // define terms
            let a_row = terms[1];
            let a_col = terms[2];
            let a_val = terms[3];
            let b_row = terms[4];
            let b_col = terms[5];
            let b_val = terms[6];
            let c_row = terms[7];
            let c_col = terms[8];
            let c_val = terms[9];
            let f = terms[10];

            // begin logic
            let a_denom = alpha_beta - beta*a_col - alpha*a_row + a_col*a_row;
            let b_denom = alpha_beta - beta*b_col - alpha*b_row + b_col*b_row;
            let c_denom = alpha_beta - beta*c_col - alpha*c_row + c_col*c_row;

            let b_poly = a_denom * b_denom * c_denom;

            let a_part_nom = eta_a_times_v_H_alpha_v_H_beta * a_val;
            let b_part_nom = eta_b_times_v_H_alpha_v_H_beta * b_val;
            let c_part_nom = eta_c_times_v_H_alpha_v_H_beta * c_val;

            let a_poly = {
                let summand_0 = a_part_nom * b_denom * c_denom;
                let summand_1 = b_part_nom * b_denom * c_denom;
                let summand_2 = c_part_nom * b_denom * c_denom;

                summand_0 + summand_1 + summand_2
            };

            a_poly - b_poly - f
        }
    };
}