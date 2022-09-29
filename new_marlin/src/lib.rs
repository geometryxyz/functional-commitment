pub mod ahp;
pub mod data_structures;
pub mod error;

use crate::ahp::{AHPForR1CS, EvaluationsProvider};
use ::zero_over_k::{
    virtual_oracle::generic_shifting_vo::{vo_term::VOTerm, GenericShiftingVO},
    vo_constant,
    zero_over_k::ZeroOverK,
};
use ahp::{
    constraint_systems::arithmetize_matrix,
    indexer::{Index, Matrix},
};
use ark_ff::{to_bytes, PrimeField, UniformRand};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{Evaluations, PCRandomness};
use ark_poly_commit::{LabeledCommitment, PCUniversalParams};
use ark_relations::r1cs::ConstraintSynthesizer;
use ark_std::rand::RngCore;
use data_structures::{Proof, ProverKey, UniversalSRS, VerifierKey};
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use new_ac_compiler::R1CSfIndex;

#[macro_use]
extern crate ark_std;

use ark_std::{
    collections::BTreeMap,
    format,
    marker::PhantomData,
    string::{String, ToString},
    vec,
    vec::Vec,
};

pub use error::*;

// pub use ahp::constraint_systems::arithmetize_matrix; //TODO: for

/// The compiled argument system.
pub struct Marlin<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng>(
    PhantomData<F>,
    PhantomData<PC>,
    PhantomData<FS>,
);

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> Marlin<F, PC, FS> {
    pub const PROTOCOL_NAME: &'static [u8] = b"INDEX-PRIVATE-MARLIN";

    /// Generate the universal prover and verifier keys for the
    /// argument system.
    pub fn universal_setup<R: RngCore>(
        index_info: &R1CSfIndex,
        rng: &mut R,
    ) -> Result<UniversalSRS<F, PC>, Error<PC::Error>> {
        let max_degree = AHPForR1CS::<F>::max_degree(index_info)?;
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

    pub fn index(
        srs: &UniversalSRS<F, PC>,
        index_info: &R1CSfIndex,
        a: Matrix<F>,
        b: Matrix<F>,
        c: Matrix<F>,
    ) -> Result<(ProverKey<F, PC>, VerifierKey<F, PC>, Vec<PC::Randomness>), Error<PC::Error>> {
        if !index_info.check_domains_sizes::<F>() {
            return Err(Error::DomainHLargerThanDomainK);
        }
        let domain_k = GeneralEvaluationDomain::<F>::new(index_info.number_of_non_zero_entries)
            .ok_or(Error::DomainTooLarge)?;
        let domain_h = GeneralEvaluationDomain::<F>::new(index_info.number_of_constraints)
            .ok_or(Error::DomainTooLarge)?;

        let supported_hiding_bound = 1;

        let max_degree = AHPForR1CS::<F>::max_degree(index_info)?;
        let degree_bounds = AHPForR1CS::<F>::get_degree_bounds(index_info);

        let (committer_key, verifier_key) = PC::trim(
            &srs,
            max_degree,
            supported_hiding_bound,
            Some(&degree_bounds),
        )
        .map_err(Error::from_pc_err)?;

        let a_arith = arithmetize_matrix(&a, domain_k, domain_h, "a", false);
        let b_arith = arithmetize_matrix(&b, domain_k, domain_h, "b", false);
        let c_arith = arithmetize_matrix(&c, domain_k, domain_h, "c", true);

        let polys = [
            a_arith.row.clone(),
            a_arith.col.clone(),
            a_arith.val.clone(),
            b_arith.row.clone(),
            b_arith.col.clone(),
            b_arith.val.clone(),
            c_arith.row.clone(),
            c_arith.col.clone(),
            c_arith.val.clone(),
        ];

        let (matrix_poly_commits, matrix_poly_rands): (_, _) =
            PC::commit(&committer_key, &polys, None).map_err(Error::from_pc_err)?;

        let matrix_poly_commits = matrix_poly_commits
            .iter()
            .map(|c| c.commitment().clone())
            .collect::<Vec<_>>();

        let vk: VerifierKey<F, PC> = VerifierKey {
            commits: matrix_poly_commits,
            verifier_key,
            index_info: index_info.clone(),
        };

        let index = Index::<F> {
            index_info: index_info.clone(),
            a_arith,
            b_arith,
            c_arith,

            a,
            b,
            c,
        };

        let pk = ProverKey {
            index: index.clone(),
            vk: vk.clone(),
            committer_key: committer_key.clone(),
        };

        Ok((pk, vk, matrix_poly_rands))
        // arith matrices
    }

    pub fn prove<R: RngCore>(
        index_pk: &ProverKey<F, PC>,
        assignment: Vec<F>,
        zk_rng: &mut R,
    ) -> Result<Proof<F, PC>, Error<PC::Error>> {
        let prover_time = start_timer!(|| "Marlin::Prover");

        let prover_init_state = AHPForR1CS::prover_init(&index_pk.index, assignment)?;
        let public_input = prover_init_state.public_input();

        let mut fs_rng =
            FS::initialize(&to_bytes![&Self::PROTOCOL_NAME, &index_pk.vk, &public_input].unwrap());

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
            AHPForR1CS::verifier_first_round(&index_pk.vk.index_info, &mut fs_rng)?;
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
            AHPForR1CS::prover_index_private_third_round(&verifier_second_msg, prover_state)?;

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
        let polynomials: Vec<_> = prover_first_oracles
            .iter()
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

        let labeled_comms: Vec<_> = first_comms
            .iter()
            .cloned()
            .chain(second_comms.iter().cloned())
            .chain(third_comms.iter().cloned())
            .collect();

        // Gather commitment randomness together.
        let comm_rands: Vec<PC::Randomness> = first_comm_rands
            .into_iter()
            .chain(second_comm_rands)
            .chain(third_comm_rands)
            .collect();

        // Compute the AHP verifier's query set.
        let (query_set, verifier_state) =
            AHPForR1CS::verifier_query_set(verifier_state, &mut fs_rng);
        let lc_s = AHPForR1CS::construct_linear_combinations(
            // &public_input,
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
            prover_third_oracles.f.clone(), // f
        ];

        let rational_sumcheck_vo = GenericShiftingVO::new(
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &vec![F::one(); 10],
            rational_sumcheck_oracle!(verifier_first_msg, verifier_second_msg, domain_k),
        )?;

        let labels = vec![
            "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
        ];
        let mut rational_sumcheck_commitments = index_pk
            .vk
            .commits
            .iter()
            .zip(labels.iter())
            .map(|(commitment, &label)| {
                LabeledCommitment::new(label.into(), commitment.clone(), None)
            })
            .collect::<Vec<_>>();

        rational_sumcheck_commitments.push(LabeledCommitment::new(
            "f".into(),
            third_comms[0].commitment().clone(),
            None,
        ));

        let empty_rands = vec![PC::Randomness::empty(); concrete_oracles.len()];

        let rational_sumcheck_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &rational_sumcheck_commitments,
            empty_rands.as_slice(),
            None,
            &rational_sumcheck_vo,
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

        let proof = Proof::new(
            commitments,
            evaluations,
            prover_messages,
            pc_proof,
            rational_sumcheck_proof,
        );
        // proof.print_size_info();
        end_timer!(prover_time);
        Ok(proof)
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied where circuit is private.
    pub fn verify<R: RngCore>(
        vk: &data_structures::VerifierKey<F, PC>,
        public_input: &Vec<F>,
        proof: Proof<F, PC>,
        rng: &mut R,
    ) -> Result<bool, Error<PC::Error>> {
        let verifier_time = start_timer!(|| "Marlin::Verify");

        // let public_input =
        //     let domain_x = GeneralEvaluationDomain::<F>::new(public_input.len()).unwrap();

        //     let mut unpadded_input = public_input.to_vec();
        //     unpadded_input.resize(
        //         core::cmp::max(public_input.len(), domain_x.size() - 1),
        //         F::zero(),
        //     );

        //     unpadded_input
        // };

        let mut fs_rng =
            FS::initialize(&to_bytes![&Self::PROTOCOL_NAME, &vk, &public_input].unwrap());

        // --------------------------------------------------------------------
        // First round

        let first_comms = &proof.commitments[0];
        fs_rng.absorb(&to_bytes![first_comms, proof.prover_messages[0]].unwrap());

        let (verifier_first_msg, verifier_state) =
            AHPForR1CS::verifier_first_round(&vk.index_info, &mut fs_rng)?;
        // --------------------------------------------------------------------

        // --------------------------------------------------------------------
        // Second round
        let second_comms = &proof.commitments[1];
        fs_rng.absorb(&to_bytes![second_comms, proof.prover_messages[1]].unwrap());

        let (verifier_second_msg, verifier_state) =
            AHPForR1CS::verifier_second_round(verifier_state, &mut fs_rng);
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
        let index_info = vk.index_info.clone();
        let degree_bounds = AHPForR1CS::<F>::prover_first_round_degree_bounds(&index_info)
            .chain(AHPForR1CS::<F>::prover_second_round_degree_bounds(
                &index_info,
            ))
            .chain(AHPForR1CS::<F>::prover_third_round_degree_bounds(
                &index_info,
            ))
            .collect::<Vec<_>>();

        // Gather commitments in one vector.
        let commitments: Vec<_> = first_comms
            .into_iter()
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

        let lc_s = AHPForR1CS::construct_linear_combinations(&evaluations, &verifier_state)?;

        let evaluations_are_correct = PC::check_combinations(
            &vk.verifier_key,
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
            rational_sumcheck_oracle!(verifier_first_msg, verifier_second_msg, domain_k),
        )?;

        let labels = vec![
            "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
        ];
        let mut rational_sumcheck_commitments = vk
            .commits
            .iter()
            .zip(labels.iter())
            .map(|(commitment, &label)| {
                LabeledCommitment::new(label.into(), commitment.clone(), None)
            })
            .collect::<Vec<_>>();

        rational_sumcheck_commitments.push(LabeledCommitment::new(
            "f".into(),
            third_comms[0].clone(),
            None,
        ));

        ZeroOverK::<F, PC, FS>::verify(
            proof.rational_sumcheck_zero_over_k_proof,
            &rational_sumcheck_commitments,
            None,
            &rational_sumcheck_vo,
            &domain_k,
            &vk.verifier_key,
        )?;

        end_timer!(verifier_time, || format!(
            " PC::Check for AHP Verifier linear equations: {}",
            evaluations_are_correct
        ));

        Ok(evaluations_are_correct)
    }
}

// TODO: we must make this nicer! :)
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
            let v_H_alpha_v_H_beta = vh_alpha * vh_beta;

            let eta_a_times_v_H_alpha_v_H_beta =
                vo_constant!($verifier_first_msg.eta_a * v_H_alpha_v_H_beta);
            let eta_b_times_v_H_alpha_v_H_beta =
                vo_constant!($verifier_first_msg.eta_b * v_H_alpha_v_H_beta);
            let eta_c_times_v_H_alpha_v_H_beta =
                vo_constant!($verifier_first_msg.eta_c * v_H_alpha_v_H_beta);

            let alpha_beta = alpha.clone() * beta.clone();

            // define terms
            let a_row = terms[1].clone();
            let a_col = terms[2].clone();
            let a_val = terms[3].clone();
            let b_row = terms[4].clone();
            let b_col = terms[5].clone();
            let b_val = terms[6].clone();
            let c_row = terms[7].clone();
            let c_col = terms[8].clone();
            let c_val = terms[9].clone();
            let f = terms[10].clone();

            // begin logic
            let a_denom =
                alpha_beta.clone() - beta.clone() * a_col.clone() - alpha.clone() * a_row.clone()
                    + a_col.clone() * a_row.clone();
            let b_denom =
                alpha_beta.clone() - beta.clone() * b_col.clone() - alpha.clone() * b_row.clone()
                    + b_col.clone() * b_row.clone();
            let c_denom = alpha_beta - beta.clone() * c_col.clone() - alpha.clone() * c_row.clone()
                + c_col.clone() * c_row.clone();

            let b_poly = a_denom.clone() * b_denom.clone() * c_denom.clone();

            let a_part_nom = eta_a_times_v_H_alpha_v_H_beta * a_val;
            let b_part_nom = eta_b_times_v_H_alpha_v_H_beta * b_val;
            let c_part_nom = eta_c_times_v_H_alpha_v_H_beta * c_val;

            let a_poly = {
                let summand_0 = a_part_nom * b_denom.clone() * c_denom.clone();
                let summand_1 = b_part_nom * a_denom.clone() * c_denom.clone();
                let summand_2 = c_part_nom * a_denom.clone() * b_denom.clone();

                summand_0 + summand_1 + summand_2
            };

            a_poly - b_poly * f
        }
    };
}
