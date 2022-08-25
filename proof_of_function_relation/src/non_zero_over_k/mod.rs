use crate::error::{to_pc_error, Error};
use crate::non_zero_over_k::{piop::PIOPforNonZeroOverK, proof::Proof};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use std::marker::PhantomData;
use zero_over_k::{
    virtual_oracle::generic_shifting_vo::{presets, GenericShiftingVO},
    zero_over_k::ZeroOverK,
};

pub mod piop;
pub mod proof;
mod tests;

pub struct NonZeroOverK<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _fs_rng: PhantomData<FS>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> NonZeroOverK<F, PC, FS> {
    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        domain: &GeneralEvaluationDomain<F>,
        f: &LabeledPolynomial<F, DensePolynomial<F>>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        f_rand: &PC::Randomness,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        //-----------------------------------------------
        // INIT PROVER
        let prover_initial_state = PIOPforNonZeroOverK::prover_init(domain, f)?;

        //-----------------------------------------------
        // FIRST ROUND
        let (_, prover_first_oracles, _prover_state) =
            PIOPforNonZeroOverK::prover_first_round(prover_initial_state, rng)?;

        //-----------------------------------------------
        // RUN SUBPROTOCOLS
        let (commitments, rands) = PC::commit(ck, &[prover_first_oracles.g.clone()], Some(rng))
            .map_err(to_pc_error::<F, PC>)?;

        let concrete_oracles = [f.clone(), prover_first_oracles.g.clone()];

        let alphas = vec![F::one(), F::one()];
        let inverse_check_oracle =
            GenericShiftingVO::new(&vec![0, 1], &alphas, presets::inverse_check)?;

        let zero_over_k_proof = ZeroOverK::<F, PC, FS>::prove(
            &concrete_oracles,
            &[f_commit.clone(), commitments[0].clone()],
            &[f_rand.clone(), rands[0].clone()],
            f.degree_bound(),
            &inverse_check_oracle,
            &domain,
            ck,
            rng,
        )?;

        let proof = Proof {
            g_commit: commitments[0].commitment().clone(),
            zero_over_k_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        domain: &GeneralEvaluationDomain<F>,
        f_commit: PC::Commitment,
        enforced_degree_bound: Option<usize>,
        proof: Proof<F, PC>,
    ) -> Result<(), Error> {
        let bounded_f_commit =
            LabeledCommitment::new(String::from("f"), f_commit, enforced_degree_bound);
        let g_commit = LabeledCommitment::new(
            String::from("g"),
            proof.g_commit.clone(),
            enforced_degree_bound,
        );

        let concrete_oracles_commitments = [bounded_f_commit.clone(), g_commit];
        let alphas = vec![F::one(), F::one()];
        let inverse_check_oracle =
            GenericShiftingVO::new(&vec![0, 1], &alphas, presets::inverse_check)?;

        ZeroOverK::<F, PC, FS>::verify(
            proof.zero_over_k_proof,
            &concrete_oracles_commitments,
            enforced_degree_bound,
            &inverse_check_oracle,
            &domain,
            &vk,
        )
        .map_err(Error::from)
    }
}
