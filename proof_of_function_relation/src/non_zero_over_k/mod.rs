use crate::non_zero_over_k::{piop::PIOPforNonZeroOverK, proof::Proof};
use crate::{
    commitment::AdditivelyHomomorphicPCS,
    error::{to_pc_error, Error},
    virtual_oracle::inverse_check_oracle::InverseCheckOracle,
    zero_over_k::ZeroOverK,
};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use std::marker::PhantomData;

pub mod piop;
pub mod proof;
mod tests;

pub struct NonZeroOverK<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, D: Digest> NonZeroOverK<F, PC, D> {
    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        domain: &GeneralEvaluationDomain<F>,
        f: &LabeledPolynomial<F, DensePolynomial<F>>,
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
        let g = prover_first_oracles.g;

        let inverse_check_oracle = InverseCheckOracle::new();
        let concrete_oracles = [f.clone(), g];
        let alphas = vec![F::one(), F::one()];
        let (commitments, rands) =
            PC::commit(ck, &concrete_oracles, None).map_err(to_pc_error::<F, PC>)?;

        let zero_over_k_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &commitments,
            &rands,
            f.degree_bound(),
            &inverse_check_oracle,
            &alphas,
            &domain,
            ck,
            rng,
        )?;

        let proof = Proof {
            g_commit: commitments[1].commitment().clone(),
            zero_over_k_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        domain: &GeneralEvaluationDomain<F>,
        f_commit: LabeledCommitment<PC::Commitment>,
        proof: Proof<F, PC>,
    ) -> Result<(), Error> {
        //TODO check g bound
        let g_commit = LabeledCommitment::new(String::from("g"), proof.g_commit.clone(), None);

        let concrete_oracles_commitments = [f_commit.clone(), g_commit];
        let zero_over_k_vo = InverseCheckOracle::new();
        let alphas = vec![F::one(), F::one()];

        ZeroOverK::<F, PC, D>::verify(
            proof.zero_over_k_proof,
            &concrete_oracles_commitments,
            f_commit.degree_bound(),
            &zero_over_k_vo,
            &domain,
            &alphas,
            &vk,
        )
    }
}
