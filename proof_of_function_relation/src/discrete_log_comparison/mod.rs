use crate::{
    discrete_log_comparison::{piop::PIOPforDLComparison, proof::Proof},
    error::{to_pc_error, Error},
    geo_seq::GeoSeqTest,
    non_zero_over_k::NonZeroOverK,
    subset_over_k::SubsetOverK,
};
use ark_ff::{to_bytes, PrimeField, SquareRootField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use ark_std::marker::PhantomData;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use std::iter;
use zero_over_k::{
    virtual_oracle::generic_shifting_vo::{
        presets::{self, square_check},
        GenericShiftingVO,
    },
    zero_over_k::ZeroOverK,
};

use fiat_shamir_rng::FiatShamirRng;

pub mod piop;
pub mod proof;
mod tests;

pub struct DLComparison<
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
> {
    _field: PhantomData<F>,
    _polynomial_commitment_scheme: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F, PC, FS> DLComparison<F, PC, FS>
where
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"Discrete-log Comparison";

    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        f: &LabeledPolynomial<F, DensePolynomial<F>>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        f_rand: &PC::Randomness,
        g: &LabeledPolynomial<F, DensePolynomial<F>>,
        g_commit: &LabeledCommitment<PC::Commitment>,
        g_rand: &PC::Randomness,
        enforced_degree_bound: Option<usize>,
        fs_rng: &mut FS,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        let prover_initial_state =
            PIOPforDLComparison::prover_init(domain_k, domain_h, f, g, enforced_degree_bound)?;

        //------------------------------------------------------------------
        // First Round
        let (_, prover_first_oracles, prover_state) =
            PIOPforDLComparison::prover_first_round(prover_initial_state, rng)?;

        //------------------------------------------------------------------
        // Commit Phase

        let one_poly = DensePolynomial::from_coefficients_vec(vec![F::one()]);
        let one_poly =
            LabeledPolynomial::new(String::from("one"), one_poly, enforced_degree_bound, None);

        // commit to s, f_prime, g_prime, s_prime, h and the constant 1 polynomial
        // order of commitments is: s, f_prime, g_prime, s_prime, h, one
        let (commitments, rands) = PC::commit(
            ck,
            prover_first_oracles.iter().chain(iter::once(&one_poly)),
            Some(rng),
        )
        .map_err(to_pc_error::<F, PC>)?;

        let fs_bytes =
            &to_bytes![Self::PROTOCOL_NAME, commitments].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let alphas = [F::one(), F::one()];
        let square_check_vo = GenericShiftingVO::new(&[0, 1], &alphas, square_check)?;

        //------------------------------------------------------------------
        // Run sub-protocols

        // Step 4a: Zero over K for f = (f')^2
        let f_prime_square_proof = ZeroOverK::<F, PC, FS>::prove(
            &[f.clone(), prover_first_oracles.f_prime.clone()],
            &[f_commit.clone(), commitments[1].clone()], // f and f'
            &[f_rand.clone(), rands[1].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            &ck,
            rng,
        )?;

        // Step 4b: Zero over K for g = (g')^2
        let g_prime_square_proof = ZeroOverK::<F, PC, FS>::prove(
            &[g.clone(), prover_first_oracles.g_prime.clone()],
            &[g_commit.clone(), commitments[2].clone()], // g and g'
            &[g_rand.clone(), rands[2].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            &ck,
            rng,
        )?;

        // Step 4c: Zero over K for s = (s')^2
        let s_prime_square_proof = ZeroOverK::<F, PC, FS>::prove(
            &[
                prover_first_oracles.s.clone(),
                prover_first_oracles.s_prime.clone(),
            ],
            &[commitments[0].clone(), commitments[3].clone()], // s and s'
            &[rands[0].clone(), rands[3].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            &ck,
            rng,
        )?;

        // // SANITY CHECK
        // for element in domain_k.elements() {
        //     let eval = prover_first_oracles.f_prime.evaluate(&element)
        //         - prover_first_oracles.s_prime.evaluate(&element)
        //             * prover_first_oracles.g_prime.evaluate(&element);
        //     assert_eq!(eval, F::zero());
        // }

        // Step 4d: Zero over K for f' = (s')*(g')
        let product_check_vo =
            GenericShiftingVO::new(&[0, 1, 2], &vec![F::one(); 3], presets::abc_product_check)?;
        let f_prime_product_proof = ZeroOverK::<F, PC, FS>::prove(
            &[
                prover_first_oracles.f_prime.clone(),
                prover_first_oracles.s_prime.clone(),
                prover_first_oracles.g_prime.clone(),
            ],
            &[
                commitments[1].clone(),
                commitments[3].clone(),
                commitments[2].clone(),
            ], // f', s' and g'
            &[rands[1].clone(), rands[3].clone(), rands[2].clone()],
            enforced_degree_bound,
            &product_check_vo,
            &domain_k,
            &ck,
            rng,
        )?;

        // Step 5: Geometric sequence test on h
        let delta = prover_state
            .delta
            .expect("Delta should be computed in the prover's first round");
        let a_s = prover_state
            .a_s
            .expect("\'a\' values should be computed in the prover's first round");
        let c_s = prover_state
            .c_s
            .expect("\'c\' values should be computed in the prover's first round");

        let h_proof = GeoSeqTest::<F, PC, FS>::prove(
            &ck,
            delta,
            &prover_first_oracles.h,
            &commitments[4].clone(),
            &rands[4],
            &a_s,
            &c_s,
            &domain_k,
            rng,
        )?;

        // Step 6a: Subset over K for f'
        let f_prime_subset_proof = SubsetOverK::<F, PC, FS>::prove();

        // Step 6b: Subset over K for g'
        let g_prime_subset_proof = SubsetOverK::<F, PC, FS>::prove();

        // Step 6c: Subset over K for s'
        let s_prime_subset_proof = SubsetOverK::<F, PC, FS>::prove();

        // Step 7a: Non-zero over K for f′
        let nzk_f_prime_proof = NonZeroOverK::<F, PC, FS>::prove(
            ck,
            domain_k,
            &prover_first_oracles.f_prime,
            &commitments[1].clone(),
            &rands[1].clone(),
            rng,
        )?;

        // Step 7b: Non-zero over K for g′
        let nzk_g_prime_proof = NonZeroOverK::<F, PC, FS>::prove(
            ck,
            domain_k,
            &prover_first_oracles.g_prime,
            &commitments[2].clone(),
            &rands[2].clone(),
            rng,
        )?;

        // Step 7c: Non-zero over K for s′
        let nzk_s_prime_proof = NonZeroOverK::<F, PC, FS>::prove(
            ck,
            domain_k,
            &prover_first_oracles.s_prime,
            &commitments[3].clone(),
            &rands[3].clone(),
            rng,
        )?;

        // Step 7d: Non-zero over K for s(X) − 1
        // The verifier is expected to derive a commitment to the constant one polynomial on their own
        let s_minus_one = prover_first_oracles.s.polynomial() - one_poly.polynomial();
        let s_minus_one = LabeledPolynomial::new(
            String::from("s_minus_one"),
            s_minus_one,
            enforced_degree_bound,
            Some(1),
        );

        let (s_minus_one_commitment, s_minus_one_rand) = PC::aggregate_commitments(
            &commitments,
            Some(rands.to_vec()),
            &PIOPforDLComparison::s_minus_one_linear_combination(),
        )
        .unwrap();

        let nzk_s_minus_one_proof = NonZeroOverK::<F, PC, FS>::prove(
            ck,
            domain_k,
            &s_minus_one,
            &s_minus_one_commitment,
            &s_minus_one_rand,
            rng,
        )?;

        let proof = Proof {
            // Commitments
            s_commit: commitments[0].commitment().clone(),
            f_prime_commit: commitments[1].commitment().clone(),
            g_prime_commit: commitments[2].commitment().clone(),
            s_prime_commit: commitments[3].commitment().clone(),
            h_commit: commitments[4].commitment().clone(),

            // Proofs
            f_prime_square_proof,
            g_prime_square_proof,
            s_prime_square_proof,
            f_prime_product_proof,
            h_proof,
            f_prime_subset_proof,
            g_prime_subset_proof,
            s_prime_subset_proof,
            nzk_f_prime_proof,
            nzk_g_prime_proof,
            nzk_s_prime_proof,
            nzk_s_minus_one_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        ck: &PC::CommitterKey,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        f_commit: &LabeledCommitment<PC::Commitment>,
        g_commit: &LabeledCommitment<PC::Commitment>,
        enforced_degree_bound: Option<usize>,
        proof: Proof<F, PC>,
        fs_rng: &mut FS,
    ) -> Result<(), Error> {
        // re-label f and g with the enforced degree bound
        let f_commit = LabeledCommitment::new(
            f_commit.label().clone(),
            f_commit.commitment().clone(),
            enforced_degree_bound,
        );
        let g_commit = LabeledCommitment::new(
            g_commit.label().clone(),
            g_commit.commitment().clone(),
            enforced_degree_bound,
        );

        let commitments = vec![
            LabeledCommitment::new(String::from("s"), proof.s_commit, enforced_degree_bound),
            LabeledCommitment::new(
                String::from("f_prime"),
                proof.f_prime_commit,
                enforced_degree_bound,
            ),
            LabeledCommitment::new(
                String::from("g_prime"),
                proof.g_prime_commit,
                enforced_degree_bound,
            ),
            LabeledCommitment::new(
                String::from("s_prime"),
                proof.s_prime_commit,
                enforced_degree_bound,
            ),
            LabeledCommitment::new(String::from("h"), proof.h_commit, enforced_degree_bound),
        ];

        let fs_bytes =
            &to_bytes![Self::PROTOCOL_NAME, commitments].map_err(|_| Error::ToBytesError)?;
        fs_rng.absorb(fs_bytes);

        let alphas = [F::one(), F::one()];
        let square_check_vo = GenericShiftingVO::new(&[0, 1], &alphas, square_check)?;

        // Zero over K for f_prime
        ZeroOverK::<F, PC, FS>::verify(
            proof.f_prime_square_proof,
            &[f_commit.clone(), commitments[1].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            vk,
        )?;

        // Zero over K for g_prime
        ZeroOverK::<F, PC, FS>::verify(
            proof.g_prime_square_proof,
            &[g_commit.clone(), commitments[2].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            vk,
        )?;

        // Zero over K for s_prime
        ZeroOverK::<F, PC, FS>::verify(
            proof.s_prime_square_proof,
            &[commitments[0].clone(), commitments[3].clone()],
            enforced_degree_bound,
            &square_check_vo,
            &domain_k,
            vk,
        )?;

        let product_check_vo =
            GenericShiftingVO::new(&[0, 1, 2], &vec![F::one(); 3], presets::abc_product_check)?;

        // Zero over K for f' = (s')*(g')
        ZeroOverK::<F, PC, FS>::verify(
            proof.f_prime_product_proof,
            &[
                commitments[1].clone(),
                commitments[3].clone(),
                commitments[2].clone(),
            ],
            enforced_degree_bound,
            &product_check_vo,
            &domain_k,
            vk,
        )?;

        // Geometric Sequence Test for h
        let omega: F = domain_h.element(1);
        let delta_r = omega.sqrt();
        if delta_r.is_none() {
            return Err(Error::OmegaSqrtError);
        }

        let delta = delta_r.unwrap();

        let mut a_s = vec![F::one()];
        let mut c_s = vec![domain_h.size()];

        let to_pad = domain_k.size() - domain_h.size();
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }

        GeoSeqTest::<F, PC, FS>::verify(
            delta,
            &mut a_s,
            &mut c_s,
            &domain_k,
            &commitments[4],
            enforced_degree_bound,
            proof.h_proof,
            &vk,
        )?;

        // Subset over K for f'
        SubsetOverK::<F, PC, FS>::verify(proof.f_prime_subset_proof)?;

        // Subset over K for g'
        SubsetOverK::<F, PC, FS>::verify(proof.g_prime_subset_proof)?;

        // Subset over K for s'
        SubsetOverK::<F, PC, FS>::verify(proof.s_prime_subset_proof)?;

        // Non-zero over K for f′
        NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            commitments[1].commitment().clone(),
            enforced_degree_bound,
            proof.nzk_f_prime_proof,
        )?;

        // Non-zero over K for g′
        NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            commitments[2].commitment().clone(),
            enforced_degree_bound,
            proof.nzk_g_prime_proof,
        )?;

        // Non-zero over K for s′
        NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            commitments[3].commitment().clone(),
            enforced_degree_bound,
            proof.nzk_s_prime_proof,
        )?;

        // Non-zero over K for s(X) − 1
        let one_poly = LabeledPolynomial::new(
            String::from("one"),
            DensePolynomial::from_coefficients_slice(&[F::one()]),
            enforced_degree_bound,
            None,
        );
        let (commit_to_one, _) = PC::commit(ck, &[one_poly], None).map_err(to_pc_error::<F, PC>)?;

        let (s_minus_one_commitment, _) = PC::aggregate_commitments(
            &[commitments[0].clone(), commit_to_one[0].clone()],
            None,
            &PIOPforDLComparison::s_minus_one_linear_combination(),
        )?;

        NonZeroOverK::<F, PC, FS>::verify(
            &vk,
            &domain_k,
            s_minus_one_commitment.commitment().clone(),
            enforced_degree_bound,
            proof.nzk_s_minus_one_proof,
        )?;

        Ok(())
    }
}
