use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use crate::util::{commit_polynomial, powers_of};
use crate::virtual_oracle::VirtualOracle;
use crate::zero_over_k::piop::PIOPforZeroOverK;
use crate::zero_over_k::proof::Proof;
use crate::{label_commitment, label_polynomial};
use ark_ff::to_bytes;
use ark_ff::PrimeField;
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, PCRandomness};
use ark_std::marker::PhantomData;
use digest::Digest;
use rand::Rng;

mod piop;
mod proof;
mod tests;

pub struct ZeroOverK<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _field: PhantomData<F>,
    _polynomial_commitment_scheme: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> ZeroOverK<F, PC, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Zero Over K";

    pub fn prove<VO: VirtualOracle<F>, R: Rng>(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        concrete_oracle_commitments: &[PC::Commitment],
        concrete_oracle_commit_rands: &[PC::Randomness],
        virtual_oracle: &VO,
        domain: GeneralEvaluationDomain<F>,
        ck: &PC::CommitterKey,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        let prover_initial_state =
            PIOPforZeroOverK::prover_init(domain, concrete_oracles, virtual_oracle)?;
        let verifier_initial_state = PIOPforZeroOverK::verifier_init(domain, virtual_oracle)?;

        // TODO: FiatShamir more stuff!

        // let mut fs_rng = FiatShamirRng::from_seed(
        //     &to_bytes![&Self::PROTOCOL_NAME, concrete_oracles, virtual_oracle.alphas()].unwrap(),
        // );
        let mut fs_rng = FiatShamirRng::<D>::from_seed(
            &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments].unwrap(),
        );

        //------------------------------------------------------------------
        // First Round
        let (_, prover_first_oracles, prover_state) =
            PIOPforZeroOverK::prover_first_round(prover_initial_state, rng)?;

        let random_polynomials = prover_first_oracles.random_polynomials;
        let masking_polynomials = prover_first_oracles.masking_polynomials;
        let q_1 = prover_first_oracles.q_1;

        // commit to the random polynomials
        let (r_commitments, _) =
            PC::commit(ck, random_polynomials.iter(), None).map_err(to_pc_error::<F, PC>)?;

        // commit to the masking polynomials
        let (m_commitments, m_rands) =
            PC::commit(ck, masking_polynomials.iter(), Some(rng)).map_err(to_pc_error::<F, PC>)?;

        // commit to q_1
        let (q1_commit, q1_rand) =
            PC::commit(ck, &[q_1.clone()], None).map_err(to_pc_error::<F, PC>)?;

        fs_rng.absorb(&to_bytes![r_commitments, m_commitments, q1_commit].unwrap());

        let (verifier_first_msg, verifier_state) =
            PIOPforZeroOverK::verifier_first_round(verifier_initial_state, &mut fs_rng)?;
        //------------------------------------------------------------------

        //------------------------------------------------------------------
        // Second Round

        let (prover_second_msg, prover_second_oracles, prover_state) =
            PIOPforZeroOverK::prover_second_round(&verifier_first_msg, prover_state, rng);

        let points_for_queries = PIOPforZeroOverK::verifier_query_set(verifier_state)?;
        //------------------------------------------------------------------

        let q_2 = prover_second_oracles.q_2;
        let beta_1 = points_for_queries.beta_1;
        let beta_2 = points_for_queries.beta_2;

        // commit to q_1
        let (q1_commit, q1_rand) = commit_polynomial::<F, PC>(ck, &q_1, Some(rng))?;

        // open q1 at beta_1
        let q1_eval = q_1.evaluate(&beta_1);
        let q1_opening = PC::open(
            ck,
            &[q_1],
            &[label_commitment!(q1_commit)],
            &beta_1,
            F::one(),
            &[q1_rand],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let r_commitments = r_commitments
            .iter()
            .map(|r_c| r_c.commitment().clone())
            .collect::<Vec<_>>();

        let q2_eval = q_2.evaluate(&beta_2);

        // let (q2_commit, q2_rand) = commit_polynomial::<F, PC>(ck, &q2)?;
        let q2_commitment = PC::multi_scalar_mul(
            &r_commitments,
            powers_of(verifier_first_msg.c)
                .take(random_polynomials.len())
                .collect::<Vec<_>>()
                .as_slice(),
        );

        //since we commit with randomnes None, linear combination of all randomness will be just empty
        let homomorphic_randomness = PC::Randomness::empty();

        let q2_opening = PC::open(
            ck,
            &[q_2],
            &[label_commitment!(q2_commitment)],
            &beta_2,
            F::one(),
            &[homomorphic_randomness.clone()],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        // compute the masked oracles
        let h_primes = concrete_oracles
            .iter()
            .zip(masking_polynomials.iter())
            .map(|(oracle, masking_poly)| oracle.polynomial() + masking_poly.polynomial())
            .collect::<Vec<_>>();

        let alphas = virtual_oracle.alphas();

        let h_prime_evals = h_primes
            .iter()
            .zip(alphas.iter())
            .map(|(h_prime, &alpha)| h_prime.evaluate(&(alpha * beta_1)))
            .collect::<Vec<_>>();

        let m_evals = masking_polynomials
            .iter()
            .zip(alphas.iter())
            .map(|(m_poly, &alpha)| m_poly.evaluate(&(alpha * beta_2)))
            .collect::<Vec<_>>();

        // commit to the concrete oracle
        let (concrete_oracles_commitments, _) =
            PC::commit(ck, concrete_oracles.iter(), None).map_err(to_pc_error::<F, PC>)?;

        let h_prime_openings = h_primes
            .iter()
            .zip(alphas.iter())
            .zip(concrete_oracles_commitments.iter())
            .zip(m_commitments.iter())
            .zip(concrete_oracle_commit_rands.iter())
            .zip(m_rands.iter())
            .map(
                |(
                    (
                        (((h_prime, &alpha_i), oracle_commitment), m_commitment),
                        concrete_oracle_rand,
                    ),
                    m_rand,
                )|
                 -> Result<PC::Proof, Error> {
                    let h_prime_commitment = PC::multi_scalar_mul(
                        &[
                            oracle_commitment.commitment().clone(),
                            m_commitment.commitment().clone(),
                        ],
                        &[F::one(), F::one()],
                    );
                    let homomorphic_randomness =
                        PC::aggregate_randomness(&[concrete_oracle_rand.clone(), m_rand.clone()]);

                    PC::open(
                        ck,
                        &[label_polynomial!(h_prime)],
                        &[label_commitment!(h_prime_commitment)],
                        &(alpha_i * beta_1),
                        F::one(),
                        &[homomorphic_randomness],
                        None,
                    )
                    .map_err(to_pc_error::<F, PC>)
                },
            )
            .collect::<Result<Vec<_>, Error>>()?;

        let m_openings = masking_polynomials
            .iter()
            .zip(m_commitments.iter())
            .zip(alphas.iter())
            .map(
                |((m_poly, m_commit), &alpha_i)| -> Result<PC::Proof, Error> {
                    PC::open(
                        ck,
                        [m_poly],
                        &[m_commit.clone()],
                        &(alpha_i * beta_2),
                        F::one(),
                        &[homomorphic_randomness.clone()],
                        None,
                    )
                    .map_err(to_pc_error::<F, PC>)
                },
            )
            .collect::<Result<Vec<_>, Error>>()?;

        let m_commitments = m_commitments
            .iter()
            .map(|m_c| m_c.commitment().clone())
            .collect::<Vec<_>>();

        let proof = Proof {
            // commitments
            m_commitments,
            r_commitments,
            q1_commit,

            // evaluations
            q1_eval,
            q2_eval,
            h_prime_evals,
            m_evals,

            // opening
            q1_opening,
            q2_opening,
            m_openings,
            h_prime_openings,
        };

        Ok(proof)
    }

    pub fn verify<VO: VirtualOracle<F>>(
        proof: Proof<F, PC>,
        concrete_oracle_commitments: &[PC::Commitment],
        virtual_oracle: &VO,
        domain: GeneralEvaluationDomain<F>,
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        let verifier_initial_state = PIOPforZeroOverK::verifier_init(domain, virtual_oracle)?;

        // let mut fs_rng = FiatShamirRng::from_seed(
        //     &to_bytes![&Self::PROTOCOL_NAME, concrete_oracles, virtual_oracle.alphas()].unwrap(),
        // );
        let mut fs_rng = FiatShamirRng::<D>::from_seed(
            &to_bytes![&Self::PROTOCOL_NAME, concrete_oracle_commitments].unwrap(),
        );

        fs_rng
            .absorb(&to_bytes![proof.r_commitments, proof.m_commitments, proof.q1_commit].unwrap());

        let (verifier_first_msg, verifier_state) =
            PIOPforZeroOverK::verifier_first_round(verifier_initial_state, &mut fs_rng)?;

        let query_set = PIOPforZeroOverK::verifier_query_set_new(&verifier_state)?;

        let points_for_queries = PIOPforZeroOverK::verifier_query_set(verifier_state)?;

        // derive commitment to h_prime through additive homomorphism
        let h_prime_commitments = concrete_oracle_commitments
            .iter()
            .zip(proof.m_commitments.iter())
            .map(|(oracle_commitment, m_commitment)| {
                PC::multi_scalar_mul(
                    &[oracle_commitment.clone(), m_commitment.clone()],
                    &[F::one(), F::one()],
                )
            })
            .collect::<Vec<_>>();

        let check = PC::check(
            vk,
            &[label_commitment!(h_prime_commitments[0])],
            &(virtual_oracle.alphas()[0] * points_for_queries.beta_1),
            [proof.h_prime_evals[0]],
            &proof.h_prime_openings[0],
            F::one(),
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        assert!(check);

        // derive commitment to q2 through additive homomorphism
        let q2_commitment = PC::multi_scalar_mul(
            &proof.r_commitments,
            powers_of(verifier_first_msg.c)
                .take(concrete_oracle_commitments.len())
                .collect::<Vec<_>>()
                .as_slice(),
        );

        // compute M(beta_2)
        let big_m_at_beta_2 = &proof
            .m_evals
            .iter()
            .zip(powers_of(verifier_first_msg.c))
            .fold(F::zero(), |acc, (&m_eval, c_power)| {
                acc + (m_eval * c_power)
            });

        // evaluate z_k(beta_1), z_k(beta_2)
        let z_k_at_beta_1 = domain.evaluate_vanishing_polynomial(points_for_queries.beta_1);
        let z_k_at_beta_2 = domain.evaluate_vanishing_polynomial(points_for_queries.beta_2);

        // compute F_prime(beta_1)
        let f_prime_eval = VO::query(&proof.h_prime_evals);

        // check batched commitment

        let mut evaluations = ark_poly_commit::Evaluations::new();
        for (i, ((&alpha, &h_i_prime_eval), &m_eval)) in virtual_oracle
            .alphas()
            .iter()
            .zip(proof.h_prime_evals.iter())
            .zip(proof.m_evals.iter())
            .enumerate()
        {
            evaluations.insert(
                (format!("h_prime_{}", i), alpha * points_for_queries.beta_1),
                h_i_prime_eval,
            );
            evaluations.insert(
                (format!("m_{}", i), alpha * points_for_queries.beta_2),
                m_eval,
            );
        }

        //let's say that we have another hash map: poly_label => poly, rand

        // evaluate poly: h_prime_0, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
        // evaluate poly: h_prime_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
        // evaluate poly: m_0, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
        // evaluate poly: m_1, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"
        // evaluate poly: q_1, in point: beta_1 with value Fp256 "(293D6AB4C074E7502C1DBD28D1151523DB175383DDBC82E5392B6EA5B70C6A7A)"
        // evaluate poly: q_2, in point: beta_2 with value Fp256 "(1A6BFDB9903D732617070351C48D21919D5715900E9204C2E3A0970EC36CE81D)"

        /// point_label = [poly_labels] for each poly label take poly from poly_label map

        /// iterate query_set
        /// for each point_label, take all poly_labels where beta_1 => (h_prime_0, h_prime_1, q_1)
        /// for each point_label, take all poly_labels where beta_2 => (m_0, m_1, q_2)
        // match PC::batch_check(
        //     &vk,
        //     &[w_commit_labeled, w_commit_shifted_labeled],
        //     &query_set,
        //     &evaluations,
        //     &self.batch_opening,
        //     u,
        //     &mut OsRng,
        // ) {
        //     Ok(true) => Ok(()),
        //     Ok(false) => Err(Error::ProofVerificationError),
        //     Err(e) => panic!("{:?}", e),
        // };

        // check that M(beta_2) - q2(beta_2)*zK(beta_2) = 0
        let check_1 = *big_m_at_beta_2 - proof.q2_eval * z_k_at_beta_2;
        if check_1 != F::zero() {
            return Err(Error::Check1Failed);
        }

        //check that F_prime(beta_1) - q1(beta_1)*zK(beta_1) = 0
        let check_2 = f_prime_eval - proof.q1_eval * z_k_at_beta_1;
        if check_2 != F::zero() {
            return Err(Error::Check2Failed);
        }

        //TODO
        // we skip checking openings until we batch them into one commitment

        Ok(())
    }
}
