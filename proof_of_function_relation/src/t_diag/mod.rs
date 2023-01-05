use crate::{
    error::{to_pc_error, Error},
    geo_seq::GeoSeqTest,
    non_zero_over_k::NonZeroOverK,
    t_diag::proof::Proof,
    util::generate_sequence,
};
use ark_ff::{PrimeField, SquareRootField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;
use rand::Rng;
use std::marker::PhantomData;
use zero_over_k::{
    virtual_oracle::generic_shifting_vo::{
        presets::{self, zero_product_check},
        GenericShiftingVO,
    },
    zero_over_k::ZeroOverK,
};

use self::piop::PIOPforTDiagTest;

pub mod piop;
pub mod proof;
mod tests;

pub struct TDiag<
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F, PC, FS> TDiag<F, PC, FS>
where
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
{
    #[allow(dead_code)]
    pub const PROTOCOL_NAME: &'static [u8] = b"t-Diagonal Test";

    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        t: usize,
        row_m: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_m: &LabeledPolynomial<F, DensePolynomial<F>>,
        val_m: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_m_commitment: &LabeledCommitment<PC::Commitment>,
        col_m_commitment: &LabeledCommitment<PC::Commitment>,
        val_m_commitment: &LabeledCommitment<PC::Commitment>,
        row_m_random: &PC::Randomness,
        col_m_random: &PC::Randomness,
        val_m_random: &PC::Randomness,
        enforced_degree_bound: Option<usize>,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        number_of_constraints: usize,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        if t > domain_h.size() {
            return Err(Error::T2Large);
        }

        // Step 1a produce h1 = w^t, w^(t+1), ..., w^(n-1), 0, 0, ..., 0
        let r_h1 = domain_h.element(1);
        let mut a_s_h1 = vec![domain_h.element(t)];
        let mut c_s_h1 = vec![number_of_constraints - t];

        let to_pad = domain_k.size() - (number_of_constraints - t);
        if to_pad > 0 {
            a_s_h1.push(F::zero());
            c_s_h1.push(to_pad);
        }

        let seq = generate_sequence::<F>(r_h1, &a_s_h1.as_slice(), &c_s_h1.as_slice());
        let h1 = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let h1 = LabeledPolynomial::new(String::from("h1"), h1, enforced_degree_bound, Some(1));

        // Step 1b produce h2 = 0, 0, ..., 0, 1, 1, ..., 1
        let r_h2 = domain_h.element(0);
        let mut a_s_h2 = vec![F::zero()];
        let mut c_s_h2 = vec![number_of_constraints - t];

        let to_pad = domain_k.size() - (number_of_constraints - t);
        if to_pad > 0 {
            a_s_h2.push(F::one());
            c_s_h2.push(to_pad);
        }

        let seq = generate_sequence::<F>(r_h2, &a_s_h2.as_slice(), &c_s_h2.as_slice());
        let h2 = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let h2 = LabeledPolynomial::new(String::from("h2"), h2, enforced_degree_bound, Some(1));

        let (h_commitments, h_rands) =
            PC::commit(ck, &[h1.clone(), h2.clone()], Some(rng)).map_err(to_pc_error::<F, PC>)?;

        // Step 2: Geometric Sequence Test on h1
        let h1_seq_proof = GeoSeqTest::<F, PC, FS>::prove(
            ck,
            r_h1,
            &h1,
            &h_commitments[0],
            &h_rands[0],
            &a_s_h1,
            &c_s_h1,
            domain_k,
            rng,
        )?;

        // Step 3: Geometric Sequence Test on h2
        let h2_seq_proof = GeoSeqTest::<F, PC, FS>::prove(
            ck,
            r_h2,
            &h2,
            &h_commitments[1],
            &h_rands[1],
            &a_s_h2,
            &c_s_h2,
            domain_k,
            rng,
        )?;

        // Step 4a: h = h1 + h2
        let h = h1.polynomial() + h2.polynomial();
        let h = LabeledPolynomial::new(String::from("h"), h, enforced_degree_bound, Some(1));
        let alphas = vec![F::one(), F::one()];

        let (h_commitment, h_rand) = PC::aggregate_commitments(
            &h_commitments,
            Some(h_rands.to_vec()),
            &PIOPforTDiagTest::generate_h_linear_combination(),
        )?;

        // Step 4b: Zero over K for h = rowM
        let eq_vo = GenericShiftingVO::new(&vec![0, 1], &alphas, presets::equality_check)?;
        let h_eq_row_m = ZeroOverK::<F, PC, FS>::prove(
            &[h.clone(), row_m.clone()],
            &[h_commitment.clone(), row_m_commitment.clone()],
            &[h_rand.clone(), row_m_random.clone()],
            enforced_degree_bound,
            &eq_vo,
            domain_k,
            ck,
            rng,
        )?;

        // Step 4c: Zero over K for rowM = colM
        let row_m_eq_col_m = ZeroOverK::<F, PC, FS>::prove(
            &[row_m.clone(), col_m.clone()],
            &[row_m_commitment.clone(), col_m_commitment.clone()],
            &[row_m_random.clone(), col_m_random.clone()],
            enforced_degree_bound,
            &eq_vo,
            domain_k,
            ck,
            rng,
        )?;

        // Step 5: Zero over K for valM * h2 = 0
        let prod_vo = GenericShiftingVO::new(&[0, 1], &[F::one(), F::one()], zero_product_check)?;
        let val_m_times_h2_proof = ZeroOverK::<F, PC, FS>::prove(
            &[val_m.clone(), h2.clone()],
            &[val_m_commitment.clone(), h_commitments[1].clone()],
            &[val_m_random.clone(), h_rands[1].clone()],
            enforced_degree_bound,
            &prod_vo,
            domain_k,
            ck,
            rng,
        )?;

        // Step 6: Non-zero over K for valM + h2 != 0
        let val_plus_h2 = val_m.polynomial() + h2.polynomial();
        let val_plus_h2 = LabeledPolynomial::new(
            String::from("val_plus_h2"),
            val_plus_h2,
            enforced_degree_bound,
            Some(1),
        );
        let (val_plus_h2_commit, val_plus_h2_rand) = PC::aggregate_commitments(
            &[h_commitments[1].clone(), val_m_commitment.clone()],
            Some(vec![h_rands[1].clone(), val_m_random.clone()]),
            &PIOPforTDiagTest::generate_valM_plus_h2_linear_combination(&String::from(
                val_m.label(),
            )),
        )?;

        let val_plus_h2_proof = NonZeroOverK::<F, PC, FS>::prove(
            ck,
            domain_k,
            &val_plus_h2,
            &val_plus_h2_commit,
            &val_plus_h2_rand,
            rng,
        )?;

        let proof = Proof {
            h1_commit: h_commitments[0].commitment().clone(),
            h2_commit: h_commitments[1].commitment().clone(),
            h1_seq_proof,
            h2_seq_proof,
            h_eq_row_m,
            row_m_eq_col_m,
            val_m_times_h2_proof,
            val_plus_h2_proof,
        };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        t: usize,
        row_m_commitment: &LabeledCommitment<PC::Commitment>,
        col_m_commitment: &LabeledCommitment<PC::Commitment>,
        val_m_commitment: &LabeledCommitment<PC::Commitment>,
        enforced_degree_bound: Option<usize>,
        domain_h: &GeneralEvaluationDomain<F>,
        domain_k: &GeneralEvaluationDomain<F>,
        number_of_constraints: usize,
        proof: Proof<F, PC>,
    ) -> Result<(), Error> {
        // re-label the oracle commitments with the enforced degree bound
        let row_m_commitment = LabeledCommitment::new(
            row_m_commitment.label().clone(),
            row_m_commitment.commitment().clone(),
            enforced_degree_bound,
        );
        let col_m_commitment = LabeledCommitment::new(
            col_m_commitment.label().clone(),
            col_m_commitment.commitment().clone(),
            enforced_degree_bound,
        );
        let val_m_commitment = LabeledCommitment::new(
            val_m_commitment.label().clone(),
            val_m_commitment.commitment().clone(),
            enforced_degree_bound,
        );

        // Step 2: Geometric Sequence Test on h1
        let r_h1 = domain_h.element(1);
        let mut a_s_h1 = vec![domain_h.element(t)];
        let mut c_s_h1 = vec![number_of_constraints - t];

        let to_pad = domain_k.size() - (number_of_constraints - t);
        if to_pad > 0 {
            a_s_h1.push(F::zero());
            c_s_h1.push(to_pad);
        }

        let h_commitments = vec![
            LabeledCommitment::new(
                String::from("h1"),
                proof.h1_commit.clone(),
                enforced_degree_bound,
            ),
            LabeledCommitment::new(
                String::from("h2"),
                proof.h2_commit.clone(),
                enforced_degree_bound,
            ),
        ];

        GeoSeqTest::<F, PC, FS>::verify(
            r_h1,
            &a_s_h1,
            &c_s_h1,
            domain_k,
            &h_commitments[0],
            enforced_degree_bound,
            proof.h1_seq_proof,
            vk,
        )?;

        // h2 params
        let r_h2 = F::one();
        let mut a_s_h2 = vec![F::zero()];
        let mut c_s_h2 = vec![number_of_constraints - t];

        let to_pad = domain_k.size() - (number_of_constraints - t);
        if to_pad > 0 {
            a_s_h2.push(F::one());
            c_s_h2.push(to_pad);
        }

        // Step 3: Geometric Sequence Test on h2
        GeoSeqTest::<F, PC, FS>::verify(
            r_h2,
            &a_s_h2,
            &c_s_h2,
            domain_k,
            &h_commitments[1],
            enforced_degree_bound,
            proof.h2_seq_proof,
            vk,
        )?;

        // Step 4a: Verifier derives a commitment to h = h1 + h2
        let alphas = [F::one(), F::one()];
        let (h_commit, _) = PC::aggregate_commitments(
            &h_commitments,
            None,
            &PIOPforTDiagTest::generate_h_linear_combination(),
        )?;

        // Step 4b: Zero over K for h = rowM
        let eq_vo = GenericShiftingVO::new(&vec![0, 1], &alphas, presets::equality_check)?;
        ZeroOverK::<F, PC, FS>::verify(
            proof.h_eq_row_m,
            vec![h_commit.clone(), row_m_commitment.clone()].as_slice(),
            enforced_degree_bound,
            &eq_vo,
            domain_k,
            vk,
        )?;

        // Step 4c: Zero over K for rowM = colM
        ZeroOverK::<F, PC, FS>::verify(
            proof.row_m_eq_col_m,
            vec![row_m_commitment.clone(), col_m_commitment.clone()].as_slice(),
            enforced_degree_bound,
            &eq_vo,
            domain_k,
            vk,
        )?;

        // Step 5: Zero over K for valM * h2 = 0
        let prod_vo = GenericShiftingVO::new(&[0, 1], &[F::one(), F::one()], zero_product_check)?;
        ZeroOverK::<F, PC, FS>::verify(
            proof.val_m_times_h2_proof,
            vec![val_m_commitment.clone(), h_commitments[1].clone()].as_slice(),
            enforced_degree_bound,
            &prod_vo,
            domain_k,
            vk,
        )?;

        // Step 6: Non-zero over K for valM + h2 != 0
        let (val_plus_h2_commit, _) = PC::aggregate_commitments(
            &[val_m_commitment.clone(), h_commitments[1].clone()],
            None,
            &PIOPforTDiagTest::generate_valM_plus_h2_linear_combination(val_m_commitment.label()),
        )?;

        NonZeroOverK::<F, PC, FS>::verify(
            vk,
            domain_k,
            val_plus_h2_commit.commitment().clone(),
            enforced_degree_bound,
            proof.val_plus_h2_proof,
        )?;

        Ok(())
    }
}
