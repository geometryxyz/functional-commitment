use crate::{
    commitment::AdditivelyHomomorphicPCS,
    error::{to_pc_error, Error},
    geo_seq::GeoSeqTest,
    label_polynomial,
    non_zero_over_k::NonZeroOverK,
    t_diag::proof::Proof,
    util::generate_sequence,
    virtual_oracle::{eq_vo::EqVO, prod_vo::ProdVO},
    zero_over_k::ZeroOverK,
};
use ark_ff::{PrimeField, SquareRootField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, PCRandomness};
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use std::marker::PhantomData;

use self::piop::PIOPforTDiagTest;

pub mod piop;
pub mod proof;
mod tests;

pub struct TDiag<F: PrimeField + SquareRootField, PC: AdditivelyHomomorphicPCS<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> TDiag<F, PC, D>
where
    F: PrimeField + SquareRootField,
    PC: AdditivelyHomomorphicPCS<F>,
    D: Digest,
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
        _row_m_random: &PC::Randomness,
        _col_m_random: &PC::Randomness,
        _val_m_random: &PC::Randomness,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        if t > domain_h.size() {
            return Err(Error::T2Large);
        }

        // Step 1a produce h1 = w^t, w^(t+1), ..., w^(n-1), 0, 0, ..., 0
        let r_h1 = domain_h.element(1);
        let mut a_s_h1 = vec![domain_h.element(t)];
        let mut c_s_h1 = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s_h1.push(F::zero());
            c_s_h1.push(to_pad);
        }

        let seq = generate_sequence::<F>(r_h1, &a_s_h1.as_slice(), &c_s_h1.as_slice());
        let h1 = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let h1 = label_polynomial!(h1);

        // Step 1b produce h2 = 0, 0, ..., 0, 1, 1, ..., 1
        let r_h2 = domain_h.element(0);
        let mut a_s_h2 = vec![F::zero()];
        let mut c_s_h2 = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s_h2.push(F::one());
            c_s_h2.push(to_pad);
        }

        let seq = generate_sequence::<F>(r_h2, &a_s_h2.as_slice(), &c_s_h2.as_slice());
        let h2 = DensePolynomial::<F>::from_coefficients_slice(&domain_k.ifft(&seq));
        let h2 = label_polynomial!(h2);

        let (h_commitments, h_rands) =
            PC::commit(ck, &[h1.clone(), h2.clone()], None).map_err(to_pc_error::<F, PC>)?;

        // Step 2: Geometric Sequence Test on h1
        let h1_seq_proof = GeoSeqTest::<F, PC, D>::prove(
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
        let h2_seq_proof = GeoSeqTest::<F, PC, D>::prove(
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
        let h = label_polynomial!(h);
        let alphas = vec![F::one(), F::one()];

        let h_commitment = PC::get_commitments_lc(
            &h_commitments,
            &PIOPforTDiagTest::generate_h_linear_combination(),
        )?;

        // Step 4b: Zero over K for h = rowM
        let eq_vo = EqVO::new();
        let h_eq_row_m = ZeroOverK::<F, PC, D>::prove(
            &[h.clone(), row_m.clone()],
            &[h_commitment.clone(), row_m_commitment.clone()],
            &[PC::Randomness::empty(), PC::Randomness::empty()],
            &eq_vo,
            &alphas,
            domain_k,
            ck,
            rng,
        )?;

        // Step 4c: Zero over K for rowM = colM
        let row_m_eq_col_m = ZeroOverK::<F, PC, D>::prove(
            &[row_m.clone(), col_m.clone()],
            &[row_m_commitment.clone(), col_m_commitment.clone()],
            &[PC::Randomness::empty(), PC::Randomness::empty()],
            &eq_vo,
            &alphas,
            domain_k,
            ck,
            rng,
        )?;

        // Step 5: Zero over K for valM * h2 = 0
        let prod_vo = ProdVO::new();
        let val_m_times_h2_proof = ZeroOverK::<F, PC, D>::prove(
            &[val_m.clone(), h2.clone()],
            &[val_m_commitment.clone(), h_commitments[1].clone()],
            &[PC::Randomness::empty(), PC::Randomness::empty()],
            &prod_vo,
            &alphas,
            domain_k,
            ck,
            rng,
        )?;

        // Step 6: Non-zero over K for valM + h2 != 0
        let val_plus_h2 = val_m.polynomial() + h2.polynomial();
        let val_plus_h2 = label_polynomial!(val_plus_h2);

        let val_plus_h2_proof = NonZeroOverK::<F, PC, D>::prove(ck, domain_k, &val_plus_h2, rng)?;

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
        domain_h: &GeneralEvaluationDomain<F>,
        domain_k: &GeneralEvaluationDomain<F>,
        proof: Proof<F, PC>,
    ) -> Result<(), Error> {
        // Step 2: Geometric Sequence Test on h1
        let r_h1 = domain_h.element(1);
        let mut a_s_h1 = vec![domain_h.element(t)];
        let mut c_s_h1 = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s_h1.push(F::zero());
            c_s_h1.push(to_pad);
        }

        let h_commitments = vec![
            LabeledCommitment::new(String::from("h1"), proof.h1_commit.clone(), None),
            LabeledCommitment::new(String::from("h2"), proof.h2_commit.clone(), None),
        ];

        GeoSeqTest::<F, PC, D>::verify(
            r_h1,
            &a_s_h1,
            &c_s_h1,
            domain_k,
            &h_commitments[0],
            proof.h1_seq_proof,
            vk,
        )?;

        // h2 params
        let r_h2 = F::one();
        let mut a_s_h2 = vec![F::zero()];
        let mut c_s_h2 = vec![domain_h.size() - t];

        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s_h2.push(F::one());
            c_s_h2.push(to_pad);
        }

        // Step 3: Geometric Sequence Test on h2
        GeoSeqTest::<F, PC, D>::verify(
            r_h2,
            &a_s_h2,
            &c_s_h2,
            domain_k,
            &h_commitments[1],
            proof.h2_seq_proof,
            vk,
        )?;

        // Step 4a: Verifier derives a commitment to h = h1 + h2
        let alphas = [F::one(), F::one()];
        let h_commit = PC::get_commitments_lc(
            &h_commitments,
            &PIOPforTDiagTest::generate_h_linear_combination(),
        )?;

        // Step 4b: Zero over K for h = rowM
        let eq_vo = EqVO::new();
        ZeroOverK::<F, PC, D>::verify(
            proof.h_eq_row_m,
            vec![h_commit.clone(), row_m_commitment.clone()].as_slice(),
            &eq_vo,
            domain_k,
            &alphas,
            vk,
        )?;

        // Step 4c: Zero over K for rowM = colM
        ZeroOverK::<F, PC, D>::verify(
            proof.row_m_eq_col_m,
            vec![row_m_commitment.clone(), col_m_commitment.clone()].as_slice(),
            &eq_vo,
            domain_k,
            &alphas,
            vk,
        )?;

        // Step 5: Zero over K for valM * h2 = 0
        let prod_vo = ProdVO::new();
        ZeroOverK::<F, PC, D>::verify(
            proof.val_m_times_h2_proof,
            vec![val_m_commitment.clone(), h_commitments[1].clone()].as_slice(),
            &prod_vo,
            domain_k,
            &alphas,
            vk,
        )?;

        // Step 6: Non-zero over K for valM + h2 != 0
        let val_plus_h2_commit = PC::get_commitments_lc(
            &[val_m_commitment.clone(), h_commitments[1].clone()],
            &PIOPforTDiagTest::generate_valM_plus_h2_linear_combination(val_m_commitment.label()),
        )?;

        NonZeroOverK::<F, PC, D>::verify(
            vk,
            domain_k,
            val_plus_h2_commit,
            proof.val_plus_h2_proof,
        )?;

        Ok(())
    }
}
