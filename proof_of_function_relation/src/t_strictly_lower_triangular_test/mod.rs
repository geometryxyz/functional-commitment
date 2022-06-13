use crate::{
    commitment::HomomorphicPolynomialCommitment,
    discrete_log_comparison::DLComparison,
    error::{to_pc_error, Error},
    t_strictly_lower_triangular_test::proof::Proof,
    virtual_oracle::{inverse_check_oracle::InverseCheckOracle, VirtualOracle},
    geo_seq::{GeoSeqTest},
    subset_over_k::SubsetOverK,
};
use ark_ff::{to_bytes, PrimeField, SquareRootField};
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial};
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use std::marker::PhantomData;

pub mod proof;
mod tests;

struct TStrictlyLowerTriangular<
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> TStrictlyLowerTriangular<F, PC, D>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"t-Strictly Lower Triangular Test";

    pub fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        row_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        col_poly: &LabeledPolynomial<F, DensePolynomial<F>>,
        row_commit: &LabeledCommitment<PC::Commitment>,
        col_commit: &LabeledCommitment<PC::Commitment>,
        fs_rng: &mut FiatShamirRng<D>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        fs_rng.absorb(&to_bytes![Self::PROTOCOL_NAME].unwrap());

        // Step 2: Geometric sequence test on h
        let r = domain_h.element(1);
        
        if t > domain_h.size() {
            return Err(Error::T2Large)
        }
        let mut a_s = vec![domain_h.element(t)];
        let mut c_s = vec![domain_h.size() - t];
        
        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }
        let geo_seq_proof = GeoSeqTest::<F, PC, D>::prove(
            ck,
            r,
            &a_s,
            &c_s,
            domain_k,
            rng,
        )?;

        // Step 3: Subset over K between row_M and h
        let subset_proof = SubsetOverK::<F, PC, D>::prove();

        // Step 4: Discrete Log Comparison between row_M and col_M
        let dl_proof = DLComparison::<F, PC, D>::prove(
            ck, domain_k, domain_h, row_poly, col_poly, row_commit, col_commit, fs_rng, rng,
        )?;

        let proof = Proof { dl_proof, geo_seq_proof, subset_proof };

        Ok(proof)
    }

    pub fn verify(
        vk: &PC::VerifierKey,
        ck: &PC::CommitterKey,
        t: usize,
        domain_k: &GeneralEvaluationDomain<F>,
        domain_h: &GeneralEvaluationDomain<F>,
        row_commit: &LabeledCommitment<PC::Commitment>,
        col_commit: &LabeledCommitment<PC::Commitment>,
        proof: Proof<F, PC>,
        fs_rng: &mut FiatShamirRng<D>,
    ) -> Result<(), Error> {

        // Step 2: Geometric sequence test on h
        let mut a_s = vec![domain_h.element(t)];
        let mut c_s = vec![domain_h.size() - t];
        
        let to_pad = domain_k.size() - (domain_h.size() - t);
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }
        GeoSeqTest::<F, PC, D>::verify(
            domain_h.element(1), 
            &a_s, 
            &c_s, 
            domain_k, 
            proof.geo_seq_proof,
            vk
        )?;

        // Step 3: Subset over K between row_M and h
        SubsetOverK::<F, PC, D>::verify(proof.subset_proof)?;

        // Step 4: Discrete Log Comparison between row_M and col_M
        DLComparison::<F, PC, D>::verify(
            vk,
            ck,
            domain_k,
            domain_h,
            row_commit,
            col_commit,
            proof.dl_proof,
            fs_rng,
        )?;

        Ok(())
    }
}
