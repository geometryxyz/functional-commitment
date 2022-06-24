use crate::{
    commitment::HomomorphicPolynomialCommitment, error::Error, subset_over_k::proof::Proof,
};
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;
use digest::Digest;

pub mod proof;

// TODO: implement SubsetOverK, currently just a placeholder
pub struct SubsetOverK<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> SubsetOverK<F, PC, D>
where
    F: PrimeField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
{
    pub fn prove() -> Proof {
        Proof {}
    }

    pub fn verify(_proof: Proof) -> Result<(), Error> {
        Ok(())
    }
}
