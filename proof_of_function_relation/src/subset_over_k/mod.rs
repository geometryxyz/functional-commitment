use crate::{error::Error, subset_over_k::proof::Proof};
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;
use digest::Digest;
use zero_over_k::commitment::AdditivelyHomomorphicPCS;

pub mod proof;

// TODO: implement SubsetOverK, currently just a placeholder
pub struct SubsetOverK<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, D: Digest> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> SubsetOverK<F, PC, D>
where
    F: PrimeField,
    PC: AdditivelyHomomorphicPCS<F>,
    D: Digest,
{
    pub fn prove() -> Proof {
        Proof {}
    }

    pub fn verify(_proof: Proof) -> Result<(), Error> {
        Ok(())
    }
}
