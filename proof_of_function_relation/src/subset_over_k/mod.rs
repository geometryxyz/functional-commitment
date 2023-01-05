use crate::{error::Error, subset_over_k::proof::Proof};
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;
use fiat_shamir_rng::FiatShamirRng;
use homomorphic_poly_commit::AdditivelyHomomorphicPCS;

pub mod proof;

// TODO: implement SubsetOverK, currently just a placeholder
pub struct SubsetOverK<F: PrimeField, PC: AdditivelyHomomorphicPCS<F>, FS: FiatShamirRng> {
    _field: PhantomData<F>,
    _pc: PhantomData<PC>,
    _fs: PhantomData<FS>,
}

impl<F, PC, FS> SubsetOverK<F, PC, FS>
where
    F: PrimeField,
    PC: AdditivelyHomomorphicPCS<F>,
    FS: FiatShamirRng,
{
    pub fn prove() -> Proof {
        Proof {}
    }

    pub fn verify(_proof: Proof) -> Result<(), Error> {
        Ok(())
    }
}
