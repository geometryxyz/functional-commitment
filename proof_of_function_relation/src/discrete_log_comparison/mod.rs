use crate::commitment::HomomorphicPolynomialCommitment;
use ark_ff::{PrimeField, SquareRootField};
use ark_std::marker::PhantomData;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`

pub mod piop;
mod tests;

pub struct DLComparison<
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
> {
    _field: PhantomData<F>,
    _polynomial_commitment_scheme: PhantomData<PC>,
    _digest: PhantomData<D>,
}

impl<F, PC, D> DLComparison<F, PC, D>
where
    F: PrimeField + SquareRootField,
    PC: HomomorphicPolynomialCommitment<F>,
    D: Digest,
{
    pub const PROTOCOL_NAME: &'static [u8] = b"Discrete-log Comparison";
    
    pub fn prove() {

    }

    pub fn verify() {
        
    }
}
            