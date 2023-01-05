use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

pub mod prover;

#[allow(dead_code)]
pub struct PIOPforNonZeroOverK<F: PrimeField + FftField> {
    _field: PhantomData<F>,
}
