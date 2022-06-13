use ark_ff::{PrimeField, SquareRootField};
use ark_std::marker::PhantomData;

pub mod prover;

#[allow(dead_code)]
pub struct PIOPforDLComparison<F: PrimeField + SquareRootField> {
    _field: PhantomData<F>,
}
