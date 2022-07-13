use ark_ff::{PrimeField, SquareRootField};
use ark_poly_commit::LinearCombination;
use ark_std::marker::PhantomData;

pub mod prover;

pub struct PIOPforDLComparison<F: PrimeField + SquareRootField> {
    _field: PhantomData<F>,
}

impl<F: PrimeField + SquareRootField> PIOPforDLComparison<F> {
    pub fn s_minus_one_linear_combination() -> LinearCombination<F> {
        LinearCombination::new("s_minus_one", vec![(F::one(), "s"), (-F::one(), "one")])
    }
}
