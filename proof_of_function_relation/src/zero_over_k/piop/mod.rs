use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;

mod prover;
mod verifier;

/// A labeled DensePolynomial with coefficients over `F`
pub type LabeledPolynomial<F> = ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>;

#[allow(dead_code)]
pub struct PIOPforZeroOverK<F: PrimeField, VO: VirtualOracle<F>> {
    _field: PhantomData<F>,
    _oracle: PhantomData<VO>,
}
