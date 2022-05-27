use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;

mod prover;
mod verifier;

/// A labeled DensePolynomial with coefficients over `F`
pub type LabeledPolynomial<F> = ark_poly_commit::LabeledPolynomial<F, DensePolynomial<F>>;

pub struct PIOPforZeroOverK<F: Field> {
    _field: PhantomData<F>,
}
