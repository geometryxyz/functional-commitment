use ark_ff::PrimeField;
use ark_poly_commit::{LinearCombination, PolynomialLabel};
use ark_std::marker::PhantomData;

pub struct PIOPforTDiagTest<F: PrimeField> {
    _field: PhantomData<F>,
}

impl<F: PrimeField> PIOPforTDiagTest<F> {
    pub fn generate_h_linear_combination() -> LinearCombination<F> {
        LinearCombination::new("h", vec![(F::one(), "h1"), (F::one(), "h2")])
    }

    #[allow(non_snake_case)]
    pub fn generate_valM_plus_h2_linear_combination(
        val_label: &PolynomialLabel,
    ) -> LinearCombination<F> {
        LinearCombination::new(
            "val_plus_h2",
            vec![(F::one(), val_label.clone()), (F::one(), "h2".to_string())],
        )
    }
}
