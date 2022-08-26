use ark_ff::PrimeField;

#[derive(Clone)]
pub enum VariableType {
    INPUT,
    WITNESS,
    OUTPUT,
}
#[derive(Clone)]
pub struct Variable<F: PrimeField> {
    pub label: String,
    pub variable_type: VariableType,
    pub value: F,
}
