use ark_ff::Field;

#[derive(Clone)]
pub enum VariableType {
    Input,
    Witness,
}
#[derive(Clone)]
pub struct Variable<F: Field> {
    pub label: String,
    pub variable_type: VariableType,
    pub value: F,
}
