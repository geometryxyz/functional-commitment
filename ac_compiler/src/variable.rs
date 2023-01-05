use ark_ff::Field;

#[derive(Clone, Copy)]
pub enum VariableType {
    Input,
    Witness,
    Output,
}
#[derive(Clone)]
pub struct Variable<F: Field> {
    pub label: String,
    pub variable_type: VariableType,
    pub value: F,
}
