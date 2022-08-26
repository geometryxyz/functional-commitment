#[derive(PartialEq, Debug)]
pub enum Error {
    VarAlreadyExists(String),
    VarMissing(String),
}
