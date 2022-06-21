use crate::{
    commitment::HomomorphicPolynomialCommitment,
};
use ark_ff::{PrimeField, SquareRootField};
use std::marker::PhantomData;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
//pub struct Proof {
pub struct Proof<F: PrimeField + SquareRootField, PC: HomomorphicPolynomialCommitment<F>> {
    pub a_slt_proof: Vec<u8>,
    pub b_slt_proof: Vec<u8>,
    pub c_diag_proof: Vec<u8>,
    pub blah: PhantomData<F>,
    pub blah2: PhantomData<PC>,
}
