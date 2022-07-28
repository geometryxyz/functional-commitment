use ark_ff::{FftField, Field};

use crate::vo_constant;

use super::vo_term::VOTerm;

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[0] and terms[1] are equal for all points in the K
pub fn equality_check<F: Field>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[0].clone() - terms[1].clone()
}

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[0] is the inverse of terms[1] for all points in K
pub fn inverse_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[0].clone() * terms[1].clone() - vo_constant!(F::one())
}

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[0] is the square-root of terms[1] for all points in K
pub fn square_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[0].clone() - terms[1].clone() * terms[1].clone()
}
