use ark_ff::{FftField, Field};

use crate::vo_constant;

use super::vo_term::VOTerm;

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[1] and terms[2] are equal for all points in the K
pub fn equality_check<F: Field>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[1].clone() - terms[2].clone()
}

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[1] is the inverse of terms[2] for all points in K
pub fn inverse_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[1].clone() * terms[2].clone() - vo_constant!(F::one())
}

/// A function to be used in a virtual oracle. Should the VO evaluate to 0 for all points in a domain K,
/// we can conclude that terms[1] is the square-root of terms[2] for all points in K
pub fn square_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[1].clone() - terms[2].clone() * terms[2].clone()
}

/// A function to be used in a virtual oracle. This function can be used to check that the product of
/// terms[1] and terms[2] is 0 over a domain K
pub fn zero_product_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[1].clone() * terms[2].clone()
}

/// A function to be used in a virtual oracle. This function can be used to check that terms[1] agrees with
/// the product of terms[2] and terms[3] over a domain K
pub fn abc_product_check<F: FftField>(terms: &[VOTerm<F>]) -> VOTerm<F> {
    terms[1].clone() - terms[2].clone() * terms[3].clone()
}

/// A closure to be used in a geometric sequence test virtual oracle.
/// This specific construction expects terms[0] to be X (as enforced by default), terms[1] to be f(X) and
/// terms[2] to be f(gamma*X); where X is an indeterminate variable, f is a polynomial for which we make
/// a geometric sequence claim and gamma is the generator of the subgroup of interest.
#[macro_export]
macro_rules! geometric_seq_check {
    ($common_ratio:expr, $sequence_lengths:expr, $domain:expr) => {
        |terms: &[VOTerm<F>]| {
            // construct (f(gamma * X) - r * f(X))
            let check_next_term = terms[2].clone() - vo_constant!($common_ratio) * terms[1].clone();

            let mut evaluation_function = check_next_term;
            let mut current_sequence_starting_index = 0;
            // construct each stitch and multiply to the final polynomial
            $sequence_lengths
                .iter()
                .for_each(|current_sequence_length| {
                    let stitch = terms[0].clone()
                        - vo_constant!($domain.element(1).pow([(current_sequence_starting_index
                            + current_sequence_length
                            - 1)
                            as u64]));
                    current_sequence_starting_index += current_sequence_length;
                    evaluation_function = evaluation_function.clone() * stitch;
                });

            evaluation_function
        }
    };
}
