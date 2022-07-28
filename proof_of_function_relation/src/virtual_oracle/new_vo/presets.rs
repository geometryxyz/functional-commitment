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

#[macro_export]
macro_rules! geometric_seq_check {
    ($common_ratio:expr, $sequence_lengths:expr, $domain:expr) => {
        |terms: &[VOTerm<F>]| {
            // construct (f(gamma * x) - r * f(x))
            let check_next_term = terms[2].clone() - vo_constant!($common_ratio) * terms[1].clone();

            let mut evaluation_function = check_next_term;
            let mut starting_index = 0;
            // construct each stitch and multiply to the final polynomial
            $sequence_lengths.iter().for_each(|sequence_length| {
                let stitch = terms[0].clone()
                    - vo_constant!($domain
                        .element(1)
                        .pow([(starting_index + sequence_length - 1) as u64]));
                starting_index += sequence_length;
                evaluation_function = evaluation_function.clone() * stitch;
            });

            evaluation_function
        }
    };
}
