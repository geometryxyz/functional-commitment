use crate::to_poly;
use crate::util::shift_dense_poly;
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_poly_commit::LabeledPolynomial;
use ark_serialize::{CanonicalSerialize, SerializationError, Write};

/// TODO: rewrite these comments
/// A abstraction which represents the description of a virtual oracle. A virtual oracle is
/// composed of an indexed list of concrete oracles, the description of a virtual oracle which
/// specifies how said concrete oracles must be composed to form the virtual oracle, as well as the
/// alpha shifting coefficients for inputs to each concrete oracle. Following the definition of
/// virtual oracles on page 12 of the Efficient Functional Commitments paper, an example of a
/// virtual oracle h(X) is defined as such: h(X) = f(αX) ⋅ g(X) + β In this case, the concrete
/// oracles are f (at index 0) and g (at index 1), and β is a constant.  α is a shifting factor for
/// f and the shifting factor for g is (implicitly) 1.
///
/// ## Implementation notes
///
/// - Only multiplication and addition operations are supported in the representation of the virtual
///   oracle.
/// - Multiplication takes precedence over addition.
/// - Brackets are not supported, as it is much simpler to only represent the fully-expanded version
///   of the virtual oracle.
/// - Constants are part of the description, and are therefore visible to the verifier.
///
/// ## Who knows/runs what
///
/// In the protocol, the prover knows both the virtual oracle's description and the exact
/// definition of each concrete oracle, but the verifier only knows the description.
///
/// If the verifier knows the evaluations of each concrete oracle at X, as well as the description,
/// they can compute the evaluation of the virtual oracle at X, without needing to know the
/// definition of any of the concrete oracles.

/// A Term is part of a Description. It is composed of the product of one or more concrete
/// oracle, as well as a constant. e.g. if a virtual oracle is defined as:
///   F(X) = b(a_1X) + [2 ⋅ b(a_2X) ⋅ c(a_3X)] + d(a_4X)
/// then there are 3 terms: b, [2 ⋅ b ⋅ c], and d.
/// The constants are [1, 2, 1].
/// Given a list of concrete oracles [b, c, d], the concrete oracle indices for each term are: 
/// - [0]
/// - [0, 1]
/// - [2]
///
/// Given a list of alpha coefficients [a_1, a_2, a_3, a_4], the alpha coefficient indices for each
/// term are:
/// - [0]
/// - [1, 2]
/// - [3]
/// 
/// Note that this means that the virtual oracle is defined in disjunctive normal form (DNF).
///
/// Doing the above makes it easy to reuse the lists of alpha coefficients and concrete oracles, as
/// well as to determine the degree of the virtual oracle.
#[derive(Debug)]
pub struct Term<F> {
    pub concrete_oracle_indices: Vec<usize>,
    pub alpha_coeff_indices: Vec<usize>,
    pub constant: F
}

impl<F: Field> Term<F> {
    fn count_concrete_oracles(&self) -> usize {
        let mut indices = HashMap::<usize, bool>::new();
        for index in self.concrete_oracle_indices.iter(){
            if !indices.contains_key(index) {
                indices.insert(*index, true);
            }
        }
        indices.len()
    }
}

/// A Description is just a Vector of Terms, and a constant that is added to the end
#[derive(Debug)]
pub struct Description<F: Field> {
    pub terms: Vec<Term<F>>,
    pub constant: F
}

use std::collections::HashMap;
impl<F: Field> Description<F> {
    /// Returns the total number of concrete oracles across all Terms in a Description.
    /// If a Description contains the following Terms:
    /// - { concrete_oracle_indices: [0, 1] }
    /// - { concrete_oracle_indices: [1, 2] }
    /// Then the number of concrete oracles should be 3.
    fn count_concrete_oracles(&self) -> usize {
        let mut indices = HashMap::<usize, bool>::new();
        for term in self.terms.iter() {
            for index in term.concrete_oracle_indices.iter(){
                if !indices.contains_key(index) {
                    indices.insert(*index, true);
                }
            }
        }

        indices.len()
    }
}

#[derive(Debug)]
pub struct VirtualOracle2<F: Field> {
    description: Description<F>
}

pub trait EvaluationsProvider<F: Field> {
    fn evaluate(&self, virtual_oracle: &VirtualOracle2<F>, point: F, alpha_coeffs: &Vec<F>)  -> Result<F, EvaluationError>;
}

#[derive(Debug)]
pub struct InvalidDescriptionError;

#[derive(Debug)]
pub struct EvaluationError;

#[derive(Debug)]
pub struct InstantiationError;

impl<F: Field> VirtualOracle2<F> {
    pub fn new(description: Description<F>) -> Result<Self, InvalidDescriptionError> {
        if description.terms.len() == 0 {
            return Err(InvalidDescriptionError);
        }

        // Ensure that the lengths of the lists of indices are the same
        for term in description.terms.iter() {
            if term.concrete_oracle_indices.len() != term.alpha_coeff_indices.len() {
                return Err(InvalidDescriptionError);
            }
        }

        Ok(Self{ description: description })
    }

    /// Output a polynomial which represents the virtual oracle according to the description, the
    /// alpha coefficients, and the given concrete oracles. The alpha coefficients and concrete
    /// oracles should be arranged according to their respective lists of indices in the
    /// description.
    /// @param concrete_oracles A vector of concrete oracles whose order follows that of the
    ///                         concrete_oracle_indices in the Terms in the description
    pub fn instantiate(
        &self,
        concrete_oracles: &Vec<DensePolynomial<F>>,
        alpha_coeffs: &Vec<F>
    ) -> Result<DensePolynomial<F>, InstantiationError> {
        // Ensure that there are enough concrete oracles and alpha coefficients to fit the
        // description
        let num_cos = self.description.count_concrete_oracles();
        if num_cos < concrete_oracles.len() || num_cos < alpha_coeffs.len() {
            return Err(InstantiationError);
        }

        // get the max index from the flattened list of concrete oracle indices
        let max_co_index = self.description.terms.iter()
            .map(
                |term| term.concrete_oracle_indices.iter().max()
            )
            .max()
            .unwrap()
            .unwrap();

        // the given vector of concrete oracles must be large enough
        if max_co_index >= &concrete_oracles.len() {
            return Err(InstantiationError);
        }

        // get the max index from the flattened list of alpha coefficient indices
        let max_alpha_index = self.description.terms.iter()
            .map(
                |term| term.alpha_coeff_indices.iter().max()
            )
            .max()
            .unwrap()  // assume that the number of terms is > 0 
            .unwrap(); // and the number of alpha coefficients is > 0

        // the given vector of concrete oracles must be large enough
        if max_alpha_index >= &concrete_oracles.len() {
            return Err(InstantiationError);
        }

        let mut poly = DensePolynomial::<F>::default();

        // loop through each term
        for term in self.description.terms.iter() {
            // term_poly is the polynomial which represents the current term. its initial value is
            // the term constant which will be chained with the other concrete oracles in the term
            // via multiplication
            let mut term_poly = DensePolynomial::<F>::from_coefficients_slice(&[term.constant]);
            for (i, co_index) in term.concrete_oracle_indices.iter().enumerate() {
                // shift_dense_poly is needed to apply the alpha coefficient to each concrete
                // oracle, since resulting polynomial only has one "input" (F(X) rather than
                // the expanded form f(a_1 * X) + g(a_2 * X)
                let alpha_coeff = alpha_coeffs[i];
                term_poly = term_poly.naive_mul(&shift_dense_poly(&concrete_oracles[*co_index], &alpha_coeff));
            }
            poly = poly + term_poly;
        }

        Ok(poly + DensePolynomial::<F>::from_coefficients_slice(&[self.description.constant]))
    }
}

impl<F: Field> EvaluationsProvider<F> for Vec::<DensePolynomial<F>> {
    /// Instantiate and evaluate the virtual oracle. Returns an error if the length of
    /// concrete_oracles differs from the number of concrete oracles in the Description.
    fn evaluate(
        &self,
        virtual_oracle: &VirtualOracle2<F>,
        point: F,
        alpha_coeffs: &Vec<F>
    ) -> Result<F, EvaluationError> {
        let expected_num_concrete_oracles = virtual_oracle.description.count_concrete_oracles();

        if self.len() != expected_num_concrete_oracles {
            return Err(EvaluationError);
        }

        let poly = virtual_oracle.instantiate(&self, alpha_coeffs).unwrap();
        return Ok(poly.evaluate(&point));
    }
}

impl<F: Field> EvaluationsProvider<F> for Vec::<F> {
    /// Return the evaluation of the virtual oracle given a list of evaluations of each concrete
    /// oracle. The length of the input vector of evaluations must equal to the number of concrete
    /// oracles in the Description.
    fn evaluate(&self, virtual_oracle: &VirtualOracle2<F>, _: F, _: &Vec<F>) -> Result<F, EvaluationError> {
        let expected_num_concrete_oracles = virtual_oracle.description.count_concrete_oracles();

        if self.len() != expected_num_concrete_oracles {
            return Err(EvaluationError);
        }

        let mut total_eval = virtual_oracle.description.constant;

        // index to keep track of concrete_oracles
        let mut i = 0 as usize;

        // For each term
        for term in virtual_oracle.description.terms.iter() {
            let mut product_of_concrete_oracles = F::from(1 as u64);

            // Calculate the evaluation of the term
            for _ in 0..term.count_concrete_oracles() {
                product_of_concrete_oracles *= &self[i];
                i += 1;
            }

            total_eval += product_of_concrete_oracles;
        }

        Ok(total_eval)
    }
}

#[cfg(test)]
mod new_tests {
    use super::{Term, Description, VirtualOracle2, EvaluationsProvider};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use ark_bn254::Fr;
    use ark_ff::PrimeField;

    fn term() {
        let term = Term {
            concrete_oracle_indices: vec![0, 1],
            alpha_coeff_indices: vec![1, 2],
            constant: Fr::from(1 as u64)
        };
        assert_eq!(term.count_concrete_oracles(), 2);
    }

    fn description_case_0() -> Description<Fr> {
        // Let the virtual oracle be:
        // F(X) = [a(X)] + [b(X) ⋅ c(X)] + d(X) + 15
        //
        // Test case 1:
        //   a(X) = 0
        //   b(X) = 1
        //   c(X) = 2
        //   d(2X) = 3
        //   F(X) = (0) + (1 ⋅ 2) + 3 + 15 = 2 + 3 + 15 = 20

        // Define the terms

        // a(X) + 5
        let term0 = Term{
            concrete_oracle_indices: vec![0],
            alpha_coeff_indices: vec![0],
            constant: Fr::from(1u64),
        };

        // b(X) ⋅ c(X) + 10
        let term1 = Term{
            concrete_oracle_indices: vec![1, 2],
            alpha_coeff_indices: vec![1, 2],
            constant: Fr::from(1u64),
        };

        // d(X)
        let term2 = Term{
            concrete_oracle_indices: vec![3],
            alpha_coeff_indices: vec![3],
            constant: Fr::from(1u64),
        };

        let desc = Description::<Fr>{ terms: vec![term0, term1, term2], constant: Fr::from(15 as u64) };

        assert_eq!(desc.count_concrete_oracles(), 4 as usize);

        desc
    }

    // As the verifier, evaluate a virtual oracle
    #[test]
    fn verifier_evaluate_case_0() {
        let description = description_case_0();

        // Query the description
        let evaluations = vec![
            Fr::from(0 as u64),
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
        ];

        let vo = VirtualOracle2::new(description).unwrap();
        let result = evaluations.evaluate(&vo, Fr::default(), &Vec::<Fr>::default());
        assert!(result.is_ok());

        assert!(result.unwrap() == Fr::from(20 as u64));

        // Test EvaluationError by passing in an invalid number of evaluations
        let bad_evaluations = vec![
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
        ];
        let result = bad_evaluations.evaluate(&vo, Fr::default(), &Vec::<Fr>::default());
        assert!(result.is_err());
    }

    // As the prover, instantiate and evaluate a virtual oracle
    #[test]
    fn prover_evaluate_case_0() {
        // Let the virtual oracle be:
        // F(x) = [a(X)] + [b(X) ⋅ c(X)] + d(X) + 15
        //
        // Let the concrete oracles be:
        // a(X) = x
        // b(X) = x^2
        // c(X) = x^3
        // d(X) = x^4
        //
        // The constant is 15.
        //
        // The virtual oracle as a polynomial should be:
        // F(x) = [x] + [x^5] + [x^4] + 15
        let ax = DensePolynomial::<Fr>::from_coefficients_slice(&[
           Fr::from(0 as u64),
           Fr::from(1 as u64),
        ]);

        let bx = DensePolynomial::<Fr>::from_coefficients_slice(&[
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(1 as u64)
        ]);

        let cx = DensePolynomial::<Fr>::from_coefficients_slice(&[
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(1 as u64)
        ]);

        let dx = DensePolynomial::<Fr>::from_coefficients_slice(&[
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(0 as u64),
           Fr::from(1 as u64)
        ]);

        let description = description_case_0();
        let vo = VirtualOracle2::new(description).unwrap();

        let concrete_oracles = vec![ax, bx, cx, dx];
        let alpha_coeffs = vec![
            Fr::from(1u64),
            Fr::from(1u64),
            Fr::from(1u64),
            Fr::from(1u64),
        ];

        let p = vo.instantiate(&concrete_oracles, &alpha_coeffs);

        // F(1) = [2] + [2^5] + [2^4] + 15 = 65
        // Evaluate the polynomial at the point 2:
        let point = Fr::from(2 as u64);
        let eval = concrete_oracles.evaluate(&vo, point, &alpha_coeffs);

        assert_eq!(eval.unwrap().into_repr(), Fr::from(65 as u64).into_repr());
    }
}



///////////////////////////////////////////////////////////////////////////////
// The original code for VirtualOracle is below.

/// Encode a virtual oracle. `alphas` allow to keep track of the shift applied to the input point for each concrete oracle.
/// The `evaluate` function encodes the function `G` defined in the paper.
pub trait VirtualOracle<F: Field> {
    fn alphas(&self) -> Vec<F>;
    fn evaluate(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
        input: &F,
    ) -> F;
    fn query(concrete_oracle_evals: &[F]) -> F;

    fn instantiate(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> DensePolynomial<F>;
}

//this virtual oracle will encode f(alpha_1*x) + 2*g(alpha_2*x) - 2
#[derive(CanonicalSerialize)]
pub struct TestVirtualOracle<F: Field> {
    pub oracles: Vec<LabeledPolynomial<F, DensePolynomial<F>>>,
    pub alphas: Vec<F>,
}

impl<F: Field> VirtualOracle<F> for TestVirtualOracle<F> {
    fn alphas(&self) -> Vec<F> {
        self.alphas.to_vec()
    }

    fn evaluate(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
        input: &F,
    ) -> F {
        concrete_oracles[0].evaluate(&(alphas[0] * input))
            + F::from(2 as u64) * concrete_oracles[1].evaluate(&(alphas[1] * input))
            - F::from(2 as u64)
    }

    fn query(concrete_oracle_evals: &[F]) -> F {
        concrete_oracle_evals[0] + F::from(2 as u64) * concrete_oracle_evals[1] - F::from(2 as u64)
    }

    fn instantiate(
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> DensePolynomial<F> {
        shift_dense_poly(&concrete_oracles[0], &alphas[0])
            + &shift_dense_poly(&concrete_oracles[1], &alphas[1]) * F::from(2 as u64)
            + to_poly!(-F::from(2 as u64))
    }
}

#[cfg(test)]
mod tests {

    use super::{TestVirtualOracle, VirtualOracle};
    use ark_bn254::Fr;
    use ark_ff::One;
    use ark_poly::{univariate::DensePolynomial, UVPolynomial};
    use ark_poly_commit::LabeledPolynomial;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn create_virtual_oracle() {
        let rng = &mut test_rng();

        let f_coeffs = vec![
            Fr::from(-3i64),
            Fr::from(1u64),
            Fr::from(-2i64),
            Fr::from(1u64),
        ];
        let f_poly = DensePolynomial::from_coefficients_slice(&f_coeffs);
        let f_poly = LabeledPolynomial::new(String::from("f"), f_poly, None, None);

        let g_coeffs = vec![
            Fr::from(-2i64),
            Fr::from(0u64),
            Fr::from(1u64),
            Fr::from(0u64),
        ];
        let g_poly = DensePolynomial::from_coefficients_slice(&g_coeffs);
        let g_poly = LabeledPolynomial::new(String::from("g"), g_poly, None, None);

        let test_point = Fr::rand(rng);

        let f_part = f_poly.evaluate(&(Fr::from(3 as u64) * test_point));
        let g_part = g_poly.evaluate(&test_point);
        let virtual_oracle_eval_by_hand = f_part + Fr::from(2 as u64) * g_part - Fr::from(2 as u64);

        let alphas = &[Fr::from(3 as u64), Fr::one()];
        let oracles = &[f_poly, g_poly];
        let virtual_oracle_eval = TestVirtualOracle::<Fr>::evaluate(oracles, alphas, &test_point);

        let virtual_oracle_query_eval = TestVirtualOracle::<Fr>::query(&[f_part, g_part]);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_eval);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_query_eval);
    }
}
