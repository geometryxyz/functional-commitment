use crate::to_poly;
use crate::util::shift_dense_poly;
use ark_bn254::Fr;
use ark_ff::{Field, One};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_std::{test_rng, UniformRand};

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
///
/// ## Next steps
///
/// The following is an UML-esque design for the structs/traits:
///
/// VirtualOracle2 (struct) TODO: rename this to VirtualOracle
///   Attributes:
///     - description: Description
///   Functions (as impl):
///     - new(Description) -> VirtualOracle
///       - returns a new VirtualOracle object
///
/// EvaluationsProvider (trait)
///   - A common interface for the prover and verifier to evaluate a virtual oracle.
///   Functions:
///     - instantiate(Polynomial[]) -> Polynomial
///       - returns a Polynomial which exactly defines the virtual oracle. Only the prover may use
///         this function.
///       - the order of the Polynomials in concrete_oracle is crucial. See below.
///
///     - query(Scalar[]) -> Scalar
///       - Given a list of evaluations of concrete oracles, returns a value that is the evaluation
///         of the virtual oracle. The verifier may use this function.
///    
///
/// Description (struct)
///   Attributes:
///     - terms: Term[]
///
///   Functions:
///     - group_concrete_oracles_by_alphas: Dict{F -> integer[]} (TBD?)
///       - Returns an associative array where the index of the concrete oracles are grouped by
///         alpha coefficients. Useful for batch polynomial commitments?
///
///  Term (struct)
///   Attributes:
///     - alpha_coeffs: Scalar[]
///     - constants: Scalar[]
///
///   Functions:
///     - evaluate_as_prover(inputs: Scalar[], concrete_oracles: Polynomial[]) -> Scalar (or Error)
///
///   Notes:
///     - the length of alpha_coeffs must be exactly the number of concrete oracles in this term
///
/// ## Important notes
///
/// The order of concrete oracles passsed to instantiate() is critical. Take for example the
/// following virtual oracle:
/// h(X) = f(αX) ⋅ g(X) + β0 + a(X) ⋅ b(X) ⋅ c(X) + β1
///
/// h(X) contains two terms:
///   - Term 0:
///     - 2 concrete oracles (f and g)
///     - 1 constant (β0)
///   - Term 1:
///     - 3 concrete oracles (a, b, and c)
///     - 1 constant (β1)
///
/// The order of concrete oracles passed to instantiate() should therefore be [f, g, a, b c]

/// A term is part of a virtual oracle. It is composed of the product of one or more concrete
/// oracle, plus the sum of some constants. e.g. if a virtual oracle is defined as:
pub struct Term<F: Field> {
    alpha_coeffs: Vec<F>,
    constants: Vec<F>
}

/// A Description is just a Vector of Terms
pub struct Description<F: Field> {
    terms: Vec<Term<F>>
}

impl<F: Field> Description<F> {
    /// Returns the total number of concrete oracles across all Terms in a Description.
    fn count_concrete_oracles(&self) -> usize {
        let r = self.terms.iter().fold(0 as usize, |sum, t| sum + t.alpha_coeffs.len());
        return r;
    }
}

pub struct VirtualOracle2<F: Field> {
    description: Description<F>
}

pub trait EvaluationsProvider<F: Field> {
    fn instantiate(&self, concrete_oracles: Vec<DensePolynomial<F>>) -> Result<DensePolynomial<F>, InstantiationError>;
    fn query(&self, evaluations: Vec<F>) -> Result<F, QueryError>;
}

pub struct InstantiationError;
pub struct QueryError;

impl<F: Field> VirtualOracle2<F> {
    fn new(description: Description<F>) -> Self {
        return Self{
            description: description
        };
    }
}

impl<F: Field> EvaluationsProvider<F> for VirtualOracle2<F> {
    /// Construct the polynomial which represents the virtual oracle. Returns an error if the
    /// length of concrete_oracles differs from the number of concrete oracles in the Description.
    fn instantiate(&self, concrete_oracles: Vec<DensePolynomial<F>>) -> Result<DensePolynomial<F>, InstantiationError> {
        let expected_num_concrete_oracles = self.description.count_concrete_oracles();

        if concrete_oracles.len() != expected_num_concrete_oracles {
            return Err(InstantiationError);
        }

        // the polynomial to construct and return
        let mut poly = DensePolynomial::<F>::default();

        // index to keep track of concrete_oracles
        let mut i = 0 as usize;

        // For each term
        for term in self.description.terms.iter() {
            // Sum the constants
            let constant_sum = term.constants.iter().fold(F::from(0 as u64), |sum, c| sum + c);

            // term_product is the 1-degree polynomial: f(X) = (constant_sum ⋅ x^0)
            let mut term_product: DensePolynomial<F>;

            if term.alpha_coeffs.len() == 0 {
                // An empty polynomial
                term_product = DensePolynomial::<F>::default();
            } else {
                // Calculate the product of all concrete oracles in the term
                term_product = concrete_oracles[i].clone();
                i += 1;

                for _i in 1..term.alpha_coeffs.len() {
                    term_product = term_product.naive_mul(&concrete_oracles[i]);
                    i += 1;
                }
            }

            poly = poly + term_product + DensePolynomial::<F>::from_coefficients_slice(&[constant_sum]);
        }

        return Ok(poly);
    }

    /// Return the evaluation of the virtual oracle given a list of evaluations of each concrete
    /// oracle. The length of the input vector of evaluations must equal to the number of concrete
    /// oracles in the Description.
    fn query(&self, evaluations: Vec<F>) -> Result<F, QueryError> {
        let expected_num_concrete_oracles = self.description.count_concrete_oracles();

        if evaluations.len() != expected_num_concrete_oracles {
            return Err(QueryError);
        }

        let mut total_eval = F::from(0 as u64);

        // index to keep track of concrete_oracles
        let mut i = 0 as usize;

        // For each term
        for term in self.description.terms.iter() {
            // Sum the constants
            let sum_of_constants = term.constants.iter().fold(F::from(0 as u64), |sum, c| sum + c);
            let mut product_of_concrete_oracles = F::from(1 as u64);

            // Calculate the evaluation of the term
            for _i in 0..term.alpha_coeffs.len() {
                product_of_concrete_oracles *= &evaluations[i];
                i += 1;
            }

            total_eval += product_of_concrete_oracles + sum_of_constants;
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

    fn description_case_0() -> Description<Fr> {
        // Let the virtual oracle be:
        // F(X) = [a(X) + 5] + [b(X) ⋅ c(X) + 10] + d(2X)
        //
        // Test case 1:
        //   a(X) = 0
        //   b(X) = 1
        //   c(X) = 2
        //   d(2X) = 3
        //   F(X) = (0 + 5) + (1 ⋅ 2 + 10) + 3 = 5 + 12 + 3 = 20

        // Define the terms

        // a(X) + 5
        let term0 = Term::<Fr>{
            alpha_coeffs: vec![Fr::from(1 as u64)],
            constants: vec![Fr::from(5 as u64)],
        };

        // b(X) ⋅ c(X) + 10
        let term1 = Term::<Fr>{
            alpha_coeffs: vec![Fr::from(1 as u64), Fr::from(1 as u64)],
            constants: vec![Fr::from(10 as u64)],
        };

        // d(2X)
        let term2 = Term::<Fr>{
            alpha_coeffs: vec![Fr::from(2 as u64)],
            constants: vec![Fr::from(0 as u64)],
        };

        // Define the description
        Description::<Fr>{ terms: vec![term0, term1, term2] }
    }

    // As the verifier, create and query a description
    #[test]
    fn verify_create_and_query_description() {
        let description = description_case_0();

        assert!(description.count_concrete_oracles() == 4 as usize);

        // Query the description
        let evaluations = vec![
            Fr::from(0 as u64),
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
        ];

        let vo = VirtualOracle2::new(description);
        let result = vo.query(evaluations);
        assert!(result.is_ok());

        // TODO: figure out why result.unwrap() returns QueryError() even though it is Ok
        //assert!(result.unwrap() == Fr::from(20 as u64));

        match result {
            Ok(r) => {
                assert!(r == Fr::from(20 as u64));
            }
            Err(e) => {
                // This shouldn't happen
                assert_eq!(true, false);
            }
        }
        
        // Test QueryError by passing in an invalid number of evaluations
        let bad_evaluations = vec![
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
        ];
        let result = vo.query(bad_evaluations);
        assert!(result.is_err());
    }

    // As the prover, instantiate and evaluate a virtual oracle
    #[test]
    fn prover_instantiate_and_evaluate() {
        // Let the virtual oracle be:
        // F(x) = [a(X) + 5] + [b(X) ⋅ c(X) + 10] + d(2X)
        //
        // Let the concrete oracles be:
        // a(X) = x
        // b(X) = x^2
        // c(X) = x^3
        // d(X) = x^4
        //
        // The virtual oracle as a polynomial should be:
        // F(x) = [x + 5] + [x^5 + 10] + [x^4]
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
        let vo = VirtualOracle2::new(description);

        let result = vo.instantiate(vec![ax, bx, cx, dx]);
        assert!(result.is_ok());

        match result {
            Ok(r) => {
                // Evaluate the polynomial at the point 1:
                // F(1) = [1 + 5] + [1 + 10] + [1] = 6 + 11 + 1 = 18
                let eval = r.evaluate(&Fr::from(1 as u64));

                assert_eq!(eval.into_repr(), Fr::from(18 as u64).into_repr());
            }
            Err(e) => {
                // This shouldn't happen
                assert_eq!(true, false);
            }
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
// The original code for VirtualOracle is below.

/// Encode a virtual oracle. `alphas` allow to keep track of the shift applied to the input point for each concrete oracle.
/// The `evaluate` function encodes the function `G` defined in the paper.
pub trait VirtualOracle<F: Field, const NUM_OF_CONCRETE_ORACLES: usize> {
    fn evaluate(
        concrete_oracles: &[DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
        input: &F,
    ) -> F;
    fn query(concrete_oracle_evals: &[F; NUM_OF_CONCRETE_ORACLES]) -> F;

    fn instantiate(
        concrete_oracles: &[DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
    ) -> DensePolynomial<F>;
}

//this virtual oracle will encode f(alpha_1*x) + 2*g(alpha_2*x) - 2
pub struct TestVirtualOracle<F: Field, const NUM_OF_CONCRETE_ORACLES: usize> {
    oracles: [DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
    alphas: [F; NUM_OF_CONCRETE_ORACLES],
}

impl<F: Field, const NUM_OF_CONCRETE_ORACLES: usize> VirtualOracle<F, NUM_OF_CONCRETE_ORACLES>
    for TestVirtualOracle<F, NUM_OF_CONCRETE_ORACLES>
{
    fn evaluate(
        concrete_oracles: &[DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
        input: &F,
    ) -> F {
        concrete_oracles[0].evaluate(&(alphas[0] * input))
            + F::from(2 as u64) * concrete_oracles[1].evaluate(&(alphas[1] * input))
            - F::from(2 as u64)
    }

    fn query(concrete_oracle_evals: &[F; NUM_OF_CONCRETE_ORACLES]) -> F {
        concrete_oracle_evals[0] + F::from(2 as u64) * concrete_oracle_evals[1] - F::from(2 as u64)
    }

    fn instantiate(
        concrete_oracles: &[DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
    ) -> DensePolynomial<F> {
        shift_dense_poly(&concrete_oracles[0], &alphas[0])
            + &shift_dense_poly(&concrete_oracles[1], &alphas[1]) * F::from(2 as u64)
            + to_poly!(-F::from(2 as u64))
    }
}

#[cfg(test)]
mod tests {

    use super::{TestVirtualOracle, VirtualOracle};
    use crate::to_poly;
    use crate::util::shift_dense_poly;
    use ark_bn254::Fr;
    use ark_ff::{Field, One};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
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

        let g_coeffs = vec![
            Fr::from(-2i64),
            Fr::from(0u64),
            Fr::from(1u64),
            Fr::from(0u64),
        ];
        let g_poly = DensePolynomial::from_coefficients_slice(&g_coeffs);

        let test_point = Fr::rand(rng);

        let f_part = f_poly.evaluate(&(Fr::from(3 as u64) * test_point));
        let g_part = g_poly.evaluate(&test_point);
        let virtual_oracle_eval_by_hand = f_part + Fr::from(2 as u64) * g_part - Fr::from(2 as u64);

        let alphas = &[Fr::from(3 as u64), Fr::one()];
        let oracles = &[f_poly, g_poly];
        let virtual_oracle_eval =
            TestVirtualOracle::<Fr, 2>::evaluate(oracles, alphas, &test_point);

        let virtual_oracle_query_eval = TestVirtualOracle::<Fr, 2>::query(&[f_part, g_part]);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_eval);
        assert_eq!(virtual_oracle_eval_by_hand, virtual_oracle_query_eval);
    }
}
