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
/// ## Questions
///
/// Question 1: Which operations should be supported? Are there any operations besides
/// multiplication and addition?
///   Answer: Only multiplication and addition
///
/// Question 2: Does multiplication take precedence over addition?
///   Answer: yes
///
/// Question 3: Should brackets be supported? (It would be simpler to not support them, such that
/// the definition of the virtual oracle is fully expanded.
///
///   e.g. the following should not be supported: [a(X) + b(X)][c(X)] but the following equivalent
///   representation should be supported: c(X) ⋅ a(X) + c(X) ⋅ b(X), even if c(X) is included twice
///   internally.
///
///   Answer: for now, just implement virtual oracles in their fully expanded form, without
///   brackets.
///
/// Question 4: Consider the following virtual oracle:
///
///   h(X) = f(αX) ⋅ g(X) + β
///
///   Should constants like β be part of the description (and therefore not hidden from the
///   verifier)?
///
///   Answer: yes.
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
/// The following is a design for the structs/traits I'll write:
///
/// UML:
///
/// EvaluableVirtualOracle (trait)
///   Attributes:
///     - description: Description
///     - concrete_oracles: Polynomial[]
///
///   Functions:
///     - new(Description) -> VirtualOracle
///       - returns a new VirtualOracle object
///
///     - instantiate(ConcreteOracle[]) -> Polynomial
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

pub struct Term<F: Field> {
    alpha_coeffs: Vec<F>,
    constants: Vec<F>
}

struct TermEvalError;

impl<F: Field> Term<F> {
    // Run by the prover
    fn evaluate_as_prover(
        &self,
        inputs: Vec<F>,
        concrete_oracles: Vec<DensePolynomial<F>>
    ) -> Result<F, TermEvalError> {

        // The length of inputs and concrete_oracles must equal the number of concrete oracles
        if inputs.len() != concrete_oracles.len() {
            return Err(TermEvalError);
        }

        // Sum the constants
        let sum_of_constants = self.constants.iter().fold(F::from(0 as u64), |sum, val| sum + val);

        // Compute the sum of the evaluations of each concrete oracle
        let mut sum_of_concrete_oracle_evaluations = F::from(0 as u64);
        for (i, c) in concrete_oracles.iter().enumerate() {
            let m = self.alpha_coeffs[i] * inputs[i];
            sum_of_concrete_oracle_evaluations += c.evaluate(&m);
        }

        return Ok(sum_of_concrete_oracle_evaluations + sum_of_constants);
    }
}

pub struct Description<F: Field> {
    terms: Vec<Term<F>>
}

pub struct InstantiationError;

pub trait EvaluableVirtualOracle<F: Field> {
    fn new(description: Description<F>) -> Self;
    fn instantiate(&self, concrete_oracles: Vec<DensePolynomial<F>>) -> Result<DensePolynomial<F>, InstantiationError>;
    fn query(&self, evaluations: Vec<F>) -> F;
}

impl<F: Field> EvaluableVirtualOracle<F> for TestVirtualOracle2<F> {
    fn new(description: Description<F>) -> Self {
        return TestVirtualOracle2::<F>{
            description: description,
            concrete_oracles: Vec::<DensePolynomial<F>>::new()
        };
    }

    fn instantiate(&self, concrete_oracles: Vec<DensePolynomial<F>>) -> Result<DensePolynomial<F>, InstantiationError> {
        let expected_num_concrete_oracles = self.description.terms.iter().fold(0 as usize, |sum, t| sum + t.alpha_coeffs.len());

        if concrete_oracles.len() != expected_num_concrete_oracles {
            return Err(InstantiationError);
        }

        // TODO WIP
        //// Construct the polynomial which represents the virtual oracle
        //let mut poly = DensePolynomial::<F>::from_coefficients_slice(&[]);
        //let mut i = 0 as usize;

        //// For each term
        //for (term_i, term) in self.description.terms.iter().enumerate() {
            //// Calculate the product of all concrete oracle in the term
            //let mut term_prod = DensePolynomial::<F>::from_coefficients_slice(&[]);

            //for _i in 0..term.alpha_coeffs.len() {
                //term_prod = term_prod.naive_mul(&concrete_oracles[i]);
            //}

            //// Add the constants
            //// i 
        //}
        return Err(InstantiationError);
    }

    fn query(&self, _evaluations: Vec<F>) -> F {
        // return a dummy value for now
        return F::from(0 as u64);
    }
}

pub struct TestVirtualOracle2<F: Field> {
    description: Description<F>,
    concrete_oracles: Vec<DensePolynomial<F>>
}

#[cfg(test)]
mod new_tests {
    use super::{TestVirtualOracle2, EvaluableVirtualOracle};
    use ark_bn254::Fr;

    #[test]
    fn create_and_evaluate_term() {
    }

    #[test]
    fn create_description() {
    }

    #[test]
    fn create_test_vo2() {
        let vo2 = TestVirtualOracle2::<Fr>::new(Fr::from(0 as u64));
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
