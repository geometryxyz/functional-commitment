use crate::to_poly;
use crate::error::Error;
use crate::util::shift_dense_poly;
use ark_ff::{PrimeField, Zero};
use ark_poly::{
    univariate::{DensePolynomial},
    Polynomial, UVPolynomial,
};
use ark_poly_commit::{LabeledPolynomial, QuerySet};
use std::collections::HashMap;

pub trait VirtualOracleTrait<F: PrimeField> {
    /// abstract function which will return instantiation of virtual oracle from concrete oracles and shifting factors
    fn instantiate(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
        alphas: &[F],
    ) -> DensePolynomial<F>;

    /// from description of self get which concrete oracles should be queried and at which points
    /// points with same value should always have same label, for example alpha*beta = beta when alpha = 1
    /// that's why we gradually build point_values_to_labels from which we construct query set
    fn get_query_set(&self, labels: &Vec<String>, alphas: &Vec<(String, F)>, x: &(String, F)) -> QuerySet<F>;
}

pub struct GeoSequenceVO<F: PrimeField> {
    pi_s: Vec<u64>,
    ci_s: Vec<u64>,
    gamma: F,
    r: F,
}

impl<F: PrimeField> GeoSequenceVO<F> {
    pub fn new(pi_s: &Vec<u64>, ci_s: &Vec<u64>, gamma: F, r: F) -> Self {
        assert_eq!(pi_s.len(), ci_s.len());

        Self {
            pi_s: pi_s.clone(),
            ci_s: ci_s.clone(),
            gamma,
            r,
        }
    }
}

impl<F: PrimeField> VirtualOracleTrait<F> for GeoSequenceVO<F> {
    fn instantiate(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
        alphas: &[F],
    ) -> DensePolynomial<F> {
        assert_eq!(concrete_oracles.len(), 1);
        assert_eq!(alphas.len(), 1);

        // construct (f(gamma * x) - r * f(x))
        let mut instantiation_poly =
            shift_dense_poly(&concrete_oracles[0], &alphas[0]) + (&concrete_oracles[0] * -self.r);

        let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
            // construct x - y^(pi + ci - 1)
            let stitch_i = &x_poly + &to_poly!(-self.gamma.pow([(pi + ci - 1) as u64]));
            instantiation_poly = &instantiation_poly * &stitch_i;
        }

        instantiation_poly
    }

    fn get_query_set(&self, labels: &Vec<String>, alphas: &Vec<(String, F)>, x: &(String, F)) -> QuerySet<F> {

        assert_eq!(labels.len(), 1);
        assert_eq!(alphas.len(), 1);
        assert_eq!(alphas[0].1, self.gamma);

        let mut point_values_to_point_labels = HashMap::new();
        point_values_to_point_labels.insert(x.1, x.0.clone());

        let mut query_set = QuerySet::<F>::new();

        query_set.insert((labels[0].clone(), (x.0.clone(), x.1))); // h_prime_0 is evaluated at (beta_1 with value beta_1)
        
        let gamma_x = alphas[0].1 * x.1;
        let label = match point_values_to_point_labels.get(&gamma_x) {
            Some(label) => label.clone(),
            None => {
                let label = format!("{}_{}", alphas[0].0.clone(), x.0.clone());
                point_values_to_point_labels.insert(gamma_x, label.clone());

                label
            }
        };

        query_set.insert((labels[0].clone(), (label, gamma_x))); // h_prime_0 is evaluated at (alpha_0_beta1 with value gamma * beta_1)

        QuerySet::new()
    }
}
// for (i, alpha) in alphas.iter().enumerate() {
//     let test_point = *alpha * beta_1;
//     let label = match point_evaluations.get(&test_point) {
//         Some(label) => label.clone(),
//         None => {
//             let label = format!("alpha_{}_beta_1", i);
//             point_evaluations.insert(test_point, label.clone());

//             label
//         }
//     };
//     query_set.insert((format!("h_prime_{}", i), (label, *alpha * beta_1)));

//     let test_point = *alpha * beta_2;
//     let label = match point_evaluations.get(&test_point) {
//         Some(label) => label.clone(),
//         None => {
//             let label = format!("alpha_{}_beta_2", i);
//             point_evaluations.insert(test_point, label.clone());

//             label
//         }
//     };

//     query_set.insert((format!("m_{}", i), (label, *alpha * beta_2)));
// }

/// co[0] should be eval at alpha[0] * beta_1

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
///
/// TODO: have constant bounds instead of arbitrary vectors, to reduce overhead
#[derive(Debug)]
pub struct Term<F: PrimeField> {
    pub concrete_oracle_indices: Vec<usize>,
    pub alpha_coeff_indices: Vec<usize>,
    pub constant: F,
}

impl<F: PrimeField> Term<F> {
    fn count_concrete_oracles(&self) -> usize {
        let mut indices = HashMap::<usize, bool>::new();
        for index in self.concrete_oracle_indices.iter() {
            if !indices.contains_key(index) {
                indices.insert(*index, true);
            }
        }
        indices.len()
    }
}

/// A Description is just a Vector of Terms, and a constant that is added to the end
#[derive(Debug)]
pub struct Description<F: PrimeField> {
    pub terms: Vec<Term<F>>,
}

impl<F: PrimeField> Description<F> {
    /// Returns the total number of concrete oracles across all Terms in a Description.
    /// If a Description contains the following Terms:
    /// - { concrete_oracle_indices: [0, 1] }
    /// - { concrete_oracle_indices: [1, 2] }
    /// Then the number of concrete oracles should be 3.
    fn count_concrete_oracles(&self) -> usize {
        let mut indices = HashMap::<usize, bool>::new();
        for term in self.terms.iter() {
            for index in term.concrete_oracle_indices.iter() {
                if !indices.contains_key(index) {
                    indices.insert(*index, true);
                }
            }
        }

        indices.len()
    }
}

#[derive(Debug)]
pub struct VirtualOracle<F: PrimeField> {
    description: Description<F>,
}

pub trait EvaluationsProvider<F: PrimeField> {
    fn evaluate(
        &self,
        virtual_oracle: &VirtualOracle<F>,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error>;
}

impl<F: PrimeField> VirtualOracle<F> {
    pub fn new(description: Description<F>) -> Result<Self, Error> {
        if description.terms.len() == 0 {
            return Err(Error::InvalidDescriptionError);
        }

        // Ensure that the lengths of the lists of indices are the same
        for term in description.terms.iter() {
            if term.concrete_oracle_indices.len() != term.alpha_coeff_indices.len() {
                return Err(Error::InvalidDescriptionError);
            }
        }

        Ok(Self {
            description: description,
        })
    }

    /// Output a polynomial which represents the virtual oracle according to the description, the
    /// alpha coefficients, and the given concrete oracles. The alpha coefficients and concrete
    /// oracles should be arranged according to their respective lists of indices in the
    /// description.
    /// @param concrete_oracles A vector of concrete oracles whose order follows that of the
    ///                         concrete_oracle_indices in the Terms in the description
    /// TODO: instead of concrete_oracles and alpha_coeffs as vector, expect slices of a specific
    /// size
    pub fn instantiate(
        &self,
        concrete_oracles: &Vec<LabeledPolynomial<F, DensePolynomial<F>>>,
        alphas: &Vec<F>,
    ) -> Result<DensePolynomial<F>, Error> {
        // TODO: can there be a way to reduce the overhead of checking the indices?
        //
        //
        // Ensure that there are enough concrete oracles and alpha coefficients to fit the
        // description
         let num_cos = self.description.count_concrete_oracles();
         if num_cos < concrete_oracles.len() || num_cos < alphas.len() {
             return Err(Error::InstantiationError);
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
             return Err(Error::InstantiationError);
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
             return Err(Error::InstantiationError);
         }

        let mut poly = DensePolynomial::<F>::zero();

        // loop through each term
        for term in self.description.terms.iter() {
            // term_poly is the polynomial which represents the current term. its initial value is
            // the term constant which will be chained with the other concrete oracles in the term
            // via multiplication
            let mut term_poly = to_poly!(term.constant);
            for (&co_index, &alpha_index) in term
                .concrete_oracle_indices
                .iter()
                .zip(term.alpha_coeff_indices.iter())
            {
                let shifted_co =
                    shift_dense_poly(&concrete_oracles[co_index], &alphas[alpha_index]);
                term_poly = &term_poly * &shifted_co;
            }
            poly = poly + term_poly;
        }

        Ok(poly)
    }
}

impl<F: PrimeField> EvaluationsProvider<F> for Vec<LabeledPolynomial<F, DensePolynomial<F>>> {
    /// Instantiate and evaluate the virtual oracle. Returns an error if the length of
    /// concrete_oracles differs from the number of concrete oracles in the Description.
    fn evaluate(
        &self,
        virtual_oracle: &VirtualOracle<F>,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error> {
        let expected_num_concrete_oracles = virtual_oracle.description.count_concrete_oracles();

        if self.len() != expected_num_concrete_oracles {
            return Err(Error::EvaluationError);
        }

        let poly = virtual_oracle.instantiate(&self, alpha_coeffs).unwrap();
        return Ok(poly.evaluate(&point));
    }
}

impl<F: PrimeField> EvaluationsProvider<F> for Vec<F> {
    /// Return the evaluation of the virtual oracle given a list of evaluations of each concrete
    /// oracle. The length of the input vector of evaluations must equal to the number of concrete
    /// oracles in the Description.
    fn evaluate(
        &self,
        virtual_oracle: &VirtualOracle<F>,
        _: F,
        _: &Vec<F>,
    ) -> Result<F, Error> {
        let expected_num_concrete_oracles = virtual_oracle.description.count_concrete_oracles();

        if self.len() != expected_num_concrete_oracles {
            return Err(Error::EvaluationError);
        }

        let mut total_eval = F::zero();

        // For each term in the virtual oracle
        for term in virtual_oracle.description.terms.iter() {
            let mut term_eval = term.constant;

            // for each concrete oracle in the term
            for &co_index in term.concrete_oracle_indices.iter() {
                term_eval *= self[co_index];
            }

            total_eval += term_eval;
        }

        Ok(total_eval)
    }
}

// co = vec![f(x), g(x)] labeled
// alphas = vec![1, gamma]
// 1) create the identity polynomial and insert it to co
// co = vec![x, f(x), g(x)]
// j = [1, 1] means h_0 = co[1] = f, h_1 = co[1] = f
/// f(gamma * x) - f(x)
// VO = shift(co[j[1]], alpha)
// V0 = shift(co[1], alpha)

/// term1 = co_indices = [1] alphas_indicies = [1] c = 1
/// term2 = co_indices = [1] alphas_indicies = [0] c = -1
///
/// query set
/// h_prime_commitments = [h_prime_commitment_0] = (commitment_f + commiement_m0)
///
///
///
/// if i want vo to give me query_set
/// co_0 shoud be evaluated at point alpha_1 * x
/// (H_PRIME_0, (alpha_{}_x_label_we_sent, alpha_1 * x))
/// alpha_3 * beta_1 == alpha_7 * beta_2
/// F

// function that says :
// co[j[]]

#[cfg(test)]
mod new_tests {
    use super::{Description, EvaluationsProvider, Term, VirtualOracle};
    use crate::label_polynomial;
    use ark_bn254::Fr;
    use ark_ff::PrimeField;
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};

    #[test]
    fn term() {
        let term = Term {
            concrete_oracle_indices: vec![0, 1],
            alpha_coeff_indices: vec![1, 2],
            constant: Fr::from(1 as u64),
        };
        assert_eq!(term.count_concrete_oracles(), 2);
    }

    fn description_case_0() -> Description<Fr> {
        // Let the virtual oracle be:
        // F(X) = [a(alpha_1X)] + [b(alpha_1X) ⋅ c(alpha_1X)] + [d(alpha_1X)] + 15
        //
        // Test case 1:
        //   a(X) = 0
        //   b(X) = 1
        //   c(X) = 2
        //   d(2X) = 3
        //   F(X) = (0) + (1 ⋅ 2) + 3 + 15 = 2 + 3 + 15 = 20

        // Define the terms

        // a(X)
        let term0 = Term {
            concrete_oracle_indices: vec![0],
            alpha_coeff_indices: vec![0],
            constant: Fr::from(1u64),
        };

        // b(X) ⋅ c(X)
        let term1 = Term {
            concrete_oracle_indices: vec![1, 2],
            alpha_coeff_indices: vec![0, 0],
            constant: Fr::from(1u64),
        };

        // d(X)
        let term2 = Term {
            concrete_oracle_indices: vec![3],
            alpha_coeff_indices: vec![0],
            constant: Fr::from(1u64),
        };

        // const term
        let term3 = Term {
            concrete_oracle_indices: vec![],
            alpha_coeff_indices: vec![],
            constant: Fr::from(15u64),
        };

        let desc = Description::<Fr> {
            terms: vec![term0, term1, term2, term3],
        };

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

        let vo = VirtualOracle::new(description).unwrap();
        let result = evaluations.evaluate(&vo, Fr::default(), &Vec::<Fr>::default());
        assert!(result.is_ok());

        assert!(result.unwrap() == Fr::from(20 as u64));

        // Test EvaluationError by passing in an invalid number of evaluations
        let bad_evaluations = vec![Fr::from(1 as u64), Fr::from(2 as u64), Fr::from(3 as u64)];
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
            Fr::from(1 as u64),
        ]);

        let cx = DensePolynomial::<Fr>::from_coefficients_slice(&[
            Fr::from(0 as u64),
            Fr::from(0 as u64),
            Fr::from(0 as u64),
            Fr::from(1 as u64),
        ]);

        let dx = DensePolynomial::<Fr>::from_coefficients_slice(&[
            Fr::from(0 as u64),
            Fr::from(0 as u64),
            Fr::from(0 as u64),
            Fr::from(0 as u64),
            Fr::from(1 as u64),
        ]);

        let description = description_case_0();
        let vo = VirtualOracle::new(description).unwrap();

        let concrete_oracles = vec![
            label_polynomial!(ax),
            label_polynomial!(bx),
            label_polynomial!(cx),
            label_polynomial!(dx),
        ];
        let alpha_coeffs = vec![Fr::from(1u64)];

        let point = Fr::from(2 as u64);

        // Test the instantiated polynomial by evaluating it at the point
        let p = vo.instantiate(&concrete_oracles, &alpha_coeffs);
        assert_eq!(p.unwrap().evaluate(&point).into_repr(), Fr::from(65 as u64).into_repr());

        // F(2) = [2] + [2^5] + [2^4] + 15 = 65
        // Evaluate the polynomial at the point:
        let eval = concrete_oracles.evaluate(&vo, point, &alpha_coeffs);

        assert_eq!(eval.unwrap().into_repr(), Fr::from(65 as u64).into_repr());
    }
}
