use crate::error::Error;
use crate::to_poly;
use crate::util::shift_dense_poly;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::{PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::LabeledPolynomial;
use std::collections::HashMap;

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
pub struct NormalizedVirtualOracle<F: PrimeField> {
    description: Description<F>,
}

impl<F: PrimeField> NormalizedVirtualOracle<F> {
    // impl<F: PrimeField> NormalizedVirtualOracle<F> {
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

    /// Return the number of concrete oracles in this virtual oracle
    pub fn num_of_oracles(&self) -> usize {
        let mut indices = HashMap::<usize, bool>::new();
        for term in self.description.terms.iter() {
            for index in term.concrete_oracle_indices.iter() {
                if !indices.contains_key(index) {
                    indices.insert(*index, true);
                }
            }
        }
        indices.len()
    }
}

impl<F: PrimeField> VirtualOracle<F> for NormalizedVirtualOracle<F> {
    /// Output a polynomial which represents the virtual oracle according to the description, the
    /// alpha coefficients, and the given concrete oracles. The alpha coefficients and concrete
    /// oracles should be arranged according to their respective lists of indices in the
    /// description.
    /// @param concrete_oracles A vector of concrete oracles whose order follows that of the
    ///                         concrete_oracle_indices in the Terms in the description
    /// TODO: instead of concrete_oracles and alpha_coeffs as vector, expect slices of a specific
    /// size
    fn instantiate(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
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
        let max_co_index = self
            .description
            .terms
            .iter()
            .map(|term| term.concrete_oracle_indices.iter().max())
            .max()
            .unwrap()
            .unwrap();

        // the given vector of concrete oracles must be large enough
        if max_co_index >= &concrete_oracles.len() {
            return Err(Error::InstantiationError);
        }

        // get the max index from the flattened list of alpha coefficient indices
        let max_alpha_index = self
            .description
            .terms
            .iter()
            .map(|term| term.alpha_coeff_indices.iter().max())
            .max()
            .unwrap() // assume that the number of terms is > 0
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

    fn query(&self, evals: &[F], point: F) -> Result<F, Error> {
        // let mut total_eval = F::zero();

        // // For each term in the virtual oracle
        // for term in virtual_oracle.description.terms.iter() {
        //     let mut term_eval = term.constant;

        //     // for each concrete oracle in the term
        //     for &co_index in term.concrete_oracle_indices.iter() {
        //         term_eval *= self[co_index];
        //     }

        //     total_eval += term_eval;
        // }

        // Ok(total_eval)
        Ok(F::default())
    }

    fn num_of_oracles(&self) -> usize {
        self.num_of_oracles()
    }

    fn mapping_vector(&self) -> Vec<usize> {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::{Description, NormalizedVirtualOracle, Term};
    use crate::label_polynomial;
    use crate::virtual_oracle::{EvaluationsProvider, VirtualOracle};
    use ark_bn254::Fr;
    use ark_ff::PrimeField;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };

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
    fn verifier_evaluate_case_0() {
        let description = description_case_0();

        // Query the description
        let evaluations = vec![
            Fr::from(0 as u64),
            Fr::from(1 as u64),
            Fr::from(2 as u64),
            Fr::from(3 as u64),
        ];

        let vo = NormalizedVirtualOracle::new(description).unwrap();
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
        let vo = NormalizedVirtualOracle::new(description).unwrap();

        let concrete_oracles = vec![
            label_polynomial!(ax),
            label_polynomial!(bx),
            label_polynomial!(cx),
            label_polynomial!(dx),
        ];
        let alpha_coeffs = vec![Fr::from(1u64)];

        let point = Fr::from(2 as u64);

        // Test the instantiated polynomial by evaluating it at the point
        let p = vo.instantiate(&concrete_oracles, &alpha_coeffs).unwrap();
        assert_eq!(
            p.evaluate(&point).into_repr(),
            Fr::from(65 as u64).into_repr()
        );

        // F(2) = [2] + [2^5] + [2^4] + 15 = 65
        // Evaluate the polynomial at the point:
        let eval = concrete_oracles.evaluate(&vo, point, &alpha_coeffs);

        assert_eq!(eval.unwrap().into_repr(), Fr::from(65 as u64).into_repr());
    }
}
