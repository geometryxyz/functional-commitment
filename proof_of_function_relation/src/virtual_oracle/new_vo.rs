use std::ops::{Add, Div, Mul, Sub};

use crate::error::Error;
use crate::util::shift_dense_poly;
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_poly_commit::{Evaluations, PolynomialLabel, QuerySet};

pub struct NewVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    mapping_vector: Vec<usize>,
    shifting_coefficients: Vec<F>,
    combine_function: T,
    minimum_oracle_length: usize,
}

impl<F, T> NewVO<F, T>
where
    F: Field,
    T: Fn(&[VOTerm<F>]) -> VOTerm<F>,
{
    /// Constructor for an input-shifting virtual oracle
    pub fn new(
        mapping_vector: Vec<usize>,
        shifting_coefficients: Vec<F>,
        combine_function: T,
    ) -> Result<Self, Error> {
        let number_of_terms = mapping_vector.len();

        if shifting_coefficients.len() != number_of_terms {
            return Err(Error::InputLengthError(String::from(
                "mapping vector and shifting coefficients do not match",
            )));
        }

        let max_index = mapping_vector
            .iter()
            .max()
            .expect("mapping vector is empty")
            .clone();

        let minimum_oracle_length = max_index + 1;

        Ok(Self {
            mapping_vector,
            shifting_coefficients,
            combine_function,
            minimum_oracle_length,
        })
    }

    pub fn number_of_internal_terms(&self) -> usize {
        self.mapping_vector.len()
    }

    pub fn get_term_labels(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
    ) -> Vec<PolynomialLabel> {
        self.mapping_vector
            .iter()
            .map(|&mapped_index| concrete_oracle_labels[mapped_index].clone())
            .collect()
    }

    /// Returns the polynomial that results from the combination of the given concrete oracles
    pub fn compute_polynomial(
        &self,
        concrete_oracles: &[DensePolynomial<F>],
    ) -> Result<DensePolynomial<F>, Error> {
        self.check_conrete_oracle_length(concrete_oracles.len())?;

        let mut terms: Vec<VOTerm<F>> = Vec::new();

        // For each item in the mapping vector, we select the corresponding concrete oracle, apply the desired
        // shift and push the resulting polynomial as a term.
        self.mapping_vector
            .iter()
            .enumerate()
            .for_each(|(term_index, &mapped_index)| {
                let shifted = shift_dense_poly(
                    &concrete_oracles[mapped_index],
                    &self.shifting_coefficients[term_index],
                );
                terms.push(VOTerm::Polynomial(shifted))
            });

        let combined = (self.combine_function)(&terms);
        match combined {
            VOTerm::Evaluation(_) => Err(Error::VOFailedToInstantiate),
            VOTerm::Polynomial(poly) => Ok(poly),
        }
    }

    pub fn query(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        labeled_point: &(String, F),
    ) -> Result<QuerySet<F>, Error> {
        self.check_conrete_oracle_length(concrete_oracle_labels.len())?;

        let mut query_set = QuerySet::new();

        self.mapping_vector
            .iter()
            .enumerate()
            .for_each(|(term_index, &mapped_index)| {
                let poly_label = concrete_oracle_labels[mapped_index].clone();
                let eval_point = self.shifting_coefficients[term_index] * labeled_point.1;
                let point_label = format!("{}_times_alpha{}", labeled_point.0, term_index);

                query_set.insert((poly_label, (point_label, eval_point)));
            });

        Ok(query_set)
    }

    /// Given evalutations of each of the concrete oracles, produce the corresponding evaluation for the virtual oracle
    pub fn evaluate_from_concrete_evals(
        &self,
        concrete_oracle_labels: &[PolynomialLabel],
        eval_point: &F,
        evaluations: &Evaluations<F, F>,
    ) -> Result<F, Error> {
        let terms: Vec<VOTerm<_>> = self
            .mapping_vector
            .iter()
            .enumerate()
            .map(|(term_index, &mapped_index)| {
                let poly_label = concrete_oracle_labels[mapped_index].clone();
                let shifted_eval_point = self.shifting_coefficients[term_index] * eval_point;
                let key = (poly_label, shifted_eval_point);

                VOTerm::Evaluation(
                    evaluations
                        .get(&key)
                        .expect("Missing a concrete oracle evaluation for VO computation")
                        .clone(),
                )
            })
            .collect();

        let combined = (self.combine_function)(&terms);
        match combined {
            VOTerm::Evaluation(eval) => Ok(eval),
            VOTerm::Polynomial(_) => Err(Error::VOFailedToCompute),
        }
    }

    /// Check that enough oracles were provided.
    fn check_conrete_oracle_length(&self, input_length: usize) -> Result<(), Error> {
        if input_length < self.minimum_oracle_length {
            return Err(Error::InputLengthError(format!(
                "Mapping vector requires {} oracles/evaluations but only {} were provided",
                self.minimum_oracle_length, input_length
            )));
        }
        Ok(())
    }
}

/// A term to be manipulated by the virtual oracle's function, this is either a (shifted) polynomial or an evaluation thereof.
#[derive(Clone)]
pub enum VOTerm<F: Field> {
    Evaluation(F),
    Polynomial(DensePolynomial<F>),
}

impl<F: Field> Add for VOTerm<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval + rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(DensePolynomial::from_coefficients_slice(&[eval]) + rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(poly + DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(poly + rhs_poly),
            },
        }
    }
}

impl<F: Field> Sub for VOTerm<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval - rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(&DensePolynomial::from_coefficients_slice(&[eval]) - &rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(&poly - &DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly - &rhs_poly),
            },
        }
    }
}

impl<F: FftField> Mul for VOTerm<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval * rhs_eval),
                Self::Polynomial(rhs_poly) => Self::Polynomial(&rhs_poly * eval),
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Polynomial(&poly * rhs_eval),
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly * &rhs_poly),
            },
        }
    }
}

impl<F: FftField> Div for VOTerm<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        match self {
            Self::Evaluation(eval) => match rhs {
                Self::Evaluation(rhs_eval) => Self::Evaluation(eval / rhs_eval),
                Self::Polynomial(rhs_poly) => {
                    Self::Polynomial(&DensePolynomial::from_coefficients_slice(&[eval]) / &rhs_poly)
                }
            },
            Self::Polynomial(poly) => match rhs {
                Self::Evaluation(rhs_eval) => {
                    Self::Polynomial(&poly / &DensePolynomial::from_coefficients_slice(&[rhs_eval]))
                }
                Self::Polynomial(rhs_poly) => Self::Polynomial(&poly / &rhs_poly),
            },
        }
    }
}

#[macro_export]
macro_rules! vo_constant {
    ($field_element:expr) => {
        VOTerm::Evaluation($field_element)
    };
}

#[cfg(test)]
mod test {
    use crate::util::sample_vector;
    use crate::virtual_oracle::geometric_sequence_vo::GeoSequenceVO;
    use crate::virtual_oracle::VirtualOracle;
    use crate::{error::Error, util::generate_sequence, util::shift_dense_poly, vo_constant};
    use ark_bn254::Fr;
    use ark_ff::{Field, One, UniformRand, Zero};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_poly_commit::{evaluate_query_set, LabeledPolynomial};
    use rand::thread_rng;
    use std::iter;

    use super::{NewVO, VOTerm};

    type F = Fr;

    pub fn simple_addition(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
        concrete_terms[0].clone() + vo_constant!(F::from(2u64)) * concrete_terms[1].clone()
    }

    pub fn simple_mul(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
        concrete_terms[0].clone() * concrete_terms[1].clone()
    }

    pub fn harder_addition(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
        concrete_terms[0].clone()
            + vo_constant!(F::from(2u64)) * concrete_terms[1].clone()
            + vo_constant!(F::from(3u64)) * concrete_terms[2].clone()
    }

    #[test]
    fn test_add_oracle() {
        let rng = &mut thread_rng();
        let a_poly = DensePolynomial::<F>::rand(4, rng);
        let b_poly = DensePolynomial::<F>::rand(4, rng);
        let c_poly = DensePolynomial::<F>::rand(4, rng);
        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let two_constant_poly = DensePolynomial::from_coefficients_vec(vec![F::from(2u64)]);
        let expected = &shifted_c + &(&two_constant_poly * &shifted_a);

        let add_oracle =
            NewVO::new(mapping_vector, shifting_coefficients, simple_addition).unwrap();

        // Check that we get the right polynomial
        let sum = add_oracle.compute_polynomial(concrete_oracles).unwrap();
        assert_eq!(expected, sum);
    }

    #[test]
    fn test_mul_oracle() {
        let rng = &mut thread_rng();
        let a_poly = DensePolynomial::<F>::rand(4, rng);
        let b_poly = DensePolynomial::<F>::rand(4, rng);
        let c_poly = DensePolynomial::<F>::rand(4, rng);
        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];

        let shifting_coefficients: Vec<F> = sample_vector(rng, 2);
        let mapping_vector = vec![2, 0];

        // Compute the expected polynomial
        let shifted_c = shift_dense_poly(&c_poly, &shifting_coefficients[0]);
        let shifted_a = shift_dense_poly(&a_poly, &shifting_coefficients[1]);
        let expected = &shifted_c * &shifted_a;

        let mul_oracle = NewVO::new(mapping_vector, shifting_coefficients, simple_mul).unwrap();

        // Check that we get the right polynomial
        let prod = mul_oracle.compute_polynomial(concrete_oracles).unwrap();
        assert_eq!(expected, prod);
    }

    #[test]
    fn test_short_input_vec() {
        // mapping vector expects there to be a concrete oracle with index 1; effectively expected at last 2 concrete oracles
        let mapping_vector = vec![1];
        let shift_coefficients = vec![F::one()];
        let add_oracle = NewVO::new(mapping_vector, shift_coefficients, simple_addition).unwrap();

        // We only provide one concrete oracle
        let err_poly = add_oracle.compute_polynomial(&vec![DensePolynomial::<F>::default()]);
        assert!(err_poly.is_err());
        assert_eq!(
            err_poly.unwrap_err(),
            Error::InputLengthError(String::from(
                "Mapping vector requires 2 oracles/evaluations but only 1 were provided"
            ))
        );
    }

    #[test]
    fn test_query_set() {
        let rng = &mut thread_rng();
        let a_poly = LabeledPolynomial::new(
            String::from("a"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let b_poly = LabeledPolynomial::new(
            String::from("b"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );
        let c_poly = LabeledPolynomial::new(
            String::from("c"),
            DensePolynomial::<F>::rand(4, rng),
            None,
            None,
        );

        let concrete_oracles = &[a_poly.clone(), b_poly.clone(), c_poly.clone()];
        let oracle_labels: Vec<_> = concrete_oracles
            .iter()
            .map(|oracle| oracle.label().clone())
            .collect();
        let oracle_polys: Vec<_> = concrete_oracles
            .iter()
            .map(|oracle| oracle.polynomial().clone())
            .collect();

        let shifting_coefficients: Vec<F> = sample_vector(rng, 3);
        let mapping_vector = vec![2, 2, 0];

        let add_oracle =
            NewVO::new(mapping_vector, shifting_coefficients, harder_addition).unwrap();

        let eval_point = (String::from("beta"), F::rand(rng));
        let query_set = add_oracle.query(&oracle_labels, &eval_point).unwrap();

        let evals = evaluate_query_set(concrete_oracles.iter(), &query_set);

        let evaluated = add_oracle
            .evaluate_from_concrete_evals(&oracle_labels, &eval_point.1, &evals)
            .unwrap();
        let eval_from_poly = add_oracle
            .compute_polynomial(&oracle_polys)
            .unwrap()
            .evaluate(&eval_point.1);

        assert_eq!(evaluated, eval_from_poly)
    }

    #[test]
    fn test_geo_seq() {
        let rng = &mut thread_rng();

        let common_ratio = Fr::rand(rng);
        let mut initial_values = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut sequence_lengths = vec![5, 3, 10, 30];

        let m = sequence_lengths.iter().sum();

        let domain = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let to_pad = domain.size() - m;
        if to_pad > 0 {
            initial_values.push(Fr::from(0u64));
            sequence_lengths.push(to_pad);
        }

        let seq = generate_sequence::<Fr>(common_ratio, &initial_values, &sequence_lengths);
        let f = DensePolynomial::<Fr>::from_coefficients_slice(&domain.ifft(&seq));

        // expected terms are terms[0] = x, terms[1] = f(x), terms[2] = f(gamma*x)
        let vo_eval_function = |terms: &[VOTerm<F>]| {
            // construct (f(gamma * x) - r * f(x))
            let check_next_term = terms[2].clone() - vo_constant!(common_ratio) * terms[1].clone();

            let mut evaluation_function = check_next_term;
            let mut starting_index = 0;
            // construct each stitch and multiply to the final polynomial
            sequence_lengths.iter().for_each(|sequence_length| {
                let stitch = terms[0].clone()
                    - vo_constant!(domain
                        .element(1)
                        .pow([(starting_index + sequence_length - 1) as u64]));
                starting_index += sequence_length;
                evaluation_function = evaluation_function.clone() * stitch;
            });

            evaluation_function
        };

        let new_geo_vo = NewVO::new(
            vec![0, 1, 1],
            vec![F::one(), F::one(), domain.element(1)],
            &vo_eval_function,
        )
        .unwrap();

        let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        let concrete_oracles = &[x_poly.clone(), f.clone()];
        let new_poly = new_geo_vo.compute_polynomial(concrete_oracles).unwrap();

        let old_vo = GeoSequenceVO::new(&sequence_lengths, domain.element(1), common_ratio);


        let labeled_f = LabeledPolynomial::new(String::from("f"), f, None, None);
        let old_poly = old_vo
            .instantiate_in_coeffs_form(
                &[labeled_f.clone(), labeled_f.clone()],
                &[F::one(), domain.element(1)],
            )
            .unwrap();

        assert_eq!(old_poly, new_poly)

    }
}
