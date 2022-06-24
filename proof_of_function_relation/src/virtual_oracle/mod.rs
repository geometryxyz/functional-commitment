use crate::error::Error;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
};
use ark_poly_commit::LabeledPolynomial;

pub mod add_vo;
pub mod eq_vo;
pub mod geometric_sequence_vo;
pub mod inverse_check_oracle;
pub mod prod_vo;
pub mod product_check_oracle;
pub mod square_check_oracle;

pub trait VirtualOracle<F: PrimeField> {
    /// returns instantiation of virtual oracle from concrete oracles and shifting factors
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error>;

    /// returns instantiation of virtual oracle in evaluatio form
    /// everything will be performed in coset since later division with vanishin poly will not be possible unless in coset
    fn instantiate_in_evals_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
        domain: &GeneralEvaluationDomain<F>,
    ) -> Result<Vec<F>, Error>;

    /// verifier always have the access to the oracle X^n through point
    fn query(&self, evals: &[F], point: F) -> Result<F, Error>;

    /// each new (f, alpha) pair should be mapped to new h (new concrete oracle)
    /// this function provides mapping from concrete oracle inidices to h indices
    fn mapping_vector(&self) -> Vec<usize>;

    fn num_of_oracles(&self) -> usize;

    /// returns degree of instantiation polynomial given the domain size
    /// for domain_size = n,
    /// in order to stay compatible with zero over k argument where each polynomial is being masked with deg(n + 1) (r * zh)
    /// each concrete oracle will be of degree (n+1)
    fn degree_bound(&self, domain_size: usize) -> usize;

    /// in general degree of vo can be bigger then domain_size. ex: deg(f(x) * g(x)) = 2 * (n + 1) = 2n + 2
    /// scaling_factor(k) is some degree of 2 such that deg(vo) - deg(vanishing_poly) can be computed using ffts
    /// for example in plonk it is 4n
    fn compute_scaling_factor(&self, domain: &GeneralEvaluationDomain<F>) -> usize {
        let n = domain.size();
        let kn = self.degree_bound(n) - n;

        // this migh seem unintuitive, but since our concrete oracles are all deg (n+1) to represent just one in evals form we must work in double domain size
        if kn <= n {
            return 2;
        }

        let kn = if kn.is_power_of_two() {
            kn
        } else {
            kn.checked_next_power_of_two().unwrap()
        };

        kn / n
    }

    fn name(&self) -> String;
}

pub trait EvaluationsProvider<F: PrimeField> {
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error>;
}

impl<F: PrimeField> EvaluationsProvider<F> for Vec<LabeledPolynomial<F, DensePolynomial<F>>> {
    /// Instantiate and evaluate the virtual oracle. Returns an error if the length of
    /// concrete_oracles differs from the number of concrete oracles in the Description.
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        alpha_coeffs: &Vec<F>,
    ) -> Result<F, Error> {
        if self.len() != virtual_oracle.num_of_oracles() {
            return Err(Error::EvaluationError);
        }

        let poly = virtual_oracle
            .instantiate_in_coeffs_form(&self, alpha_coeffs)
            .unwrap();
        return Ok(poly.evaluate(&point));
    }
}

//1. add query funct to oracle vec<f>.evaluate => vo.query(vec<f>)
//2. return function from oracle and then return func(self) of vec<f>
impl<F: PrimeField> EvaluationsProvider<F> for Vec<F> {
    /// Return the evaluation of the virtual oracle given a list of evaluations of each concrete
    /// oracle. The length of the input vector of evaluations must equal to the number of concrete
    /// oracles in the Description.
    fn evaluate<VO: VirtualOracle<F>>(
        &self,
        virtual_oracle: &VO,
        point: F,
        _: &Vec<F>,
    ) -> Result<F, Error> {
        virtual_oracle.query(&self, point)
    }
}
