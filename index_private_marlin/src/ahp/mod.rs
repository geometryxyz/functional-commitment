use std::{borrow::Borrow, marker::PhantomData};

use ac_compiler::R1CSfIndex;
use ark_ff::{Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::{LCTerm, LinearCombination};
use ark_std::cfg_iter_mut;

use self::constraint_systems::LabeledPolynomial;

pub(crate) mod constraint_systems;

pub mod indexer;
pub mod prover;
pub mod verifier;

pub struct AHPForR1CS<F: Field> {
    field: PhantomData<F>,
}

impl<F: PrimeField> AHPForR1CS<F> {
    /// The labels for the polynomials output by the AHP index private prover.
    #[rustfmt::skip]
    pub const PROVER_POLYNOMIALS: [&'static str; 10] = [
        // First sumcheck
        "z", "z_a", "z_b", "mask_poly", "inner_mask_poly", "t", "g_1", "h_1",
        // Second sumcheck
        "f", "g_2", //TODO: REMOVE h2, it's calculated inside zero over K
    ];

    /// THe linear combinations that are statically known to evaluate to zero in index private version.
    pub const LC_WITH_ZERO_EVAL: [&'static str; 2] = ["outer_sumcheck", "f_sumcheck"];

    pub const MATRIX_POLY_LABELS: [&'static str; 9] = [
        "a_row", "a_col", "a_val", "b_row", "b_col", "b_val", "c_row", "c_col", "c_val",
    ];

    pub const WELL_FORMATION_LABELS: [&'static str; 4] =
        ["pi_lde", "vh_gt_x", "output_lde", "vh_lt_y"];

    pub(crate) fn polynomial_labels() -> impl Iterator<Item = String> {
        Self::PROVER_POLYNOMIALS.iter().map(|s| s.to_string())
    }

    pub(crate) fn matrix_poly_labels() -> impl Iterator<Item = String> {
        Self::MATRIX_POLY_LABELS.iter().map(|s| s.to_string())
    }

    pub(crate) fn well_formation_labels() -> impl Iterator<Item = String> {
        Self::WELL_FORMATION_LABELS.iter().map(|s| s.to_string())
    }

    /// The maximum degree of polynomials produced by the indexer and prover
    /// of this protocol.
    pub fn max_degree(info: &R1CSfIndex) -> Result<usize, Error> {
        let zk_bound = 1;
        let domain_h_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_constraints)
                .ok_or(Error::DomainTooLarge)?;
        let domain_k_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_non_zero_entries)
                .ok_or(Error::DomainTooLarge)?;

        /*
            From Aurora paper we know that deg(mask_poly) should be same degree of oracle we are proving sumcheck for: TODO -> add link

            For outer sumcheck we have following bounds:
                deg(t) = |H|
                deg(z) = |H| + zk_bound
                deg(r(alpha, X)) = |H| - 1, => deg of vanishing is |H| and it's derivative is |H| * X^(|H| - 1)
                deg(z_a) = |H| + zk_bound
                deg(z_b) = |H| + zk_bound
                !!! zc = za * zb !!! deg(z_c) = 2|H| + 2zk_bound

                deg(vanishing) = |H|
                deg(g1) = |H| - 2
                deg(h1) = largest degree - deg(vanishing)

            so the largest degree in outer sumcheck is r(alpha, X) * eta_z * z_c which gives:
            deg(mask_poly) = |H| - 1 + 2|H| + 2*zk_bound = 3|H| + 2*zk_bound - 1 !!!TODO: check why -2 in ark_marlin diagram.pdf

            For inner sumcheck part we have following bounds:
                deg(val) = |K| - 1
                deg(row) = |K| - 1
                deg(col) = |K| - 1
                deg(row_col) = 2|K| - 2
                deg(a) = deg(val) + 2 * deg(row_col) = |K| - 1 + 2 * (2|K| - 2) = |K| - 1 + 4|K| - 4 = 5|K| - 5
                deg(b) = 3 * deg(row_col) = 6|K| - 6
                deg(f) = |K| - 1
                deg(mask_poly) = |K| - 1
                deg(g2) = |K| - 2
                deg(h2) = max(a, f*b) - vanishing = 7|K| - 7 - |K| = 6|K| - 7

                !!!
                    when working with zero over K, our oracles are masked with polynomials of degree |K| + 1
                    when proving zero over k for rational sumcheck we have to commit to polynomial h2

                    that's why degrees of concrete oracles will be |K| + 1 instead of |K| - 1
                    deg(val) = |K| + 1
                    deg(row) = |K| + 1
                    deg(col) = |K| + 1
                    deg(row_col) = 2|K| + 2
                    deg(a) = deg(val) + 2 * deg(row_col) = |K| + 1 + 2 * (2|K| + 2) = |K| + 1 + 4|K| + 4 = 5|K| + 5
                    deg(b) = 3 * deg(row_col) = 6|K| + 6
                    deg(h2) = max(a, f*b) - vanishing = 7|K| + 7 - |K| = 6|K| + 7
                !!!
        */

        Ok(
            *[3 * domain_h_size + 2 * zk_bound - 1, 6 * domain_k_size + 7]
                .iter()
                .max()
                .unwrap(),
        )
    }

    /// Get all the strict degree bounds enforced in the AHP.
    pub fn get_degree_bounds(info: &R1CSfIndex) -> [usize; 4] {
        let mut degree_bounds = [0usize; 4];

        let h_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_constraints)
                .unwrap();
        let k_size =
            GeneralEvaluationDomain::<F>::compute_size_of_domain(info.number_of_non_zero_entries)
                .unwrap();

        degree_bounds[0] = h_size - 2; // for g1
        degree_bounds[1] = k_size - 2; // for g2
        degree_bounds[2] = 2; // mask poly in zero over k
        degree_bounds[3] = k_size + 1; // this is needed in tft proof TODO: create better index function abstraction
        degree_bounds
    }

    /// Construct the linear combinations that are checked by the AHP.
    #[allow(non_snake_case)]
    pub fn construct_linear_combinations<E>(
        evals: &E,
        state: &verifier::VerifierState<F>,
    ) -> Result<Vec<LinearCombination<F>>, Error>
    where
        E: EvaluationsProvider<F>,
    {
        let domain_h = state.domain_h;
        let domain_k = state.domain_k;
        let k_size = domain_k.size_as_field_element();

        let first_round_msg = state.first_round_msg.unwrap();
        let alpha = first_round_msg.alpha;
        let eta_a = first_round_msg.eta_a;
        let eta_b = first_round_msg.eta_b;
        let eta_c = first_round_msg.eta_c;

        let beta = state.second_round_msg.unwrap().beta;
        let gamma = state.gamma.unwrap();

        let mut linear_combinations = Vec::new();

        // Constants
        let r_alpha_at_beta = domain_h.eval_unnormalized_bivariate_lagrange_poly(alpha, beta);
        let _v_H_at_alpha = domain_h.evaluate_vanishing_polynomial(alpha);
        let v_H_at_beta = domain_h.evaluate_vanishing_polynomial(beta);

        // Outer sumcheck:
        let z_b = LinearCombination::new("z_b", vec![(F::one(), "z_b")]);
        let g_1 = LinearCombination::new("g_1", vec![(F::one(), "g_1")]);
        let t = LinearCombination::new("t", vec![(F::one(), "t")]);

        let z_b_at_beta = evals.get_lc_eval(&z_b, beta)?;
        let t_at_beta = evals.get_lc_eval(&t, beta)?;
        let g_1_at_beta = evals.get_lc_eval(&g_1, beta)?;

        #[rustfmt::skip]
        let outer_sumcheck = LinearCombination::new(
            "outer_sumcheck",
            vec![
                (F::one(), "mask_poly".into()),

                (r_alpha_at_beta * (eta_a + eta_c * z_b_at_beta), "z_a".into()),
                (r_alpha_at_beta * eta_b * z_b_at_beta, LCTerm::One),

                (-t_at_beta, "z".into()),

                (-v_H_at_beta, "h_1".into()),
                (-beta * g_1_at_beta, LCTerm::One),
            ],
        );
        debug_assert!(evals.get_lc_eval(&outer_sumcheck, beta)?.is_zero());

        linear_combinations.push(z_b);
        linear_combinations.push(g_1);
        linear_combinations.push(t);
        linear_combinations.push(outer_sumcheck);

        // f sumcheck test

        let g_2 = LinearCombination::new("g_2", vec![(F::one(), "g_2")]);
        let g_2_at_gamma = evals.get_lc_eval(&g_2, gamma)?;

        let f_sumcheck = LinearCombination::<F>::new(
            "f_sumcheck",
            vec![
                (F::one(), "inner_mask_poly".into()),
                (F::one(), "f".into()),
                (-gamma * g_2_at_gamma, LCTerm::One),
                (-(t_at_beta / k_size), LCTerm::One),
            ],
        );

        debug_assert!(evals.get_lc_eval(&f_sumcheck, gamma)?.is_zero());
        linear_combinations.push(g_2);
        linear_combinations.push(f_sumcheck);

        linear_combinations.sort_by(|a, b| a.label.cmp(&b.label));
        Ok(linear_combinations)
    }
}

/// Describes the failure modes of the AHP scheme.
#[derive(Debug)]
pub enum Error {
    /// During verification, a required evaluation is missing
    MissingEval(String),
    /// The instance generated during proving does not match that in the index.
    InstanceDoesNotMatchIndex,
    DomainTooLarge,
    /// Number of constraints is larger than number of non zero elements (we don't allow this because Discrete-log Comparison in proof of function fails)
    DomainHLargerThanDomainK,
}

/// Abstraction that provides evaluations of (linear combinations of) polynomials
///
/// Intended to provide a common interface for both the prover and the verifier
/// when constructing linear combinations via `AHPForR1CS::construct_linear_combinations`.
pub trait EvaluationsProvider<F>
where
    F: Field,
{
    /// Get the evaluation of linear combination `lc` at `point`.
    fn get_lc_eval(&self, lc: &LinearCombination<F>, point: F) -> Result<F, Error>;
}

impl<'a, F: Field> EvaluationsProvider<F> for ark_poly_commit::Evaluations<F, F> {
    fn get_lc_eval(&self, lc: &LinearCombination<F>, point: F) -> Result<F, Error> {
        let key = (lc.label.clone(), point);
        self.get(&key)
            .map(|v| *v)
            .ok_or(Error::MissingEval(lc.label.clone()))
    }
}

impl<F: Field, T: Borrow<LabeledPolynomial<F>>> EvaluationsProvider<F> for Vec<T> {
    fn get_lc_eval(&self, lc: &LinearCombination<F>, point: F) -> Result<F, Error> {
        let mut eval = F::zero();
        for (coeff, term) in lc.iter() {
            let value = if let LCTerm::PolyLabel(label) = term {
                self.iter()
                    .find(|p| {
                        let p: &LabeledPolynomial<F> = (*p).borrow();
                        p.label() == label
                    })
                    .ok_or(Error::MissingEval(format!(
                        "Missing {} for {}",
                        label, lc.label
                    )))?
                    .borrow()
                    .evaluate(&point)
            } else {
                assert!(term.is_one());
                F::one()
            };
            eval += *coeff * value
        }
        Ok(eval)
    }
}

/// The derivative of the vanishing polynomial
pub trait UnnormalizedBivariateLagrangePoly<F: ark_ff::FftField> {
    /// Evaluate the polynomial
    fn eval_unnormalized_bivariate_lagrange_poly(&self, x: F, y: F) -> F;

    /// Evaluate over a batch of inputs
    fn batch_eval_unnormalized_bivariate_lagrange_poly_with_diff_inputs(&self, x: F) -> Vec<F>;

    /// Evaluate the magic polynomial over `self`
    fn batch_eval_unnormalized_bivariate_lagrange_poly_with_same_inputs(&self) -> Vec<F>;
}

impl<F: PrimeField> UnnormalizedBivariateLagrangePoly<F> for GeneralEvaluationDomain<F> {
    fn eval_unnormalized_bivariate_lagrange_poly(&self, x: F, y: F) -> F {
        if x != y {
            (self.evaluate_vanishing_polynomial(x) - self.evaluate_vanishing_polynomial(y))
                / (x - y)
        } else {
            self.size_as_field_element() * x.pow(&[(self.size() - 1) as u64])
        }
    }

    fn batch_eval_unnormalized_bivariate_lagrange_poly_with_diff_inputs(&self, x: F) -> Vec<F> {
        let vanish_x = self.evaluate_vanishing_polynomial(x);
        let mut inverses: Vec<F> = self.elements().map(|y| x - y).collect();
        ark_ff::batch_inversion(&mut inverses);

        cfg_iter_mut!(inverses).for_each(|denominator| *denominator *= vanish_x);
        inverses
    }

    fn batch_eval_unnormalized_bivariate_lagrange_poly_with_same_inputs(&self) -> Vec<F> {
        let mut elems: Vec<F> = self
            .elements()
            .map(|e| e * self.size_as_field_element())
            .collect();
        elems[1..].reverse();
        elems
    }
}
