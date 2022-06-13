#![allow(dead_code)]

use crate::discrete_log_comparison::piop::PIOPforDLComparison;
use crate::error::Error;
use crate::label_polynomial;
use crate::util::*;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::{One, PrimeField, SquareRootField, Zero};
use ark_marlin::ahp::prover::ProverMsg;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::LabeledPolynomial;
use ark_std::rand::Rng;
use std::iter;

pub struct ProverState<'a, F: PrimeField + SquareRootField> {
    domain_k: &'a GeneralEvaluationDomain<F>,

    domain_h: &'a GeneralEvaluationDomain<F>,

    f: &'a LabeledPolynomial<F, DensePolynomial<F>>,

    g: &'a LabeledPolynomial<F, DensePolynomial<F>>,

    first_oracles: Option<ProverFirstOracles<F>>,

    pub a_s: Option<Vec<F>>,

    pub c_s: Option<Vec<usize>>,

    pub delta: Option<F>,
}

/// The first set of prover oracles
#[derive(Clone)]
pub struct ProverFirstOracles<F: PrimeField + SquareRootField> {
    pub s: LabeledPolynomial<F, DensePolynomial<F>>,

    pub f_prime: LabeledPolynomial<F, DensePolynomial<F>>,

    pub g_prime: LabeledPolynomial<F, DensePolynomial<F>>,

    pub s_prime: LabeledPolynomial<F, DensePolynomial<F>>,

    pub h: LabeledPolynomial<F, DensePolynomial<F>>,
}

impl<F: PrimeField + SquareRootField> ProverFirstOracles<F> {
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F, DensePolynomial<F>>> {
        vec![
            &self.s,
            &self.f_prime,
            &self.g_prime,
            &self.s_prime,
            &self.h,
        ]
        .into_iter()
    }
}

#[allow(dead_code)]
impl<F: PrimeField + SquareRootField> PIOPforDLComparison<F> {
    pub fn prover_init<'a>(
        domain_k: &'a GeneralEvaluationDomain<F>,
        domain_h: &'a GeneralEvaluationDomain<F>,
        f: &'a LabeledPolynomial<F, DensePolynomial<F>>,
        g: &'a LabeledPolynomial<F, DensePolynomial<F>>,
    ) -> Result<ProverState<'a, F>, Error> {
        Ok(ProverState {
            domain_k,
            domain_h,
            f,
            g,
            first_oracles: None,
            a_s: None,
            c_s: None,
            delta: None,
        })
    }

    pub fn prover_first_round<'a, R: Rng>(
        mut state: ProverState<'a, F>,
        rng: &mut R,
    ) -> Result<(ProverMsg<F>, ProverFirstOracles<F>, ProverState<'a, F>), Error> {
        let m = state.domain_k.size();
        let n = state.domain_h.size();
        let delta = state.domain_h.element(1).sqrt().unwrap();

        // Compute s
        let f_evals = state.domain_k.fft(state.f.polynomial().coeffs());
        let g_evals = state.domain_k.fft(state.g.polynomial().coeffs());

        let s_evals: Vec<F> = f_evals
            .iter()
            .zip(g_evals.iter())
            .map(|(&f_eval, g_eval)| f_eval * g_eval.inverse().unwrap())
            .collect();

        let s = DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&s_evals));
        let s = label_polynomial!(s);

        // For b in {f, g, s}, compute b_prime
        let omegas = state.domain_h.elements();
        let omega_powers_mapping = omegas
            .enumerate()
            .map(|(power, omega)| (omega, power))
            .collect::<std::collections::HashMap<_, _>>();

        let f_prime_evals = f_evals
            .iter()
            .map(|f| {
                let power = omega_powers_mapping.get(f).expect("F evals are wrong");
                delta.pow(&[(*power) as u64])
            })
            .collect::<Vec<_>>();

        let f_prime =
            DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&f_prime_evals));
        let f_prime = label_polynomial!(f_prime);

        let g_prime_evals = g_evals
            .iter()
            .map(|g| {
                let power = omega_powers_mapping.get(g).expect("G evals are wrong");
                delta.pow(&[(*power) as u64])
            })
            .collect::<Vec<_>>();

        let g_prime =
            DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&g_prime_evals));
        let g_prime = label_polynomial!(g_prime);

        let s_prime_evals = s_evals
            .iter()
            .map(|s| {
                let power = omega_powers_mapping.get(s).expect("S evals are wrong");
                delta.pow(&[(*power) as u64])
            })
            .collect::<Vec<_>>();

        let s_prime =
            DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&s_prime_evals));
        let s_prime = label_polynomial!(s_prime);

        // Compute the sequence h
        let mut a_s = vec![F::one()];
        let mut c_s = vec![n];

        let to_pad = m - n;
        if to_pad > 0 {
            a_s.push(F::zero());
            c_s.push(to_pad);
        }

        let seq = generate_sequence(delta, &a_s, &c_s);
        let h = DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&seq));
        let h = label_polynomial!(h);

        // create ProverFirstOracles struct
        let prover_oracles = ProverFirstOracles {
            s,
            f_prime,
            g_prime,
            s_prime,
            h,
        };

        // Prover message
        let msg = ProverMsg::EmptyMessage;

        // Update Prover state
        state.first_oracles = Some(prover_oracles.clone());
        state.a_s = Some(a_s);
        state.c_s = Some(c_s);
        state.delta = Some(delta);

        Ok((msg, prover_oracles, state))
    }
}
