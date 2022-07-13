use crate::error::Error;
use crate::non_zero_over_k::piop::PIOPforNonZeroOverK;
use ark_ff::{FftField, PrimeField};
use ark_marlin::ahp::prover::ProverMsg;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::LabeledPolynomial;
use ark_std::rand::Rng;

pub struct ProverState<'a, F: PrimeField + FftField> {
    domain_k: &'a GeneralEvaluationDomain<F>,

    f: &'a LabeledPolynomial<F, DensePolynomial<F>>,

    first_oracles: Option<ProverFirstOracles<F>>,
}

/// The first set of prover oracles
#[derive(Clone)]
pub struct ProverFirstOracles<F: PrimeField + FftField> {
    pub g: LabeledPolynomial<F, DensePolynomial<F>>,
}

impl<F: PrimeField + FftField> ProverFirstOracles<F> {
    pub fn iter(&self) -> impl Iterator<Item = &LabeledPolynomial<F, DensePolynomial<F>>> {
        vec![&self.g].into_iter()
    }
}

impl<F: PrimeField + FftField> PIOPforNonZeroOverK<F> {
    pub fn prover_init<'a>(
        domain_k: &'a GeneralEvaluationDomain<F>,
        f: &'a LabeledPolynomial<F, DensePolynomial<F>>,
    ) -> Result<ProverState<'a, F>, Error> {
        Ok(ProverState {
            domain_k,
            f,
            first_oracles: None,
        })
    }

    pub fn prover_first_round<'a, R: Rng>(
        mut state: ProverState<'a, F>,
        _rng: &mut R,
    ) -> Result<(ProverMsg<F>, ProverFirstOracles<F>, ProverState<'a, F>), Error> {
        let f_evals = state.domain_k.fft(state.f.polynomial());

        // Check that all the f_evals are nonzero; otherwise, .inverse() will return None and
        // .unwrap() will panic
        for f in f_evals.iter() {
            if f.inverse().is_none() {
                return Err(Error::FEvalIsZero);
            }
        }

        let g_evals = f_evals
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        let g = DensePolynomial::<F>::from_coefficients_slice(&state.domain_k.ifft(&g_evals));
        let g = LabeledPolynomial::new(
            String::from("g"),
            g.clone(),
            state.f.degree_bound(),
            Some(1),
        );

        // create ProverFirstOracles struct
        let prover_oracles = ProverFirstOracles { g };

        // Prover message
        let msg = ProverMsg::EmptyMessage;

        // Update Prover state
        state.first_oracles = Some(prover_oracles.clone());

        Ok((msg, prover_oracles, state))
    }
}
