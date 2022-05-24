use crate::{
    commitment::HomomorphicPolynomialCommitment,
    error::{to_pc_error, Error},
    label_polynomial, 
    util::{shift_dense_poly, commit_polynomial},
    virtual_oracle::VirtualOracle,
};
use ark_ff::{PrimeField, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::LabeledPolynomial;
use ark_std::marker::PhantomData;
use rand::Rng;

struct ZeroOverK<F, PC>
where
    F: PrimeField,
    PC: HomomorphicPolynomialCommitment<F>,
{
    _field: PhantomData<F>,
    _commitment_scheme: PhantomData<PC>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>> ZeroOverK<F, PC> {
    pub fn prove<
        const NUM_OF_CONCRETE_ORACLES: usize,
        VO: VirtualOracle<F, NUM_OF_CONCRETE_ORACLES>,
        R: Rng,
    >(
        concrete_oracles: &[DensePolynomial<F>; NUM_OF_CONCRETE_ORACLES],
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
        domain: &GeneralEvaluationDomain<F>,
        ck: &PC::CommitterKey,
        rng: &mut R,
    ) -> Result<(), Error> {
        // generate the masking polynomials and keep a record of the random polynomials that were used
        let (random_polynomials, masking_polynomials) = Self::compute_maskings(domain, alphas, rng);

        // commit to the random polynomials
        let (r_commits, _) =
            PC::commit(ck, random_polynomials.iter(), None).map_err(to_pc_error::<F, PC>)?;

        // commit to the masking polynomials
        let (m_commits, _) =
            PC::commit(ck, masking_polynomials.iter(), None).map_err(to_pc_error::<F, PC>)?;

        // compute the masked oracles
        let h_primes = concrete_oracles
            .iter()
            .zip(masking_polynomials.iter())
            .map(|(oracle, masking_poly)| oracle + masking_poly.polynomial())
            .collect::<Vec<_>>();

        // compute the blinded virtual oracle F'[X]
        let f_prime = VO::instantiate(h_primes.as_slice().try_into().unwrap(), alphas);
        let (q1, r) = DenseOrSparsePolynomial::from(&f_prime)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(
                &domain.vanishing_polynomial(),
            ))
            .unwrap();
        assert_eq!(r, DensePolynomial::<F>::zero());

        // commit to q_1
        let (q1_commit, _) = commit_polynomial::<F, PC>(ck, &q1)?;

        // let mut oracle_div_z_h = oracle_evals.as_slice().clone();
        // domain.divide_by_vanishing_poly_on_coset_in_place(&oracle_div_z_h);
        //we divide with z_h, we want to divide in evaluation form,

        Ok(())
    }

    /// computes array of m_i = ri(alpha_i^-1) * zk(alpha_i^-1)
    fn compute_maskings<R: Rng, const NUM_OF_CONCRETE_ORACLES: usize>(
        domain: &GeneralEvaluationDomain<F>,
        alphas: &[F; NUM_OF_CONCRETE_ORACLES],
        rng: &mut R,
    ) -> (
        Vec<LabeledPolynomial<F, DensePolynomial<F>>>,
        Vec<LabeledPolynomial<F, DensePolynomial<F>>>,
    ) {
        //r is defined as polynomial degree < 2
        let degree = 1;

        let mut random_polynomials = Vec::with_capacity(NUM_OF_CONCRETE_ORACLES);
        let mut masking_polynomials = Vec::with_capacity(NUM_OF_CONCRETE_ORACLES);

        for i in 0..NUM_OF_CONCRETE_ORACLES {
            let r = DensePolynomial::<F>::rand(degree, rng);
            let shifting_factor = alphas[i].inverse().unwrap();
            let r_shifted = shift_dense_poly(&r, &shifting_factor);
            let vanishing_shifted =
                shift_dense_poly(&domain.vanishing_polynomial().into(), &shifting_factor);

            random_polynomials[i] = label_polynomial!(r_shifted.clone());
            masking_polynomials[i] = label_polynomial!(&r_shifted * &vanishing_shifted);
        }

        (random_polynomials, masking_polynomials)
    }

    pub fn verify() -> () {
        // let a_coeffs = vec![F::from(5), F::from(0), F::from(0), F::from(0)];
        // let a_poly = DensePolynomial::from_coefficients_slice(&a_coeffs);
        // let p: P = q.into();
    }
}
