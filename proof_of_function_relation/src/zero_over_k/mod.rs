use crate::{
    commitment::HomomorphicPolynomialCommitment,
    error::{to_pc_error, Error},
    label_commitment, label_polynomial,
    transcript::TranscriptProtocol,
    util::{commit_polynomial, powers_of, shift_dense_poly},
    virtual_oracle::VirtualOracle,
};
use ark_ff::{PrimeField, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_poly_commit::{LabeledPolynomial, PCRandomness};
use ark_std::marker::PhantomData;
use merlin::Transcript;
use rand::Rng;

mod proof;

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
        let mut transcript = Transcript::new(b"zero_over_k");

        // commit to the random polynomials
        let labeled_concrete_oracles = concrete_oracles
            .iter()
            .map(|oracle| label_polynomial!(oracle))
            .collect::<Vec<_>>();
        let (concrete_oracles_commitments, _) =
            PC::commit(ck, labeled_concrete_oracles.iter(), None).map_err(to_pc_error::<F, PC>)?;

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

        // let h_prime_commitments =

        // compute the blinded virtual oracle F'[X]
        let f_prime = VO::instantiate(h_primes.as_slice().try_into().unwrap(), alphas);
        let (q1, r) = DenseOrSparsePolynomial::from(&f_prime)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(
                &domain.vanishing_polynomial(),
            ))
            .unwrap();
        assert_eq!(r, DensePolynomial::<F>::zero());

        // commit to q_1
        let (q1_commit, q1_rand) = commit_polynomial::<F, PC>(ck, &q1)?;

        // round 1
        let r_commits = r_commits
            .iter()
            .map(|r_c| r_c.commitment().clone())
            .collect::<Vec<_>>();
        let m_commits = m_commits
            .iter()
            .map(|m_c| m_c.commitment().clone())
            .collect::<Vec<_>>();
        let concrete_oracles_commitments = concrete_oracles_commitments
            .iter()
            .map(|oracle_c| oracle_c.commitment().clone())
            .collect::<Vec<_>>();

        transcript.append(b"oracles", &concrete_oracles_commitments.to_vec());
        transcript.append(b"alphas", &alphas.to_vec());
        transcript.append(b"rs", &r_commits);
        transcript.append(b"ms", &m_commits);
        transcript.append(b"q", &q1_commit);

        let beta_1: F = transcript.challenge_scalar(b"beta_1");
        let beta_2: F = transcript.challenge_scalar(b"beta_2");
        let c: F = transcript.challenge_scalar(b"c");

        //check that beta_1, beta_2, c in F* \ K
        assert_ne!(beta_1, F::zero());
        assert_ne!(domain.evaluate_vanishing_polynomial(beta_1), F::zero());

        assert_ne!(beta_2, F::zero());
        assert_ne!(domain.evaluate_vanishing_polynomial(beta_2), F::zero());

        assert_ne!(c, F::zero());
        assert_ne!(domain.evaluate_vanishing_polynomial(c), F::zero());

        // open q1 at beta_1
        let q1_eval = q1.evaluate(&beta_1);
        let q1_open = PC::open(
            ck,
            &[label_polynomial!(q1)],
            &[label_commitment!(q1_commit)],
            &beta_1,
            F::one(),
            &[q1_rand],
            None,
        );

        // q_2 is defined as r1 + c*r2 + c^2r3 + ...
        let q2 = random_polynomials
            .iter()
            .zip(powers_of(c))
            .fold(DensePolynomial::<F>::zero(), |acc_poly, (r, c_power)| {
                acc_poly + (r.polynomial() * c_power)
            });

        let q2_eval = q2.evaluate(&beta_2);

        //since we commit with randomnes None, linear combination of all randomness will be just empty
        let homomorphic_randomness = PC::Randomness::empty();

        // let (q2_commit, q2_rand) = commit_polynomial::<F, PC>(ck, &q2)?;
        let q2_commitment = PC::multi_scalar_mul(
            &r_commits,
            powers_of(c)
                .take(random_polynomials.len())
                .collect::<Vec<_>>()
                .as_slice(),
        );
        let q2_open = PC::open(
            ck,
            &[label_polynomial!(q2)],
            &[label_commitment!(q2_commitment)],
            &beta_2,
            F::one(),
            &[homomorphic_randomness.clone()],
            None,
        );

        let one = F::one();
        let h_openings = h_primes
            .iter()
            .zip(alphas.iter())
            .zip(concrete_oracles_commitments.iter())
            .zip(m_commits.iter())
            .map(|(((h_prime, &alpha_i), oracle_commitment), m_commitment)| {
                let h_prime_commitment = PC::multi_scalar_mul(
                    &[oracle_commitment.clone(), m_commitment.clone()],
                    &[one, one],
                );
                PC::open(
                    ck,
                    &[label_polynomial!(h_prime)],
                    &[label_commitment!(h_prime_commitment)],
                    &(alpha_i * beta_1),
                    F::one(),
                    &[homomorphic_randomness.clone()],
                    None,
                )
            })
            .collect::<Vec<_>>();

        let m_openings = masking_polynomials
            .iter()
            .zip(m_commits.iter())
            .zip(alphas.iter())
            .map(|((m_poly, m_commit), &alpha_i)| {
                PC::open(
                    ck,
                    [m_poly],
                    &[label_commitment!(m_commit)],
                    &(alpha_i * beta_2),
                    F::one(),
                    &[homomorphic_randomness.clone()],
                    None,
                )
            })
            .collect::<Vec<_>>();

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
