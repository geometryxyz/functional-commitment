//! Useful commitment stuff
use ark_ec::{msm::VariableBaseMSM, PairingEngine};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{
    sonic_pc::SonicKZG10, LabeledCommitment, PCRandomness, PolynomialCommitment,
};

/// A homomorphic polynomial commitment
pub trait HomomorphicPolynomialCommitment<F>: PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
    Self::VerifierKey: core::fmt::Debug,
{
    // type HomomorphicRandomness: Add;

    /// Combine a linear combination of homomorphic commitments
    fn multi_scalar_mul(
        commitments: &[LabeledCommitment<Self::Commitment>],
        scalars: &[F],
    ) -> Self::Commitment;
    fn aggregate_randomness(rands: &[Self::Randomness]) -> Self::Randomness;
}

/// The Default KZG-style commitment scheme
pub type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;

/// A single KZG10 commitment
pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Commitment;

pub type KZGRandomness<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Randomness;

impl<E> HomomorphicPolynomialCommitment<E::Fr> for KZG10<E>
where
    E: PairingEngine,
{
    fn multi_scalar_mul(
        commitments: &[LabeledCommitment<KZG10Commitment<E>>],
        scalars: &[E::Fr],
    ) -> KZG10Commitment<E> {
        let scalars_repr = scalars
            .iter()
            .map(<E::Fr as PrimeField>::into_repr)
            .collect::<Vec<_>>();

        let points_repr = commitments
            .iter()
            .map(|c| c.commitment().0)
            .collect::<Vec<_>>();

        ark_poly_commit::kzg10::Commitment::<E>(
            VariableBaseMSM::multi_scalar_mul(&points_repr, &scalars_repr).into(),
        )
    }

    fn aggregate_randomness(rands: &[KZGRandomness<E>]) -> KZGRandomness<E> {
        if rands.len() == 0 {
            return KZGRandomness::<E>::empty();
        }
        let mut acc = rands[0].clone();
        for rand in rands.iter().skip(1) {
            acc += rand;
        }

        acc
    }
}

#[cfg(test)]
mod test {
    use crate::commitment::{HomomorphicPolynomialCommitment, KZG10};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::One;
    use ark_ff::UniformRand;
    use ark_marlin::ahp::EvaluationsProvider;
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
    };
    use ark_poly_commit::LinearCombination;
    use ark_poly_commit::QuerySet;
    use ark_poly_commit::{
        Evaluations as CommitEvaluations, LabeledCommitment, LabeledPolynomial,
        PolynomialCommitment,
    };
    use ark_std::rand::thread_rng;
    use rand_core::OsRng;
    use std::iter;

    type F = Fr;
    type PC = KZG10<Bn254>;

    fn generate_test_polynomials_with_degree_bounds() -> (
        GeneralEvaluationDomain<F>,
        LabeledPolynomial<F, DensePolynomial<F>>,
        LabeledPolynomial<F, DensePolynomial<F>>,
        LabeledPolynomial<F, DensePolynomial<F>>,
    ) {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        //test oracle to be zero at roots of unity
        let a_evals = vec![F::from(2u64), F::from(4u64)];
        let a_poly = Evaluations::from_vec_and_domain(a_evals, domain).interpolate();
        let a_poly = LabeledPolynomial::new(String::from("a"), a_poly, Some(4), Some(1));

        let b_evals = vec![
            F::from(5u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let b_poly = Evaluations::from_vec_and_domain(b_evals, domain).interpolate();
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, Some(4), Some(1));

        let a_plus_b_poly = a_poly.polynomial() + b_poly.polynomial();
        let a_plus_b_poly =
            LabeledPolynomial::new(String::from("a_plus_b"), a_plus_b_poly, Some(4), Some(1));

        (domain, a_poly, b_poly, a_plus_b_poly)
    }

    fn generate_test_polynomials_no_degree_bounds() -> (
        GeneralEvaluationDomain<F>,
        LabeledPolynomial<F, DensePolynomial<F>>,
        LabeledPolynomial<F, DensePolynomial<F>>,
        LabeledPolynomial<F, DensePolynomial<F>>,
    ) {
        let n = 4;
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();

        //test oracle to be zero at roots of unity
        let a_evals = vec![F::from(2u64), F::from(4u64)];
        let a_poly = Evaluations::from_vec_and_domain(a_evals, domain).interpolate();
        let a_poly = LabeledPolynomial::new(String::from("a"), a_poly, None, Some(1));

        let b_evals = vec![
            F::from(5u64),
            F::from(4 as u64),
            F::from(6 as u64),
            F::from(8 as u64),
        ];
        let b_poly = Evaluations::from_vec_and_domain(b_evals, domain).interpolate();
        let b_poly = LabeledPolynomial::new(String::from("b"), b_poly, None, Some(1));

        let a_plus_b_poly = a_poly.polynomial() + b_poly.polynomial();
        let a_plus_b_poly =
            LabeledPolynomial::new(String::from("a_plus_b"), a_plus_b_poly, None, Some(1));

        (domain, a_poly, b_poly, a_plus_b_poly)
    }

    #[test]
    fn test_commit_with_bounds() {
        let rng = &mut thread_rng();
        let maximum_degree: usize = 16;
        let hiding_bound = 1;
        let enforced_degree_bounds = [4];

        let (_domain, a_poly, b_poly, a_plus_b_poly) =
            generate_test_polynomials_with_degree_bounds();

        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            maximum_degree,
            hiding_bound,
            Some(&enforced_degree_bounds),
        )
        .unwrap();

        let (commitments, rands) = PC::commit(&ck, &[a_poly, b_poly], Some(rng)).unwrap();
        let coefficients = vec![Fr::one(); commitments.len()];

        let a_plus_b_commit = PC::multi_scalar_mul(&commitments, &coefficients);
        let a_plus_b_commit =
            LabeledCommitment::new(String::from("a_plus_b"), a_plus_b_commit, Some(4));

        let challenge = Fr::rand(rng);

        let evaluation_point = Fr::rand(rng);
        let eval = a_plus_b_poly.evaluate(&evaluation_point);

        let homomorphic_randomness = PC::aggregate_randomness(&rands);

        let opening_proof = PC::open(
            &ck,
            &[a_plus_b_poly.clone()],
            &[a_plus_b_commit.clone()],
            &evaluation_point,
            challenge,
            &[homomorphic_randomness],
            None,
        )
        .unwrap();

        let res = PC::check(
            &vk,
            &[a_plus_b_commit],
            &evaluation_point,
            iter::once(eval),
            &opening_proof,
            challenge,
            None,
        );

        println!("{:?}", res);
    }

    #[test]
    // Note that we cannot use LinearCombination with polynomials that have degree bounds
    fn test_linear_combinations() {
        let rng = &mut thread_rng();

        // Parameters
        let maximum_degree: usize = 16;
        let hiding_bound = 1;
        let enforced_degree_bounds = [4];

        // Setup the commitment scheme
        let pp = PC::setup(maximum_degree, None, &mut OsRng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            maximum_degree,
            hiding_bound,
            Some(&enforced_degree_bounds),
        )
        .unwrap();

        // Define polynomials and a linear combination
        let (_domain, a_poly, b_poly, a_plus_b_poly) = generate_test_polynomials_no_degree_bounds();
        let polynomials = vec![a_poly, b_poly];
        let linear_combination =
            LinearCombination::new("a_plus_b", vec![(F::one(), "a"), (F::one(), "b")]);
        let linear_combinations = vec![linear_combination];

        // Commit Phase
        let (commitments, rands) = PC::commit(&ck, &polynomials, Some(rng)).unwrap();

        // Derive evaluation point and generate a query set
        let evaluation_point = Fr::rand(rng);
        let mut query_set = QuerySet::new();
        query_set.insert((
            linear_combinations[0].label().clone(),
            (String::from("challenge"), evaluation_point),
        ));

        for (poly_label, (point_label, point)) in &query_set {
            println!(
                "Query poly {} at point {}: {}",
                poly_label, point_label, point
            );
        }

        // Evaluation Phase, here we only output the evaluation of the linear combination
        let manual_eval = a_plus_b_poly.evaluate(&evaluation_point);
        let lc_eval = polynomials
            .get_lc_eval(&linear_combinations[0], evaluation_point)
            .unwrap();
        assert_eq!(manual_eval, lc_eval);

        let mut evaluations = CommitEvaluations::new();
        evaluations.insert(
            (linear_combinations[0].label().clone(), evaluation_point),
            lc_eval,
        );
        for ((poly_label, point), eval) in &evaluations {
            println!(
                "Evaluated poly {} at point {}. Value is {}",
                poly_label, point, eval
            );
        }

        //Generate opening proof
        let opening_challenge = Fr::rand(rng);
        let opening_proof = PC::open_combinations(
            &ck,
            &linear_combinations,
            &polynomials,
            &commitments,
            &query_set,
            opening_challenge,
            &rands,
            Some(rng),
        )
        .unwrap();

        // Verify opening proof
        assert_eq!(
            true,
            PC::check_combinations(
                &vk,
                &linear_combinations,
                &commitments,
                &query_set,
                &evaluations,
                &opening_proof,
                opening_challenge,
                rng
            )
            .is_ok()
        )
    }
}
