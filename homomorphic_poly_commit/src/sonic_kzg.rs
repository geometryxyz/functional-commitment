use std::collections::BTreeMap;

use ark_ec::PairingEngine;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{
    sonic_pc::SonicKZG10, LCTerm, LabeledCommitment, LinearCombination, PCCommitment, PCRandomness,
};

use crate::{error::Error, AdditivelyHomomorphicPCS};

/// The Default KZG-style commitment scheme
pub type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;

/// A single KZG10 commitment
// pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
//     <E as PairingEngine>::Fr,
//     DensePolynomial<<E as PairingEngine>::Fr>,
// >>::Commitment;

// pub type KZGRandomness<E> = <KZG10<E> as PolynomialCommitment<
//     <E as PairingEngine>::Fr,
//     DensePolynomial<<E as PairingEngine>::Fr>,
// >>::Randomness;

impl<E: PairingEngine> AdditivelyHomomorphicPCS<E::Fr> for SonicKZG10<E, DensePolynomial<E::Fr>> {
    fn aggregate_commitments(
        commitments: &[LabeledCommitment<Self::Commitment>],
        randomness: Option<Vec<Self::Randomness>>,
        lc: &LinearCombination<E::Fr>,
    ) -> Result<(LabeledCommitment<Self::Commitment>, Self::Randomness), Error> {
        let degree_bound = commitments[0].degree_bound();

        let randomness = randomness.map_or(
            vec![Self::Randomness::empty(); commitments.len()],
            |rands| rands,
        );

        // create mapping of label -> commitment and fail if all degree bounds are not the same
        let label_comm_mapping = commitments
            .iter()
            .zip(randomness.iter())
            .map(|(comm, rand)| {
                if comm.degree_bound() != degree_bound {
                    // Can only accumulate commitments that have the same degree bound
                    return Err(Error::MismatchedDegreeBounds(format!(
                        "{} has degree bound {:?}, but {} has degree bound {:?}",
                        commitments[0].label(),
                        degree_bound,
                        comm.label(),
                        comm.degree_bound()
                    )));
                }
                Ok((
                    comm.label().clone(),
                    (comm.commitment().clone(), rand.clone()),
                ))
            })
            .collect::<Result<BTreeMap<_, (_, _)>, Error>>()?;

        // initial values
        let mut aggregate_commitment = Self::Commitment::empty();
        let mut aggregate_randomness = Self::Randomness::empty();

        for (coef, term) in lc.iter() {
            match term {
                // No support for constant terms
                LCTerm::One => return Err(Error::ConstantTermInAggregation),

                // Find the corresponding commitment and randomness in our map; aggregate.
                LCTerm::PolyLabel(label) => match label_comm_mapping.get(label) {
                    Some((comm, rand)) => {
                        aggregate_commitment += (*coef, &comm);
                        aggregate_randomness += (*coef, rand);
                    }
                    None => {
                        return Err(Error::MissingCommitment(format!(
                            "Could not find object with label '{}' when computing '{}'",
                            label,
                            lc.label()
                        )))
                    }
                },
            }
        }

        Ok((
            LabeledCommitment::new(lc.label().clone(), aggregate_commitment, degree_bound),
            aggregate_randomness,
        ))
    }
}

#[cfg(test)]
mod test {
    use crate::{sonic_kzg::KZG10, AdditivelyHomomorphicPCS};
    use ark_bn254::{Bn254, Fr};
    use ark_ff::One;
    use ark_ff::UniformRand;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::UVPolynomial;
    use ark_poly_commit::LinearCombination;
    use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
    use ark_std::rand::thread_rng;
    use rand_core::OsRng;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_aggregate_comm_with_rand() {
        // Parameters
        let rng = &mut thread_rng();
        let maximum_degree: usize = 16;
        let hiding_bound = 1;
        let enforced_degree_bounds = [10];

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
        let a_unlabeled: DensePolynomial<F> = DensePolynomial::rand(7, rng);
        let a_poly = LabeledPolynomial::new(String::from("a"), a_unlabeled, Some(10), Some(1));

        let b_unlabeled: DensePolynomial<F> = DensePolynomial::rand(5, rng);
        let b_poly = LabeledPolynomial::new(String::from("b"), b_unlabeled, Some(10), Some(1));

        let a_plus_2b_poly = a_poly.polynomial().clone() + (b_poly.polynomial() * F::from(2u64));
        let a_plus_2b_poly =
            LabeledPolynomial::new(String::from("a_plus_2b"), a_plus_2b_poly, Some(10), Some(1));
        let polynomials = vec![a_poly.clone(), b_poly.clone()];
        let linear_combination =
            LinearCombination::new("a_plus_2b", vec![(F::one(), "a"), (F::from(2u64), "b")]);

        // Commit Phase
        let (commitments, rands) = PC::commit(&ck, &polynomials, Some(rng)).unwrap();
        let (test_commitment, test_rand) =
            PC::aggregate_commitments(&commitments, Some(rands.to_vec()), &linear_combination)
                .unwrap();

        // Derive evaluation point and generate a query set
        let evaluation_point = Fr::rand(rng);

        // Evaluation Phase, here we only output the evaluation of the linear combination
        let manual_eval = a_plus_2b_poly.evaluate(&evaluation_point);

        // Opening phase
        let opening_challenge = F::rand(rng);
        let lc_opening_proof = PC::open(
            &ck,
            &[a_plus_2b_poly],
            &[test_commitment.clone()],
            &evaluation_point,
            opening_challenge,
            &[test_rand],
            Some(rng),
        )
        .unwrap();

        // Verify
        let res = PC::check(
            &vk,
            &[test_commitment],
            &evaluation_point,
            vec![manual_eval],
            &lc_opening_proof,
            opening_challenge,
            Some(rng),
        )
        .unwrap();

        assert_eq!(true, res)
    }
}
