use std::collections::BTreeMap;

use ark_ec::PairingEngine;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{
    kzg10, marlin_pc::MarlinKZG10, LCTerm, LabeledCommitment, LinearCombination, PCCommitment,
    PCRandomness, PolynomialCommitment,
};

use crate::{error::Error, AdditivelyHomomorphicPCS};

/// The Default KZG-style commitment scheme
pub type KZG10<E> = MarlinKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;

/// A single KZG10 commitment
// pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
//     <E as PairingEngine>::Fr,
//     DensePolynomial<<E as PairingEngine>::Fr>,
// >>::Commitment;

pub type KZGRandomness<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Randomness;

impl<E: PairingEngine> AdditivelyHomomorphicPCS<E::Fr> for MarlinKZG10<E, DensePolynomial<E::Fr>> {
    fn aggregate_commitments(
        commitments: &[LabeledCommitment<Self::Commitment>],
        randomness: Option<Vec<Self::Randomness>>,
        lc: &LinearCombination<E::Fr>,
    ) -> Result<(LabeledCommitment<Self::Commitment>, Self::Randomness), Error> {
        let randomness = randomness.map_or(
            vec![
                Self::Randomness {
                    rand: kzg10::Randomness::empty(),
                    shifted_rand: Some(kzg10::Randomness::empty())
                };
                commitments.len()
            ],
            |rands| rands,
        );

        let degree_bound = commitments[0].degree_bound();
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
        let mut aggregate_commitment = kzg10::Commitment::empty();
        let mut aggregate_shifted_commitment = kzg10::Commitment::empty();
        let mut aggregate_randomness = kzg10::Randomness::empty();
        let mut aggregate_shifted_randomness = kzg10::Randomness::empty();

        for (coef, term) in lc.iter() {
            match term {
                // No support for constant terms
                LCTerm::One => return Err(Error::ConstantTermInAggregation),

                // Find the corresponding commitment and randomness in our map; aggregate.
                LCTerm::PolyLabel(label) => match label_comm_mapping.get(label) {
                    Some((comm, rand)) => {
                        if degree_bound.is_some() {
                            aggregate_shifted_commitment += (
                                *coef,
                                &comm.shifted_comm.expect(
                                    "Degree bounded polynomial must have shifted commitment",
                                ),
                            );
                            aggregate_shifted_randomness += (
                                *coef,
                                &rand.shifted_rand.clone().expect(
                                    "Degree bounded polynomial must have shifted commitment",
                                ),
                            );
                        }

                        aggregate_commitment += (*coef, &comm.comm);
                        aggregate_randomness += (*coef, &rand.rand);
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

        let (shifted_comm, shifted_rand) = match degree_bound.is_some() {
            true => (
                Some(aggregate_shifted_commitment),
                Some(aggregate_shifted_randomness),
            ),
            false => (None, None),
        };

        let commitment = Self::Commitment {
            comm: aggregate_commitment,
            shifted_comm,
        };

        let randomness = Self::Randomness {
            rand: aggregate_randomness,
            shifted_rand,
        };

        Ok((
            LabeledCommitment::new(lc.label().clone(), commitment, degree_bound),
            randomness,
        ))
    }
}

#[cfg(test)]
mod test {
    use crate::{marlin_kzg::KZG10, AdditivelyHomomorphicPCS};
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
        let degree_bound = 10;

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
        let a_poly = LabeledPolynomial::new(
            String::from("a"),
            a_unlabeled,
            Some(degree_bound),
            Some(hiding_bound),
        );

        let b_unlabeled: DensePolynomial<F> = DensePolynomial::rand(5, rng);
        let b_poly = LabeledPolynomial::new(
            String::from("b"),
            b_unlabeled,
            Some(degree_bound),
            Some(hiding_bound),
        );

        let a_plus_2b_poly = a_poly.polynomial().clone() + (b_poly.polynomial() * F::from(2u64));
        let a_plus_2b_poly = LabeledPolynomial::new(
            String::from("a_plus_2b_poly"),
            a_plus_2b_poly,
            Some(degree_bound),
            Some(hiding_bound),
        );

        let polynomials = vec![a_poly.clone(), b_poly.clone()];
        let linear_combination =
            LinearCombination::new("a_plus_2b", vec![(F::one(), "a"), (F::from(2u64), "b")]);

        // Commit Phase
        let (commitments, rands) = PC::commit(&ck, &polynomials, Some(rng)).unwrap();
        let (test_commitment, test_rand) =
            PC::aggregate_commitments(&commitments, Some(rands.to_vec()), &linear_combination)
                .unwrap();

        for comm in &commitments {
            println!("{}: {:?}", comm.label(), comm.commitment().shifted_comm);
        }
        println!(
            "{}: {:?}",
            test_commitment.label(),
            test_commitment.commitment().shifted_comm
        );

        let mut manual_commitment = commitments[0].commitment().comm;
        manual_commitment += (F::from(2u64), &commitments[1].commitment().comm);

        let mut manual_rand = rands[0].clone().rand;
        manual_rand += (F::from(2u64), &rands[1].clone().rand);

        assert_eq!(test_commitment.commitment().comm, manual_commitment);
        assert_eq!(test_rand.rand, manual_rand);

        println!("PASSED SANITY");

        // Derive evaluation point and generate a query set
        let evaluation_point = Fr::rand(rng);

        // Evaluation Phase, here we only output the evaluation of the linear combination
        let manual_eval = a_plus_2b_poly.evaluate(&evaluation_point);

        let opening_challenge = F::rand(rng);

        // Opening phase
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
