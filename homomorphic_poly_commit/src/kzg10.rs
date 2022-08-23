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
    fn get_commitments_lc(
        commitments: &[LabeledCommitment<Self::Commitment>],
        lc: &LinearCombination<E::Fr>,
    ) -> Result<LabeledCommitment<Self::Commitment>, Error> {
        let mut aggregate_commitment = Self::Commitment::empty();

        let degree_bound = commitments[0].degree_bound();
        for comm in commitments {
            if comm.degree_bound() != degree_bound {
                // Can only accumulate polynomials and commitments that have the same degree bound
                return Err(Error::MismatchedDegreeBounds(format!(
                    "{} has degree bound {:?}, but {} has degree bound {:?}",
                    commitments[0].label(),
                    degree_bound,
                    comm.label(),
                    comm.degree_bound()
                )));
            }
        }

        for (coef, term) in lc.iter() {
            let commitment = if let LCTerm::PolyLabel(label) = term {
                commitments
                    .iter()
                    .find(|&c| c.label() == label)
                    .ok_or(Error::MissingCommitment(format!(
                        "Could not find object with label '{}' when computing '{}'",
                        label,
                        lc.label()
                    )))?
                    .commitment()
                    .clone()
            } else {
                Self::Commitment::empty()
            };
            aggregate_commitment += (*coef, &commitment);
        }

        Ok(LabeledCommitment::new(
            lc.label().clone(),
            aggregate_commitment,
            degree_bound,
        ))
    }

    fn get_commitments_lc_with_rands(
        commitments: &[LabeledCommitment<Self::Commitment>],
        hiding_rands: &[Self::Randomness],
        lc: &LinearCombination<E::Fr>,
    ) -> Result<(LabeledCommitment<Self::Commitment>, Self::Randomness), Error> {
        if commitments.len() != hiding_rands.len() {
            return Err(Error::InputLengthError(format!(
                "There are {} commitments and {} randomness values",
                commitments.len(),
                hiding_rands.len()
            )));
        }

        let degree_bound = commitments[0].degree_bound();
        for comm in commitments {
            if comm.degree_bound() != degree_bound {
                // Can only accumulate polynomials and commitments that have the same degree bound
                return Err(Error::MismatchedDegreeBounds(format!(
                    "{} has degree bound {:?}, but {} has degree bound {:?}",
                    commitments[0].label(),
                    degree_bound,
                    comm.label(),
                    comm.degree_bound()
                )));
            }
        }

        let mut aggregate_commitment = Self::Commitment::empty();
        let mut aggregate_randomness = Self::Randomness::empty();

        for (coef, term) in lc.iter() {
            let (comm, rand) = if let LCTerm::PolyLabel(label) = term {
                let current_pair = commitments
                    .iter()
                    .zip(hiding_rands.iter())
                    .find(|&c| c.0.label() == label)
                    .ok_or(Error::MissingCommitment(format!(
                        "Could not find object with label '{}' when computing '{}'",
                        label,
                        lc.label()
                    )))?;
                (current_pair.0.commitment().clone(), current_pair.1.clone())
            } else {
                println!("I'm here");
                (Self::Commitment::empty(), Self::Randomness::empty())
            };
            aggregate_commitment += (*coef, &comm);
            aggregate_randomness += (*coef, &rand);
        }

        Ok((
            LabeledCommitment::new(lc.label().clone(), aggregate_commitment, degree_bound),
            aggregate_randomness,
        ))
    }
}

#[cfg(test)]
mod test {
    use crate::{kzg10::KZG10, AdditivelyHomomorphicPCS};
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
            PC::get_commitments_lc_with_rands(&commitments, &rands, &linear_combination).unwrap();

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
