use crate::commitment::HomomorphicPolynomialCommitment;
use crate::error::{to_pc_error, Error};
use crate::geo_seq::proof::Proof;
use crate::label_polynomial;
use crate::util::generate_sequence;
use crate::virtual_oracle::geometric_sequence_vo::GeoSequenceVO;
use crate::virtual_oracle::normalized_vo::{Description, Term};
use crate::virtual_oracle::{EvaluationsProvider, VirtualOracle};
use crate::zero_over_k::proof::Proof as Z_Proof;
use crate::zero_over_k::ZeroOverK;
use ark_bn254::Fr;
use ark_ff::{to_bytes, PrimeField};
use ark_marlin::rng::FiatShamirRng;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};
use ark_poly_commit::{Evaluations, LabeledCommitment, LabeledPolynomial, QuerySet};
use ark_std::rand::thread_rng;
use digest::Digest; // Note that in the latest Marlin commit, Digest has been replaced by an arkworks trait `FiatShamirRng`
use rand::Rng;
use rand_core::OsRng;

mod proof;
mod tests;

pub struct GeoSeqTest<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> {
    _zero_over_k: ZeroOverK<F, PC, D>,
}

impl<F: PrimeField, PC: HomomorphicPolynomialCommitment<F>, D: Digest> GeoSeqTest<F, PC, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"Geometric Sequence Test";
    // TODO: for both prove() and verify:
    // TODO: degree bounds!
    // TODO: have an assertion that domain is large enough given m
    // TODO: move the padding outside and the check that the length is correct
    // TODO: verifier should check that the size of the sequence is correct given the domain
    fn prove<R: Rng>(
        ck: &PC::CommitterKey,
        r: F,
        a_s: &Vec<F>,
        c_s: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        rng: &mut R,
    ) -> Result<Proof<F, PC>, Error> {
        let seq = generate_sequence::<F>(r, &a_s.as_slice(), &c_s.as_slice());

        // Generate f() such that f(w^n) = a_i*r^n
        let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&seq));
        let f = label_polynomial!(f);
        // Generate the GeoSequenceVO virtual oracle
        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);

        let concrete_oracles = [f.clone()];
        let alphas = [F::from(1u64), domain.element(1)];

        let (concrete_oracle_commitments, concrete_oracle_rands) =
            PC::commit(ck, &concrete_oracles, None).unwrap();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(
            &to_bytes![
                &Self::PROTOCOL_NAME,
                a_s,
                c_s.iter().map(|&x| x as u64).collect::<Vec<_>>(),
                r,
                &concrete_oracle_commitments.to_vec(),
                &alphas.to_vec()
            ]
            .unwrap(),
        );

        let mut query_set = QuerySet::new();
        let pi_s = geo_seq_vo.get_pi_s();
        for (i, &pi) in pi_s.iter().enumerate() {
            query_set.insert((
                f.label().clone(),
                (
                    format!("gamma_pi_{}", i),
                    domain.element(1).pow([pi as u64]),
                ),
            ));
        }

        let separation_challenge = F::rand(&mut fs_rng);
        let opening_proof = PC::batch_open(
            ck,
            &concrete_oracles,
            &concrete_oracle_commitments,
            &query_set,
            separation_challenge,
            &concrete_oracle_rands,
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let z_proof = ZeroOverK::<F, PC, D>::prove(
            &concrete_oracles,
            &concrete_oracle_commitments,
            &concrete_oracle_rands,
            &geo_seq_vo,
            &alphas.to_vec(),
            &domain,
            &ck,
            rng,
        )?;

        let proof = Proof::<F, PC> {
            z_proof,
            f_commit: concrete_oracle_commitments[0].commitment().clone(),
            opening_proof,
        };
        Ok(proof)
    }

    pub fn verify(
        r: F,
        a_s: &Vec<F>,
        c_s: &Vec<usize>,
        domain: &GeneralEvaluationDomain<F>,
        proof: Proof<F, PC>,
        vk: &PC::VerifierKey,
    ) -> Result<(), Error> {
        let alphas = [F::from(1u64), domain.element(1)];
        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);
        let concrete_oracle_commitments = [LabeledCommitment::new(
            String::from("f"),
            proof.f_commit,
            None,
        )];

        // Test that for all i in n, check that f(gamma^p_i) = a_i
        let pi_s = geo_seq_vo.get_pi_s();
        let points = pi_s
            .iter()
            .map(|&pi| domain.element(1).pow([pi as u64]))
            .collect::<Vec<_>>();

        let mut fs_rng = FiatShamirRng::<D>::from_seed(
            &to_bytes![
                &Self::PROTOCOL_NAME,
                a_s,
                c_s.iter().map(|&x| x as u64).collect::<Vec<_>>(),
                r,
                &concrete_oracle_commitments.to_vec(),
                &alphas.to_vec()
            ]
            .unwrap(),
        );

        let mut query_set = QuerySet::new();
        for (i, &point_i) in points.iter().enumerate() {
            query_set.insert((String::from("f"), (format!("gamma_pi_{}", i), point_i)));
        }
        let mut evaluations = ark_poly_commit::Evaluations::new();
        for (&point_i, &a_i) in points.iter().zip(a_s.iter()) {
            evaluations.insert((String::from("f"), point_i), a_i);
        }

        let separation_challenge = F::rand(&mut fs_rng);
        match PC::batch_check(
            vk,
            &concrete_oracle_commitments,
            &query_set,
            &evaluations,
            &proof.opening_proof,
            separation_challenge,
            &mut OsRng,
        ) {
            Ok(true) => Ok(()),
            Ok(false) => Err(Error::ProofVerificationError),
            Err(e) => panic!("{:?}", e),
        }?;

        // let seq = generate_sequence::<F>(r, &a_s.as_slice(), &c_s.as_slice());
        // let f = DensePolynomial::<F>::from_coefficients_slice(&domain.ifft(&seq));
        ZeroOverK::<F, PC, D>::verify(
            proof.z_proof,
            &concrete_oracle_commitments,
            &geo_seq_vo,
            &domain,
            &alphas,
            vk,
        )
    }

    //pub fn prove2(
    //seq: &Vec<F>,
    //r: F,
    //a_s: &[F],
    //c_s: &[usize],
    //) -> Result<Proof<F, PC>, Error>{
    //// Check that seq is valid; otherwise, return an error
    //if !GeoSeqTest::<F, PC, D>::naive_verify(&seq, r, a_s, c_s) {
    //return Err(Error::InvalidGeoSeq);
    //}

    //// Generate f() such that f(w^n) = a_i*r^n
    //let m: usize = c_s.iter().sum();
    //let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
    //let to_pad: usize = domain.size() - m;

    //// Generate the GeoSequenceVO virtual oracle
    //let coeffs = domain.ifft(seq);
    //let coeffs_len = coeffs.len();
    //let f = DensePolynomial::<F>::from_coefficients_vec(coeffs);
    //let mut a_vec: Vec<F> = a_s.to_vec();
    //let mut c_vec = c_s.to_vec();

    //if to_pad > 0 {
    //a_vec.push(F::from(0u64));
    //c_vec.push(to_pad);
    //}

    //let concrete_oracles = [label_polynomial!(f)];
    //let alphas = [F::from(1u64), domain.element(1)];

    //let geo_seq_vo = GeoSequenceVO::new(&c_vec, domain.element(1), r);

    ////let instantiatied_geo_seq_vo = geo_seq_vo
    ////.instantiate(
    ////&concrete_oracles,
    ////&[F::from(1u64), domain.element(1)],
    ////)
    ////.unwrap();
    ////for root_of_unity in domain.elements() {
    ////let eval = instantiatied_geo_seq_vo.evaluate(&root_of_unity);
    ////assert_eq!(eval, F::from(0u64));
    ////}

    //let mut rng = thread_rng();

    //// perform zero over k
    //let max_degree = 80;
    //let pp = PC::setup(max_degree, None, &mut thread_rng()).unwrap();
    //let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();

    //let (concrete_oracle_commitments, concrete_oracle_rands) =
    //PC::commit(&ck, &concrete_oracles, None).unwrap();

    //let z_proof = ZeroOverK::<F, PC, D>::prove(
    //&concrete_oracles,
    //&concrete_oracle_commitments,
    //&concrete_oracle_rands,
    //&geo_seq_vo,
    //&alphas.to_vec(),
    //&domain,
    //&ck,
    //&mut rng,
    //);

    //let proof = Proof::<F, PC> {
    //z_proof: z_proof.unwrap(),
    //concrete_oracle_commitments,
    //seq: seq.clone(),
    //r: r,
    //a_s: a_vec,
    //c_s: c_vec,
    //};
    //Ok(proof)
    //}

    /// Inefficiently verify that the sequence is valid
    pub fn naive_verify(seq: &Vec<F>, r: F, a_s: &[F], c_s: &[usize]) -> bool {
        if a_s.len() != c_s.len() {
            return false;
        }

        if c_s.iter().fold(0, |x, y| x + y) != seq.len() {
            return false;
        }

        let mut i = 0;
        for (a_i, a) in a_s.iter().enumerate() {
            let mut r_pow = F::from(1 as u64);
            for _ in 0..c_s[a_i] {
                let expected = *a * r_pow;

                if expected != seq[i] {
                    return false;
                }

                r_pow = r * r_pow;
                i += 1;
            }
        }

        return true;
    }
}
