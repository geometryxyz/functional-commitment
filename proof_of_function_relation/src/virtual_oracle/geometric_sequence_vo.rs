use crate::error::Error;
use crate::to_poly;
use crate::util::shift_dense_poly;
use crate::virtual_oracle::VirtualOracle;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, UVPolynomial, GeneralEvaluationDomain, EvaluationDomain};
use ark_poly_commit::LabeledPolynomial;
use std::iter;

pub struct GeoSequenceVO<F: PrimeField> {
    pub pi_s: Vec<usize>,
    pub ci_s: Vec<usize>,
    gamma: F,
    r: F,
}

impl<F: PrimeField> GeoSequenceVO<F> {
    pub fn new(ci_s: &Vec<usize>, gamma: F, r: F) -> Self {
        let pi_s = iter::once(0)
            .chain(ci_s.iter().scan(0, |st, elem| {
                *st += elem;
                Some(*st)
            }))
            .collect::<Vec<_>>();

        Self {
            pi_s: pi_s.clone(),
            ci_s: ci_s.clone(),
            gamma,
            r,
        }
    }

    pub fn get_pi_s(&self) -> Vec<usize> {
        self.pi_s.clone()
    }
}

impl<F: PrimeField> VirtualOracle<F> for GeoSequenceVO<F> {
    fn instantiate_in_coeffs_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
    ) -> Result<DensePolynomial<F>, Error> {
        if concrete_oracles.len() != 2 || alphas.len() != 2 {
            return Err(Error::InstantiationError);
        }

        // construct (f(gamma * x) - r * f(x))
        let mut instantiation_poly = shift_dense_poly(&concrete_oracles[1], &alphas[1])
            + (&shift_dense_poly(&concrete_oracles[0], &alphas[0]) * -self.r);

        let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
            // construct x - y^(pi + ci - 1)
            let stitch_i = &x_poly + &to_poly!(-self.gamma.pow([(pi + ci - 1) as u64]));
            // let stitch_i = x_poly.clone();
            instantiation_poly = &instantiation_poly * &stitch_i;
        }

        Ok(instantiation_poly)
    }

    fn instantiate_in_evals_form(
        &self,
        concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
        alphas: &[F],
        domain: &GeneralEvaluationDomain<F>,
    ) -> Result<Vec<F>, Error> {
        if concrete_oracles.len() != 2 || alphas.len() != 2 {
            return Err(Error::InstantiationError);
        }

        let n = domain.size();
        let k = self.compute_scaling_factor(domain);

        let domain_kn = GeneralEvaluationDomain::<F>::new(k * n).unwrap();

        //TODO REMOVE THIS
        // let mut instantiation_poly = shift_dense_poly(&concrete_oracles[1], &alphas[1])
        //     + (&shift_dense_poly(&concrete_oracles[0], &alphas[0]) * -self.r);

        // let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        // for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
        //     // construct x - y^(pi + ci - 1)
        //     let stitch_i = &x_poly + &to_poly!(-self.gamma.pow([(pi + ci - 1) as u64]));
        //     // let stitch_i = x_poly.clone();
        //     instantiation_poly = &instantiation_poly * &stitch_i;
        // }

        // println!("DEG: {}", instantiation_poly.degree());
        /////////

        // let f_evals = domain_kn.coset_fft(concrete_oracles[0].polynomial());
        let f_evals = domain_kn.coset_fft(&shift_dense_poly(concrete_oracles[0].polynomial(), &alphas[0]));
        let f_sh_evals = domain_kn.coset_fft(&shift_dense_poly(concrete_oracles[1].polynomial(), &alphas[1]));


        let x_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero(), F::one()]);
        let x_evals = domain_kn.coset_fft(&x_poly);

        let vo_evals = (0..domain_kn.size())
            .map(|i| {
                let mut seq_part = f_sh_evals[i] - self.r * f_evals[i];

                for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
                    // construct x - y^(pi + ci - 1)
                    let stitch_i = x_evals[i] - self.gamma.pow([(pi + ci - 1) as u64]);
                    // let stitch_i = domain.element(i);
                    seq_part *= stitch_i;
                }

                seq_part

            })
            .collect::<Vec<_>>();

        Ok(vo_evals)
    }

    fn num_of_oracles(&self) -> usize {
        return 2;
    }

    fn query(&self, evals: &[F], point: F) -> Result<F, Error> {
        if evals.len() != 2 {
            return Err(Error::EvaluationError);
        }

        // evals = [shifted_f, f]
        let mut eval = evals[1] - self.r * evals[0];
        for (&pi, &ci) in self.pi_s.iter().zip(self.ci_s.iter()) {
            // construct x - y^(pi + ci - 1)
            let stitch_i = point - self.gamma.pow([(pi + ci - 1) as u64]);
            eval *= stitch_i;
        }

        Ok(eval)
    }

    /// this map encodes at which concrete oracle should h_i point
    fn mapping_vector(&self) -> Vec<usize> {
        // in geo sequence test we have just one concrete oracle f0, hence
        // h0 = f0, h1 = f0
        Vec::from([0, 0])
    }

    fn degree_bound(&self, domain_size: usize) -> usize {
        let n = self.ci_s.len();
        domain_size + 1 + n
    }

    // fn compute_scaling_factor(&self, _domain: &GeneralEvaluationDomain<F>) -> usize {
    //     2
    // }

    fn name(&self) -> String {
        String::from("geo_seq")
    }
}

#[cfg(test)]
mod test {

    use super::{GeoSequenceVO, VirtualOracle};
    use crate::{label_polynomial, util::{generate_sequence, compute_vanishing_poly_over_coset}};
    use ark_bn254::Fr;
    use ark_poly::{
        univariate::{DensePolynomial, DenseOrSparsePolynomial}, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_ff::{Field, Zero};

    /// @param a_s The initial value of each sequence

    #[test]
    fn test_geo_seq_instantiation() {
        let r = Fr::from(2u64);
        let mut a_s = vec![
            Fr::from(2u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(11u64),
        ];
        let mut c_s = vec![5, 3, 10, 30];

        let m = c_s.iter().sum();

        // domain is a subset of the scalar field F.
        // it's significant because domain.element(1) = gamma, domain.element(2) = gamma^2, etc
        // m is also called the order of some subgroup (TBD)
        // m is the length that we want (it must be a power of 2 if you want to ffts and iffts; otherwise
        // you don't need it to be so).
        let domain = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let to_pad = domain.size() - m;
        if to_pad > 0 {
            a_s.push(Fr::from(0u64));
            c_s.push(to_pad);
        }

        let seq = generate_sequence::<Fr>(r, &a_s.as_slice(), &c_s.as_slice());
        let f = DensePolynomial::<Fr>::from_coefficients_slice(&domain.ifft(&seq));

        let geo_seq_vo = GeoSequenceVO::new(&c_s, domain.element(1), r);

        let k = geo_seq_vo.compute_scaling_factor(&domain);
        let domain_kn = GeneralEvaluationDomain::<Fr>::new(k * domain.size()).unwrap();
        let vh = compute_vanishing_poly_over_coset(&domain_kn, domain.size() as u64);

        let geo_seq_vo_evals = geo_seq_vo
            .instantiate_in_evals_form(
                &[label_polynomial!(f), label_polynomial!(f)],
                &[Fr::from(1u64), domain.element(1)],
                &domain
            )
            .unwrap();

        let quotient_evals = geo_seq_vo_evals
            .iter()
            .zip(vh.evals.iter())
            .map(|(&nominator_eval, &denominator_eval)| {
                nominator_eval * denominator_eval.inverse().unwrap()
            })
            .collect::<Vec<_>>();

        let quotient =
            DensePolynomial::from_coefficients_slice(&domain_kn.coset_ifft(&quotient_evals));

        
        let geo_seq_vo_coeffs = geo_seq_vo
            .instantiate_in_coeffs_form(
                &[label_polynomial!(f), label_polynomial!(f)],
                &[Fr::from(1u64), domain.element(1)],
            )
            .unwrap();

        let (q_1, _r) = DenseOrSparsePolynomial::from(&geo_seq_vo_coeffs)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(
                &domain.vanishing_polynomial(),
            ))
            .unwrap();

        assert_eq!(_r, DensePolynomial::<Fr>::zero());
        println!("div deg: {}, eval deg: {}", q_1.degree(), quotient.degree());
        assert_eq!(q_1, quotient);
    }
}
