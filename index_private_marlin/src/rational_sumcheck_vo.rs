use ::zero_over_k::virtual_oracle::VirtualOracle;

use crate::{virtual_oracle::Error, virtual_oracle::VirtualOracle};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
use ark_poly_commit::LabeledPolynomial;

pub struct RationalSumcheckVO<F: PrimeField> {
    pub eta_a: F, 
    pub eta_b: F, 
    pub eta_c: F, 

    pub alpha: F, 
    pub beta: F,

    pub vh_alpha: F, 
    pub vh_beta: F
}

impl<F: PrimeField> RationalSumcheckVO<F> {
    pub fn new(
        eta_a: F, 
        eta_b: F,
        eta_c: F,
        alpha: F,
        beta: F,
        vh_alpha: F,
        vh_beta: F,
    ) -> Self {
        Self { 
            eta_a, 
            eta_b, 
            eta_c, 
            alpha, 
            beta, 
            vh_alpha, 
            vh_beta 
        }
    }
}

// Oracles are in order [row_a, col_a, val_a, row_b, col_b, val_b, row_c, col_c, val_c, f]
//                         0     1      2       3      4      5      6      7      8    9
// impl<F: PrimeField> VirtualOracle<F> for RationalSumcheckVO<F> {
//     fn mapping_vector(&self) -> Vec<usize> {
//         vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
//     }

//     fn shifting_coefficients(&self) -> Vec<F> {
//         vec![F::One(); 10]
//     }

//     fn instantiate_in_coeffs_form(
//         &self,
//         concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
//         _alphas: &[F],
//     ) -> Result<DensePolynomial<F>, Error> {
//         if concrete_oracles.len() != 10 {
//             return Err(Error::InstantiationError);
//         }

//         let v_H_alpha_v_H_beta = self.vh_alpha * self.vh_beta;
//         let eta_a_times_v_H_alpha_v_H_beta = self.eta_a * v_H_alpha_v_H_beta;
//         let eta_b_times_v_H_alpha_v_H_beta = self.eta_b * v_H_alpha_v_H_beta;
//         let eta_c_times_v_H_alpha_v_H_beta = self.eta_c * v_H_alpha_v_H_beta;

//         let a_row = concrete_oracles[0].polynomial();
//         let a_col = concrete_oracles[1].polynomial();
//         let a_val = concrete_oracles[2].polynomial();
//         let b_row = concrete_oracles[3].polynomial();
//         let b_col = concrete_oracles[4].polynomial();
//         let b_val = concrete_oracles[5].polynomial();
//         let c_row = concrete_oracles[6].polynomial();
//         let c_col = concrete_oracles[7].polynomial();
//         let c_val = concrete_oracles[8].polynomial();
//         let f = concrete_oracles[9].polynomial();

//         let a_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[self.beta])) - a_row)
//             * &(&(DensePolynomial::from_coefficients_slice(&[self.alpha])) - a_col);

//         let b_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[self.beta])) - b_row)
//             * &(&(DensePolynomial::from_coefficients_slice(&[self.alpha])) - b_col);

//         let c_part_denom = &(&(DensePolynomial::from_coefficients_slice(&[self.beta])) - c_row)
//             * &(&(DensePolynomial::from_coefficients_slice(&[self.alpha])) - c_col);

//         let b_poly = &(&a_part_denom * &b_part_denom) * &c_part_denom;

//         let a_part_nom =
//             &(DensePolynomial::from_coefficients_slice(&[eta_a_times_v_H_alpha_v_H_beta])) * a_val;
//         let b_part_nom =
//             &(DensePolynomial::from_coefficients_slice(&[eta_b_times_v_H_alpha_v_H_beta])) * b_val;
//         let c_part_nom =
//             &(DensePolynomial::from_coefficients_slice(&[eta_c_times_v_H_alpha_v_H_beta])) * c_val;

//         let a_poly = {
//             let summand_0 = &a_part_nom * &(&b_part_denom * &c_part_denom);
//             let summand_1 = &b_part_nom * &(&a_part_denom * &c_part_denom);
//             let summand_2 = &c_part_nom * &(&a_part_denom * &b_part_denom);

//             summand_0 + summand_1 + summand_2
//         };

//         Ok(&a_poly - &(&b_poly - f))
//     }

//     fn instantiate_in_evals_form(
//         &self,
//         concrete_oracles: &[LabeledPolynomial<F, DensePolynomial<F>>],
//         _alphas: &[F],
//         domain: &GeneralEvaluationDomain<F>,
//     ) -> Result<Vec<F>, Error> {
//         if concrete_oracles.len() != 10 {
//             return Err(Error::InstantiationError);
//         }

//         let n = domain.size();
//         let k = self.compute_scaling_factor(domain);

//         let domain_kn = GeneralEvaluationDomain::<F>::new(k * n).unwrap();

//         let a_row_evals = domain_kn.coset_fft(concrete_oracles[0].polynomial());
//         let a_col_evals = domain_kn.coset_fft(concrete_oracles[1].polynomial());
//         let a_val_evals = domain_kn.coset_fft(concrete_oracles[2].polynomial());
//         let b_row_evals = domain_kn.coset_fft(concrete_oracles[3].polynomial());
//         let b_col_evals = domain_kn.coset_fft(concrete_oracles[4].polynomial());
//         let b_val_evals = domain_kn.coset_fft(concrete_oracles[5].polynomial());
//         let c_row_evals = domain_kn.coset_fft(concrete_oracles[6].polynomial());
//         let c_col_evals = domain_kn.coset_fft(concrete_oracles[7].polynomial());
//         let c_val_evals = domain_kn.coset_fft(concrete_oracles[8].polynomial());
//         let f_evals = domain_kn.coset_fft(concrete_oracles[9].polynomial());

//         let vh_alpha_beta = self.vh_alpha * self.vh_beta;
//         let alpha_beta = self.alpha * self.beta;

//         let vo_evals = (0..domain_kn.size())
//             .map(|i| {
//                 let a_denom = alpha_beta - self.beta*a_col_evals[i] - self.alpha * a_row_evals[i] + a_col_evals[i]*a_row_evals[i];
//                 let b_denom = alpha_beta - self.beta*b_col_evals[i] - self.alpha * b_row_evals[i] + b_col_evals[i]*b_row_evals[i];
//                 let c_denom = alpha_beta - self.beta*c_col_evals[i] - self.alpha * c_row_evals[i] + c_col_evals[i]*c_row_evals[i];

//                 let a_nom = self.eta_a * vh_alpha_beta * a_val_evals[i];
//                 let b_nom = self.eta_b * vh_alpha_beta * b_val_evals[i];
//                 let c_nom = self.eta_c * vh_alpha_beta * c_val_evals[i];

//                 let b_eval = a_denom * b_denom * c_denom;

//                 let a_eval = {
//                     let summand_0 = a_nom * b_denom * c_denom;
//                     let summand_1 = b_nom * a_denom * b_denom;
//                     let summand_2 = c_nom * a_denom * c_denom;
        
//                     summand_0 + summand_1 + summand_2
//                 };

//                 a_eval - b_eval*f_evals[i]
//             })
//             .collect::<Vec<_>>();

//         Ok(vo_evals)
//     }

//     fn query(&self, evals: &[F], _point: F) -> Result<F, Error> {
//         if evals.len() != 10 {
//             return Err(Error::EvaluationError);
//         }

//         let vh_alpha_beta = self.vh_alpha * self.vh_beta;
//         let alpha_beta = self.alpha * self.beta;

//         let a_denom = alpha_beta - self.beta*evals[1] - self.alpha * evals[0] + evals[0]*evals[1];
//         let b_denom = alpha_beta - self.beta*evals[4] - self.alpha * evals[3] + evals[3]*evals[4];
//         let c_denom = alpha_beta - self.beta*evals[7] - self.alpha * evals[6] + evals[6]*evals[7];

//         let a_nom = self.eta_a * vh_alpha_beta * evals[2];
//         let b_nom = self.eta_b * vh_alpha_beta * evals[5];
//         let c_nom = self.eta_c * vh_alpha_beta * evals[8];

//         let b_eval = a_denom * b_denom * c_denom;

//         let a_eval = {
//             let summand_0 = a_nom * b_denom * c_denom;
//             let summand_1 = b_nom * a_denom * b_denom;
//             let summand_2 = c_nom * a_denom * c_denom;

//             summand_0 + summand_1 + summand_2
//         };

//         Ok(a_eval - b_eval * evals[9])
//     }

//     fn num_of_oracles(&self) -> usize {
//         return 10;
//     }

//     /// this map encodes at which concrete oracle should h_i point
//     fn mapping_vector(&self) -> Vec<usize> {
//         // h0 = f0, h1 = f1
//         Vec::from([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
//     }

//     // fn degree_bound(&self, domain_size: usize) -> usize {
//     //     6 * domain_size + 7
//     // }

//     // fn name(&self) -> String {
//     //     String::from("rational_sumcheck")
//     // }
// }