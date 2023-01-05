use ark_ff::{FftField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};

// given the x coords construct Li polynomials
pub fn construct_lagrange_basis<F: FftField>(evaulation_domain: &[F]) -> Vec<DensePolynomial<F>> {
    let mut bases = Vec::with_capacity(evaulation_domain.len());
    for i in 0..evaulation_domain.len() {
        let mut l_i = DensePolynomial::from_coefficients_slice(&[F::one()]);
        let x_i = evaulation_domain[i];
        for j in 0..evaulation_domain.len() {
            if j != i {
                let nom =
                    DensePolynomial::from_coefficients_slice(&[-evaulation_domain[j], F::one()]);
                let denom = x_i - evaulation_domain[j];

                l_i = &l_i * &(&nom * denom.inverse().unwrap());
            }
        }

        bases.push(l_i);
    }

    bases
}

pub fn construct_vanishing<F: FftField>(evaulation_domain: &[F]) -> DensePolynomial<F> {
    let mut v_h = DensePolynomial::from_coefficients_slice(&[F::one()]);
    for point in evaulation_domain {
        v_h = &v_h * &DensePolynomial::from_coefficients_slice(&[-*point, F::one()]);
    }

    v_h
}

// pub fn well_formation<F: PrimeField>(concrete_terms: &[VOTerm<F>]) -> VOTerm<F> {
//     concrete_terms[0].clone()
//         + vo_constant!(F::from(2u64)) * concrete_terms[2].clone()
//         + vo_constant!(F::from(3u64)) * concrete_terms[3].clone()
// }

#[cfg(test)]
mod tests {
    use super::{construct_lagrange_basis, construct_vanishing};
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, Evaluations,
        GeneralEvaluationDomain, MixedRadixEvaluationDomain, Polynomial, UVPolynomial,
    };

    type F = Fr;

    #[test]
    fn test_assignment_formation() {
        let instance_len = 4;
        let wnts_len = 4;

        let h_len = instance_len + wnts_len;

        let domain_h = GeneralEvaluationDomain::<F>::new(h_len).unwrap();
        let _domain_x = GeneralEvaluationDomain::<F>::new(instance_len).unwrap();

        let mut instance_evals = vec![F::from(5u64), F::from(4u64), F::from(3u64), F::from(6u64)];
        let _wtns_evals = vec![
            F::from(12u64),
            F::from(15u64),
            F::from(13u64),
            F::from(56u64),
        ];

        // start with x - 1 -> there is alway at least one PI
        let mut x_vanish = DensePolynomial::from_coefficients_slice(&[-F::one(), F::one()]);
        for i in 1..instance_evals.len() {
            x_vanish = &x_vanish
                * &DensePolynomial::from_coefficients_slice(&[-domain_h.element(i), F::one()]);
        }

        instance_evals.extend(vec![F::zero(); domain_h.size() - instance_len]);
        for x in &instance_evals {
            println!("x: {}", x);
        }

        let _x_poly =
            Evaluations::from_vec_and_domain(instance_evals.clone(), domain_h).interpolate();
        // let x_evals = domain_h.fft(&x_poly);

        // println!("x deg: {}", x_poly.degree());

        let mixed_domain = MixedRadixEvaluationDomain::<F>::new(4).unwrap();
        println!("{}", mixed_domain.size());
    }

    #[test]
    fn encode_assignment() {
        let instance_len = 4;
        let wnts_len = 4;

        let h_len = instance_len + wnts_len;

        let domain_h = GeneralEvaluationDomain::<F>::new(h_len).unwrap();
        let domain_x = GeneralEvaluationDomain::<F>::new(instance_len).unwrap();

        let instance_evals = vec![F::from(5u64), F::from(4u64), F::from(3u64), F::from(6u64)];
        let wtns_evals = vec![
            F::from(12u64),
            F::from(15u64),
            F::from(13u64),
            F::from(56u64),
        ];

        let x_poly =
            Evaluations::from_vec_and_domain(instance_evals.clone(), domain_x).interpolate();
        let x_evals = domain_h.fft(&x_poly);

        // for x in x_evals {
        //     println!("x: {}", x);
        // }

        let ratio = domain_h.size() / domain_x.size();

        let w_poly_evals = (0..domain_h.size())
            .map(|k| {
                if k % ratio == 0 {
                    F::zero()
                } else {
                    wtns_evals[k - (k / ratio) - 1] - instance_evals[k - (k / ratio) - 1]
                }
            })
            .collect::<Vec<F>>();

        let w_poly_evals_like_marlin = (0..domain_h.size())
            .map(|k| {
                if k % ratio == 0 {
                    F::zero()
                } else {
                    wtns_evals[k - (k / ratio) - 1] - x_evals[k]
                }
            })
            .collect::<Vec<F>>();

        // assert_eq!(w_poly_evals, w_poly_evals_like_marlin);

        // for (w, w_m) in w_poly_evals.iter().zip(w_poly_evals_like_marlin.iter()) {
        //     println!("w: {}, w_m: {}", w, w_m);
        // }

        let w_poly = DensePolynomial::<F>::from_coefficients_slice(&domain_h.ifft(&w_poly_evals));
        let (_w_poly, remainder) = w_poly.divide_by_vanishing_poly(domain_x).unwrap();
        assert!(remainder.is_zero());

        let w_poly_marlin = DensePolynomial::<F>::from_coefficients_slice(
            &domain_h.ifft(&w_poly_evals_like_marlin),
        );
        let (w_poly_marlin, remainder) = w_poly_marlin.divide_by_vanishing_poly(domain_x).unwrap();
        assert!(remainder.is_zero());

        let mut assignment = vec![];
        assignment.extend(instance_evals.iter());
        assignment.extend(wtns_evals.iter());

        for a in &assignment {
            println!("{}", a);
        }

        let z_poly = DensePolynomial::<F>::from_coefficients_slice(&assignment);

        let w_vh_x = w_poly_marlin.mul_by_vanishing_poly(domain_x);
        let z_poly_by_hand = w_vh_x + x_poly;

        assert_eq!(z_poly, z_poly_by_hand);

        // for h in domain_x.elements() {
        //     let rhs = w_poly.evaluate(&h) * domain_x.evaluate_vanishing_polynomial(h) + x_poly.evaluate(&h);
        //     assert_eq!(
        //         z_poly.evaluate(&h),
        //         rhs
        //     );
        // }
    }

    #[test]
    fn test_artificial_domain() {}

    #[test]
    fn test_lagrange_basis() {
        let instance_len = 4;
        let wnts_len = 4;

        let h_len = instance_len + wnts_len;
        let domain_h = GeneralEvaluationDomain::<F>::new(h_len).unwrap();

        let instance_evals = vec![F::from(5u64), F::from(4u64), F::from(3u64), F::from(6u64)];
        let wtns_evals = vec![
            F::from(12u64),
            F::from(15u64),
            F::from(13u64),
            F::from(56u64),
        ];

        let elems: Vec<F> = domain_h.elements().collect();
        let pi_roots_of_unity: Vec<F> = elems.iter().take(instance_len).map(|&x| x).collect();
        let wtns_roots_of_unity: Vec<F> = elems.iter().skip(instance_len).map(|&x| x).collect();

        let assginment_evals: Vec<F> = instance_evals
            .iter()
            .chain(wtns_evals.iter())
            .map(|&x| x)
            .collect();
        let assignment =
            DensePolynomial::from_coefficients_slice(&domain_h.ifft(&assginment_evals));

        let bases = construct_lagrange_basis(&pi_roots_of_unity);

        let mut x_poly = DensePolynomial::<F>::zero();
        for (eval, base) in bases.iter().zip(instance_evals.iter()) {
            x_poly += &(eval * *base);
        }

        let x_evals = domain_h.fft(&x_poly);
        for x_eval in &x_evals {
            println!("{}", x_eval);
        }

        println!("x deg: {}", x_poly.degree());

        let v_h = construct_vanishing(&wtns_roots_of_unity);

        for w in &elems {
            let assignment_eval = assignment.evaluate(w);
            let v_h_eval = v_h.evaluate(w);
            let x_eval = x_poly.evaluate(w);

            assert_eq!(F::zero(), v_h_eval * (assignment_eval - x_eval));
        }
    }
}
