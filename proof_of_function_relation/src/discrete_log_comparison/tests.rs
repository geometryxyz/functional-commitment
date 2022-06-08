#[cfg(test)]
mod tests {
    use ark_bn254::Fr;
    use ark_ff::{Field, SquareRootField};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::rand::thread_rng;
    use ark_std::UniformRand;
    use rand::Rng;

    type F = Fr;

    #[test]
    fn test_square_root_friendly() {
        let m = 6;
        let n = 4;

        let domain = GeneralEvaluationDomain::new(m).unwrap();
        let domain2 = GeneralEvaluationDomain::new(n).unwrap();

        let generator: F = domain.element(1);
        let g2: F = domain2.element(1);

        assert_eq!(g2.pow(&[n as u64]), F::from(1u64));

        // assert_eq!(generator, g2);

        let sq_root = generator.sqrt().unwrap();
        let sq_root2 = g2.sqrt().unwrap();

        assert_eq!(sq_root2.pow(&[2 * n as u64]), F::from(1u64));

        println!("{}", sq_root);

        assert_eq!(generator, sq_root * sq_root);
        assert_eq!(g2, sq_root2 * sq_root2);
    }
}
