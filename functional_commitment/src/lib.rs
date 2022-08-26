#[cfg(test)]
mod tests {
    use ac2tft::{sample_matrices, SparseMatrices, printmatrix};
    use ark_poly_commit::PolynomialCommitment;
    use ark_std::rand::thread_rng;
    use ark_ff::{PrimeField, to_bytes};
    use ark_bn254::{Fr, Bn254};
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use index_private_marlin::AHPForR1CS;
    use proof_of_function_relation::{t_strictly_lower_triangular_test::TStrictlyLowerTriangular};
    use proof_of_function_relation::{t_diag::TDiag};
    use homomorphic_poly_commit::{AdditivelyHomomorphicPCS, marlin_kzg::KZG10};
    use fiat_shamir_rng::{FiatShamirRng, SimpleHashFiatShamirRng};
    use rand_chacha::ChaChaRng;
    use blake2::Blake2s;

    type FS = SimpleHashFiatShamirRng<Blake2s, ChaChaRng>;

    type F = Fr;
    type PC = KZG10<Bn254>;

    #[test]
    fn test_matrix_arithmetization() {
        let rng = &mut thread_rng();        
        let matrices: SparseMatrices<F> = sample_matrices::<F>();

        let number_of_public_variables = 8;

        let index = AHPForR1CS::index_from_functional_triple(
            matrices.0, 
            matrices.1, 
            matrices.2, 
            number_of_public_variables
        ).unwrap();

        let domain_h = GeneralEvaluationDomain::<F>::new(index.index_info.num_constraints).unwrap();
        let domain_k = GeneralEvaluationDomain::<F>::new(index.index_info.num_non_zero).unwrap();

        let enforced_degree_bound = domain_k.size() + 1;
        let enforced_hiding_bound = 1;

        let max_degree = 20;
        let pp = PC::setup(max_degree, None, rng).unwrap();
        let (ck, vk) = PC::trim(
            &pp,
            max_degree,
            enforced_hiding_bound,
            Some(&[2, enforced_degree_bound]),
        )
        .unwrap();

        let (commits, rands) = PC::commit(&ck, index.iter_individual_matrices(), None).unwrap();
        
        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        let row_evals = index.a_arith.evals_on_K.row;
        let col_evals = index.a_arith.evals_on_K.col;
        
        row_evals.evals.iter().zip(col_evals.evals.iter()).enumerate().for_each(|(i, (row, col))| if row == col {println!("equal at gamma_{}", i)});

        let a_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
            &ck,
            number_of_public_variables,
            &domain_k,
            &domain_h, 
            &index.a_arith.col, 
            &commits[1],
            &rands[1], 
            &index.a_arith.row, 
            &commits[0], 
            &rands[0], 
            None, 
            &mut fs_rng,
            rng
        ).unwrap();

        let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

        assert_eq!(
            true,
            TStrictlyLowerTriangular::verify(
                &vk,
                &ck,
                number_of_public_variables,
                &domain_k,
                &domain_h,
                &commits[1],
                &commits[0],
                None,
                a_proof,
                &mut fs_rng
            ).is_ok());


            let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());

            let b_proof = TStrictlyLowerTriangular::<F, PC, FS>::prove(
                &ck,
                number_of_public_variables,
                &domain_k,
                &domain_h, 
                &index.b_arith.col, 
                &commits[4],
                &rands[4], 
                &index.b_arith.row, 
                &commits[3], 
                &rands[3], 
                None, 
                &mut fs_rng,
                rng
            ).unwrap();
    
            let mut fs_rng = FS::initialize(&to_bytes!(b"Testing :)").unwrap());
            assert_eq!(
                true,
                TStrictlyLowerTriangular::verify(
                    &vk,
                    &ck,
                    number_of_public_variables,
                    &domain_k,
                    &domain_h,
                    &commits[4],
                    &commits[3],
                    None,
                    b_proof,
                    &mut fs_rng
                ).is_ok());


            // let c_proof = TDiag::<F, PC, FS>::prove(
            //     &ck,
            //     number_of_public_variables,
            //     &index.c_arith.row,
            //     &index.c_arith.col,
            //     &index.c_arith.val,
            //     &commits[6],
            //     &commits[7],
            //     &commits[8],
            //     &rands[6],
            //     &rands[7],
            //     &rands[8],
            //     None,
            //     &domain_k,
            //     &domain_h, 
            //     rng
            // ).unwrap();

            // let is_valid = TDiag::<F, PC, FS>::verify(
            //     &vk,
            //     number_of_public_variables,
            //     &commits[6],
            //     &commits[7],
            //     &commits[8],
            //     None,
            //     &domain_h,
            //     &domain_k,
            //     c_proof,
            // );
    
            // assert!(is_valid.is_ok());

    }

    #[test]
    fn reindex_test() {
        let domain_big = GeneralEvaluationDomain::<F>::new(8).unwrap();
        let domain_small = GeneralEvaluationDomain::<F>::new(1).unwrap();

        for i in 0..8 {
            let reindexed = domain_big.reindex_by_subdomain(domain_small, i);
            println!("for index: {}, reindexed is : {}", i, reindexed);
        }

    }

    #[test]
    fn read_matrices() {
        let matrices: SparseMatrices<F> = sample_matrices::<F>();

        printmatrix!(matrices.1)
    }
}
