// #[cfg(test)]
// pub mod tests {
//     // use holographic_lincheck::matrix::MatrixIndexer;
//     use ark_poly::univariate::DensePolynomial;
//     use ark_poly::{
//         EvaluationDomain,
//         GeneralEvaluationDomain,
//     };
//     use ark_poly_commit::PolynomialCommitment;
//     use ark_bn254::{Bn254, Fr};
//     use ark_ff::{Field, to_bytes};
//     use ark_std::rand::thread_rng;
//     use ark_marlin::rng::FiatShamirRng;
//     use ark_relations::{
//         lc,
//         r1cs::{
//             Variable,
//             ConstraintSystem,
//             OptimizationGoal,
//             ConstraintSynthesizer,
//             ConstraintSystemRef,
//             ConstraintMatrices,
//             SynthesisError,
//             Matrix,
//         },
//     };
//     //use ark_marlin::constraint_systems::make_matrices_square_for_indexer;
//     use blake2::Blake2s;
//     use crate::{SparseMatrices, GateType, GateInput, Gate, gates_to_sparse_matrices};
//     use crate::error::Error;
//     use proof_of_function_relation::commitment::KZG10;
//     use proof_of_function_relation::t_strictly_lower_triangular_test::TStrictlyLowerTriangular;
//     use proof_of_function_relation::t_strictly_lower_triangular_test::proof::Proof as TSltProof;
//     use proof_of_function_relation::t_diag::TDiag;
//     use proof_of_function_relation::t_diag::proof::Proof as TDiagProof;

//     type F = Fr;
//     type PC = KZG10<Bn254>;
//     type D = Blake2s;

//     #[test]
//     fn test_gate_formatting() {
//         // g0: (1, 1, +)
//         let gate = Gate {
//             left: GateInput::Constant(F::from(1u64)),
//             right: GateInput::Constant(F::from(1u64)),
//             symbol: GateType::Add,
//             label: String::from("g0"),
//         };

//         assert_eq!(
//             format!("{}", gate),
//             "g0: (Fp256 \"(0000000000000000000000000000000000000000000000000000000000000001)\", \
//             Fp256 \"(0000000000000000000000000000000000000000000000000000000000000001)\", +)"
//         );
//     }

//     pub fn sample_gates_0() -> Vec<Gate<F>> {
//         // Encodes x^3 + 2x + 5
//         // g0: (x, x, *)
//         // g1: (g0, x, *)
//         // g2: (x, 2, *)
//         // g3: (g1, g2, +)
//         // g4: (g3, 5, +)

//         // s: [*, *, *, +, +]
//         // l: [x, g0, g1, g3]
//         // r: [x, 2, g2, 5]

//         let g0 = Gate::<Fr> {
//             left: GateInput::Input(String::from("x")),
//             right: GateInput::Input(String::from("x")),
//             symbol: GateType::Mul,
//             label: String::from("g0"),
//         };

//         let g1 = Gate::<Fr> {
//             left: GateInput::Gate(Box::new(g0.clone())),
//             right: GateInput::Input(String::from("x")),
//             symbol: GateType::Mul,
//             label: String::from("g1"),
//         };

//         let g2 = Gate::<Fr> {
//             left: GateInput::Input(String::from("x")),
//             right: GateInput::Constant(F::from(2u64)),
//             symbol: GateType::Mul,
//             label: String::from("g2"),
//         };

//         let g3 = Gate::<Fr> {
//             left: GateInput::Gate(Box::new(g1.clone())),
//             right: GateInput::Gate(Box::new(g2.clone())),
//             symbol: GateType::Add,
//             label: String::from("g3"),
//         };

//         let g4 = Gate::<Fr> {
//             left: GateInput::Gate(Box::new(g3.clone())),
//             right: GateInput::Constant(F::from(5u64)),
//             symbol: GateType::Add,
//             label: String::from("g4"),
//         };

//         return vec![g0, g1, g2, g3, g4];
//     }

//     #[test]
//     fn test_gate_input_eq() {
//         let gates = sample_gates_0();
//         let g0 = &gates[0];
//         let g1 = &gates[1];
//         let g2 = &gates[2];
//         assert_eq!(g0.left, g0.right);
//         assert_eq!(g0.left, g1.right);
//         assert_eq!(g2.right, GateInput::Constant(F::from(2u64)));
//     }

//     #[test]
//     fn test_into_gate() {
//         let gates = sample_gates_0();
//         let g0 = &gates[0];
//         let g1 = &gates[1];
//         assert_eq!(g1.left.into_gate().unwrap(), *g0);
//         assert_eq!(g1.right.into_gate(), Err(Error::GateInputNotGate));
//     }

//     #[derive(Copy, Clone)]
//     pub struct SampleCircuit0<F: Field> {
//         x: Option<F>,
//         out: Option<F>,
//     }

//     impl<F: Field> ConstraintSynthesizer<F> for SampleCircuit0<F> {
//         fn generate_constraints(
//             self,
//             cs: ConstraintSystemRef<F>,
//         ) -> Result<(), SynthesisError> {
//             // Encodes x^3 + 2x + 5
//             let x = cs.new_input_variable(|| self.x.ok_or(SynthesisError::AssignmentMissing))?;
//             let out = cs.new_input_variable(|| self.out.ok_or(SynthesisError::AssignmentMissing))?;

//             // x^3
//             let x2 = cs.new_witness_variable(|| {
//                 let mut x = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 x.square_in_place();
//                 Ok(x)
//             })?;

//             let x3 = cs.new_witness_variable(|| {
//                 let x = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 let mut x2 = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 let mut x3 = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 x2.mul_assign(&x);
//                 x3.mul_assign(&x2);
//                 Ok(x3)
//             })?;

//             // Enforce that x3 = x2 * x
//             cs.enforce_constraint(lc!() + x, lc!() + x2, lc!() + x3)?;

//             // x_double = 2x
//             let x_double = cs.new_witness_variable(|| {
//                 let x = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 Ok(x * F::from(2u64))
//             })?;

//             // Enforce that x_double = 2 * x
//             cs.enforce_constraint(
//                 lc!() + (F::from(2u64), Variable::One),
//                 lc!() + x,
//                 lc!() + x_double
//             )?;

//             // x_double_plus_5
//             let x_double_plus_5 = cs.new_witness_variable(|| {
//                 let x = self.x.ok_or(SynthesisError::AssignmentMissing)?;
//                 Ok(x * F::from(2u64) + F::from(5u64))
//             })?;

//             // Enforce that x_double_plus_5 = x_double_plus_5_w * 1
//             cs.enforce_constraint(
//                 lc!() + (F::from(1u64), Variable::One),
//                 lc!() + (F::from(1u64), x_double) + (F::from(5u64), Variable::One),
//                 lc!() + x_double_plus_5
//             )?;

//             // out = x^3 + 2x + 5
//             cs.enforce_constraint(
//                 lc!() + (F::from(1u64), Variable::One),
//                 lc!() + (F::from(1u64), x3) + (F::from(1u64), x_double_plus_5),
//                 lc!() + out
//             )?;

//             Ok(())
//         }
//     }

//     #[test]
//     pub fn test_sample_circuit_0() {
//         // Encodes x^3 + 2x + 5
//         let circ = SampleCircuit0 {
//             x: Some(Fr::from(1u64)),
//             out: Some(Fr::from(9u64))
//         };
//         let pcs = ConstraintSystem::new_ref();
//         pcs.set_optimization_goal(OptimizationGoal::Weight);
//         pcs.set_mode(ark_relations::r1cs::SynthesisMode::Prove {
//             construct_matrices: true,
//         });
//         let _ = circ.generate_constraints(pcs.clone());

//         // Requires importing a Marlin fork
//         //make_matrices_square_for_indexer(pcs.clone());

//         let constraint_matrices: ConstraintMatrices<Fr> = pcs.to_matrices().unwrap();

//         //print_matrix("sample_circuit_0 matrix a", constraint_matrices.a);
//         //print_matrix("sample_circuit_0 matrix b", constraint_matrices.b);
//         //print_matrix("sample_circuit_0 matrix c", constraint_matrices.c);

//         //let gates = sample_gates_0();
//         //let matrices = gates_to_sparse_matrices(gates);
//         //print_matrix("sample_gates_0 matrix a", matrices.0);
//         //print_matrix("sample_gates_0 matrix b", matrices.1);
//         //print_matrix("sample_gates_0 matrix c", matrices.2);
//     }

//     // #[test]
//     // pub fn test_sample_gates_0_to_polys() {
//     //     let gates = sample_gates_0();
//     //     let matrices = gates_to_sparse_matrices(gates);

//     //     //print_polynomial("gates_matrix_a rows", &result.row_poly);
//     //     //print_polynomial("gates_matrix_a cols", &result.col_poly);
//     //     //print_polynomial("gates_matrix_a vals", &result.vals_poly);
//     //     //print_polynomial("gates_matrix_a row_col", &result.row_col_poly);

//     //     let num_constraints = matrices.0.len();

//     //     // TODO: this may be wrong, but setting the size of domain_h to num_constraints and
//     //     // domain_k to t causes an error in DLComparision
//     //     let domain_h = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     let domain_k = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     let indexed_a = MatrixIndexer::index(&domain_k, &domain_h, &matrices.0);
//     //     do_t_slt_proof(indexed_a, domain_k, domain_h, calc_top_non_zero(matrices.0));

//     //     let num_constraints = matrices.1.len();
//     //     let domain_h = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     let domain_k = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     let indexed_b = MatrixIndexer::index(&domain_k, &domain_h, &matrices.1);
//     //     do_t_slt_proof(indexed_b, domain_k, domain_h, calc_top_non_zero(matrices.1));

//     //     // TODO: note that t-diag fails
//     //     //let num_constraints = matrices.2.len();
//     //     //let domain_h = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     //let domain_k = GeneralEvaluationDomain::<F>::new(num_constraints).unwrap();
//     //     //let indexed_c = MatrixIndexer::index(&domain_k, &domain_h, &matrices.2);
//     //     //do_t_diag_proof(indexed_c, domain_k, domain_h, calc_top_non_zero(matrices.2));
//     // }

//     pub fn calc_top_non_zero(
//         matrix: Matrix<F>
//     ) -> usize {
//         let mut t = 0;
//         for row in matrix.iter() {
//             if row.len() == 0 {
//                 t += 1;
//             } else {
//                 break
//             }
//         }
//         t
//     }

//     // pub fn do_t_slt_proof(
//     //     indexed: MatrixIndexer<F>,
//     //     domain_k: GeneralEvaluationDomain<F>,
//     //     domain_h: GeneralEvaluationDomain<F>,
//     //     t: usize,
//     // ) -> TSltProof<F, PC> {
//     //     let row_poly = indexed.row_poly;
//     //     let col_poly = indexed.col_poly;

//     //     let mut rng = thread_rng();
//     //     let max_degree = 20;
//     //     let pp = PC::setup(max_degree, None, &mut rng).unwrap();
//     //     let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();
//     //     let (commitments, _) =
//     //         PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(&mut rng)).unwrap();

//     //     let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());

//     //     let proof = TStrictlyLowerTriangular::<F, PC, D>::prove(
//     //         &ck,
//     //         t,
//     //         &domain_k,
//     //         &domain_h,

//     //         &col_poly,
//     //         &row_poly,
//     //         &commitments[1].clone(),
//     //         &commitments[0].clone(),

//     //         //&row_poly,
//     //         //&col_poly,
//     //         //&commitments[0].clone(),
//     //         //&commitments[1].clone(),

//     //         &mut fs_rng,
//     //         &mut rng,
//     //     ).unwrap();

//     //     proof
//     // }

//     // pub fn do_t_diag_proof(
//     //     indexed: MatrixIndexer<F>,
//     //     domain_k: GeneralEvaluationDomain<F>,
//     //     domain_h: GeneralEvaluationDomain<F>,
//     //     t: usize,
//     // ) -> TDiagProof<F, PC> {
//     //     let row_poly = indexed.row_poly;
//     //     let col_poly = indexed.col_poly;
//     //     let val_poly = indexed.vals_poly;

//     //     let mut rng = thread_rng();
//     //     let max_degree = 20;
//     //     let pp = PC::setup(max_degree, None, &mut rng).unwrap();
//     //     let (ck, vk) = PC::trim(&pp, max_degree, 0, None).unwrap();
//     //     let (commitments, rands) =
//     //         PC::commit(&ck, &[row_poly.clone(), col_poly.clone()], Some(&mut rng)).unwrap();

//     //     let mut fs_rng = FiatShamirRng::<D>::from_seed(&to_bytes!(b"Testing :)").unwrap());
//     //     let proof = TDiag::<F, PC, D>::prove(
//     //         &ck,
//     //         t,
//     //         &row_poly,
//     //         &col_poly,
//     //         &val_poly,
//     //         &commitments[0],
//     //         &commitments[1],
//     //         &commitments[2],
//     //         &rands[0],
//     //         &rands[1],
//     //         &rands[2],
//     //         &domain_k,
//     //         &domain_h,
//     //         &mut rng,
//     //     )
//     //     .unwrap();

//     //     proof
//     // }

//     fn print_matrix(
//         label: &str,
//         matrix: Matrix<F>,
//     ) {
//         println!("{}:", label);
//         for (i, i_vec) in matrix.iter().enumerate() {
//             print!("{}: ", i);
//             for j_tuple in i_vec.iter() {
//                 print!("({}, {}), ", j_tuple.0, j_tuple.1);
//             }
//             println!("");
//         }
//         println!("");
//     }

//     fn print_joint_matrix(
//         label: &str,
//         matrix: Vec<Vec<usize>>,
//     ) {
//         println!("{}:", label);
//         println!("{:?}:", matrix);
//         println!("");
//     }

//     fn print_polynomial(
//         label: &str,
//         poly: &DensePolynomial<Fr>,
//     ) {
//         println!("{} coeffs: ", label,);
//         for c in poly.coeffs.iter() {
//             println!("{}", c);
//         }
//         println!("");
//     }
// }
