#[cfg(test)]
mod tests {
    use crate::{example_circuits::sample_circuit_2, printmatrix, vanilla_ac2tft, R1CSfIndex};
    use ark_bn254::Fr;

    type F = Fr;

    #[test]
    fn manual_check_simple_circuit() {
        let circuit = sample_circuit_2();

        let index: R1CSfIndex<F> = vanilla_ac2tft(circuit);

        index.iter_matrices().for_each(|matrix| {
            printmatrix!(matrix);
            println!("")
        })
    }
}
