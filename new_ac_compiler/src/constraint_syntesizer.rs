pub trait ConstraintSynthesizer {}

/*
    Circuit

    impl ConstraintSynthesizer for circuit

        fn build_constraints(
            constraint_builder: &mut ConstraintBuilder
        ) {
            let x = constraint_builder.new_input(x, 7)?;
            let twp = constraint_builder.new_input(two, 2)?;
            let five = constraint_builder.new_input(five, 5)?;

            x_sq = cs.enforce_constraint(x, x, mul)?;
            x_qubed = cs.enforce_constraint(x_sq, x, mul)?;

            two_x = cs.enforce_constraint(x, two, mul)?;

            x_qubed_plus_2x = cs.enforce_constraint(x_qubed, two_x, plus)?;
            output = cs.enforce_constraint(x_qubed_plus_2x, five, plus)?;
        }
*/
