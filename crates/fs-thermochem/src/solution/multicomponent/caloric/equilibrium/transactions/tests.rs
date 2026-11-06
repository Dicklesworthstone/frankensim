use super::*;

#[test]
fn g1_adiabatic_homogeneous_demixing_matches_binary_calorics_without_a_heat_sink() {
    for datum in [100.0, 1e20] {
        let model = model(2, 1.0, 800.0, datum, false);
        let amounts = [3.5, 6.5];
        let initial = model
            .homogeneous(Temperature::new(400.0), &amounts, || false)
            .unwrap();
        let retained = initial.clone();
        let branch = state(&model, 400.0, &amounts, &seed(2, 0.95));
        let control = controls(350.0, 450.0);
        let result = model
            .relax_homogeneous(&initial, &branch, 0.0, control, || false)
            .unwrap();
        let binary = BinaryLiquidCaloricModelV1::new(
            BinaryRegularSolutionV1::new(
                [
                    model.thermodynamics().species()[0].clone(),
                    model.thermodynamics().species()[1].clone(),
                ],
                1.0,
                800.0,
                Temperature::new(100.0),
                Temperature::new(1000.0),
            )
            .unwrap(),
            [
                model.heat_capacities()[0].clone(),
                model.heat_capacities()[1].clone(),
            ],
        )
        .unwrap();
        let feed = binary.homogeneous(initial.temperature(), amounts).unwrap();
        let oracle = binary
            .add_heat(
                &feed,
                0.0,
                BinaryCaloricSolveV1 {
                    temperature_bracket: control.thermal.temperature_bracket,
                    temperature_tolerance_k: control.thermal.temperature_tolerance_k,
                    energy_tolerance_j: control.thermal.energy_tolerance_j,
                    maximum_iterations: control.thermal.maximum_iterations,
                    coexistence: BinaryCoexistenceOptionsV1::default(),
                },
                || false,
            )
            .unwrap();
        near(
            result.final_state().temperature().value(),
            oracle.final_state().temperature().value(),
            2e-8,
            0.0,
        );
        near(
            result.entropy_change_j_per_k(),
            oracle.entropy_change_j_per_k(),
            3e-8,
            2e-10,
        );
        near(
            result.final_state().isobaric_heat_capacity_j_per_k(),
            oracle
                .final_state()
                .isobaric_heat_capacity_j_per_k()
                .unwrap(),
            1e-7,
            2e-10,
        );
        assert!(result.final_state().temperature().value() > initial.temperature().value());
        assert!(result.entropy_change_j_per_k() > 0.0);
        assert!(result.residual_j().abs() <= control.thermal.energy_tolerance_j);
        assert_eq!(result.final_state().candidate().phases().len(), 2);
        assert!(result.final_state().redistribution_heat_capacity_j_per_k() > 0.0);
        conservation(result.final_state(), &amounts);
        if datum == 100.0 {
            near(
                result.final_state().enthalpy_j(),
                initial.enthalpy_j(),
                2e-6,
                2e-14,
            );
        }
        assert_eq!(initial, retained);
        assert_eq!(branch.temperature().value(), 400.0);
    }
}

#[test]
fn g1_unequal_temperature_mixing_preserves_every_initial_and_final_phase() {
    let model = model(3, 0.0, 1600.0, 100.0, false);
    let left = state(&model, 350.0, &[2.0, 3.0, 5.0], &seed(3, 0.95));
    let right = model
        .homogeneous(Temperature::new(500.0), &[1.0; 3], || false)
        .unwrap();
    let branch = state(&model, 400.0, &[3.0, 4.0, 6.0], &seed(3, 0.95));
    let control = controls(330.0, 500.0);
    let result = model
        .mix_equilibrium(&left, &right, &branch, 0.0, control, || false)
        .unwrap();
    let initial_mixing_h = 10.0 * ternary_mixing(350.0)[0] + R * 1600.0;
    let mut interval: [f64; 2] = [350.0, 500.0];
    for _ in 0..100 {
        let t = 0.5 * interval[0] + 0.5 * interval[1];
        let balance = 1320.0 * (t - 350.0) + 360.0 * (t - 500.0) + 13.0 * ternary_mixing(t)[0]
            - initial_mixing_h;
        if balance > 0.0 {
            interval[1] = t;
        } else {
            interval[0] = t;
        }
    }
    let expected = 0.5 * interval[0] + 0.5 * interval[1];
    let entropy = 1320.0 * (expected / 350.0).ln()
        + 360.0 * (expected / 500.0).ln()
        + 13.0 * ternary_mixing(expected)[1]
        - 10.0 * ternary_mixing(350.0)[1]
        - 3.0 * R * 3.0_f64.ln();
    near(
        result.final_state().temperature().value(),
        expected,
        2e-8,
        0.0,
    );
    near(result.entropy_change_j_per_k(), entropy, 3e-8, 2e-10);
    assert!(entropy > 0.0);
    near(
        result.final_state().enthalpy_j(),
        left.enthalpy_j() + right.enthalpy_j(),
        2e-6,
        2e-14,
    );
    conservation(result.final_state(), &[3.0, 4.0, 6.0]);
    let reverse = model
        .mix_equilibrium(&right, &left, &branch, 0.0, control, || false)
        .unwrap();
    near(
        reverse.final_state().temperature().value(),
        expected,
        2e-8,
        0.0,
    );
    near(reverse.entropy_change_j_per_k(), entropy, 3e-8, 2e-10);
    // Both inputs can already contain several phases; no homogeneous recast of
    // either finite initial inventory is used to determine their stored heat.
    let separated_right = state(&model, 450.0, &[1.0; 3], &seed(3, 0.95));
    let separated = model
        .mix_equilibrium(&left, &separated_right, &branch, 0.0, control, || false)
        .unwrap();
    near(
        separated.final_state().enthalpy_j(),
        left.enthalpy_j() + separated_right.enthalpy_j(),
        2e-6,
        2e-14,
    );
    assert!(separated.entropy_change_j_per_k() > 0.0);
    conservation(separated.final_state(), &[3.0, 4.0, 6.0]);
}

fn explicit_mixing(model: &MulticomponentLiquidCaloricModelV1, x: &[f64]) -> [f64; 2] {
    let n = x.len();
    let mut h = 0.0;
    let mut s = -R
        * x.iter()
            .filter(|&&x| x > 0.0)
            .map(|x| x * x.ln())
            .sum::<f64>();
    for i in 0..n {
        for j in i + 1..n {
            h += R * model.thermodynamics().chi_inverse_temperature_k()[i * n + j] * x[i] * x[j];
            s -= R * model.thermodynamics().chi_constants()[i * n + j] * x[i] * x[j];
        }
    }
    [h, s]
}

#[test]
fn g1_asymmetric_solvent_partition_uses_actual_phase_weights_and_signed_heat() {
    // Same unequal solvent affinities as the existing two-phase fixture at
    // 400 K, now with authored nonzero B_ij and distinct per-phase enthalpies.
    let reference = model(3, 0.0, 0.0, 100.0, false);
    let potential = MulticomponentRegularSolutionV1::new(
        reference.thermodynamics().species().to_vec(),
        vec![0.0, 1.0, 0.25, 1.0, 0.0, 0.325, 0.25, 0.325, 0.0],
        vec![0.0, 1200.0, 100.0, 1200.0, 0.0, 350.0, 100.0, 350.0, 0.0],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let model =
        MulticomponentLiquidCaloricModelV1::new(potential, reference.heat_capacities().to_vec())
            .unwrap();
    let guess = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.75, 0.05, 0.2], vec![0.05, 0.75, 0.2]],
        phase_weights: vec![0.55, 0.45],
    };
    let amounts = [4.5, 3.5, 2.0];
    let initial = model
        .homogeneous(Temperature::new(400.0), &amounts, || false)
        .unwrap();
    let branch = state(&model, 400.0, &amounts, &guess);
    let target = state(&model, 410.0, &amounts, &guess);
    let phase_mixing: Vec<_> = target
        .candidate()
        .phases()
        .iter()
        .map(|phase| explicit_mixing(&model, phase.thermodynamics().mole_fractions()))
        .collect();
    assert!((phase_mixing[0][0] - phase_mixing[1][0]).abs() > 10.0);
    let total = |index| {
        target
            .candidate()
            .phases()
            .iter()
            .zip(&phase_mixing)
            .map(|(phase, mixing)| phase.total_moles() * mixing[index])
            .sum::<f64>()
    };
    let weighted_h = total(0);
    let wrong_equal_weights = 5.0 * (phase_mixing[0][0] + phase_mixing[1][0]);
    assert!((weighted_h - wrong_equal_weights).abs() > 10.0);
    let initial_mixing = explicit_mixing(&model, &[0.45, 0.35, 0.2]);
    let heat = 1100.0 * 10.0 + weighted_h - 10.0 * initial_mixing[0];
    let entropy = 1100.0 * (410.0_f64 / 400.0).ln() + total(1) - 10.0 * initial_mixing[1];
    let result = model
        .relax_homogeneous(&initial, &branch, heat, controls(380.0, 425.0), || false)
        .unwrap();
    near(result.final_state().temperature().value(), 410.0, 2e-8, 0.0);
    near(result.enthalpy_change_j(), heat, 1e-6, 0.0);
    near(result.entropy_change_j_per_k(), entropy, 3e-8, 2e-10);
    near(
        target.enthalpy_j(),
        1100.0 * 310.0 + weighted_h,
        1e-7,
        2e-14,
    );
    conservation(result.final_state(), &amounts);
}

#[test]
fn g1_isolated_mixing_cannot_force_a_higher_gibbs_local_branch() {
    let reference = model(3, 0.5, 1000.0, 100.0, false);
    let capacities = reference
        .thermodynamics()
        .species()
        .iter()
        .map(|id| {
            PureLiquidHeatCapacityV1::new(
                id.clone(),
                Temperature::new(100.0),
                vec![
                    (Temperature::new(100.0), MolarThermalQuantityV1::new(1000.0)),
                    (
                        Temperature::new(1000.0),
                        MolarThermalQuantityV1::new(1000.0),
                    ),
                ],
            )
            .unwrap()
        })
        .collect();
    let model =
        MulticomponentLiquidCaloricModelV1::new(reference.thermodynamics().clone(), capacities)
            .unwrap();
    let left = state(&model, 420.0, &[0.5; 3], &seed(3, 0.9));
    let right = left.clone();
    let homogeneous_branch = state(
        &model,
        420.0,
        &[1.0; 3],
        &MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![vec![1.0 / 3.0; 3]],
            phase_weights: vec![1.0],
        },
    );
    assert!(matches!(
        model.mix_equilibrium(
            &left,
            &right,
            &homogeneous_branch,
            0.0,
            controls(410.0, 425.0),
            || false
        ),
        Err(MulticomponentEquilibriumCaloricErrorV1::Caloric(
            MulticomponentCaloricErrorV1::EntropyDecrease { .. }
        ))
    ));
    assert_eq!(left, right);
}

#[test]
fn g3_input_inventory_source_branch_work_and_cancellation_refuse_atomically() {
    let model = model(3, 0.0, 1600.0, 100.0, false);
    let left = state(&model, 350.0, &[2.0, 3.0, 5.0], &seed(3, 0.95));
    let right = model
        .homogeneous(Temperature::new(500.0), &[1.0; 3], || false)
        .unwrap();
    let branch = state(&model, 400.0, &[3.0, 4.0, 6.0], &seed(3, 0.95));
    let retained = (left.clone(), right.clone(), branch.clone());
    let control = controls(330.0, 500.0);
    assert!(matches!(
        model.mix_equilibrium(&left, &right, &left, 0.0, control, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
    ));
    let foreign = super::model(3, 0.0, 1600.0, 350.0, false);
    assert!(matches!(
        foreign.mix_equilibrium(&left, &right, &branch, 0.0, control, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
    ));
    let exhausted = MulticomponentEquilibriumCaloricSolveV1 {
        maximum_phase_solves: 1,
        ..control
    };
    assert!(matches!(
        model.mix_equilibrium(&left, &right, &branch, 0.0, exhausted, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::PhaseSolveBudget { .. })
    ));
    assert!(
        model
            .mix_equilibrium(&left, &right, &branch, 1e10, control, || false)
            .is_err()
    );
    let mut outside = control;
    outside.thermal.temperature_bracket[1] = Temperature::new(1100.0);
    assert!(matches!(
        model.mix_equilibrium(&left, &right, &branch, 0.0, outside, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
    ));
    let mut count = 0;
    model
        .mix_equilibrium(&left, &right, &branch, 0.0, control, || {
            count += 1;
            false
        })
        .unwrap();
    for cut in [1, 250, count] {
        let mut polls = 0;
        let result = model.mix_equilibrium(&left, &right, &branch, 0.0, control, || {
            polls += 1;
            polls == cut
        });
        assert!(
            matches!(
                result,
                Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
            ),
            "{result:?}"
        );
        assert_eq!(polls, cut);
    }
    assert_eq!((left, right, branch), retained);
}
