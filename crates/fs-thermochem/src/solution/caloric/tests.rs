use super::*;

fn source(species: &str, cp: f64) -> PureLiquidHeatCapacityV1 {
    PureLiquidHeatCapacityV1::new(
        SpeciesId::new(species).unwrap(),
        Temperature::new(100.0),
        vec![
            (Temperature::new(100.0), MolarThermalQuantityV1::new(cp)),
            (Temperature::new(1000.0), MolarThermalQuantityV1::new(cp)),
        ],
    )
    .unwrap()
}
fn model(a: f64, b: f64, cp: [f64; 2]) -> BinaryLiquidCaloricModelV1 {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    BinaryLiquidCaloricModelV1::new(pair, [source("A", cp[0]), source("B", cp[1])]).unwrap()
}
fn controls() -> BinaryCaloricSolveV1 {
    BinaryCaloricSolveV1 {
        temperature_bracket: [Temperature::new(200.0), Temperature::new(800.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
        coexistence: BinaryCoexistenceOptionsV1::default(),
    }
}
fn near(x: f64, y: f64, atol: f64) {
    assert!(
        (x - y).abs() <= atol,
        "{x:.17e} != {y:.17e}, tolerance {atol:e}"
    );
}
fn pure(model: &BinaryLiquidCaloricModelV1, t: f64, n: [f64; 2]) -> BinaryCaloricStateV1 {
    model
        .equilibrate(Temperature::new(t), n, controls().coexistence, || false)
        .unwrap()
}

#[test]
fn g1_variable_cp_integrates_each_source_segment_and_preserves_tiny_increments() {
    let cp = PureLiquidHeatCapacityV1::new(
        SpeciesId::new("A").unwrap(),
        Temperature::new(100.0),
        [(100.0, 10.0), (300.0, 30.0), (500.0, 20.0)]
            .into_iter()
            .map(|(t, c)| (Temperature::new(t), MolarThermalQuantityV1::new(c)))
            .collect(),
    )
    .unwrap();
    let forward = cp
        .enthalpy_change_j_per_mol(Temperature::new(200.0), Temperature::new(400.0))
        .unwrap();
    near(forward, 100.0 * 25.0 + 100.0 * 27.5, 1e-12);
    assert_eq!(
        -forward,
        cp.enthalpy_change_j_per_mol(Temperature::new(400.0), Temperature::new(200.0))
            .unwrap()
    );
    let t0 = 350.0_f64;
    let t1 = t0.next_up();
    let tiny = cp
        .enthalpy_change_j_per_mol(Temperature::new(t0), Temperature::new(t1))
        .unwrap();
    near(tiny, 27.5 * (t1 - t0), 1e-26);
    assert!(
        cp.enthalpy_change_j_per_mol(Temperature::new(99.0), Temperature::new(200.0))
            .is_err()
    );

    // Positive capacities may be subnormal; their crossed-segment integral
    // remains representable over a sufficiently large temperature interval.
    let smallest = f64::from_bits(1);
    let cp = source("A", smallest);
    assert_eq!(
        cp.enthalpy_change_j_per_mol(Temperature::new(300.0), Temperature::new(600.0))
            .unwrap(),
        smallest * 300.0,
    );
    assert!(matches!(
        cp.enthalpy_change_j_per_mol(
            Temperature::new(300.0),
            Temperature::new(300.0_f64.next_up()),
        ),
        Err(BinaryCaloricErrorV1::Arithmetic(_)),
    ));
    let model = model(0.0, 0.0, [smallest; 2]);
    let initial = pure(&model, 500.0, [1.0, 0.0]);
    let unchanged = model.add_heat(&initial, 0.0, controls(), || false).unwrap();
    assert_eq!(unchanged.final_state().temperature(), initial.temperature());
}

#[test]
fn g1_ideal_unequal_temperature_mixing_uses_species_heat_capacities() {
    let model = model(0.0, 0.0, [100.0, 200.0]);
    let left = pure(&model, 300.0, [1.0, 0.0]);
    let right = pure(&model, 500.0, [0.0, 3.0]);
    let result = model.mix(&left, &right, 0.0, controls(), || false).unwrap();
    near(
        result.final_state().temperature().value(),
        (100.0 * 300.0 + 600.0 * 500.0) / 700.0,
        1e-9,
    );
    assert_eq!(result.final_state().component_moles(), [1.0, 3.0]);
    assert_eq!(
        result.final_state().equilibrium().unwrap().phases().len(),
        1
    );
    near(result.energy_residual_j(), 0.0, 1e-7);
    let reverse = model.mix(&right, &left, 0.0, controls(), || false).unwrap();
    assert_eq!(result, reverse);
}

#[test]
fn g1_nonideal_heat_of_mixing_changes_temperature_with_the_correct_sign() {
    for b in [100.0, -100.0] {
        let model = model(0.0, b, [100.0; 2]);
        let left = pure(&model, 300.0, [1.0, 0.0]);
        let right = pure(&model, 300.0, [0.0, 1.0]);
        let result = model.mix(&left, &right, 0.0, controls(), || false).unwrap();
        near(
            result.final_state().temperature().value(),
            300.0 - R * b / (4.0 * 100.0),
            1e-9,
        );
        near(
            result.final_state().enthalpy_j(),
            left.enthalpy_j() + right.enthalpy_j(),
            1e-7,
        );
    }
}

#[test]
fn g1_coexistence_enthalpy_matches_independent_high_precision_thermal_root() {
    let model = model(0.0, 1000.0, [100.0; 2]);
    let left = pure(&model, 400.0, [1.0, 0.0]);
    let right = pure(&model, 400.0, [0.0, 1.0]);
    let result = model.mix(&left, &right, 0.0, controls(), || false).unwrap();
    // Independent 70-digit decimal bisection of ln(a/(1-a)) + (1000/T)(1-2a)=0
    // and 100(T-400)+R*1000*a(1-a)=0. No production solver supplied these values.
    near(
        result.final_state().temperature().value(),
        390.497_105_202_694_2,
        1e-9,
    );
    let equilibrium = result.final_state().equilibrium().unwrap();
    assert_eq!(equilibrium.phases().len(), 2);
    near(
        equilibrium.phases()[0].mole_fraction_a(),
        0.131_616_441_327_145_56,
        2e-12,
    );
    near(
        equilibrium.mixing_enthalpy_j_per_mol(),
        950.289_479_730_579,
        2e-8,
    );
    near(result.energy_residual_j(), 0.0, 1e-7);
}

#[test]
fn g3_unstable_homogeneous_parcel_demixes_and_releases_excess_enthalpy() {
    let model = model(0.0, 900.0, [100.0; 2]);
    let initial = model
        .homogeneous(Temperature::new(300.0), [2.0, 1.0])
        .unwrap();
    assert!(initial.equilibrium().is_none());
    let result = model.add_heat(&initial, 0.0, controls(), || false).unwrap();
    assert!(result.final_state().temperature().value() > 300.0);
    assert_eq!(
        result.final_state().equilibrium().unwrap().phases().len(),
        2
    );
    near(
        result.final_state().enthalpy_j(),
        initial.enthalpy_j(),
        1e-7,
    );
    assert_eq!(
        result.final_state().component_moles(),
        initial.component_moles()
    );
    assert!(initial.equilibrium().is_none());
}

#[test]
fn g3_heating_and_cooling_cross_the_miscibility_gap_without_losing_inventory() {
    let model = model(0.0, 900.0, [100.0; 2]);
    let cold = pure(&model, 300.0, [1.0, 1.0]);
    // chi=3 endpoint computed independently from its scalar common tangent.
    let a = 0.070_720_181_679_944_8;
    let heat = 2.0 * (100.0 * 300.0 + R * 900.0 * (0.25 - a * (1.0 - a)));
    let hot = model.add_heat(&cold, heat, controls(), || false).unwrap();
    near(hot.final_state().temperature().value(), 600.0, 1e-9);
    assert_eq!(hot.final_state().equilibrium().unwrap().phases().len(), 1);
    let recovered = model
        .add_heat(hot.final_state(), -heat, controls(), || false)
        .unwrap();
    near(recovered.final_state().temperature().value(), 300.0, 2e-9);
    assert_eq!(
        recovered
            .final_state()
            .equilibrium()
            .unwrap()
            .phases()
            .len(),
        2
    );
    assert_eq!(
        recovered.final_state().component_moles(),
        cold.component_moles()
    );
}

#[test]
fn g1_equilibrium_heat_capacity_includes_moving_coexistence_compositions() {
    let model = model(0.0, 900.0, [100.0; 2]);
    let t = 350.0;
    let step = 1e-3;
    let before = pure(&model, t - step, [0.3, 0.7]);
    let after = pure(&model, t + step, [0.3, 0.7]);
    let midpoint = pure(&model, t, [0.3, 0.7]);
    let a = midpoint.equilibrium().unwrap().phases()[0].mole_fraction_a();
    let curvature = 1.0 / (a * (1.0 - a)) - 2.0 * 900.0 / t;
    let expected = 100.0 + R * 900.0 * 900.0 * (1.0 - 2.0 * a).powi(2) / (t * t * curvature);
    near(
        (after.enthalpy_j() - before.enthalpy_j()) / (2.0 * step),
        expected,
        2e-6,
    );
    assert!(expected > 100.0);
}

#[test]
fn g0_heat_balance_retains_a_small_input_against_a_large_reference_energy() {
    let model = model(0.0, 0.0, [100.0; 2]);
    let initial = pure(&model, 700.0, [1e15, 1e15]);
    let mut policy = controls();
    policy.temperature_bracket = [Temperature::new(699.0), Temperature::new(701.0)];
    policy.temperature_tolerance_k = 1e-12;
    policy.energy_tolerance_j = 1.0;
    // One representable temperature step has resolvable heat, even when the
    // corresponding addition to the extensive reference energy rounds poorly.
    let target = 700.0_f64.next_up();
    let heat = (target - 700.0) * 2e17;
    let result = model.add_heat(&initial, heat, policy, || false).unwrap();
    assert_eq!(result.final_state().temperature().value(), target);
    near(result.energy_residual_j(), 0.0, 1.0);
    let increment = model
        .enthalpy_change_j(&initial, result.final_state())
        .unwrap();
    near(increment, heat, 1.0);
    assert_eq!(
        model
            .enthalpy_change_j(result.final_state(), &initial)
            .unwrap(),
        -increment
    );
    let other_inventory = pure(&model, 700.0, [2e15, 1e15]);
    assert!(model.enthalpy_change_j(&initial, &other_inventory).is_err());
}

#[test]
fn g0_source_ownership_domains_budgets_and_cancellation_are_atomic() {
    let model = model(0.0, 900.0, [100.0; 2]);
    let initial = pure(&model, 300.0, [1.0, 1.0]);
    let snapshot = initial.clone();
    let other = super::tests::model(0.0, 900.0, [101.0, 100.0]);
    assert!(other.add_heat(&initial, 0.0, controls(), || false).is_err());
    assert!(matches!(
        model.add_heat(&initial, 1e8, controls(), || false),
        Err(BinaryCaloricErrorV1::HeatOutsideBracket { .. })
    ));
    let mut short = controls();
    short.maximum_iterations = 1;
    assert!(matches!(
        model.add_heat(&initial, 123.0, short, || false),
        Err(BinaryCaloricErrorV1::IterationLimit)
    ));
    let mut calls = 0;
    assert!(
        model
            .add_heat(&initial, 100.0, controls(), || {
                calls += 1;
                calls > 7
            })
            .is_err()
    );
    assert_eq!(initial, snapshot);
    let mut outside = controls();
    outside.temperature_bracket[0] = Temperature::new(99.0);
    assert!(model.add_heat(&initial, 0.0, outside, || false).is_err());
    assert!(
        BinaryLiquidCaloricModelV1::new(
            model.thermodynamics().clone(),
            [source("B", 100.0), source("A", 100.0)]
        )
        .is_err()
    );
}

#[test]
fn g0_equilibrium_caloric_totals_obey_the_same_gibbs_identity() {
    let model = model(0.7, 900.0, [100.0; 2]);
    for t in [250.0, 400.0, 800.0] {
        let state = pure(&model, t, [3.0, 7.0]);
        let equilibrium = state.equilibrium().unwrap();
        near(
            equilibrium.mixing_gibbs_j_per_mol(),
            equilibrium.mixing_enthalpy_j_per_mol() - t * equilibrium.mixing_entropy_j_per_mol_k(),
            2e-12,
        );
    }
}
