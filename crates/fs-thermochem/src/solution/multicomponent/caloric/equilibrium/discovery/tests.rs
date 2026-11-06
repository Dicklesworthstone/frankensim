use super::super::discovery::MulticomponentDiscoveredCaloricSolveV1;
use super::*;
use crate::solution::caloric::BinaryCaloricErrorV1;
use crate::solution::multicomponent::discovery::{
    MulticomponentPhaseDiscoveryControlV1, MulticomponentPhaseDiscoveryErrorV1,
};

fn control(lo: f64, hi: f64) -> MulticomponentDiscoveredCaloricSolveV1 {
    MulticomponentDiscoveredCaloricSolveV1 {
        thermal: controls(lo, hi).thermal,
        discovery: MulticomponentPhaseDiscoveryControlV1 {
            equilibrium: MulticomponentEquilibriumControlV1 {
                chemical_potential_tolerance: 2e-13,
                ..Default::default()
            },
            ..Default::default()
        },
        response: Default::default(),
        maximum_candidate_solves: 16_384,
    }
}

fn binary(model: &MulticomponentLiquidCaloricModelV1) -> BinaryLiquidCaloricModelV1 {
    BinaryLiquidCaloricModelV1::new(
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
    .unwrap()
}

fn binary_controls(control: MulticomponentDiscoveredCaloricSolveV1) -> BinaryCaloricSolveV1 {
    BinaryCaloricSolveV1 {
        temperature_bracket: control.thermal.temperature_bracket,
        temperature_tolerance_k: control.thermal.temperature_tolerance_k,
        energy_tolerance_j: control.thermal.energy_tolerance_j,
        maximum_iterations: control.thermal.maximum_iterations,
        coexistence: Default::default(),
    }
}

#[test]
fn g1_heating_removes_and_cooling_creates_a_phase_with_original_binary_h_s() {
    for datum in [100.0, 1e20] {
        let model = model(2, 1.0, 800.0, datum, false);
        let amounts = [3.5, 6.5];
        let cold = state(&model, 400.0, &amounts, &seed(2, 0.95));
        let hot = model
            .homogeneous(Temperature::new(900.0), &amounts, || false)
            .unwrap();
        let retained = (cold.clone(), hot.clone());
        let binary = binary(&model);
        let binary_cold = binary
            .equilibrate(cold.temperature(), amounts, Default::default(), || false)
            .unwrap();
        let binary_hot = binary
            .equilibrate(hot.temperature(), amounts, Default::default(), || false)
            .unwrap();
        let heat = binary.enthalpy_change_j(&binary_cold, &binary_hot).unwrap();
        let entropy = binary
            .entropy_change_j_per_k(&binary_cold, &binary_hot)
            .unwrap();
        let control = control(350.0, 980.0);
        let heated = model
            .add_heat_discovered(&cold, heat, control, || false)
            .unwrap();
        let cooled = model
            .add_heat_discovered(&hot, -heat, control, || false)
            .unwrap();
        for (result, t, phases, q, ds) in [
            (&heated, 900.0, 1, heat, entropy),
            (&cooled, 400.0, 2, -heat, -entropy),
        ] {
            near(result.final_state().temperature().value(), t, 2e-8, 0.0);
            assert_eq!(result.final_state().candidate().phases().len(), phases);
            assert_eq!(result.discovery_report().selected_phase_count(), phases);
            near(result.balance().enthalpy_change_j(), q, 1e-6, 0.0);
            near(result.balance().entropy_change_j_per_k(), ds, 3e-8, 2e-10);
            assert!(result.balance().residual_j().abs() <= control.thermal.energy_tolerance_j);
            assert!(result.completed_searches() >= 3);
            assert!(result.balance().completed_phase_solves() >= result.completed_searches());
            assert_eq!(result.controls(), control);
            conservation(result.final_state(), &amounts);
        }
        near(
            heated.final_state().isobaric_heat_capacity_j_per_k(),
            1060.0,
            1e-10,
            0.0,
        );
        assert!(cooled.final_state().redistribution_heat_capacity_j_per_k() > 0.0);
        assert_eq!((cold, hot), retained);
    }
}

#[test]
fn g1_adiabatic_demixing_and_unequal_temperature_phase_mixing_conserve_full_storage() {
    let model = model(2, 1.0, 800.0, 100.0, true);
    let binary = binary(&model);
    let control = control(350.0, 980.0);
    let initial = model
        .homogeneous(Temperature::new(400.0), &[3.5, 6.5], || false)
        .unwrap();
    let binary_initial = binary
        .homogeneous(initial.temperature(), [3.5, 6.5])
        .unwrap();
    let relaxed = model
        .add_heat_discovered(&initial, 0.0, control, || false)
        .unwrap();
    let oracle = binary
        .add_heat(&binary_initial, 0.0, binary_controls(control), || false)
        .unwrap();
    near(
        relaxed.final_state().temperature().value(),
        oracle.final_state().temperature().value(),
        2e-8,
        0.0,
    );
    near(
        relaxed.balance().entropy_change_j_per_k(),
        oracle.entropy_change_j_per_k(),
        3e-8,
        2e-10,
    );
    assert_eq!(relaxed.final_state().candidate().phases().len(), 2);
    assert!(relaxed.balance().entropy_change_j_per_k() > 0.0);
    near(
        relaxed.final_state().enthalpy_j(),
        initial.enthalpy_j(),
        1e-6,
        0.0,
    );

    let left = state(&model, 400.0, &[2.0, 3.0], &seed(2, 0.95));
    let right = model
        .homogeneous(Temperature::new(900.0), &[1.5, 3.5], || false)
        .unwrap();
    let binary_left = binary
        .equilibrate(left.temperature(), [2.0, 3.0], Default::default(), || false)
        .unwrap();
    let binary_right = binary.homogeneous(right.temperature(), [1.5, 3.5]).unwrap();
    let mixed = model
        .mix_discovered(&left, &right, 0.0, control, || false)
        .unwrap();
    let oracle = binary
        .mix(
            &binary_left,
            &binary_right,
            0.0,
            binary_controls(control),
            || false,
        )
        .unwrap();
    near(
        mixed.final_state().temperature().value(),
        oracle.final_state().temperature().value(),
        2e-8,
        0.0,
    );
    near(
        mixed.balance().entropy_change_j_per_k(),
        oracle.entropy_change_j_per_k(),
        3e-8,
        2e-10,
    );
    near(
        mixed.final_state().enthalpy_j(),
        left.enthalpy_j() + right.enthalpy_j(),
        1e-6,
        0.0,
    );
    conservation(mixed.final_state(), &[3.5, 6.5]);
    let reversed = model
        .mix_discovered(&right, &left, 0.0, control, || false)
        .unwrap();
    near(
        mixed.final_state().temperature().value(),
        reversed.final_state().temperature().value(),
        1e-9,
        0.0,
    );
    assert!(mixed.balance().entropy_change_j_per_k() > 0.0);
}

#[test]
fn g1_three_phase_disappearance_retains_exactly_absent_species() {
    let model = model(4, 0.0, 1600.0, 100.0, false);
    let amounts = [2.0, 3.0, 5.0, 0.0];
    let initial = state(
        &model,
        400.0,
        &amounts,
        &MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![
                vec![0.9, 0.05, 0.05, 0.0],
                vec![0.05, 0.9, 0.05, 0.0],
                vec![0.05, 0.05, 0.9, 0.0],
            ],
            phase_weights: vec![1.0; 3],
        },
    );
    let expected = model
        .homogeneous(Temperature::new(700.0), &amounts, || false)
        .unwrap();
    let heat = expected.enthalpy_j() - initial.enthalpy_j();
    let result = model
        .add_heat_discovered(&initial, heat, control(380.0, 750.0), || false)
        .unwrap();
    assert_eq!(initial.candidate().phases().len(), 3);
    assert_eq!(result.final_state().candidate().phases().len(), 1);
    near(result.final_state().temperature().value(), 700.0, 2e-8, 0.0);
    near(result.balance().enthalpy_change_j(), heat, 1e-6, 0.0);
    conservation(result.final_state(), &amounts);
    assert_eq!(
        result.final_state().candidate().phases()[0].component_moles()[3],
        0.0
    );
}

#[test]
fn g3_admitted_metastable_homogeneous_search_cannot_decrease_adiabatic_entropy() {
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
    let initial = state(&model, 420.0, &[1.0; 3], &seed(3, 0.9));
    let mut control = control(410.0, 425.0);
    control.discovery.maximum_phases = 1;
    control.discovery.maximum_seed_trials = 1;
    assert!(matches!(
        model.add_heat_discovered(&initial, 0.0, control, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Caloric(
            MulticomponentCaloricErrorV1::EntropyDecrease { .. }
        ))
    ));
}

#[test]
fn g3_first_order_enthalpy_gap_is_refused_without_a_fabricated_latent_plateau() {
    // Symmetric ternary regular solution: the r=2/3, q=1/6 three-phase
    // state and uniform state share a tangent at chi=4 ln(2). Their mixing
    // enthalpies are n R B/4 and n R B/3. Intermediate heat requires an
    // additional coexistence owner, not interpolation across a search switch.
    let model = model(3, 0.5, 1000.0, 100.0, false);
    let initial = state(&model, 420.0, &[1.0; 3], &seed(3, 0.9));
    let transition = 1000.0 / (4.0 * 2.0_f64.ln() - 0.5);
    let target_enthalpy =
        360.0 * (transition - 100.0) + 3.0 * R * 1000.0 * (0.25 + 1.0 / 3.0) / 2.0;
    let result = model.add_heat_discovered(
        &initial,
        target_enthalpy - initial.enthalpy_j(),
        control(transition - 4.0, transition + 4.0),
        || false,
    );
    assert!(
        matches!(
            result,
            Err(MulticomponentEquilibriumCaloricErrorV1::Caloric(
                MulticomponentCaloricErrorV1::Caloric(BinaryCaloricErrorV1::IterationLimit)
            ))
        ),
        "{result:?}"
    );
}

#[test]
fn g3_search_budget_unstable_feed_source_and_cancellation_refuse_atomically() {
    let model = model(2, 1.0, 800.0, 100.0, false);
    let initial = model
        .homogeneous(Temperature::new(400.0), &[3.5, 6.5], || false)
        .unwrap();
    let retained = initial.clone();
    let mut control = control(350.0, 980.0);
    control.maximum_candidate_solves = control.discovery.maximum_seed_trials;
    assert!(matches!(
        model.add_heat_discovered(&initial, 0.0, control, || false),
        Err(
            MulticomponentEquilibriumCaloricErrorV1::DiscoverySolveBudget {
                completed_candidate_solves: 5,
                required_trial_allowance: 64,
            }
        )
    ));
    control.maximum_candidate_solves = 16_384;
    control.discovery.maximum_seed_trials = 1;
    assert!(matches!(
        model.add_heat_discovered(&initial, 0.0, control, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Discovery(
            MulticomponentPhaseDiscoveryErrorV1::NoAdmissibleCandidate {
                attempted_candidates: 1
            }
        ))
    ));
    control.discovery.maximum_seed_trials = 64;
    let foreign = super::model(2, 1.0, 800.0, 350.0, false);
    assert!(
        model
            .add_heat_discovered(
                &foreign
                    .homogeneous(initial.temperature(), initial.component_moles(), || false)
                    .unwrap(),
                0.0,
                control,
                || false,
            )
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        model.add_heat_discovered(&initial, 0.0, control, || {
            polls += 1;
            polls > 100
        }),
        Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
    );
    assert!(polls > 100);
    assert_eq!(initial, retained);
}
