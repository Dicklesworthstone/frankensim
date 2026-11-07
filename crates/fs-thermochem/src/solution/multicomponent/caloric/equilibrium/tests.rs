use super::*;
use crate::solution::caloric::{BinaryCaloricSolveV1, BinaryLiquidCaloricModelV1};
use crate::solution::multicomponent::MulticomponentRegularSolutionV1;
use crate::solution::multicomponent::equilibrium::MulticomponentEquilibriumControlV1;
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use crate::{MolarThermalQuantityV1, SpeciesId};

#[path = "coexistence/tests.rs"]
mod coexistence;
#[path = "discovery/tests.rs"]
mod discovery;
#[path = "transactions/tests.rs"]
mod transactions;

fn potential(n: usize, a: f64, b: f64) -> MulticomponentRegularSolutionV1 {
    let matrix = |coefficient| {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| if i == j { 0.0 } else { coefficient }))
            .collect()
    };
    MulticomponentRegularSolutionV1::new(
        (0..n)
            .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
            .collect(),
        matrix(a),
        matrix(b),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}

fn model(
    n: usize,
    a: f64,
    b: f64,
    datum: f64,
    variable: bool,
) -> MulticomponentLiquidCaloricModelV1 {
    let law = potential(n, a, b);
    let sources = law
        .species()
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let cp = 80.0 + 40.0 * i as f64;
            let knots = if variable {
                vec![
                    (100.0, 0.5 * cp),
                    (350.0, cp),
                    (450.0, 1.5 * cp),
                    (1000.0, 2.0 * cp),
                ]
            } else {
                vec![(100.0, cp), (1e20, cp)]
            };
            PureLiquidHeatCapacityV1::new(
                id.clone(),
                Temperature::new(datum),
                knots
                    .into_iter()
                    .map(|(t, value)| (Temperature::new(t), MolarThermalQuantityV1::new(value)))
                    .collect(),
            )
            .unwrap()
        })
        .collect();
    MulticomponentLiquidCaloricModelV1::new(law, sources).unwrap()
}

fn seed(n: usize, rich: f64) -> MulticomponentEquilibriumSeedV1 {
    let dilute = (1.0 - rich) / (n - 1) as f64;
    MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: (0..n)
            .map(|i| (0..n).map(|j| if i == j { rich } else { dilute }).collect())
            .collect(),
        phase_weights: vec![1.0; n],
    }
}

fn state(
    model: &MulticomponentLiquidCaloricModelV1,
    t: f64,
    amounts: &[f64],
    seed: &MulticomponentEquilibriumSeedV1,
) -> MulticomponentEquilibriumCaloricStateV1 {
    let candidate = model
        .thermodynamics()
        .equilibrium_candidate(
            Temperature::new(t),
            amounts,
            seed,
            MulticomponentEquilibriumControlV1 {
                chemical_potential_tolerance: 2e-13,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    model
        .equilibrium_caloric_state(candidate, Default::default(), || false)
        .unwrap()
}

fn controls(lo: f64, hi: f64) -> MulticomponentEquilibriumCaloricSolveV1 {
    MulticomponentEquilibriumCaloricSolveV1 {
        thermal: MulticomponentCaloricSolveV1 {
            temperature_bracket: [Temperature::new(lo), Temperature::new(hi)],
            temperature_tolerance_k: 1e-8,
            energy_tolerance_j: 1e-6,
            entropy_tolerance_j_per_k: 1e-8,
            maximum_iterations: 128,
        },
        response: Default::default(),
        maximum_temperature_step_k: 12.0,
        maximum_phase_solves: 4096,
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn conservation(state: &MulticomponentEquilibriumCaloricStateV1, offered: &[f64]) {
    assert_eq!(state.component_moles(), offered);
    for (i, &amount) in offered.iter().enumerate() {
        let recovered: f64 = state
            .candidate()
            .phases()
            .iter()
            .map(|phase| phase.component_moles()[i])
            .sum();
        near(recovered, amount, 0.0, 64.0 * f64::EPSILON);
    }
}

#[test]
fn g1_binary_finite_heating_and_cooling_match_original_equilibrium_calorics() {
    for (a, b, initial_t, final_t, bracket) in [
        (0.0, 1000.0, 280.0, 400.0, [260.0, 430.0]),
        (5.0, -900.0, 400.0, 600.0, [380.0, 650.0]),
    ] {
        let model = model(2, a, b, 100.0, true);
        let original = state(&model, initial_t, &[3.5, 6.5], &seed(2, 0.95));
        let binary = BinaryLiquidCaloricModelV1::new(
            BinaryRegularSolutionV1::new(
                [
                    model.thermodynamics().species()[0].clone(),
                    model.thermodynamics().species()[1].clone(),
                ],
                a,
                b,
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
        let exact = |t| {
            binary
                .equilibrate(Temperature::new(t), [3.5, 6.5], Default::default(), || {
                    false
                })
                .unwrap()
        };
        let initial_binary = exact(initial_t);
        let final_binary = exact(final_t);
        let heat = binary
            .enthalpy_change_j(&initial_binary, &final_binary)
            .unwrap();
        let entropy = binary
            .entropy_change_j_per_k(&initial_binary, &final_binary)
            .unwrap();
        let controls = controls(bracket[0], bracket[1]);
        let result = model
            .add_heat_equilibrium(&original, heat, controls, || false)
            .unwrap();
        near(
            result.final_state().temperature().value(),
            final_t,
            2e-8,
            0.0,
        );
        near(result.enthalpy_change_j(), heat, 1e-6, 0.0);
        near(result.entropy_change_j_per_k(), entropy, 2e-8, 2e-11);
        near(
            result.final_state().isobaric_heat_capacity_j_per_k(),
            final_binary.isobaric_heat_capacity_j_per_k().unwrap(),
            1e-7,
            2e-10,
        );
        near(
            result.residual_j(),
            0.0,
            controls.thermal.energy_tolerance_j,
            0.0,
        );
        assert_eq!(result.final_state().candidate().phases().len(), 2);
        assert!(result.completed_phase_solves() > 10);
        assert!(result.completed_phase_solves() <= controls.maximum_phase_solves);
        conservation(result.final_state(), original.component_moles());
        assert_eq!(original.temperature().value(), initial_t);
        let cooled = model
            .add_heat_equilibrium(result.final_state(), -heat, controls, || false)
            .unwrap();
        near(
            cooled.final_state().temperature().value(),
            initial_t,
            3e-8,
            0.0,
        );
        near(cooled.entropy_change_j_per_k(), -entropy, 3e-8, 2e-11);
        near(
            model
                .equilibrium_enthalpy_change_j(&original, result.final_state())
                .unwrap(),
            heat,
            1e-6,
            0.0,
        );
        near(
            model
                .equilibrium_entropy_change_j_per_k(&original, result.final_state())
                .unwrap(),
            entropy,
            2e-8,
            2e-11,
        );
        // This large thermal transaction also agrees with the independent
        // shipped binary owner's full nonlinear heat inverse.
        let reference = binary
            .add_heat(
                &initial_binary,
                heat,
                BinaryCaloricSolveV1 {
                    temperature_bracket: controls.thermal.temperature_bracket,
                    temperature_tolerance_k: controls.thermal.temperature_tolerance_k,
                    energy_tolerance_j: controls.thermal.energy_tolerance_j,
                    maximum_iterations: controls.thermal.maximum_iterations,
                    coexistence: BinaryCoexistenceOptionsV1::default(),
                },
                || false,
            )
            .unwrap();
        near(
            result.final_state().temperature().value(),
            reference.final_state().temperature().value(),
            2e-8,
            0.0,
        );
    }
}

// Independent scalar common tangent of three permutations of (r,d,d).
fn ternary_rich(chi: f64) -> f64 {
    let mut bracket = [2.0 / 3.0, 1.0 - 1e-14];
    let residual = |r: f64| (2.0 * r / (1.0 - r)).ln() - 0.5 * chi * (3.0 * r - 1.0);
    assert!(residual(bracket[0]) < 0.0 && residual(bracket[1]) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * bracket[0] + 0.5 * bracket[1];
        if residual(middle) > 0.0 {
            bracket[1] = middle;
        } else {
            bracket[0] = middle;
        }
    }
    0.5 * bracket[0] + 0.5 * bracket[1]
}

fn ternary_mixing(t: f64) -> [f64; 2] {
    let rich = ternary_rich(1600.0 / t);
    let dilute = 0.5 * (1.0 - rich);
    [
        R * 1600.0 * (2.0 * rich * dilute + dilute * dilute),
        -R * (rich * rich.ln() + 2.0 * dilute * dilute.ln()),
    ]
}

#[test]
fn g1_ternary_large_heat_closes_independent_energy_entropy_and_phase_amounts() {
    let amounts = [2.0, 3.0, 5.0];
    let initial_t = 330.0;
    let target_t = 480.0;
    let initial_mix = ternary_mixing(initial_t);
    let final_mix = ternary_mixing(target_t);
    let heat = 1320.0 * (target_t - initial_t) + 10.0 * (final_mix[0] - initial_mix[0]);
    let entropy = 1320.0 * (target_t / initial_t).ln() + 10.0 * (final_mix[1] - initial_mix[1]);
    for datum in [100.0, 1e20] {
        let model = model(3, 0.0, 1600.0, datum, false);
        let original = state(&model, initial_t, &amounts, &seed(3, 0.97));
        let result = model
            .add_heat_equilibrium(&original, heat, controls(310.0, 500.0), || false)
            .unwrap();
        let final_state = result.final_state();
        near(final_state.temperature().value(), target_t, 2e-8, 0.0);
        near(result.entropy_change_j_per_k(), entropy, 3e-8, 3e-11);
        near(result.enthalpy_change_j(), heat, 1e-6, 0.0);
        assert!(final_state.redistribution_heat_capacity_j_per_k() > 0.0);
        conservation(final_state, &amounts);
        let rich = ternary_rich(1600.0 / target_t);
        let dilute = 0.5 * (1.0 - rich);
        for (i, phase) in final_state.candidate().phases().iter().enumerate() {
            near(
                phase.total_moles(),
                (amounts[i] - 10.0 * dilute) / (rich - dilute),
                2e-10,
                2e-10,
            );
            near(
                phase.thermodynamics().mole_fractions()[i],
                rich,
                2e-11,
                2e-11,
            );
        }
        let lower_entropy = heat / target_t;
        let upper_entropy = heat / initial_t;
        assert!(result.entropy_change_j_per_k() > lower_entropy);
        assert!(result.entropy_change_j_per_k() < upper_entropy);
        if datum == 1e20 {
            // The entire 150 K transaction is below one ulp of the arbitrary
            // absolute reference H; direct increments still recover it.
            assert_eq!(original.enthalpy_j(), final_state.enthalpy_j());
            assert!(result.enthalpy_change_j() > 200_000.0);
        }
        let reverse_h = model
            .equilibrium_enthalpy_change_j(final_state, &original)
            .unwrap();
        let reverse_s = model
            .equilibrium_entropy_change_j_per_k(final_state, &original)
            .unwrap();
        assert_eq!(reverse_h, -result.enthalpy_change_j());
        assert_eq!(reverse_s, -result.entropy_change_j_per_k());
    }
}

#[test]
fn g1_finite_heat_preserves_inventory_scaling_and_phase_permutation() {
    let model = model(3, 0.0, 1600.0, 300.0, false);
    let mix0 = ternary_mixing(350.0);
    let mix1 = ternary_mixing(440.0);
    let heat = 1320.0 * 90.0 + 10.0 * (mix1[0] - mix0[0]);
    for scale in [1e-60, 1.0, 1e60] {
        let amounts = [2.0 * scale, 3.0 * scale, 5.0 * scale];
        let mut guess = seed(3, 0.95);
        guess.phase_mole_fractions.swap(0, 2);
        let initial = state(&model, 350.0, &amounts, &guess);
        let mut control = controls(330.0, 470.0);
        control.thermal.energy_tolerance_j *= scale;
        control.thermal.entropy_tolerance_j_per_k *= scale;
        let result = model
            .add_heat_equilibrium(&initial, heat * scale, control, || false)
            .unwrap();
        near(result.final_state().temperature().value(), 440.0, 2e-8, 0.0);
        near(result.enthalpy_change_j() / scale, heat, 1e-6, 0.0);
        conservation(result.final_state(), &amounts);
        let expected = ternary_rich(1600.0 / 440.0);
        near(
            result.final_state().candidate().phases()[0]
                .thermodynamics()
                .mole_fractions()[2],
            expected,
            2e-11,
            2e-11,
        );
    }
}

#[test]
fn g0_one_phase_uses_original_piecewise_cp_and_retains_absent_species() {
    let model = model(3, 0.0, 0.0, 100.0, true);
    let amounts = [1.0, 1e-24, 0.0];
    let initial = state(
        &model,
        250.0,
        &amounts,
        &MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![amounts.to_vec()],
            phase_weights: vec![1.0],
        },
    );
    let homogeneous = model
        .homogeneous(Temperature::new(250.0), &amounts, || false)
        .unwrap();
    let target = model
        .homogeneous(Temperature::new(550.0), &amounts, || false)
        .unwrap();
    let heat = model.enthalpy_change_j(&homogeneous, &target).unwrap();
    let result = model
        .add_heat_equilibrium(&initial, heat, controls(200.0, 600.0), || false)
        .unwrap();
    near(result.final_state().temperature().value(), 550.0, 2e-8, 0.0);
    assert_eq!(
        result.final_state().redistribution_heat_capacity_j_per_k(),
        0.0
    );
    assert_eq!(
        result.final_state().candidate().phases()[0].component_moles(),
        amounts
    );
    near(
        result.entropy_change_j_per_k(),
        model.entropy_change_j_per_k(&homogeneous, &target).unwrap(),
        1e-8,
        2e-12,
    );
}

#[test]
fn g1_two_phase_heating_and_cooling_move_real_trace_amounts_without_a_floor() {
    let reference = model(3, 0.0, 0.0, 300.0, false);
    let law = MulticomponentRegularSolutionV1::new(
        reference.thermodynamics().species().to_vec(),
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 800.0, 0.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let model =
        MulticomponentLiquidCaloricModelV1::new(law, reference.heat_capacities().to_vec()).unwrap();
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
    let before = binary
        .equilibrate(
            Temperature::new(350.0),
            [3.5, 6.5],
            Default::default(),
            || false,
        )
        .unwrap();
    let after = binary
        .equilibrate(
            Temperature::new(480.0),
            [3.5, 6.5],
            Default::default(),
            || false,
        )
        .unwrap();
    let heat = binary.enthalpy_change_j(&before, &after).unwrap();
    let entropy = binary.entropy_change_j_per_k(&before, &after).unwrap();
    for trace in [0.0, 1e-24] {
        let amounts = [3.5, 6.5, trace];
        let initial = state(
            &model,
            350.0,
            &amounts,
            &MulticomponentEquilibriumSeedV1 {
                phase_mole_fractions: vec![
                    vec![0.95, 0.05, trace / 10.0],
                    vec![0.05, 0.95, trace / 10.0],
                ],
                phase_weights: vec![1.0; 2],
            },
        );
        let control = controls(330.0, 500.0);
        let heated = model
            .add_heat_equilibrium(&initial, heat, control, || false)
            .unwrap();
        near(heated.final_state().temperature().value(), 480.0, 2e-8, 0.0);
        near(heated.entropy_change_j_per_k(), entropy, 3e-8, 2e-11);
        let cooled = model
            .add_heat_equilibrium(heated.final_state(), -heat, control, || false)
            .unwrap();
        near(cooled.final_state().temperature().value(), 350.0, 2e-8, 0.0);
        for endpoint in [heated.final_state(), cooled.final_state()] {
            conservation(endpoint, &amounts);
            for phase in endpoint.candidate().phases() {
                if trace == 0.0 {
                    assert_eq!(phase.component_moles()[2], 0.0);
                } else {
                    assert!(phase.component_moles()[2] > 0.0);
                    // The symmetric inert tracer shares the same concentration
                    // in both phases and follows their changing finite amounts.
                    near(
                        phase.component_moles()[2] / phase.total_moles(),
                        trace / 10.0,
                        0.0,
                        2e-10,
                    );
                }
            }
        }
        if trace != 0.0 {
            assert_ne!(
                initial.candidate().phases()[0].component_moles()[2],
                heated.final_state().candidate().phases()[0].component_moles()[2]
            );
        }
    }
}

#[test]
fn g3_original_owner_controls_budget_and_cancellation_refuse_atomically() {
    let model = model(3, 0.0, 1600.0, 300.0, false);
    let initial = state(&model, 350.0, &[2.0, 3.0, 5.0], &seed(3, 0.95));
    let retained = initial.clone();
    let control = controls(330.0, 470.0);
    let exhausted = MulticomponentEquilibriumCaloricSolveV1 {
        maximum_phase_solves: 1,
        ..control
    };
    assert!(matches!(
        model.add_heat_equilibrium(&initial, 20_000.0, exhausted, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::PhaseSolveBudget {
            completed_phase_solves: 1
        })
    ));
    assert!(matches!(
        model.add_heat_equilibrium(&initial, 20_000.0, control, || true),
        Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        model.add_heat_equilibrium(&initial, 20_000.0, control, || {
            polls += 1;
            polls == 250
        }),
        Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
    ));
    assert_eq!(polls, 250);
    for value in [0.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            model.add_heat_equilibrium(
                &initial,
                20_000.0,
                MulticomponentEquilibriumCaloricSolveV1 {
                    maximum_temperature_step_k: value,
                    ..control
                },
                || false
            ),
            Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
        ));
    }
    let foreign = super::tests::model(3, 0.0, 1600.0, 350.0, false);
    assert!(matches!(
        foreign.add_heat_equilibrium(&initial, 20_000.0, control, || false),
        Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
    ));
    let foreign_potential = super::tests::model(3, 0.1, 1600.0, 300.0, false);
    assert!(matches!(
        foreign_potential.equilibrium_caloric_state(
            initial.candidate().clone(),
            Default::default(),
            || false
        ),
        Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(_))
    ));
    let unchanged = model
        .add_heat_equilibrium(&initial, 0.0, control, || false)
        .unwrap();
    assert_eq!(unchanged.final_state(), &initial);
    assert_eq!(unchanged.completed_phase_solves(), 0);
    assert_eq!(unchanged.enthalpy_change_j(), 0.0);
    assert_eq!(unchanged.entropy_change_j_per_k(), 0.0);
    assert_eq!(initial, retained);
}

#[test]
fn g3_a_thermal_bracket_crossing_phase_loss_publishes_no_frozen_cp_endpoint() {
    let model = model(2, 0.0, 1000.0, 300.0, false);
    let initial = state(&model, 350.0, &[2.5, 7.5], &seed(2, 0.95));
    let retained = initial.clone();
    let error = model
        .add_heat_equilibrium(&initial, 200_000.0, controls(330.0, 650.0), || false)
        .unwrap_err();
    assert!(
        matches!(
            error,
            MulticomponentEquilibriumCaloricErrorV1::Equilibrium(_)
                | MulticomponentEquilibriumCaloricErrorV1::Response(_)
                | MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(_)
                | MulticomponentEquilibriumCaloricErrorV1::PhaseSolveBudget { .. }
        ),
        "unexpected refusal: {error:?}"
    );
    assert_eq!(initial, retained);
    assert_eq!(initial.candidate().phases().len(), 2);
}
