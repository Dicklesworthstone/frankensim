//! Physical time error includes latent storage and charges discarded proposals.
use super::{fixture::*, implicit_time::PhysicalBalance};
use fs_blake3::ContentHash;
use fs_couple::finite_gas_thermal::wall_film::{
    implicit::{WallFilmImplicitControls, WallFilmImplicitError},
    window::{
        WallFilmImplicitWindowControls, WallFilmWindowError, WallFilmWindowSolve,
        adaptive::{WallFilmAdaptiveControls, WallFilmAdaptiveError, WallFilmAdaptiveReport},
    },
};
use fs_material::{
    evaporation::gas::{GasLiquidBinding, wall_film::GasWallFilmInventory},
    phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStep,
        UniformEnthalpyStepInput,
    },
    thermochemical_gas::energy::parcel::ClosedGasParcel,
};

fn policy(span: f64) -> WallFilmAdaptiveControls {
    let integration = WallFilmImplicitWindowControls {
        window: WallFilmWindowSolve {
            maximum_step_s: span,
            minimum_step_s: 1e-7,
            maximum_attempts: 20_000,
            energy_tolerance_j: 1e-6,
            mass_tolerance_kg: 1e-13,
        },
        solve: WallFilmImplicitControls {
            maximum_iterations: 10_000,
            maximum_evaluations: 100_000,
            ..WallFilmImplicitControls::default()
        },
        maximum_gas_iterations: 100_000,
    };
    let mut controls = WallFilmAdaptiveControls::new(integration);
    controls.relative_tolerance = 0.0;
    controls.enthalpy_absolute_tolerance_j_kg = 1e4;
    controls.temperature_absolute_tolerance_k = 1.0;
    controls.phase_fraction_absolute_tolerance = 1.0;
    controls.mass_absolute_tolerance_kg = 1e-6;
    controls
}

fn balance(
    input: UniformEnthalpyStepInput<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &GasLiquidBinding,
    gas_heat: f64,
    next: &UniformEnthalpyStep<WallFilmAdaptiveReport>,
) {
    let report = &next.report.window.physical;
    let wall_delta = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    let gas_delta = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_delta = film_energy(report.next_film, binding) - film_energy(film, binding);
    near(
        wall_delta + gas_delta + film_delta,
        input.internal_heat_j + gas_heat + report.radiation_heat_j,
        1e-6,
    );
    near(report.internal_body_heat_j, input.internal_heat_j, 1e-14);
    near(report.external_gas_heat_j, gas_heat, 1e-14);
    near(
        report.next_film.mass_kg() + report.next_gas.species_mass_kg(&water()).unwrap(),
        film.mass_kg() + gas.species_mass_kg(&water()).unwrap(),
        1e-13,
    );
    near(
        report.next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-17,
    );
    assert!(report.sum_absolute_energy_residual_j <= 1e-6);
    assert!(next.report.window.sum_absolute_mass_residual_kg <= 1e-13);
    assert!(next.report.maximum_accepted_normalized_error <= 1.0);
    assert_eq!(
        report.accepted_intervals,
        2 * next.report.accepted_comparisons
    );
    assert_eq!(report.samples.len(), report.accepted_intervals + 1);
    assert_eq!(report.samples.last().unwrap().time_s, input.duration_s);
    assert_eq!(
        report.attempted_intervals,
        report.accepted_intervals + report.rejected_intervals
    );
}

#[test]
fn tighter_physical_tolerances_improve_wall_film_and_gas_against_continuous_balance() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let input = input(&curve, 300.0, 0.1);
    let film = GasWallFilmInventory::Wet {
        mass_kg: 1e-6,
        specific_enthalpy_j_kg: 0.0,
    };
    let environment = environment(&gas, &source, film, &binding, &body);
    let reference = PhysicalBalance {
        wall_capacity: input.mass_kg * CP_WALL,
        carrier_mass: gas.species_mass_kg(&carrier()).unwrap(),
        total_volatile_mass: film.mass_kg() + gas.species_mass_kg(&water()).unwrap(),
        sensible_conductance: body.convection_w_per_m2_k * AREA,
        phase_exposure: environment.boundary.mass_transfer_coefficient_m_s * AREA,
    }
    .reference([300.0, film.mass_kg(), 320.0], input.duration_s);
    let loose = policy(input.duration_s);
    let mut tight = loose;
    tight.enthalpy_absolute_tolerance_j_kg = 0.8;
    tight.temperature_absolute_tolerance_k = 1e-3;
    tight.mass_absolute_tolerance_kg = 1e-11;
    let (coarse, fine) = scope(|cx, _| {
        (
            environment
                .advance_implicit_adaptive(cx, input, loose)
                .unwrap(),
            environment
                .advance_implicit_adaptive(cx, input, tight)
                .unwrap(),
        )
    });
    let state = |next: &UniformEnthalpyStep<WallFilmAdaptiveReport>| {
        [
            next.state.temperature_k(),
            next.report.window.physical.next_film.mass_kg(),
            next.report
                .window
                .physical
                .next_gas
                .state()
                .state()
                .temperature,
        ]
    };
    for (i, (coarse_value, fine_value)) in state(&coarse).into_iter().zip(state(&fine)).enumerate()
    {
        let loose_error = (coarse_value - reference[i]).abs();
        let tight_error = (fine_value - reference[i]).abs();
        assert!(
            tight_error < 0.2 * loose_error,
            "component {i}: loose error {loose_error:e}, tight error {tight_error:e}, reference {reference:?}"
        );
    }
    assert_eq!(coarse.report.accepted_comparisons, 1);
    assert!(fine.report.accuracy_rejections > 0);
    assert!(fine.report.accepted_comparisons > coarse.report.accepted_comparisons);
    for next in [&coarse, &fine] {
        balance(input, &gas, film, &binding, 0.0, next);
        assert_eq!(next.report.admission_rejections, 0);
        assert_eq!(
            next.report.window.physical.attempted_intervals,
            3 * next.report.comparisons
        );
        assert_eq!(
            next.report.comparisons,
            next.report.accepted_comparisons + next.report.accuracy_rejections
        );
    }
}

#[test]
fn wall_phase_fraction_alone_detects_latent_time_error_at_constant_temperature() {
    const LATENT_WALL: f64 = 100_000.0;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x86; 32]),
        [
            (280.0, -20.0 * CP_WALL, 0.0),
            (300.0, 0.0, 0.0),
            (300.0, LATENT_WALL, 1.0),
            (340.0, LATENT_WALL + 40.0 * CP_WALL, 1.0),
        ]
        .map(
            |(temperature_k, specific_enthalpy_j_kg, liquid_mass_fraction)| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
        )
        .to_vec(),
    )
    .unwrap();
    let input = UniformEnthalpyStepInput {
        initial: curve.state_at_specific_enthalpy(0.5 * LATENT_WALL).unwrap(),
        ..input(&curve, 300.0, 0.125)
    };
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let binding = binding(&source);
    let film = GasWallFilmInventory::Wet {
        mass_kg: 1e-6,
        specific_enthalpy_j_kg: 0.0,
    };
    let body = body(0.0);
    let environment = environment(&gas, &source, film, &binding, &body);
    let mut loose = policy(input.duration_s);
    loose.enthalpy_absolute_tolerance_j_kg = 1e12;
    loose.temperature_absolute_tolerance_k = 1e12;
    loose.mass_absolute_tolerance_kg = 1.0;
    let mut tight = loose;
    tight.phase_fraction_absolute_tolerance = 1e-6;
    let (coarse, fine) = scope(|cx, _| {
        (
            environment
                .advance_implicit_adaptive(cx, input, loose)
                .unwrap(),
            environment
                .advance_implicit_adaptive(cx, input, tight)
                .unwrap(),
        )
    });
    let first_error = coarse.report.last_error.unwrap();
    assert_eq!(first_error.wall_temperature_difference_k, 0.0);
    assert!(first_error.wall_phase_fraction_difference > tight.phase_fraction_absolute_tolerance);
    assert_eq!(coarse.report.accepted_comparisons, 1);
    assert!(fine.report.accuracy_rejections > 0);
    assert!(fine.report.accepted_comparisons > 1);
    for next in [&coarse, &fine] {
        let report = &next.report.window.physical;
        assert!(
            report
                .samples
                .iter()
                .all(|sample| sample.wall_temperature_k == 300.0)
        );
        assert!(report.next_film.mass_kg() > 0.0);
        assert_eq!(report.condensed_mass_kg, 0.0);
        near(
            next.state.liquid_mass_fraction(),
            0.5 - report.evaporated_mass_kg * LATENT / (input.mass_kg * LATENT_WALL),
            1e-10,
        );
        balance(input, &gas, film, &binding, 0.0, next);
    }
    assert!(
        fine.report
            .last_error
            .unwrap()
            .wall_phase_fraction_difference
            <= tight.phase_fraction_absolute_tolerance
    );
}

#[test]
fn only_fine_sources_and_states_commit_while_coarse_and_rejected_work_spend_shared_budgets() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(10.0);
    let mut input = input(&curve, 300.0, 0.125);
    input.internal_heat_j = 0.5;
    let mut environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    environment.gas_heat_j = 0.125;
    let controls = policy(input.duration_s);
    let (accepted, coarse, first, fine) = scope(|cx, _| {
        let accepted = environment
            .advance_implicit_adaptive(cx, input, controls)
            .unwrap();
        let coarse = environment
            .advance_implicit(cx, input, controls.integration.solve)
            .unwrap();
        let half_input = UniformEnthalpyStepInput {
            duration_s: 0.5 * input.duration_s,
            internal_heat_j: 0.5 * input.internal_heat_j,
            ..input
        };
        let mut half_environment = environment.clone();
        half_environment.gas_heat_j *= 0.5;
        let first = half_environment
            .advance_implicit(cx, half_input, controls.integration.solve)
            .unwrap();
        let mut last_environment = half_environment.clone();
        last_environment.gas = &first.report.next_gas;
        last_environment.film = first.report.next_film;
        let fine = last_environment
            .advance_implicit(
                cx,
                UniformEnthalpyStepInput {
                    initial: first.state,
                    ..half_input
                },
                controls.integration.solve,
            )
            .unwrap();
        (accepted, coarse, first, fine)
    });
    let report = &accepted.report.window;
    let physical = &report.physical;
    assert_eq!(accepted.report.comparisons, 1);
    assert_eq!(physical.accepted_intervals, 2);
    assert_eq!(physical.attempted_intervals, 3);
    assert_eq!(physical.rejected_intervals, 1);
    near(
        accepted.state.specific_enthalpy_j_kg(),
        fine.state.specific_enthalpy_j_kg(),
        1e-6,
    );
    near(
        physical.next_gas.state().state().temperature,
        fine.report.next_gas.state().state().temperature,
        1e-7,
    );
    near(
        physical.next_film.mass_kg(),
        fine.report.next_film.mass_kg(),
        1e-14,
    );
    near(
        physical.vapor_energy_into_gas_j,
        first.report.phase.energy_into_gas_j + fine.report.phase.energy_into_gas_j,
        1e-8,
    );
    near(
        physical.latent_energy_j,
        first.report.phase.latent_energy_j + fine.report.phase.latent_energy_j,
        1e-8,
    );
    assert_eq!(
        report.residual_evaluations,
        coarse.report.residual_evaluations
            + first.report.residual_evaluations
            + fine.report.residual_evaluations
    );
    assert_eq!(
        report.gas_iterations,
        coarse.report.gas_iterations + first.report.gas_iterations + fine.report.gas_iterations
    );
    balance(
        input,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        environment.gas_heat_j,
        &accepted,
    );

    // A budget sufficient for the retained physical pair must still pay for
    // the discarded coarse proposal. It cannot be renewed after comparison.
    let mut limited = controls;
    limited.integration.solve.maximum_evaluations =
        first.report.residual_evaluations + fine.report.residual_evaluations;
    assert!(matches!(
        scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::Window(
            WallFilmWindowError::Implicit(WallFilmImplicitError::EvaluationBudget)
        ))
    ));
    limited = controls;
    limited.integration.maximum_gas_iterations =
        first.report.gas_iterations + fine.report.gas_iterations;
    assert!(scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)).is_err());
    limited = controls;
    limited.integration.window.maximum_attempts = 2;
    assert!(matches!(
        scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::Window(
            WallFilmWindowError::WorkLimit {
                accepted_time_s: 0.0,
                attempted_intervals: 2
            }
        ))
    ));

    limited = controls;
    limited.mass_absolute_tolerance_kg = 1e-14;
    limited.maximum_comparisons = 1;
    assert!(matches!(
        scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::ComparisonBudget {
            accepted_time_s: 0.0,
            comparisons: 1
        })
    ));
    limited.maximum_reductions = 0;
    assert!(
        matches!(scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::TimeAccuracy { accepted_time_s: 0.0, error, .. }) if error.normalized > 1.0)
    );
    scope(|cx, gate| {
        gate.request();
        assert!(matches!(
            environment.advance_implicit_adaptive(cx, input, controls),
            Err(WallFilmAdaptiveError::Window(
                WallFilmWindowError::Cancelled
            ))
        ));
    });
    assert_eq!(gas, original_gas);
    assert_eq!(environment.film, GasWallFilmInventory::Dry);
    assert_eq!(input.initial.specific_enthalpy_j_kg(), 0.0);
}

#[test]
fn actual_donor_admission_retries_keep_heat_totals_and_spend_shared_work() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let mut input = input(&curve, 300.0, 0.2);
    input.internal_heat_j = 0.5;
    let mut environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    environment.gas_heat_j = 0.2;
    environment.solve.phase.maximum_inventory_fraction = 0.01;
    let controls = policy(input.duration_s);
    let accepted = scope(|cx, _| {
        environment
            .advance_implicit_adaptive(cx, input, controls)
            .unwrap()
    });
    assert!(accepted.report.admission_rejections > 0);
    assert!(accepted.report.accepted_comparisons > 1);
    assert!(accepted.report.window.physical.next_film.mass_kg() > 0.0);
    assert!(accepted.report.window.physical.condensed_mass_kg > 0.0);
    balance(
        input,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        environment.gas_heat_j,
        &accepted,
    );

    let mut limited = controls;
    limited.integration.window.maximum_attempts = 1;
    assert!(matches!(
        scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::Window(
            WallFilmWindowError::WorkLimit {
                accepted_time_s: 0.0,
                attempted_intervals: 1,
            }
        ))
    ));
    limited = controls;
    limited.integration.solve.maximum_evaluations = 1;
    assert!(matches!(
        scope(|cx, _| environment.advance_implicit_adaptive(cx, input, limited)),
        Err(WallFilmAdaptiveError::Window(
            WallFilmWindowError::Implicit(WallFilmImplicitError::EvaluationBudget)
        ))
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(environment.film, GasWallFilmInventory::Dry);
    assert_eq!(input.initial.specific_enthalpy_j_kg(), 0.0);
}
