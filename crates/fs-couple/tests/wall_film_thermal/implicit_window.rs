//! Whole physical windows retain the joint wall/film/gas endpoint and clock.
use super::fixture::*;
use std::cell::Cell;

use fs_blake3::ContentHash;
use fs_couple::finite_gas_thermal::wall_film::{
    WallFilmThermalError,
    implicit::{WallFilmImplicitControls, WallFilmImplicitError},
    window::{
        WallFilmImplicitWindowControls, WallFilmImplicitWindowReport, WallFilmWindowError,
        WallFilmWindowSolve,
    },
};
use fs_material::{
    evaporation::gas::{GasEvaporationError, GasLiquidBinding, wall_film::GasWallFilmInventory},
    phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStep,
        UniformEnthalpyStepInput,
    },
    thermochemical_gas::energy::parcel::ClosedGasParcel,
};
use fs_thermal::evaporation::EvaporationError;

fn window(step: f64) -> WallFilmImplicitWindowControls {
    WallFilmImplicitWindowControls {
        window: WallFilmWindowSolve {
            maximum_step_s: step,
            minimum_step_s: 1e-7,
            maximum_attempts: 2000,
            energy_tolerance_j: 1e-6,
            mass_tolerance_kg: 1e-13,
        },
        solve: WallFilmImplicitControls {
            maximum_iterations: 4096,
            maximum_evaluations: 100_000,
            ..WallFilmImplicitControls::default()
        },
        maximum_gas_iterations: 100_000,
    }
}

fn balances(
    input: UniformEnthalpyStepInput<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &GasLiquidBinding,
    gas_heat: f64,
    next: &UniformEnthalpyStep<WallFilmImplicitWindowReport>,
) {
    let report = &next.report.physical;
    let wall_change = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    let gas_change = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_change = film_energy(report.next_film, binding) - film_energy(film, binding);
    near(
        wall_change + gas_change + film_change,
        input.internal_heat_j + gas_heat + report.radiation_heat_j,
        1e-6,
    );
    near(
        wall_change,
        input.internal_heat_j + next.external_heat_j,
        1e-6,
    );
    near(report.internal_body_heat_j, input.internal_heat_j, 1e-14);
    near(report.external_gas_heat_j, gas_heat, 1e-14);
    near(report.body_heat_apportionment_residual_j, 0.0, 1e-14);
    near(report.gas_heat_apportionment_residual_j, 0.0, 1e-14);
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
    near(
        report.next_film.mass_kg() + report.next_gas.mass_kg(),
        film.mass_kg() + gas.mass_kg(),
        1e-13,
    );
    assert!(report.sum_absolute_energy_residual_j <= 1e-6);
    assert!(next.report.sum_absolute_mass_residual_kg <= 1e-13);
    assert!(next.report.charged_gas_iterations >= next.report.gas_iterations);
    assert_eq!(report.samples.first().unwrap().time_s, 0.0);
    assert_eq!(report.samples.last().unwrap().time_s, input.duration_s);
    assert_eq!(report.samples.len(), report.accepted_intervals + 1);
    assert_eq!(
        report.attempted_intervals,
        report.accepted_intervals + report.rejected_intervals
    );
    assert!(report.samples.windows(2).all(|s| s[0].time_s < s[1].time_s));
    assert_eq!(report.next_gas.volume_m3(), gas.volume_m3());
    assert_eq!(
        report.next_gas.point().axes().get("wall-film-lot"),
        Some(&13.0)
    );
}

#[test]
fn simultaneous_source_window_matches_manual_endpoint_composition_on_the_same_clock() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(10.0);
    let mut input = input(&curve, 300.0, 0.3125);
    input.internal_heat_j = 0.5;
    let mut environment = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    environment.gas_heat_j = 0.125;
    let controls = window(0.125);
    let next = scope(|cx, _| {
        environment
            .advance_implicit_window(cx, input, controls)
            .unwrap()
    });
    let report = &next.report.physical;
    assert_eq!(report.accepted_intervals, 3);
    assert_eq!(report.rejected_intervals, 0);
    assert!(report.condensed_mass_kg > 0.0);
    assert!(report.solid_to_gas_heat_j < 0.0);
    assert_eq!(
        report.samples.iter().map(|s| s.time_s).collect::<Vec<_>>(),
        [0.0, 0.125, 0.25, 0.3125]
    );

    let (manual_state, manual_gas, manual_film, transferred) = scope(|cx, _| {
        let mut state = input.initial;
        let mut gas = gas.clone();
        let mut film = GasWallFilmInventory::Dry;
        let mut previous = 0.0;
        let mut phase = [0.0; 3];
        for end in [0.125, 0.25, input.duration_s] {
            // Independent cumulative source endpoints: no heat is renewed by
            // the final half-length interval or by the nonlinear iterations.
            let fraction = end / input.duration_s - previous / input.duration_s;
            let mut interval = environment.clone();
            interval.gas = &gas;
            interval.film = film;
            interval.gas_heat_j = environment.gas_heat_j * fraction;
            let accepted = interval
                .advance_implicit(
                    cx,
                    UniformEnthalpyStepInput {
                        initial: state,
                        duration_s: end - previous,
                        internal_heat_j: input.internal_heat_j * fraction,
                        ..input
                    },
                    controls.solve,
                )
                .unwrap();
            phase[0] += accepted.report.phase.energy_into_gas_j;
            phase[1] += accepted.report.phase.latent_energy_j;
            phase[2] += accepted.report.solid_to_gas_heat_j;
            state = accepted.state;
            gas = accepted.report.next_gas;
            film = accepted.report.next_film;
            previous = end;
        }
        (state, gas, film, phase)
    });
    near(
        next.state.specific_enthalpy_j_kg(),
        manual_state.specific_enthalpy_j_kg(),
        1e-6,
    );
    near(
        report.next_gas.state().state().temperature,
        manual_gas.state().state().temperature,
        1e-7,
    );
    near(report.next_film.mass_kg(), manual_film.mass_kg(), 1e-14);
    near(
        report.next_film.specific_enthalpy_j_kg().unwrap(),
        manual_film.specific_enthalpy_j_kg().unwrap(),
        1e-5,
    );
    near(report.vapor_energy_into_gas_j, transferred[0], 1e-8);
    near(report.latent_energy_j, transferred[1], 1e-8);
    near(report.solid_to_gas_heat_j, transferred[2], 1e-8);
    balances(
        input,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        environment.gas_heat_j,
        &next,
    );
}

#[test]
fn one_finite_inventory_survives_coupled_onset_dryout_and_recondensation() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let cold_input = input(&curve, 300.0, 0.1);
    let (cold, hot, recooled) = scope(|cx, _| {
        let cold = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body)
            .advance_implicit_window(cx, cold_input, window(0.025))
            .unwrap();
        let cold_report = &cold.report.physical;
        assert!(cold_report.next_film.mass_kg() > 0.0);
        assert!(cold_report.samples.iter().any(|s| s.film_formed));
        balances(
            cold_input,
            &gas,
            GasWallFilmInventory::Dry,
            &binding,
            0.0,
            &cold,
        );

        let warm_input = UniformEnthalpyStepInput {
            initial: cold.state,
            duration_s: 2.0,
            internal_heat_j: 40.0,
            ..cold_input
        };
        let hot = environment(
            &cold_report.next_gas,
            &source,
            cold_report.next_film,
            &binding,
            &body,
        )
        .advance_implicit_window(cx, warm_input, window(0.125))
        .unwrap();
        let hot_report = &hot.report.physical;
        assert_eq!(hot_report.next_film, GasWallFilmInventory::Dry);
        assert!(hot_report.evaporated_mass_kg >= cold_report.next_film.mass_kg());
        let event_index = hot_report
            .samples
            .iter()
            .position(|s| s.dry_out_time_s.is_some())
            .unwrap();
        let exhaustion = hot_report.samples[event_index].dry_out_time_s.unwrap();
        // This is the accepted endpoint-rate exhaustion estimate on the
        // whole-window clock, not an exact continuous event-time oracle.
        assert!(exhaustion > hot_report.samples[event_index - 1].time_s);
        assert!(exhaustion <= hot_report.samples[event_index].time_s);
        assert!(
            hot_report.samples[event_index..]
                .iter()
                .all(|s| s.film_mass_kg == 0.0)
        );
        balances(
            warm_input,
            &cold_report.next_gas,
            cold_report.next_film,
            &binding,
            0.0,
            &hot,
        );

        let cooling_input = UniformEnthalpyStepInput {
            initial: hot.state,
            duration_s: 2.0,
            internal_heat_j: -40.0,
            ..cold_input
        };
        let recooled = environment(
            &hot_report.next_gas,
            &source,
            hot_report.next_film,
            &binding,
            &body,
        )
        .advance_implicit_window(cx, cooling_input, window(0.125))
        .unwrap();
        assert!(recooled.report.physical.next_film.mass_kg() > 0.0);
        assert!(
            recooled
                .report
                .physical
                .samples
                .iter()
                .any(|s| s.film_formed)
        );
        assert!(recooled.state.temperature_k() < hot.state.temperature_k());
        balances(
            cooling_input,
            &hot_report.next_gas,
            hot_report.next_film,
            &binding,
            0.0,
            &recooled,
        );
        (cold, hot, recooled)
    });
    for next in [&cold, &hot, &recooled] {
        let report = &next.report.physical;
        near(
            report.next_gas.species_mass_kg(&water()).unwrap() + report.next_film.mass_kg(),
            gas.species_mass_kg(&water()).unwrap(),
            1e-13,
        );
    }
    near(
        hot.report.physical.internal_body_heat_j + recooled.report.physical.internal_body_heat_j,
        0.0,
        1e-14,
    );
    assert_eq!(gas, original_gas);
}

#[test]
fn native_wall_plateau_pays_repeated_evaporation_and_retains_exact_dry_state() {
    const WALL_LATENT: f64 = 100_000.0;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x85; 32]),
        [
            (280.0, -20.0 * CP_WALL, 0.0),
            (300.0, 0.0, 0.0),
            (300.0, WALL_LATENT, 1.0),
            (340.0, WALL_LATENT + 40.0 * CP_WALL, 1.0),
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
        initial: curve.state_at_specific_enthalpy(0.5 * WALL_LATENT).unwrap(),
        ..input(&curve, 300.0, 0.5)
    };
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let binding = binding(&source);
    let film_mass = 1e-6;
    let film = GasWallFilmInventory::Wet {
        mass_kg: film_mass,
        specific_enthalpy_j_kg: 0.0,
    };
    let body = body(0.0);
    let next = scope(|cx, _| {
        environment(&gas, &source, film, &binding, &body)
            .advance_implicit_window(cx, input, window(0.0625))
            .unwrap()
    });
    let report = &next.report.physical;
    assert_eq!(report.accepted_intervals, 8);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.next_film, GasWallFilmInventory::Dry);
    assert!(report.samples[1].film_mass_kg > 0.0);
    assert!(report.samples.iter().all(|s| s.wall_temperature_k == 300.0));
    near(report.evaporated_mass_kg, film_mass, 1e-14);
    assert_eq!(report.condensed_mass_kg, 0.0);
    let released = film_mass * LATENT;
    near(report.latent_energy_j, released, 1e-8);
    near(report.wall_heat_into_phase_system_j, released, 1e-8);
    near(next.external_heat_j, -released, 1e-8);
    near(
        next.state.liquid_mass_fraction(),
        0.5 - released / (input.mass_kg * WALL_LATENT),
        1e-10,
    );
    assert_eq!(next.state.phase_curve_identity(), curve.identity());

    // Every vapor increment leaves at exactly 300 K. Integrating the real
    // NASA9 Cp/R=3.5 stream and Cv/R=2.5 gas is independent of subdivision.
    let gas_capacity = 2.5 * 100_000.0 * GAS_VOLUME / 290.0;
    let gas_t = 290.0
        + film_mass * (R * 290.0 + 3.5 * R * (300.0 - 290.0))
            / (gas_capacity + film_mass * 2.5 * R);
    near(report.next_gas.state().state().temperature, gas_t, 2e-7);
    near(
        report.vapor_energy_into_gas_j,
        film_mass * (-13.4e6 + 3.5 * R * 300.0),
        1e-8,
    );
    assert_eq!(report.solid_to_gas_heat_j, 0.0);
    assert_eq!(report.radiation_heat_j, 0.0);
    balances(input, &gas, film, &binding, 0.0, &next);
}

#[test]
fn late_whole_window_work_source_and_cancellation_refusals_publish_no_partial_state() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let input = input(&curve, 300.0, 0.375);
    let initial = input.initial;
    let base = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    let controls = window(0.125);
    let first = scope(|cx, _| {
        base.advance_implicit(
            cx,
            UniformEnthalpyStepInput {
                duration_s: 0.125,
                ..input
            },
            controls.solve,
        )
        .unwrap()
    });
    let mut one = controls;
    one.window.maximum_attempts = 1;
    assert!(matches!(
        scope(|cx, _| base.advance_implicit_window(cx, input, one)),
        Err(WallFilmWindowError::WorkLimit {
            accepted_time_s: 0.125,
            ..
        })
    ));
    let mut evaluations = controls;
    evaluations.solve.maximum_evaluations = first.report.residual_evaluations;
    assert!(scope(|cx, _| base.advance_implicit_window(cx, input, evaluations)).is_err());
    let mut gas_work = controls;
    gas_work.maximum_gas_iterations = first.report.gas_iterations;
    assert!(scope(|cx, _| base.advance_implicit_window(cx, input, gas_work)).is_err());

    let queries = Cell::new(0usize);
    let counted_h = |t| {
        queries.set(queries.get() + 1);
        liquid_h(t)
    };
    let mut counted = base.clone();
    counted.boundary.liquid_enthalpy_at_temperature = &counted_h;
    let accepted = scope(|cx, _| {
        counted
            .advance_implicit_window(cx, input, controls)
            .unwrap()
    });
    assert_eq!(accepted.report.physical.accepted_intervals, 3);
    let final_query = queries.get();
    assert!(final_query > 4);

    // Trigger at the last original-source observation of a previously
    // successful complete window, after earlier speculative intervals exist.
    queries.set(0);
    let refuse_h = |t| {
        queries.set(queries.get() + 1);
        if queries.get() == final_query {
            Err(EvaporationError::Invalid(
                "late implicit window liquid source",
            ))
        } else {
            liquid_h(t)
        }
    };
    let mut refused = base.clone();
    refused.boundary.liquid_enthalpy_at_temperature = &refuse_h;
    assert!(matches!(
        scope(|cx, _| refused.advance_implicit_window(cx, input, controls)),
        Err(WallFilmWindowError::Implicit(
            WallFilmImplicitError::Physical(WallFilmThermalError::Phase(
                GasEvaporationError::Surface(EvaporationError::Invalid(
                    "late implicit window liquid source"
                ))
            ))
        ))
    ));
    assert_eq!(queries.get(), final_query);
    scope(|cx, gate| {
        queries.set(0);
        let cancel_h = |t| {
            queries.set(queries.get() + 1);
            if queries.get() == final_query {
                gate.request();
            }
            liquid_h(t)
        };
        let mut cancelled = base.clone();
        cancelled.boundary.liquid_enthalpy_at_temperature = &cancel_h;
        assert!(matches!(
            cancelled.advance_implicit_window(cx, input, controls),
            Err(WallFilmWindowError::Cancelled)
        ));
        assert_eq!(queries.get(), final_query);
    });
    assert_eq!(gas, original_gas);
    assert_eq!(input.initial, initial);
    assert_eq!(base.film, GasWallFilmInventory::Dry);
}
