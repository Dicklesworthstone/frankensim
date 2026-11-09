//! G0/G1 live film transport on the actual coupled thermal/species endpoint.
use super::*;
use fs_thermal::evaporation::GasFilm;

fn transport(surface_t: f64, gas_t: f64, pressure: f64) -> (f64, f64) {
    let film_t = f64::midpoint(surface_t, gas_t);
    // Authored binary-diffusion scaling used as a manufactured transport law.
    // Source-backed correlations are supplied by fs-airflow, not this fixture.
    let coefficient = 0.03 * (film_t / 300.0).powf(1.75) * (100_000.0 / pressure);
    (film_t, coefficient)
}

fn live_films(
    thermal: &DepositionSimulation<'_>,
    gas: &ClosedGasParcel,
    films: &mut BTreeMap<CellId, GasFilm>,
) -> Result<(), ImplicitGasError> {
    for (&id, film) in films {
        let (t, k) = transport(
            thermal.temperature_k(id)?,
            gas.state().state().temperature,
            gas.state().state().pressure,
        );
        film.temperature_k = t;
        film.mass_transfer_coefficient_m_s = k;
    }
    Ok(())
}

#[test]
fn unchanged_films_preserve_the_existing_endpoint_and_work_exactly() {
    let source = catalog(-1e7);
    let gas = gas(&source, 310.0, 1000.0, 0.01);
    let initial = thermal(vec![cell(1, 300.0, 0.0, 1e-3)], None);
    let binding = binding(-1e7);
    let surfaces = [surface(1, &binding, 0.02, 10.0)];
    let old = solve_deposition_with_gas(
        &initial,
        &surfaces,
        &gas,
        &source,
        0.1,
        controls(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let calls = Cell::new(0usize);
    let next = solve_deposition_with_gas_and_films(
        &initial,
        &surfaces,
        &gas,
        &source,
        0.1,
        controls(),
        empty_heat,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(())
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(next.thermal().cells(), old.thermal().cells());
    assert_eq!(next.gas(), old.gas());
    assert_eq!(next.report(), old.report());
    // Initial source admission is followed by each actual residual, including
    // both fresh final checks. Cached residuals do not rerun physical sources.
    assert_eq!(calls.get(), next.report().residual_evaluations + 1);
}

#[test]
fn live_films_close_both_signed_rates_against_the_returned_gas_and_liquid() {
    let source = catalog(-1e7);
    let gas = gas(&source, 310.0, 1500.0, 0.01);
    let saved_gas = gas.clone();
    let initial = thermal(
        vec![cell(1, 305.0, 2.0, 1e-3), cell(2, 290.0, 0.0, 1e-3)],
        Some((1e-4, 0.2)),
    );
    let saved = initial.clone();
    let binding = binding(-1e7);
    // Zero authored coefficient proves that every actual rate uses the live
    // callback, including the initial reference and final gas reconstruction.
    let surfaces = [
        surface(1, &binding, 0.0, 10.0),
        surface(2, &binding, 0.0, 10.0),
    ];
    let target = 0.1;
    let result = solve_deposition_with_gas_and_films(
        &initial,
        &surfaces,
        &gas,
        &source,
        target,
        controls(),
        empty_heat,
        live_films,
        || Ok(()),
    )
    .unwrap();
    let final_gas = result.gas();
    let p = final_gas.state().state().pressure;
    let pv = final_gas.species_partial_pressure_pa(&vapor()).unwrap();
    for (id, sign) in [(CellId(1), 1.0), (CellId(2), -1.0)] {
        let t = result.thermal().temperature_k(id).unwrap();
        let (tf, k) = transport(t, final_gas.state().state().temperature, p);
        let expected = target * 9e-4 * k * unit_film_rate(t, tf, p, pv);
        let phase = &result.report().thermal.thermal.surfaces[&id];
        assert!(sign * phase.transferred_mass_kg > 0.0);
        near(phase.transferred_mass_kg, expected, 1e-14);
        near(
            initial.cells()[&id].mass_kg - result.thermal().cells()[&id].mass_kg,
            expected,
            1e-14,
        );
    }
    near(
        result.report().energy_residual_j,
        0.0,
        controls().phase.gas.energy_tolerance_j,
    );
    near(
        result.report().gas_mass_residual_kg,
        0.0,
        controls().phase.gas.mass_tolerance_kg,
    );
    near(
        final_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-14,
    );
    unchanged(&initial, &saved, &gas, &saved_gas);
}

#[test]
fn film_sources_cannot_replace_gas_pressures_or_the_active_surface_roster() {
    let source = catalog(0.0);
    let gas = gas(&source, 310.0, 1000.0, 0.01);
    let saved_gas = gas.clone();
    let initial = thermal(vec![cell(1, 300.0, 0.0, 1e-3)], None);
    let saved = initial.clone();
    let binding = binding(0.0);
    let surfaces = [surface(1, &binding, 0.02, 0.0)];
    for change in 0..5 {
        let result = solve_deposition_with_gas_and_films(
            &initial,
            &surfaces,
            &gas,
            &source,
            0.1,
            controls(),
            empty_heat,
            |_, _, films| {
                match change {
                    0 => films.get_mut(&CellId(1)).unwrap().total_pressure_pa += 1.0,
                    1 => films.get_mut(&CellId(1)).unwrap().vapor_partial_pressure_pa += 1.0,
                    2 => {
                        films.remove(&CellId(1));
                    }
                    3 => {
                        films.insert(CellId(2), films[&CellId(1)]);
                    }
                    _ => {
                        films
                            .get_mut(&CellId(1))
                            .unwrap()
                            .mass_transfer_coefficient_m_s = -1.0
                    }
                }
                Ok(())
            },
            || Ok(()),
        );
        assert!(result.is_err());
        unchanged(&initial, &saved, &gas, &saved_gas);
    }
}

#[test]
fn a_film_source_failure_in_the_reconstructed_gas_recheck_discards_both_states() {
    let source = catalog(0.0);
    let gas = gas(&source, 310.0, 1000.0, 0.01);
    let saved_gas = gas.clone();
    let initial = thermal(vec![cell(1, 300.0, 0.0, 1e-3)], None);
    let saved = initial.clone();
    let binding = binding(0.0);
    let surfaces = [surface(1, &binding, 0.02, 0.0)];
    let baseline = solve_deposition_with_gas(
        &initial,
        &surfaces,
        &gas,
        &source,
        0.1,
        controls(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let final_call = baseline.report().residual_evaluations + 1;
    let mut calls = 0;
    let result = solve_deposition_with_gas_and_films(
        &initial,
        &surfaces,
        &gas,
        &source,
        0.1,
        controls(),
        empty_heat,
        |_, _, _| {
            calls += 1;
            if calls == final_call {
                Err(ImplicitGasError::Invalid(
                    "original film source refused final gas state",
                ))
            } else {
                Ok(())
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Invalid(
            "original film source refused final gas state"
        ))
    ));
    assert_eq!(calls, final_call);
    unchanged(&initial, &saved, &gas, &saved_gas);
}

#[test]
fn adaptive_and_event_timelines_re_resolve_transport_only_for_active_liquid() {
    use crate::evaporation::gas::implicit::adaptive::{
        AdaptiveImplicitGasControls, advance_deposition_with_gas_adaptive_and_films,
    };
    use crate::evaporation::gas::implicit::timeline::advance_deposition_with_gas_timeline_and_films;

    let source = catalog(0.0);
    let gas = gas(&source, 300.0, 1000.0, 0.01);
    let binding = binding(0.0);
    let surfaces = [surface(1, &binding, 0.02, 0.0)];
    let mut policy = AdaptiveImplicitGasControls::new(controls());
    policy.initial_dt_s = 0.1;
    policy.maximum_dt_s = 0.1;
    policy.solve.thermal.maximum_iterations = 1000;
    policy.solve.thermal.maximum_evaluations = 10000;
    let original = thermal(vec![cell(1, 300.0, 0.0, 1e-3)], None);
    let calls = Cell::new(0usize);
    let seal = |thermal: &DepositionSimulation<'_>,
                _: &ClosedGasParcel,
                films: &mut BTreeMap<CellId, GasFilm>| {
        calls.set(calls.get() + 1);
        assert_eq!(films.len(), 1);
        assert!(thermal.cells().contains_key(&CellId(1)));
        films
            .get_mut(&CellId(1))
            .unwrap()
            .mass_transfer_coefficient_m_s = 0.0;
        Ok(())
    };
    let adaptive = advance_deposition_with_gas_adaptive_and_films(
        &original,
        &surfaces,
        &gas,
        &source,
        0.1,
        policy,
        empty_heat,
        seal,
        || Ok(()),
    )
    .unwrap();
    assert!(calls.get() > 3);
    assert_eq!(adaptive.gas(), &gas);
    assert_eq!(adaptive.thermal().cells(), original.cells());
    assert_eq!(
        adaptive.report().surface_totals[&CellId(1)].transferred_mass_kg,
        0.0
    );

    let future = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.05,
            event: Event::Birth(cell(1, 300.0, 0.0, 1e-3)),
        }],
        &[(MaterialId(1), &ORIGINAL_LIQUID)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 32,
        },
    )
    .unwrap();
    let calls = Cell::new(0usize);
    let timeline = advance_deposition_with_gas_timeline_and_films(
        &future,
        &surfaces,
        &gas,
        &source,
        0.15,
        policy,
        empty_heat,
        |thermal, _, films| {
            calls.set(calls.get() + 1);
            assert!(thermal.time_s() >= 0.05);
            assert_eq!(films.len(), 1);
            assert!(thermal.cells().contains_key(&CellId(1)));
            films
                .get_mut(&CellId(1))
                .unwrap()
                .mass_transfer_coefficient_m_s = 0.0;
            Ok(())
        },
        || Ok(()),
    )
    .unwrap();
    assert!(calls.get() > 3);
    assert_eq!(timeline.report().events, 1);
    assert_eq!(timeline.thermal().time_s(), 0.15);
    assert_eq!(timeline.gas(), &gas);
    assert_eq!(
        timeline.report().surface_totals[&CellId(1)].transferred_mass_kg,
        0.0
    );
    assert!(future.cells().is_empty());
    assert_eq!(future.time_s(), 0.0);
}
