//! G0 live-source admission/atomicity and G1 independently closed film feedback.
use super::*;
use crate::enthalpy::evaporation::implicit::{
    EvaporatingCellBinding, EvaporatingCellSource, EvaporatingNetworkBackwardEulerStep,
    EvaporationImplicitError, EvaporationNetworkImplicitControls,
};
use std::cell::Cell;

fn source(id: u64, law: &dyn EvaporationLaw) -> EvaporatingCellSource<'_> {
    EvaporatingCellSource {
        id: CellId(id),
        law,
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    }
}

fn one_film(value: GasFilm) -> BTreeMap<CellId, GasFilm> {
    BTreeMap::from([(CellId(1), value)])
}

fn empty_heat(
    _: &DepositionSimulation<'_>,
    _: &BTreeMap<CellId, LiquidSurfaceState>,
    _: &BTreeMap<CellId, GasFilm>,
) -> Result<BoundaryExchangeStepCoefficients<()>, ThermalEvaporationError> {
    Ok(coefficients(()))
}

fn unchanged(actual: &DepositionSimulation<'_>, expected: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s(), expected.time_s());
    assert_eq!(actual.cursor, expected.cursor);
    assert_eq!(actual.cells(), expected.cells());
    assert_eq!(actual.links(), expected.links());
    assert_eq!(actual.covered_areas, expected.covered_areas);
    assert_eq!(actual.schedule, expected.schedule);
    assert_eq!(actual.budget, expected.budget);
}

// Deliberately analytic environmental feedback, not a finite gas conservation
// model or an empirical material/film correlation. Every film input depends on
// the thermal endpoint, including the non-evaporating neighbor's temperature.
fn feedback_film(a: f64, b: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 95_000.0 + 400.0 * (a - 300.0) + 250.0 * (b - 300.0),
        vapor_partial_pressure_pa: 1500.0 + 50.0 * (a - 300.0) + 100.0 * (b - 300.0),
        temperature_k: 295.0 + 0.6 * (a - 300.0) + 0.4 * (b - 300.0),
        mass_transfer_coefficient_m_s: 0.008 * (1.0 + 0.02 * (a - 300.0) + 0.01 * (b - 300.0)),
    }
}

fn feedback_films(
    thermal: &DepositionSimulation<'_>,
) -> Result<BTreeMap<CellId, GasFilm>, ThermalEvaporationError> {
    assert_eq!(thermal.cells()[&CellId(1)].mass_kg, 0.1);
    Ok(one_film(feedback_film(
        thermal.temperature_k(CellId(1))?,
        thermal.temperature_k(CellId(2))?,
    )))
}

fn feedback_heat(
    thermal: &DepositionSimulation<'_>,
    liquids: &BTreeMap<CellId, LiquidSurfaceState>,
    films: &BTreeMap<CellId, GasFilm>,
) -> Result<BoundaryExchangeStepCoefficients<(f64, f64, f64, GasFilm)>, ThermalEvaporationError> {
    let liquid = liquids[&CellId(1)];
    let mut result = coefficients((
        thermal.temperature_k(CellId(1))?,
        thermal.temperature_k(CellId(2))?,
        liquid.mass_kg,
        films[&CellId(1)],
    ));
    result
        .transport
        .heat
        .bulk
        .sources_w
        .insert(CellId(1), 60.0 * (liquid.mass_kg - 0.1));
    Ok(result)
}

#[test]
fn endpoint_films_and_mass_dependent_heat_close_the_same_original_network() {
    let dt = 500.0;
    let target = [307.0, 304.0];
    let endpoint_film = feedback_film(target[0], target[1]);
    let raw_dm = transfer(&VARYING, target[0], endpoint_film, 0.008, dt);
    let mass = 0.1 - raw_dm;
    let dm = 0.1 - mass;
    let latent = LATENT - 2000.0 * (target[0] - 300.0);
    let required = [0.1 * CP * 7.0 + dm * latent, -0.2 * CP * 6.0];
    let contact_w = 0.4 * (target[1] - target[0]);
    let authored = [
        required[0] / dt - contact_w - 60.0 * (mass - 0.1),
        required[1] / dt + contact_w,
    ];
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, authored[0]),
            cell(2, 0.2, 310.0, 0.0, authored[1]),
        ],
        Some(0.4),
        vec![],
    );
    let snapshot = old.clone();
    let initial_film = feedback_film(300.0, 310.0);
    assert_ne!(initial_film, endpoint_film);
    assert!(dm > 0.0);
    let binding = EvaporatingCellBinding {
        id: CellId(1),
        law: &VARYING,
        film: initial_film,
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    };
    let step = EvaporatingNetworkBackwardEulerStep::new(&old, &[binding], dt, || Ok(())).unwrap();
    // The energy coordinates are derived from the chosen temperatures; neither
    // the low-level residual nor Newton supplies the reference mass/heat balance.
    let trial = step
        .evaluate_with_films(
            &[0.1 * CP * 7.0, -0.2 * CP * 6.0],
            feedback_films,
            feedback_heat,
            || Ok(()),
        )
        .unwrap();
    for &residual in trial.residual_j() {
        close(residual, 0.0, 1e-8);
    }
    close(
        trial.audit().surfaces[&CellId(1)].transferred_mass_kg,
        dm,
        1e-15,
    );
    assert_eq!(trial.payload().3, endpoint_film);

    let (next, payload, report) = old
        .solve_implicit_evaporating_network_step_with_films(
            &[source(1, &VARYING)],
            dt,
            EvaporationNetworkImplicitControls::default(),
            feedback_films,
            feedback_heat,
            || Ok(()),
        )
        .unwrap();
    close(payload.0, target[0], 1e-9);
    close(payload.1, target[1], 1e-9);
    close(payload.2, mass, 1e-12);
    assert_eq!(payload.3, feedback_film(payload.0, payload.1));
    let surface = report.physical.surfaces[&CellId(1)];
    close(next.temperature_k(CellId(1)).unwrap(), target[0], 1e-9);
    close(next.temperature_k(CellId(2)).unwrap(), target[1], 1e-9);
    close(next.cells()[&CellId(1)].mass_kg, mass, 1e-12);
    close(surface.liquid.area_m2, 0.008, 1e-17);
    close(surface.transferred_mass_kg, dm, 1e-12);
    close(surface.latent_energy_j, dm * latent, 1e-6);
    close(surface.supplied_heat_j, required[0], 1e-6);
    close(
        report.physical.external_heat_j,
        required[0] + required[1],
        1e-6,
    );
    close(
        next.energy_j() - old.energy_j() + surface.vapor_energy_j,
        required[0] + required[1],
        1e-6,
    );
    assert_eq!(next.links(), old.links());
    assert_eq!(next.cells()[&CellId(2)].mass_kg, 0.2);
    unchanged(&old, &snapshot);
}

#[test]
fn endpoint_humidity_can_reverse_initial_evaporation_into_condensation() {
    let dt = 500.0;
    let raw_dm = transfer(&CONSTANT, 300.0, film(6000.0), 0.01, dt);
    let dm = 0.1 - (0.1 - raw_dm);
    assert!(transfer(&CONSTANT, 300.0, film(1000.0), 0.01, dt) > 0.0);
    assert!(dm < 0.0);
    // Remove precisely the endpoint condensation heat: the mass changes even
    // though the accepted temperature is unchanged. A frozen initial film
    // would evaporate and could not close this independently authored balance.
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, dm * LATENT / dt)],
        None,
        vec![],
    );
    let mut queried_times = Vec::new();
    let (next, endpoint, report) = old
        .solve_implicit_evaporating_network_step_with_films(
            &[source(1, &CONSTANT)],
            dt,
            EvaporationNetworkImplicitControls::default(),
            |thermal| {
                queried_times.push(thermal.time_s());
                Ok(one_film(film(1000.0 + 10.0 * thermal.time_s())))
            },
            |_, _, films| Ok(coefficients(films[&CellId(1)])),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(queried_times[0], 0.0);
    assert!(queried_times[1..].iter().all(|&time| time == dt));
    assert_eq!(endpoint, film(6000.0));
    close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-11);
    close(next.cells()[&CellId(1)].mass_kg, 0.1 - dm, 1e-14);
    let surface = report.physical.surfaces[&CellId(1)];
    assert!(surface.transferred_mass_kg < 0.0);
    assert!(surface.latent_energy_j < 0.0);
    assert!(surface.vapor_energy_j < 0.0);
    close(surface.transferred_mass_kg, dm, 1e-14);
    close(report.physical.external_heat_j, dm * LATENT, 1e-8);
    close(
        next.energy_j() - old.energy_j() + surface.vapor_energy_j,
        dm * LATENT,
        1e-8,
    );
}

#[test]
fn fixed_film_resolver_preserves_exact_multi_surface_results_and_work() {
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, 12.0),
            cell(2, 0.2, 305.0, 0.0, -3.0),
        ],
        Some(0.4),
        vec![],
    );
    let films = BTreeMap::from([(CellId(1), film(0.0)), (CellId(2), film(8000.0))]);
    let sources = [source(2, &VARYING), source(1, &VARYING)];
    let bindings: Vec<_> = sources
        .iter()
        .map(|source| EvaporatingCellBinding {
            id: source.id,
            law: source.law,
            film: films[&source.id],
            pure_liquid_enthalpy_range_j_kg: source.pure_liquid_enthalpy_range_j_kg,
        })
        .collect();
    let controls = EvaporationNetworkImplicitControls::default();
    let (fixed, fixed_payload, fixed_report) = old
        .solve_implicit_evaporating_network_step(
            &bindings,
            250.0,
            controls,
            |_, liquids| Ok(coefficients(liquids.clone())),
            || Ok(()),
        )
        .unwrap();
    let calls = Cell::new(0usize);
    let (resolved, payload, report) = old
        .solve_implicit_evaporating_network_step_with_films(
            &sources,
            250.0,
            controls,
            |_| {
                calls.set(calls.get() + 1);
                Ok(films.clone())
            },
            |_, liquids, current_films| {
                assert_eq!(current_films, &films);
                Ok(coefficients(liquids.clone()))
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(resolved.cells(), fixed.cells());
    assert_eq!(resolved.links(), fixed.links());
    assert_eq!(resolved.time_s(), fixed.time_s());
    assert_eq!(payload, fixed_payload);
    assert_eq!(report, fixed_report);
    assert_eq!(calls.get(), report.residual_evaluations + 1);
    assert!(report.physical.surfaces[&CellId(1)].transferred_mass_kg > 0.0);
    assert!(report.physical.surfaces[&CellId(2)].transferred_mass_kg < 0.0);
}

#[test]
fn every_film_map_and_source_refusal_is_validated_before_heat_or_publication() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 15.0)], None, vec![]);
    let snapshot = old.clone();
    let sources = [source(1, &VARYING)];
    let controls = EvaporationNetworkImplicitControls::default();
    for fail_at_endpoint in [false, true] {
        for invalid in 0..4 {
            let film_calls = Cell::new(0usize);
            let heat_calls = Cell::new(0usize);
            let result = old.solve_implicit_evaporating_network_step_with_films(
                &sources,
                1.0,
                controls,
                |thermal| {
                    film_calls.set(film_calls.get() + 1);
                    if fail_at_endpoint && thermal.time_s() == 0.0 {
                        return Ok(one_film(film(0.0)));
                    }
                    Ok(match invalid {
                        0 => BTreeMap::new(),
                        1 => BTreeMap::from([(CellId(1), film(0.0)), (CellId(2), film(0.0))]),
                        2 => BTreeMap::from([(CellId(2), film(0.0))]),
                        _ => one_film(GasFilm {
                            vapor_partial_pressure_pa: 100_000.0,
                            ..film(0.0)
                        }),
                    })
                },
                |_, _, _| {
                    heat_calls.set(heat_calls.get() + 1);
                    Ok(coefficients(()))
                },
                || Ok(()),
            );
            assert!(matches!(result, Err(EvaporationImplicitError::Physical(_))));
            assert_eq!(film_calls.get(), if fail_at_endpoint { 2 } else { 1 });
            assert_eq!(heat_calls.get(), 0);
            unchanged(&old, &snapshot);
        }
    }
    for refuse_film in [true, false] {
        let heat_calls = Cell::new(0usize);
        let result = old.solve_implicit_evaporating_network_step_with_films(
            &sources,
            1.0,
            controls,
            |thermal| {
                if refuse_film && thermal.time_s() > 0.0 {
                    Err(DepositionError::InvalidInput("live film source refused").into())
                } else {
                    Ok(one_film(film(0.0)))
                }
            },
            |_, _, _| {
                heat_calls.set(heat_calls.get() + 1);
                Err::<BoundaryExchangeStepCoefficients<()>, _>(
                    DepositionError::InvalidInput("live heat source refused").into(),
                )
            },
            || Ok(()),
        );
        let expected = if refuse_film {
            "live film source refused"
        } else {
            "live heat source refused"
        };
        assert!(matches!(
            result,
            Err(EvaporationImplicitError::Physical(
                ThermalEvaporationError::Thermal(DepositionError::InvalidInput(
                    message
                ))
            )) if message == expected
        ));
        assert_eq!(heat_calls.get(), if refuse_film { 0 } else { 1 });
        unchanged(&old, &snapshot);
    }

    let no_sources: [EvaporatingCellSource<'_>; 0] = [];
    let (heat_only, (), report) = old
        .solve_implicit_evaporating_network_step_with_films(
            &no_sources,
            1.0,
            controls,
            |_| Ok(BTreeMap::new()),
            empty_heat,
            || Ok(()),
        )
        .unwrap();
    close(
        heat_only.temperature_k(CellId(1)).unwrap(),
        300.0 + 15.0 / 400.0,
        1e-10,
    );
    assert!(report.physical.surfaces.is_empty());
    assert_eq!(heat_only.mass_kg(), old.mass_kg());
    assert!(matches!(
        old.solve_implicit_evaporating_network_step_with_films(
            &no_sources,
            1.0,
            controls,
            |_| Ok(one_film(film(0.0))),
            empty_heat,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::Physical(_))
    ));
    unchanged(&old, &snapshot);
}

#[test]
fn final_films_and_heat_are_fresh_funded_and_cancellable() {
    let raw_dm = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0);
    let dm = 0.1 - (0.1 - raw_dm);
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, dm * LATENT)],
        None,
        vec![],
    );
    let snapshot = old.clone();
    let sources = [source(1, &CONSTANT)];
    let controls = EvaporationNetworkImplicitControls::default();
    let film_calls = Cell::new(0usize);
    let heat_calls = Cell::new(0usize);
    let mut polls = 0usize;
    let (next, payload, report) = old
        .solve_implicit_evaporating_network_step_with_films(
            &sources,
            1.0,
            controls,
            |_| {
                film_calls.set(film_calls.get() + 1);
                Ok(one_film(film(0.0)))
            },
            |_, _, _| {
                heat_calls.set(heat_calls.get() + 1);
                Ok(coefficients((film_calls.get(), heat_calls.get())))
            },
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    // The exact temperature is already the initial Newton guess. The separate
    // final physical recheck still queries both sources and spends work.
    close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-12);
    assert_eq!(report.newton_iterations, 0);
    assert_eq!(report.residual_evaluations, 2);
    assert_eq!(payload, (3, 2));
    assert_eq!(film_calls.get(), report.residual_evaluations + 1);
    assert_eq!(heat_calls.get(), report.residual_evaluations);

    film_calls.set(0);
    heat_calls.set(0);
    let limited = EvaporationNetworkImplicitControls {
        maximum_evaluations: 1,
        ..controls
    };
    let result = old.solve_implicit_evaporating_network_step_with_films(
        &sources,
        1.0,
        limited,
        |_| {
            film_calls.set(film_calls.get() + 1);
            Ok(one_film(film(0.0)))
        },
        |_, _, _| {
            heat_calls.set(heat_calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::Budget(
                "implicit evaporation residual evaluations"
            ))
        ))
    ));
    assert_eq!(film_calls.get(), 2);
    assert_eq!(heat_calls.get(), 1);

    film_calls.set(0);
    heat_calls.set(0);
    let result = old.solve_implicit_evaporating_network_step_with_films(
        &sources,
        1.0,
        controls,
        |_| {
            film_calls.set(film_calls.get() + 1);
            if film_calls.get() == 3 {
                Err(DepositionError::InvalidInput("final environmental recheck refused").into())
            } else {
                Ok(one_film(film(0.0)))
            }
        },
        |_, _, _| {
            heat_calls.set(heat_calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::InvalidInput(
                "final environmental recheck refused"
            ))
        ))
    ));
    assert_eq!(film_calls.get(), 3);
    assert_eq!(heat_calls.get(), 1);

    let mut count = 0usize;
    let result = old.solve_implicit_evaporating_network_step_with_films(
        &sources,
        1.0,
        controls,
        |_| Ok(one_film(film(0.0))),
        empty_heat,
        || {
            count += 1;
            if count == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(count, polls);
    unchanged(&old, &snapshot);
}
