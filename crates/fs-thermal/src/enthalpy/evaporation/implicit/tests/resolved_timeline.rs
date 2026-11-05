//! Endpoint film sources on the one conservative adaptive event timeline.
use super::*;
use std::cell::Cell;

fn source(id: u64) -> EvaporatingCellSource<'static, Interface> {
    EvaporatingCellSource {
        id: CellId(id),
        law: &CONSTANT,
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    }
}

fn policy(span: f64, mass_tolerance: f64) -> AdaptiveEvaporationControls {
    let mut controls =
        AdaptiveEvaporationControls::new(EvaporationNetworkImplicitControls::default());
    controls.solve.maximum_iterations = 10_000;
    controls.solve.maximum_evaluations = 100_000;
    controls.initial_dt_s = span;
    controls.maximum_dt_s = span;
    controls.minimum_dt_s = 1e-5;
    controls.mass_absolute_tolerance_kg = mass_tolerance;
    controls.relative_tolerance = 0.0;
    controls
}

fn balance<T>(
    payload: T,
    liquids: &BTreeMap<CellId, LiquidSurfaceState>,
    films: &BTreeMap<CellId, GasFilm>,
) -> BoundaryExchangeStepCoefficients<T> {
    let mut heat = coefficients(payload);
    for (&id, liquid) in liquids {
        heat.transport.heat.bulk.sources_w.insert(
            id,
            transfer(&CONSTANT, 300.0, films[&id], liquid.area_m2, 1.0) * LATENT,
        );
    }
    heat
}

fn unchanged(actual: &DepositionSimulation<'_>, original: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s, original.time_s);
    assert_eq!(actual.cursor, original.cursor);
    assert_eq!(actual.cells, original.cells);
    assert_eq!(actual.links, original.links);
    assert_eq!(actual.covered_areas, original.covered_areas);
    assert_eq!(actual.schedule, original.schedule);
}

#[test]
fn humidity_from_another_cells_temperature_converges_to_the_integrated_mass() {
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, 0.0),
            cell(2, 1.0, 300.0, 0.0, 400.0),
        ],
        None,
        vec![],
    );
    // The second original cell heats at exactly 0.1 K/s. The authored humidity
    // source therefore gives pv(t)=1000+100t Pa, without evolving a gas parcel.
    // Integrate the Stefan logarithm analytically to obtain the exact transfer.
    let duration = 20.0;
    let primitive = |pressure: f64| pressure * ((pressure / 96_000.0).ln() - 1.0);
    let prefactor = 0.01 * 0.01 * 0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0);
    let exact_transfer = prefactor * (primitive(99_000.0) - primitive(97_000.0)) / 100.0;
    let run = |mass_tolerance| {
        let film_calls = Cell::new(0usize);
        let heat_calls = Cell::new(0usize);
        let (next, _, audit) = old
            .advance_adaptive_evaporating_network_with_films_to(
                &[source(1)],
                duration,
                policy(duration, mass_tolerance),
                |thermal| {
                    film_calls.set(film_calls.get() + 1);
                    let humidity = 1000.0 + 1000.0 * (thermal.temperature_k(CellId(2))? - 300.0);
                    Ok(BTreeMap::from([(CellId(1), film(humidity))]))
                },
                |_, liquids, films| {
                    heat_calls.set(heat_calls.get() + 1);
                    Ok(balance((), liquids, films))
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(heat_calls.get(), audit.residual_evaluations);
        assert_eq!(
            film_calls.get(),
            audit.residual_evaluations + 3 * audit.comparisons
        );
        close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-10);
        close(next.temperature_k(CellId(2)).unwrap(), 302.0, 1e-10);
        close(
            audit.physical.external_heat_j,
            400.0 * duration + audit.physical.transferred_mass_kg * LATENT,
            1e-6,
        );
        close(audit.physical.energy_residual_j, 0.0, 1e-6);
        let error = (audit.physical.transferred_mass_kg - exact_transfer).abs();
        (error, audit)
    };
    let (loose_error, loose) = run(1e-5);
    let (tight_error, tight) = run(1e-8);
    assert_eq!(loose.accepted_comparisons, 1);
    assert!(tight.accuracy_rejections > 0);
    assert!(tight.accepted_comparisons > loose.accepted_comparisons);
    assert!(
        tight_error < loose_error * 0.15,
        "{tight_error:e} vs {loose_error:e}"
    );
    assert!(tight_error < 5e-7);
}

#[test]
fn delayed_interfaces_resolve_empty_prefixes_and_use_left_limit_event_exposure() {
    let incoming = cell(2, 0.1, 300.0, 0.0, 0.0);
    let mut final_boundary = incoming.boundary;
    final_boundary.reservoir = Some(Reservoir {
        area_m2: 0.01,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.5, 300.0, 0.0, 40.0)],
        None,
        vec![
            TimedEvent {
                time_s: 1.0,
                event: Event::Birth(incoming),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.004,
                    conductance_w_k: 0.0,
                }),
            },
            TimedEvent {
                time_s: 2.0,
                event: Event::Separation(9),
            },
            TimedEvent {
                time_s: 3.0,
                event: Event::Boundary(CellId(2), final_boundary),
            },
        ],
    );
    let mut empty_initial = 0;
    let mut empty_heat = 0;
    let (next, payload, audit) = old
        .advance_adaptive_evaporating_network_with_films_to(
            &[source(2)],
            3.0,
            policy(1.0, 1e-4),
            |thermal| {
                if thermal.cells().contains_key(&CellId(2)) {
                    let humidity = 1000.0 + 500.0 * (thermal.time_s() - 1.0);
                    Ok(BTreeMap::from([(CellId(2), film(humidity))]))
                } else {
                    empty_initial += 1;
                    assert!(thermal.time_s() <= 1.0);
                    Ok(BTreeMap::new())
                }
            },
            |thermal, liquids, films| {
                assert_eq!(liquids.len(), films.len());
                let payload = if liquids.is_empty() {
                    empty_heat += 1;
                    assert!(thermal.time_s() <= 1.0);
                    None
                } else {
                    let area = if thermal.time_s() <= 2.0 { 0.006 } else { 0.01 };
                    close(liquids[&CellId(2)].area_m2, area, 1e-17);
                    Some((
                        thermal.time_s(),
                        area,
                        films[&CellId(2)].vapor_partial_pressure_pa,
                    ))
                };
                Ok(balance(payload, liquids, films))
            },
            || Ok(()),
        )
        .unwrap();
    assert!(empty_initial > 0 && empty_heat > 0);
    let expected = [(1.5, 0.006), (2.0, 0.006), (2.5, 0.01), (3.0, 0.01)]
        .into_iter()
        .map(|(t, area)| {
            transfer(
                &CONSTANT,
                300.0,
                film(1000.0 + 500.0 * (t - 1.0)),
                area,
                0.5,
            )
        })
        .sum::<f64>();
    close(audit.physical.transferred_mass_kg, expected, 1e-15);
    close(next.cells()[&CellId(2)].mass_kg, 0.1 - expected, 1e-15);
    close(
        audit.physical.external_heat_j,
        120.0 + expected * LATENT,
        1e-7,
    );
    assert_eq!(audit.surface_totals[&CellId(2)].intervals, 4);
    assert_eq!(audit.physical.events, 4);
    assert_eq!(audit.physical.substeps, 6);
    assert_eq!(payload.unwrap().unwrap(), (3.0, 0.01, 2000.0));
    assert_eq!(next.exposed_area_m2(CellId(2)).unwrap(), 0.0);
    assert_eq!(audit.last_surfaces[&CellId(2)].liquid.area_m2, 0.01);
}

#[test]
fn each_fine_interval_resolves_its_own_inventory_history_and_fixed_path_is_unchanged() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let duration = 100.0;
    let rate = transfer(&CONSTANT, 300.0, film(1000.0), 0.01, 1.0);
    let half_mass = 0.1 - 0.5 * duration * rate;
    let fine_mass = half_mass - 0.5 * duration * rate * (half_mass / 0.1);
    let saw_coarse = Cell::new(false);
    let saw_second_fine = Cell::new(false);
    let (next, payload, audit) = old
        .advance_adaptive_evaporating_network_with_films_to(
            &[source(1)],
            duration,
            policy(duration, 1e-3),
            |thermal| {
                let retained_mass = thermal.cells()[&CellId(1)].mass_kg;
                let mut local = film(1000.0);
                local.mass_transfer_coefficient_m_s *= retained_mass / 0.1;
                if thermal.time_s() == duration {
                    if retained_mass == 0.1 {
                        saw_coarse.set(true);
                    } else {
                        close(retained_mass, half_mass, 1e-15);
                        saw_second_fine.set(true);
                    }
                }
                Ok(BTreeMap::from([(CellId(1), local)]))
            },
            |thermal, liquids, films| {
                Ok(balance(
                    (
                        thermal.cells()[&CellId(1)].mass_kg,
                        liquids[&CellId(1)].mass_kg,
                        films[&CellId(1)].mass_transfer_coefficient_m_s,
                    ),
                    liquids,
                    films,
                ))
            },
            || Ok(()),
        )
        .unwrap();
    assert!(saw_coarse.get() && saw_second_fine.get());
    assert_eq!(audit.accepted_comparisons, 1);
    let (retained, proposed, coefficient) = payload.unwrap();
    close(retained, half_mass, 1e-15);
    close(proposed, fine_mass, 1e-15);
    close(coefficient, 0.01 * half_mass / 0.1, 1e-17);
    close(next.cells()[&CellId(1)].mass_kg, fine_mass, 1e-15);
    close(audit.physical.transferred_mass_kg, 0.1 - fine_mass, 1e-15);
    close(
        audit.physical.external_heat_j,
        (0.1 - fine_mass) * LATENT,
        1e-7,
    );
    close(
        audit.last_surfaces[&CellId(1)].transferred_mass_kg,
        half_mass - fine_mass,
        1e-15,
    );
    assert_eq!(audit.surface_totals[&CellId(1)].intervals, 2);
    assert!((fine_mass - (0.1 - duration * rate)).abs() > 1e-8);

    let fixed_calls = Cell::new(0usize);
    let resolved_calls = Cell::new(0usize);
    let constant_films = BTreeMap::from([(CellId(1), film(1000.0))]);
    let fixed_binding = EvaporatingCellBinding {
        id: CellId(1),
        law: &CONSTANT,
        film: film(1000.0),
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    };
    let (fixed, fixed_payload, fixed_audit) = old
        .advance_adaptive_evaporating_network_to(
            &[fixed_binding],
            duration,
            policy(duration, 1e-3),
            |_, liquids| {
                fixed_calls.set(fixed_calls.get() + 1);
                Ok(balance((), liquids, &constant_films))
            },
            || Ok(()),
        )
        .unwrap();
    let (resolved, resolved_payload, resolved_audit) = old
        .advance_adaptive_evaporating_network_with_films_to(
            &[source(1)],
            duration,
            policy(duration, 1e-3),
            |_| Ok(constant_films.clone()),
            |_, liquids, films| {
                resolved_calls.set(resolved_calls.get() + 1);
                Ok(balance((), liquids, films))
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(fixed.cells(), resolved.cells());
    assert_eq!(fixed_payload, resolved_payload);
    assert_eq!(fixed_audit, resolved_audit);
    assert_eq!(fixed_calls.get(), resolved_calls.get());
    assert_eq!(fixed_calls.get(), fixed_audit.residual_evaluations);
}

#[test]
fn film_source_errors_and_cancellation_after_accuracy_rejection_are_atomic() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let before = old.clone();
    for source_error in [true, false] {
        let saw_second_half = Cell::new(false);
        let reached_retry = Cell::new(false);
        let cancel = Cell::new(false);
        let result = old.advance_adaptive_evaporating_network_with_films_to(
            &[source(1)],
            20.0,
            policy(20.0, 1e-8),
            |thermal| {
                let t = thermal.time_s();
                if t == 10.0 && thermal.cells()[&CellId(1)].mass_kg < 0.1 {
                    saw_second_half.set(true);
                }
                // First comparison has completed both fine intervals. Its
                // time-varying humidity fails the mass-accuracy allowance, so
                // the next source evaluation restarts at the unchanged origin.
                if t == 0.0 && saw_second_half.get() {
                    reached_retry.set(true);
                    if source_error {
                        return Err(
                            DepositionError::InvalidInput("humidity source unavailable").into()
                        );
                    }
                    cancel.set(true);
                }
                Ok(BTreeMap::from([(CellId(1), film(1000.0 + 100.0 * t))]))
            },
            |_, liquids, films| Ok(balance((), liquids, films)),
            || {
                if cancel.get() {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(reached_retry.get());
        if source_error {
            assert!(matches!(
                result,
                Err(AdaptiveEvaporationError::Implicit(
                    EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                        DepositionError::InvalidInput("humidity source unavailable")
                    ))
                ))
            ));
        } else {
            assert!(matches!(
                result,
                Err(AdaptiveEvaporationError::Implicit(
                    EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                        DepositionError::Cancelled
                    ))
                ))
            ));
        }
        unchanged(&old, &before);
    }
}
