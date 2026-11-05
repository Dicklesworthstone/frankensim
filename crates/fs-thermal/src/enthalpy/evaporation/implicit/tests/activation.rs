//! Selected liquid inventories appear on the original deposition event clock.
use super::*;
use std::cell::Cell;

fn binding(id: u64) -> EvaporatingCellBinding<'static, Interface> {
    EvaporatingCellBinding {
        id: CellId(id),
        law: &CONSTANT,
        film: local_film(CellId(id)),
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    }
}

fn local_film(id: CellId) -> GasFilm {
    film(if id == CellId(3) { 6000.0 } else { 0.0 })
}

fn controls() -> AdaptiveEvaporationControls {
    let mut controls =
        AdaptiveEvaporationControls::new(EvaporationNetworkImplicitControls::default());
    controls.solve.maximum_iterations = 1000;
    controls.solve.maximum_evaluations = 20_000;
    controls.minimum_dt_s = 1e-4;
    controls.mass_absolute_tolerance_kg = 1e-12;
    controls.relative_tolerance = 0.0;
    controls
}

// Balance each actual interface's latent power, so its exact temperature stays
// 300 K while its signed mass follows its exposed area and activation time.
fn balanced(
    _: &DepositionSimulation<'_>,
    liquids: &BTreeMap<CellId, LiquidSurfaceState>,
) -> Result<BoundaryExchangeStepCoefficients<Vec<(CellId, f64)>>, ThermalEvaporationError> {
    let mut result = coefficients(
        liquids
            .iter()
            .map(|(&id, state)| (id, state.area_m2))
            .collect(),
    );
    for (&id, state) in liquids {
        result.transport.heat.bulk.sources_w.insert(
            id,
            transfer(&CONSTANT, 300.0, local_film(id), state.area_m2, 1.0) * LATENT,
        );
    }
    Ok(result)
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
fn existing_solid_heats_before_liquid_birth_and_endpoint_birth_has_no_prior_flux() {
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.5, 300.0, 0.0, 40.0)],
        None,
        vec![
            TimedEvent {
                time_s: 1.0,
                event: Event::Birth(cell(2, 0.1, 300.0, 0.0, 0.0)),
            },
            TimedEvent {
                time_s: 3.0,
                event: Event::Birth(cell(3, 0.2, 300.0, 0.0, 0.0)),
            },
            TimedEvent {
                time_s: 4.0,
                event: Event::Birth(cell(4, 0.3, 290.0, 0.0, 0.0)),
            },
        ],
    );
    let before = old.clone();
    // This last birth is outside the requested window and its pure-liquid
    // admission must wait until it actually arrives.
    let future = EvaporatingCellBinding {
        pure_liquid_enthalpy_range_j_kg: [0.0, CP * 30.0],
        ..binding(4)
    };
    let selected = [future, binding(3), binding(2)];
    let empty: [EvaporatingCellBinding<'_, Interface>; 0] = [];
    let (prefix, prefix_payload, prefix_audit) = old
        .advance_adaptive_evaporating_network_to(&selected, 0.5, controls(), balanced, || Ok(()))
        .unwrap();
    let (heat_only, heat_payload, heat_audit) = old
        .advance_adaptive_evaporating_network_to(&empty, 0.5, controls(), balanced, || Ok(()))
        .unwrap();
    assert_eq!(prefix.cells(), heat_only.cells());
    assert_eq!(prefix_payload, heat_payload);
    assert_eq!(prefix_audit, heat_audit);
    assert!(prefix_audit.surface_totals.is_empty());

    let mut saw_heat_only = false;
    let mut saw_liquid = false;
    let (next, payload, audit) = old
        .advance_adaptive_evaporating_network_to(
            &selected,
            3.0,
            controls(),
            |thermal, liquids| {
                if thermal.time_s() <= 1.0 {
                    assert!(liquids.is_empty());
                    assert_eq!(thermal.cells().len(), 1);
                    saw_heat_only = true;
                } else {
                    assert_eq!(liquids.keys().copied().collect::<Vec<_>>(), [CellId(2)]);
                    assert_eq!(thermal.cells().len(), 2);
                    saw_liquid = true;
                }
                balanced(thermal, liquids)
            },
            || Ok(()),
        )
        .unwrap();
    assert!(saw_heat_only && saw_liquid);
    assert_eq!(payload.unwrap(), [(CellId(2), 0.01)]);
    close(next.temperature_k(CellId(1)).unwrap(), 300.06, 1e-11);
    close(next.temperature_k(CellId(2)).unwrap(), 300.0, 1e-11);
    let dm = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 2.0);
    close(next.cells()[&CellId(2)].mass_kg, 0.1 - dm, 1e-15);
    assert_eq!(next.cells()[&CellId(3)].mass_kg, 0.2);
    assert!(!next.cells().contains_key(&CellId(4)));
    assert_eq!(
        audit.surface_totals.keys().copied().collect::<Vec<_>>(),
        [CellId(2)]
    );
    assert_eq!(audit.surface_totals[&CellId(2)].intervals, 4);
    assert_eq!(audit.last_surfaces.len(), 1);
    assert_eq!(audit.physical.events, 2);
    assert_eq!(audit.physical.substeps, 6);
    close(audit.physical.incoming_mass_kg, 0.3, 1e-15);
    close(audit.physical.transferred_mass_kg, dm, 1e-15);
    close(audit.physical.external_heat_j, 120.0 + dm * LATENT, 1e-7);
    close(audit.physical.energy_residual_j, 0.0, 1e-7);
    close(audit.physical.mass_balance_residual_kg, 0.0, 1e-14);
    unchanged(&old, &before);
}

#[test]
fn successive_births_use_event_exposure_and_keep_opposing_surface_accounts() {
    let mut final_boundary = cell(2, 0.1, 300.0, 0.0, 0.0).boundary;
    final_boundary.reservoir = Some(Reservoir {
        area_m2: 0.008,
        temperature_k: 300.0,
        conductance_w_k: 0.0,
    });
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.5, 300.0, 0.0, 0.0)],
        None,
        vec![
            TimedEvent {
                time_s: 1.0,
                event: Event::Birth(cell(2, 0.1, 300.0, 0.0, 0.0)),
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
                event: Event::Birth(cell(3, 0.2, 300.0, 0.0, 0.0)),
            },
            TimedEvent {
                time_s: 2.0,
                event: Event::Contact(ThermalLink {
                    id: 10,
                    a: CellId(2),
                    b: CellId(3),
                    area_m2: 0.002,
                    conductance_w_k: 0.0,
                }),
            },
            TimedEvent {
                time_s: 3.0,
                event: Event::Separation(9),
            },
            TimedEvent {
                time_s: 4.0,
                event: Event::Boundary(CellId(2), final_boundary),
            },
        ],
    );
    let mut seen = [false; 4];
    let (next, payload, audit) = old
        .advance_adaptive_evaporating_network_to(
            &[binding(3), binding(2)],
            4.0,
            controls(),
            |thermal, liquids| {
                let t = thermal.time_s();
                let (part, expected): (_, &[(CellId, f64)]) = if t <= 1.0 {
                    (0, &[])
                } else if t <= 2.0 {
                    (1, &[(CellId(2), 0.006)])
                } else if t <= 3.0 {
                    (2, &[(CellId(2), 0.004), (CellId(3), 0.008)])
                } else {
                    (3, &[(CellId(2), 0.008), (CellId(3), 0.008)])
                };
                seen[part] = true;
                assert_eq!(liquids.len(), expected.len());
                for &(id, area) in expected {
                    close(liquids[&id].area_m2, area, 1e-17);
                }
                balanced(thermal, liquids)
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(seen, [true; 4]);
    assert_eq!(audit.physical.events, 6);
    assert_eq!(audit.physical.substeps, 8);
    for (id, exposure_time, mass, intervals) in [
        (CellId(2), 0.006 + 0.004 + 0.008, 0.1, 6),
        (CellId(3), 0.008 + 0.008, 0.2, 4),
    ] {
        let dm = transfer(&CONSTANT, 300.0, local_film(id), exposure_time, 1.0);
        close(next.cells()[&id].mass_kg, mass - dm, 1e-15);
        close(next.temperature_k(id).unwrap(), 300.0, 1e-11);
        close(audit.surface_totals[&id].transferred_mass_kg, dm, 1e-15);
        close(audit.surface_totals[&id].latent_energy_j, dm * LATENT, 1e-7);
        assert_eq!(audit.surface_totals[&id].intervals, intervals);
    }
    assert!(audit.surface_totals[&CellId(2)].transferred_mass_kg > 0.0);
    assert!(audit.surface_totals[&CellId(3)].transferred_mass_kg < 0.0);
    assert_eq!(payload.unwrap(), [(CellId(2), 0.008), (CellId(3), 0.008)]);
    close(next.exposed_area_m2(CellId(2)).unwrap(), 0.0, 1e-17);
    close(audit.last_surfaces[&CellId(2)].liquid.area_m2, 0.008, 1e-17);
    close(audit.physical.energy_residual_j, 0.0, 1e-7);
    close(audit.physical.mass_balance_residual_kg, 0.0, 1e-14);
}

#[test]
fn unscheduled_ids_refuse_before_callbacks_and_invalid_liquid_birth_is_atomic() {
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.5, 300.0, 0.0, 40.0)],
        None,
        vec![TimedEvent {
            time_s: 1.0,
            event: Event::Birth(cell(2, 0.1, 290.0, 0.0, 0.0)),
        }],
    );
    let before = old.clone();
    let mut calls = 0;
    let result = old.advance_adaptive_evaporating_network_to(
        &[binding(2), binding(3)],
        2.0,
        controls(),
        |thermal, liquids| {
            calls += 1;
            balanced(thermal, liquids)
        },
        || Ok(()),
    );
    assert!(matches!(result, Err(AdaptiveEvaporationError::Invalid(_))));
    assert_eq!(calls, 0);
    unchanged(&old, &before);

    let liquid = EvaporatingCellBinding {
        pure_liquid_enthalpy_range_j_kg: [0.0, 30.0 * CP],
        ..binding(2)
    };
    let result = old.advance_adaptive_evaporating_network_to(
        &[liquid],
        2.0,
        controls(),
        |thermal, liquids| {
            calls += 1;
            assert!(liquids.is_empty());
            assert!(thermal.time_s() <= 1.0);
            balanced(thermal, liquids)
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(AdaptiveEvaporationError::Implicit(
            EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                DepositionError::MaterialDomain
            ))
        ))
    ));
    assert!(calls > 0);
    unchanged(&old, &before);
}

#[test]
fn source_failure_and_cancellation_after_activation_discard_the_whole_candidate() {
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.5, 300.0, 0.0, 40.0)],
        None,
        vec![TimedEvent {
            time_s: 1.0,
            event: Event::Birth(cell(2, 0.1, 300.0, 0.0, 0.0)),
        }],
    );
    let before = old.clone();
    let mut reached_birth = false;
    let result = old.advance_adaptive_evaporating_network_to(
        &[binding(2)],
        2.0,
        controls(),
        |thermal, liquids| {
            if !liquids.is_empty() {
                reached_birth = true;
                return Err(DepositionError::InvalidInput("original activated heat source").into());
            }
            balanced(thermal, liquids)
        },
        || Ok(()),
    );
    assert!(reached_birth);
    assert!(matches!(
        result,
        Err(AdaptiveEvaporationError::Implicit(
            EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                DepositionError::InvalidInput("original activated heat source")
            ))
        ))
    ));
    unchanged(&old, &before);

    let cancel = Cell::new(false);
    let result = old.advance_adaptive_evaporating_network_to(
        &[binding(2)],
        2.0,
        controls(),
        |thermal, liquids| {
            if !liquids.is_empty() {
                cancel.set(true);
            }
            balanced(thermal, liquids)
        },
        || {
            if cancel.get() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(cancel.get());
    assert!(matches!(
        result,
        Err(AdaptiveEvaporationError::Implicit(
            EvaporationImplicitError::Physical(ThermalEvaporationError::Thermal(
                DepositionError::Cancelled
            ))
        ))
    ));
    unchanged(&old, &before);
}
