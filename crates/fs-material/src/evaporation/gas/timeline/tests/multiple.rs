//! Shared-vapor films use one original heat clock and preserve actual inventories.
use super::super::multiple::advance_deposition_with_gas as advance_multiple;
use super::*;
use crate::evaporation::gas::chamber::GasLiquidChamberSolve;

fn film_birth(id: u64, time: f64, temperature: f64, power: f64) -> TimedEvent {
    let mut event = birth(time, power);
    if let Event::Birth(cell) = &mut event.event {
        cell.id = CellId(id);
        cell.specific_h_j_kg = h(temperature);
    }
    event
}

fn two_films(events: Vec<TimedEvent>, steps: usize) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &LIQUID)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 4,
            max_substeps: steps,
        },
    )
    .unwrap()
}

fn selected<'a>(id: u64, binding: &'a GasLiquidBinding) -> GasEvaporationTimelineSurface<'a> {
    GasEvaporationTimelineSurface {
        cell_id: CellId(id),
        ..surface(binding)
    }
}

fn joint_controls() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 2,
    }
}

#[test]
fn cold_and_warm_films_share_gas_heat_clock_and_datum_independent_balance() {
    let mut endpoints = Vec::new();
    for reference in [0.0, 1e8] {
        let catalog = catalog(reference);
        let humid = catalog
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![0.985, 0.015]).unwrap(),
            )
            .unwrap();
        let initial_gas = ClosedGasParcel::new(humid, point(300.0, 100000.0), 1e-3).unwrap();
        let original = two_films(
            vec![
                film_birth(1, 0.0, 280.0, 0.5),
                film_birth(2, 0.0, 320.0, 1.5),
            ],
            1000,
        )
        .advance_to(0.0, heat(0.01), || Ok(()))
        .unwrap()
        .0;
        let binding = GasLiquidBinding {
            species: id(),
            liquid_enthalpy_offset_j_kg: reference,
        };
        let surfaces = [selected(2, &binding), selected(1, &binding)];
        let run = |surfaces: &[GasEvaporationTimelineSurface<'_>]| {
            advance_multiple(
                &original,
                surfaces,
                &initial_gas,
                &catalog,
                0.2,
                heat(0.01),
                joint_controls(),
                || Ok(()),
            )
            .unwrap()
        };
        let result = run(&surfaces);
        let reversed = run(&[surfaces[1], surfaces[0]]);
        assert_eq!(result.gas(), reversed.gas());
        assert_eq!(result.thermal().cells(), reversed.thermal().cells());
        assert_eq!(result.audit(), reversed.audit());
        let audit = result.audit();
        assert!(audit.cells[&CellId(1)].condensed_mass_kg > 0.0);
        assert_eq!(audit.cells[&CellId(1)].evaporated_mass_kg, 0.0);
        assert!(audit.cells[&CellId(2)].evaporated_mass_kg > 0.0);
        assert_eq!(audit.cells[&CellId(2)].condensed_mass_kg, 0.0);
        assert_eq!(result.thermal().time_s(), 0.2);
        near(audit.thermal.external_heat_j, 0.4, 1e-12);
        assert_eq!(audit.phase_steps, 20);
        assert_eq!(audit.cells[&CellId(1)].phase_steps, audit.phase_steps);
        assert_eq!(audit.cells[&CellId(2)].phase_steps, audit.phase_steps);
        near(audit.mass_residual_kg, 0.0, 1e-14);
        assert!(audit.maximum_species_mass_residual_kg <= 1e-14);
        assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
        near(
            audit.thermal_energy_change_j + audit.gas_energy_change_j,
            0.4,
            1e-6,
        );
        let liquid_change = (result.thermal().cells()[&CellId(1)].mass_kg
            - original.cells()[&CellId(1)].mass_kg)
            + (result.thermal().cells()[&CellId(2)].mass_kg - original.cells()[&CellId(2)].mass_kg);
        near(
            liquid_change + result.gas().species_mass_kg(&id()).unwrap()
                - initial_gas.species_mass_kg(&id()).unwrap(),
            0.0,
            1e-14,
        );
        if reference == 0.0 {
            near(
                result.thermal().energy_j() - original.energy_j()
                    + result.gas().mass_kg() * result.gas().state().specific_internal_energy_j_kg()
                    - initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg(),
                0.4,
                1e-6,
            );
        }
        endpoints.push((
            result.gas().state().state().temperature,
            result.thermal().temperature_k(CellId(1)).unwrap(),
            result.thermal().temperature_k(CellId(2)).unwrap(),
            audit.transferred_mass_kg,
        ));
    }
    near(endpoints[0].0, endpoints[1].0, 1e-8);
    near(endpoints[0].1, endpoints[1].1, 1e-10);
    near(endpoints[0].2, endpoints[1].2, 1e-10);
    near(endpoints[0].3, endpoints[1].3, 1e-14);
}

#[test]
fn distinct_births_and_covering_keep_shared_interval_exposure_and_birth_energy() {
    let catalog = catalog(1e8);
    let initial_gas = gas(&catalog, 1e-3);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 1e8,
    };
    let mut covered = boundary(2.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 1e-3,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let original = two_films(
        vec![
            film_birth(1, 0.1, 300.0, 1.0),
            film_birth(2, 0.2, 300.0, 2.0),
            TimedEvent {
                time_s: 0.3,
                event: Event::Boundary(CellId(2), covered),
            },
        ],
        1000,
    );
    let surfaces = [selected(1, &binding), selected(2, &binding)];
    let born = advance_multiple(
        &original,
        &surfaces,
        &initial_gas,
        &catalog,
        0.1,
        heat(0.1),
        joint_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(born.audit().phase_steps, 0);
    assert_eq!(born.gas(), &initial_gas);
    let result = advance_multiple(
        &original,
        &surfaces,
        &initial_gas,
        &catalog,
        0.4,
        heat(0.1),
        joint_controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = result.audit();
    assert_eq!(audit.thermal.events, 3);
    assert!(audit.cells[&CellId(2)].transferred_mass_kg > 0.0);
    assert_eq!(result.thermal().exposed_area_m2(CellId(2)), Ok(0.0));
    assert_eq!(audit.last_step.as_ref().unwrap().duration_s, 0.4 - 0.3);
    near(audit.thermal.external_heat_j, 0.7, 1e-12);
    near(audit.incoming_energy_j, 0.02 * (h(300.0) + 1e8), 1e-8);
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.incoming_energy_j,
        0.7,
        1e-6,
    );
    assert!(audit.sum_absolute_energy_residual_j <= 1e-6);
    assert!(original.cells().is_empty());
}

#[test]
fn shared_inventory_retry_and_late_cancellation_do_not_replay_heat_or_publish_inputs() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-5);
    let before_gas = initial_gas.clone();
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let surfaces = [selected(1, &binding), selected(2, &binding)];
    let original = two_films(
        vec![
            film_birth(1, 0.0, 300.0, 1.0),
            film_birth(2, 0.0, 300.0, 1.0),
        ],
        1000,
    );
    let mut polls = 0;
    let result = advance_multiple(
        &original,
        &surfaces,
        &initial_gas,
        &catalog,
        2.0,
        heat(2.0),
        joint_controls(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let audit = result.audit();
    assert!(audit.rejected_intervals > 0);
    assert!(audit.attempted_thermal_substeps >= audit.thermal.substeps + audit.rejected_intervals);
    assert!(audit.attempted_thermal_substeps <= 1000);
    assert_eq!(audit.thermal.events, 2);
    near(audit.thermal.external_heat_j, 4.0, 1e-12);
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.incoming_energy_j,
        4.0,
        1e-6,
    );
    let mut seen = 0;
    assert!(matches!(
        advance_multiple(
            &original,
            &surfaces,
            &initial_gas,
            &catalog,
            2.0,
            heat(2.0),
            joint_controls(),
            || {
                seen += 1;
                if seen == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    let limited = two_films(
        vec![
            film_birth(1, 0.0, 300.0, 1.0),
            film_birth(2, 0.0, 300.0, 1.0),
        ],
        1,
    );
    assert!(matches!(
        advance_multiple(
            &limited,
            &surfaces,
            &initial_gas,
            &catalog,
            2.0,
            heat(2.0),
            joint_controls(),
            || Ok(())
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert!(
        advance_multiple(
            &original,
            &[surfaces[0], surfaces[0]],
            &initial_gas,
            &catalog,
            0.0,
            heat(0.1),
            joint_controls(),
            || Ok(())
        )
        .is_err()
    );
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(initial_gas, before_gas);
}

#[test]
fn singleton_retains_original_endpoint_and_complete_multi_film_clock_refines() {
    let catalog = catalog(0.0);
    let initial_gas = gas(&catalog, 1e-4);
    let binding = GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let original = thermal(vec![birth(0.0, 2.0)], 1000);
    let old = super::advance_deposition_with_gas(
        &original,
        surface(&binding),
        &initial_gas,
        &catalog,
        0.2,
        heat(0.05),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let joint = advance_multiple(
        &original,
        &[surface(&binding)],
        &initial_gas,
        &catalog,
        0.2,
        heat(0.05),
        joint_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(old.gas(), joint.gas());
    assert_eq!(old.thermal().cells(), joint.thermal().cells());
    let original = two_films(
        vec![
            film_birth(1, 0.0, 300.0, 1.0),
            film_birth(2, 0.0, 310.0, 2.0),
        ],
        10000,
    );
    let surfaces = [selected(1, &binding), selected(2, &binding)];
    let run = |dt| {
        advance_multiple(
            &original,
            &surfaces,
            &initial_gas,
            &catalog,
            1.0,
            heat(dt),
            joint_controls(),
            || Ok(()),
        )
        .unwrap()
        .gas()
        .species_partial_pressure_pa(&id())
        .unwrap()
    };
    let reference = run(0.003125);
    let a = (run(0.1) - reference).abs();
    let b = (run(0.05) - reference).abs();
    let c = (run(0.025) - reference).abs();
    assert!(a > 1e-6);
    assert!(b < 0.75 * a, "{a}, {b}");
    assert!(c < 0.75 * b, "{b}, {c}");
}
