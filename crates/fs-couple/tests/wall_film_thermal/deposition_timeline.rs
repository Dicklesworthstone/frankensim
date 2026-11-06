//! Native events, original birth energy and one continuing wall-film inventory.
use super::{deposition::*, fixture::*};
use std::cell::Cell;

use fs_blake3::ContentHash;
use fs_couple::finite_gas_thermal::wall_film::{
    deposition::{
        DepositionWallFilmBinding, advance_deposition_with_wall_film,
        timeline::{
            DepositionWallFilmTimelineControls, advance_deposition_with_wall_film_timeline,
        },
    },
    window::WallFilmWindowSolve,
};
use fs_material::{
    evaporation::gas::{implicit::ImplicitGasError, wall_film::GasWallFilmInventory},
    phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve},
};
use fs_thermal::enthalpy::{
    Boundary, CellId, DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId,
    StepConfig, TimedEvent,
};

fn timeline(step: f64) -> DepositionWallFilmTimelineControls {
    let mut controls = DepositionWallFilmTimelineControls::new(
        policy(),
        WallFilmWindowSolve {
            maximum_step_s: step,
            minimum_step_s: 1e-7,
            maximum_attempts: 2000,
            energy_tolerance_j: 1e-6,
            mass_tolerance_kg: 1e-13,
        },
    );
    controls.step.implicit.thermal.maximum_iterations = 4096;
    controls.step.implicit.thermal.maximum_evaluations = 100_000;
    controls
}
fn plateau(t: f64, datum: f64, identity: u8) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([identity; 32]),
        [
            (280.0, datum + (280.0 - t) * CP_WALL, 0.0),
            (t, datum, 0.0),
            (t, datum + 100_000.0, 1.0),
            (340.0, datum + 100_000.0 + (340.0 - t) * CP_WALL, 1.0),
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
    .unwrap()
}

#[test]
fn native_contact_and_boundary_events_match_manual_endpoint_composition_and_heat_clock() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let mut wall = cell(1, 1, 0.0, AREA + CONTACT_AREA);
    wall.boundary.htc_w_m2_k = 10.0;
    wall.boundary.source_w = 0.5;
    let changed = Boundary {
        ambient_k: 285.0,
        htc_w_m2_k: 20.0,
        source_w: 2.0,
        ..wall.boundary
    };
    let final_boundary = Boundary {
        source_w: 5.0,
        ..changed
    };
    let mut touch = contact(0.2);
    touch.time_s = 0.0625;
    let thermal = network(
        vec![wall, cell(2, 1, 10.0 * CP_WALL, AREA + CONTACT_AREA)],
        &[(MaterialId(1), &curve)],
        vec![
            touch,
            TimedEvent {
                time_s: 0.0625,
                event: Event::Boundary(CellId(1), changed),
            },
            TimedEvent {
                time_s: 0.375,
                event: Event::Separation(LINK),
            },
            TimedEvent {
                time_s: 0.375,
                event: Event::Boundary(CellId(1), final_boundary),
            },
        ],
    );
    let interface = patch(&binding, GasWallFilmInventory::Dry, 0.5 * AREA, 10.0);
    let controls = timeline(0.125);
    let target = 0.375;
    let total_gas_heat = 0.1875;
    let resolve =
        |thermal: &DepositionSimulation<'_>,
         gas: &fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel| {
            Ok::<_, ImplicitGasError>(coefficients([
                thermal.time_s(),
                thermal.links().len() as f64,
                thermal.cells()[&CellId(1)].boundary.source_w,
                gas.state().state().temperature,
            ]))
        };
    let next = advance_deposition_with_wall_film_timeline(
        &thermal,
        interface,
        &gas,
        &source,
        target,
        total_gas_heat,
        controls,
        resolve,
        || Ok(()),
    )
    .unwrap();
    let mut manual = thermal.clone();
    let mut manual_gas = gas.clone();
    let mut manual_film = GasWallFilmInventory::Dry;
    let mut external = 0.0;
    let mut vapor = 0.0;
    let mut last_payload = None;
    for end in [0.0625, 0.1875, 0.3125, target] {
        let start = manual.time_s();
        let heat = total_gas_heat * (end / target) - total_gas_heat * (start / target);
        let candidate = advance_deposition_with_wall_film(
            &manual,
            DepositionWallFilmBinding {
                film: manual_film,
                ..interface
            },
            &manual_gas,
            &source,
            end,
            heat,
            controls.step,
            resolve,
            || Ok(()),
        )
        .unwrap();
        let cell = &candidate.thermal().cells()[&CellId(1)];
        let unoccupied =
            candidate.thermal().exposed_area_m2(CellId(1)).unwrap() - interface.boundary.area_m2;
        let ambient = unoccupied
            * cell.boundary.htc_w_m2_k
            * (candidate.thermal().temperature_k(CellId(1)).unwrap() - cell.boundary.ambient_k);
        external += (end - start) * (cell.boundary.source_w - ambient);
        vapor += candidate.report().phase.energy_into_gas_j;
        let (advanced, film, gas, _, payload, _) = candidate.into_parts();
        let (after_events, audit) = advanced
            .advance_to(
                end,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-6,
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(audit.substeps, 0);
        manual = after_events;
        manual_film = film;
        manual_gas = gas;
        last_payload = Some(payload);
    }
    assert_eq!(next.report().accepted_intervals, 4);
    assert_eq!(next.report().events, 4);
    assert_eq!(next.report().rejected_intervals, 0);
    assert_eq!(next.thermal().time_s(), target);
    assert!(next.thermal().links().is_empty());
    assert_eq!(next.thermal().cells()[&CellId(1)].boundary, final_boundary);
    assert_eq!(next.payload().unwrap()[0], target);
    assert_eq!(next.payload().unwrap()[1], 1.0);
    assert_eq!(next.payload().unwrap()[2], 2.0);
    for id in [CellId(1), CellId(2)] {
        near(
            next.thermal().cells()[&id].specific_h_j_kg,
            manual.cells()[&id].specific_h_j_kg,
            1e-6,
        );
    }
    near(
        next.gas().state().state().temperature,
        manual_gas.state().state().temperature,
        1e-7,
    );
    near(next.film().mass_kg(), manual_film.mass_kg(), 1e-14);
    near(next.payload().unwrap()[3], last_payload.unwrap()[3], 1e-7);
    near(next.report().external_thermal_heat_j, external, 1e-8);
    near(next.report().external_gas_heat_j, total_gas_heat, 1e-14);
    near(next.report().vapor_energy_into_gas_j, vapor, 1e-8);
    near(next.report().gas_heat_apportionment_residual_j, 0.0, 1e-14);
    assert!(next.report().sum_absolute_energy_residual_j <= controls.window.energy_tolerance_j);
    balances(
        &thermal,
        &gas,
        GasWallFilmInventory::Dry,
        &binding,
        (next.thermal(), next.gas(), next.film()),
        external + total_gas_heat,
    );
}

#[test]
fn a_born_latent_neighbor_keeps_large_birth_datum_out_of_continuing_heat_and_phase() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let wall_law = plateau(300.0, 0.0, 0x88);
    let mut reference: Option<[f64; 5]> = None;
    for datum in [0.0, 2.0_f64.powi(40)] {
        let neighbor_law = plateau(310.0, datum, if datum == 0.0 { 0x89 } else { 0x8a });
        let born = cell(2, 2, datum + 50_000.0, AREA + CONTACT_AREA);
        let mut touch = contact(0.125);
        touch.time_s = 0.125;
        let thermal = network(
            vec![cell(1, 1, 50_000.0, AREA + CONTACT_AREA)],
            &[(MaterialId(1), &wall_law), (MaterialId(2), &neighbor_law)],
            vec![
                TimedEvent {
                    time_s: 0.125,
                    event: Event::Birth(born.clone()),
                },
                touch,
            ],
        );
        let next = advance_deposition_with_wall_film_timeline(
            &thermal,
            patch(&binding, GasWallFilmInventory::Dry, 0.5 * AREA, 0.0),
            &gas,
            &source,
            0.375,
            0.0,
            timeline(0.125),
            empty_heat,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(next.report().events, 2);
        assert_eq!(next.report().accepted_intervals, 3);
        assert_eq!(next.report().incoming_mass_kg, MASS);
        assert_eq!(next.report().incoming_energy_j, MASS * (datum + 50_000.0));
        assert_eq!(next.thermal().temperature_k(CellId(1)).unwrap(), 300.0);
        assert_eq!(next.thermal().temperature_k(CellId(2)).unwrap(), 310.0);
        let wall_h = next.thermal().cells()[&CellId(1)].specific_h_j_kg;
        let neighbor_delta =
            MASS * (next.thermal().cells()[&CellId(2)].specific_h_j_kg - born.specific_h_j_kg);
        // Exactly two 0.125s contacts, G=0.125W/K and deltaT=10K.
        near(neighbor_delta, -0.3125, 1e-10);
        let wall_delta = MASS * (wall_h - 50_000.0) + neighbor_delta;
        near(next.report().wall_storage_change_j, wall_delta, 1e-8);
        assert!(next.report().wall_storage_change_j.abs() < 10.0);
        near(
            MASS * (wall_h - 50_000.0),
            0.3125 - next.report().wall_heat_into_phase_system_j,
            1e-8,
        );
        let gas_delta = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
            - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
        near(
            wall_delta + gas_delta + film_energy(next.film(), &binding),
            0.0,
            1e-7,
        );
        near(next.report().combined_energy_residual_j, 0.0, 1e-7);
        near(
            next.gas().species_mass_kg(&water()).unwrap() + next.film().mass_kg(),
            gas.species_mass_kg(&water()).unwrap(),
            1e-13,
        );
        let observed = [
            wall_h,
            next.gas().state().state().temperature,
            next.film().mass_kg(),
            next.report().vapor_energy_into_gas_j,
            wall_delta,
        ];
        if let Some(reference) = reference {
            for ((actual, expected), tolerance) in observed
                .into_iter()
                .zip(reference)
                .zip([1e-5, 1e-7, 1e-14, 1e-8, 1e-7])
            {
                near(actual, expected, tolerance);
            }
        } else {
            reference = Some(observed);
        }
    }
}

#[test]
fn late_events_sources_cancellation_and_shared_work_never_publish_a_partial_timeline() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let cells = vec![
        cell(1, 1, 0.0, AREA + CONTACT_AREA),
        cell(2, 1, 10.0 * CP_WALL, AREA + CONTACT_AREA),
    ];
    let thermal = network(cells.clone(), &[(MaterialId(1), &curve)], vec![]);
    let interface = patch(&binding, GasWallFilmInventory::Dry, AREA, 0.0);
    let controls = timeline(0.0625);
    let first = advance_deposition_with_wall_film(
        &thermal,
        interface,
        &gas,
        &source,
        0.0625,
        0.0,
        controls.step,
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let queries = Cell::new(0usize);
    let polls = Cell::new(0usize);
    let accepted = advance_deposition_with_wall_film_timeline(
        &thermal,
        interface,
        &gas,
        &source,
        0.25,
        0.0,
        controls,
        |_, _| {
            queries.set(queries.get() + 1);
            Ok(coefficients(()))
        },
        || {
            polls.set(polls.get() + 1);
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(accepted.report().accepted_intervals, 4);

    let mut limited = controls;
    limited.step.implicit.thermal.maximum_evaluations = first.report().residual_evaluations;
    assert!(
        advance_deposition_with_wall_film_timeline(
            &thermal,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            limited,
            empty_heat,
            || Ok(())
        )
        .is_err()
    );
    limited = controls;
    limited.maximum_gas_iterations = first.report().gas_iterations;
    assert!(
        advance_deposition_with_wall_film_timeline(
            &thermal,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            limited,
            empty_heat,
            || Ok(())
        )
        .is_err()
    );
    limited = controls;
    limited.window.maximum_attempts = 1;
    assert!(
        advance_deposition_with_wall_film_timeline(
            &thermal,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            limited,
            empty_heat,
            || Ok(())
        )
        .is_err()
    );

    let native_limited = DepositionSimulation::new(
        cells
            .iter()
            .cloned()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect(),
        &[(MaterialId(1), &curve)],
        DepositionBudget {
            max_substeps: 1,
            ..thermal.work_budget()
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let seen = Cell::new(0.0_f64);
    assert!(
        advance_deposition_with_wall_film_timeline(
            &native_limited,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            controls,
            |thermal, _| {
                seen.set(seen.get().max(thermal.time_s()));
                Ok(coefficients(()))
            },
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(seen.get(), 0.0625);

    let mut cover = contact(0.1);
    cover.time_s = 0.125;
    if let Event::Contact(ref mut link) = cover.event {
        link.area_m2 = 2.0 * CONTACT_AREA;
    }
    let covered = network(cells, &[(MaterialId(1), &curve)], vec![cover]);
    seen.set(0.0);
    assert!(matches!(
        advance_deposition_with_wall_film_timeline(
            &covered,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            controls,
            |thermal, _| {
                seen.set(seen.get().max(thermal.time_s()));
                Ok(coefficients(()))
            },
            || Ok(())
        ),
        Err(ImplicitGasError::Invalid(_))
    ));
    assert_eq!(seen.get(), 0.125);
    assert_eq!(covered.pending_events().len(), 1);
    assert!(covered.links().is_empty());

    let last_query = queries.get();
    queries.set(0);
    assert!(matches!(
        advance_deposition_with_wall_film_timeline(
            &thermal,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            controls,
            |_, _| {
                queries.set(queries.get() + 1);
                if queries.get() == last_query {
                    Err(ImplicitGasError::Invalid("late timeline source"))
                } else {
                    Ok(coefficients(()))
                }
            },
            || Ok(())
        ),
        Err(ImplicitGasError::Invalid("late timeline source"))
    ));
    assert_eq!(queries.get(), last_query);
    let last_poll = polls.get();
    polls.set(0);
    assert!(
        advance_deposition_with_wall_film_timeline(
            &thermal,
            interface,
            &gas,
            &source,
            0.25,
            0.0,
            controls,
            empty_heat,
            || {
                polls.set(polls.get() + 1);
                if polls.get() == last_poll {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        )
        .is_err()
    );
    assert_eq!(polls.get(), last_poll);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(gas, original_gas);
    assert_eq!(interface.film, GasWallFilmInventory::Dry);
}
