//! Completion uses the same physical time for native heat, vapor and all cells.
use super::*;
use crate::evaporation::gas::{
    equilibrium::{
        heat::boundary::{FreezingHeatSolve, advance_to_freezing_with_heat},
        tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
    },
    tests::{Phase, binding, controls, gas},
    timeline::equilibrium::{
        advance_deposition_with_gas,
        tests::{birth, cell, heat, schedule, selection},
    },
};
use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{Boundary, CellId, DepositionError, Event, ThermalLink, TimedEvent};

#[test]
fn no_event_and_zero_window_preserve_the_existing_timeline_endpoint_and_accounts() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, bulk(&chart, 0.5).specific_enthalpy_j_kg, 2.0),
        )],
        100,
    );
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    for target in [0.0, 0.1] {
        let old = advance_deposition_with_gas(
            &original,
            selection(&chart, &binding, TRANSFER),
            &atmosphere,
            &source,
            target,
            heat(0.1),
            controls(),
            || Ok(()),
        )
        .unwrap();
        let result = advance_deposition_to_freezing(
            &original,
            selection(&chart, &binding, TRANSFER),
            &atmosphere,
            &source,
            target,
            heat(0.1),
            controls(),
            || Ok(()),
        )
        .unwrap();
        assert!(!result.boundary().freezing_completed);
        assert_eq!(result.boundary().elapsed_duration_s, target);
        assert_eq!(result.boundary().remaining_duration_s, 0.0);
        assert_eq!(result.exchange().thermal().cells(), old.thermal().cells());
        assert_eq!(
            result.exchange().thermal().pending_events(),
            old.thermal().pending_events()
        );
        assert_eq!(result.exchange().thermal().time_s(), old.thermal().time_s());
        assert_eq!(result.exchange().gas(), old.gas());
        assert_eq!(result.exchange().audit(), old.audit());
    }
    assert!(original.cells().is_empty());
}

#[test]
fn vapor_and_native_heat_replay_to_one_event_time_including_the_neighbor() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 1e-4);
    let cooling_w = -0.02;
    let neighbor_w = 1.0;
    let original = schedule(
        &chart,
        vec![
            birth(
                0.0,
                cell(&chart, 1, initial.specific_enthalpy_j_kg, cooling_w),
            ),
            birth(
                0.0,
                cell(
                    &chart,
                    2,
                    bulk(&chart, 0.5).specific_enthalpy_j_kg,
                    neighbor_w,
                ),
            ),
        ],
        256,
    )
    .advance_to(0.0, heat(2.0), || Ok(()))
    .unwrap()
    .0;
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    let mut phase_controls = controls();
    // At completion, one third of the pre-vapor liquid evaporates while the
    // other two thirds freeze; this event must fit the caller's donor policy.
    phase_controls.maximum_inventory_fraction = 0.5;
    let result = advance_deposition_to_freezing(
        &original,
        selection(&chart, &binding, TRANSFER),
        &atmosphere,
        &source,
        2.0,
        heat(2.0),
        phase_controls,
        || Ok(()),
    )
    .unwrap();
    let endpoint = result.exchange();
    let audit = endpoint.audit();
    let phase = audit.phase_inventory.as_ref().unwrap();
    let solid_h = chart.knots()[1].specific_enthalpy_j_kg;
    // Independent frozen-film/total-enthalpy root. The phase is isothermal
    // throughout this one split, so the initial Stefan rate remains constant.
    let r = initial.area_m2 * TRANSFER * 100_000.0 / ((0.25 * CP) * MELTING)
        * (100_000.0_f64 / 98_500.0).ln();
    let headroom = initial.mass_kg * (initial.specific_enthalpy_j_kg - solid_h);
    let expected_time = headroom / (r * (VAPORIZATION + FUSION) - cooling_w);
    let elapsed = result.boundary().elapsed_duration_s;
    near(elapsed, expected_time, 2e-5);
    near(audit.transferred_mass_kg, r * elapsed, 1e-13);
    assert!(result.boundary().freezing_completed);
    assert_eq!(endpoint.thermal().time_s(), elapsed);
    assert_eq!(result.boundary().remaining_duration_s, 2.0 - elapsed);
    near(audit.last_step.as_ref().unwrap().duration_s, elapsed, 1e-12);
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].specific_h_j_kg,
        solid_h
    );
    assert_eq!(phase.final_state.liquid_mass_kg, 0.0);
    assert!(phase.final_state.equilibrium.is_fully_solid());
    assert_eq!(
        phase.final_state.equilibrium.phase_curve_identity(),
        chart.identity()
    );
    near(
        audit.thermal.external_heat_j,
        (cooling_w + neighbor_w) * elapsed,
        1e-10,
    );
    let neighbor = &endpoint.thermal().cells()[&CellId(2)];
    let old_neighbor = &original.cells()[&CellId(2)];
    near(
        neighbor.mass_kg * (neighbor.specific_h_j_kg - old_neighbor.specific_h_j_kg),
        neighbor_w * elapsed,
        1e-9,
    );
    near(
        phase.fusion_energy_j,
        audit.transferred_mass_kg * VAPORIZATION - cooling_w * elapsed,
        2e-6,
    );
    near(
        endpoint.thermal().mass_kg() + endpoint.gas().mass_kg(),
        original.mass_kg() + atmosphere.mass_kg(),
        phase_controls.gas.mass_tolerance_kg,
    );
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j - audit.thermal.external_heat_j,
        0.0,
        2e-6,
    );
    assert!(audit.attempted_thermal_substeps > audit.thermal.substeps);
    assert!(audit.rejected_intervals > 0);
}

#[test]
fn fully_covered_cell_still_completes_native_freezing_and_preserves_the_gas() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.001);
    let original = schedule(
        &chart,
        vec![
            birth(0.0, cell(&chart, 1, initial.specific_enthalpy_j_kg, -4.0)),
            birth(
                0.0,
                cell(&chart, 2, bulk(&chart, 0.5).specific_enthalpy_j_kg, 2.0),
            ),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: initial.area_m2,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        256,
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    assert_eq!(original.exposed_area_m2(CellId(1)), Ok(0.0));
    let result = advance_deposition_to_freezing(
        &original,
        selection(&chart, &binding, TRANSFER),
        &atmosphere,
        &source,
        1.0,
        heat(1.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let endpoint = result.exchange();
    let audit = endpoint.audit();
    let phase = audit.phase_inventory.as_ref().unwrap();
    // One original joule of latent headroom, removed at four watts. Coverage
    // suppresses vapor transfer, but cannot suppress bulk thermal completion.
    near(result.boundary().elapsed_duration_s, 0.25, 3e-7);
    near(result.boundary().remaining_duration_s, 0.75, 3e-7);
    assert!(result.boundary().freezing_completed);
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].specific_h_j_kg,
        chart.knots()[1].specific_enthalpy_j_kg
    );
    assert_eq!(phase.final_state.liquid_mass_kg, 0.0);
    near(phase.fusion_energy_j, 1.0, 1e-9);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(endpoint.gas(), &atmosphere);
    near(
        audit.thermal.external_heat_j,
        -2.0 * endpoint.thermal().time_s(),
        1e-10,
    );
    near(
        endpoint.thermal().energy_j() - original.energy_j(),
        -0.5,
        2e-6,
    );
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].geometry,
        original.cells()[&CellId(1)].geometry
    );
    assert_eq!(endpoint.thermal().mass_kg(), original.mass_kg());
}

#[test]
fn exact_saturation_allows_native_freezing_without_bypassing_original_source_admission() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.001);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, initial.specific_enthalpy_j_kg, -4.0),
        )],
        256,
    );
    let atmosphere = gas(&source, MELTING, 1500.0, 0.01);
    let result = advance_deposition_to_freezing(
        &original,
        selection(&chart, &binding, TRANSFER),
        &atmosphere,
        &source,
        1.0,
        heat(1.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let endpoint = result.exchange();
    assert!(result.boundary().freezing_completed);
    near(endpoint.thermal().time_s(), 0.25, 3e-7);
    assert_eq!(
        endpoint.thermal().exposed_area_m2(CellId(1)),
        Ok(initial.area_m2)
    );
    assert_eq!(endpoint.audit().transferred_mass_kg, 0.0);
    assert_eq!(endpoint.gas(), &atmosphere);
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].specific_h_j_kg,
        chart.knots()[1].specific_enthalpy_j_kg
    );
    assert_eq!(
        endpoint
            .audit()
            .phase_inventory
            .as_ref()
            .unwrap()
            .final_state
            .liquid_mass_kg,
        0.0
    );
    let wrong_binding = crate::evaporation::gas::GasLiquidBinding {
        liquid_enthalpy_offset_j_kg: 1000.0,
        ..binding.clone()
    };
    assert!(matches!(
        advance_deposition_to_freezing(
            &original,
            selection(&chart, &wrong_binding, TRANSFER),
            &atmosphere,
            &source,
            1.0,
            heat(1.0),
            controls(),
            || Ok(()),
        ),
        Err(GasEvaporationTimelineError::Exchange(
            crate::evaporation::gas::GasEvaporationError::CaloricMismatch { .. }
        ))
    ));
    assert!(original.cells().is_empty());
}

#[test]
fn delayed_birth_and_same_time_boundary_events_keep_the_original_order_at_completion() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.001);
    let original = schedule(
        &chart,
        vec![
            birth(
                0.0,
                cell(&chart, 2, bulk(&chart, 0.5).specific_enthalpy_j_kg, 0.0),
            ),
            birth(0.125, cell(&chart, 1, initial.specific_enthalpy_j_kg, -4.0)),
            TimedEvent {
                time_s: 0.25,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: initial.area_m2,
                    conductance_w_k: 0.0,
                }),
            },
            TimedEvent {
                time_s: 0.375,
                event: Event::Boundary(
                    CellId(1),
                    Boundary {
                        ambient_k: MELTING,
                        htc_w_m2_k: 0.0,
                        source_w: 4.0,
                        reservoir: None,
                    },
                ),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Separation(7),
            },
        ],
        256,
    );
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    let result = advance_deposition_to_freezing(
        &original,
        selection(&chart, &binding, 0.0),
        &atmosphere,
        &source,
        1.0,
        heat(1.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let endpoint = result.exchange();
    assert!(result.boundary().freezing_completed);
    assert_eq!(endpoint.thermal().time_s(), 0.375);
    assert_eq!(result.boundary().elapsed_duration_s, 0.375);
    assert_eq!(result.boundary().remaining_duration_s, 0.625);
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].boundary.source_w,
        4.0
    );
    assert_eq!(endpoint.thermal().exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(endpoint.thermal().pending_events().len(), 1);
    assert_eq!(endpoint.thermal().pending_events()[0].time_s, 0.5);
    assert_eq!(endpoint.audit().thermal.events, 4);
    near(endpoint.audit().thermal.external_heat_j, -1.0, 1e-12);
    near(
        endpoint
            .audit()
            .phase_inventory
            .as_ref()
            .unwrap()
            .fusion_energy_j,
        1.0,
        1e-9,
    );
    assert_eq!(endpoint.gas(), &atmosphere);
    assert!(original.cells().is_empty());
    assert_eq!(original.pending_events().len(), 5);
}

#[test]
fn discarded_boundary_probes_spend_the_original_native_work_budget() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, bulk(&chart, 0.001).specific_enthalpy_j_kg, -4.0),
        )],
        1,
    );
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    let saved_gas = atmosphere.clone();
    assert!(matches!(
        advance_deposition_to_freezing(
            &original,
            selection(&chart, &binding, 0.0),
            &atmosphere,
            &source,
            1.0,
            heat(1.0),
            controls(),
            || Ok(()),
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.pending_events().len(), 1);
    assert_eq!(atmosphere, saved_gas);
}

#[test]
fn original_source_and_initial_liquid_admission_survive_noops_and_late_cancellation() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, bulk(&chart, 0.001).specific_enthalpy_j_kg, -4.0),
        )],
        256,
    );
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    let other =
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x7e; 32]), chart.knots().to_vec())
            .unwrap();
    // A future/zero-window birth still binds the exact registered chart.
    assert!(matches!(
        advance_deposition_to_freezing(
            &original,
            selection(&other, &binding, 0.0),
            &atmosphere,
            &source,
            0.0,
            heat(1.0),
            controls(),
            || Ok(()),
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let already_solid = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, chart.knots()[1].specific_enthalpy_j_kg, 0.0),
        )],
        256,
    );
    assert!(
        advance_deposition_to_freezing(
            &already_solid,
            selection(&chart, &binding, 0.0),
            &atmosphere,
            &source,
            0.0,
            heat(1.0),
            controls(),
            || Ok(()),
        )
        .is_err()
    );
    let mut calls = 0;
    advance_deposition_to_freezing(
        &original,
        selection(&chart, &binding, 0.0),
        &atmosphere,
        &source,
        1.0,
        heat(1.0),
        controls(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let final_gate = calls;
    let mut calls = 0;
    assert!(matches!(
        advance_deposition_to_freezing(
            &original,
            selection(&chart, &binding, 0.0),
            &atmosphere,
            &source,
            1.0,
            heat(1.0),
            controls(),
            || {
                calls += 1;
                if calls == final_gate {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.pending_events().len(), 1);
}

#[test]
fn native_heat_gas_heat_and_vapor_complete_freezing_on_the_neighbor_clock() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 1e-4);
    let cooling_w = -0.02;
    let neighbor_w = 1.0;
    let htc = 1.0;
    let target = 2.0;
    let original = schedule(
        &chart,
        vec![
            birth(
                0.0,
                cell(&chart, 1, initial.specific_enthalpy_j_kg, cooling_w),
            ),
            birth(
                0.0,
                cell(
                    &chart,
                    2,
                    bulk(&chart, 0.5).specific_enthalpy_j_kg,
                    neighbor_w,
                ),
            ),
        ],
        256,
    )
    .advance_to(0.0, heat(target), || Ok(()))
    .unwrap()
    .0;
    let atmosphere = gas(&source, 280.0, 0.0, 0.01);
    let mut phase_controls = controls();
    phase_controls.maximum_inventory_fraction = 0.5;
    let gas_heat = EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [
            chart.knots()[0].specific_enthalpy_j_kg,
            chart.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    };
    let surface = selection(&chart, &binding, TRANSFER);
    let result = advance_deposition_to_freezing_with_heat(
        &original,
        surface,
        &atmosphere,
        &source,
        target,
        heat(target),
        phase_controls,
        htc,
        gas_heat,
        || Ok(()),
    )
    .unwrap();
    let endpoint = result.exchange();
    let audit = endpoint.audit();
    let phase = audit.phase_inventory.as_ref().unwrap();
    let elapsed = result.boundary().elapsed_duration_s;
    let solid_h = chart.knots()[1].specific_enthalpy_j_kg;
    assert!(result.boundary().freezing_completed && elapsed > 0.0 && elapsed < target);
    assert_eq!(endpoint.thermal().time_s(), elapsed);
    assert_eq!(result.boundary().remaining_duration_s, target - elapsed);
    assert_eq!(audit.last_step.as_ref().unwrap().duration_s, elapsed);
    assert_eq!(
        endpoint.thermal().cells()[&CellId(1)].specific_h_j_kg,
        solid_h
    );
    assert_eq!(phase.final_state.liquid_mass_kg, 0.0);
    assert!(audit.transferred_mass_kg > 0.0 && audit.sensible_heat_into_liquid_j < 0.0);

    // Independent finite-gas Cv and isothermal backward-Euler heat balance.
    // The original dry-gas Cv/R is 2.5; native heat keeps the body at melting.
    let capacity = 100_000.0 * 0.01 / 280.0 * 2.5;
    let exposure = htc * initial.area_m2 * elapsed;
    let expected_heat = exposure * (280.0 - MELTING) / (1.0 + exposure / capacity);
    let heated_gas_t = 280.0 - expected_heat / capacity;
    let heated_pressure = 100_000.0 * heated_gas_t / 280.0;
    let expected_mass = initial.area_m2 * TRANSFER * heated_pressure / ((0.25 * CP) * heated_gas_t)
        * (heated_pressure / (heated_pressure - 1500.0)).ln()
        * elapsed;
    let headroom = initial.mass_kg * (initial.specific_enthalpy_j_kg - solid_h);
    near(audit.sensible_heat_into_liquid_j, expected_heat, 1e-8);
    near(audit.transferred_mass_kg, expected_mass, 1e-13);
    near(
        headroom + cooling_w * elapsed + expected_heat - expected_mass * (FUSION + VAPORIZATION),
        0.0,
        2e-6,
    );
    near(
        phase.fusion_energy_j,
        audit.transferred_mass_kg * VAPORIZATION
            - cooling_w * elapsed
            - audit.sensible_heat_into_liquid_j,
        2e-6,
    );

    // Every native cell shares the terminal duration, including the neighbor
    // which has no vapor interface or reciprocal gas-heat source of its own.
    let neighbor = &endpoint.thermal().cells()[&CellId(2)];
    let old_neighbor = &original.cells()[&CellId(2)];
    near(
        neighbor.mass_kg * (neighbor.specific_h_j_kg - old_neighbor.specific_h_j_kg),
        neighbor_w * elapsed,
        1e-9,
    );
    near(
        audit.thermal.external_heat_j,
        (cooling_w + neighbor_w) * elapsed,
        1e-10,
    );
    let before = original.energy_j()
        + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
    let after = endpoint.thermal().energy_j()
        + endpoint.gas().mass_kg() * endpoint.gas().state().specific_internal_energy_j_kg();
    near(after - before, (cooling_w + neighbor_w) * elapsed, 2e-6);
    near(
        endpoint.thermal().mass_kg() + endpoint.gas().mass_kg(),
        original.mass_kg() + atmosphere.mass_kg(),
        phase_controls.gas.mass_tolerance_kg,
    );

    // Replay actual native heat to the accepted time, then ask the original
    // scalar gas-heat/phase owner to consume exactly that same interval.
    let native = original
        .advance_to(elapsed, heat(target), || Ok(()))
        .unwrap()
        .0;
    let native_cell = &native.cells()[&CellId(1)];
    let mut post_native = initial;
    post_native.mass_kg = native_cell.mass_kg;
    post_native.specific_enthalpy_j_kg = native_cell.specific_h_j_kg;
    post_native.area_m2 = native.exposed_area_m2(CellId(1)).unwrap();
    let replay = advance_to_freezing_with_heat(
        surface.inventory,
        &Phase,
        &binding,
        post_native,
        &atmosphere,
        &source,
        TRANSFER,
        htc,
        elapsed,
        phase_controls,
        gas_heat,
        FreezingHeatSolve {
            maximum_iterations: 128,
        },
        || false,
    )
    .unwrap();
    assert!(replay.boundary().freezing_completed);
    assert_eq!(replay.boundary().elapsed_duration_s, elapsed);
    assert_eq!(replay.boundary().remaining_duration_s, 0.0);
    assert_eq!(
        replay.bulk().mass_kg,
        endpoint.thermal().cells()[&CellId(1)].mass_kg
    );
    assert_eq!(replay.bulk().specific_enthalpy_j_kg, solid_h);
    assert_eq!(replay.gas(), endpoint.gas());
    assert_eq!(replay.step(), audit.last_step.as_ref().unwrap());
    assert_eq!(
        replay.heat().heat_into_liquid_j,
        audit.sensible_heat_into_liquid_j
    );
    assert_eq!(
        native.cells()[&CellId(2)],
        endpoint.thermal().cells()[&CellId(2)]
    );
    assert_eq!(
        original.cells()[&CellId(1)].specific_h_j_kg,
        initial.specific_enthalpy_j_kg
    );
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn zero_gas_heat_exactly_delegates_completion_despite_unused_invalid_heat_controls() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, bulk(&chart, 1e-4).specific_enthalpy_j_kg, -0.02),
        )],
        256,
    );
    let atmosphere = gas(&source, MELTING, 0.0, 0.01);
    let mut phase_controls = controls();
    phase_controls.maximum_inventory_fraction = 0.5;
    let unused = EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [f64::NAN; 2],
        temperature_tolerance_k: f64::NAN,
        maximum_iterations: 0,
    };
    for target in [0.0, 2.0] {
        let old = advance_deposition_to_freezing(
            &original,
            selection(&chart, &binding, TRANSFER),
            &atmosphere,
            &source,
            target,
            heat(2.0),
            phase_controls,
            || Ok(()),
        )
        .unwrap();
        let new = advance_deposition_to_freezing_with_heat(
            &original,
            selection(&chart, &binding, TRANSFER),
            &atmosphere,
            &source,
            target,
            heat(2.0),
            phase_controls,
            0.0,
            unused,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(new.boundary(), old.boundary());
        assert_eq!(
            new.exchange().thermal().time_s(),
            old.exchange().thermal().time_s()
        );
        assert_eq!(
            new.exchange().thermal().cells(),
            old.exchange().thermal().cells()
        );
        assert_eq!(
            new.exchange().thermal().pending_events(),
            old.exchange().thermal().pending_events()
        );
        assert_eq!(new.exchange().gas(), old.exchange().gas());
        assert_eq!(new.exchange().audit(), old.exchange().audit());
        if target > 0.0 {
            assert!(new.boundary().freezing_completed);
        }
    }
    assert!(original.cells().is_empty());
    assert_eq!(original.time_s(), 0.0);
}
