//! G1 original-clock fusion balances; G0 identity, exposure and rollback gates.
use super::*;
use crate::evaporation::gas::{
    equilibrium::{
        heat::advance_equilibrium_surface_with_heat,
        tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
    },
    tests::{Phase, binding, controls, gas},
};
use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositionBudget, MaterialId, ThermalLink, TimedEvent,
};

fn boundary(power_w: f64) -> Boundary {
    Boundary {
        ambient_k: MELTING,
        htc_w_m2_k: 0.0,
        source_w: power_w,
        reservoir: None,
    }
}

fn cell(chart: &EquilibriumEnthalpyPhaseCurve, id: u64, h: f64, power_w: f64) -> DepositedCell {
    let initial = bulk(chart, 0.5);
    let state = chart.values_at_specific_enthalpy(h).unwrap();
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: initial.mass_kg,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: initial.mass_kg / state.bulk_density_kg_m3,
            surface_area_m2: initial.area_m2,
        },
        boundary: boundary(power_w),
    }
}

fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}

fn schedule(
    chart: &EquilibriumEnthalpyPhaseCurve,
    events: Vec<TimedEvent>,
    maximum_steps: usize,
) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), chart)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: maximum_steps,
        },
    )
    .unwrap()
}

fn selection<'a>(
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    binding: &'a GasLiquidBinding,
    coefficient: f64,
) -> EquilibriumTimelineSurface<'a> {
    EquilibriumTimelineSurface {
        cell_id: CellId(1),
        law: &Phase,
        binding,
        inventory: EquilibriumEvaporation::new(chart).unwrap(),
        mass_transfer_coefficient_m_s: coefficient,
    }
}

fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}

fn film_heat(chart: &EquilibriumEnthalpyPhaseCurve) -> EquilibriumGasHeatSolve {
    EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [
            chart.knots()[0].specific_enthalpy_j_kg,
            chart.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

#[test]
fn birth_contact_and_heat_after_covering_share_one_original_fusion_inventory() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.5);
    let initial_gas = gas(&source, 330.0, 0.0, 0.01);
    let born = cell(&chart, 1, initial.specific_enthalpy_j_kg, 2.0);
    let geometry = born.geometry;
    let original = schedule(
        &chart,
        vec![
            birth(0.0, cell(&chart, 2, initial.specific_enthalpy_j_kg, 0.0)),
            birth(0.1, born),
            TimedEvent {
                time_s: 0.2,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: initial.area_m2,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        100,
    );
    let result = advance_deposition_with_gas_and_heat(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        0.3,
        heat(0.1),
        controls(),
        60.0,
        film_heat(&chart),
        || Ok(()),
    )
    .unwrap();
    // Only [.1,.2] has a liquid interface: birth receives no earlier flux,
    // and the contact at .2 cannot erase that interval's starting exposure.
    let post_network_heat = LiquidSurfaceState {
        specific_enthalpy_j_kg: initial.specific_enthalpy_j_kg + 2.0 * 0.1 / initial.mass_kg,
        ..initial
    };
    let one_exposed_interval = advance_equilibrium_surface_with_heat(
        EquilibriumEvaporation::new(&chart).unwrap(),
        &Phase,
        &binding,
        post_network_heat,
        &initial_gas,
        &source,
        TRANSFER,
        60.0,
        0.1,
        controls(),
        film_heat(&chart),
        || false,
    )
    .unwrap();
    let audit = result.audit();
    let phase = audit.phase_inventory.as_ref().unwrap();
    near(
        audit.transferred_mass_kg,
        one_exposed_interval.step().transferred_mass_kg,
        1e-18,
    );
    near(
        audit.sensible_heat_into_liquid_j,
        one_exposed_interval.heat().heat_into_liquid_j,
        1e-10,
    );
    assert!(audit.transferred_mass_kg > 0.0 && audit.sensible_heat_into_liquid_j > 0.0);
    near(audit.thermal.external_heat_j, 0.4, 1e-12);
    // Whole-window fusion includes native heat AFTER covering, which cannot
    // appear in a sum of the vapor-transfer substeps' phase reports.
    near(
        phase.fusion_energy_j,
        audit.transferred_mass_kg * VAPORIZATION
            - audit.thermal.external_heat_j
            - audit.sensible_heat_into_liquid_j,
        1e-8,
    );
    near(phase.initial.liquid_mass_kg, 0.005, 1e-16);
    assert_eq!(phase.final_state.equilibrium.temperature_k(), MELTING);
    assert_eq!(
        phase.final_state.equilibrium.phase_curve_identity(),
        chart.identity()
    );
    assert_eq!(result.thermal().exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(audit.last_step.as_ref().unwrap().transferred_mass_kg, 0.0);
    assert_eq!(result.thermal().cells()[&CellId(1)].geometry, geometry);
    near(
        audit.mass_residual_kg,
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    near(audit.energy_residual_j, 0.0, heat(0.1).energy_tolerance_j);
    near(
        audit.thermal_energy_change_j + audit.gas_energy_change_j
            - audit.incoming_energy_j
            - audit.thermal.external_heat_j,
        0.0,
        1e-6,
    );
    near(
        result.thermal().mass_kg() + result.gas().mass_kg() - initial_gas.mass_kg(),
        audit.thermal.incoming_mass_kg,
        1e-14,
    );
    assert!(std::ptr::eq(result.catalog(), &source));
    assert!(original.cells().is_empty());
}

#[test]
fn native_contact_conduction_freezes_pcm_without_vapor_transfer_or_external_heat() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.5);
    let solid_cp = 0.5 * CP;
    let cold_h = chart.knots()[1].specific_enthalpy_j_kg - 50.0 * solid_cp;
    let original = schedule(
        &chart,
        vec![
            birth(0.0, cell(&chart, 1, initial.specific_enthalpy_j_kg, 0.0)),
            birth(0.0, cell(&chart, 2, cold_h, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-4,
                    conductance_w_k: 0.05,
                }),
            },
        ],
        100,
    )
    .advance_to(0.0, heat(0.01), || Ok(()))
    .unwrap()
    .0;
    let initial_gas = gas(&source, MELTING, 1000.0, 0.01);
    let result = advance_deposition_with_gas(
        &original,
        selection(&chart, &binding, 0.0),
        &initial_gas,
        &source,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    // One original explicit heat substep: Q = G*(300-250)*dt. The exposed
    // cell stays isothermal while its neighbor gains the SAME extensive heat.
    let q = 0.05 * 50.0 * 0.01;
    let audit = result.audit();
    let phase = audit.phase_inventory.as_ref().unwrap();
    near(phase.fusion_energy_j, q, 1e-9);
    near(phase.freezing_mass_kg, q / FUSION, 1e-14);
    assert_eq!(phase.final_state.equilibrium.temperature_k(), MELTING);
    near(
        result.thermal().temperature_k(CellId(2)).unwrap(),
        250.0 + q / (0.01 * solid_cp),
        1e-10,
    );
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.sensible_heat_into_liquid_j, 0.0);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(result.gas(), &initial_gas);
    assert_eq!(result.thermal().mass_kg(), original.mass_kg());
    near(result.thermal().energy_j(), original.energy_j(), 1e-8);
    near(audit.energy_residual_j, 0.0, 1e-8);
}

#[test]
fn reciprocal_gas_heat_changes_phase_and_zero_coefficient_preserves_the_phase_only_lane() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.5);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, initial.specific_enthalpy_j_kg, 0.0),
        )],
        100,
    )
    .advance_to(0.0, heat(2.0), || Ok(()))
    .unwrap()
    .0;
    for gas_t in [340.0, 260.0] {
        let initial_gas = gas(&source, gas_t, 1000.0, 0.01);
        let result = advance_deposition_with_gas_and_heat(
            &original,
            selection(&chart, &binding, 0.0),
            &initial_gas,
            &source,
            2.0,
            heat(2.0),
            controls(),
            60.0,
            film_heat(&chart),
            || Ok(()),
        )
        .unwrap();
        // Independent constant-Cv gas fixture, with 1 mole percent vapor.
        let gas_capacity = 100_000.0 * 0.01 / gas_t * (0.99 * 2.5 + 0.01 * 3.0);
        let exposure = 60.0 * initial.area_m2 * 2.0;
        let q = exposure * (gas_t - MELTING) / (1.0 + exposure / gas_capacity);
        let audit = result.audit();
        let phase = audit.phase_inventory.as_ref().unwrap();
        near(audit.sensible_heat_into_liquid_j, q, 5e-9);
        near(phase.fusion_energy_j, -q, 1e-8);
        near(
            result.gas().state().state().temperature,
            gas_t - q / gas_capacity,
            1e-8,
        );
        assert_eq!(phase.final_state.equilibrium.temperature_k(), MELTING);
        assert_eq!(audit.transferred_mass_kg, 0.0);
        assert_eq!(audit.thermal.external_heat_j, 0.0);
        assert_eq!(result.thermal().mass_kg(), original.mass_kg());
        assert_eq!(result.gas().mass_kg(), initial_gas.mass_kg());
        near(
            audit.thermal_energy_change_j + audit.gas_energy_change_j,
            0.0,
            1e-7,
        );
    }
    let initial_gas = gas(&source, 310.0, 1000.0, 0.01);
    let plain = advance_deposition_with_gas(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        0.2,
        heat(0.1),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let disabled = advance_deposition_with_gas_and_heat(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        0.2,
        heat(0.1),
        controls(),
        0.0,
        EquilibriumGasHeatSolve {
            bulk_enthalpy_range_j_kg: [f64::NAN; 2],
            temperature_tolerance_k: f64::NAN,
            maximum_iterations: 0,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(disabled.thermal().cells(), plain.thermal().cells());
    assert_eq!(disabled.gas(), plain.gas());
    assert_eq!(disabled.audit(), plain.audit());
}

#[test]
fn original_source_identity_zero_window_and_late_cancellation_are_transactional() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.5);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, initial.specific_enthalpy_j_kg, 0.0),
        )],
        100,
    )
    .advance_to(0.0, heat(0.1), || Ok(()))
    .unwrap()
    .0;
    let initial_gas = gas(&source, MELTING, 0.0, 0.01);
    let saved_cells = original.cells().clone();
    let saved_gas = initial_gas.clone();
    let idle = advance_deposition_with_gas(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        0.0,
        heat(0.1),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(idle.thermal().cells(), &saved_cells);
    assert_eq!(idle.gas(), &initial_gas);
    let phase = idle.audit().phase_inventory.as_ref().unwrap();
    assert_eq!(phase.initial, phase.final_state);
    assert_eq!(phase.freezing_mass_kg, 0.0);
    assert_eq!(idle.audit().phase_steps, 0);
    // Equal numbers from a different material identity cannot stand in for
    // the timeline's registered enthalpy law, even in a zero-duration call.
    let other =
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x7d; 32]), chart.knots().to_vec())
            .unwrap();
    assert_ne!(other.identity(), chart.identity());
    assert!(matches!(
        advance_deposition_with_gas(
            &original,
            selection(&other, &binding, TRANSFER),
            &initial_gas,
            &source,
            0.0,
            heat(0.1),
            controls(),
            || Ok(()),
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let mut polls = 0;
    advance_deposition_with_gas(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        0.1,
        heat(0.1),
        controls(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert!(polls > 1);
    let mut seen = 0;
    assert!(matches!(
        advance_deposition_with_gas(
            &original,
            selection(&chart, &binding, TRANSFER),
            &initial_gas,
            &source,
            0.1,
            heat(0.1),
            controls(),
            || {
                seen += 1;
                if seen == polls {
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
    assert_eq!(original.cells(), &saved_cells);
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(initial_gas, saved_gas);
}

#[test]
fn finite_atmosphere_step_limits_retry_the_complete_interval_and_spend_native_work() {
    let source = catalog(0.0);
    let chart = curve();
    let binding = binding(0.0);
    let initial = bulk(&chart, 0.5);
    let initial_gas = gas(&source, MELTING, 0.0, 1e-5);
    let saved_gas = initial_gas.clone();
    let limited = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, initial.specific_enthalpy_j_kg, 2.0),
        )],
        1,
    );
    assert!(matches!(
        advance_deposition_with_gas(
            &limited,
            selection(&chart, &binding, TRANSFER),
            &initial_gas,
            &source,
            2.0,
            heat(2.0),
            controls(),
            || Ok(()),
        ),
        Err(GasEvaporationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert!(limited.cells().is_empty());
    assert_eq!(limited.time_s(), 0.0);
    let original = schedule(
        &chart,
        vec![birth(
            0.0,
            cell(&chart, 1, initial.specific_enthalpy_j_kg, 2.0),
        )],
        1000,
    );
    let accepted = advance_deposition_with_gas(
        &original,
        selection(&chart, &binding, TRANSFER),
        &initial_gas,
        &source,
        2.0,
        heat(2.0),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = accepted.audit();
    assert!(audit.rejected_intervals > 0);
    assert!(audit.attempted_thermal_substeps > audit.thermal.substeps);
    near(audit.thermal.external_heat_j, 4.0, 1e-12);
    let phase = audit.phase_inventory.as_ref().unwrap();
    near(
        phase.fusion_energy_j,
        audit.transferred_mass_kg * VAPORIZATION - 4.0,
        1e-7,
    );
    assert_eq!(phase.final_state.equilibrium.temperature_k(), MELTING);
    near(
        audit.mass_residual_kg,
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    near(audit.energy_residual_j, 0.0, heat(2.0).energy_tolerance_j);
    assert!(original.cells().is_empty());
    assert_eq!(initial_gas, saved_gas);
}
