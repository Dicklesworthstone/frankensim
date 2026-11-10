//! Plastic heat enters original finite cells while native transport evolves.

use super::super::super::MASS;
use super::super::{CRITICAL_SHIFT, binding, close, controls, events, loaded, source};
use fs_contact::normal_patch::adhesive::SourceJkrSphereContact;
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{
        JkrRollingSnapHeatError,
        dissipation::{JkrRollingDissipationHeat, timed::JkrTimedRollingHeatControls},
    },
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLink, TimedEvent,
};
use fs_tribo::HeatPartition;

const STEP: f64 = 1.0 / 16.0;

fn timed_controls() -> JkrTimedRollingHeatControls {
    JkrTimedRollingHeatControls {
        heat: controls(),
        transport: StepConfig {
            max_dt_s: STEP,
            energy_tolerance_j: 2e-22,
        },
        maximum_transport_substeps: 64,
        total_energy_tolerance_j: 8e-22,
    }
}

// Manufactured constant heat capacity with the original two material-card
// identities. The small capacity resolves source-owned adhesion heat directly.
fn linear_curves(
    normal: &SourceJkrSphereContact,
    heat_capacity_j_k: f64,
) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let elastic = normal.normal_state().elastic();
    [
        elastic.surface_a_state().resolved().card_identity(),
        elastic.surface_b_state().resolved().card_identity(),
    ]
    .map(|card| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card,
            fs_material::phase::SolidLiquidPhase::Solid,
            [0.0, 1000.0]
                .map(|rise| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: rise * heat_capacity_j_k / MASS,
                    temperature_k: 290.0 + rise,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}

fn boundary(conductance_w_k: f64) -> Boundary {
    Boundary {
        ambient_k: 290.0,
        htc_w_m2_k: 0.0,
        reservoir: Some(Reservoir {
            area_m2: 1e-6,
            conductance_w_k,
            temperature_k: 290.0,
        }),
        source_w: 0.0,
    }
}

fn timeline(
    caloric: &[EquilibriumEnthalpyPhaseCurve; 2],
    reservoir_w_k: f64,
    link_w_k: f64,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'_> {
    let mut schedule: Vec<_> = (0..2)
        .map(|side| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side),
                segment_id: side,
                body_id: side,
                material_id: MaterialId(side),
                mass_kg: MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [side as f64, 0.0, 0.0],
                    end_m: [side as f64 + 1e-3, 0.0, 0.0],
                    volume_m3: 1e-9,
                    surface_area_m2: 6e-6,
                },
                boundary: boundary(reservoir_w_k),
            }),
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 0,
            a: CellId(0),
            b: CellId(1),
            area_m2: 1e-6,
            conductance_w_k: link_w_k,
        }),
    });
    schedule.extend(later);
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, timed_controls().transport, || Ok(()))
    .unwrap()
    .0
}

#[test]
fn plastic_power_heats_during_native_conduction_and_reservoir_cooling() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_rolling_loss_j();
    let angle = 2.0 * CRITICAL_SHIFT / normal.effective_radius_m();
    let original = events.advance_rolling(&initial, angle).unwrap();
    let loss = original.report().irreversible_rolling_loss_j;
    assert!(loss > 0.0);
    let caloric = linear_curves(&normal, loss);
    // C = loss J/K, reservoir G = C/s, intercell G = C/s. The sum
    // and difference modes therefore decay at 1/s and 3/s respectively.
    let mut thermal = timeline(&caloric, loss, loss, vec![]);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, initial, &caloric, &thermal),
        HeatPartition::new(0.5, 0.25, 0.25).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let report = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            0.0,
            angle,
            1.0,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(thermal.time_s(), 1.0);
    assert_eq!(report.transport.substeps, 16);
    assert_eq!(heat.state(), original.state());
    assert_eq!(report.normal.work.dissipated_work_j(), 0.0);
    close(heat.state().accumulated_rolling_loss_j(), prior_loss + loss);
    close(heat.rolling_loss_ledger().dissipated_work_j(), loss);
    close(report.injected_surface_heat_j[0], 0.5 * loss);
    close(report.injected_surface_heat_j[1], 0.25 * loss);
    let sum = 0.75 * loss * (1.0 - (1.0 - STEP).powi(16));
    let difference = 0.25 * loss * (1.0 - (1.0 - 3.0 * STEP).powi(16)) / 3.0;
    for (side, expected) in [(0, (sum + difference) / 2.0), (1, (sum - difference) / 2.0)] {
        let actual = MASS * thermal.cells()[&CellId(side)].specific_h_j_kg;
        close(actual, expected);
    }
    close(report.external_transport_heat_j, sum - 0.75 * loss);
    assert!(report.external_transport_heat_j < 0.0);
    assert!(sum < 0.7 * 0.75 * loss); // An endpoint impulse cannot produce this result.
    assert!(report.absolute_stage_energy_residual_j <= timed_controls().total_energy_tolerance_j);
    // A second accepted interval uses the existing source history and clock.
    let next = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            0.0,
            0.0,
            0.5,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.initial_time_s, 1.0);
    assert_eq!(thermal.time_s(), 1.5);
    assert_eq!(next.rolling_loss.dissipated_work_j(), 0.0);
    close(heat.rolling_loss_ledger().dissipated_work_j(), loss);
}

#[test]
fn original_snap_is_at_start_and_scheduled_cooling_is_at_its_physical_time() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = events.initial_detached_state(-1e-6).unwrap();
    let release = events
        .advance_normal(&initial, 0.0)
        .unwrap()
        .report()
        .normal_released_energy_j;
    let caloric = linear_curves(&normal, release);
    let changes = (0..2)
        .map(|side| TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(CellId(side), boundary(2.0 * release)),
        })
        .collect();
    let mut thermal = timeline(&caloric, 0.0, 0.0, changes);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, initial, &caloric, &thermal),
        HeatPartition::new(0.5, 0.25, 0.25).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let report = heat
        .apply_interval(
            &mut thermal,
            heat.stamp(),
            0.0,
            0.0,
            1.0,
            timed_controls(),
            || Ok(()),
        )
        .unwrap();
    close(report.normal.normal_release.dissipated_work_j(), release);
    close(report.normal.thermal.stored_energy_change_j, 0.75 * release);
    assert_eq!(report.transport.events, 2);
    assert_eq!(report.rolling_loss.dissipated_work_j(), 0.0);
    assert_eq!(report.injected_surface_heat_j, [0.0; 2]);
    for (side, fraction) in [(0, 0.25), (1, 0.5)] {
        close(
            MASS * thermal.cells()[&CellId(side)].specific_h_j_kg,
            fraction * release * (1.0 - 2.0 * STEP).powi(8),
        );
    }
    close(
        report.thermal.external_heat_j,
        report.thermal.stored_energy_change_j,
    );
    assert_eq!(
        report.thermal.events,
        report.initial_events.events + report.transport.events
    );
    assert_eq!(thermal.time_s(), 1.0);
}

#[test]
fn budget_final_cancellation_and_stale_clock_leave_every_owner_unchanged() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = events.initial_detached_state(-1e-6).unwrap();
    let release = events
        .advance_normal(&initial, 0.0)
        .unwrap()
        .report()
        .normal_released_energy_j;
    let caloric = linear_curves(&normal, release);
    let mut thermal = timeline(&caloric, release, 0.0, vec![]);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, initial, &caloric, &thermal),
        HeatPartition::new(0.5, 0.25, 0.25).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let cells = thermal.cells().clone();
    let angle = 3.0 * CRITICAL_SHIFT / normal.effective_radius_m();
    let mut bounded = timed_controls();
    bounded.maximum_transport_substeps = 1;
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp.clone(),
            0.0,
            angle,
            1.0,
            bounded,
            || Ok(()),
        ),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heat.normal_release_ledger().dissipated_work_j(), 0.0);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    let mut proposal_calls = 0usize;
    heat.propose_interval(
        &thermal,
        stamp.clone(),
        0.0,
        angle,
        1.0,
        timed_controls(),
        || {
            proposal_calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut calls = 0usize;
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp.clone(),
            0.0,
            angle,
            1.0,
            timed_controls(),
            || {
                calls += 1;
                if calls == proposal_calls + 1 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    // Moving only the thermal owner cannot silently resynchronize contact history.
    thermal = thermal
        .advance_to(STEP, timed_controls().transport, || Ok(()))
        .unwrap()
        .0;
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp.clone(),
            0.0,
            angle,
            1.0,
            timed_controls(),
            || Ok(()),
        ),
        Err(JkrRollingSnapHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), STEP);
}

#[test]
fn finite_caloric_refusal_discards_partial_transport_and_plastic_history() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let angle = 2.0 * CRITICAL_SHIFT / normal.effective_radius_m();
    let loss = events
        .advance_rolling(&initial, angle)
        .unwrap()
        .report()
        .irreversible_rolling_loss_j;
    let caloric = linear_curves(&normal, loss * 1e-6);
    let mut thermal = timeline(&caloric, 0.0, 0.0, vec![]);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, initial, &caloric, &thermal),
        HeatPartition::new(0.5, 0.25, 0.25).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let cells = thermal.cells().clone();
    assert!(matches!(
        heat.apply_interval(
            &mut thermal,
            stamp.clone(),
            0.0,
            angle,
            1.0,
            timed_controls(),
            || Ok(()),
        ),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &cells);
}
