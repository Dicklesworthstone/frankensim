//! Sealed normal events drive original finite-cell transport at physical times.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{motion::*, *};
use fs_couple::friction_heat::adhesive::snap::{timed::*, *};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_thermal::enthalpy::*;
use fs_tribo::{DissipationStep, HeatPartition};
use std::collections::BTreeMap;

const RELATIVE_MASS: f64 = 1.2e-5;
const CELL_MASS: f64 = 1e-6;

fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn scales(source: &SourceJkrSphereContact) -> (f64, f64, f64, f64) {
    let h = source.indentation_scale_m();
    let force = PI * 0.05 * 0.01 / 6.0;
    let speed = (force * h / RELATIVE_MASS).sqrt();
    (h, force, speed, h / speed)
}
fn motion(source: &SourceJkrSphereContact) -> JkrNormalMotion<'_> {
    let (h, force, speed, time) = scales(source);
    JkrNormalMotion::new(
        source,
        RELATIVE_MASS,
        JkrNormalMotionBudget {
            maximum_step_s: 10.0 * time,
            maximum_substep_s: 0.01 * time,
            maximum_substeps: 4096,
            maximum_root_iterations: 48,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -100.0 * h,
            maximum_accumulated_travel_m: 200.0 * h,
            maximum_normal_speed_m_s: 100.0 * speed,
            maximum_force_n: 1e4 * force,
            coordinate_absolute_tolerance_m: 1e-14 * h,
            coordinate_relative_tolerance: 1e-13,
            energy_absolute_tolerance_j: 1e-11 * force * h,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-12 * RELATIVE_MASS * speed,
            impulse_relative_tolerance: 1e-11,
        },
    )
    .unwrap()
}
fn collision(source: &SourceJkrSphereContact) -> JkrNormalMotionStep {
    let (h, _, v, time) = scales(source);
    let motion = motion(source);
    let initial = motion
        .initialize(source.detached_state(-0.25 * h).unwrap(), 10.0 * v)
        .unwrap();
    let step = motion
        .advance(initial, 0.0, 8.0 * time, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(step.event_count(), 2);
    step
}
fn release(source: &SourceJkrSphereContact) -> f64 {
    -source.evaluate_attached(0.0).unwrap().reversible_energy_j
}
fn curves(source: &SourceJkrSphereContact, widths: [f64; 2]) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let e = source.normal_state().elastic();
    let ids = [
        e.surface_a_state().resolved().card_identity(),
        e.surface_b_state().resolved().card_identity(),
    ];
    std::array::from_fn(|side| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ids[side],
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 20.0 * release(source) / CELL_MASS * widths[side],
                    temperature_k: 320.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
            ],
        )
        .unwrap()
    })
}
fn boundary(source: &SourceJkrSphereContact, side: usize) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: (side == 1).then_some(Reservoir {
            area_m2: 1e-7,
            conductance_w_k: 0.5 * release(source) / scales(source).3,
            temperature_k: 300.0,
        }),
    }
}
fn cell(source: &SourceJkrSphereContact, id: u64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id % 2),
        mass_kg: CELL_MASS,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 1e-3, 0.0, 0.0],
            volume_m3: CELL_MASS / 1000.0,
            surface_area_m2: 6e-6,
        },
        boundary: boundary(source, id as usize),
    }
}
fn topology(
    source: &SourceJkrSphereContact,
    step: &JkrNormalMotionStep,
    origin: f64,
) -> Vec<TimedEvent> {
    // Explicit independently supplied thermal closure. Only the caller's
    // authored schedule creates/removes this link; the heat adapter does not.
    step.events()
        .map(|e| TimedEvent {
            time_s: origin + step.initial_state().time_s() + e.offset_s,
            event: if e.jump.kind == JkrJumpKind::Attach {
                Event::Contact(ThermalLink {
                    id: 10,
                    a: CellId(0),
                    b: CellId(1),
                    area_m2: 1e-7,
                    conductance_w_k: release(source) / scales(source).3,
                })
            } else {
                Event::Separation(10)
            },
        })
        .collect()
}
fn thermal<'a>(
    source: &SourceJkrSphereContact,
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    origin: f64,
    mut extra: Vec<TimedEvent>,
    source_w: f64,
) -> DepositionSimulation<'a> {
    let mut a = cell(source, 0);
    a.boundary.source_w = source_w;
    let mut schedule = vec![
        TimedEvent {
            time_s: origin,
            event: Event::Birth(a),
        },
        TimedEvent {
            time_s: origin,
            event: Event::Birth(cell(source, 1)),
        },
    ];
    schedule.append(&mut extra);
    schedule.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &curves[0]), (MaterialId(1), &curves[1])],
        DepositionBudget {
            max_cells: 4,
            max_links: 2,
            max_events: 16,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(origin, controls(source).transport, || Ok(()))
    .unwrap()
    .0
}
fn partition() -> HeatPartition {
    HeatPartition::new(0.6, 0.2, 0.2).unwrap()
}
fn binding<'s>(
    source: &'s SourceJkrSphereContact,
    curves: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
    initial: JkrState,
) -> JkrSnapHeat<'s> {
    JkrSnapHeat::new(
        source,
        initial,
        thermal,
        [CellId(0), CellId(1)],
        [&curves[0], &curves[1]],
        partition(),
        Some("resolved elastic waves".into()),
    )
    .unwrap()
}
fn controls(source: &SourceJkrSphereContact) -> JkrTimedSnapHeatControls {
    JkrTimedSnapHeatControls {
        heat: JkrSnapHeatControls {
            maximum_segments: MAX_JKR_NORMAL_EVENTS + 1,
            cell_energy_tolerance_j: 1e-22,
            thermal_energy_tolerance_j: 2e-22,
            transfer_energy_tolerance_j: 4e-22,
        },
        transport: StepConfig {
            max_dt_s: 0.02 * scales(source).3,
            energy_tolerance_j: 1e-20,
        },
        maximum_transport_substeps: 5000,
        total_energy_tolerance_j: 5e-20,
    }
}
fn manual<'a>(
    source: &SourceJkrSphereContact,
    original: &DepositionSimulation<'a>,
    motion: &JkrNormalMotionStep,
    origin: f64,
) -> DepositionSimulation<'a> {
    let controls = controls(source);
    let mut next = original.clone();
    for e in motion.events() {
        next = next
            .advance_to(
                origin + (motion.initial_state().time_s() + e.offset_s),
                controls.transport,
                || Ok(()),
            )
            .unwrap()
            .0;
        let work = DissipationStep::from_work(e.jump.released_energy_j, partition()).unwrap();
        next = next
            .apply_heat_j(
                &BTreeMap::from([
                    (CellId(0), work.surface_a_heat_j()),
                    (CellId(1), work.surface_b_heat_j()),
                ]),
                controls.heat.cell_energy_tolerance_j,
                controls.heat.thermal_energy_tolerance_j,
                || Ok(()),
            )
            .unwrap()
            .0;
    }
    next.advance_to(origin + motion.state().time_s(), controls.transport, || {
        Ok(())
    })
    .unwrap()
    .0
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scaled error {:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn attachment_heat_conducts_and_cools_before_detachment_and_matches_native_manual_splitting() {
    let source = source();
    let step = collision(&source);
    let caloric = curves(&source, [1.0; 2]);
    let thermal = thermal(&source, &caloric, 0.0, topology(&source, &step, 0.0), 0.0);
    let heat = binding(
        &source,
        &caloric,
        &thermal,
        step.initial_state().contact_state(),
    );
    let expected = manual(&source, &thermal, &step, 0.0);
    let mut edited = step;
    edited.released_energy_j = 1e100;
    edited.normal_work_j = 1e100;
    edited.contact.reversible_energy_j = -1e100;
    let result = heat
        .propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            &edited,
            0.0,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(result.thermal().cells(), expected.cells());
    assert_eq!(result.thermal().links(), expected.links());
    assert_eq!(result.thermal().time_s(), expected.time_s());
    let report = result.report();
    assert_eq!(report.events.len(), 2);
    assert!(report.events[0].thermal_time_s < report.events[1].thermal_time_s);
    assert!(report.events[1].transport.external_heat_j < 0.0);
    assert!(report.transport.external_heat_j < 0.0);
    close(
        report.snap.work.dissipated_work_j(),
        step.released_energy_j,
        release(&source),
        1e-14,
    );
    close(
        result.thermal().energy_j() + report.snap.work.other_work_j()
            - report.transport.external_heat_j,
        step.released_energy_j,
        release(&source),
        1e-10,
    );
    assert!(report.absolute_stage_energy_residual_j <= controls(&source).total_energy_tolerance_j);
    let marched = thermal
        .advance_to(
            step.state().time_s(),
            controls(&source).transport,
            || Ok(()),
        )
        .unwrap()
        .0;
    let end_lumped = heat
        .propose_normal_motion(
            &marched,
            heat.stamp(),
            &step,
            controls(&source).heat,
            || Ok(()),
        )
        .unwrap();
    assert!(
        end_lumped.thermal().temperature_k(CellId(0)).unwrap()
            - result.thermal().temperature_k(CellId(0)).unwrap()
            > 0.1
    );
}

#[test]
fn exact_time_schedule_order_birth_energy_and_explicit_clock_origin_are_preserved() {
    let source = source();
    let step = collision(&source);
    let caloric = curves(&source, [1.0; 2]);
    let origin = 2.0;
    let first = origin + (step.initial_state().time_s() + step.events().next().unwrap().offset_s);
    let mut schedule = topology(&source, &step, origin);
    let mut changed = boundary(&source, 0);
    changed.reservoir = Some(Reservoir {
        area_m2: 1e-7,
        conductance_w_k: release(&source) / scales(&source).3,
        temperature_k: 300.0,
    });
    schedule.push(TimedEvent {
        time_s: first,
        event: Event::Boundary(CellId(0), changed),
    });
    let mut incoming = cell(&source, 2);
    incoming.specific_h_j_kg = release(&source) / CELL_MASS;
    schedule.push(TimedEvent {
        time_s: first,
        event: Event::Birth(incoming),
    });
    let thermal = thermal(&source, &caloric, origin, schedule, 0.0);
    let heat = binding(
        &source,
        &caloric,
        &thermal,
        step.initial_state().contact_state(),
    );
    assert!(matches!(
        heat.propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            &step,
            0.0,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Invalid(_))
    ));
    let result = heat
        .propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            &step,
            origin,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        result.thermal().cells(),
        manual(&source, &thermal, &step, origin).cells()
    );
    assert_eq!(result.report().events[0].thermal_time_s, first);
    assert_eq!(result.report().events[0].transport.events, 3); // link, boundary, birth before heat
    assert_eq!(result.report().transport.events, 4); // final separation
    close(
        result.report().transport.incoming_mass_kg,
        CELL_MASS,
        CELL_MASS,
        1e-14,
    );
    close(
        result.report().transport.incoming_energy_j,
        release(&source),
        release(&source),
        1e-14,
    );
    close(
        result.thermal().energy_j() + result.report().snap.work.other_work_j()
            - result.report().transport.external_heat_j
            - result.report().transport.incoming_energy_j,
        step.released_energy_j,
        release(&source),
        1e-10,
    );
    // The interval end still advances, but the earlier positive event offset
    // collapses on this absolute clock. It must refuse, never move the snap.
    let coarse_origin = 2.0_f64.powi(40);
    let coarse = self::thermal(&source, &caloric, coarse_origin, vec![], 0.0);
    let coarse_heat = binding(
        &source,
        &caloric,
        &coarse,
        step.initial_state().contact_state(),
    );
    assert!(coarse_origin + step.state().time_s() > coarse_origin);
    assert!(matches!(
        coarse_heat.propose_timed_normal_motion(
            &coarse,
            coarse_heat.stamp(),
            &step,
            coarse_origin,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Invalid(_))
    ));
}

#[test]
fn zero_event_intervals_still_transport_heat_and_continue_sealed_history() {
    let source = source();
    let equilibrium = source.evaluate_attached(0.0).unwrap();
    let motion = motion(&source);
    let initial = motion.initialize(equilibrium.state(), 0.0).unwrap();
    let step = motion
        .advance(
            initial,
            equilibrium.normal_force_n,
            scales(&source).3,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(step.event_count(), 0);
    let caloric = curves(&source, [1.0; 2]);
    let thermal = thermal(
        &source,
        &caloric,
        0.0,
        vec![],
        0.1 * release(&source) / scales(&source).3,
    );
    let heat = binding(&source, &caloric, &thermal, initial.contact_state());
    let expected = thermal
        .advance_to(
            step.state().time_s(),
            controls(&source).transport,
            || Ok(()),
        )
        .unwrap()
        .0;
    let result = heat
        .propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            &step,
            0.0,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(result.thermal().cells(), expected.cells());
    assert!(result.report().events.is_empty());
    assert_eq!(result.report().snap.work.dissipated_work_j(), 0.0);
    assert_eq!(result.report().snap.thermal.external_heat_j, 0.0);
    assert!(result.report().transport.external_heat_j > 0.0);
    assert_eq!(result.report().transport, result.report().final_transport);
    let (accepted, thermal, _) = result.into_parts();
    assert!(matches!(
        accepted.propose_timed_normal_motion(
            &thermal,
            accepted.stamp(),
            &step,
            0.0,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
    let next = motion
        .advance(
            step.state(),
            equilibrium.normal_force_n,
            scales(&source).3,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(
        accepted
            .propose_timed_normal_motion(
                &thermal,
                accepted.stamp(),
                &next,
                0.0,
                controls(&source),
                || Ok(())
            )
            .is_ok()
    );
}

#[test]
fn work_limits_cancellation_and_late_caloric_failures_preserve_both_owners() {
    let source = source();
    let step = collision(&source);
    let caloric = curves(&source, [1.0; 2]);
    let thermal = thermal(&source, &caloric, 0.0, topology(&source, &step, 0.0), 0.0);
    let heat = binding(
        &source,
        &caloric,
        &thermal,
        step.initial_state().contact_state(),
    );
    let original_cells = thermal.cells().clone();
    let stamp = heat.stamp();
    let mut polls = 0;
    let reference = heat
        .propose_timed_normal_motion(&thermal, stamp, &step, 0.0, controls(&source), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, polls / 2, polls - 1, polls] {
        let mut count = 0;
        assert!(matches!(
            heat.propose_timed_normal_motion(
                &thermal,
                stamp,
                &step,
                0.0,
                controls(&source),
                || {
                    count += 1;
                    if count == stop {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                }
            ),
            Err(JkrSnapHeatError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(thermal.cells(), &original_cells);
    }
    let mut limited = controls(&source);
    limited.maximum_transport_substeps = reference.report().transport.substeps - 1;
    assert!(matches!(
        heat.propose_timed_normal_motion(&thermal, stamp, &step, 0.0, limited, || Ok(())),
        Err(JkrSnapHeatError::Thermal(DepositionError::Budget(_)))
    ));
    let retry = heat
        .propose_timed_normal_motion(&thermal, stamp, &step, 0.0, controls(&source), || Ok(()))
        .unwrap();
    assert_eq!(retry.report(), reference.report());
    assert_eq!(retry.thermal().cells(), reference.thermal().cells());
    // Refuse during transport AFTER the original attachment deposit. The cold
    // boundary is scheduled exactly at attachment, so no earlier cooling occurs.
    let mut cold_schedule = topology(&source, &step, 0.0);
    let mut cold = boundary(&source, 0);
    cold.reservoir = Some(Reservoir {
        area_m2: 1e-7,
        conductance_w_k: release(&source) / scales(&source).3,
        temperature_k: 250.0,
    });
    cold_schedule.push(TimedEvent {
        time_s: step.events().next().unwrap().offset_s,
        event: Event::Boundary(CellId(0), cold),
    });
    let cold_thermal = self::thermal(&source, &caloric, 0.0, cold_schedule, 0.0);
    let cold_heat = binding(
        &source,
        &caloric,
        &cold_thermal,
        step.initial_state().contact_state(),
    );
    assert!(matches!(
        cold_heat.propose_timed_normal_motion(
            &cold_thermal,
            cold_heat.stamp(),
            &step,
            0.0,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(cold_thermal.time_s(), 0.0);
    assert_eq!(cold_heat.ledger().dissipated_work_j(), 0.0);
    assert_eq!(cold_thermal.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    // A distinct failure inside the snap itself cannot publish earlier transport.
    let narrow = curves(&source, [1.0, 1e-4]);
    let narrow_thermal = self::thermal(&source, &narrow, 0.0, topology(&source, &step, 0.0), 0.0);
    let narrow_heat = binding(
        &source,
        &narrow,
        &narrow_thermal,
        step.initial_state().contact_state(),
    );
    assert!(matches!(
        narrow_heat.propose_timed_normal_motion(
            &narrow_thermal,
            narrow_heat.stamp(),
            &step,
            0.0,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(narrow_thermal.time_s(), 0.0);
    assert_eq!(narrow_heat.ledger().dissipated_work_j(), 0.0);
}
