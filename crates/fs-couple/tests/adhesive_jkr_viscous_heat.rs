//! Original viscous support work becomes native heat over its physical intervals.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{motion::*, *};
use fs_couple::friction_heat::adhesive::snap::{timed::*, viscous::*, *};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_thermal::enthalpy::*;
use fs_tribo::{DissipationStep, HeatPartition};
use std::collections::BTreeMap;

const MASS: f64 = 1.2e-5;
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
    let velocity = (force * h / MASS).sqrt();
    (h, force, velocity, h / velocity)
}
fn release(source: &SourceJkrSphereContact) -> f64 {
    -source.evaluate_attached(0.0).unwrap().reversible_energy_j
}
fn motion(source: &SourceJkrSphereContact, damping: f64) -> JkrNormalMotion<'_> {
    let (h, force, velocity, time) = scales(source);
    JkrNormalMotion::new_with_damped_support(
        source,
        MASS,
        JkrNormalMotionBudget {
            maximum_step_s: 10.0 * time,
            maximum_substep_s: 0.01 * time,
            maximum_substeps: 4096,
            maximum_root_iterations: 48,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -100.0 * h,
            maximum_accumulated_travel_m: 200.0 * h,
            maximum_normal_speed_m_s: 100.0 * velocity,
            maximum_force_n: 1e4 * force,
            coordinate_absolute_tolerance_m: 1e-14 * h,
            coordinate_relative_tolerance: 1e-13,
            energy_absolute_tolerance_j: 1e-11 * force * h,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-12 * MASS * velocity,
            impulse_relative_tolerance: 1e-11,
        },
        0.0,
        0.0,
        damping,
    )
    .unwrap()
}
fn constant_motion(
    source: &SourceJkrSphereContact,
) -> (JkrNormalMotion<'_>, JkrNormalMotionTrace, f64) {
    let (h, _, velocity, time) = scales(source);
    let velocity = 2.0 * velocity;
    let damping = release(source) / (velocity * velocity * time);
    let law = motion(source, damping);
    let state = law
        .initialize(source.detached_state(-20.0 * h).unwrap(), velocity)
        .unwrap();
    let load = damping * velocity;
    let trace = law
        .advance_with_dissipation_trace(state, load, time, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(trace.step().event_count(), 0);
    (law, trace, load)
}
fn collision(source: &SourceJkrSphereContact) -> JkrNormalMotionTrace {
    let (h, _, velocity, time) = scales(source);
    let law = motion(source, 0.01 * MASS / time);
    let state = law
        .initialize(source.detached_state(-0.25 * h).unwrap(), 10.0 * velocity)
        .unwrap();
    let trace = law
        .advance_with_dissipation_trace(state, 0.0, 8.0 * time, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(trace.step().event_count(), 2);
    trace
}
fn curves(source: &SourceJkrSphereContact, width: f64) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let e = source.normal_state().elastic();
    [
        e.surface_a_state().resolved().card_identity(),
        e.surface_b_state().resolved().card_identity(),
    ]
    .map(|card| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card,
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 20.0 * release(source) * width / CELL_MASS,
                    temperature_k: 320.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                },
            ],
        )
        .unwrap()
    })
}
fn boundary(source: &SourceJkrSphereContact, watts: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w: watts,
        reservoir: Some(Reservoir {
            area_m2: 1e-7,
            conductance_w_k: release(source) / scales(source).3,
            temperature_k: 300.0,
        }),
    }
}
fn thermal<'m>(
    source: &SourceJkrSphereContact,
    caloric: &'m [EquilibriumEnthalpyPhaseCurve; 2],
    origin: f64,
    mut schedule: Vec<TimedEvent>,
    watts: f64,
) -> DepositionSimulation<'m> {
    for id in 0..2 {
        let mut b = boundary(source, if id == 0 { watts } else { 0.0 });
        if id == 1 {
            b.reservoir = None;
        }
        schedule.push(TimedEvent {
            time_s: origin,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(id),
                mass_kg: CELL_MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 1e-3, 0.0, 0.0],
                    volume_m3: CELL_MASS / 1000.0,
                    surface_area_m2: 6e-6,
                },
                boundary: b,
            }),
        });
    }
    schedule.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(origin, controls(source).transport, || Ok(()))
    .unwrap()
    .0
}
fn snap_partition() -> HeatPartition {
    HeatPartition::new(0.6, 0.2, 0.2).unwrap()
}
fn route() -> JkrViscousHeatRoute {
    JkrViscousHeatRoute::new(
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("support damping reservoir".into()),
    )
    .unwrap()
}
fn binding<'s>(
    source: &'s SourceJkrSphereContact,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
    trace: &JkrNormalMotionTrace,
) -> JkrSnapHeat<'s> {
    JkrSnapHeat::new(
        source,
        trace.step().initial_state().contact_state(),
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        snap_partition(),
        Some("surface snap waves".into()),
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
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scaled error {:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn constant_support_power_preserves_authored_sources_and_converges_to_the_caloric_reservoir_solution()
 {
    let source = source();
    let (_, trace, _) = constant_motion(&source);
    let caloric = curves(&source, 1.0);
    let time = scales(&source).3;
    let power = release(&source) / time;
    let changed = boundary(&source, 0.2 * power);
    let thermal = thermal(
        &source,
        &caloric,
        0.0,
        vec![TimedEvent {
            time_s: 0.4 * time,
            event: Event::Boundary(CellId(0), changed),
        }],
        0.1 * power,
    );
    let heat = binding(&source, &caloric, &thermal, &trace);
    let result = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            heat.stamp(),
            &trace,
            0.0,
            &route(),
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    let report = result.report();
    close(
        trace.step().damping_dissipation_j,
        release(&source),
        release(&source),
        1e-11,
    );
    close(
        report.support.injected_surface_heat_j[0],
        0.25 * release(&source),
        release(&source),
        1e-11,
    );
    close(
        report.support.injected_surface_heat_j[1],
        0.5 * release(&source),
        release(&source),
        1e-11,
    );
    assert_eq!(report.timed.snap.work.dissipated_work_j(), 0.0);
    assert_eq!(result.thermal().cells()[&CellId(0)].boundary, changed);
    assert_eq!(result.thermal().cells()[&CellId(1)].boundary.source_w, 0.0);
    close(
        result.thermal().temperature_k(CellId(1)).unwrap(),
        300.5,
        1.0,
        1e-11,
    );
    let first = 0.35 * (1.0 - (-0.4_f64).exp());
    let exact = 300.0 + 0.45 + (first - 0.45) * (-0.6_f64).exp();
    let error = (result.thermal().temperature_k(CellId(0)).unwrap() - exact).abs();
    assert!(error < 0.002);
    let mut refined = controls(&source);
    refined.transport.max_dt_s = 0.005 * time;
    let fine = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            heat.stamp(),
            &trace,
            0.0,
            &route(),
            refined,
            || Ok(()),
        )
        .unwrap();
    let fine_error = (fine.thermal().temperature_k(CellId(0)).unwrap() - exact).abs();
    assert!(fine_error < 0.6 * error, "{fine_error:e} versus {error:e}");
    close(
        result.thermal().energy_j() + report.support.work.other_work_j()
            - report.external_transport_heat_j,
        trace.step().damping_dissipation_j,
        release(&source),
        1e-10,
    );
}

#[derive(Clone)]
struct Powers([f64; 2]);
impl DepositionObserver for Powers {
    type Error = DepositionError;
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        _: f64,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([
            (CellId(0), self.0[0]),
            (CellId(1), self.0[1]),
        ]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}
fn manual<'m>(
    source: &SourceJkrSphereContact,
    thermal: &DepositionSimulation<'m>,
    trace: &JkrNormalMotionTrace,
) -> DepositionSimulation<'m> {
    let mut cuts = vec![0.0, trace.step().state().time_s()];
    for segment in trace.segments() {
        cuts.extend([segment.start_offset_s(), segment.end_offset_s()]);
    }
    cuts.extend(trace.step().events().map(|event| event.offset_s));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut next = thermal.clone();
    let mut index = 0;
    let config = controls(source);
    for end in cuts {
        while trace
            .segments()
            .get(index)
            .is_some_and(|segment| segment.end_offset_s() <= next.time_s())
        {
            index += 1;
        }
        let power = trace
            .segments()
            .get(index)
            .filter(|segment| segment.start_offset_s() <= next.time_s())
            .map_or(0.0, |segment| segment.power_w());
        next = next
            .advance_observed(
                &Powers([0.25 * power, 0.5 * power]),
                end,
                config.transport,
                || Ok(()),
            )
            .unwrap()
            .0;
        for event in trace.step().events().filter(|event| event.offset_s == end) {
            let work =
                DissipationStep::from_work(event.jump.released_energy_j, snap_partition()).unwrap();
            next = next
                .apply_heat_j(
                    &BTreeMap::from([
                        (CellId(0), work.surface_a_heat_j()),
                        (CellId(1), work.surface_b_heat_j()),
                    ]),
                    config.heat.cell_energy_tolerance_j,
                    config.heat.thermal_energy_tolerance_j,
                    || Ok(()),
                )
                .unwrap()
                .0;
        }
    }
    next
}

#[test]
fn original_smooth_losses_and_snap_heat_match_manual_native_transport_with_separate_destinations() {
    let source = source();
    let trace = collision(&source);
    let caloric = curves(&source, 1.0);
    let schedule = trace
        .step()
        .events()
        .map(|event| TimedEvent {
            time_s: event.offset_s,
            event: if event.jump.kind == JkrJumpKind::Attach {
                Event::Contact(ThermalLink {
                    id: 10,
                    a: CellId(0),
                    b: CellId(1),
                    area_m2: 1e-7,
                    conductance_w_k: release(&source) / scales(&source).3,
                })
            } else {
                Event::Separation(10)
            },
        })
        .collect();
    let thermal = thermal(&source, &caloric, 0.0, schedule, 0.0);
    let heat = binding(&source, &caloric, &thermal, &trace);
    let result = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            heat.stamp(),
            &trace,
            0.0,
            &route(),
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    let expected = manual(&source, &thermal, &trace);
    assert_eq!(result.thermal().cells(), expected.cells());
    assert_eq!(result.thermal().links(), expected.links());
    let report = result.report();
    assert_eq!(report.timed.events.len(), 2);
    assert!(report.support.work.dissipated_work_j() > 0.0);
    assert_eq!(
        report.support.outside_destination.as_deref(),
        Some("support damping reservoir")
    );
    assert_eq!(
        report.timed.snap.outside_destination.as_deref(),
        Some("surface snap waves")
    );
    close(
        report.support.work.dissipated_work_j(),
        trace.step().damping_dissipation_j,
        release(&source),
        1e-13,
    );
    close(
        report.timed.snap.work.dissipated_work_j(),
        trace.step().released_energy_j,
        release(&source),
        1e-13,
    );
    close(
        result.thermal().energy_j()
            + report.support.work.other_work_j()
            + report.timed.snap.work.other_work_j()
            - report.external_transport_heat_j,
        trace.step().damping_dissipation_j + trace.step().released_energy_j,
        release(&source),
        1e-10,
    );
    let snap_only = heat
        .propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            trace.step(),
            0.0,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    let support =
        DissipationStep::from_work(trace.step().damping_dissipation_j, route().partition())
            .unwrap();
    let lumped = snap_only
        .thermal()
        .apply_heat_j(
            &BTreeMap::from([
                (CellId(0), support.surface_a_heat_j()),
                (CellId(1), support.surface_b_heat_j()),
            ]),
            controls(&source).heat.cell_energy_tolerance_j,
            controls(&source).heat.thermal_energy_tolerance_j,
            || Ok(()),
        )
        .unwrap()
        .0;
    assert!(
        lumped.temperature_k(CellId(0)).unwrap()
            - result.thermal().temperature_k(CellId(0)).unwrap()
            > 1e-4
    );
}

#[test]
fn explicit_outside_support_loss_preserves_thermal_stepping_and_full_trace_continuation() {
    let source = source();
    let (motion, trace, load) = constant_motion(&source);
    let caloric = curves(&source, 1.0);
    let thermal = thermal(
        &source,
        &caloric,
        0.0,
        vec![],
        0.1 * release(&source) / scales(&source).3,
    );
    let heat = binding(&source, &caloric, &thermal, &trace);
    let partition = HeatPartition::new(0.0, 0.0, 1.0).unwrap();
    assert!(JkrViscousHeatRoute::new(partition, None).is_err());
    assert!(JkrViscousHeatRoute::new(partition, Some("  ".into())).is_err());
    let route =
        JkrViscousHeatRoute::new(partition, Some("unresolved support bath".into())).unwrap();
    let baseline = heat
        .propose_timed_normal_motion(
            &thermal,
            heat.stamp(),
            trace.step(),
            0.0,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    let result = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            heat.stamp(),
            &trace,
            0.0,
            &route,
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(result.thermal().cells(), baseline.thermal().cells());
    assert_eq!(&result.report().timed, baseline.report());
    assert_eq!(result.report().support.injected_surface_heat_j, [0.0; 2]);
    assert_eq!(
        result.report().support.work.other_work_j(),
        trace.step().damping_dissipation_j
    );
    let (accepted, thermal, _) = result.into_parts();
    assert!(matches!(
        accepted.propose_timed_damped_normal_motion(
            &thermal,
            accepted.stamp(),
            &trace,
            0.0,
            &route,
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
    let next = motion
        .advance_with_dissipation_trace(
            trace.step().state(),
            load,
            scales(&source).3,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(
        accepted
            .propose_timed_damped_normal_motion(
                &thermal,
                accepted.stamp(),
                &next,
                0.0,
                &route,
                controls(&source),
                || Ok(())
            )
            .is_ok()
    );
}

#[test]
fn cancellation_late_source_refusal_and_unrepresentable_source_times_publish_nothing() {
    let source = source();
    let (_, trace, _) = constant_motion(&source);
    let caloric = curves(&source, 1.0);
    let thermal = thermal(&source, &caloric, 0.0, vec![], 0.0);
    let heat = binding(&source, &caloric, &thermal, &trace);
    let stamp = heat.stamp();
    let mut polls = 0;
    let reference = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            stamp,
            &trace,
            0.0,
            &route(),
            controls(&source),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls - 1, polls] {
        let mut count = 0;
        assert!(matches!(
            heat.propose_timed_damped_normal_motion(
                &thermal,
                stamp,
                &trace,
                0.0,
                &route(),
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
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    }
    let retry = heat
        .propose_timed_damped_normal_motion(
            &thermal,
            stamp,
            &trace,
            0.0,
            &route(),
            controls(&source),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(retry.report(), reference.report());
    assert_eq!(retry.thermal().cells(), reference.thermal().cells());
    let narrow = curves(&source, 1e-3);
    let limited = self::thermal(&source, &narrow, 0.0, vec![], 0.0);
    let limited_heat = binding(&source, &narrow, &limited, &trace);
    assert!(matches!(
        limited_heat.propose_timed_damped_normal_motion(
            &limited,
            limited_heat.stamp(),
            &trace,
            0.0,
            &route(),
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(limited.time_s(), 0.0);
    assert_eq!(limited_heat.ledger().dissipated_work_j(), 0.0);
    let origin = 2.0_f64.powi(40);
    let coarse = self::thermal(&source, &caloric, origin, vec![], 0.0);
    let coarse_heat = binding(&source, &caloric, &coarse, &trace);
    assert!(origin + trace.step().state().time_s() > origin);
    assert!(matches!(
        coarse_heat.propose_timed_damped_normal_motion(
            &coarse,
            coarse_heat.stamp(),
            &trace,
            origin,
            &route(),
            controls(&source),
            || Ok(())
        ),
        Err(JkrSnapHeatError::Invalid(_))
    ));
}
