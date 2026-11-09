//! G1/G3/G4: actual finite-body twist and calorics share physical time.

#[path = "native/temperature.rs"]
mod temperature;

use super::{
    CRITICAL_TWIST, binding, controls as heat_controls, curves, events, ledgers, source, thermal,
};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrJumpKind, torsion::events::SourceJkrTorsionEvents,
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    torsion::{
        JkrTorsionHeat, JkrTorsionHeatError, JkrTorsionHeatStamp,
        timed::{
            JkrTimedTorsionHeatControls,
            native::{
                RigidJkrTorsionIntervalHeatControls, RigidJkrTorsionIntervalHeatError,
                RigidJkrTorsionIntervalHeatReport,
            },
        },
    },
};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_mbd::tangential_friction::{
    RigidContactBinding,
    adhesive::torsion::{
        RigidAdhesiveTorsionControls,
        events::{
            RigidAdhesiveTorsionEventContact, RigidAdhesiveTorsionEventState,
            apply_rigid_adhesive_torsion_events,
            prepare_rigid_adhesive_torsion_events_with_wrenches,
        },
    },
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId, StepConfig,
    TimedEvent,
};
use fs_tribo::partial_slip::TangentFrame;

const DT: f64 = 1.0 / 1024.0;

fn controls() -> RigidJkrTorsionIntervalHeatControls {
    RigidJkrTorsionIntervalHeatControls {
        mechanical: RigidAdhesiveTorsionControls {
            maximum_bodies: 2,
            maximum_step_s: 0.01,
            maximum_angle_increment_rad: 1e-3,
            maximum_iterations: 96,
            angle_tolerance_rad: 1e-18,
            geometry_tolerance_m: 1e-15,
            angular_velocity_tolerance_rad_s: 1e-10,
            linear_momentum_tolerance_ns: 1e-25,
            angular_momentum_tolerance_nms: 1e-25,
            energy_tolerance_j: 1e-20,
        },
        heat: JkrTimedTorsionHeatControls {
            heat: heat_controls(),
            transport: StepConfig {
                max_dt_s: DT / 4.0,
                energy_tolerance_j: 2e-22,
            },
            maximum_transport_substeps: 64,
            total_energy_tolerance_j: 8e-22,
        },
    }
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
    state: RigidAdhesiveTorsionEventState,
}
impl Fixture {
    fn new(source: &SourceJkrTorsionEvents, delta: f64, rate: f64) -> Self {
        let properties = [
            (2e-5, Vec3::new(0.8e-9, 1e-9, 1.2e-9)),
            (3e-5, Vec3::new(1.5e-9, 2e-9, 2.5e-9)),
        ]
        .map(|(mass, inertia)| {
            MassProperties::from_principal_frame(
                mass,
                Vec3::ZERO,
                inertia,
                UnitQuaternion::IDENTITY,
            )
            .unwrap()
        });
        let arms = [Vec3::new(0.0, 0.0, 0.01), Vec3::new(0.0, 0.0, -0.01)];
        let bodies = core::array::from_fn(|i| {
            let omega = Vec3::new(0.0, 0.0, if i == 0 { rate } else { 0.0 });
            RigidBodyState::new(
                Pose::new(arms[i].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                properties[i].angular_momentum_body_checked(omega).unwrap(),
            )
            .unwrap()
        });
        let binding = RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        };
        let frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
        let state = RigidAdhesiveTorsionEventContact {
            binding,
            frame,
            source,
        }
        .initial_attached_state(
            &bodies,
            &properties,
            delta,
            controls().mechanical,
            &CancelGate::new(),
        )
        .unwrap();
        Self {
            bodies,
            properties,
            binding,
            frame,
            state,
        }
    }

    fn contact<'a>(
        &self,
        source: &'a SourceJkrTorsionEvents,
    ) -> RigidAdhesiveTorsionEventContact<'a> {
        RigidAdhesiveTorsionEventContact {
            binding: self.binding,
            frame: self.frame,
            source,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn apply(
        &mut self,
        heat: &mut JkrTorsionHeat<'_>,
        thermal: &mut DepositionSimulation<'_>,
        stamp: JkrTorsionHeatStamp,
        delta: f64,
        wrenches: [Wrench; 2],
        limits: RigidJkrTorsionIntervalHeatControls,
        gate: &CancelGate,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrTorsionIntervalHeatReport, RigidJkrTorsionIntervalHeatError> {
        heat.apply_native_interval(
            &mut self.bodies,
            &self.properties,
            &mut self.state,
            self.binding,
            self.frame,
            thermal,
            stamp,
            delta,
            wrenches,
            DT,
            limits,
            gate,
            checkpoint,
        )
    }
}

fn loaded_native(events: &SourceJkrTorsionEvents) -> Fixture {
    let delta = events
        .normal_source()
        .evaluate_force_controlled(0.0)
        .unwrap()
        .indentation_m;
    let mut f = Fixture::new(events, delta, 3.0 * CRITICAL_TWIST / DT);
    let contact = f.contact(events);
    apply_rigid_adhesive_torsion_events(
        &mut f.bodies,
        &f.properties,
        &mut f.state,
        contact,
        delta,
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(f.state.source_state().torsion_storage_j() > 0.0);
    assert!(f.state.source_state().accumulated_torsion_loss_j() > 0.0);
    f
}

fn powered_thermal(
    caloric: &[EquilibriumEnthalpyPhaseCurve; 2],
    powers: [f64; 2],
) -> DepositionSimulation<'_> {
    let initial = thermal(caloric, 1.0);
    let schedule = initial
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.boundary.source_w = powers[cell.id.0 as usize];
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_substeps: 64,
            ..initial.budget()
        },
    )
    .unwrap()
    .advance_to(0.0, controls().heat.transport, || Ok(()))
    .unwrap()
    .0
}

fn cell_change(before: &DepositionSimulation<'_>, after: &DepositionSimulation<'_>) -> f64 {
    after
        .cells()
        .iter()
        .map(|(id, cell)| {
            cell.mass_kg * (cell.specific_h_j_kg - before.cells()[id].specific_h_j_kg)
        })
        .sum()
}

fn body_changes(
    before: [RigidBodyState; 2],
    after: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    wrenches: [Wrench; 2],
) -> (f64, f64) {
    let mut kinetic = 0.0;
    let mut work = 0.0;
    for i in 0..2 {
        kinetic += before[i]
            .kinetic_energy_change_to(properties[i], after[i])
            .unwrap();
        let v0 = before[i]
            .center_of_mass_velocity_world(properties[i])
            .unwrap();
        let v1 = after[i]
            .center_of_mass_velocity_world(properties[i])
            .unwrap();
        let omega0 = properties[i]
            .angular_velocity_body_checked(before[i].angular_momentum_body())
            .unwrap();
        let omega1 = properties[i]
            .angular_velocity_body_checked(after[i].angular_momentum_body())
            .unwrap();
        work += 0.5
            * DT
            * (wrenches[i].force_world.dot(v0.add(v1))
                + wrenches[i].torque_body.dot(omega0.add(omega1)));
    }
    (kinetic, work)
}

#[test]
fn g1_driven_twist_and_background_heat_close_each_shared_physical_interval() {
    let normal = source(0.05);
    let events = events(&normal);
    let mut f = Fixture::new(&events, 0.0, 0.0);
    let torque = 10.0
        * f.state
            .source_state()
            .attached()
            .unwrap()
            .0
            .law()
            .critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, -2e-8, 0.0),
            torque_body: Vec3::new(0.0, 0.0, torque),
        },
        Wrench::ZERO,
    ];
    let release = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let powers = [release / (100.0 * DT), release / (200.0 * DT)];
    let background = DT * (powers[0] + powers[1]);
    let caloric = curves(&normal, 1.0);
    let mut thermal = powered_thermal(&caloric, powers);
    let initial_thermal = thermal.clone();
    let initial_bodies = f.bodies;
    let initial_state = f.state.clone();
    let mut heat = binding(
        events.clone(),
        f.state.source_state().clone(),
        &caloric,
        &thermal,
    );
    let gate = CancelGate::new();
    let mut supplied = 0.0;
    let mut outside = 0.0;
    let mut normal_work = 0.0;
    for step in 0..8 {
        let delta = if step < 4 { 0.0 } else { 1e-8 };
        let before = f.bodies;
        let prior = f.state.source_state().clone();
        let old_thermal = thermal.clone();
        let direct = prepare_rigid_adhesive_torsion_events_with_wrenches(
            &f.bodies,
            &f.properties,
            &f.state,
            f.contact(&events),
            delta,
            wrenches,
            DT,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let stamp = heat.stamp();
        let report = f
            .apply(
                &mut heat,
                &mut thermal,
                stamp,
                delta,
                wrenches,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(f.bodies, direct.bodies_after);
        assert_eq!(f.state, direct.state_after);
        assert_eq!(report.mechanics, direct.report);
        assert_eq!(heat.state(), f.state.source_state());
        assert_eq!(report.heat.initial_time_s, step as f64 * DT);
        assert_eq!(thermal.time_s(), (step + 1) as f64 * DT);
        assert_eq!(report.heat.transport.substeps, 4);
        assert_eq!(report.heat.thermal.incoming_energy_j, 0.0);
        assert_eq!(report.heat.normal.work.dissipated_work_j(), 0.0);
        assert_eq!(
            report.heat.torsion_loss.dissipated_work_j(),
            report.mechanics.irreversible_torsion_loss_j
        );
        let (kinetic, work) = body_changes(before, f.bodies, f.properties, wrenches);
        let contact = f.state.source_state().normal_response().reversible_energy_j
            - prior.normal_response().reversible_energy_j
            + f.state.source_state().torsion_storage_j()
            - prior.torsion_storage_j();
        let prescribed = report.mechanics.normal.normal_contact_work_j
            + report.mechanics.normal.torsion_parameter_work_j;
        let escaped =
            report.heat.normal.work.other_work_j() + report.heat.torsion_loss.other_work_j();
        let residual = kinetic + contact + cell_change(&old_thermal, &thermal) + escaped
            - work
            - prescribed
            - background;
        assert!(
            residual.abs() <= 2e-20,
            "independent interval energy: {residual:e}"
        );
        assert!((report.mechanics.external_work_j - work).abs() <= 1e-20);
        assert!((report.heat.external_transport_heat_j - background).abs() <= 1e-22);
        assert!(report.combined_energy_residual_j.abs() <= 2e-20);
        supplied += work + prescribed + background;
        outside += escaped;
        normal_work += prescribed;
    }
    let kinetic: f64 = (0..2)
        .map(|i| {
            initial_bodies[i]
                .kinetic_energy_change_to(f.properties[i], f.bodies[i])
                .unwrap()
        })
        .sum();
    let contact = f.state.source_state().normal_response().reversible_energy_j
        - initial_state
            .source_state()
            .normal_response()
            .reversible_energy_j
        + f.state.source_state().torsion_storage_j()
        - initial_state.source_state().torsion_storage_j();
    assert!(
        (kinetic + contact + cell_change(&initial_thermal, &thermal) + outside - supplied).abs()
            <= 2e-20
    );
    assert!(normal_work != 0.0);
    assert!(heat.torsion_loss_ledger().dissipated_work_j() > 0.0);
    assert!(outside > 0.0);
    assert_eq!(thermal.time_s(), 8.0 * DT);
}

#[test]
fn g3_fold_release_detached_work_and_fresh_twist_follow_the_thermal_clock() {
    let normal = source(0.05);
    let events = events(&normal);
    let mut f = loaded_native(&events);
    let expected_fold_twist = f.state.source_state().torsion_storage_j() / 36.0;
    let prior_loss = f.state.source_state().accumulated_torsion_loss_j();
    let torque = 2.0
        * f.state
            .source_state()
            .attached()
            .unwrap()
            .0
            .law()
            .critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(2e-8, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, torque),
        },
        Wrench::ZERO,
    ];
    let caloric = curves(&normal, 1.0);
    let powers = [1e-12, 2e-12];
    let mut thermal = powered_thermal(&caloric, powers);
    let mut heat = binding(
        events.clone(),
        f.state.source_state().clone(),
        &caloric,
        &thermal,
    );
    let gate = CancelGate::new();
    let fold = normal.minimum_attached_indentation_m();
    for (phase, delta) in [fold, fold - 1e-8, 0.0].into_iter().enumerate() {
        let before = f.bodies;
        let prior = f.state.source_state().clone();
        let old_thermal = thermal.clone();
        let normal_stage = events.advance_normal(&prior, delta).unwrap();
        let direct = prepare_rigid_adhesive_torsion_events_with_wrenches(
            &f.bodies,
            &f.properties,
            &f.state,
            f.contact(&events),
            delta,
            wrenches,
            DT,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let stamp = heat.stamp();
        let report = f
            .apply(
                &mut heat,
                &mut thermal,
                stamp,
                delta,
                wrenches,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(report.mechanics, direct.report);
        assert_eq!(f.bodies, direct.bodies_after);
        assert_eq!(f.state, direct.state_after);
        assert_eq!(heat.state(), f.state.source_state());
        assert_eq!(report.heat.initial_time_s, phase as f64 * DT);
        assert_eq!(thermal.time_s(), (phase + 1) as f64 * DT);
        assert_eq!(report.heat.normal.segments[0], *normal_stage.report());
        assert_eq!(report.heat.normal.thermal.substeps, 0);
        assert_eq!(report.heat.transport.substeps, 4);
        let (kinetic, work) = body_changes(before, f.bodies, f.properties, wrenches);
        let contact = f.state.source_state().normal_response().reversible_energy_j
            - prior.normal_response().reversible_energy_j
            + f.state.source_state().torsion_storage_j()
            - prior.torsion_storage_j();
        let escaped =
            report.heat.normal.work.other_work_j() + report.heat.torsion_loss.other_work_j();
        let prescribed = report.mechanics.normal.normal_contact_work_j
            + report.mechanics.normal.torsion_parameter_work_j;
        assert!(
            (kinetic + contact + cell_change(&old_thermal, &thermal) + escaped
                - work
                - prescribed
                - DT * (powers[0] + powers[1]))
                .abs()
                <= 2e-20
        );
        if phase == 0 {
            let jump = report
                .mechanics
                .normal
                .normal_transition
                .unwrap()
                .jump
                .unwrap();
            assert_eq!(jump.kind, JkrJumpKind::Detach);
            assert_eq!(jump.indentation_m, fold);
            super::close(
                report.heat.normal.torsion_release.dissipated_work_j(),
                expected_fold_twist,
            );
            assert!(report.heat.normal.normal_release.dissipated_work_j() > 0.0);
        } else if phase == 1 {
            assert_eq!(report.heat.normal.work.dissipated_work_j(), 0.0);
        } else {
            let mechanical = report.mechanics.mechanical.as_ref().unwrap();
            assert_eq!(mechanical.contact.torsion().storage_before_j, 0.0);
            assert!(report.heat.normal.normal_release.dissipated_work_j() > 0.0);
            assert_eq!(report.heat.normal.torsion_release.dissipated_work_j(), 0.0);
            assert!(report.heat.torsion_loss.dissipated_work_j() > 0.0);
        }
        super::close(
            heat.torsion_release_ledger().dissipated_work_j(),
            expected_fold_twist,
        );
        if phase < 2 {
            assert_eq!(
                f.state.source_state().normal_response().branch,
                JkrBranch::Detached
            );
            assert!(report.mechanics.mechanical.is_none());
            assert!(report.mechanics.detached_impulses.is_some());
            assert_eq!(report.heat.torsion_loss.dissipated_work_j(), 0.0);
            assert_eq!(report.heat.injected_surface_heat_j, [0.0; 2]);
            assert_eq!(
                f.state.source_state().accumulated_torsion_loss_j(),
                prior_loss
            );
            assert!((kinetic - work).abs() <= 1e-20);
            // The snap itself adds no angular impulse; only actual body loads do.
            for i in 0..2 {
                assert_eq!(
                    f.bodies[i].angular_momentum_body(),
                    before[i]
                        .angular_momentum_body()
                        .add(wrenches[i].torque_body.scale(DT))
                );
            }
        }
    }
    super::close(
        f.state.source_state().accumulated_torsion_loss_j(),
        prior_loss + heat.torsion_loss_ledger().dissipated_work_j(),
    );
}

#[test]
fn g4_late_thermal_refusal_cancellation_and_stale_inputs_preserve_every_owner() {
    let normal = source(0.05);
    let events = events(&normal);
    let mut f = loaded_native(&events);
    let original_bodies = f.bodies;
    let original_state = f.state.clone();
    let caloric = curves(&normal, 1.0);
    let mut thermal = powered_thermal(&caloric, [1e-12, 2e-12]);
    let original_thermal = thermal.clone();
    let mut heat = binding(events, f.state.source_state().clone(), &caloric, &thermal);
    let stamp = heat.stamp();
    let original_ledgers = ledgers(&heat);
    let target = normal.minimum_attached_indentation_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 1e-12),
        },
        Wrench::ZERO,
    ];
    let gate = CancelGate::new();
    let mut bounded = controls();
    bounded.heat.maximum_transport_substeps = 1;
    // Release and a transport slice are private when the remaining work refuses.
    assert!(matches!(
        f.apply(
            &mut heat,
            &mut thermal,
            stamp.clone(),
            target,
            wrenches,
            bounded,
            &gate,
            || Ok(())
        ),
        Err(RigidJkrTorsionIntervalHeatError::Heat(
            JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(DepositionError::Budget(_)))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(f.state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    assert_eq!(thermal.cells(), original_thermal.cells());
    assert_eq!(thermal.links(), original_thermal.links());
    assert_eq!(thermal.time_s(), 0.0);

    let mut proposal_calls = 0;
    let proposal = heat
        .propose_native_interval(
            &f.bodies,
            &f.properties,
            &f.state,
            f.binding,
            f.frame,
            &thermal,
            stamp.clone(),
            target,
            wrenches,
            DT,
            controls(),
            &gate,
            || {
                proposal_calls += 1;
                Ok(())
            },
        )
        .unwrap();
    assert!(
        proposal
            .report()
            .heat
            .normal
            .torsion_release
            .dissipated_work_j()
            > 0.0
    );
    assert_eq!(proposal.thermal().time_s(), DT);
    for request_gate in [false, true] {
        let gate = CancelGate::new();
        let mut calls = 0;
        assert!(matches!(
            f.apply(
                &mut heat,
                &mut thermal,
                stamp.clone(),
                target,
                wrenches,
                controls(),
                &gate,
                || {
                    calls += 1;
                    if calls == proposal_calls + 1 {
                        if request_gate {
                            gate.request();
                        } else {
                            return Err(DepositionError::Cancelled);
                        }
                    }
                    Ok(())
                }
            ),
            Err(RigidJkrTorsionIntervalHeatError::Heat(
                JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(DepositionError::Cancelled))
            ))
        ));
        assert_eq!(calls, proposal_calls + 1);
        assert_eq!(f.bodies, original_bodies);
        assert_eq!(f.state, original_state);
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(ledgers(&heat), original_ledgers);
        assert_eq!(thermal.cells(), original_thermal.cells());
        assert_eq!(thermal.links(), original_thermal.links());
        assert_eq!(thermal.time_s(), 0.0);
    }
    let mut advanced = thermal
        .advance_to(DT / 4.0, controls().heat.transport, || Ok(()))
        .unwrap()
        .0;
    let advanced_cells = advanced.cells().clone();
    assert!(matches!(
        f.apply(
            &mut heat,
            &mut advanced,
            stamp.clone(),
            target,
            wrenches,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrTorsionIntervalHeatError::Heat(
            JkrTorsionHeatError::Heat(JkrSnapHeatError::HistoryMismatch)
        ))
    ));
    assert_eq!(advanced.cells(), &advanced_cells);
    assert_eq!(advanced.time_s(), DT / 4.0);
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(f.state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    let report = f
        .apply(
            &mut heat,
            &mut thermal,
            stamp.clone(),
            target,
            wrenches,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(&report, proposal.report());
    assert_eq!(&f.bodies, proposal.bodies_after());
    assert_eq!(&f.state, proposal.native_state());
    assert_eq!(heat.stamp(), proposal.history().stamp());
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    assert_eq!(thermal.links(), proposal.thermal().links());
    assert_eq!(thermal.time_s(), proposal.thermal().time_s());
    let accepted_stamp = heat.stamp();
    assert!(matches!(
        f.apply(
            &mut heat,
            &mut thermal,
            stamp,
            target,
            wrenches,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrTorsionIntervalHeatError::HistoryMismatch)
    ));
    assert_eq!(&f.bodies, proposal.bodies_after());
    assert_eq!(&f.state, proposal.native_state());
    assert_eq!(heat.stamp(), accepted_stamp);
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    assert_eq!(thermal.time_s(), DT);
}
