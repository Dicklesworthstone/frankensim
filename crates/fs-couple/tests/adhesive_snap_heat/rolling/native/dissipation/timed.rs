//! G1/G3/G4: native rolling, calorics and transport share physical time.

use super::{binding, controls as same_time_controls, curves, events, fixture, source, thermal};
use fs_contact::normal_patch::adhesive::{JkrBranch, JkrJumpKind};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{
        JkrRollingSnapHeatError,
        dissipation::{
            JkrRollingDissipationHeat,
            timed::{
                JkrTimedRollingHeatControls,
                native::{RigidJkrRollingIntervalHeatControls, RigidJkrRollingSnapHeatError},
            },
        },
    },
};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_mbd::tangential_friction::adhesive::rolling::events::prepare_rigid_adhesive_rolling_events_with_wrenches;
use fs_mbd::{MassProperties, RigidBodyState, Vec3, Wrench};
use fs_thermal::enthalpy::{
    DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId, StepConfig,
    TimedEvent,
};
use fs_tribo::{HeatPartition, WorkLedger};

// Dyadic endpoints make the four native transport slices exact clock intervals.
const INTERVAL: f64 = 1.0 / 1024.0;

fn controls() -> RigidJkrRollingIntervalHeatControls {
    let same_time = same_time_controls();
    RigidJkrRollingIntervalHeatControls {
        mechanical: same_time.mechanical,
        heat: JkrTimedRollingHeatControls {
            heat: same_time.heat,
            transport: StepConfig {
                max_dt_s: INTERVAL / 4.0,
                energy_tolerance_j: 2e-22,
            },
            maximum_transport_substeps: 64,
            total_energy_tolerance_j: 8e-22,
        },
    }
}

// Keep the ancestor's original material calorics, masses and cell geometry;
// supply explicit ordinary powers independently of the contact heat channel.
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
            * INTERVAL
            * (wrenches[i].force_world.dot(v0.add(v1))
                + wrenches[i].torque_body.dot(omega0.add(omega1)));
    }
    (kinetic, work)
}

fn ledgers(heat: &JkrRollingDissipationHeat<'_>) -> [WorkLedger; 3] {
    [
        *heat.normal_release_ledger(),
        *heat.rolling_release_ledger(),
        *heat.rolling_loss_ledger(),
    ]
}

#[test]
fn g1_sustained_applied_torque_and_background_heat_share_each_native_interval() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, _) = fixture(&events, true);
    f.bodies = f
        .bodies
        .map(|body| RigidBodyState::new(body.pose(), Vec3::ZERO, Vec3::ZERO).unwrap());
    let gate = CancelGate::new();
    let mut state = f
        .contact(&events)
        .initial_attached_state(&f.bodies, &f.properties, 0.0, controls().mechanical, &gate)
        .unwrap();
    let torque = 10.0
        * state
            .source_state()
            .attached()
            .unwrap()
            .0
            .law()
            .critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, -2e-8, 0.0),
            torque_body: Vec3::new(torque, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let release = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let powers = [release / (100.0 * INTERVAL), release / (200.0 * INTERVAL)];
    let background = INTERVAL * (powers[0] + powers[1]);
    let caloric = curves(&normal, 1.0);
    let mut thermal = powered_thermal(&caloric, powers);
    let initial_thermal = thermal.clone();
    let initial_state = state.clone();
    let initial_bodies = f.bodies;
    let mut heat = JkrRollingDissipationHeat::new(
        binding(
            events.clone(),
            state.source_state().clone(),
            &caloric,
            &thermal,
        ),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let mut supplied = 0.0;
    let mut outside = 0.0;
    let mut normal_work = 0.0;
    for step in 0..8 {
        let delta = if step < 4 { 0.0 } else { 1e-8 };
        let before = f.bodies;
        let prior = state.source_state().clone();
        let old_thermal = thermal.clone();
        let expected = prepare_rigid_adhesive_rolling_events_with_wrenches(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&events),
            delta,
            wrenches,
            INTERVAL,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let report = heat
            .apply_native_interval(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                delta,
                wrenches,
                INTERVAL,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(f.bodies, expected.bodies_after);
        assert_eq!(state, expected.state_after);
        assert_eq!(report.mechanics, expected.report);
        assert_eq!(heat.state(), state.source_state());
        assert_eq!(report.heat.initial_time_s, step as f64 * INTERVAL);
        assert_eq!(thermal.time_s(), (step + 1) as f64 * INTERVAL);
        assert_eq!(report.heat.final_time_s, thermal.time_s());
        assert_eq!(report.heat.transport.substeps, 4);
        assert_eq!(report.heat.thermal.incoming_energy_j, 0.0);
        let mechanical = report.mechanics.mechanical.as_ref().unwrap();
        assert!(mechanical.rolling_rate_after_rad_s > mechanical.rolling_rate_before_rad_s);
        assert_eq!(
            report.heat.rolling_loss.dissipated_work_j(),
            report.mechanics.irreversible_rolling_loss_j
        );
        assert_eq!(report.heat.normal.work.dissipated_work_j(), 0.0);
        let (kinetic, body_work) = body_changes(before, f.bodies, f.properties, wrenches);
        let contact = state.source_state().normal_response().reversible_energy_j
            - prior.normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j()
            - prior.rolling_storage_j();
        let prescribed = report.mechanics.normal.normal_contact_work_j
            + report.mechanics.normal.rolling_parameter_work_j;
        let escaped =
            report.heat.normal.work.other_work_j() + report.heat.rolling_loss.other_work_j();
        let residual = kinetic + contact + cell_change(&old_thermal, &thermal) + escaped
            - body_work
            - prescribed
            - background;
        assert!(
            residual.abs() <= 2e-20,
            "independent interval energy: {residual:e}"
        );
        assert!((report.mechanics.external_work_j - body_work).abs() <= 1e-20);
        assert!((report.heat.external_transport_heat_j - background).abs() <= 1e-22);
        assert!(report.combined_energy_residual_j.abs() <= 2e-20);
        supplied += body_work + prescribed + background;
        normal_work += prescribed;
        outside += escaped;
    }
    let kinetic: f64 = (0..2)
        .map(|i| {
            initial_bodies[i]
                .kinetic_energy_change_to(f.properties[i], f.bodies[i])
                .unwrap()
        })
        .sum();
    let contact = state.source_state().normal_response().reversible_energy_j
        - initial_state
            .source_state()
            .normal_response()
            .reversible_energy_j
        + state.source_state().rolling_storage_j()
        - initial_state.source_state().rolling_storage_j();
    assert!(
        (kinetic + contact + cell_change(&initial_thermal, &thermal) + outside - supplied).abs()
            <= 2e-20
    );
    assert!(normal_work != 0.0);
    assert!(heat.rolling_loss_ledger().dissipated_work_j() > 0.0);
    assert!(outside > 0.0);
    assert_eq!(thermal.time_s(), 8.0 * INTERVAL);
}

#[test]
fn g3_separation_and_fresh_recontact_preserve_release_channels_as_real_time_advances() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let torque = 2.0
        * state
            .source_state()
            .attached()
            .unwrap()
            .0
            .law()
            .critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(2e-8, 0.0, 0.0),
            torque_body: Vec3::new(torque, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let release = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let powers = [release / (100.0 * INTERVAL), 0.0];
    let caloric = curves(&normal, 1.0);
    let mut thermal = powered_thermal(&caloric, powers);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(
            events.clone(),
            state.source_state().clone(),
            &caloric,
            &thermal,
        ),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let initial_loss = state.source_state().accumulated_rolling_loss_j();
    let gate = CancelGate::new();
    let fold = normal.minimum_attached_indentation_m();
    let mut spring_release = 0.0;
    for (phase, delta) in [fold, fold - 1e-8, 0.0].into_iter().enumerate() {
        let before = f.bodies;
        let prior = state.source_state().clone();
        let old_thermal = thermal.clone();
        let normal_stage = events.advance_normal(&prior, delta).unwrap();
        let expected = prepare_rigid_adhesive_rolling_events_with_wrenches(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&events),
            delta,
            wrenches,
            INTERVAL,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let report = heat
            .apply_native_interval(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                delta,
                wrenches,
                INTERVAL,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(report.mechanics, expected.report);
        assert_eq!(f.bodies, expected.bodies_after);
        assert_eq!(state, expected.state_after);
        assert_eq!(heat.state(), state.source_state());
        assert_eq!(report.heat.initial_time_s, phase as f64 * INTERVAL);
        assert_eq!(thermal.time_s(), (phase + 1) as f64 * INTERVAL);
        assert_eq!(report.heat.normal.segments[0], *normal_stage.report());
        assert_eq!(report.heat.normal.thermal.substeps, 0);
        assert_eq!(report.heat.transport.substeps, 4);
        assert_eq!(
            report.heat.normal.normal_release.dissipated_work_j(),
            normal_stage.report().normal_released_energy_j
        );
        assert_eq!(
            report.heat.normal.rolling_release.dissipated_work_j(),
            normal_stage.report().rolling_released_energy_j
        );
        let (kinetic, work) = body_changes(before, f.bodies, f.properties, wrenches);
        let contact = state.source_state().normal_response().reversible_energy_j
            - prior.normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j()
            - prior.rolling_storage_j();
        let escaped =
            report.heat.normal.work.other_work_j() + report.heat.rolling_loss.other_work_j();
        let prescribed = report.mechanics.normal.normal_contact_work_j
            + report.mechanics.normal.rolling_parameter_work_j;
        let residual = kinetic + contact + cell_change(&old_thermal, &thermal) + escaped
            - work
            - prescribed
            - INTERVAL * powers[0];
        assert!(residual.abs() <= 2e-20);
        assert_ne!(before, f.bodies);
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
            assert!(report.heat.normal.normal_release.dissipated_work_j() > 0.0);
            spring_release = report.heat.normal.rolling_release.dissipated_work_j();
            assert!(spring_release > 0.0);
            assert!(report.heat.normal.thermal.stored_energy_change_j > 0.0);
        } else if phase == 1 {
            assert_eq!(report.heat.normal.work.dissipated_work_j(), 0.0);
            assert!((kinetic - work).abs() <= 1e-20);
            assert_eq!(
                heat.rolling_release_ledger().dissipated_work_j(),
                spring_release
            );
        } else {
            let mechanical = report.mechanics.mechanical.as_ref().unwrap();
            assert_eq!(mechanical.contact.rolling().storage_before_j, 0.0);
            assert!(report.heat.normal.normal_release.dissipated_work_j() > 0.0);
            assert_eq!(report.heat.normal.rolling_release.dissipated_work_j(), 0.0);
            assert_eq!(
                heat.rolling_release_ledger().dissipated_work_j(),
                spring_release
            );
            assert!(report.heat.rolling_loss.dissipated_work_j() > 0.0);
        }
        if phase < 2 {
            assert_eq!(
                state.source_state().normal_response().branch,
                JkrBranch::Detached
            );
            assert!(report.mechanics.mechanical.is_none());
            assert!(report.mechanics.detached_impulses.is_some());
            assert_eq!(report.heat.rolling_loss.dissipated_work_j(), 0.0);
            assert_eq!(report.heat.injected_surface_heat_j, [0.0; 2]);
            assert_eq!(
                state.source_state().accumulated_rolling_loss_j(),
                initial_loss
            );
        }
    }
}

#[test]
fn g4_late_transport_cancellation_and_stale_inputs_leave_all_owners_for_exact_retry() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let original_bodies = f.bodies;
    let original_state = state.clone();
    let caloric = curves(&normal, 1.0);
    let mut thermal = powered_thermal(&caloric, [1e-12, 2e-12]);
    let original_thermal = thermal.clone();
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let original_ledgers = ledgers(&heat);
    let target = normal.minimum_attached_indentation_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, 0.0, 0.0),
            torque_body: Vec3::new(1e-12, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let gate = CancelGate::new();
    let mut bounded = controls();
    bounded.heat.maximum_transport_substeps = 1;
    // Original normal/spring heat and one transport slice have been proposed
    // before the remaining three slices exceed this native transport budget.
    assert!(matches!(
        heat.apply_native_interval(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            target,
            wrenches,
            INTERVAL,
            bounded,
            &gate,
            || Ok(()),
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(DepositionError::Budget(_)))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    assert_eq!(thermal.cells(), original_thermal.cells());
    assert_eq!(thermal.links(), original_thermal.links());
    assert_eq!(thermal.time_s(), original_thermal.time_s());
    let mut proposal_calls = 0;
    let proposal = heat
        .propose_native_interval(
            &f.bodies,
            &f.properties,
            &state,
            f.binding,
            f.frame,
            &thermal,
            stamp.clone(),
            target,
            wrenches,
            INTERVAL,
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
            .rolling_release
            .dissipated_work_j()
            > 0.0
    );
    assert_eq!(proposal.thermal().time_s(), INTERVAL);
    for request_gate in [false, true] {
        let gate = CancelGate::new();
        let mut calls = 0;
        assert!(matches!(
            heat.apply_native_interval(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                stamp.clone(),
                target,
                wrenches,
                INTERVAL,
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
                },
            ),
            Err(RigidJkrRollingSnapHeatError::Heat(
                JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
                    DepositionError::Cancelled
                ))
            ))
        ));
        assert_eq!(calls, proposal_calls + 1);
        assert_eq!(f.bodies, original_bodies);
        assert_eq!(state, original_state);
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(ledgers(&heat), original_ledgers);
        assert_eq!(thermal.cells(), original_thermal.cells());
        assert_eq!(thermal.links(), original_thermal.links());
        assert_eq!(thermal.time_s(), 0.0);
    }
    let mut advanced = thermal
        .advance_to(INTERVAL / 4.0, controls().heat.transport, || Ok(()))
        .unwrap()
        .0;
    let advanced_cells = advanced.cells().clone();
    assert!(matches!(
        heat.apply_native_interval(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut advanced,
            stamp.clone(),
            target,
            wrenches,
            INTERVAL,
            controls(),
            &gate,
            || Ok(()),
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::HistoryMismatch)
        ))
    ));
    assert_eq!(advanced.cells(), &advanced_cells);
    assert_eq!(advanced.time_s(), INTERVAL / 4.0);
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), original_ledgers);
    let report = heat
        .apply_native_interval(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            target,
            wrenches,
            INTERVAL,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(&report, proposal.report());
    assert_eq!(&f.bodies, proposal.bodies_after());
    assert_eq!(&state, proposal.native_state());
    assert_eq!(heat.stamp(), proposal.history().stamp());
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    assert_eq!(thermal.links(), proposal.thermal().links());
    assert_eq!(thermal.time_s(), proposal.thermal().time_s());
    let accepted_stamp = heat.stamp();
    assert!(matches!(
        heat.apply_native_interval(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp,
            target,
            wrenches,
            INTERVAL,
            controls(),
            &gate,
            || Ok(()),
        ),
        Err(RigidJkrRollingSnapHeatError::HistoryMismatch)
    ));
    assert_eq!(&f.bodies, proposal.bodies_after());
    assert_eq!(&state, proposal.native_state());
    assert_eq!(heat.stamp(), accepted_stamp);
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    assert_eq!(thermal.links(), proposal.thermal().links());
    assert_eq!(thermal.time_s(), INTERVAL);
}
