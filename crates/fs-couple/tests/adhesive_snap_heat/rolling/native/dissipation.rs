//! G1/G3/G4: continuous plastic heat, release and all-owner rollback.

#[path = "dissipation/timed.rs"]
mod timed;

#[path = "dissipation/temperature.rs"]
mod temperature;

use super::{DT, binding, controls, curves, events, fixture, source, thermal};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{
        JkrRollingSnapHeatError,
        dissipation::{JkrRollingDissipationHeat, native::RigidJkrRollingSnapHeatError},
    },
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::adhesive::rolling::events::prepare_rigid_adhesive_rolling_events;
use fs_mbd::{RigidBodyState, Vec3, Wrench};
use fs_thermal::enthalpy::{CellId, DepositionError};
use fs_tribo::HeatPartition;

#[test]
fn actual_native_plastic_heat_and_snap_release_close_one_energy_account() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let snap = binding(
        events.clone(),
        state.source_state().clone(),
        &caloric,
        &thermal,
    );
    let prior_loss = state.source_state().accumulated_rolling_loss_j();
    let mut heat = JkrRollingDissipationHeat::new(
        snap,
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let gate = CancelGate::new();
    let mut converted_loss = 0.0;
    for (phase, delta) in [0.0, normal.minimum_attached_indentation_m(), 0.0]
        .into_iter()
        .enumerate()
    {
        let before = f.bodies;
        let cells = thermal.cells().clone();
        let old_storage = state.source_state().normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j();
        let expected = prepare_rigid_adhesive_rolling_events(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&events),
            delta,
            DT,
            controls().mechanical,
            &gate,
        )
        .unwrap();
        let report = heat
            .apply_native(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                delta,
                DT,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(f.bodies, expected.bodies_after);
        assert_eq!(state, expected.state_after);
        assert_eq!(report.mechanics, expected.report);
        assert_eq!(heat.state(), state.source_state());
        assert_eq!(
            report.heat.rolling_loss.dissipated_work_j(),
            report.mechanics.irreversible_rolling_loss_j
        );
        converted_loss += report.heat.rolling_loss.dissipated_work_j();
        let kinetic: f64 = (0..2)
            .map(|i| {
                before[i]
                    .kinetic_energy_change_to(f.properties[i], f.bodies[i])
                    .unwrap()
            })
            .sum();
        let contact = state.source_state().normal_response().reversible_energy_j
            + state.source_state().rolling_storage_j()
            - old_storage;
        let actual_heat: f64 = [CellId(0), CellId(1)]
            .into_iter()
            .map(|id| {
                let cell = &thermal.cells()[&id];
                cell.mass_kg * (cell.specific_h_j_kg - cells[&id].specific_h_j_kg)
            })
            .sum();
        // Plastic loss has already become caloric/outside energy. Adding it
        // again here would double-count the very dissipation being transferred.
        let actual_residual = kinetic + contact + actual_heat + report.heat.work.other_work_j()
            - report.mechanics.applied_work_j;
        assert!(actual_residual.abs() <= 2e-20);
        assert!(report.combined_energy_residual_j.abs() <= 2e-20);
        assert_eq!(thermal.time_s(), 0.0);
        if phase == 0 {
            assert!(report.heat.rolling_loss.dissipated_work_j() > 0.0);
            assert_eq!(report.heat.normal_release.dissipated_work_j(), 0.0);
            assert_eq!(report.heat.rolling_release.dissipated_work_j(), 0.0);
            assert!(actual_heat > 0.0);
        } else if phase == 1 {
            assert_eq!(f.bodies, before);
            assert_eq!(report.heat.rolling_loss.dissipated_work_j(), 0.0);
            assert!(report.heat.rolling_release.dissipated_work_j() > 0.0);
            assert!(actual_heat > 0.0);
        } else {
            assert!(report.heat.normal_release.dissipated_work_j() > 0.0);
            assert_eq!(
                report
                    .mechanics
                    .mechanical
                    .unwrap()
                    .contact
                    .rolling()
                    .storage_before_j,
                0.0
            );
        }
    }
    assert_eq!(
        heat.rolling_loss_ledger().dissipated_work_j(),
        converted_loss
    );
    super::super::close(
        heat.state().accumulated_rolling_loss_j(),
        prior_loss + converted_loss,
    );
}

#[test]
fn plastic_heat_refusal_and_final_cancellation_roll_back_bodies_and_all_ledgers() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, false);
    let original_bodies = f.bodies;
    let original_state = state.clone();
    let gate = CancelGate::new();
    let mechanical = prepare_rigid_adhesive_rolling_events(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&events),
        0.0,
        DT,
        controls().mechanical,
        &gate,
    )
    .unwrap();
    assert!(mechanical.report.irreversible_rolling_loss_j > 0.0);
    assert_eq!(mechanical.report.released_energy_j, 0.0);
    let attachment_release = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let narrow = curves(
        &normal,
        mechanical.report.irreversible_rolling_loss_j / (100.0 * attachment_release),
    );
    let mut refused_thermal = thermal(&narrow, 1.0);
    let refused_cells = refused_thermal.cells().clone();
    let mut refused = JkrRollingDissipationHeat::new(
        binding(
            events.clone(),
            state.source_state().clone(),
            &narrow,
            &refused_thermal,
        ),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let refused_stamp = refused.stamp();
    assert!(matches!(
        refused.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut refused_thermal,
            refused_stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(_))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(refused.stamp(), refused_stamp);
    assert_eq!(refused.rolling_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(refused_thermal.cells(), &refused_cells);

    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let mut calls = 0;
    let proposal = heat
        .propose_native(
            &f.bodies,
            &f.properties,
            &state,
            f.binding,
            f.frame,
            &thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    assert_ne!(*proposal.bodies_after(), original_bodies);
    let final_check = calls + 1;
    calls = 0;
    assert!(matches!(
        heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || {
                calls += 1;
                if calls == final_check {
                    gate.request();
                }
                Ok(())
            }
        ),
        Err(RigidJkrRollingSnapHeatError::Heat(
            JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(DepositionError::Cancelled))
        ))
    ));
    assert_eq!(f.bodies, original_bodies);
    assert_eq!(state, original_state);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &cells);
    let gate = CancelGate::new();
    let accepted = heat
        .apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp.clone(),
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(&accepted, proposal.report());
    assert_eq!(&f.bodies, proposal.bodies_after());
    let accepted_cells = thermal.cells().clone();
    assert!(matches!(
        heat.apply_native(
            &mut f.bodies,
            &f.properties,
            &mut state,
            f.binding,
            f.frame,
            &mut thermal,
            stamp,
            0.0,
            DT,
            controls(),
            &gate,
            || Ok(())
        ),
        Err(RigidJkrRollingSnapHeatError::HistoryMismatch)
    ));
    assert_eq!(thermal.cells(), &accepted_cells);
}

#[test]
fn applied_torque_drives_rolling_from_rest_and_heats_finite_cells() {
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
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(torque, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let mut input = 0.0;
    for _ in 0..8 {
        let report = heat
            .apply_native_with_wrenches(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                0.0,
                wrenches,
                1e-3,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        let mechanical = report.mechanics.mechanical.as_ref().unwrap();
        assert!(mechanical.rolling_rate_after_rad_s > mechanical.rolling_rate_before_rad_s);
        assert!(report.mechanics.external_work_j > 0.0);
        assert_eq!(report.mechanics.external_wrenches, wrenches);
        assert_eq!(report.heat.normal_release.dissipated_work_j(), 0.0);
        assert_eq!(report.heat.rolling_release.dissipated_work_j(), 0.0);
        assert_eq!(
            report.heat.rolling_loss.dissipated_work_j(),
            report.mechanics.irreversible_rolling_loss_j
        );
        input += report.mechanics.external_work_j;
        assert!(report.combined_energy_residual_j.abs() <= 2e-20);
    }
    let kinetic: f64 = (0..2)
        .map(|i| f.bodies[i].kinetic_energy(f.properties[i]).unwrap())
        .sum();
    let actual_heat: f64 = [CellId(0), CellId(1)]
        .into_iter()
        .map(|id| {
            let c = &thermal.cells()[&id];
            c.mass_kg * (c.specific_h_j_kg - cells[&id].specific_h_j_kg)
        })
        .sum();
    assert!(actual_heat > 0.0);
    assert!(heat.rolling_loss_ledger().dissipated_work_j() > 0.0);
    let total = kinetic
        + heat.state().rolling_storage_j()
        + actual_heat
        + heat.rolling_loss_ledger().other_work_j();
    assert!((total - input).abs() <= 2e-20);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn detached_applied_work_stays_kinetic_until_fresh_recontact() {
    let normal = source(0.05);
    let events = events(&normal);
    let (mut f, mut state) = fixture(&events, true);
    let attached = events.initial_attached_state(0.0).unwrap();
    let torque = 10.0 * attached.attached().unwrap().0.law().critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(2e-8, 0.0, 0.0),
            torque_body: Vec3::new(torque, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let mut heat = JkrRollingDissipationHeat::new(
        binding(events, state.source_state().clone(), &caloric, &thermal),
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let gate = CancelGate::new();
    for delta in [-1e-6, 0.0] {
        let bodies_before = f.bodies;
        let cells = thermal.cells().clone();
        let stored_before =
            heat.state().normal_response().reversible_energy_j + heat.state().rolling_storage_j();
        let report = heat
            .apply_native_with_wrenches(
                &mut f.bodies,
                &f.properties,
                &mut state,
                f.binding,
                f.frame,
                &mut thermal,
                heat.stamp(),
                delta,
                wrenches,
                DT,
                controls(),
                &gate,
                || Ok(()),
            )
            .unwrap();
        let kinetic: f64 = (0..2)
            .map(|i| {
                bodies_before[i]
                    .kinetic_energy_change_to(f.properties[i], f.bodies[i])
                    .unwrap()
            })
            .sum();
        let actual_heat: f64 = [CellId(0), CellId(1)]
            .into_iter()
            .map(|id| {
                let c = &thermal.cells()[&id];
                c.mass_kg * (c.specific_h_j_kg - cells[&id].specific_h_j_kg)
            })
            .sum();
        let storage = heat.state().normal_response().reversible_energy_j
            + heat.state().rolling_storage_j()
            - stored_before;
        assert!(
            (kinetic + storage + actual_heat + report.heat.work.other_work_j()
                - report.mechanics.applied_work_j)
                .abs()
                <= 2e-20
        );
        assert!(report.combined_energy_residual_j.abs() <= 2e-20);
        if delta < 0.0 {
            assert_ne!(f.bodies, bodies_before);
            assert!(report.mechanics.external_work_j > 1e-20);
            assert!((kinetic - report.mechanics.external_work_j).abs() <= 1e-20);
            assert!(report.mechanics.detached_impulses.is_some());
            assert_eq!(report.heat.work.dissipated_work_j(), 0.0);
            assert_eq!(thermal.cells(), &cells);
            assert!(heat.state().attached().is_none());
            assert_eq!(heat.state().rolling_storage_j(), 0.0);
        } else {
            assert!(report.heat.normal_release.dissipated_work_j() > 0.0);
            assert!(actual_heat > 0.0);
            assert_eq!(
                report
                    .mechanics
                    .mechanical
                    .unwrap()
                    .contact
                    .rolling()
                    .storage_before_j,
                0.0
            );
        }
    }
}
