//! G0/G1/G3: original rolling release, finite calorics and atomic recontact.

#[path = "rolling/native.rs"]
mod native;

#[path = "rolling/dissipation.rs"]
mod dissipation;

use super::{close, controls, curves, source, thermal};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrJumpKind, SourceJkrSphereContact,
    rolling::{
        JkrRollingControls,
        events::{JkrRollingEventState, SourceJkrRollingEvents},
    },
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{JkrRollingSnapHeat, JkrRollingSnapHeatError, JkrRollingSnapSegment},
};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::HeatPartition;

const CRITICAL_SHIFT: f64 = 1e-8;

fn events(normal: &SourceJkrSphereContact) -> SourceJkrRollingEvents {
    SourceJkrRollingEvents::new(
        normal.clone(),
        JkrRollingControls {
            critical_rolling_displacement_m: CRITICAL_SHIFT,
            maximum_critical_displacement_to_contact_radius: 0.01,
        },
    )
    .unwrap()
}

fn loaded(events: &SourceJkrRollingEvents) -> JkrRollingEventState {
    let initial = events.initial_attached_state(0.0).unwrap();
    events
        .advance_rolling(
            &initial,
            3.0 * CRITICAL_SHIFT / events.normal_source().effective_radius_m(),
        )
        .unwrap()
        .into_parts()
        .0
}

fn binding<'s>(
    events: SourceJkrRollingEvents,
    state: JkrRollingEventState,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
) -> JkrRollingSnapHeat<'s> {
    JkrRollingSnapHeat::new(
        events,
        state,
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("resolved elastic waves".into()),
    )
    .unwrap()
}

#[test]
fn detachment_heats_original_normal_and_retained_rolling_release_once() {
    let normal = source(0.05);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let events = events(&normal);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_rolling_loss_j();
    let old_rolling_storage = initial.rolling_storage_j();
    assert!(prior_loss > 0.0);
    let fold = normal.minimum_attached_indentation_m();
    let canonical = normal
        .advance_monotone(initial.normal_response().state(), fold)
        .unwrap();
    let mut heat = binding(events, initial.clone(), &caloric, &thermal);
    let stamp = heat.stamp();
    let report = heat
        .apply(
            &mut thermal,
            stamp.clone(),
            &[JkrRollingSnapSegment::NormalIndentation(fold)],
            controls(),
            || Ok(()),
        )
        .unwrap();
    let mechanical = report.segments[0];
    assert_eq!(mechanical.normal_transition, Some(canonical));
    assert_eq!(canonical.jump.unwrap().kind, JkrJumpKind::Detach);
    close(
        report.normal_release.dissipated_work_j(),
        canonical.irreversible_work_j,
    );
    // At a_s^3 = pi W R^2/(8E*), the independent D&T stiffness is pi W.
    let rolling_release = 0.5
        * core::f64::consts::PI
        * normal.normal_state().adhesion_energy_j_per_m2()
        * CRITICAL_SHIFT.powi(2);
    close(report.rolling_release.dissipated_work_j(), rolling_release);
    close(
        mechanical.rolling_parameter_work_j,
        rolling_release - old_rolling_storage,
    );
    assert!(mechanical.rolling_parameter_work_j < 0.0);
    assert_eq!(mechanical.irreversible_rolling_loss_j, 0.0);
    close(
        report.work.dissipated_work_j(),
        canonical.irreversible_work_j + rolling_release,
    );
    close(
        report.thermal.stored_energy_change_j,
        0.75 * report.work.dissipated_work_j(),
    );
    let actual_storage_change = report.final_state.normal_response().reversible_energy_j
        + report.final_state.rolling_storage_j()
        - initial.normal_response().reversible_energy_j
        - initial.rolling_storage_j();
    let combined =
        actual_storage_change + report.thermal.stored_energy_change_j + report.work.other_work_j()
            - mechanical.applied_work_j;
    assert!(combined.abs() <= controls().transfer_energy_tolerance_j);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heat.state().normal_response().branch, JkrBranch::Detached);
    assert!(heat.state().attached().is_none());
    assert_eq!(heat.state().rolling_storage_j(), 0.0);
    assert_eq!(heat.state().accumulated_rolling_loss_j(), prior_loss);
    let accepted_cells = thermal.cells().clone();
    let accepted_state = heat.state().clone();
    assert_eq!(
        heat.apply(
            &mut thermal,
            stamp,
            &[JkrRollingSnapSegment::NormalIndentation(fold)],
            controls(),
            || Ok(())
        ),
        Err(JkrRollingSnapHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    );
    assert_eq!(heat.state(), &accepted_state);
    assert_eq!(thermal.cells(), &accepted_cells);
    let no_jump = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[JkrRollingSnapSegment::NormalIndentation(fold)],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(no_jump.work.dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &accepted_cells);
}

#[test]
fn recontact_starts_unstrained_and_preserves_separate_loss_and_release_ledgers() {
    let normal = source(0.05);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let events = events(&normal);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_rolling_loss_j();
    let fold = normal.minimum_attached_indentation_m();
    let mut heat = binding(events, initial, &caloric, &thermal);
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[
                JkrRollingSnapSegment::NormalIndentation(fold),
                JkrRollingSnapSegment::RollingAngleIncrement(0.4),
                JkrRollingSnapSegment::NormalIndentation(0.0),
            ],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.segments[1].applied_work_j, 0.0);
    assert_eq!(report.segments[1].irreversible_rolling_loss_j, 0.0);
    assert_eq!(
        report.segments[2]
            .normal_transition
            .unwrap()
            .jump
            .unwrap()
            .kind,
        JkrJumpKind::Attach
    );
    assert_eq!(report.segments[2].rolling_released_energy_j, 0.0);
    let (_, fresh) = heat.state().attached().unwrap();
    assert_eq!(fresh.elastic_shift_m(), 0.0);
    assert_eq!(heat.state().rolling_storage_j(), 0.0);
    assert_eq!(heat.state().accumulated_rolling_loss_j(), prior_loss);
    assert!(heat.rolling_release_ledger().dissipated_work_j() > 0.0);
    let release_before = *heat.rolling_release_ledger();
    let cells_before = thermal.cells().clone();
    let rolling = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[JkrRollingSnapSegment::RollingAngleIncrement(
                2.0 * CRITICAL_SHIFT / normal.effective_radius_m(),
            )],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(rolling.segments[0].irreversible_rolling_loss_j > 0.0);
    assert!(heat.state().accumulated_rolling_loss_j() > prior_loss);
    assert_eq!(rolling.work.dissipated_work_j(), 0.0);
    assert_eq!(*heat.rolling_release_ledger(), release_before);
    assert_eq!(thermal.cells(), &cells_before);
}

#[test]
fn prescribed_area_work_changes_the_spring_without_generating_snap_heat() {
    let normal = source(0.05);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let events = events(&normal);
    let state = loaded(&events);
    let old_storage = state.rolling_storage_j();
    let mut heat = binding(events, state, &caloric, &thermal);
    let cells = thermal.cells().clone();
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[JkrRollingSnapSegment::NormalIndentation(1e-6)],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(report.segments[0].rolling_parameter_work_j > 0.0);
    assert!(heat.state().rolling_storage_j() > old_storage);
    assert_eq!(report.work.dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &cells);
}

#[test]
fn source_caloric_and_late_cancellation_refusals_leave_cycle_and_cells_retryable() {
    let normal = source(0.05);
    let caloric = curves(&normal, 1.0);
    let original = thermal(&caloric, 1.0);
    let events = events(&normal);
    let initial = loaded(&events);
    let fold = normal.minimum_attached_indentation_m();
    let segments = [JkrRollingSnapSegment::NormalIndentation(fold)];
    let reference = binding(events.clone(), initial.clone(), &caloric, &original);
    let mut calls = 0;
    let accepted = reference
        .propose(&original, reference.stamp(), &segments, controls(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut thermal = original.clone();
    let mut heat = binding(events.clone(), initial.clone(), &caloric, &thermal);
    let stamp = heat.stamp();
    let mut polled = 0;
    assert_eq!(
        heat.apply(&mut thermal, stamp.clone(), &segments, controls(), || {
            polled += 1;
            if polled == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    );
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.state(), &initial);
    assert_eq!(thermal.cells(), original.cells());
    let retry = heat
        .apply(&mut thermal, stamp, &segments, controls(), || Ok(()))
        .unwrap();
    assert_eq!(&retry, accepted.report());
    assert_eq!(thermal.cells(), accepted.thermal().cells());

    // The first finite receiver admits heat before the second rejects it.
    let mut narrow = curves(&normal, 1e-5);
    narrow[0] = caloric[0].clone();
    let mut refused_thermal = self::thermal(&narrow, 1.0);
    let cells = refused_thermal.cells().clone();
    let mut heat = binding(events, initial.clone(), &narrow, &refused_thermal);
    let stamp = heat.stamp();
    assert_eq!(
        heat.apply(
            &mut refused_thermal,
            stamp.clone(),
            &segments,
            controls(),
            || Ok(())
        ),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    );
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.state(), &initial);
    assert_eq!(heat.normal_release_ledger().dissipated_work_j(), 0.0);
    assert_eq!(heat.rolling_release_ledger().dissipated_work_j(), 0.0);
    assert_eq!(refused_thermal.cells(), &cells);
    assert!(
        heat.apply(
            &mut refused_thermal,
            heat.stamp(),
            &[JkrRollingSnapSegment::NormalIndentation(f64::NAN)],
            controls(),
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(heat.state(), &initial);
    assert_eq!(refused_thermal.cells(), &cells);
}
