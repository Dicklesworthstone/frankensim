//! G1/G3: plastic rolling heat, separate spring release and finite calorics.

#[path = "dissipation/timed.rs"]
mod timed;

use super::{CRITICAL_SHIFT, binding, close, controls, curves, events, loaded, source, thermal};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    rolling::{
        JkrRollingSnapHeat, JkrRollingSnapHeatError, JkrRollingSnapSegment as Segment,
        dissipation::JkrRollingDissipationHeat,
    },
};
use fs_thermal::enthalpy::CellId;
use fs_tribo::HeatPartition;

#[test]
fn elastic_rolling_and_prescribed_area_work_do_not_heat_the_receivers() {
    let normal = source(0.05);
    let events = events(&normal);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let initial = events.initial_attached_state(0.0).unwrap();
    let snap = binding(events, initial, &caloric, &thermal);
    let mut heat = JkrRollingDissipationHeat::new(
        snap,
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let theta = CRITICAL_SHIFT / normal.effective_radius_m();
    let cells = thermal.cells().clone();
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[
                Segment::RollingAngleIncrement(0.5 * theta),
                Segment::NormalIndentation(1e-6),
                Segment::RollingAngleIncrement(-0.5 * theta),
            ],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(report.segments[0].storage_change_j > 0.0);
    assert!(report.segments[1].rolling_parameter_work_j > 0.0);
    assert!(report.segments[2].storage_change_j < 0.0);
    assert_eq!(report.rolling_loss.dissipated_work_j(), 0.0);
    assert_eq!(report.work.dissipated_work_j(), 0.0);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &cells);
}

#[test]
fn plastic_loss_release_and_recontact_keep_distinct_destinations_and_history() {
    let normal = source(0.05);
    let events = events(&normal);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_rolling_loss_j();
    assert!(prior_loss > 0.0);
    let old_storage = initial.normal_response().reversible_energy_j + initial.rolling_storage_j();
    let snap = binding(events, initial, &caloric, &thermal);
    let mut heat = JkrRollingDissipationHeat::new(
        snap,
        HeatPartition::new(0.8, 0.1, 0.1).unwrap(),
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    assert_eq!(heat.state().accumulated_rolling_loss_j(), prior_loss);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    let theta = CRITICAL_SHIFT / normal.effective_radius_m();
    let path = [
        Segment::RollingAngleIncrement(2.0 * theta),
        Segment::NormalIndentation(normal.minimum_attached_indentation_m()),
        Segment::RollingAngleIncrement(0.3),
        Segment::NormalIndentation(0.0),
        Segment::RollingAngleIncrement(3.0 * theta),
    ];
    let cells = thermal.cells().clone();
    let stamp = heat.stamp();
    let report = heat
        .apply(&mut thermal, stamp.clone(), &path, controls(), || Ok(()))
        .unwrap();
    let plastic = report.rolling_loss.dissipated_work_j();
    let release =
        report.normal_release.dissipated_work_j() + report.rolling_release.dissipated_work_j();
    assert!(plastic > 0.0);
    assert!(report.normal_release.dissipated_work_j() > 0.0);
    assert!(report.rolling_release.dissipated_work_j() > 0.0);
    assert_eq!(report.segments[2].applied_work_j, 0.0);
    assert_eq!(report.segments[2].irreversible_rolling_loss_j, 0.0);
    assert_eq!(
        report.segments[4]
            .rolling
            .unwrap()
            .rolling()
            .storage_before_j,
        0.0
    );
    close(
        heat.state().accumulated_rolling_loss_j(),
        prior_loss + plastic,
    );
    close(heat.rolling_loss_ledger().dissipated_work_j(), plastic);
    close(
        report.work.surface_a_heat_j(),
        0.25 * release + 0.8 * plastic,
    );
    close(
        report.work.surface_b_heat_j(),
        0.5 * release + 0.1 * plastic,
    );
    close(report.work.other_work_j(), 0.25 * release + 0.1 * plastic);
    assert_eq!(
        report.release_outside_destination.as_deref(),
        Some("resolved elastic waves")
    );
    assert_eq!(
        report.rolling_outside_destination.as_deref(),
        Some("rolling microstructural rearrangement")
    );
    let actual_heat: f64 = [CellId(0), CellId(1)]
        .into_iter()
        .map(|id| {
            let c = &thermal.cells()[&id];
            c.mass_kg * (c.specific_h_j_kg - cells[&id].specific_h_j_kg)
        })
        .sum();
    close(actual_heat, 0.75 * release + 0.9 * plastic);
    let storage = heat.state().normal_response().reversible_energy_j
        + heat.state().rolling_storage_j()
        - old_storage;
    let work: f64 = report.segments.iter().map(|r| r.applied_work_j).sum();
    assert!((storage + actual_heat + report.work.other_work_j() - work).abs() < 1e-20);
    assert!(report.transfer_residual_j.abs() <= controls().transfer_energy_tolerance_j);
    assert_eq!(thermal.time_s(), 0.0);
    let accepted = heat.stamp();
    let accepted_cells = thermal.cells().clone();
    assert_eq!(
        heat.apply(&mut thermal, stamp, &path, controls(), || Ok(())),
        Err(JkrRollingSnapHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    );
    assert_eq!(heat.stamp(), accepted);
    assert_eq!(thermal.cells(), &accepted_cells);
}

#[test]
fn rolling_conversion_requires_a_destination_and_can_leave_all_loss_outside() {
    let normal = source(0.05);
    let events = events(&normal);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let initial = loaded(&events);
    let snap = binding(events, initial, &caloric, &thermal);
    let partition = HeatPartition::new(0.0, 0.0, 1.0).unwrap();
    assert!(JkrRollingDissipationHeat::new(snap.clone(), partition, None).is_err());
    let mut heat = JkrRollingDissipationHeat::new(
        snap,
        partition,
        Some("unresolved deformation processes".into()),
    )
    .unwrap();
    let cells = thermal.cells().clone();
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[Segment::RollingAngleIncrement(
                2.0 * CRITICAL_SHIFT / normal.effective_radius_m(),
            )],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(report.rolling_loss.dissipated_work_j() > 0.0);
    assert_eq!(
        report.work.other_work_j(),
        report.rolling_loss.dissipated_work_j()
    );
    assert_eq!(thermal.cells(), &cells);
}

#[test]
fn transfer_budget_preserves_weak_plastic_loss_beside_a_large_snap() {
    let normal = source(0.05);
    let events = events(&normal);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let initial = events.initial_detached_state(-1e-6).unwrap();
    let attached = events.advance_normal(&initial, 0.0).unwrap();
    let angle = (CRITICAL_SHIFT / normal.effective_radius_m()) * (1.0 + 8.0 * f64::EPSILON);
    let rolling = events.advance_rolling(attached.state(), angle).unwrap();
    let weak_loss = rolling.report().irreversible_rolling_loss_j;
    let snap_energy = attached.report().normal_released_energy_j;
    assert!(weak_loss > 0.0 && weak_loss < 0.125 * f64::EPSILON * snap_energy);
    assert_eq!(snap_energy + weak_loss, snap_energy);
    let outside = HeatPartition::new(0.0, 0.0, 1.0).unwrap();
    let snap = JkrRollingSnapHeat::new(
        events,
        initial,
        &thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        outside,
        Some("resolved elastic waves".into()),
    )
    .unwrap();
    let mut heat = JkrRollingDissipationHeat::new(
        snap,
        outside,
        Some("rolling microstructural rearrangement".into()),
    )
    .unwrap();
    let stamp = heat.stamp();
    let cells = thermal.cells().clone();
    let mut exact = controls();
    exact.transfer_energy_tolerance_j = 0.0;
    assert!(matches!(
        heat.apply(
            &mut thermal,
            stamp.clone(),
            &[Segment::NormalIndentation(0.0), Segment::RollingAngleIncrement(angle)],
            exact,
            || Ok(())
        ),
        Err(JkrRollingSnapHeatError::Heat(JkrSnapHeatError::TransferBalance(defect)))
            if defect < 0.0 && (defect + weak_loss).abs() < 1e-12 * weak_loss
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.rolling_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &cells);
}
