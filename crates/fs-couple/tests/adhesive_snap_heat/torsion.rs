//! G1/G3/G4: retained twist release, plastic loss and original finite calorics.

#[path = "torsion/native.rs"]
mod native;
#[path = "torsion/timed.rs"]
mod timed;

use super::{MASS, close, controls, curves, source, thermal};
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrJumpKind, SourceJkrSphereContact,
    torsion::{
        JkrTorsionControls,
        events::{JkrTorsionEventState, SourceJkrTorsionEvents},
    },
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    torsion::{JkrTorsionHeat, JkrTorsionHeatError, JkrTorsionHeatSegment as Segment},
};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::{HeatPartition, WorkLedger};

const CRITICAL_TWIST: f64 = 1e-5;

fn events(normal: &SourceJkrSphereContact) -> SourceJkrTorsionEvents {
    SourceJkrTorsionEvents::new(
        normal.clone(),
        JkrTorsionControls {
            critical_twist_angle_rad: CRITICAL_TWIST,
            maximum_critical_twist_angle_rad: 0.01,
            maximum_mode_iii_energy_release_to_adhesion: 0.01,
        },
    )
    .unwrap()
}

fn loaded(events: &SourceJkrTorsionEvents) -> JkrTorsionEventState {
    let delta = events
        .normal_source()
        .evaluate_force_controlled(0.0)
        .unwrap()
        .indentation_m;
    events
        .advance_torsion(
            &events.initial_attached_state(delta).unwrap(),
            2.5 * CRITICAL_TWIST,
        )
        .unwrap()
        .into_parts()
        .0
}

fn binding<'s>(
    events: SourceJkrTorsionEvents,
    state: JkrTorsionEventState,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
) -> JkrTorsionHeat<'s> {
    JkrTorsionHeat::new(
        events,
        state,
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("resolved elastic waves".into()),
        HeatPartition::new(0.5, 0.25, 0.25).unwrap(),
        Some("torsional microstructural rearrangement".into()),
    )
    .unwrap()
}

fn ledgers(heat: &JkrTorsionHeat<'_>) -> [WorkLedger; 3] {
    [
        *heat.normal_release_ledger(),
        *heat.torsion_release_ledger(),
        *heat.torsion_loss_ledger(),
    ]
}

#[test]
fn original_fold_twist_energy_heats_once_and_recontact_starts_unstrained() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let prior_loss = initial.accumulated_torsion_loss_j();
    assert!(prior_loss > 0.0);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let original_cells = thermal.cells().clone();
    let mut heat = binding(events.clone(), initial.clone(), &caloric, &thermal);
    assert_eq!(ledgers(&heat), [WorkLedger::default(); 3]);
    assert_eq!(thermal.cells(), &original_cells);
    let elastic = normal.normal_state().elastic();
    let g = [elastic.surface_a_state(), elastic.surface_b_state()]
        .map(|s| s.elastic_state().shear_modulus_pa());
    let reduced_g = (g[0].recip() + g[1].recip()).recip();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let e = elastic.reduced_modulus_pa();
    let phi = initial.attached().unwrap().1.elastic_twist_rad();
    let a_s = (PI * w * r * r / (8.0 * e)).cbrt();
    // Two circular-punch compliances and the original mode-I fold imply
    // U_twist,fold = pi/3 (G_bar/E*) W R^2 phi_e^2 = U_zero_force/36.
    let twist_release = (PI / 3.0) * (reduced_g / e) * w * r * r * phi * phi;
    let normal_release = 0.4 * PI * w * a_s * a_s;
    let fold = normal.minimum_attached_indentation_m();
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[Segment::NormalIndentation(fold)],
            controls(),
            || Ok(()),
        )
        .unwrap();
    close(report.torsion_release.dissipated_work_j(), twist_release);
    close(twist_release, initial.torsion_storage_j() / 36.0);
    close(report.normal_release.dissipated_work_j(), normal_release);
    close(
        report.segments[0].torsion_parameter_work_j,
        twist_release - initial.torsion_storage_j(),
    );
    assert_eq!(report.torsion_loss.dissipated_work_j(), 0.0);
    close(
        report.work.dissipated_work_j(),
        normal_release + twist_release,
    );
    close(
        report.thermal.stored_energy_change_j,
        0.75 * (normal_release + twist_release),
    );
    for (id, fraction) in [(CellId(0), 0.25), (CellId(1), 0.5)] {
        close(
            MASS * thermal.cells()[&id].specific_h_j_kg,
            fraction * (normal_release + twist_release),
        );
    }
    let storage = -initial.normal_response().reversible_energy_j - initial.torsion_storage_j();
    assert!(
        (storage + report.thermal.stored_energy_change_j + report.work.other_work_j()
            - report.segments[0].applied_work_j)
            .abs()
            <= controls().transfer_energy_tolerance_j
    );
    assert_eq!(heat.state().normal_response().branch, JkrBranch::Detached);
    assert_eq!(heat.state().accumulated_torsion_loss_j(), prior_loss);
    let recontact = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[
                Segment::TorsionAngleIncrement(0.4),
                Segment::NormalIndentation(0.0),
            ],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(recontact.segments[0].applied_work_j, 0.0);
    assert_eq!(recontact.segments[0].irreversible_torsion_loss_j, 0.0);
    assert_eq!(
        recontact.segments[1]
            .normal_transition
            .unwrap()
            .jump
            .unwrap()
            .kind,
        JkrJumpKind::Attach
    );
    assert_eq!(heat.state().attached().unwrap().1.elastic_twist_rad(), 0.0);
    assert_eq!(heat.state().torsion_storage_j(), 0.0);
    close(
        heat.torsion_release_ledger().dissipated_work_j(),
        twist_release,
    );
    assert_eq!(heat.torsion_loss_ledger().dissipated_work_j(), 0.0);
    assert_eq!(heat.state().accumulated_torsion_loss_j(), prior_loss);
    assert_eq!(thermal.time_s(), 0.0);
    // Binding a marginal attached source is admission, not another snap event.
    let marginal = events.initial_attached_state(fold).unwrap();
    assert_eq!(
        binding(events, marginal.clone(), &caloric, &thermal).state(),
        &marginal
    );
}

#[test]
fn elastic_and_area_work_stay_mechanical_while_original_yield_heats_its_declared_shares() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = events.initial_attached_state(0.0).unwrap();
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = binding(events, initial, &caloric, &thermal);
    let elastic = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[
                Segment::TorsionAngleIncrement(0.25 * CRITICAL_TWIST),
                Segment::NormalIndentation(1e-8),
                Segment::TorsionAngleIncrement(-0.25 * CRITICAL_TWIST),
            ],
            controls(),
            || Ok(()),
        )
        .unwrap();
    assert!(elastic.segments[0].storage_change_j > 0.0);
    assert!(elastic.segments[1].torsion_parameter_work_j > 0.0);
    assert!(elastic.segments[2].storage_change_j < 0.0);
    assert_eq!(elastic.work.dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &cells);
    let moment = heat
        .state()
        .attached()
        .unwrap()
        .0
        .law()
        .critical_moment_n_m();
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[Segment::TorsionAngleIncrement(3.0 * CRITICAL_TWIST)],
            controls(),
            || Ok(()),
        )
        .unwrap();
    let loss = 2.0 * CRITICAL_TWIST * moment;
    close(report.torsion_loss.dissipated_work_j(), loss);
    close(report.torsion_loss.surface_a_heat_j(), 0.5 * loss);
    close(report.torsion_loss.surface_b_heat_j(), 0.25 * loss);
    close(report.torsion_loss.other_work_j(), 0.25 * loss);
    assert_eq!(report.normal_release.dissipated_work_j(), 0.0);
    assert_eq!(report.torsion_release.dissipated_work_j(), 0.0);
    close(report.thermal.stored_energy_change_j, 0.75 * loss);
    close(
        report.segments[0].storage_change_j
            + report.thermal.stored_energy_change_j
            + report.work.other_work_j(),
        report.segments[0].applied_work_j,
    );
    assert_eq!(
        report.release_outside_destination.as_deref(),
        Some("resolved elastic waves")
    );
    assert_eq!(
        report.torsion_outside_destination.as_deref(),
        Some("torsional microstructural rearrangement")
    );
    assert!(report.absolute_transfer_residual_j <= controls().transfer_energy_tolerance_j);
}

#[test]
fn finite_caloric_refusal_preserves_the_loaded_source_and_all_loss_ledgers() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let loss = events
        .advance_torsion(&initial, 2.0 * CRITICAL_TWIST)
        .unwrap()
        .report()
        .irreversible_torsion_loss_j;
    assert!(loss > 0.0);
    let release = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let caloric = curves(&normal, loss / (100.0 * release));
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = binding(events, initial.clone(), &caloric, &thermal);
    let stamp = heat.stamp();
    assert!(matches!(
        heat.apply(
            &mut thermal,
            stamp.clone(),
            &[Segment::TorsionAngleIncrement(2.0 * CRITICAL_TWIST)],
            controls(),
            || Ok(()),
        ),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(heat.state(), &initial);
    assert_eq!(ledgers(&heat), [WorkLedger::default(); 3]);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn budget_and_final_callback_refusal_leave_a_deterministic_retry_without_redeposit() {
    let normal = source(0.05);
    let events = events(&normal);
    let initial = loaded(&events);
    let caloric = curves(&normal, 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let original = thermal.clone();
    let mut heat = binding(events, initial, &caloric, &thermal);
    let stamp = heat.stamp();
    let path = [
        Segment::TorsionAngleIncrement(CRITICAL_TWIST),
        Segment::NormalIndentation(normal.minimum_attached_indentation_m()),
    ];
    let mut bounded = controls();
    bounded.maximum_segments = 1;
    assert!(
        heat.apply(&mut thermal, stamp.clone(), &path, bounded, || Ok(()))
            .is_err()
    );
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(thermal.cells(), original.cells());
    let mut proposal_calls = 0;
    let proposal = heat
        .propose(&thermal, stamp.clone(), &path, controls(), || {
            proposal_calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(proposal.report().torsion_loss.dissipated_work_j() > 0.0);
    assert!(proposal.report().torsion_release.dissipated_work_j() > 0.0);
    let mut calls = 0;
    assert!(matches!(
        heat.apply(&mut thermal, stamp.clone(), &path, controls(), || {
            calls += 1;
            if calls == proposal_calls + 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(calls, proposal_calls + 1);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(ledgers(&heat), [WorkLedger::default(); 3]);
    assert_eq!(thermal.cells(), original.cells());
    assert_eq!(thermal.time_s(), original.time_s());
    let report = heat
        .apply(&mut thermal, stamp.clone(), &path, controls(), || Ok(()))
        .unwrap();
    assert_eq!(&report, proposal.report());
    assert_eq!(heat.stamp(), proposal.history().stamp());
    assert_eq!(ledgers(&heat), ledgers(proposal.history()));
    assert_eq!(thermal.cells(), proposal.thermal().cells());
    let accepted = heat.stamp();
    assert!(matches!(
        heat.apply(&mut thermal, stamp, &path, controls(), || Ok(())),
        Err(JkrTorsionHeatError::Heat(JkrSnapHeatError::HistoryMismatch))
    ));
    assert_eq!(heat.stamp(), accepted);
    assert_eq!(thermal.cells(), proposal.thermal().cells());
}
