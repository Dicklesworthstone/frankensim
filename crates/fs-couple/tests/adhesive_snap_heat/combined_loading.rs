//! Shared shear/twist failure, explicit release, and actual finite caloric storage.

mod temperature;
mod timed;

use super::{MASS, close, curves, source, thermal};
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::{
    JkrMixedModeEnvelope, SourceJkrMixedMode,
    combined::{JkrCombinedBoundary, loading::*},
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeatError,
    combined_loading::{
        JkrCombinedHeat, JkrCombinedHeatControls, JkrCombinedHeatError,
        JkrCombinedHeatSegment as Segment,
    },
};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::{HeatPartition, WorkLedger};

fn mixed_source() -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(
        source(0.05),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}

fn controls() -> JkrCombinedHeatControls {
    JkrCombinedHeatControls {
        heat: super::controls(),
        total_energy_tolerance_j: 1e-20,
    }
}

// Dimensional directions make the original stationarity equation
// (y²-d)² + lambda² (b y²+t) - 16y = 0. For d=-1,b=4,t=8,
// y=1,lambda=1 is its joint fold; no combined limit solver enters this input.
fn driver<'s>(
    source: &'s SourceJkrMixedMode,
    d: f64,
    b: f64,
    t: f64,
    gate: &CancelGate,
) -> (JkrCombinedLoading<'s>, f64) {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let a = (PI * w * r * r / (8.0 * e)).cbrt();
    let depth = a * a / r;
    let scale = e * depth * depth;
    let shear = (t * scale / (4.0 * source.reduced_shear_modulus_pa())).sqrt();
    let twist = (b * scale / (8.0 * source.torsional_reduced_shear_modulus_pa() * a * a)).sqrt();
    (
        JkrCombinedLoading::new(
            source,
            d * depth,
            [0.6 * shear, 0.8 * shear],
            twist,
            JkrCombinedLoadingControls {
                maximum_iterations: 128,
                maximum_evaluations: 192,
                energy_absolute_tolerance_j: 0.0,
                energy_relative_tolerance: 3e-12,
            },
            gate,
        )
        .unwrap(),
        a,
    )
}

fn binding<'s>(
    loading: JkrCombinedLoading<'s>,
    initial: JkrCombinedLoadingState,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    thermal: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> JkrCombinedHeat<'s> {
    JkrCombinedHeat::new(
        loading,
        initial,
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("resolved elastic waves".into()),
        gate,
    )
    .unwrap()
}

#[test]
fn loading_unloading_and_sign_reversal_return_mechanical_work_without_heating() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.2, &gate).unwrap();
    let caloric = curves(source.normal_source(), 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[
                Segment::LoadFactor(0.7),
                Segment::LoadFactor(-0.4),
                Segment::LoadFactor(0.2),
            ],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert!(report.segments[0].prescribed_work_j > 0.0);
    assert!(report.segments[1].prescribed_work_j < 0.0);
    assert!(report.segments[2].prescribed_work_j < 0.0);
    assert_eq!(heat.state().contact(), initial.contact());
    assert_eq!(heat.state().steps(), initial.steps() + 3);
    assert!(report.prescribed_work_j.abs() <= controls().total_energy_tolerance_j);
    assert!(report.storage_change_j.abs() <= controls().total_energy_tolerance_j);
    assert_eq!(report.released_energy_j, 0.0);
    assert_eq!(report.work, WorkLedger::default());
    assert_eq!(*heat.ledger(), WorkLedger::default());
    assert_eq!(report.thermal.stored_energy_change_j, 0.0);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn shared_radius_fold_deposits_the_complete_signed_potential_once() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, a) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let limit = loading.limit().unwrap();
    assert_eq!(limit.kind, JkrCombinedBoundary::RadialFold);
    close(limit.contact_radius_m, a);
    close(limit.load_factor, 1.0);
    let caloric = curves(source.normal_source(), 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let original_cells = thermal.cells().clone();
    let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
    let stamp = heat.stamp();
    // Reaching the fold is a loaded contact, not permission to deposit heat.
    let at_fold = heat
        .propose(
            &thermal,
            stamp,
            &[Segment::LoadFactor(2.0 * limit.load_factor)],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    let fold = at_fold.history().state().contact().unwrap();
    assert_eq!(fold.boundary, Some(JkrCombinedBoundary::RadialFold));
    assert!(fold.tangential_elastic_energy_j > 0.0);
    assert!(fold.torsional_elastic_energy_j > 0.0);
    assert!(fold.surface_energy_j < 0.0);
    assert_eq!(at_fold.report().work, WorkLedger::default());
    assert_eq!(at_fold.thermal().cells(), &original_cells);
    assert_eq!(heat.stamp(), stamp);

    // Independent elimination of U_a=U_aa=0 from the four-store potential:
    // the original mode-I snap release alone cannot supply this value.
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let depth = a * a / normal.effective_radius_m();
    let expected_release =
        e * a * (6.0 * depth * depth / 5.0 - 2.0 * fold.coordinates[0] * depth / 3.0)
            + (8.0 / 3.0)
                * source.torsional_reduced_shear_modulus_pa()
                * a.powi(3)
                * fold.coordinates[3].powi(2);
    let report = heat
        .apply(
            &mut thermal,
            stamp,
            &[
                Segment::LoadFactor(2.0 * limit.load_factor),
                Segment::DetachAtFold,
            ],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    close(report.released_energy_j, expected_release);
    close(report.work.dissipated_work_j(), expected_release);
    close(report.work.surface_a_heat_j(), 0.25 * expected_release);
    close(report.work.surface_b_heat_j(), 0.5 * expected_release);
    close(report.work.other_work_j(), 0.25 * expected_release);
    close(
        report.thermal.stored_energy_change_j,
        0.75 * expected_release,
    );
    assert_eq!(report.segments[0].accepted_load_factor, limit.load_factor);
    assert_eq!(
        report.segments[0].untraversed_load_factor,
        limit.load_factor
    );
    assert_eq!(report.segments[0].released_energy_j, 0.0);
    assert_eq!(report.segments[1].prescribed_work_j, 0.0);
    assert_eq!(heat.state().contact(), None);
    assert_eq!(report.final_state, heat.state());
    close(heat.state().released_energy_j(), expected_release);
    close(report.storage_change_j, -initial.stored_energy_j());
    assert!(
        (-initial.stored_energy_j()
            + report.thermal.stored_energy_change_j
            + report.work.other_work_j()
            - report.prescribed_work_j)
            .abs()
            <= controls().total_energy_tolerance_j
    );
    assert!(report.transfer_residual_j.abs() <= controls().heat.transfer_energy_tolerance_j);
    assert!(report.absolute_transfer_residual_j <= controls().heat.transfer_energy_tolerance_j);
    assert!(report.combined_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
    assert!(report.absolute_stage_energy_residual_j <= controls().total_energy_tolerance_j);
    assert_eq!(
        report.outside_destination.as_deref(),
        Some("resolved elastic waves")
    );
    for (id, share) in [(CellId(0), 0.25), (CellId(1), 0.5)] {
        close(
            MASS * thermal.cells()[&id].specific_h_j_kg,
            share * expected_release,
        );
        assert_eq!(thermal.cells()[&id].geometry, original_cells[&id].geometry);
        assert_eq!(thermal.cells()[&id].mass_kg, original_cells[&id].mass_kg);
    }
    let accepted_cells = thermal.cells().clone();
    let accepted_stamp = heat.stamp();
    assert!(matches!(
        heat.apply(
            &mut thermal,
            stamp,
            &[Segment::DetachAtFold],
            controls(),
            &gate,
            || Ok(())
        ),
        Err(JkrCombinedHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    ));
    assert!(matches!(
        heat.apply(
            &mut thermal,
            accepted_stamp,
            &[Segment::DetachAtFold],
            controls(),
            &gate,
            || Ok(())
        ),
        Err(JkrCombinedHeatError::Source(
            JkrCombinedLoadingError::NotAtSeparatingFold
        ))
    ));
    assert_eq!(heat.stamp(), accepted_stamp);
    assert_eq!(thermal.cells(), &accepted_cells);
    let free = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[Segment::LoadFactor(-2.0)],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(free.work, WorkLedger::default());
    assert_eq!(free.prescribed_work_j, 0.0);
    assert_eq!(heat.state().contact(), None);
    assert_eq!(thermal.cells(), &accepted_cells);
    close(heat.ledger().dissipated_work_j(), expected_release);
}

#[test]
fn imported_separated_history_does_not_deposit_prior_release_or_rebond() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let fold = loading
        .advance_to_load_factor(initial, 2.0, &gate)
        .unwrap()
        .state();
    let separated = loading.detach_at_fold(fold, &gate).unwrap().state();
    assert!(separated.released_energy_j() > 0.0);
    let caloric = curves(source.normal_source(), 1.0);
    let mut thermal = thermal(&caloric, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = binding(loading, separated, &caloric, &thermal, &gate);
    assert_eq!(heat.state(), separated);
    assert_eq!(*heat.ledger(), WorkLedger::default());
    let report = heat
        .apply(
            &mut thermal,
            heat.stamp(),
            &[Segment::LoadFactor(0.0)],
            controls(),
            &gate,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.released_energy_j, 0.0);
    assert_eq!(report.work, WorkLedger::default());
    assert_eq!(heat.state().contact(), None);
    assert_eq!(
        heat.state().released_energy_j(),
        separated.released_energy_j()
    );
    assert_eq!(thermal.cells(), &cells);
}

#[test]
fn nonseparating_boundaries_refuse_the_entire_loading_and_heat_transaction() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let caloric = curves(source.normal_source(), 1.0);
    for (b, t, kind) in [
        (8.0, 2.0, JkrCombinedBoundary::OpeningEdge),
        (2.0, 8.0, JkrCombinedBoundary::RadialFold),
    ] {
        let (loading, _) = driver(&source, 1.0, b, t, &gate);
        let limit = loading.limit().unwrap();
        assert_eq!(limit.kind, kind);
        let initial = loading.initial_state(0.0, &gate).unwrap();
        let mut thermal = thermal(&caloric, 1.0);
        let cells = thermal.cells().clone();
        let mut heat = binding(loading, initial, &caloric, &thermal, &gate);
        let stamp = heat.stamp();
        assert!(matches!(
            heat.apply(
                &mut thermal,
                stamp,
                &[
                    Segment::LoadFactor(2.0 * limit.load_factor),
                    Segment::DetachAtFold
                ],
                controls(),
                &gate,
                || Ok(())
            ),
            Err(JkrCombinedHeatError::Source(
                JkrCombinedLoadingError::NotAtSeparatingFold
            ))
        ));
        assert_eq!(heat.stamp(), stamp);
        assert_eq!(heat.state(), initial);
        assert_eq!(*heat.ledger(), WorkLedger::default());
        assert_eq!(thermal.cells(), &cells);
    }
}

#[test]
fn source_caloric_identity_and_finite_receiver_ownership_are_not_substitutable() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let (foreign, _) = driver(&source, -1.0, 4.0, 9.0, &gate);
    let foreign_state = foreign.initial_state(0.0, &gate).unwrap();
    // Even an identical zero-load physical point is not the same loading family.
    assert_eq!(foreign_state.contact(), initial.contact());
    let caloric = curves(source.normal_source(), 1.0);
    let thermal = thermal(&caloric, 1.0);
    assert!(matches!(
        JkrCombinedHeat::new(
            loading,
            foreign_state,
            &thermal,
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            None,
            &gate
        ),
        Err(JkrCombinedHeatError::Source(
            JkrCombinedLoadingError::StateMismatch
        ))
    ));
    let copy = caloric.clone();
    for selected in [[&copy[0], &copy[1]], [&caloric[1], &caloric[0]]] {
        assert!(
            JkrCombinedHeat::new(
                loading,
                initial,
                &thermal,
                [CellId(0), CellId(1)],
                selected,
                HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
                None,
                &gate
            )
            .is_err()
        );
    }
    assert!(
        JkrCombinedHeat::new(
            loading,
            initial,
            &thermal,
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
            None,
            &gate
        )
        .is_err()
    );
    let heat = binding(loading, initial, &caloric, &thermal, &gate);
    let changed_mass = super::thermal(&caloric, 2.0);
    assert!(matches!(
        heat.propose(
            &changed_mass,
            heat.stamp(),
            &[Segment::LoadFactor(0.2)],
            controls(),
            &gate,
            || Ok(())
        ),
        Err(JkrCombinedHeatError::Heat(JkrSnapHeatError::Invalid(_)))
    ));
    assert_eq!(heat.state(), initial);
}

#[test]
fn late_cancellation_and_second_cell_caloric_failure_preserve_both_owners() {
    let source = mixed_source();
    let gate = CancelGate::new_clock_free();
    let (loading, _) = driver(&source, -1.0, 4.0, 8.0, &gate);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let caloric = curves(source.normal_source(), 1.0);
    let original = thermal(&caloric, 1.0);
    let segments = [Segment::LoadFactor(2.0), Segment::DetachAtFold];
    let mut reference_thermal = original.clone();
    let mut reference = binding(loading, initial, &caloric, &original, &gate);
    let mut calls = 0;
    let expected = reference
        .apply(
            &mut reference_thermal,
            reference.stamp(),
            &segments,
            controls(),
            &gate,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    assert!(calls > 1);
    let mut heat = binding(loading, initial, &caloric, &original, &gate);
    let mut thermal = original.clone();
    let stamp = heat.stamp();
    let mut polls = 0;
    assert!(matches!(
        heat.apply(&mut thermal, stamp, &segments, controls(), &gate, || {
            polls += 1;
            if polls == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(JkrCombinedHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(polls, calls);
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(*heat.ledger(), WorkLedger::default());
    assert_eq!(thermal.cells(), original.cells());
    let retry = heat
        .apply(&mut thermal, stamp, &segments, controls(), &gate, || Ok(()))
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(thermal.cells(), reference_thermal.cells());

    // The first receiver has enough caloric range; the second does not. No
    // privately staged first-cell heating or release history may escape.
    let mut narrow = curves(source.normal_source(), 1e-6);
    narrow[0] = caloric[0].clone();
    let mut thermal = super::thermal(&narrow, 1.0);
    let cells = thermal.cells().clone();
    let mut heat = binding(loading, initial, &narrow, &thermal, &gate);
    let stamp = heat.stamp();
    assert!(matches!(
        heat.apply(&mut thermal, stamp, &segments, controls(), &gate, || Ok(())),
        Err(JkrCombinedHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(heat.stamp(), stamp);
    assert_eq!(*heat.ledger(), WorkLedger::default());
    assert_eq!(thermal.cells(), &cells);
}
