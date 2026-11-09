//! Bounded phase-field refinement and whole-interval accounting oracles.

use super::*;
use super::super::tests::{
    assert_inventories, bulk_energy, field, interface_energy, near, ring,
};

fn accuracy(step: f64) -> MulticomponentPhaseFieldTimeControl {
    MulticomponentPhaseFieldTimeControl {
        composition: MulticomponentMiscibilityTimeControl {
            initial_step_s: step,
            minimum_step_s: 1e-9,
            maximum_step_s: step,
            absolute_mole_fraction_tolerance: 1e-9,
            relative_mole_fraction_tolerance: 0.0,
            maximum_attempts: 2_000,
        },
        absolute_interfacial_energy_tolerance_j: 1e-3,
        relative_interfacial_energy_tolerance: 0.0,
    }
}

fn edge_energies(field: &MulticomponentPhaseField, kappa: f64) -> Vec<f64> {
    field
        .network()
        .faces()
        .iter()
        .map(|face| {
            let left = field.cell_mole_fractions(face.cells[0]).unwrap();
            let right = field.cell_mole_fractions(face.cells[1]).unwrap();
            kappa * face.area_m2 / (4.0 * face.distance_m)
                * left
                    .iter()
                    .zip(right)
                    .map(|(a, b)| (b - a) * (b - a))
                    .sum::<f64>()
        })
        .collect()
}

#[test]
fn large_initial_interval_refines_and_matches_a_fine_full_component_reference() {
    let (initial, chi, kappa) = ring(5, 0.01);
    let original = initial.clone();
    let controls = MulticomponentMiscibilitySolve::default();
    let candidate = initial
        .advance_candidate(1.0, accuracy(1.0), controls, || false)
        .unwrap();
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    assert_eq!(candidate.field().elapsed_s(), 1.0);
    assert_eq!(initial, original);

    // Binary-representable subintervals reach the exact endpoint. This is a
    // same-law fine reference, independent of the adaptive error estimator.
    let mut reference = initial.clone();
    for _ in 0..2_048 {
        reference = reference
            .step_candidate(1.0 / 2_048.0, controls, || false)
            .unwrap()
            .into_parts()
            .0;
    }
    assert_eq!(reference.elapsed_s(), 1.0);
    for (&actual, &expected) in candidate
        .field()
        .mole_fractions()
        .iter()
        .zip(reference.mole_fractions())
    {
        near(actual, expected, 2e-7);
    }
    let actual_bulk = bulk_energy(candidate.field(), &chi) - bulk_energy(&initial, &chi);
    let actual_interface =
        interface_energy(candidate.field(), kappa) - interface_energy(&initial, kappa);
    let reference_energy = bulk_energy(&reference, &chi) - bulk_energy(&initial, &chi)
        + interface_energy(&reference, kappa) - interface_energy(&initial, kappa);
    near(candidate.advance().bulk_free_energy_change_j, actual_bulk, 3e-7);
    near(candidate.advance().interfacial_energy_change_j, actual_interface, 1e-8);
    near(
        candidate.advance().free_energy_change_j,
        actual_bulk + actual_interface,
        3e-7,
    );
    near(candidate.advance().free_energy_change_j, reference_energy, 0.05);
    assert!(candidate.advance().free_energy_change_j < 0.0);
    assert_inventories(&initial, candidate.field());

    let report = candidate.advance();
    let faces = initial.network().faces().len();
    assert_eq!(report.attempted_intervals, report.accepted_intervals + report.rejected_intervals);
    assert!(report.completed_steps >= 3 * report.accepted_intervals);
    assert!(report.completed_steps <= 3 * report.attempted_intervals);
    assert!(
        report.completed_maxwell_stefan_evaluations
            <= faces * (report.accepted_intervals + report.attempted_intervals)
    );
    assert!(report.completed_maxwell_stefan_evaluations < 2 * faces * report.attempted_intervals);
}

#[test]
fn interface_energy_l1_prevents_opposite_face_errors_from_hiding_behind_total_energy() {
    let (base, _, kappa) = ring(2, 0.01);
    let cells = base.network().cells().len();
    let initial = field(
        base.network(),
        (0..cells)
            .map(|cell| {
                let mut x = base.cell_mole_fractions(cell).unwrap().to_vec();
                let short_wave =
                    0.0025 * (2.0 * std::f64::consts::PI * 5.0 * cell as f64 / cells as f64).cos();
                x[0] += short_wave;
                x[1] -= short_wave;
                x
            })
            .collect(),
    );
    let snapshot = initial.clone();
    let controls = MulticomponentMiscibilitySolve::default();
    let coarse = initial.step_candidate(1.0, controls, || false).unwrap();
    let first = initial.step_candidate(0.5, controls, || false).unwrap();
    let fine = first.field().step_candidate(0.5, controls, || false).unwrap();
    let errors = edge_energies(fine.field(), kappa)
        .into_iter()
        .zip(edge_energies(coarse.field(), kappa))
        .map(|(fine, coarse)| fine - coarse)
        .collect::<Vec<_>>();
    assert!(errors.iter().any(|error| *error > 0.0));
    assert!(errors.iter().any(|error| *error < 0.0));
    let signed_error = errors.iter().sum::<f64>().abs();
    let l1_error = errors.iter().map(|error| error.abs()).sum::<f64>();
    assert!(l1_error > 1.1 * signed_error);

    // A total-energy-only estimator would accept this trial. The authored
    // absolute allowance lies strictly between |sum(error)| and sum(|error|).
    let mut time = accuracy(1.0);
    time.composition.absolute_mole_fraction_tolerance = 1.0;
    time.composition.maximum_attempts = 1;
    time.absolute_interfacial_energy_tolerance_j = 0.5 * (signed_error + l1_error);
    assert!(matches!(
        initial.advance_candidate(1.0, time, controls, || false),
        Err(MulticomponentPhaseFieldAdvanceError::WorkLimit {
            attempted_intervals: 1,
            accepted_intervals: 0,
            elapsed_s,
        }) if elapsed_s == 0.0
    ));
    assert_eq!(initial, snapshot);
    time.composition.maximum_attempts = 64;
    let candidate = initial.advance_candidate(1.0, time, controls, || false).unwrap();
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    assert_eq!(candidate.field().elapsed_s(), 1.0);
    assert_inventories(&initial, candidate.field());
}

#[test]
fn complete_window_uses_one_energy_allowance_and_late_cancellation_is_atomic() {
    let (initial, chi, kappa) = ring(8, 1e-4);
    let snapshot = initial.clone();
    // Deliberately large authored energy allowances make individual unstable
    // explicit steps inspectable. This tests accounting, not physical accuracy.
    // The alternating mode has lambda=-4/15; lambda*duration=-7 makes the
    // two-half path grow more in energy than any one of its three trial steps.
    let duration = 26.25;
    let loose_controls = MulticomponentMiscibilitySolve {
        free_energy_tolerance_j: 1e6,
        ..MulticomponentMiscibilitySolve::default()
    };
    let coarse = initial.step_candidate(duration, loose_controls, || false).unwrap();
    let first = initial.step_candidate(0.5 * duration, loose_controls, || false).unwrap();
    let fine = first
        .field()
        .step_candidate(0.5 * duration, loose_controls, || false)
        .unwrap();
    let energy = |state: &MulticomponentPhaseField| {
        bulk_energy(state, &chi) + interface_energy(state, kappa)
    };
    let initial_energy = energy(&initial);
    let first_energy = energy(first.field());
    let final_energy = energy(fine.field());
    let maximum_trial_increase = (energy(coarse.field()) - initial_energy)
        .max(first_energy - initial_energy)
        .max(final_energy - first_energy);
    let allowance = 1.01 * maximum_trial_increase;
    assert!(allowance > 0.0);
    assert!(final_energy - initial_energy > 1.01 * allowance);
    let mut time = accuracy(duration);
    time.composition.minimum_step_s = duration;
    time.composition.absolute_mole_fraction_tolerance = 1.0;
    time.composition.maximum_attempts = 4;
    time.absolute_interfacial_energy_tolerance_j = 1e6;

    let controls = MulticomponentMiscibilitySolve {
        free_energy_tolerance_j: allowance,
        ..loose_controls
    };
    // All full/half candidates fit individually, but the complete physical
    // window exceeds the one declared allowance and must remain unpublished.
    initial.step_candidate(duration, controls, || false).unwrap();
    initial.step_candidate(0.5 * duration, controls, || false).unwrap();
    first.field().step_candidate(0.5 * duration, controls, || false).unwrap();
    assert!(matches!(
        initial.advance_candidate(duration, time, controls, || false),
        Err(MulticomponentPhaseFieldAdvanceError::Solver(
            MulticomponentMiscibilityError::FreeEnergy { change_j, tolerance_j }
        )) if change_j > tolerance_j && tolerance_j == allowance
    ));
    assert_eq!(initial, snapshot);

    let mut completed_polls = 0;
    let accepted = initial
        .advance_candidate(duration, time, loose_controls, || {
            completed_polls += 1;
            false
        })
        .unwrap();
    assert_eq!(accepted.field(), fine.field());
    assert_eq!(accepted.advance().attempted_intervals, 1);
    assert_eq!(accepted.advance().accepted_intervals, 1);
    assert_eq!(accepted.advance().rejected_intervals, 0);
    assert_eq!(accepted.advance().completed_steps, 3);
    assert_eq!(
        accepted.advance().completed_maxwell_stefan_evaluations,
        2 * initial.network().faces().len(),
    );
    near(
        accepted.advance().free_energy_change_j,
        final_energy - initial_energy,
        3e-7,
    );
    for face in 0..initial.network().faces().len() {
        for species in 0..initial.network().species_count() {
            near(
                accepted.advance().transferred_mol[face][species],
                first.step().transferred_mol[face][species]
                    + fine.step().transferred_mol[face][species],
                1e-12,
            );
        }
    }
    let mut polls = 0;
    assert!(matches!(
        initial.advance_candidate(2.0 * duration, time, loose_controls, || {
            polls += 1;
            polls > completed_polls + 8
        }),
        Err(MulticomponentPhaseFieldAdvanceError::Solver(
            MulticomponentMiscibilityError::Cancelled
        ))
    ));
    assert!(polls > completed_polls);
    assert_eq!(initial, snapshot);
}
