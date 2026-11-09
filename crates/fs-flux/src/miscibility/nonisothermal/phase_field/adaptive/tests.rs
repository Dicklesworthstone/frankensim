//! Time refinement, physical receipts and interface-sensitive demixing control.

use super::*;
use crate::miscibility::adaptive::BinaryMiscibilityTimeControl;
use crate::miscibility::nonisothermal::phase_field::tests::{Law, fixture, link, near};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn time(step: f64) -> BinaryNonisothermalPhaseFieldTimeControl {
    BinaryNonisothermalPhaseFieldTimeControl {
        thermal: BinaryNonisothermalTimeControl {
            composition: BinaryMiscibilityTimeControl {
                initial_step_s: step,
                minimum_step_s: 1e-9,
                maximum_step_s: step,
                absolute_mole_fraction_tolerance: 1e-7,
                relative_mole_fraction_tolerance: 0.0,
                maximum_attempts: 20_000,
            },
            absolute_temperature_tolerance_k: 1e-4,
            relative_temperature_tolerance: 0.0,
        },
        absolute_interfacial_energy_tolerance_j: 1e-6,
        relative_interfacial_energy_tolerance: 0.0,
    }
}

fn spinodal() -> BinaryNonisothermalPhaseField {
    fixture(
        Law {
            a: 0.0,
            b: 1200.0,
            ..Law::default()
        },
        &[[0.001, 0.49, 400.0], [0.001, 0.51, 400.0]],
        &[link(0, 1, 1e-4, 0.0)],
    )
}

#[test]
fn g1_interface_energy_controls_demixing_with_loose_bulk_tolerances() {
    let initial = spinodal();
    let mut loose = time(2.0);
    loose.thermal.composition.absolute_mole_fraction_tolerance = 1.0;
    loose.thermal.absolute_temperature_tolerance_k = 100.0;
    loose.absolute_interfacial_energy_tolerance_j = 100.0;
    let coarse = initial
        .advance_candidate(2.0, loose, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(coarse.advance().accepted_intervals, 1);
    let mut controlled = loose;
    controlled.absolute_interfacial_energy_tolerance_j = 1e-7;
    let fine = initial
        .advance_candidate(2.0, controlled, BinaryNonisothermalSolve::default(), || {
            false
        })
        .unwrap();
    let mut reference = controlled;
    reference.thermal.composition.initial_step_s = 0.001;
    reference.thermal.composition.maximum_step_s = 0.001;
    reference.absolute_interfacial_energy_tolerance_j = 1e-9;
    let reference = initial
        .advance_candidate(2.0, reference, BinaryNonisothermalSolve::default(), || {
            false
        })
        .unwrap();
    assert!(fine.advance().rejected_intervals > 0);
    assert!(fine.advance().accepted_intervals > 1);
    assert!(fine.advance().maximum_accepted_error_ratio <= 1.0);
    let expected = reference.field().mole_fractions()[0][0];
    let coarse_error = (coarse.field().mole_fractions()[0][0] - expected).abs();
    let fine_error = (fine.field().mole_fractions()[0][0] - expected).abs();
    assert!(
        fine_error < 0.1 * coarse_error,
        "{fine_error:e} versus {coarse_error:e}"
    );
}

#[test]
fn g0_equal_total_interface_energy_cannot_hide_spatial_energy_error() {
    let links = [link(0, 1, 0.0, 0.0), link(1, 2, 0.0, 0.0)];
    let a = fixture(
        Law::default(),
        &[
            [0.001, 0.3, 400.0],
            [0.001, 0.4, 400.0],
            [0.001, 0.8, 400.0],
        ],
        &links,
    );
    let b = fixture(
        Law::default(),
        &[
            [0.001, 0.2, 400.0],
            [0.001, 0.6, 400.0],
            [0.001, 0.7, 400.0],
        ],
        &links,
    );
    // Both profiles contain 1.5 mol of each species and 17 J of interface
    // energy, but the edge allocations [1,16] and [16,1] differ by 30 J in L1.
    let mut control = time(1.0);
    control.absolute_interfacial_energy_tolerance_j = 1.0;
    near(
        interface_error_ratio(&a, &b, control, &mut || false).unwrap(),
        30.0,
        2e-13,
    );
}

#[test]
fn g1_temperature_error_refines_uniform_composition_fourier_relaxation() {
    let initial = fixture(
        Law::default(),
        &[[0.001, 0.5, 350.0], [0.001, 0.5, 550.0]],
        &[link(0, 1, 0.0, 2.0)],
    );
    let mut accuracy = time(4.0);
    accuracy
        .thermal
        .composition
        .absolute_mole_fraction_tolerance = 1.0;
    accuracy.thermal.absolute_temperature_tolerance_k = 1e-3;
    accuracy.absolute_interfacial_energy_tolerance_j = 1.0;
    let result = initial
        .advance_candidate(4.0, accuracy, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let exact_left = 450.0 - 100.0 * (-0.1_f64 * 4.0).exp();
    assert!(result.advance().rejected_intervals > 0);
    assert!(result.advance().accepted_intervals > 1);
    assert_eq!(result.field().mole_fractions(), initial.mole_fractions());
    near(result.field().temperatures()[0].value(), exact_left, 0.15);
    near(
        result.field().temperatures()[0].value() + result.field().temperatures()[1].value(),
        900.0,
        1e-8,
    );
}

#[test]
fn g0_only_accepted_demixing_transfers_enter_species_and_local_energy_receipts() {
    let initial = spinodal();
    let result = initial
        .advance_candidate(1.0, time(1.0), BinaryNonisothermalSolve::default(), || {
            false
        })
        .unwrap();
    assert!(result.advance().rejected_intervals > 0);
    assert!(result.advance().completed_step_candidates >= 3 * result.advance().accepted_intervals);
    let amount = result.advance().transferred_a_mol[0];
    near(
        initial.component_moles()[0][0] - result.field().component_moles()[0][0],
        amount,
        2e-14,
    );
    near(
        result.field().component_moles()[1][0] - initial.component_moles()[1][0],
        amount,
        2e-14,
    );
    let cell_energy = |field: &BinaryNonisothermalPhaseField, cell: usize| {
        let n = field.component_moles()[cell];
        let total = n[0] + n[1];
        let dx = field.mole_fractions()[0][0] - field.mole_fractions()[1][0];
        total * 40.0 * (field.temperatures()[cell].value() - 300.0)
            + R * 1200.0 * n[0] * n[1] / total
            + 0.25 * 200.0 * dx * dx
    };
    for cell in 0..2 {
        let change = cell_energy(result.field(), cell) - cell_energy(&initial, cell);
        let expected =
            if cell == 0 { -1.0 } else { 1.0 } * result.advance().transferred_energy_j[0];
        near(change, expected, 3e-8);
    }
    near(
        result.advance().homogeneous_enthalpy_change_j
            + result.advance().interfacial_energy_change_j,
        result.advance().energy_residual_j,
        2e-13,
    );
    assert!(
        result.advance().energy_residual_j.abs()
            <= BinaryNonisothermalSolve::default().energy_tolerance_j
    );
    assert!(result.advance().entropy_change_j_per_k >= 0.0);
}

#[test]
fn g0_work_minimum_step_and_late_cancellation_leave_the_field_unchanged() {
    let initial = spinodal();
    let saved = initial.clone();
    let mut accuracy = time(2.0);
    accuracy.thermal.composition.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_candidate(2.0, accuracy, BinaryNonisothermalSolve::default(), || false),
        Err(BinaryNonisothermalAdvanceError::WorkLimit { .. })
    ));
    accuracy.thermal.composition.minimum_step_s = 2.0;
    assert!(matches!(
        initial.advance_candidate(2.0, accuracy, BinaryNonisothermalSolve::default(), || false),
        Err(BinaryNonisothermalAdvanceError::TemporalAccuracy { .. })
    ));
    let mut polls = 0usize;
    let expected = initial
        .advance_candidate(0.2, time(0.2), BinaryNonisothermalSolve::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0usize;
    assert!(matches!(
        initial.advance_candidate(0.2, time(0.2), BinaryNonisothermalSolve::default(), || {
            seen += 1;
            seen == polls
        }),
        Err(BinaryNonisothermalAdvanceError::Solver(
            BinaryNonisothermalError::Cancelled
        ))
    ));
    assert_eq!(initial, saved);
    assert_eq!(
        initial
            .advance_candidate(0.2, time(0.2), BinaryNonisothermalSolve::default(), || {
                false
            })
            .unwrap(),
        expected
    );
}

#[test]
fn g0_continuation_reaches_the_exact_clock_with_a_short_final_remainder() {
    let initial = spinodal();
    let mut accuracy = time(0.1);
    accuracy.thermal.composition.minimum_step_s = 0.1;
    accuracy
        .thermal
        .composition
        .absolute_mole_fraction_tolerance = 1.0;
    accuracy.thermal.absolute_temperature_tolerance_k = 10.0;
    accuracy.absolute_interfacial_energy_tolerance_j = 10.0;
    let first = initial
        .advance_candidate(0.2, accuracy, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let second = first
        .field()
        .advance_candidate(0.07, accuracy, BinaryNonisothermalSolve::default(), || {
            false
        })
        .unwrap();
    assert_eq!(first.field().elapsed_s(), initial.elapsed_s() + 0.2);
    assert_eq!(second.field().elapsed_s(), first.field().elapsed_s() + 0.07);
    assert_eq!(second.advance().accepted_intervals, 1);
    assert_eq!(second.advance().duration_s, 0.07);
}
