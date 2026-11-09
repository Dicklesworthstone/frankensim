//! Independent thermal/composition modes and accepted-path conservation.

use super::*;
use crate::miscibility::adaptive::BinaryMiscibilityTimeControl;
use crate::miscibility::nonisothermal::advection::tests::{enthalpy, fixture, near, opposing_flow};

fn time(composition_tolerance: f64, temperature_tolerance: f64) -> BinaryNonisothermalTimeControl {
    BinaryNonisothermalTimeControl {
        composition: BinaryMiscibilityTimeControl {
            initial_step_s: 1.0,
            minimum_step_s: 1e-7,
            maximum_step_s: 1.0,
            absolute_mole_fraction_tolerance: composition_tolerance,
            relative_mole_fraction_tolerance: 0.0,
            maximum_attempts: 20_000,
        },
        absolute_temperature_tolerance_k: temperature_tolerance,
        relative_temperature_tolerance: 0.0,
    }
}

fn advance(
    initial: &BinaryNonisothermalField,
    flow: &BinaryMiscibilityAdvection,
    duration: f64,
    time: BinaryNonisothermalTimeControl,
) -> BinaryNonisothermalTransportAdvanceCandidate {
    initial
        .advance_transport_candidate(
            flow,
            duration,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            time,
            || false,
        )
        .unwrap()
}

#[test]
fn g1_temperature_error_drives_refinement_even_with_uniform_composition() {
    let initial = fixture(
        [20.0, 60.0],
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.4, 340.0], [0.001, 0.4, 560.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0005);
    // Composition is uniform, so the composition-only estimate is zero.
    // With equal capacities in the two cells, T_left=450-110*exp(-t).
    let exact = 450.0 - 110.0 * (-1.0_f64).exp();
    let coarse = advance(&initial, &flow, 1.0, time(1.0, 0.1));
    let fine = advance(&initial, &flow, 1.0, time(1.0, 0.001));
    let coarse_error = (coarse.field().temperatures()[0].value() - exact).abs();
    let fine_error = (fine.field().temperatures()[0].value() - exact).abs();
    assert!(
        fine_error < coarse_error / 5.0,
        "{fine_error:e} vs {coarse_error:e}"
    );
    assert!(fine_error < 0.1, "thermal error {fine_error:e} K");
    assert!(fine.advance().accepted_intervals > coarse.advance().accepted_intervals);
    assert!(fine.advance().rejected_intervals > 0);
    assert!(fine.advance().maximum_accepted_error_ratio <= 1.0);
    for fractions in fine.field().mole_fractions() {
        near(fractions[0], 0.4, 2e-13);
    }
    near(
        fine.field().temperatures()[0].value() + fine.field().temperatures()[1].value(),
        900.0,
        1e-8,
    );
    assert_eq!(fine.field().elapsed_s(), 1.0);
}

#[test]
fn g1_single_donor_trials_prevent_cfl_subcycling_from_hiding_mixing_error() {
    let initial = fixture(
        [40.0; 2],
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.9, 400.0], [0.001, 0.1, 400.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.001);
    let step = |field: &BinaryNonisothermalField, duration| {
        field
            .transport_candidate(
                &flow,
                duration,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryNonisothermalSolve::default(),
                || false,
            )
            .unwrap()
    };
    // Both naive paths use the identical two CFL=.5 donor steps. They report
    // zero error after artificial instantaneous homogenization, not exp decay.
    let naive_coarse = step(&initial, 1.0);
    let naive_half = step(&initial, 0.5);
    let naive_fine = step(naive_half.field(), 0.5);
    assert_eq!(naive_coarse.field(), naive_fine.field());
    near(
        naive_coarse.field().mole_fractions()[0][0],
        0.5,
        2e-15,
    );
    let result = advance(&initial, &flow, 1.0, time(1e-6, 1.0));
    let exact = 0.5 + 0.4 * (-2.0_f64).exp();
    near(result.field().mole_fractions()[0][0], exact, 2e-4);
    assert!(result.advance().accepted_intervals > 1);
    assert!(result.advance().maximum_accepted_step_s <= 0.5);
    assert!(result.advance().maximum_advection_courant <= 0.5);
    assert!(result.advance().maximum_accepted_error_ratio <= 1.0);
    near(result.field().temperatures()[0].value(), 400.0, 1e-8);
}

#[test]
fn g0_only_accepted_half_steps_contribute_to_species_and_energy_receipts() {
    let cp = [20.0, 60.0];
    let b = -100.0;
    let initial = fixture(
        cp,
        b,
        1e-5,
        0.03,
        &[[0.001, 0.3, 340.0], [0.002, 0.7, 520.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0003);
    let result = advance(&initial, &flow, 0.3, time(1e-5, 0.001));
    let report = result.advance();
    assert!(report.rejected_intervals > 0);
    for (cell, sign) in [(0, -1.0), (1, 1.0)] {
        for species in 0..2 {
            near(
                result.field().component_moles()[cell][species]
                    - initial.component_moles()[cell][species],
                sign * report.transferred_component_moles[0][species],
                3e-13,
            );
        }
        let actual_heat = enthalpy(
            result.field().component_moles()[cell],
            result.field().temperatures()[cell].value(),
            cp,
            b,
        ) - enthalpy(
            initial.component_moles()[cell],
            initial.temperatures()[cell].value(),
            cp,
            b,
        );
        near(actual_heat, sign * report.transferred_energy_j[0], 1e-8);
    }
    assert!(
        report
            .component_residual_mol
            .iter()
            .all(|r| r.abs() <= 1e-12)
    );
    assert!(report.maximum_capacity_residual_mol <= 1e-12);
    near(report.energy_residual_j, 0.0, 1e-8);
    assert!(report.entropy_change_j_per_k > 0.0);
    assert!(report.completed_transport_steps >= 3 * report.accepted_intervals);
    assert_eq!(result.field().elapsed_s(), 0.3);
}

#[test]
fn g0_exhausted_work_and_an_impossible_minimum_cfl_step_are_atomic() {
    let initial = fixture(
        [40.0; 2],
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.2, 400.0], [0.001, 0.8, 400.0]],
        &[[0, 1]],
    );
    let original = initial.clone();
    let flow = opposing_flow(&initial, 0.001);
    let mut limited = time(1e-12, 1.0);
    limited.composition.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_transport_candidate(
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            limited,
            || false,
        ),
        Err(BinaryNonisothermalAdvanceError::WorkLimit { .. })
    ));
    let mut impossible = time(1e-3, 1.0);
    impossible.composition.minimum_step_s = 0.6;
    // Valid ordered time controls, but a single CFL-admitted trial is <=.5 s.
    assert!(
        initial
            .advance_transport_candidate(
                &flow,
                1.0,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryNonisothermalSolve::default(),
                impossible,
                || false,
            )
            .is_err()
    );
    assert_eq!(initial, original);
}

#[test]
fn g4_late_cancellation_publishes_nothing_and_replays_identically() {
    let initial = fixture(
        [20.0, 60.0],
        -100.0,
        0.0,
        0.0,
        &[[0.001, 0.3, 350.0], [0.001, 0.7, 500.0]],
        &[[0, 1]],
    );
    let original = initial.clone();
    let flow = opposing_flow(&initial, 0.0004);
    let controls = time(1e-4, 0.02);
    let mut polls = 0;
    let expected = initial
        .advance_transport_candidate(
            &flow,
            0.2,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            controls,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        initial.advance_transport_candidate(
            &flow,
            0.2,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            controls,
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(BinaryNonisothermalAdvanceError::Solver(
            BinaryNonisothermalError::Cancelled
        ))
    ));
    assert_eq!(initial, original);
    assert_eq!(advance(&initial, &flow, 0.2, controls), expected);
}

#[test]
fn g0_continuation_accepts_short_final_remainders_with_one_exact_clock() {
    let initial = fixture(
        [20.0, 60.0],
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.3, 350.0], [0.001, 0.7, 550.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0002);
    let mut controls = time(1.0, 100.0);
    controls.composition.initial_step_s = 0.125;
    controls.composition.minimum_step_s = 0.125;
    controls.composition.maximum_step_s = 0.125;
    let first = advance(&initial, &flow, 0.3, controls);
    assert_eq!(first.advance().accepted_intervals, 3);
    near(first.advance().minimum_accepted_step_s, 0.05, 2e-15);
    assert_eq!(first.advance().maximum_accepted_step_s, 0.125);
    assert_eq!(first.field().elapsed_s(), 0.3);
    let second = advance(first.field(), &flow, 0.07, controls);
    assert_eq!(second.advance().accepted_intervals, 1);
    assert_eq!(second.advance().duration_s, 0.07);
    assert_eq!(
        second.field().elapsed_s(),
        first.field().elapsed_s() + 0.07
    );
    assert_eq!(advance(first.field(), &flow, 0.07, controls), second);
}
