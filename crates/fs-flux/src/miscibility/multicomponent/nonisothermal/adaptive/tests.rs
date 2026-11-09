use super::super::tests::{Case, near};
use super::*;

fn time(step: f64, temperature_tolerance: f64) -> MulticomponentNonisothermalTimeControl {
    MulticomponentNonisothermalTimeControl {
        composition: MulticomponentMiscibilityTimeControl {
            initial_step_s: step,
            minimum_step_s: 1e-10,
            maximum_step_s: step,
            absolute_mole_fraction_tolerance: 1e-8,
            relative_mole_fraction_tolerance: 1e-6,
            maximum_attempts: 10000,
        },
        absolute_temperature_tolerance_k: temperature_tolerance,
        relative_temperature_tolerance: 0.0,
    }
}

#[test]
fn heat_only_temperature_accuracy_controls_refinement_against_exact_solution() {
    let case = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        conductivity: 12.0,
        ..Case::default()
    };
    let initial = case.field();
    let controls = MulticomponentNonisothermalSolve::default();
    let loose = initial
        .advance_candidate(0.5, time(0.5, 0.2), controls, || false)
        .unwrap();
    let tight = initial
        .advance_candidate(0.5, time(0.5, 1e-5), controls, || false)
        .unwrap();
    let capacities = [56.0, 112.0];
    let equilibrium = 500.0;
    let contrast = 300.0 * (-12.0_f64 * (1.0 / capacities[0] + 1.0 / capacities[1]) * 0.5).exp();
    let exact = [
        equilibrium - 2.0 * contrast / 3.0,
        equilibrium + contrast / 3.0,
    ];
    let error = |field: &MulticomponentNonisothermalField| {
        field
            .temperatures()
            .iter()
            .zip(exact)
            .map(|(t, expected)| (t.value() - expected).abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(tight.advance().accepted_intervals > loose.advance().accepted_intervals);
    assert!(error(tight.field()) < 0.02 * error(loose.field()));
    assert!(error(tight.field()) < 0.02);
    assert_eq!(tight.field().component_moles(), initial.component_moles());
    assert_eq!(tight.field().elapsed_s(), 0.5);
    assert!(tight.advance().maximum_accepted_error_ratio <= 1.0);
    assert!(tight.advance().entropy_change_j_per_k > 0.0);
    assert!(tight.advance().energy_residual_j.abs() <= controls.energy_tolerance_j);
    near(
        tight.advance().transferred_energy_j[0],
        -capacities[0] * (tight.field().temperatures()[0].value() - 300.0),
        1e-7,
    );
}

#[test]
fn complete_coupled_half_steps_and_work_receipts_are_reused_exactly() {
    let initial = Case {
        conductivity: 2.0,
        heat: [300.0, -200.0, 0.0],
        ..Case::default()
    }
    .field();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut bounds = time(0.001, 100.0);
    bounds.composition.absolute_mole_fraction_tolerance = 1.0;
    let candidate = initial
        .advance_candidate(0.001, bounds, controls, || false)
        .unwrap();
    let first = initial.step_candidate(0.0005, controls, || false).unwrap();
    let second = first
        .field()
        .step_candidate(0.0005, controls, || false)
        .unwrap();
    assert_eq!(candidate.field(), second.field());
    assert_eq!(candidate.advance().attempted_intervals, 1);
    assert_eq!(candidate.advance().accepted_intervals, 1);
    assert_eq!(candidate.advance().completed_steps, 3);
    assert_eq!(candidate.advance().completed_maxwell_stefan_evaluations, 2);
    for species in 0..3 {
        near(
            candidate.advance().transferred_mol[0][species],
            first.step().transferred_mol[0][species] + second.step().transferred_mol[0][species],
            1e-20,
        );
    }
    near(
        candidate.advance().transferred_energy_j[0],
        first.step().transferred_energy_j[0] + second.step().transferred_energy_j[0],
        1e-15,
    );
}

#[test]
fn coupled_refinement_keeps_all_species_energy_and_physical_receipts() {
    let initial = Case {
        conductivity: 4.0,
        heat: [1200.0, -600.0, 0.0],
        ..Case::default()
    }
    .field();
    let controls = MulticomponentNonisothermalSolve::default();
    let candidate = initial
        .advance_candidate(0.4, time(0.4, 1e-4), controls, || false)
        .unwrap();
    assert!(candidate.advance().accepted_intervals > 1);
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    assert!(
        candidate
            .advance()
            .component_residual_mol
            .iter()
            .all(|x| x.abs() <= controls.amount_balance_tolerance_mol)
    );
    assert!(
        candidate.advance().maximum_capacity_residual_mol <= controls.amount_balance_tolerance_mol
    );
    assert!(candidate.advance().energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(candidate.advance().entropy_change_j_per_k >= 0.0);
    for species in 0..3 {
        near(
            candidate.field().component_moles()[species] - initial.component_moles()[species],
            -candidate.advance().transferred_mol[0][species],
            2e-15,
        );
        near(
            candidate.field().component_moles()[3 + species]
                - initial.component_moles()[3 + species],
            candidate.advance().transferred_mol[0][species],
            2e-15,
        );
    }
    let mut reference = initial.clone();
    for _ in 0..400 {
        reference = reference
            .step_candidate(0.001, controls, || false)
            .unwrap()
            .into_parts()
            .0;
    }
    for (a, b) in candidate
        .field()
        .temperatures()
        .iter()
        .zip(reference.temperatures())
    {
        near(a.value(), b.value(), 0.015);
    }
    for (&a, &b) in candidate
        .field()
        .mole_fractions()
        .iter()
        .zip(reference.mole_fractions())
    {
        near(a, b, 1e-5);
    }
}

#[test]
fn failed_trials_and_late_cancellation_never_publish_partial_heat_or_species() {
    let initial = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        volumes: vec![0.001; 2],
        conductivity: 100.0,
        ..Case::default()
    }
    .field();
    let before = initial.clone();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut limited = time(0.6, 1e-5);
    limited.composition.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_candidate(0.6, limited, controls, || false),
        Err(MulticomponentNonisothermalAdvanceError::WorkLimit {
            accepted_intervals: 0,
            ..
        })
    ));
    assert_eq!(initial, before);
    let bounds = time(0.01, 0.01);
    let mut calls = 0usize;
    let expected = initial
        .advance_candidate(0.01, bounds, controls, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut seen = 0usize;
    assert!(matches!(
        initial.advance_candidate(0.01, bounds, controls, || {
            seen += 1;
            seen == calls
        }),
        Err(MulticomponentNonisothermalAdvanceError::Solver(
            MulticomponentNonisothermalError::Cancelled
        ))
    ));
    assert_eq!(initial, before);
    assert_eq!(
        initial
            .advance_candidate(0.01, bounds, controls, || false)
            .unwrap(),
        expected
    );
}

#[test]
fn final_short_interval_and_continued_clock_share_exact_endpoints() {
    let initial = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        temperatures: vec![300.0; 2],
        ..Case::default()
    }
    .field();
    let controls = MulticomponentNonisothermalSolve::default();
    let continued = initial
        .step_candidate(8.0, controls, || false)
        .unwrap()
        .into_parts()
        .0;
    let mut bounds = time(0.1, 1e-6);
    bounds.composition.minimum_step_s = 0.05;
    let candidate = continued
        .advance_candidate(0.13, bounds, controls, || false)
        .unwrap();
    assert_eq!(candidate.field().elapsed_s(), 8.0 + 0.13);
    assert_eq!(
        candidate.field().component_moles(),
        continued.component_moles()
    );
    assert_eq!(candidate.field().temperatures(), continued.temperatures());
    assert_eq!(candidate.advance().accepted_intervals, 2);
    assert!(candidate.advance().minimum_accepted_step_s < bounds.composition.minimum_step_s);
    assert_eq!(candidate.advance().energy_residual_j, 0.0);
}

#[test]
fn temporal_errors_include_trace_low_words_and_zero_kelvin_allowances_refuse() {
    let mut a = Case::default().field();
    let mut b = a.clone();
    a.state.amounts[0] = 1e-250;
    b.state.amounts[0] = 1e-250;
    a.state.fractions[0] = 1e-250;
    b.state.fractions[0] = 1e-250;
    b.state.amount_remainders[0] = 1e-270;
    let mut bounds = time(0.01, 1e-6);
    bounds.composition.absolute_mole_fraction_tolerance = 0.0;
    bounds.composition.relative_mole_fraction_tolerance = 1e-25;
    assert!(error_ratio(&a, &b, bounds, &mut || false).unwrap() > 10000.0);
    bounds.absolute_temperature_tolerance_k = 0.0;
    bounds.relative_temperature_tolerance = 0.0;
    assert!(
        a.advance_candidate(
            0.01,
            bounds,
            MulticomponentNonisothermalSolve::default(),
            || false
        )
        .is_err()
    );
}
