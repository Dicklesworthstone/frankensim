use super::*;
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;

fn field(amounts: Vec<Vec<f64>>) -> MulticomponentMiscibilityField {
    let n = amounts[0].len();
    let species: Vec<_> = (0..n)
        .map(|i| SpeciesId::new(format!("species_{i}")).unwrap())
        .collect();
    let model = MulticomponentRegularSolutionV1::new(
        species.clone(),
        vec![0.0; n * n],
        vec![0.0; n * n],
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let transport = MaxwellStefanV1::new(
        species
            .into_iter()
            .enumerate()
            .map(|(i, id)| (id, MolarMass::new(0.018 + 0.01 * i as f64)))
            .collect(),
    )
    .unwrap();
    let mut diffusion = vec![1.0; n * n];
    for i in 0..n {
        diffusion[i * n + i] = 0.0;
    }
    let cells = amounts
        .iter()
        .map(|cell| MulticomponentMiscibilityCell {
            volume_m3: cell.iter().sum(),
        })
        .collect();
    MulticomponentMiscibilityNetwork::new(
        model,
        transport,
        Temperature::new(300.0),
        AmountConcentration::new(1.0),
        cells,
        vec![MulticomponentMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_matrix: 0,
        }],
        vec![diffusion],
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
    .bind(amounts, || false)
    .unwrap()
}

fn time(step: f64) -> MulticomponentMiscibilityTimeControl {
    MulticomponentMiscibilityTimeControl {
        initial_step_s: step,
        minimum_step_s: 1e-9,
        maximum_step_s: step,
        absolute_mole_fraction_tolerance: 1e-8,
        relative_mole_fraction_tolerance: 1e-7,
        maximum_attempts: 20_000,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

#[test]
fn small_binary_perturbation_matches_independent_diffusive_relaxation() {
    // For ideal binary liquids, the exact linearization of this finite-volume
    // graph is d(delta)/dt=-2*D*delta. Finite-jump nonlinear terms are O(delta^3).
    let perturbation = 1e-4;
    let initial = field(vec![
        vec![0.5 + perturbation, 0.5 - perturbation],
        vec![0.5 - perturbation, 0.5 + perturbation],
    ]);
    let mut accuracy = time(0.4);
    accuracy.absolute_mole_fraction_tolerance = 1e-13;
    accuracy.relative_mole_fraction_tolerance = 1e-13;
    let candidate = initial
        .advance_candidate(
            0.4,
            accuracy,
            MulticomponentMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    let expected = 0.5 + perturbation * (-0.8_f64).exp();
    near(candidate.field().component_moles()[0], expected, 3e-9);
    assert!(candidate.advance().rejected_intervals > 0);
    assert!(candidate.advance().maximum_accepted_error_ratio <= 1.0);
    assert!(candidate.advance().free_energy_change_j < 0.0);
    assert_eq!(candidate.field().elapsed_s(), 0.4);
    assert_eq!(initial.elapsed_s(), 0.0);
}

#[test]
fn ternary_receipts_include_only_the_two_accepted_half_steps() {
    let initial = field(vec![vec![0.55, 0.3, 0.15], vec![0.2, 0.35, 0.45]]);
    let controls = MulticomponentMiscibilitySolve::default();
    let duration = 0.03;
    let first = initial
        .step_candidate(duration * 0.5, controls, || false)
        .unwrap();
    let second = first
        .field()
        .step_candidate(duration * 0.5, controls, || false)
        .unwrap();
    let mut accuracy = time(duration);
    accuracy.absolute_mole_fraction_tolerance = 1.0;
    accuracy.maximum_attempts = 1;
    let candidate = initial
        .advance_candidate(duration, accuracy, controls, || false)
        .unwrap();
    assert_eq!(candidate.field(), second.field());
    assert_eq!(candidate.advance().completed_steps, 3);
    assert_eq!(candidate.advance().completed_maxwell_stefan_evaluations, 2);
    assert_eq!(candidate.advance().accepted_intervals, 1);
    for species in 0..3 {
        let transfer =
            first.step().transferred_mol[0][species] + second.step().transferred_mol[0][species];
        near(
            candidate.advance().transferred_mol[0][species],
            transfer,
            1e-17,
        );
        near(
            initial.component_moles()[species] - candidate.field().component_moles()[species],
            transfer,
            1e-16,
        );
        near(
            candidate.advance().component_residual_mol[species],
            0.0,
            1e-15,
        );
    }
}

#[test]
fn oversize_physical_interval_refines_without_clipping_species() {
    let initial = field(vec![vec![0.8, 0.15, 0.05], vec![0.05, 0.25, 0.7]]);
    let controls = MulticomponentMiscibilitySolve::default();
    assert!(matches!(
        initial.step_candidate(3.0, controls, || false),
        Err(MulticomponentMiscibilityError::Positivity { .. })
    ));
    let mut accuracy = time(3.0);
    accuracy.absolute_mole_fraction_tolerance = 1e-6;
    accuracy.relative_mole_fraction_tolerance = 1e-5;
    let candidate = initial
        .advance_candidate(3.0, accuracy, controls, || false)
        .unwrap();
    let replay = initial
        .advance_candidate(3.0, accuracy, controls, || false)
        .unwrap();
    assert_eq!(candidate, replay);
    assert!(candidate.advance().rejected_intervals > 0);
    // A rejected timestep keeps the same original field and its cached rates.
    // Even a completed rejected fine trial requires only a new midpoint solve.
    assert!(
        candidate.advance().completed_maxwell_stefan_evaluations
            <= candidate.advance().accepted_intervals + candidate.advance().attempted_intervals
    );
    assert!(
        candidate.advance().completed_maxwell_stefan_evaluations
            < 2 * candidate.advance().attempted_intervals
    );
    assert!(
        candidate
            .field()
            .component_moles()
            .iter()
            .all(|&amount| amount > 0.0)
    );
    assert_eq!(candidate.field().elapsed_s(), 3.0);
    for species in 0..3 {
        let expected = initial.component_moles()[species] + initial.component_moles()[3 + species];
        let actual = candidate.field().component_moles()[species]
            + candidate.field().component_moles()[3 + species];
        near(actual, expected, 2e-14);
        near(
            candidate.advance().component_residual_mol[species],
            0.0,
            1e-14,
        );
    }
    assert!(candidate.advance().free_energy_change_j < 0.0);
}

#[test]
fn work_accuracy_and_late_cancellation_leave_original_field_intact() {
    let initial = field(vec![vec![0.7, 0.2, 0.1], vec![0.1, 0.2, 0.7]]);
    let snapshot = initial.clone();
    let controls = MulticomponentMiscibilitySolve::default();
    let mut accuracy = time(0.2);
    accuracy.absolute_mole_fraction_tolerance = 1e-15;
    accuracy.relative_mole_fraction_tolerance = 1e-15;
    accuracy.maximum_attempts = 1;
    assert!(matches!(
        initial.advance_candidate(0.2, accuracy, controls, || false),
        Err(MulticomponentMiscibilityAdvanceError::WorkLimit { .. })
    ));
    accuracy.minimum_step_s = 0.2;
    assert!(matches!(
        initial.advance_candidate(0.2, accuracy, controls, || false),
        Err(MulticomponentMiscibilityAdvanceError::TemporalAccuracy { .. })
    ));
    let mut polls = 0usize;
    let cancelled = initial.advance_candidate(0.2, time(0.2), controls, || {
        polls += 1;
        polls > 40
    });
    assert!(matches!(
        cancelled,
        Err(MulticomponentMiscibilityAdvanceError::Solver(
            MulticomponentMiscibilityError::Cancelled
        ))
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn independent_trace_low_words_cannot_hide_from_temporal_error() {
    let a = field(vec![vec![1e-300, 0.4, 0.6], vec![1e-300, 0.4, 0.6]]);
    let mut b = a.clone();
    b.amount_remainders[0] = 1e-310;
    let mut accuracy = time(0.1);
    accuracy.absolute_mole_fraction_tolerance = 0.0;
    accuracy.relative_mole_fraction_tolerance = 1e-100;
    assert_eq!(error_ratio(&a, &a, accuracy, &mut || false).unwrap(), 0.0);
    assert_eq!(
        error_ratio(&a, &b, accuracy, &mut || false).unwrap(),
        f64::INFINITY
    );
    b.amount_remainders[0] = 0.0;
    b.amount_remainders[2] = 1e-17;
    accuracy.relative_mole_fraction_tolerance = 1e-18;
    assert!(error_ratio(&a, &b, accuracy, &mut || false).unwrap() > 1.0);
}

#[test]
fn a_short_final_remainder_still_reaches_the_exact_endpoint() {
    let initial = field(vec![vec![0.2, 0.3, 0.5], vec![0.2, 0.3, 0.5]]);
    let mut accuracy = time(0.06);
    accuracy.minimum_step_s = 0.05;
    let result = initial
        .advance_candidate(
            0.07,
            accuracy,
            MulticomponentMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(result.field().elapsed_s(), 0.07);
    assert_eq!(result.advance().accepted_intervals, 2);
    assert!(result.advance().minimum_accepted_step_s < accuracy.minimum_step_s);
    assert_eq!(result.field().component_moles(), initial.component_moles());
}
