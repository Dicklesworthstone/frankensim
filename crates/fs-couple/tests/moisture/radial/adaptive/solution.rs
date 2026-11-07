//! Physical finite-bath feedback and real unequal-shell temporal error scales.

use super::*;
use crate::{atmosphere as air, solution as bath};
use fs_couple::moisture::solution::FiniteMoistureSolution;

const WATER_DENSITY: f64 = 0.018 / bath::WATER_VOLUME;

fn check_account<Step>(
    original: &MoistureRadial,
    bath: &FiniteMoistureSolution,
    next: &MoistureRadial,
    next_bath: &FiniteMoistureSolution,
    report: &MoistureAdvanceReport<Step, 2>,
) {
    near(
        next_bath.solvent_mass_kg() + next.water_mass_kg(),
        bath.solvent_mass_kg() + original.water_mass_kg(),
        1e-10,
    );
    near(report.stored_water_kg[0], original.water_mass_kg(), 1e-13);
    near(report.stored_water_kg[1], next.water_mass_kg(), 1e-13);
    assert_eq!(
        report.reservoir_water_kg.unwrap(),
        [bath.solvent_mass_kg(), next_bath.solvent_mass_kg()],
    );
    assert_eq!(
        report.activity.unwrap(),
        [
            bath.solution().water_activity().unwrap(),
            next_bath.solution().water_activity().unwrap(),
        ],
    );
    let changed = original
        .moisture()
        .iter()
        .zip(next.moisture())
        .zip(original.dry_masses_kg())
        .map(|((&old, &new), &mass)| mass * (new - old))
        .sum::<f64>();
    near(report.water_change_kg, changed, 1e-12);
    near(report.boundary_water_kg[0], 0.0, 1e-13);
    near(report.boundary_water_kg[1], changed, 1e-10);
    near(report.transport_balance_residual_kg, 0.0, 1e-10);
    near(report.water_balance_residual_kg, 0.0, 1e-10);
    assert_eq!(
        next_bath.solution().solute_moles(),
        bath.solution().solute_moles()
    );
    assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
}

fn finite_bath_reference(
    symmetry: RadialMoistureSymmetry,
    duration: f64,
    steps: usize,
) -> ([f64; 2], f64) {
    let (mass, _) = two_shell_parameters(symmetry);
    let (interior, boundary) = match symmetry {
        RadialMoistureSymmetry::Sphere => (
            500.0 * 0.001 * 4.0 * core::f64::consts::PI / (1.0 / 0.0275 - 1.0 / 0.0675),
            500.0 * 0.001 * 4.0 * core::f64::consts::PI / (1.0 / 0.0675 - 1.0 / 0.1),
        ),
        RadialMoistureSymmetry::Cylinder { length_m } => (
            500.0 * 0.001 * 2.0 * core::f64::consts::PI * length_m
                / (0.0675_f64 / 0.0275).ln(),
            500.0 * 0.001 * 2.0 * core::f64::consts::PI * length_m
                / (0.1_f64 / 0.0675).ln(),
        ),
    };
    let activity = |m: [f64; 2]| {
        let volume = 1e-4
            - (mass[0] * (m[0] - 0.22) + mass[1] * (m[1] - 0.08))
                / WATER_DENSITY;
        (-bath::WATER_VOLUME * 0.8 / volume).exp()
    };
    // Independent continuous water equations. Eliminate actual solvent volume
    // with the closed mass balance, then integrate the resulting two-state ODE.
    let rate = |m: [f64; 2]| {
        let exchange = interior * (m[1] - m[0]);
        [
            exchange / mass[0],
            (-exchange + boundary * (0.25 * activity(m) - m[1])) / mass[1],
        ]
    };
    let h = duration / steps as f64;
    let mut m = [0.22, 0.08];
    for _ in 0..steps {
        let a = rate(m);
        let b = rate(std::array::from_fn(|i| m[i] + 0.5 * h * a[i]));
        let c = rate(std::array::from_fn(|i| m[i] + 0.5 * h * b[i]));
        let d = rate(std::array::from_fn(|i| m[i] + h * c[i]));
        m = std::array::from_fn(|i| {
            m[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    (m, activity(m))
}

#[test]
fn radial_solution_adaptive_refines_actual_curved_fields_and_finite_bath_activity() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            let original = specimen(
                symmetry,
                vec![0.02, 0.035, 0.1],
                vec![0.22, 0.08],
                [0.001; 2],
                cx,
            );
            let bath = bath::reservoir(1e-4, vec![0.8]);
            let saved = bath.solution().clone();
            let faces = [None, Some(air::sorption(0.25))];
            let duration = 0.1;
            let (reference, activity) = finite_bath_reference(symmetry, duration, 4096);
            let (checked, checked_activity) = finite_bath_reference(symmetry, duration, 8192);
            for i in 0..2 {
                near(reference[i], checked[i], 5e-13);
            }
            near(activity, checked_activity, 5e-13);
            for second_order in [false, true] {
                let mut errors = Vec::new();
                let mut work = Vec::new();
                for tolerance in [1e-4, 1e-6] {
                    let time = controls(duration, tolerance);
                    let (next_bath, next, evaluations) = if second_order {
                        let (r, m, report) = bath
                            .advance_radial_sdirk2_adaptive_to(
                                &original,
                                duration,
                                &faces,
                                config(),
                                nonlinear(),
                                bath::controls(1e-4),
                                time,
                                cx,
                            )
                            .unwrap();
                        check_account(&original, &bath, &m, &r, &report);
                        (r, m, report.evaluations)
                    } else {
                        let (r, m, report) = bath
                            .advance_radial_adaptive_to(
                                &original,
                                duration,
                                &faces,
                                config(),
                                nonlinear(),
                                bath::controls(1e-4),
                                time,
                                cx,
                            )
                            .unwrap();
                        check_account(&original, &bath, &m, &r, &report);
                        (r, m, report.evaluations)
                    };
                    let moisture_error = next
                        .moisture()
                        .iter()
                        .zip(reference)
                        .map(|(&a, b)| (a - b).abs())
                        .fold(0.0_f64, f64::max);
                    errors.push(
                        moisture_error.max(
                            (next_bath.solution().water_activity().unwrap() - activity).abs(),
                        ),
                    );
                    work.push(evaluations);
                    assert_eq!(next.time_s(), duration);
                    assert!(next_bath.solvent_mass_kg() < bath.solvent_mass_kg());
                    assert!(
                        next_bath.solution().water_activity().unwrap()
                            < bath.solution().water_activity().unwrap()
                    );
                }
                assert!(
                    errors[1] < 0.3 * errors[0],
                    "{symmetry:?}, second_order={second_order}: {errors:?}"
                );
                assert!(
                    errors[1] < if second_order { 8e-6 } else { 8e-5 },
                    "{errors:?}"
                );
                assert!(work[1] > work[0]);
            }
            assert_eq!(bath.solution(), &saved);
            assert_eq!(original.moisture(), &[0.22, 0.08]);
            assert_eq!(original.time_s(), 0.0);
        }
    });
}

fn decay(z: f64, second_order: bool) -> f64 {
    if second_order {
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
    } else {
        1.0 / (1.0 + z)
    }
}

#[test]
fn sealed_solution_radial_error_keeps_opposing_shell_water_discrepancies() {
    with_cx(|cx| {
        let symmetry = RadialMoistureSymmetry::Sphere;
        let original = specimen(
            symmetry,
            vec![0.02, 0.035, 0.1],
            vec![0.22, 0.08],
            [0.001; 2],
            cx,
        );
        let (mass, rate) = two_shell_parameters(symmetry);
        assert!(mass[1] > 25.0 * mass[0]);
        let bath = bath::reservoir(1e-4, vec![0.8]);
        let duration = 0.05;
        for second_order in [false, true] {
            let coarse = two_shell_mode(mass, decay(rate * duration, second_order));
            let fine = two_shell_mode(
                mass,
                decay(0.5 * rate * duration, second_order).powi(2),
            );
            let factor = if second_order { 1.0 / 3.0 } else { 1.0 };
            let discrepancy =
                mass[0] * (coarse[0] - fine[0]).abs()
                    + mass[1] * (coarse[1] - fine[1]).abs();
            let mut time = controls(duration, 1.0);
            time.reservoir_water_absolute_tolerance_kg = 2.0 * factor * discrepancy;
            if second_order {
                let (r, m, report) = bath
                    .advance_radial_sdirk2_adaptive_to(
                        &original,
                        duration,
                        &[None, None],
                        config(),
                        nonlinear(),
                        bath::controls(1e-4),
                        time,
                        cx,
                    )
                    .unwrap();
                check_account(&original, &bath, &m, &r, &report);
                near(report.maximum_accepted_error_ratio, 0.5, 3e-8);
                assert_eq!(report.accepted_intervals, 1);
                assert_eq!(report.evaluations, 6);
                near(report.water_change_kg, 0.0, 1e-12);
            } else {
                let (r, m, report) = bath
                    .advance_radial_adaptive_to(
                        &original,
                        duration,
                        &[None, None],
                        config(),
                        nonlinear(),
                        bath::controls(1e-4),
                        time,
                        cx,
                    )
                    .unwrap();
                check_account(&original, &bath, &m, &r, &report);
                near(report.maximum_accepted_error_ratio, 0.5, 3e-8);
                assert_eq!(report.accepted_intervals, 1);
                assert_eq!(report.evaluations, 3);
                near(report.water_change_kg, 0.0, 1e-12);
            }
        }
    });
}

#[test]
fn pure_solution_radial_tracks_finite_water_when_activity_cannot_change() {
    with_cx(|cx| {
        let radius: f64 = 0.05;
        let rate = 3.0 * 0.001 / radius.powi(2);
        let duration = 0.4 / rate;
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.1],
            [0.001; 2],
            cx,
        );
        let bath = bath::reservoir(5e-4, vec![]);
        let faces = [None, Some(air::sorption(0.25))];
        for second_order in [false, true] {
            let coarse = 0.25 + (0.1 - 0.25) * decay(rate * duration, second_order);
            let fine =
                0.25 + (0.1 - 0.25) * decay(0.5 * rate * duration, second_order).powi(2);
            let factor = if second_order { 1.0 / 3.0 } else { 1.0 };
            let mut time = controls(duration, 1.0);
            time.reservoir_water_absolute_tolerance_kg =
                5.0 * factor * original.dry_mass_kg() * (coarse - fine).abs();
            let inspect = |r: &FiniteMoistureSolution, m: &MoistureRadial| {
                near(m.moisture()[0], fine, 3e-12);
                assert_eq!(r.solution().water_activity().unwrap(), 1.0);
                assert!(r.solvent_mass_kg() < bath.solvent_mass_kg());
            };
            if second_order {
                let (r, m, report) = bath
                    .advance_radial_sdirk2_adaptive_to(
                        &original,
                        duration,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(5e-4),
                        time,
                        cx,
                    )
                    .unwrap();
                check_account(&original, &bath, &m, &r, &report);
                inspect(&r, &m);
                near(report.maximum_accepted_error_ratio, 0.2, 3e-7);
                assert_eq!(report.accepted_intervals, 1);
                assert_eq!(report.evaluations, 6);
            } else {
                let (r, m, report) = bath
                    .advance_radial_adaptive_to(
                        &original,
                        duration,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(5e-4),
                        time,
                        cx,
                    )
                    .unwrap();
                check_account(&original, &bath, &m, &r, &report);
                inspect(&r, &m);
                near(report.maximum_accepted_error_ratio, 0.2, 3e-7);
                assert_eq!(report.accepted_intervals, 1);
                assert_eq!(report.evaluations, 3);
            }
        }
    });
}

#[test]
fn radial_solution_adaptive_refusal_discards_an_accepted_prefix_and_cancellation_is_atomic() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.05],
            vec![0.1],
            [0.001; 2],
            cx,
        );
        let bath = bath::reservoir(1e-4, vec![0.2]);
        let saved = bath.solution().clone();
        let faces = [None, Some(air::sorption(0.25))];
        let step = 1.0 / 256.0;
        for second_order in [false, true] {
            let budget = if second_order { 6 } else { 3 };
            let mut time = controls(step, 1.0);
            time.maximum_evaluations = budget;
            if second_order {
                assert!(matches!(
                    bath.advance_radial_sdirk2_adaptive_to(
                        &original, 2.0 * step, &faces, config(), nonlinear(),
                        bath::controls(1e-4), time, cx,
                    ),
                    Err(MoistureAdvanceError::WorkLimit {
                        time_s, accepted_intervals: 1, evaluations
                    }) if time_s == step && evaluations == budget
                ));
            } else {
                assert!(matches!(
                    bath.advance_radial_adaptive_to(
                        &original, 2.0 * step, &faces, config(), nonlinear(),
                        bath::controls(1e-4), time, cx,
                    ),
                    Err(MoistureAdvanceError::WorkLimit {
                        time_s, accepted_intervals: 1, evaluations
                    }) if time_s == step && evaluations == budget
                ));
            }
            assert_eq!(bath.solution(), &saved);
            assert_eq!(original.moisture(), &[0.1]);
            assert_eq!(original.time_s(), 0.0);
        }
    });
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.05],
            vec![0.1],
            [0.001; 2],
            &cx,
        );
        let bath = bath::reservoir(1e-4, vec![0.2]);
        let saved = bath.solution().clone();
        gate.request();
        for second_order in [false, true] {
            let faces = [None, Some(air::sorption(0.25))];
            if second_order {
                assert!(matches!(
                    bath.advance_radial_sdirk2_adaptive_to(
                        &original, 1.0, &faces, config(), nonlinear(),
                        bath::controls(1e-4), controls(1.0, 1e-5), &cx,
                    ),
                    Err(MoistureAdvanceError::Cancelled)
                ));
            } else {
                assert!(matches!(
                    bath.advance_radial_adaptive_to(
                        &original, 1.0, &faces, config(), nonlinear(),
                        bath::controls(1e-4), controls(1.0, 1e-5), &cx,
                    ),
                    Err(MoistureAdvanceError::Cancelled)
                ));
            }
        }
        assert_eq!(bath.solution(), &saved);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
