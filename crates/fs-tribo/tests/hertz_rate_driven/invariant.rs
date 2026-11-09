//! Independent analytic pressure-moment events, original static continuation,
//! clipped phases, weak-slope limits and bounded exact-path refusal.
use super::*;

// Test-only closed forms use the platform exp/log independently of production's
// deterministic decay/series kernels. Fixtures stay away from their cancellation
// regimes; those limits have separate constant/linear analytical references.
fn moments(initial: f64, a: f64, b: f64, t: f64) -> (f64, f64, f64) {
    if b == 0.0 {
        return (
            initial + a * t,
            initial * t + 0.5 * a * t * t,
            initial * initial * t + initial * a * t * t + a * a * t * t * t / 3.0,
        );
    }
    let equilibrium = a / b;
    let delta = initial - equilibrium;
    let e = (-b * t).exp();
    (
        equilibrium + delta * e,
        equilibrium * t + delta * (1.0 - e) / b,
        equilibrium * equilibrium * t
            + 2.0 * equilibrium * delta * (1.0 - e) / b
            + delta * delta * (1.0 - e * e) / (2.0 * b),
    )
}

#[test]
fn original_affine_translation_and_spin_stop_then_hold_with_exact_heat_and_impulses() {
    let (mu0, slope, initial) = (0.35, 0.4, 0.7);
    for axis in [0, 2] {
        let (dry, viscous) = if axis == 0 {
            (N * mu0, N * slope)
        } else {
            (N * K * mu0, 0.4 * N * slope)
        };
        let mobility = diagonal()[axis][axis];
        let demand = 0.1;
        let a = mobility * (demand - dry);
        let b = mobility * viscous;
        let stop = (1.0 + b * initial / (-a)).ln() / b;
        let (_, distance, square) = moments(initial, a, b, stop);
        for sign in [-1.0, 1.0] {
            let r = run(
                affine(mu0, slope),
                scalar(axis, sign * initial, sign * demand, stop + 0.4),
            );
            assert!(r.exact_invariant);
            assert!(!r.adaptive_error_estimates_available);
            assert_eq!(r.trajectory.mode, HertzDrivenMode::Resting);
            assert_eq!(q(&r), [0.0; 3]);
            assert!(r.trajectory.arrest_time_is_exact);
            let actual_stop = r.trajectory.first_arrest_time_s.unwrap();
            close(actual_stop, stop, 3e-14);
            assert_eq!(r.stopping_time_bracket_s, Some([actual_stop; 2]));
            close(travel(&r)[axis], sign * distance, 4e-14);
            close(
                impulse(&r)[axis],
                -sign * (dry * stop + viscous * distance + demand * 0.4),
                5e-14,
            );
            close(
                r.trajectory.dissipated_energy_j,
                dry * distance + viscous * square,
                8e-14,
            );
            close(r.trajectory.drive_work_j, demand * distance, 8e-14);
            close(
                r.trajectory.kinetic_energy_change_j,
                -initial * initial / (2.0 * mobility),
                8e-14,
            );
            close(r.trajectory.energy_residual_j, 0.0, 8e-14);
        }
    }
}

#[test]
fn over_static_drive_reverses_after_localized_zero_and_phase_splitting_preserves_work() {
    let (mu0, slope, initial) = (0.35, 0.4, 0.6);
    for axis in [0, 2] {
        let mobility = diagonal()[axis][axis];
        let (dry, viscous, demand) = if axis == 0 {
            (N * mu0, N * slope, -2.0)
        } else {
            (N * K * mu0, 0.4 * N * slope, -1.2)
        };
        let b = mobility * viscous;
        let before = mobility * (demand - dry);
        let stop = (1.0 + b * initial / (-before)).ln() / b;
        let (_, x_before, q2_before) = moments(initial, before, b, stop);
        let (reversed, x_after, q2_after) = moments(0.0, mobility * (-demand - dry), b, 0.7);
        let input = scalar(axis, initial, demand, stop + 0.7);
        let patch = HertzFrictionPatch::new(N, A, affine(mu0, slope)).unwrap();
        let r = patch
            .advance_driven(&input, config(), grid(1), controls(), || false)
            .unwrap();
        assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
        assert!(r.exact_invariant && r.trajectory.arrest_time_is_exact);
        close(r.trajectory.first_arrest_time_s.unwrap(), stop, 3e-14);
        close(q(&r)[axis], -reversed, 5e-14);
        close(travel(&r)[axis], x_before - x_after, 5e-14);
        close(
            r.trajectory.dissipated_energy_j,
            dry * (x_before + x_after) + viscous * (q2_before + q2_after),
            2e-13,
        );
        close(
            r.trajectory.drive_work_j,
            demand * (x_before - x_after),
            2e-13,
        );
        close(r.trajectory.energy_residual_j, 0.0, 2e-13);
        let split_time = 0.5 * stop;
        let first = patch
            .advance_driven(
                &HertzDrivenInput {
                    dt_s: split_time,
                    ..input
                },
                config(),
                grid(1),
                controls(),
                || false,
            )
            .unwrap();
        let continuation = HertzDrivenInput {
            velocity_m_s: first.trajectory.velocity_after_m_s,
            spin_rate_rad_s: first.trajectory.spin_rate_after_rad_s,
            dt_s: input.dt_s - split_time,
            mode: first.trajectory.mode,
            ..input
        };
        let second = patch
            .advance_driven(&continuation, config(), grid(1), controls(), || false)
            .unwrap();
        close(q(&second)[axis], q(&r)[axis], 8e-14);
        close(
            split_time + second.trajectory.first_arrest_time_s.unwrap(),
            stop,
            8e-14,
        );
        close(
            first.trajectory.dissipated_energy_j + second.trajectory.dissipated_energy_j,
            r.trajectory.dissipated_energy_j,
            2e-13,
        );
        close(
            first.trajectory.drive_work_j + second.trajectory.drive_work_j,
            r.trajectory.drive_work_j,
            2e-13,
        );
    }
}

#[test]
fn clipped_translation_crosses_zero_reverses_and_exits_the_active_region() {
    let (mu0, slope, initial, demand) = (0.4, -1.0, 0.6, -2.0);
    let mobility = diagonal()[0][0];
    let cutoff = 0.4;
    let first_coast = (initial - cutoff) / (-mobility * demand);
    let b = mobility * N * slope;
    let a = mobility * (demand - N * mu0);
    let braking = (1.0 + b * cutoff / (-a)).ln() / b;
    let (_, x_brake, q2_brake) = moments(cutoff, a, b, braking);
    let reverse_a = mobility * (-demand - N * mu0);
    let active_reverse = (1.0 + (-b) * cutoff / reverse_a).ln() / (-b);
    let (_, x_reverse, q2_reverse) = moments(0.0, reverse_a, b, active_reverse);
    let dt = 1.0;
    let last_coast = dt - first_coast - braking - active_reverse;
    assert!(last_coast > 0.0);
    let expected = -(cutoff + (-mobility * demand) * last_coast);
    let distance = 0.5 * (initial + cutoff) * first_coast + x_brake
        - x_reverse
        - cutoff * last_coast
        - 0.5 * (-mobility * demand) * last_coast * last_coast;
    let heat = N * mu0 * (x_brake + x_reverse) + N * slope * (q2_brake + q2_reverse);
    let r = run(affine(mu0, slope), scalar(0, initial, demand, dt));
    assert!(r.exact_invariant);
    assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
    close(
        r.trajectory.first_arrest_time_s.unwrap(),
        first_coast + braking,
        8e-14,
    );
    close(q(&r)[0], expected, 2e-13);
    close(travel(&r)[0], distance, 2e-13);
    close(r.trajectory.dissipated_energy_j, heat, 2e-13);
    close(r.trajectory.drive_work_j, demand * distance, 3e-13);
    close(r.trajectory.energy_residual_j, 0.0, 3e-13);
    assert!(r.trajectory.accepted_steps <= 6);
}

#[test]
fn equilibrium_asymptote_and_weak_slope_limits_do_not_manufacture_arrest() {
    let mu0 = 0.6;
    let demand = 0.7;
    let slope = -1.0;
    let equilibrium = (N * mu0 - demand) / (N * (-slope));
    let steady = run(affine(mu0, slope), scalar(0, equilibrium, demand, 2.0));
    close(q(&steady)[0], equilibrium, 2e-13);
    assert_eq!(steady.trajectory.first_arrest_time_s, None);
    let viscous = run(affine(0.0, 0.3), scalar(2, 0.7, 0.0, 12.0));
    let expected = 0.7 * (-diagonal()[2][2] * 0.4 * N * 0.3 * 12.0).exp();
    close(q(&viscous)[2], expected, 2e-14);
    assert_eq!(viscous.trajectory.first_arrest_time_s, None);
    assert_eq!(viscous.trajectory.mode, HertzDrivenMode::Sliding);
    assert!(q(&viscous)[2] > 0.0);
    let patch = HertzFrictionPatch::new(N, A, affine(0.0, 0.3)).unwrap();
    assert!(matches!(
        patch.advance_driven(
            &scalar(0, 0.7, 0.0, 1000.0),
            config(),
            grid(1),
            controls(),
            || false
        ),
        Err(HertzRateDrivenError::UnresolvedDrivenEvent { .. })
    ));
    // Tiny beta*t must tend smoothly to constant acceleration and its exact
    // triangular travel/heat; no cancellation in t-phi1 may erase them.
    let initial = 0.7;
    let dry = N * 0.4;
    let lambda = diagonal()[0][0];
    for slope in [-1e-14, 1e-14] {
        let stop = initial / (lambda * dry);
        let r = run(affine(0.4, slope), scalar(0, initial, 0.0, 2.0 * stop));
        close(r.trajectory.first_arrest_time_s.unwrap(), stop, 5e-14);
        close(travel(&r)[0], 0.5 * initial * stop, 5e-14);
        close(
            r.trajectory.dissipated_energy_j,
            initial * initial / (2.0 * lambda),
            5e-14,
        );
    }
}

#[test]
fn original_affine_classification_and_exact_path_limits_are_preserved() {
    for law in [
        FrictionLaw::Stribeck {
            static_mu: 0.4,
            kinetic_mu: 0.4,
            characteristic_speed: 0.2,
            viscous_per_speed: 0.3,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 0.4,
            kinetic_mu: 0.4,
            characteristic_speed: 0.2,
        },
        FrictionLaw::ShiftedPower {
            static_mu: 0.4,
            kinetic_mu: 0.4,
            characteristic_speed: 0.2,
            exponent: 0.7,
        },
    ] {
        let r = run(law, scalar(0, 0.7, 0.0, 2.0));
        assert!(r.exact_invariant);
        assert_eq!(r.trajectory.mode, HertzDrivenMode::Resting);
        assert!(r.stopping_time_bracket_s.is_some());
    }
    let patch = HertzFrictionPatch::new(N, A, affine(0.4, 0.3)).unwrap();
    let input = scalar(0, 0.7, 0.0, 2.0);
    let snapshot = input;
    let c = HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_evaluations: 1,
            ..controls().integration
        },
        ..controls()
    };
    assert_eq!(
        patch.advance_driven(&input, config(), grid(1), c, || false),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::EvaluationBudget
        ))
    );
    let c = HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_attempts: 1,
            ..controls().integration
        },
        ..controls()
    };
    assert_eq!(
        patch.advance_driven(&input, config(), grid(1), c, || false),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::AttemptBudget
        ))
    );
    let mut polls = 0;
    assert_eq!(
        patch.advance_driven(&input, config(), grid(1), controls(), || {
            polls += 1;
            polls > 3
        }),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::Cancelled
        ))
    );
    assert_eq!(
        patch.advance_driven(
            &input,
            HertzWrenchConfig {
                maximum_force_n: 0.01,
                ..config()
            },
            grid(1),
            controls(),
            || false
        ),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::EffortBudget
        ))
    );
    assert_eq!(input, snapshot);
    // One representable duration before the computed linear stop rounds q to
    // zero. It must refuse rather than publish a fictitious early arrest.
    let patch = HertzFrictionPatch::new(N, A, affine(0.35, 0.0)).unwrap();
    let stop = 0.3_f64 / 0.7;
    let before = f64::from_bits(stop.to_bits() - 1);
    let input = make_input(
        [0.3, 0.0, 0.0],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0; 3],
        before,
        HertzDrivenMode::Sliding,
    );
    assert!(matches!(
        patch.advance_driven(&input, config(), grid(1), controls(), || false),
        Err(HertzRateDrivenError::UnresolvedDrivenEvent { .. })
    ));
}

#[test]
fn long_near_collinear_phase_cannot_hide_transverse_motion_in_final_cancellation() {
    let epsilon = 128.0 * f64::EPSILON;
    let matrix = [[1.0, 0.0, epsilon], [0.0, 1.0, 0.0], [epsilon, 0.0, 1.0]];
    let (initial, duration, drive) = (2e6, 4e6, 0.4);
    let (dry, viscous) = (N * 0.4, N * (-5e-8));
    let (end, _, _) = moments(initial, drive - dry, viscous, duration);
    let contact_impulse = end - initial - drive * duration;
    let transverse_drive = -epsilon * contact_impulse / duration;
    // Both relative collinearity tests pass, and the unperturbed scalar
    // prediction has exactly cancelling final transverse contributions.
    assert!(epsilon < 256.0 * f64::EPSILON);
    assert!(transverse_drive.abs() < 256.0 * f64::EPSILON * drive);
    close(
        transverse_drive * duration + epsilon * contact_impulse,
        0.0,
        1e-20,
    );
    let halfway = moments(initial, drive - dry, viscous, 0.5 * duration).0;
    let half_impulse = halfway - initial - drive * 0.5 * duration;
    let interior = transverse_drive * 0.5 * duration + epsilon * half_impulse;
    assert!(interior.abs() > 1e-9);
    let input = make_input(
        [initial, 0.0, 0.0],
        matrix,
        [drive, 0.0, transverse_drive],
        duration,
        HertzDrivenMode::Sliding,
    );
    let patch = HertzFrictionPatch::new(N, A, affine(0.4, -5e-8)).unwrap();
    let controls = HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_step_s: duration,
            maximum_attempts: 1,
            scaled_rate_tolerance_m_s: 1e-10,
            ..controls().integration
        },
        ..controls()
    };
    // The analytic path must decline before consuming an evaluation budget.
    // One bounded full-matrix RK attempt cannot resolve this long horizon.
    assert_eq!(
        patch.advance_driven(&input, config(), grid(1), controls, || false),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::AttemptBudget
        ))
    );
}
