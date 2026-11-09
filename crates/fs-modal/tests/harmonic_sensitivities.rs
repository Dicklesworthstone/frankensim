//! Continuous-tone transient derivatives against independent variational ODEs.

use fs_modal::time::{HarmonicForceTransition, LinearForceTransition, ModalAcousticState};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

// Independent RK4 solves the original motion and physical variational ODEs,
// not the production entire series, convolution, or their algebraic derivatives.
// Coordinates are [q,v,q_k,v_k,q_c,v_c]; the prescribed force clock is fixed.
fn reference(omega: f64, c: f64, forcing: f64, duration: f64, fc: f64, fs: f64) -> [f64; 6] {
    let rhs = |t: f64, x: [f64; 6]| {
        let force = fc * (forcing * t).cos() + fs * (forcing * t).sin();
        [
            x[1],
            force - omega * omega * x[0] - c * x[1],
            x[3],
            -x[0] - omega * omega * x[2] - c * x[3],
            x[5],
            -x[1] - omega * omega * x[4] - c * x[5],
        ]
    };
    let add = |a: [f64; 6], b: [f64; 6], scale: f64| std::array::from_fn(|i| a[i] + scale * b[i]);
    let n = 16384;
    let h = duration / n as f64;
    let mut x = [0.3, -0.2, 0.0, 0.0, 0.0, 0.0];
    for j in 0..n {
        let t = j as f64 * h;
        let a = rhs(t, x);
        let b = rhs(t + h * 0.5, add(x, a, h * 0.5));
        let d = rhs(t + h * 0.5, add(x, b, h * 0.5));
        let e = rhs(t + h, add(x, d, h));
        x = std::array::from_fn(|i| x[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]));
    }
    x
}

#[test]
fn startup_and_driven_partials_match_independent_variational_dynamics() {
    for (omega, c, forcing, h) in [
        (0.0, 0.0, 3.0, 0.5),
        (0.0, 3.0, 2.0, 0.7),
        (2.0, 0.0, 2.0, 1.7),
        (2.0, 0.2, 3.1, 0.7),
        (2.0, 4.0, 2.0, 0.7),
        (2.0, 3.999999, 2.0, 0.7),
        (2.0, 4.000001, 2.0, 0.7),
        (2.0, 20.0, 1.1, 0.7),
        (30.0, 1e-3, 30.0, 0.4),
        (2.0, 0.2, 0.1, 0.3),
    ] {
        let motion =
            HarmonicForceTransition::try_new_with_damping_rate(omega, c, forcing, h).unwrap();
        let prepared = motion.prepare_sensitivities().unwrap();
        let frame = prepared.advance(state(0.3, -0.2), 0.7, -0.4).unwrap();
        assert_eq!(
            frame.state,
            motion.advance(state(0.3, -0.2), 0.7, -0.4).unwrap()
        );
        assert_eq!(
            frame,
            prepared.advance(state(0.3, -0.2), 0.7, -0.4).unwrap()
        );
        let p = frame.parameters;
        let actual = [
            frame.state.displacement_m_sqrt_kg,
            frame.state.velocity_m_sqrt_kg_per_s,
            p.displacement_by_stiffness,
            p.velocity_by_stiffness,
            p.displacement_by_damping_rate,
            p.velocity_by_damping_rate,
        ];
        for (actual, expected) in actual
            .into_iter()
            .zip(reference(omega, c, forcing, h, 0.7, -0.4))
        {
            close(actual, expected, 2e-11);
        }
    }
}

#[test]
fn exact_resonance_keeps_frequency_derivative_of_secular_growth() {
    let omega = 2.0_f64;
    let h = 1.7_f64;
    let phase = omega * h;
    let motion = HarmonicForceTransition::try_new_with_damping_rate(omega, 0.0, omega, h).unwrap();
    let frame = motion
        .prepare_sensitivities()
        .unwrap()
        .advance(state(0.0, 0.0), 1.0, 0.0)
        .unwrap();
    close(
        frame.state.displacement_m_sqrt_kg,
        h * phase.sin() / (2.0 * omega),
        2e-13,
    );
    // Differentiate (cos(Omega*h)-cos(sqrt(k)*h))/(k-Omega^2) at fixed Omega.
    // Its removable resonant limit contains the secular h^2 contribution.
    close(
        frame.parameters.displacement_by_stiffness,
        h * h * phase.cos() / (8.0 * omega * omega) - h * phase.sin() / (8.0 * omega.powi(3)),
        2e-13,
    );
    close(
        frame.parameters.velocity_by_stiffness,
        h * phase.cos() / (8.0 * omega * omega)
            - h * h * phase.sin() / (8.0 * omega)
            - phase.sin() / (8.0 * omega.powi(3)),
        2e-13,
    );
}

#[test]
fn dc_matches_held_force_and_static_balance_has_no_damping_sensitivity() {
    for (omega, c, h) in [
        (0.0, 0.0, 0.5),
        (0.0, 3.0, 0.7),
        (2.0, 4.0, 0.3),
        (1.0, 2.0, 40.0),
    ] {
        let harmonic = HarmonicForceTransition::try_new_with_damping_rate(omega, c, 0.0, h)
            .unwrap()
            .prepare_sensitivities()
            .unwrap();
        let held = LinearForceTransition::try_new_with_damping_rate(omega, c, h)
            .unwrap()
            .prepare_sensitivities()
            .unwrap();
        let a = harmonic
            .advance(state(0.3, -0.2), 0.7, -4.0)
            .unwrap()
            .parameters;
        let b = held.advance(state(0.3, -0.2), 0.7, 0.7).unwrap().parameters;
        for (a, b) in [
            (a.displacement_by_stiffness, b.displacement_by_stiffness),
            (a.velocity_by_stiffness, b.velocity_by_stiffness),
            (
                a.displacement_by_damping_rate,
                b.displacement_by_damping_rate,
            ),
            (a.velocity_by_damping_rate, b.velocity_by_damping_rate),
        ] {
            close(a, b, 3e-14);
        }
        assert_eq!(
            harmonic
                .pullback(state(0.3, -0.2), 0.7, -4.0, 0.2, 0.8)
                .unwrap()
                .sine_force,
            0.0
        );
    }
    let prepared = HarmonicForceTransition::try_new_with_damping_rate(1.0, 2.0, 0.0, 40.0)
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
    let equilibrium = prepared
        .advance(state(1.0, 0.0), 1.0, 0.0)
        .unwrap()
        .parameters;
    assert_eq!(equilibrium.displacement_by_damping_rate, 0.0);
    assert_eq!(equilibrium.velocity_by_damping_rate, 0.0);
    let driven = prepared
        .advance(state(0.0, 0.0), 1.0, 0.0)
        .unwrap()
        .parameters;
    let expected = -40.0_f64.powi(3) * (-40.0_f64).exp() / 6.0;
    assert!((driven.displacement_by_damping_rate / expected - 1.0).abs() < 2e-11);
}

#[test]
fn reverse_sweep_sums_all_transient_parameter_and_force_contributions() {
    let count = 9;
    let h = 0.13;
    let (omega, c, forcing) = (4.0, 0.6, 6.0);
    let motion = HarmonicForceTransition::try_new_with_damping_rate(omega, c, forcing, h).unwrap();
    let prepared = motion.prepare_sensitivities().unwrap();
    let mut history = Vec::new();
    let mut x = state(0.3, -0.2);
    let mut source = (0.7, -0.4);
    for _ in 0..count {
        history.push((x, source));
        x = motion.advance(x, source.0, source.1).unwrap();
        source = motion.advance_force_phase(source.0, source.1).unwrap();
    }
    let lambda = (0.2, -0.8);
    let mut cotangent = lambda;
    let (mut k_bar, mut c_bar, mut source_gain_bar) = (0.0, 0.0, 0.0);
    for (state, (fc, fs)) in history.into_iter().rev() {
        let pb = prepared
            .pullback(state, fc, fs, cotangent.0, cotangent.1)
            .unwrap();
        k_bar += pb.stiffness;
        c_bar += pb.damping_rate;
        source_gain_bar += pb.cosine_force * fc + pb.sine_force * fs;
        cotangent = (pb.initial_displacement, pb.initial_velocity);
    }
    let complete =
        HarmonicForceTransition::try_new_with_damping_rate(omega, c, forcing, h * count as f64)
            .unwrap();
    let direct = complete
        .prepare_sensitivities()
        .unwrap()
        .pullback(state(0.3, -0.2), 0.7, -0.4, lambda.0, lambda.1)
        .unwrap();
    close(k_bar, direct.stiffness, 2e-13);
    close(c_bar, direct.damping_rate, 2e-13);
    close(cotangent.0, direct.initial_displacement, 2e-13);
    close(cotangent.1, direct.initial_velocity, 2e-13);
    close(
        source_gain_bar,
        0.7 * direct.cosine_force - 0.4 * direct.sine_force,
        2e-13,
    );
    let oracle = reference(omega, c, forcing, h * count as f64, 0.7, -0.4);
    close(k_bar, lambda.0 * oracle[2] + lambda.1 * oracle[3], 2e-11);
    close(c_bar, lambda.0 * oracle[4] + lambda.1 * oracle[5], 2e-11);
}

#[test]
fn input_transpose_retains_both_phase_quadratures_and_refuses_invalid_inputs() {
    for forcing in [0.0, 0.1, 3.0, 17.0] {
        let motion =
            HarmonicForceTransition::try_new_with_damping_rate(3.0, 0.2, forcing, 0.7).unwrap();
        let prepared = motion.prepare_sensitivities().unwrap();
        let pb = prepared
            .pullback(state(0.3, -0.2), 0.7, -0.4, 0.2, 0.8)
            .unwrap();
        let direction = motion.advance(state(0.4, -0.7), 0.3, 0.9).unwrap();
        close(
            0.2 * direction.displacement_m_sqrt_kg + 0.8 * direction.velocity_m_sqrt_kg_per_s,
            0.4 * pb.initial_displacement - 0.7 * pb.initial_velocity
                + 0.3 * pb.cosine_force
                + 0.9 * pb.sine_force,
            2e-14,
        );
        assert!(prepared.advance(state(f64::NAN, 0.0), 0.0, 0.0).is_err());
        assert!(
            prepared
                .advance(state(0.0, 0.0), f64::INFINITY, 0.0)
                .is_err()
        );
        assert!(
            prepared
                .pullback(state(0.0, 0.0), 0.0, 0.0, f64::NAN, 0.0)
                .is_err()
        );
    }
}
