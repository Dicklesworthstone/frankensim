//! G1/G3: physical source-clock derivatives, including transient phase history.

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

// Independent RK4 of [q,v,q_Omega,v_Omega,q_phi,v_phi]. The physical initial
// state, natural frequency, damping, duration and initial phase stay fixed.
fn reference(omega: f64, c: f64, forcing: f64, duration: f64, fc: f64, fs: f64) -> [f64; 6] {
    let rhs = |t: f64, x: [f64; 6]| {
        let (sine, cosine) = (forcing * t).sin_cos();
        let phase_force = fs * cosine - fc * sine;
        [
            x[1],
            fc * cosine + fs * sine - omega * omega * x[0] - c * x[1],
            x[3],
            t * phase_force - omega * omega * x[2] - c * x[3],
            x[5],
            phase_force - omega * omega * x[4] - c * x[5],
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
fn source_clock_and_phase_match_independent_transient_variational_dynamics() {
    for (omega, c, forcing, h) in [
        (0.0, 0.0, 0.0, 0.5),
        (0.0, 3.0, 2.0, 0.7),
        (2.0, 0.0, 2.0, 1.7),
        (2.0, 0.2, 3.1, 0.7),
        (2.0, 4.0, 2.0, 0.7),
        (2.0, 3.999999, 2.0, 0.7),
        (2.0, 4.000001, 2.0, 0.7),
        (2.0, 20.0, 1.1, 0.7),
        (30.0, 1e-3, 30.0, 0.4),
        (2.0, 0.2, 0.001, 0.3),
        (2.0, 0.2, 9.0, 3.0),
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
            p.displacement_by_forcing_angular_frequency,
            p.velocity_by_forcing_angular_frequency,
            p.displacement_by_initial_phase,
            p.velocity_by_initial_phase,
        ];
        for (a, b) in actual
            .into_iter()
            .zip(reference(omega, c, forcing, h, 0.7, -0.4))
        {
            close(a, b, 2e-11);
        }
        // The clock changes the source, not the homogeneous state transition.
        let other = prepared
            .advance(state(-4.0, 9.0), 0.7, -0.4)
            .unwrap()
            .parameters;
        assert_eq!(
            p.displacement_by_forcing_angular_frequency,
            other.displacement_by_forcing_angular_frequency
        );
        assert_eq!(p.velocity_by_initial_phase, other.velocity_by_initial_phase);
    }
}

#[test]
fn exact_dc_retains_the_latent_sine_clock_and_phase_derivatives() {
    let h = 0.7_f64;
    let (fc, fs) = (0.6, -0.4);
    let free = HarmonicForceTransition::try_new_with_damping_rate(0.0, 0.0, 0.0, h)
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
    let p = free.advance(state(0.3, -0.2), fc, fs).unwrap().parameters;
    close(
        p.displacement_by_forcing_angular_frequency,
        fs * h.powi(3) / 6.0,
        2e-15,
    );
    close(
        p.velocity_by_forcing_angular_frequency,
        fs * h.powi(2) / 2.0,
        2e-15,
    );
    close(p.displacement_by_initial_phase, fs * h.powi(2) / 2.0, 2e-15);
    close(p.velocity_by_initial_phase, fs * h, 2e-15);
    let cosine = free.advance(state(0.3, -0.2), fc, 0.0).unwrap().parameters;
    assert_eq!(cosine.displacement_by_forcing_angular_frequency, 0.0);
    assert_eq!(cosine.velocity_by_initial_phase, 0.0);
    for (omega, c) in [(0.0, 3.0), (2.0, 4.0), (2.0, 0.2)] {
        let p = HarmonicForceTransition::try_new_with_damping_rate(omega, c, 0.0, h)
            .unwrap()
            .prepare_sensitivities()
            .unwrap()
            .advance(state(0.3, -0.2), fc, fs)
            .unwrap()
            .parameters;
        let ramp = LinearForceTransition::try_new_with_damping_rate(omega, c, h).unwrap();
        let frequency = ramp.advance(state(0.0, 0.0), 0.0, h * fs).unwrap();
        let phase = ramp.advance(state(0.0, 0.0), fs, fs).unwrap();
        close(
            p.displacement_by_forcing_angular_frequency,
            frequency.displacement_m_sqrt_kg,
            3e-14,
        );
        close(
            p.velocity_by_forcing_angular_frequency,
            frequency.velocity_m_sqrt_kg_per_s,
            3e-14,
        );
        close(
            p.displacement_by_initial_phase,
            phase.displacement_m_sqrt_kg,
            3e-14,
        );
        close(
            p.velocity_by_initial_phase,
            phase.velocity_m_sqrt_kg_per_s,
            3e-14,
        );
    }
}

#[test]
fn exact_resonance_keeps_the_driving_frequency_secular_derivative() {
    let omega = 2.0_f64;
    let h = 1.7_f64;
    let (sine, cosine) = (omega * h).sin_cos();
    let p = HarmonicForceTransition::try_new_with_damping_rate(omega, 0.0, omega, h)
        .unwrap()
        .prepare_sensitivities()
        .unwrap()
        .advance(state(0.0, 0.0), 1.0, 0.0)
        .unwrap()
        .parameters;
    // Removable limit of the original finite-time forced response at Omega=omega,
    // with the oscillator's omega held fixed. No steady-state particular exists.
    close(
        p.displacement_by_forcing_angular_frequency,
        h * h * cosine / (4.0 * omega) - h * sine / (4.0 * omega * omega),
        3e-14,
    );
    close(
        p.velocity_by_forcing_angular_frequency,
        h * cosine / (4.0 * omega) - h * h * sine / 4.0 - sine / (4.0 * omega * omega),
        3e-14,
    );
    close(
        p.displacement_by_initial_phase,
        h * cosine / (2.0 * omega) - sine / (2.0 * omega * omega),
        3e-14,
    );
    close(p.velocity_by_initial_phase, -h * sine / 2.0, 3e-14);
}

#[test]
fn reverse_trajectory_includes_source_phase_clock_dependence_at_every_interval() {
    let intervals = [0.03, 0.2, 0.07, 0.19, 0.11, 0.13];
    let duration: f64 = intervals.iter().sum();
    let (omega, c, forcing) = (4.0, 0.6, 6.0);
    let initial = state(0.3, -0.2);
    let initial_source = (0.7, -0.4);
    let mut x = initial;
    let mut source = initial_source;
    let mut history = Vec::new();
    for h in intervals {
        let motion =
            HarmonicForceTransition::try_new_with_damping_rate(omega, c, forcing, h).unwrap();
        let prepared = motion.prepare_sensitivities().unwrap();
        history.push((motion, prepared, x, source));
        x = motion.advance(x, source.0, source.1).unwrap();
        source = motion.advance_force_phase(source.0, source.1).unwrap();
    }
    // Include a terminal source observation as well as the mechanical endpoint.
    let endpoint = (0.2, -0.8);
    let endpoint_source = (0.3, 0.7);
    let mut state_bar = endpoint;
    let mut source_bar = endpoint_source;
    let (mut frequency_bar, mut k_bar, mut c_bar) = (0.0, 0.0, 0.0);
    for (motion, prepared, x, (fc, fs)) in history.into_iter().rev() {
        let body = prepared
            .pullback(x, fc, fs, state_bar.0, state_bar.1)
            .unwrap();
        let clock = motion
            .pullback_force_phase(fc, fs, source_bar.0, source_bar.1)
            .unwrap();
        frequency_bar += body.forcing_angular_frequency + clock.forcing_angular_frequency;
        k_bar += body.stiffness;
        c_bar += body.damping_rate;
        state_bar = (body.initial_displacement, body.initial_velocity);
        source_bar = (
            body.cosine_force + clock.cosine_force,
            body.sine_force + clock.sine_force,
        );
    }
    // Apply the initial phase chain rule ONCE to the original source covector;
    // also summing body.initial_phase would count the same dependence twice.
    let phase_bar = initial_source.1 * source_bar.0 - initial_source.0 * source_bar.1;
    let complete =
        HarmonicForceTransition::try_new_with_damping_rate(omega, c, forcing, duration).unwrap();
    let body = complete
        .prepare_sensitivities()
        .unwrap()
        .pullback(
            initial,
            initial_source.0,
            initial_source.1,
            endpoint.0,
            endpoint.1,
        )
        .unwrap();
    let clock = complete
        .pullback_force_phase(
            initial_source.0,
            initial_source.1,
            endpoint_source.0,
            endpoint_source.1,
        )
        .unwrap();
    close(
        frequency_bar,
        body.forcing_angular_frequency + clock.forcing_angular_frequency,
        3e-13,
    );
    close(k_bar, body.stiffness, 3e-13);
    close(c_bar, body.damping_rate, 3e-13);
    close(state_bar.0, body.initial_displacement, 3e-13);
    close(state_bar.1, body.initial_velocity, 3e-13);
    close(source_bar.0, body.cosine_force + clock.cosine_force, 3e-13);
    close(source_bar.1, body.sine_force + clock.sine_force, 3e-13);
    let whole = complete
        .advance(initial, initial_source.0, initial_source.1)
        .unwrap();
    close(
        x.displacement_m_sqrt_kg,
        whole.displacement_m_sqrt_kg,
        3e-13,
    );
    close(
        x.velocity_m_sqrt_kg_per_s,
        whole.velocity_m_sqrt_kg_per_s,
        3e-13,
    );
    let oracle = reference(
        omega,
        c,
        forcing,
        duration,
        initial_source.0,
        initial_source.1,
    );
    let source_phase = endpoint_source.0 * source.1 - endpoint_source.1 * source.0;
    close(
        frequency_bar,
        endpoint.0 * oracle[2] + endpoint.1 * oracle[3] + duration * source_phase,
        2e-11,
    );
    close(
        phase_bar,
        endpoint.0 * oracle[4] + endpoint.1 * oracle[5] + source_phase,
        2e-11,
    );
}

#[test]
fn phase_pullback_transposes_both_quadratures_and_refusals_leave_reuse_unchanged() {
    for forcing in [0.0, 0.1, 3.0, 17.0] {
        let motion =
            HarmonicForceTransition::try_new_with_damping_rate(3.0, 0.2, forcing, 0.7).unwrap();
        let prepared = motion.prepare_sensitivities().unwrap();
        let expected = prepared
            .pullback(state(0.3, -0.2), 0.7, -0.4, 0.2, 0.8)
            .unwrap();
        close(
            expected.initial_phase,
            -0.4 * expected.cosine_force - 0.7 * expected.sine_force,
            3e-14,
        );
        let phase = motion.pullback_force_phase(0.7, -0.4, 0.2, 0.8).unwrap();
        let direction = motion.advance_force_phase(0.3, 0.9).unwrap();
        close(
            0.2 * direction.0 + 0.8 * direction.1,
            0.3 * phase.cosine_force + 0.9 * phase.sine_force,
            3e-14,
        );
        assert!(
            motion
                .pullback_force_phase(f64::NAN, 0.0, 0.2, 0.8)
                .is_err()
        );
        assert!(
            motion
                .pullback_force_phase(0.7, -0.4, f64::INFINITY, 0.8)
                .is_err()
        );
        assert!(prepared.advance(state(0.0, 0.0), f64::NAN, 0.0).is_err());
        assert_eq!(
            expected,
            prepared
                .pullback(state(0.3, -0.2), 0.7, -0.4, 0.2, 0.8)
                .unwrap()
        );
    }
}
