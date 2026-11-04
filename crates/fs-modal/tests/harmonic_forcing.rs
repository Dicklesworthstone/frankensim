//! G0/G1/G3 continuous harmonic loads, resonance and Fourier superposition.

use fs_math::c64::C64;
use fs_modal::time::{
    HarmonicForceTransition, ModalAcousticMode, ModalAcousticState, advance_exact_zoh,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 8.0e-12 * scale.max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn g1_free_viscous_drag_under_a_tone_matches_first_order_velocity_solution() {
    let (rate, drive, h, fc, fs) = (2.5_f64, 3.0, 0.7, 1.2, -0.4);
    let initial = state(-0.2, 0.3);
    // Solve v'+c*v=Fc*cos(Omega*t)+Fs*sin(Omega*t) directly, then integrate
    // v=a*cos+b*sin+d*exp(-c*t) for displacement. No oscillator formula.
    let denominator = rate * rate + drive * drive;
    let a = (rate * fc - drive * fs) / denominator;
    let b = (drive * fc + rate * fs) / denominator;
    let d = initial.velocity_m_sqrt_kg_per_s - a;
    let (sine, cosine) = (drive * h).sin_cos();
    let expected_v = a * cosine + b * sine + d * (-rate * h).exp();
    let expected_q = initial.displacement_m_sqrt_kg + a * sine / drive
        + b * (1.0 - cosine) / drive - d * (-rate * h).exp_m1() / rate;
    let actual = HarmonicForceTransition::try_new_with_damping_rate(0.0, rate, drive, h)
        .unwrap().advance(initial, fc, fs).unwrap();
    close(actual.displacement_m_sqrt_kg, expected_q, expected_q.abs());
    close(actual.velocity_m_sqrt_kg_per_s, expected_v, expected_v.abs());
    let half = HarmonicForceTransition::try_new_with_damping_rate(0.0, rate, drive, h * 0.5).unwrap();
    let first = half.advance(initial, fc, fs).unwrap();
    let (s, c) = (drive * h * 0.5).sin_cos();
    let second = half.advance(first, fc * c + fs * s, fs * c - fc * s).unwrap();
    close(second.displacement_m_sqrt_kg, expected_q, expected_q.abs());
    close(second.velocity_m_sqrt_kg_per_s, expected_v, expected_v.abs());
}

#[test]
fn g1_viscous_free_dc_tone_reaches_terminal_velocity() {
    for (rate, h) in [(2.0_f64, 0.3), (1.0e6, 0.4)] {
        let force = 1.3;
        let initial = state(-0.2, 0.3);
        let terminal = force / rate;
        let transient = initial.velocity_m_sqrt_kg_per_s - terminal;
        let expected_v = terminal + transient * (-rate * h).exp();
        let expected_q = initial.displacement_m_sqrt_kg + terminal * h
            - transient * (-rate * h).exp_m1() / rate;
        let actual = HarmonicForceTransition::try_new_with_damping_rate(0.0, rate, 0.0, h)
            .unwrap().advance(initial, force, 0.0).unwrap();
        close(actual.displacement_m_sqrt_kg, expected_q, expected_q.abs());
        close(actual.velocity_m_sqrt_kg_per_s, expected_v, expected_v.abs());
    }
}

#[test]
fn g0_explicit_harmonic_damping_rate_refuses_invalid_parameters() {
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(HarmonicForceTransition::try_new_with_damping_rate(0.0, rate, 3.0, 0.1).is_err());
    }
    assert!(HarmonicForceTransition::try_new_with_damping_rate(0.0, f64::MAX, 3.0, 2.0).is_err());
    assert!(HarmonicForceTransition::try_new(0.0, 1.0, 3.0, 0.1).is_err());
}

#[test]
fn g1_matches_independent_high_precision_dimensional_exponentials() {
    // 115-digit Decimal reference: exponentiate h*G for
    // G=[[0,1,0,0],[-w^2,-2*zeta*w,1,0],[0,0,0,Omega],[0,0,-Omega,0]]
    // and apply it to [q0,v0,Fc,Fs]. Every input is its exact f64 value.
    // The independent full 4x4 exponential used an induced infinity norm
    // below 1/8 and a 1e-108 Taylor tail before squaring back. It neither
    // uses the scalar transition nor decomposes cosine into a DC remainder.
    let cases = [
        (
            9.0,
            0.08,
            2.4,
            0.3,
            1.3,
            -0.4,
            0.2,
            -0.1,
            -0.120_581_881_328_301_67,
            -0.546_453_151_939_771_7,
        ),
        (
            5.0,
            0.999_999_995,
            5.0,
            0.4,
            3.0,
            -1.2,
            0.03,
            -0.2,
            0.019_939_299_975_225_363,
            -0.174_412_217_939_457_85,
        ),
        (
            5.0,
            1.0,
            5.0,
            0.4,
            3.0,
            -1.2,
            0.03,
            -0.2,
            0.019_939_299_983_345_97,
            -0.174_412_217_559_114_96,
        ),
        (
            5.0,
            1.000_000_005,
            5.0,
            0.4,
            3.0,
            -1.2,
            0.03,
            -0.2,
            0.019_939_299_991_466_58,
            -0.174_412_217_178_772_06,
        ),
        (
            100.0,
            1.0e12,
            3.0,
            0.1,
            2.0e14,
            -1.0e14,
            0.0,
            0.0,
            0.091_062_817_074_474_47,
            0.807_576_385_790_394_7,
        ),
        (
            1.0e-6,
            0.05,
            30_000.0,
            1.0 / 48_000.0,
            100.0,
            300.0,
            0.0,
            0.0,
            3.430_500_685_257_776e-8,
            0.003_840_693_048_079_264,
        ),
        (
            1.0,
            0.0,
            1.0e-8,
            0.1,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            -1.665_833_531_718_477_2e-20,
        ),
        (
            1.0e-160, 0.0, 1.0e-160, 1.0e-4, 1.0, 3.0, 0.0, 0.0, 5.0e-9, 0.0001,
        ),
    ];
    for (w, zeta, drive, h, fc, fs, q, v, expected_q, expected_v) in cases {
        let step = HarmonicForceTransition::try_new(w, zeta, drive, h).unwrap();
        let actual = step.advance(state(q, v), fc, fs).unwrap();
        close(actual.displacement_m_sqrt_kg, expected_q, expected_q.abs());
        close(
            actual.velocity_m_sqrt_kg_per_s,
            expected_v,
            expected_v.abs(),
        );
    }
}

#[test]
fn g1_undamped_resonance_retains_secular_growth_and_both_quadratures() {
    let (fc, fs) = (1.3, -0.7);
    for (w, h) in [
        (11.0, 0.003),
        (11.0, 0.3),
        (11.0, 10.0),
        (core::f64::consts::TAU * 440.0, 1.0),
    ] {
        let angle = w * h;
        let (s, c) = angle.sin_cos();
        // Independent resonant solution with q(0)=v(0)=0. These t*sin and
        // t*cos terms must survive exactly coincident force/natural poles.
        let q = fc * h * s / (2.0 * w) + fs * (s - angle * c) / (2.0 * w * w);
        let v = fc * (s + angle * c) / (2.0 * w) + fs * h * s / 2.0;
        let actual = HarmonicForceTransition::try_new(w, 0.0, w, h)
            .unwrap()
            .forced_increment(fc, fs)
            .unwrap();
        close(
            actual.displacement_m_sqrt_kg,
            q,
            (fc.abs() + fs.abs()) * h / w,
        );
        close(
            actual.velocity_m_sqrt_kg_per_s,
            v,
            (fc.abs() + fs.abs()) * h,
        );
    }
    let exact = HarmonicForceTransition::try_new(11.0, 0.0, 11.0, 0.3)
        .unwrap()
        .forced_increment(fc, fs)
        .unwrap();
    for relative_offset in [-1.0e-12, 1.0e-12] {
        let nearby =
            HarmonicForceTransition::try_new(11.0, 0.0, 11.0 * (1.0 + relative_offset), 0.3)
                .unwrap()
                .forced_increment(fc, fs)
                .unwrap();
        close(
            nearby.displacement_m_sqrt_kg,
            exact.displacement_m_sqrt_kg,
            0.1,
        );
        close(
            nearby.velocity_m_sqrt_kg_per_s,
            exact.velocity_m_sqrt_kg_per_s,
            1.0,
        );
    }
}

#[test]
fn g1_high_harmonics_follow_continuous_complex_response_without_hold_error() {
    let (w, zeta, h, fc, fs) = (8.0, 0.07, 0.31, 1.4, -0.6);
    for ratio in [0.5, 1.0, 2.0, 7.0, 31.0] {
        let drive = ratio * w;
        let phasor = C64::new(fc, fs) / C64::new(w * w - drive * drive, -2.0 * zeta * w * drive);
        let initial = state(phasor.re, drive * phasor.im);
        let transition = HarmonicForceTransition::try_new(w, zeta, drive, h).unwrap();
        let actual = transition.advance(initial, fc, fs).unwrap();
        let (s, c) = (drive * h).sin_cos();
        close(
            actual.displacement_m_sqrt_kg,
            phasor.re * c + phasor.im * s,
            phasor.abs(),
        );
        close(
            actual.velocity_m_sqrt_kg_per_s,
            drive * (-phasor.re * s + phasor.im * c),
            drive * phasor.abs(),
        );
        let shifted = transition.advance_force_phase(fc, fs).unwrap();
        close(shifted.0, fc * c + fs * s, fc.abs() + fs.abs());
        close(shifted.1, fs * c - fc * s, fc.abs() + fs.abs());
    }
}

#[test]
fn g1_dc_force_recovers_existing_zoh_in_every_damping_regime() {
    for (w, zeta) in [
        (0.0, 0.0),
        (7.0, 0.0),
        (7.0, 0.02),
        (7.0, 1.0),
        (7.0, 1.000_000_01),
        (7.0, 30.0),
    ] {
        for h in [1.0e-5, 0.2, 5.0] {
            let initial = state(0.03, -0.7);
            let mode = ModalAcousticMode {
                angular_frequency_rad_s: w,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            };
            let expected = advance_exact_zoh(mode, initial, 2.0, h);
            let actual = HarmonicForceTransition::try_new(w, zeta, 0.0, h)
                .unwrap()
                .advance(initial, 2.0, 13.0)
                .unwrap();
            close(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                0.03 + expected.displacement_m_sqrt_kg.abs(),
            );
            close(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
                0.7 + expected.velocity_m_sqrt_kg_per_s.abs(),
            );
        }
    }
}

#[test]
fn g1_free_coordinate_has_the_exact_sine_and_cosine_impulse() {
    let (drive, h, fc, fs) = (3.0, 0.7, 1.2, -0.4);
    let initial = state(0.1, -0.3);
    let angle: f64 = drive * h;
    let (s, c) = angle.sin_cos();
    let actual = HarmonicForceTransition::try_new(0.0, 0.0, drive, h)
        .unwrap()
        .advance(initial, fc, fs)
        .unwrap();
    let q = 0.1 - 0.3 * h + (fc * (1.0 - c) + fs * (angle - s)) / (drive * drive);
    let v = -0.3 + (fc * s + fs * (1.0 - c)) / drive;
    close(actual.displacement_m_sqrt_kg, q, 1.0);
    close(actual.velocity_m_sqrt_kg_per_s, v, 1.0);
}

fn multi_step(
    transitions: &[HarmonicForceTransition],
    initial: ModalAcousticState,
    forces: &mut [(f64, f64)],
) -> ModalAcousticState {
    let mut output = transitions[0].advance(initial, 0.0, 0.0).unwrap();
    for (transition, (fc, fs)) in transitions.iter().zip(forces) {
        let part = transition.forced_increment(*fc, *fs).unwrap();
        output.displacement_m_sqrt_kg += part.displacement_m_sqrt_kg;
        output.velocity_m_sqrt_kg_per_s += part.velocity_m_sqrt_kg_per_s;
        (*fc, *fs) = transition.advance_force_phase(*fc, *fs).unwrap();
    }
    output
}

#[test]
fn g3_fourier_superposition_subdivision_and_checkpoint_replay() {
    for zeta in [0.0, 0.04, 1.0, 3.0, 1.0e8] {
        let (w, h) = (8.0, 0.38);
        let drives = [0.0, 4.0, 24.0];
        let forces = [(1.0, 0.0), (0.7, 0.1), (-0.2, 0.5)];
        let whole: Vec<_> = drives
            .iter()
            .map(|&drive| HarmonicForceTransition::try_new(w, zeta, drive, h).unwrap())
            .collect();
        let parts: Vec<_> = drives
            .iter()
            .map(|&drive| HarmonicForceTransition::try_new(w, zeta, drive, h / 17.0).unwrap())
            .collect();
        let initial = state(0.12, -0.3);
        let mut whole_forces = forces;
        let expected = multi_step(&whole, initial, &mut whole_forces);
        let mut phased = forces;
        let mut actual = initial;
        for _ in 0..8 {
            actual = multi_step(&parts, actual, &mut phased);
        }
        let checkpoint = (actual, phased);
        for _ in 8..17 {
            actual = multi_step(&parts, actual, &mut phased);
        }
        close(
            actual.displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            0.3,
        );
        close(
            actual.velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            1.0,
        );
        let (mut resumed, mut resumed_forces) = checkpoint;
        for _ in 8..17 {
            resumed = multi_step(&parts, resumed, &mut resumed_forces);
        }
        assert_eq!(
            resumed.displacement_m_sqrt_kg.to_bits(),
            actual.displacement_m_sqrt_kg.to_bits()
        );
        assert_eq!(
            resumed.velocity_m_sqrt_kg_per_s.to_bits(),
            actual.velocity_m_sqrt_kg_per_s.to_bits()
        );
        assert_eq!(resumed_forces, phased);
    }
}

#[test]
fn g0_invalid_inputs_and_unrepresentable_trials_refuse_without_hidden_state() {
    for (w, zeta, drive, h) in [
        (-1.0, 0.0, 1.0, 0.1),
        (1.0, -0.1, 1.0, 0.1),
        (1.0, 0.0, -1.0, 0.1),
        (1.0, 0.0, 1.0, 0.0),
        (0.0, 0.1, 1.0, 0.1),
        (f64::NAN, 0.0, 1.0, 0.1),
        (1.0, 0.0, f64::INFINITY, 0.1),
        (1.0e300, 0.0, 1.0, 1.0e300),
    ] {
        assert!(HarmonicForceTransition::try_new(w, zeta, drive, h).is_err());
    }
    let transition = HarmonicForceTransition::try_new(4.0, 0.1, 4.0, 0.2).unwrap();
    let initial = state(0.2, -0.3);
    let first = transition.advance(initial, 1.0, -0.4).unwrap();
    assert!(transition.advance(state(f64::NAN, 0.0), 1.0, -0.4).is_err());
    assert!(transition.forced_increment(f64::INFINITY, 1.0).is_err());
    assert!(transition.advance_force_phase(f64::NAN, 1.0).is_err());
    assert_eq!(transition.advance(initial, 1.0, -0.4).unwrap(), first);
    let long = HarmonicForceTransition::try_new(0.0, 0.0, 0.0, 1.0e200).unwrap();
    assert!(long.forced_increment(1.0e200, 0.0).is_err());
}
