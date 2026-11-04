//! Physical scalar harmonic work/heat with explicit viscous drag.

use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAcousticWorkspace, PreparedHarmonicForce, PreparedModalHarmonic,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-11 * scale.max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn free_viscous_tone_has_analytic_motion_work_and_positive_heat() {
    let (drive, h, fc, fs) = (3.0_f64, 0.7, 1.2, -0.4);
    let initial = state(-0.2, 0.3);
    for rate in [0.01_f64, 2.5, 1.0e6] {
        // Solve v'+c*v=Fc*cos(Omega*t)+Fs*sin(Omega*t) independently:
        // v=a*cos(Omega*t)+b*sin(Omega*t)+d*exp(-c*t).
        // Integrate F*v and c*v^2 in elementary closed form, retaining the
        // cosine/sine and homogeneous/forced cross terms in both powers.
        let denominator = rate * rate + drive * drive;
        let a = (rate * fc - drive * fs) / denominator;
        let b = (drive * fc + rate * fs) / denominator;
        let d = initial.velocity_m_sqrt_kg_per_s - a;
        let (sine, cosine) = (drive * h).sin_cos();
        let decay = (-rate * h).exp();
        let cc = 0.5 * h + (2.0 * drive * h).sin() / (4.0 * drive);
        let ss = 0.5 * h - (2.0 * drive * h).sin() / (4.0 * drive);
        let cs = sine * sine / (2.0 * drive);
        let ec = (rate + decay * (-rate * cosine + drive * sine)) / denominator;
        let es = (drive + decay * (-rate * sine - drive * cosine)) / denominator;
        let ee = -(-2.0 * rate * h).exp_m1() / (2.0 * rate);
        let expected_work =
            fc * a * cc + (fc * b + fs * a) * cs + fs * b * ss + d * (fc * ec + fs * es);
        let expected_heat = rate
            * (a * a * cc
                + 2.0 * a * b * cs
                + b * b * ss
                + 2.0 * d * (a * ec + b * es)
                + d * d * ee);
        let expected_q =
            initial.displacement_m_sqrt_kg + a * sine / drive + b * (1.0 - cosine) / drive
                - d * (-rate * h).exp_m1() / rate;
        let expected_v = a * cosine + b * sine + d * decay;
        let prepared =
            PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, drive, h).unwrap();
        let frame = prepared.advance_with_energy(initial, fc, fs).unwrap();
        near(
            frame.state.displacement_m_sqrt_kg,
            expected_q,
            expected_q.abs(),
        );
        near(
            frame.state.velocity_m_sqrt_kg_per_s,
            expected_v,
            expected_v.abs(),
        );
        near(frame.input_work_j, expected_work, expected_work.abs());
        near(frame.viscous_dissipation_j, expected_heat, expected_heat);
        assert!(frame.viscous_dissipation_j > 0.0);
        near(
            frame.input_work_j - frame.viscous_dissipation_j,
            0.5 * (expected_v.powi(2) - initial.velocity_m_sqrt_kg_per_s.powi(2)),
            expected_heat.max(expected_work.abs()),
        );
        assert_eq!(
            frame.state,
            prepared.motion().advance(initial, fc, fs).unwrap()
        );
    }
}

#[test]
fn dc_drag_returns_negative_source_work_and_independent_heat() {
    let (rate, h, force) = (4.0_f64, 0.2, -0.3);
    // Absolute-position subtraction loses the whole displacement increment.
    // Work must still retain the real force/velocity integral.
    let initial = state(1.0e20, 0.7);
    let terminal = force / rate;
    let d = initial.velocity_m_sqrt_kg_per_s - terminal;
    let i0 = -(-rate * h).exp_m1() / rate;
    let i2 = -(-2.0 * rate * h).exp_m1() / (2.0 * rate);
    let work = force * (terminal * h + d * i0);
    let heat = rate * (terminal * terminal * h + 2.0 * terminal * d * i0 + d * d * i2);
    let frame = PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, 0.0, h)
        .unwrap()
        .advance_with_energy(initial, force, 1.0e6)
        .unwrap();
    // At DC the sine quadrature vanishes, even for a nonzero supplied value.
    assert_eq!(
        frame.state.displacement_m_sqrt_kg,
        initial.displacement_m_sqrt_kg
    );
    assert!(frame.input_work_j < 0.0);
    near(frame.input_work_j, work, work.abs());
    near(frame.viscous_dissipation_j, heat, heat);
}

#[test]
fn sub_roundoff_drag_heat_survives_and_zero_drag_is_exactly_lossless() {
    let (rate, h) = (1.0e-18_f64, 0.7);
    let initial = state(0.0, 0.7);
    let frame = PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, 3.0, h)
        .unwrap()
        .advance_with_energy(initial, 0.0, 0.0)
        .unwrap();
    let expected = -0.5 * initial.velocity_m_sqrt_kg_per_s.powi(2) * (-2.0 * rate * h).exp_m1();
    assert_eq!(
        frame.state.velocity_m_sqrt_kg_per_s,
        initial.velocity_m_sqrt_kg_per_s
    );
    assert!(frame.viscous_dissipation_j > 0.0);
    near(frame.viscous_dissipation_j, expected, expected);
    assert_eq!(frame.input_work_j, 0.0);
    let lossless = PreparedHarmonicForce::try_new_with_damping_rate(0.0, 0.0, 3.0, h)
        .unwrap()
        .advance_with_energy(initial, 1.2, -0.4)
        .unwrap();
    assert_eq!(lossless.viscous_dissipation_j, 0.0);
    assert_eq!(lossless.dissipation_roundoff_tolerance_j, 0.0);
}

#[test]
fn force_phase_subdivision_preserves_work_heat_and_replay() {
    let (rate, drive, h, fc, fs) = (2.5_f64, 3.0, 0.7, 1.2, -0.4);
    let initial = state(-0.2, 0.3);
    let whole = PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, drive, h).unwrap();
    let half = PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, drive, 0.5 * h).unwrap();
    let one = whole.advance_with_energy(initial, fc, fs).unwrap();
    let first = half.advance_with_energy(initial, fc, fs).unwrap();
    let (next_fc, next_fs) = half.motion().advance_force_phase(fc, fs).unwrap();
    let second = half
        .advance_with_energy(first.state, next_fc, next_fs)
        .unwrap();
    near(
        second.state.displacement_m_sqrt_kg,
        one.state.displacement_m_sqrt_kg,
        one.state.displacement_m_sqrt_kg.abs(),
    );
    near(
        second.state.velocity_m_sqrt_kg_per_s,
        one.state.velocity_m_sqrt_kg_per_s,
        one.state.velocity_m_sqrt_kg_per_s.abs(),
    );
    near(
        first.input_work_j + second.input_work_j,
        one.input_work_j,
        one.input_work_j.abs(),
    );
    near(
        first.viscous_dissipation_j + second.viscous_dissipation_j,
        one.viscous_dissipation_j,
        one.viscous_dissipation_j,
    );
    assert_eq!(
        whole.clone().advance_with_energy(initial, fc, fs).unwrap(),
        one
    );
}

#[test]
fn ratio_rate_and_transactional_bank_share_identical_power_arithmetic() {
    let (omega, drive, h, fc, fs) = (4.0, 3.0, 0.25, 1.2, -0.4);
    let initial = state(0.1, -0.2);
    for ratio in [0.0, 0.125, 1.0, 32.0] {
        let scalar = PreparedHarmonicForce::try_new(omega, ratio, drive, h)
            .unwrap()
            .advance_with_energy(initial, fc, fs)
            .unwrap();
        let explicit =
            PreparedHarmonicForce::try_new_with_damping_rate(omega, 2.0 * ratio * omega, drive, h)
                .unwrap()
                .advance_with_energy(initial, fc, fs)
                .unwrap();
        assert_eq!(scalar, explicit);
        let mut model = ModalAcousticTimeModel::try_new(
            48_000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: ratio,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        model.restore_states(&[initial]).unwrap();
        let prepared = PreparedModalHarmonic::try_new(&model, drive, h).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&model);
        let frame = model
            .step_harmonic_into(&prepared, &[fc], &[fs], &mut workspace)
            .unwrap();
        assert_eq!(scalar.input_work_j, frame.input_work_j);
        assert_eq!(scalar.viscous_dissipation_j, frame.viscous_dissipation_j);
        assert_eq!(
            scalar.dissipation_roundoff_tolerance_j,
            frame.modal_viscous_losses[0].roundoff_tolerance_j
        );
        assert_eq!(scalar.state, model.states()[0]);
    }
}

#[test]
fn invalid_trials_refuse_without_mutating_preparation_or_state() {
    for rate in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, 3.0, 0.7).is_err());
    }
    assert!(PreparedHarmonicForce::try_new_with_damping_rate(0.0, f64::MAX, 3.0, 2.0).is_err());
    assert!(PreparedHarmonicForce::try_new(0.0, 1.0, 3.0, 0.7).is_err());
    let initial = state(-0.2, 0.3);
    let prepared = PreparedHarmonicForce::try_new_with_damping_rate(0.0, 2.5, 3.0, 0.7).unwrap();
    let before = prepared.advance_with_energy(initial, 1.2, -0.4).unwrap();
    for (trial, fc, fs) in [
        (state(f64::NAN, 0.3), 1.2, -0.4),
        (initial, f64::INFINITY, -0.4),
        (initial, 1.2, f64::NAN),
        (state(0.0, 1.0e200), 0.0, 0.0),
    ] {
        assert!(prepared.advance_with_energy(trial, fc, fs).is_err());
        assert_eq!(
            prepared.advance_with_energy(initial, 1.2, -0.4).unwrap(),
            before
        );
    }
}
