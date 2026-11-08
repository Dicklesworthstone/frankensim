use super::super::{CompactSourceRecordingBudget, CompactSourceWelchStream};
use super::*;
use fs_fft::welch::{WelchDetrend, WelchWindow};
use fs_math::det;

const DT: f64 = 1.0 / 8000.0;
const BIN: usize = 8;
const ORIGIN: [[f64; 3]; 1] = [[0.0; 3]];

fn acoustic() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 343.0).unwrap()
}
fn clock() -> SourceSampleClock {
    SourceSampleClock::new(0.25, DT, 0.0).unwrap()
}
fn config() -> WelchConfig {
    WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    }
}
fn budget() -> JointSourceRecordingBudget {
    JointSourceRecordingBudget {
        maximum_compact_sources: 2,
        maximum_quadrupole_sources: 2,
        welch: WelchBudget {
            maximum_channels: 20,
            maximum_samples_per_channel: 128,
            maximum_transforms: 64,
            maximum_workspace_bytes: 1 << 20,
        },
    }
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn phase(n: usize) -> (f64, f64) {
    det::sin_cos(core::f64::consts::TAU * BIN as f64 * n as f64 / 64.0)
}

#[test]
fn analytic_compact_stress_interference_cancels_pressure_but_preserves_velocity() {
    // Independent negative-time monopole + xx stress formula on the x axis:
    // a=2/r²-2ik/r-k², Q=-i*a/(rho*omega), G=exp(ikr)/(4pi*r).
    // p=0 but ux=(2k/r²+4i/r³)*G/(rho*omega). This is a generic
    // disjoint-mechanism fixture, not overlapping FW-H surface/volume geometry.
    let radius = 1.7;
    let omega = core::f64::consts::TAU * 1000.0;
    let k = omega / acoustic().sound_speed_m_s();
    let rho_omega = acoustic().density_kg_m3() * omega;
    let a = C64::new(2.0 / (radius * radius) - k * k, -2.0 * k / radius);
    let q = C64::new(0.0, -1.0 / rho_omega) * a;
    let (sin, cos) = det::sin_cos(k * radius);
    let green = C64::new(cos, sin).scale(1.0 / (2.0 * core::f64::consts::TAU * radius));
    let ux = (C64::new(
        2.0 * k / (radius * radius),
        4.0 / (radius * radius * radius),
    ) * green)
        .scale(1.0 / rho_omega);
    let independent_pressure_square = (a * green).norm_sq();
    assert!(independent_pressure_square > 200.0);
    let expected_velocity_square = 0.5 * ux.norm_sq();
    assert!(expected_velocity_square > 1e-9);
    let mut recording = JointSourceWelchStream::new(
        acoustic(),
        &ORIGIN,
        &ORIGIN,
        clock(),
        128,
        config(),
        budget(),
        || false,
    )
    .unwrap();
    for n in 0..128 {
        let (s, c) = phase(n);
        recording
            .push_sample(
                &ORIGIN,
                &ORIGIN,
                clock().time_at(n).unwrap(),
                &[[q.re * c + q.im * s, 0.0, 0.0, 0.0]],
                &[[c, 0.0, 0.0, 0.0, 0.0, 0.0]],
                || false,
            )
            .unwrap();
    }
    let recorded = recording.finish().unwrap();
    let (compact, stress) = recorded
        .prepare_observers([radius, 0.0, 0.0], BIN, || false)
        .unwrap();
    let actual = recorded
        .evaluate_bin(
            &compact,
            &stress,
            BIN,
            CoherenceBudget {
                max_modes: 3,
                max_source_mode_pairs: 6,
            },
            || false,
        )
        .unwrap();
    assert!(actual.mean_square_pressure_pa2 < 1e-22 * independent_pressure_square);
    close(
        actual.velocity_covariance_m2_s2[0][0].re,
        expected_velocity_square,
        expected_velocity_square,
    );
    assert_eq!(actual.velocity_covariance_m2_s2[1][1].re, 0.0);
    assert_eq!(actual.velocity_covariance_m2_s2[2][2].re, 0.0);
    let flux_scale = det::sqrt(independent_pressure_square * expected_velocity_square);
    for axis in 0..3 {
        close(actual.active_intensity_w_m2[axis], 0.0, flux_scale);
        close(actual.reactive_intensity_w_m2[axis], 0.0, flux_scale);
    }
}

#[test]
fn zero_stress_preserves_the_existing_compact_only_fields_and_normalization() {
    let positions = [[0.03, -0.01, 0.02], [-0.06, 0.04, -0.03]];
    let stress_positions = [[-0.1, 0.2, 0.0]];
    let mut joint = JointSourceWelchStream::new(
        acoustic(),
        &positions,
        &stress_positions,
        clock(),
        128,
        config(),
        budget(),
        || false,
    )
    .unwrap();
    let mut compact = CompactSourceWelchStream::new(
        acoustic(),
        &positions,
        clock(),
        128,
        config(),
        CompactSourceRecordingBudget {
            maximum_sources: 2,
            welch: budget().welch,
        },
        || false,
    )
    .unwrap();
    for n in 0..128 {
        let (s, c) = phase(n);
        let rows = [
            [0.001 * c + 0.0002 * s, 0.4 * s, -0.1 * c, 0.2 * c],
            [-0.0004 * s, 0.1 * c, 0.3 * c, -0.4 * s],
        ];
        let time = clock().time_at(n).unwrap();
        joint
            .push_sample(
                &positions,
                &stress_positions,
                time,
                &rows,
                &[[0.0; 6]],
                || false,
            )
            .unwrap();
        compact
            .push_sample(&positions, time, &rows, || false)
            .unwrap();
    }
    let joint = joint.finish().unwrap();
    let compact = compact.finish().unwrap();
    let receiver = [0.3, -0.2, 0.4];
    let pair = joint.prepare_observers(receiver, BIN, || false).unwrap();
    let compact_observer = compact.prepare_observer(receiver, BIN, || false).unwrap();
    let work = CoherenceBudget {
        max_modes: 3,
        max_source_mode_pairs: 9,
    };
    let actual = joint
        .evaluate_bin(&pair.0, &pair.1, BIN, work, || false)
        .unwrap();
    let expected = compact
        .evaluate_bin(&compact_observer, BIN, work, || false)
        .unwrap();
    assert_eq!(
        actual.mean_square_pressure_pa2,
        expected.mean_square_pressure_pa2
    );
    assert_eq!(
        actual.velocity_covariance_m2_s2,
        expected.velocity_covariance_m2_s2
    );
    assert_eq!(actual.active_intensity_w_m2, expected.active_intensity_w_m2);
    assert_eq!(
        actual.reactive_intensity_w_m2,
        expected.reactive_intensity_w_m2
    );
}

fn frame(n: usize) -> ([f64; 4], [f64; 6]) {
    let (s, c) = phase(n);
    (
        [0.001 * c + 0.0002 * s, c, -0.3 * s, 0.5 * c],
        [c, -0.4 * s, 0.3 * c, 0.2 * s, -0.7 * c, 0.6 * s],
    )
}
fn prefix() -> JointSourceWelchStream {
    let mut recording = JointSourceWelchStream::new(
        acoustic(),
        &ORIGIN,
        &ORIGIN,
        clock(),
        64,
        config(),
        budget(),
        || false,
    )
    .unwrap();
    for n in 0..63 {
        let (compact, stress) = frame(n);
        recording
            .push_sample(
                &ORIGIN,
                &ORIGIN,
                clock().time_at(n).unwrap(),
                &[compact],
                &[stress],
                || false,
            )
            .unwrap();
    }
    recording
}

#[test]
fn joint_layout_work_and_atomic_retry_remain_bound_to_original_sources() {
    let mut recording = prefix();
    let time = clock().time_at(63).unwrap();
    let (compact, stress) = frame(63);
    assert!(
        recording
            .push_sample(
                &ORIGIN,
                &[[0.1, 0.0, 0.0]],
                time,
                &[compact],
                &[stress],
                || false
            )
            .is_err()
    );
    assert!(
        recording
            .push_sample(&ORIGIN, &ORIGIN, time + DT, &[compact], &[stress], || false)
            .is_err()
    );
    assert!(
        recording
            .push_sample(&ORIGIN, &ORIGIN, time, &[compact], &[], || false)
            .is_err()
    );
    let mut nonfinite = stress;
    nonfinite[5] = f64::NAN;
    assert!(
        recording
            .push_sample(&ORIGIN, &ORIGIN, time, &[compact], &[nonfinite], || false)
            .is_err()
    );
    let mut reference = prefix();
    let mut polls = 0;
    reference
        .push_sample(&ORIGIN, &ORIGIN, time, &[compact], &[stress], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut retry_polls = 0;
    assert!(matches!(
        recording.push_sample(&ORIGIN, &ORIGIN, time, &[compact], &[stress], || {
            retry_polls += 1;
            retry_polls == polls
        }),
        Err(AeroacError::Welch(fs_fft::welch::WelchError::Cancelled))
    ));
    assert_eq!(recording.accepted_sample_count(), 63);
    assert_eq!(recording.completed_segment_count(), 0);
    recording
        .push_sample(&ORIGIN, &ORIGIN, time, &[compact], &[stress], || false)
        .unwrap();
    let recorded = recording.finish().unwrap();
    let reference = reference.finish().unwrap();
    for bin in 0..=32 {
        assert_eq!(
            recorded
                .spectrum()
                .factors(bin)
                .unwrap()
                .collect::<Vec<_>>(),
            reference
                .spectrum()
                .factors(bin)
                .unwrap()
                .collect::<Vec<_>>()
        );
    }
    let receiver = [0.3, -0.2, 0.4];
    let pair = recorded.prepare_observers(receiver, BIN, || false).unwrap();
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    assert!(
        recorded
            .evaluate_bin(&pair.0, &pair.1, BIN + 1, work, || false)
            .is_err()
    );
    assert!(
        recorded
            .evaluate_bin(
                &pair.0,
                &pair.1,
                BIN,
                CoherenceBudget {
                    max_source_mode_pairs: 1,
                    ..work
                },
                || false
            )
            .is_err()
    );
    let other = PreparedQuadrupoleObserver::new(
        acoustic(),
        pair.0.frequency_hz(),
        receiver,
        &[[0.1, 0.0, 0.0]],
        1,
    )
    .unwrap();
    assert!(
        recorded
            .evaluate_bin(&pair.0, &other, BIN, work, || false)
            .is_err()
    );
    let wrong_medium = FluidAcousticState::try_new(1.3, 343.0).unwrap();
    let other = PreparedCompactSourceObserver::new(
        wrong_medium,
        pair.0.frequency_hz(),
        receiver,
        &ORIGIN,
        1,
    )
    .unwrap();
    assert!(
        recorded
            .evaluate_bin(&other, &pair.1, BIN, work, || false)
            .is_err()
    );
    let mut limits = budget();
    limits.welch.maximum_channels = 9;
    assert!(
        JointSourceWelchStream::new(
            acoustic(),
            &ORIGIN,
            &ORIGIN,
            clock(),
            64,
            config(),
            limits,
            || false
        )
        .is_err()
    );
    limits = budget();
    limits.welch.maximum_workspace_bytes = size_of::<[f64; 3]>();
    assert!(
        JointSourceWelchStream::new(
            acoustic(),
            &ORIGIN,
            &ORIGIN,
            clock(),
            64,
            config(),
            limits,
            || false
        )
        .is_err()
    );
}
