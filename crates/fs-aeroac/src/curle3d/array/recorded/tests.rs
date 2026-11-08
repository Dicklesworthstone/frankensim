use super::*;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor};
use fs_fft::welch::{WelchDetrend, WelchWindow};
use fs_math::det;

const DT: f64 = 1.0 / 8192.0;
const EPOCH: f64 = 0.125;
const POSITIONS: [[f64; 3]; 2] = [[0.03, -0.01, 0.02], [-0.06, 0.04, -0.03]];
const RECEIVER: [f64; 3] = [0.3, -0.2, 0.4];
const BIN: usize = 5;

fn acoustic() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 343.0).unwrap()
}

fn clock() -> SourceSampleClock {
    SourceSampleClock::new(EPOCH, DT, 0.0).unwrap()
}

fn config() -> WelchConfig {
    WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    }
}

fn budget() -> CompactSourceRecordingBudget {
    CompactSourceRecordingBudget {
        maximum_sources: 2,
        welch: WelchBudget {
            maximum_channels: 8,
            maximum_samples_per_channel: 256,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
    }
}

fn phasors() -> [CompactSourcePhasor; 2] {
    [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(0.003, 0.0005),
            force_on_fluid_n: [
                C64::new(0.3, 0.2),
                C64::new(-0.1, 0.4),
                C64::new(0.15, -0.2),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-0.001, 0.002),
            force_on_fluid_n: [
                C64::new(-0.2, 0.1),
                C64::new(0.25, -0.08),
                C64::new(0.35, 0.12),
            ],
        },
    ]
}

fn values(sample: usize) -> [[f64; 4]; 2] {
    let (s, c) = det::sin_cos(core::f64::consts::TAU * BIN as f64 * sample as f64 / 64.0);
    // Re[P exp(-i omega t)] = Re(P) cos + Im(P) sin. The original
    // Welch owner conjugates its forward FFT exactly once for this convention.
    phasors().map(|source| source.components().map(|z| z.re * c + z.im * s))
}

fn stream(samples: usize, config: WelchConfig) -> CompactSourceWelchStream {
    CompactSourceWelchStream::new(
        acoustic(),
        &POSITIONS,
        clock(),
        samples,
        config,
        budget(),
        || false,
    )
    .unwrap()
}

fn append(recording: &mut CompactSourceWelchStream, first: usize, end: usize) {
    for n in first..end {
        recording
            .push_sample(&POSITIONS, clock().time_at(n).unwrap(), &values(n), || {
                false
            })
            .unwrap();
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * scale.max(1e-24),
        "{actual:e} != {expected:e}"
    );
}

fn field_matches_tone(actual: &CompactFieldStatistics, tone: &CoherentCompactField) {
    let pressure_square = 0.5 * tone.pressure_pa.norm_sq();
    close(
        actual.mean_square_pressure_pa2,
        pressure_square,
        pressure_square,
    );
    for i in 0..3 {
        let intensity_scale = 0.5 * tone.pressure_pa.abs() * tone.particle_velocity_m_s[i].abs();
        close(
            actual.active_intensity_w_m2[i],
            tone.active_intensity_w_m2[i],
            intensity_scale,
        );
        close(
            actual.reactive_intensity_w_m2[i],
            tone.reactive_intensity_w_m2[i],
            intensity_scale,
        );
        for j in 0..3 {
            let expected =
                (tone.particle_velocity_m_s[i] * tone.particle_velocity_m_s[j].conj()).scale(0.5);
            close(
                actual.velocity_covariance_m2_s2[i][j].re,
                expected.re,
                expected.abs(),
            );
            close(
                actual.velocity_covariance_m2_s2[i][j].im,
                expected.im,
                expected.abs(),
            );
        }
    }
}

#[test]
fn g1_phased_volume_force_records_match_complete_harmonic_fields_with_overlap() {
    let mut acquisition = stream(128, config());
    append(&mut acquisition, 0, 128);
    let recorded = acquisition.finish().unwrap();
    assert_eq!(recorded.clock(), clock());
    assert_eq!(recorded.acoustic(), acoustic());
    assert_eq!(recorded.source_positions_m(), &POSITIONS);
    let observer = recorded.prepare_observer(RECEIVER, BIN, || false).unwrap();
    let actual = recorded
        .evaluate_bin(
            &observer,
            BIN,
            CoherenceBudget {
                max_modes: 3,
                max_source_mode_pairs: 6,
            },
            || false,
        )
        .unwrap();
    field_matches_tone(&actual, &observer.evaluate(&phasors()).unwrap());
    let separate_power: f64 = phasors()
        .into_iter()
        .enumerate()
        .map(|(i, source)| {
            let mut independent = [CompactSourcePhasor::default(); 2];
            independent[i] = source;
            0.5 * observer
                .evaluate(&independent)
                .unwrap()
                .pressure_pa
                .norm_sq()
        })
        .sum();
    assert!((actual.mean_square_pressure_pa2 - separate_power).abs() > 0.01 * separate_power);
}

#[test]
fn g3_pressure_null_keeps_the_recorded_particle_velocity_field() {
    let positions = [[-0.1, 0.0, 0.0], [0.1, 0.0, 0.0]];
    let mut acquisition = CompactSourceWelchStream::new(
        acoustic(),
        &positions,
        clock(),
        64,
        config(),
        budget(),
        || false,
    )
    .unwrap();
    for n in 0..64 {
        let q = 0.001 * det::cos(core::f64::consts::TAU * BIN as f64 * n as f64 / 64.0);
        acquisition
            .push_sample(
                &positions,
                clock().time_at(n).unwrap(),
                &[[q, 0.0, 0.0, 0.0], [-q, 0.0, 0.0, 0.0]],
                || false,
            )
            .unwrap();
    }
    let recorded = acquisition.finish().unwrap();
    let observer = recorded
        .prepare_observer([0.0, 0.0, 0.2], BIN, || false)
        .unwrap();
    let actual = recorded
        .evaluate_bin(
            &observer,
            BIN,
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: 2,
            },
            || false,
        )
        .unwrap();
    assert_eq!(actual.mean_square_pressure_pa2, 0.0);
    assert!(actual.velocity_covariance_m2_s2[0][0].re > 1e-8);
    assert_eq!(actual.active_intensity_w_m2, [0.0; 3]);
    assert_eq!(actual.reactive_intensity_w_m2, [0.0; 3]);
}

fn assert_batch_factors(
    recorded: &RecordedCompactSourceSpectrum,
    samples: usize,
    config: WelchConfig,
) {
    let batch = WelchCrossSpectrum::estimate_with(
        8,
        samples,
        |channel, n| values(n)[channel / 4][channel % 4],
        DT,
        config,
        budget().welch,
        || false,
    )
    .unwrap();
    for bin in 0..=config.segment_samples / 2 {
        assert_eq!(
            recorded
                .spectrum()
                .factors(bin)
                .unwrap()
                .collect::<Vec<_>>(),
            batch.factors(bin).unwrap().collect::<Vec<_>>()
        );
    }
}

#[test]
fn g5_streaming_hann_factors_match_the_original_joint_transform_with_a_partial_tail() {
    let config = WelchConfig {
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
        ..config()
    };
    let mut acquisition = stream(101, config);
    append(&mut acquisition, 0, 101);
    assert_eq!(acquisition.accepted_sample_count(), 101);
    assert_eq!(acquisition.completed_segment_count(), 2);
    assert_batch_factors(&acquisition.finish().unwrap(), 101, config);
}

#[test]
fn g4_bad_samples_and_mid_transform_cancellation_preserve_the_accepted_record() {
    let mut acquisition = stream(64, config());
    append(&mut acquisition, 0, 63);
    let time = clock().time_at(63).unwrap();
    let sample = values(63);
    assert!(
        acquisition
            .push_sample(&POSITIONS, time + DT, &sample, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&POSITIONS, time - DT, &sample, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&[POSITIONS[1], POSITIONS[0]], time, &sample, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&POSITIONS, time, &sample[..1], || false)
            .is_err()
    );
    let mut bad = sample;
    bad[1][3] = f64::NAN;
    assert!(
        acquisition
            .push_sample(&POSITIONS, time, &bad, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        acquisition
            .push_sample(&POSITIONS, time, &sample, || {
                polls += 1;
                polls > 20
            })
            .is_err()
    );
    assert!(polls > 20);
    assert_eq!(acquisition.accepted_sample_count(), 63);
    assert_eq!(acquisition.completed_segment_count(), 0);
    acquisition
        .push_sample(&POSITIONS, time, &sample, || false)
        .unwrap();
    assert!(
        acquisition
            .push_sample(&POSITIONS, time + DT, &sample, || false)
            .is_err()
    );
    assert_batch_factors(&acquisition.finish().unwrap(), 64, config());
}

#[test]
fn g0_recording_and_field_admission_keep_medium_geometry_clock_and_work_limits() {
    let build = |positions: &[[f64; 3]], limits| {
        CompactSourceWelchStream::new(acoustic(), positions, clock(), 64, config(), limits, || {
            false
        })
    };
    assert!(build(&[], budget()).is_err());
    assert!(
        build(
            &POSITIONS,
            CompactSourceRecordingBudget {
                maximum_sources: 1,
                ..budget()
            }
        )
        .is_err()
    );
    let mut limits = budget();
    limits.welch.maximum_channels = 7;
    assert!(build(&POSITIONS, limits).is_err());
    limits = budget();
    limits.welch.maximum_workspace_bytes = size_of::<[f64; 3]>();
    assert!(build(&POSITIONS, limits).is_err());
    assert!(build(&[[f64::NAN, 0.0, 0.0]], budget()).is_err());
    assert!(stream(64, config()).finish().is_err());
    let mut acquisition = stream(64, config());
    append(&mut acquisition, 0, 64);
    let recorded = acquisition.finish().unwrap();
    assert!(recorded.prepare_observer(RECEIVER, 0, || false).is_err());
    assert!(recorded.prepare_observer(RECEIVER, 32, || false).is_err());
    assert!(recorded.prepare_observer(RECEIVER, BIN, || true).is_err());
    let observer = recorded.prepare_observer(RECEIVER, BIN, || false).unwrap();
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    assert!(
        recorded
            .evaluate_bin(&observer, BIN + 1, work, || false)
            .is_err()
    );
    assert!(
        recorded
            .evaluate_bin(
                &observer,
                BIN,
                CoherenceBudget {
                    max_source_mode_pairs: 1,
                    ..work
                },
                || false
            )
            .is_err()
    );
    for (medium, positions) in [
        (FluidAcousticState::try_new(1.3, 343.0).unwrap(), POSITIONS),
        (acoustic(), [POSITIONS[1], POSITIONS[0]]),
    ] {
        let other = PreparedCompactSourceObserver::new(
            medium,
            observer.frequency_hz(),
            RECEIVER,
            &positions,
            2,
        )
        .unwrap();
        assert!(recorded.evaluate_bin(&other, BIN, work, || false).is_err());
    }
}
