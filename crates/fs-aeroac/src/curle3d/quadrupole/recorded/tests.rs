use super::*;
use crate::curle3d::quadrupole::resolved_reynolds_moment_n_m;
use fs_fft::welch::{WelchDetrend, WelchWindow};

const POSITIONS: [[f64; 3]; 2] = [[0.01, 0.02, 0.03], [0.04, -0.02, 0.01]];
const DT: f64 = 1.0 / 8192.0;
const EPOCH: f64 = 0.125;

fn config() -> WelchConfig {
    WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
    }
}

fn budget() -> QuadrupoleRecordingBudget {
    QuadrupoleRecordingBudget {
        maximum_sources: 2,
        welch: WelchBudget {
            maximum_channels: 12,
            maximum_samples_per_channel: 256,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
    }
}

fn clock() -> QuadrupoleSampleClock {
    QuadrupoleSampleClock::new(EPOCH, DT, 0.0).unwrap()
}

fn moment(j: usize, cell: usize) -> [f64; 6] {
    let theta = core::f64::consts::TAU * j as f64 / 64.0;
    let u = [
        0.4 + (3.0 * theta + cell as f64 * 0.2).cos(),
        (2.0 * theta).sin() * (1.0 + cell as f64),
        0.2 * (7.0 * theta).cos(),
    ];
    resolved_reynolds_moment_n_m(1.2, 0.001, u).unwrap()
}

fn stream(samples: usize) -> QuadrupoleWelchStream {
    QuadrupoleWelchStream::new(&POSITIONS, clock(), samples, config(), budget(), || false).unwrap()
}

fn append(stream: &mut QuadrupoleWelchStream, first: usize, end: usize) {
    for j in first..end {
        stream
            .push_sample(
                &POSITIONS,
                clock().time_at(j).unwrap(),
                &[moment(j, 0), moment(j, 1)],
                || false,
            )
            .unwrap();
    }
}

fn assert_matches_batch(recorded: &RecordedQuadrupoleSpectrum, samples: usize) {
    let batch = WelchCrossSpectrum::estimate_with(
        12,
        samples,
        |channel, j| moment(j, channel / 6)[channel % 6],
        DT,
        config(),
        budget().welch,
        || false,
    )
    .unwrap();
    for bin in 0..=config().segment_samples / 2 {
        assert_eq!(
            recorded.spectrum().frequency_hz(bin),
            batch.frequency_hz(bin)
        );
        let actual: Vec<_> = recorded.spectrum().factors(bin).unwrap().collect();
        let expected: Vec<_> = batch.factors(bin).unwrap().collect();
        assert_eq!(actual, expected, "bin {bin}: source-phase factors changed");
    }
}

#[test]
fn g5_streaming_physical_tensors_match_original_joint_welch_with_overlap_and_tail() {
    let mut acquisition = stream(101);
    append(&mut acquisition, 0, 101);
    assert_eq!(acquisition.accepted_sample_count(), 101);
    assert_eq!(acquisition.completed_segment_count(), 2);
    let recorded = acquisition.finish().unwrap();
    assert_eq!(recorded.clock(), clock());
    assert_eq!(recorded.source_positions_m(), &POSITIONS);
    assert_matches_batch(&recorded, 101);
}

#[test]
fn g4_sample_clock_geometry_and_mid_transform_refusals_preserve_the_accepted_record() {
    let mut acquisition = stream(64);
    append(&mut acquisition, 0, 63);
    let time = clock().time_at(63).unwrap();
    let moments = [moment(63, 0), moment(63, 1)];
    assert!(
        acquisition
            .push_sample(&POSITIONS, time + DT, &moments, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&POSITIONS, time - DT, &moments, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&[POSITIONS[1], POSITIONS[0]], time, &moments, || false)
            .is_err()
    );
    assert!(
        acquisition
            .push_sample(&POSITIONS, time, &moments[..1], || false)
            .is_err()
    );
    let mut bad = moments;
    bad[1][5] = f64::NAN;
    assert!(matches!(
        acquisition.push_sample(&POSITIONS, time, &bad, || false),
        Err(AeroacError::Welch(_))
    ));
    let mut polls = 0;
    assert!(
        acquisition
            .push_sample(&POSITIONS, time, &moments, || {
                polls += 1;
                polls > 20
            })
            .is_err()
    );
    assert!(polls > 20);
    assert_eq!(acquisition.accepted_sample_count(), 63);
    assert_eq!(acquisition.completed_segment_count(), 0);
    acquisition
        .push_sample(&POSITIONS, time, &moments, || false)
        .unwrap();
    assert!(
        acquisition
            .push_sample(&POSITIONS, time + DT, &moments, || false)
            .is_err()
    );
    assert_matches_batch(&acquisition.finish().unwrap(), 64);
}

#[test]
fn g1_recorded_geometry_binds_the_original_pressure_and_vector_velocity_observer() {
    let mut acquisition = stream(64);
    append(&mut acquisition, 0, 64);
    let recorded = acquisition.finish().unwrap();
    let acoustic = FluidAcousticState::try_new(1.2, 343.0).unwrap();
    let observer = recorded
        .prepare_observer(acoustic, [0.3, 0.2, 0.1], 6, || false)
        .unwrap();
    let work = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    let field = recorded.evaluate_bin(&observer, 6, work, || false).unwrap();
    assert!(field.mean_square_pressure_pa2 > 0.0);
    assert!(
        field
            .velocity_covariance_m2_s2
            .iter()
            .flatten()
            .all(|v| v.re.is_finite() && v.im.is_finite())
    );
    let other = PreparedQuadrupoleObserver::new(
        acoustic,
        observer.frequency_hz(),
        [0.3, 0.2, 0.1],
        &[POSITIONS[1], POSITIONS[0]],
        2,
    )
    .unwrap();
    assert!(recorded.evaluate_bin(&other, 6, work, || false).is_err());
    assert!(recorded.evaluate_bin(&observer, 7, work, || false).is_err());
    assert!(recorded.evaluate_bin(&observer, 6, work, || true).is_err());
    assert!(
        recorded
            .prepare_observer(acoustic, [0.3, 0.2, 0.1], 0, || false)
            .is_err()
    );
    assert!(
        recorded
            .prepare_observer(acoustic, [0.3, 0.2, 0.1], 32, || false)
            .is_err()
    );
}

#[test]
fn g0_source_clock_budget_and_incomplete_record_admission_are_explicit() {
    assert!(QuadrupoleSampleClock::new(0.0, 0.0, 0.0).is_err());
    assert!(QuadrupoleSampleClock::new(0.0, DT, DT * 0.25).is_err());
    // dt=0.75 ULP can pass one adjacent pair yet later repeat a timestamp.
    assert!(QuadrupoleSampleClock::new(2.0f64.powi(50), 0.1875, 0.0).is_err());
    assert!(QuadrupoleSampleClock::new(1e20, DT, 0.0).is_err());
    let tolerant = QuadrupoleSampleClock::new(0.0, DT, DT / 16.0).unwrap();
    assert!(
        tolerant
            .admit(10, tolerant.time_at(10).unwrap() + DT / 32.0)
            .is_ok()
    );
    assert!(
        tolerant
            .admit(10, tolerant.time_at(10).unwrap() + DT / 8.0)
            .is_err()
    );
    let mut tiny = budget();
    tiny.maximum_sources = 1;
    assert!(QuadrupoleWelchStream::new(&POSITIONS, clock(), 64, config(), tiny, || false).is_err());
    tiny = budget();
    tiny.welch.maximum_workspace_bytes = 47;
    assert!(QuadrupoleWelchStream::new(&POSITIONS, clock(), 64, config(), tiny, || false).is_err());
    tiny = budget();
    tiny.welch.maximum_channels = 11;
    assert!(QuadrupoleWelchStream::new(&POSITIONS, clock(), 64, config(), tiny, || false).is_err());
    assert!(
        QuadrupoleWelchStream::new(&POSITIONS, clock(), 64, config(), budget(), || true).is_err()
    );
    assert!(stream(64).finish().is_err());
}
