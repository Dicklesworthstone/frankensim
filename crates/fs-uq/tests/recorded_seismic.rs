use fs_uq::seismic::{LinearSdofControls, LinearSdofError};
use fs_uq::{
    AccelerationUnit, BilinearIdaModel, GroundMotionRecord, GroundMotionRecordError,
    GroundMotionSampling, IdaIntensity, KanaiTajimi, ida_fragility_from_records,
    ida_fragility_with_sampling,
};

fn record(dt: f64, samples: Vec<f64>) -> GroundMotionRecord {
    GroundMotionRecord::try_new(dt, samples, AccelerationUnit::MetresPerSecondSquared).unwrap()
}

fn model() -> BilinearIdaModel {
    BilinearIdaModel {
        angular_frequency_rad_s: 12.0,
        damping_ratio: 0.05,
        yield_displacement_m: 0.005,
        post_yield_stiffness_ratio: 0.1,
        ductility_limit: 3.0,
    }
}

const PEER: &str = "PEER NGA STRONG MOTION DATABASE RECORD
SYNTHETIC FORMAT FIXTURE, station, component 90
ACCELERATION TIME SERIES IN UNITS OF G
NPTS= 5, DT= .0200 SEC,
-1.0E-1  5.0D-2 0
2.5e-2 -5.0d-4
";

#[test]
fn peer_and_legacy_headers_preserve_signed_samples_clock_and_units() {
    let parsed = GroundMotionRecord::from_peer_at2(PEER, 5).unwrap();
    assert_eq!(parsed.dt_s(), 0.02);
    assert_eq!(parsed.duration_s(), 0.08);
    assert_eq!(
        parsed.description(),
        Some("SYNTHETIC FORMAT FIXTURE, station, component 90")
    );
    for (actual, value_g) in parsed
        .acceleration_m_s2()
        .iter()
        .zip([-0.1, 0.05, 0.0, 0.025, -0.0005])
    {
        assert!((actual - 9.80665 * value_g).abs() < 2e-16);
    }
    let legacy = PEER
        .replace(
            "PEER NGA STRONG MOTION DATABASE RECORD",
            "PACIFIC ENGINEERING AND ANALYSIS STRONG-MOTION DATA",
        )
        .replace("TIME SERIES", "TIME HISTORY")
        .replace(
            "UNITS OF G",
            "UNITS OF G. FILTER POINTS: HP=0.1 Hz LP=20 Hz",
        );
    assert_eq!(
        GroundMotionRecord::from_peer_at2(&legacy, 5).unwrap(),
        parsed
    );
    let cm = GroundMotionRecord::try_new(
        0.02,
        vec![-98.0665, 49.03325, 0.0, 24.516625, -0.4903325],
        AccelerationUnit::CentimetresPerSecondSquared,
    )
    .unwrap();
    for (a, b) in cm
        .acceleration_m_s2()
        .iter()
        .zip(parsed.acceleration_m_s2())
    {
        assert!((a - b).abs() < 3e-16);
    }
}

#[test]
fn ingestion_refuses_bad_units_counts_nonfinite_values_and_sample_budget() {
    assert_eq!(
        GroundMotionRecord::from_peer_at2(PEER, 4),
        Err(GroundMotionRecordError::SampleLimit {
            requested: 5,
            maximum: 4
        })
    );
    for text in [
        PEER.replace("ACCELERATION", "VELOCITY"),
        PEER.replace("UNITS OF G", "UNITS OF CM/SEC"),
        PEER.replace("DT= .0200", "DT= 0"),
        PEER.replace("SEC,", "MSEC,"),
        PEER.replace("NPTS= 5", "NPTS= 5 NPTS= 5"),
    ] {
        assert!(GroundMotionRecord::from_peer_at2(&text, 10).is_err());
    }
    assert_eq!(
        GroundMotionRecord::from_peer_at2(&PEER.replace("NPTS= 5", "NPTS= 6"), 10),
        Err(GroundMotionRecordError::SampleCount {
            expected: 6,
            actual: 5
        })
    );
    assert_eq!(
        GroundMotionRecord::from_peer_at2(&PEER.replace("NPTS= 5", "NPTS= 4"), 10),
        Err(GroundMotionRecordError::SampleCount {
            expected: 4,
            actual: 5
        })
    );
    for value in ["NaN", "inf", "1e308", "not-a-number"] {
        assert_eq!(
            GroundMotionRecord::from_peer_at2(&PEER.replace("-1.0E-1", value), 5),
            Err(GroundMotionRecordError::Sample { index: 0 })
        );
    }
    assert!(
        GroundMotionRecord::try_new(0.01, vec![0.0], AccelerationUnit::MetresPerSecondSquared)
            .is_err()
    );
}

#[test]
fn distinct_record_clocks_match_independent_one_interval_solutions() {
    let clocks = [0.1, 0.2, 0.075];
    let amplitudes = [0.5, 1.0, 3.0];
    let suite: Vec<_> = clocks
        .into_iter()
        .zip(amplitudes)
        .map(|(dt, a)| record(dt, vec![0.0, -a]))
        .collect();
    let m = BilinearIdaModel {
        angular_frequency_rad_s: 10.0,
        damping_ratio: 0.2,
        yield_displacement_m: 0.1,
        post_yield_stiffness_ratio: 1.0,
        ductility_limit: 0.02,
    };
    let levels = [0.0, 0.5, 2.0];
    let report = ida_fragility_from_records(&suite, IdaIntensity::ScaleFactor, &levels, m);
    for (row, scale) in levels.into_iter().enumerate() {
        for (column, (dt, a)) in clocks.into_iter().zip(amplitudes).enumerate() {
            // M*a1 + C*v1 + K*u1 = scale*A, u1=h²*a1/4, v1=h*a1/2.
            let expected = scale * a * dt * dt / 4.0 / (1.0 + 2.0 * dt + 25.0 * dt * dt) / 0.1;
            assert!((report.peak_ductilities[row][column] - expected).abs() < 2e-16);
            assert_eq!(report.record_scales[row][column], scale);
        }
    }
    assert_eq!(
        report.fragility.iter().map(|p| p.p).collect::<Vec<_>>(),
        vec![0.0, 1.0 / 3.0, 2.0 / 3.0]
    );
    assert!(report.fragility.iter().all(|p| p.n == 3));
}

#[test]
fn pga_normalization_preserves_nonlinear_waveform_and_polarity_equivalence() {
    let base: Vec<_> = (0..=600)
        .map(|i| {
            let t = f64::from(i) * 0.01;
            (7.0 * t).sin() + 0.35 * (19.0 * t).sin()
        })
        .collect();
    let suite: Vec<_> = [1.0, 3.5, -2.0]
        .into_iter()
        .map(|factor| record(0.01, base.iter().map(|a| factor * a).collect()))
        .collect();
    let levels = [0.0, 0.3, 1.0, 3.0];
    let report = ida_fragility_from_records(
        &suite,
        IdaIntensity::PeakGroundAcceleration,
        &levels,
        model(),
    );
    for (row, &target) in levels.iter().enumerate() {
        for column in 0..suite.len() {
            let actual = report.peak_ductilities[row][column];
            let first = report.peak_ductilities[row][0];
            assert!((actual - first).abs() < 2e-10 * first.max(1.0));
            assert!(
                (report.record_scales[row][column] * suite[column].peak_ground_acceleration()
                    - target)
                    .abs()
                    < 2e-15
            );
        }
    }
    assert!(
        report.peak_ductilities.last().unwrap()[0] > 5.0,
        "the normalization regression must exercise substantial yielding"
    );
}

#[test]
fn psa_targets_recover_the_same_physical_linear_demand_across_records() {
    let suite: Vec<_> = [(0.01, 7.0), (0.015, 17.0)]
        .into_iter()
        .map(|(dt, frequency)| {
            record(
                dt,
                (0..=300)
                    .map(|i| {
                        let t = f64::from(i) * dt;
                        0.4 * (frequency * t).sin() + 0.1 * (3.0 * t).cos()
                    })
                    .collect(),
            )
        })
        .collect();
    let m = BilinearIdaModel {
        post_yield_stiffness_ratio: 1.0,
        ductility_limit: 1.0,
        ..model()
    };
    let levels = [0.0, 0.2, 1.0];
    let measure = IdaIntensity::PseudoSpectralAcceleration {
        period_s: core::f64::consts::TAU / m.angular_frequency_rad_s,
        damping_ratio: m.damping_ratio,
    };
    let report = ida_fragility_from_records(&suite, measure, &levels, m);
    assert!((report.reference_intensities[0] - report.reference_intensities[1]).abs() > 0.01);
    for (row, target) in levels.into_iter().enumerate() {
        // For this elastic oscillator Sa = wn²*peak_displacement. The expected
        // demand follows directly from the physical intensity definition,
        // independently of waveform, sample interval, and amplitude.
        let expected = target / (m.angular_frequency_rad_s.powi(2) * m.yield_displacement_m);
        for &actual in &report.peak_ductilities[row] {
            assert!((actual - expected).abs() < 2e-12 * expected.max(1.0));
        }
    }
    assert_eq!(
        report.fragility.iter().map(|p| p.p).collect::<Vec<_>>(),
        vec![0.0, 0.0, 1.0]
    );
}

#[test]
fn supplied_synthetic_histories_reproduce_existing_ida_exactly() {
    let kt = KanaiTajimi {
        s0: 0.03,
        wg: 12.0,
        zg: 0.6,
    };
    let sampling = GroundMotionSampling {
        dt_s: 0.01,
        n_steps: 320,
        n_freq: 48,
        omega_max: 70.0,
    };
    let seeds = [17, 48, 119, 322];
    let suite: Vec<_> = seeds
        .into_iter()
        .map(|seed| record(sampling.dt_s, kt.synthesize_with_sampling(seed, sampling)))
        .collect();
    let levels = [0.1, 0.5, 1.0, 2.0];
    let m = model();
    let legacy = ida_fragility_with_sampling(
        &kt,
        sampling,
        &seeds,
        &levels,
        m.angular_frequency_rad_s,
        m.damping_ratio,
        m.yield_displacement_m,
        m.post_yield_stiffness_ratio,
        m.ductility_limit,
    );
    let supplied = ida_fragility_from_records(&suite, IdaIntensity::ScaleFactor, &levels, m);
    assert_eq!(supplied.fragility, legacy);
}

#[test]
fn peer_file_drives_the_original_yielded_equilibrium_in_si_units() {
    let text = "PEER NGA STRONG MOTION DATABASE RECORD
SYNTHETIC TWO-SAMPLE YIELD FIXTURE
ACCELERATION TIME HISTORY IN UNITS OF G
NPTS=2, DT=.2 SEC
0 -1
";
    let parsed = GroundMotionRecord::from_peer_at2(text, 2).unwrap();
    let m = BilinearIdaModel {
        angular_frequency_rad_s: 10.0,
        damping_ratio: 0.0,
        yield_displacement_m: 0.001,
        post_yield_stiffness_ratio: 0.1,
        ductility_limit: 10.0,
    };
    let report = ida_fragility_from_records(&[parsed], IdaIntensity::ScaleFactor, &[1.0], m);
    // k_post=10, plastic yield force=.09, h²/4=.01, applied acceleration=g.
    let expected = 0.01 * (9.80665 - 0.09) / 1.1 / 0.001;
    assert!((report.peak_ductilities[0][0] - expected).abs() < 3e-13);
    assert_eq!(report.fragility[0].p, 1.0);
}

#[test]
fn zero_records_remain_valid_without_physical_intensity_normalization() {
    let suite = [record(0.02, vec![0.0; 100])];
    let result =
        ida_fragility_from_records(&suite, IdaIntensity::ScaleFactor, &[0.0, 3.0], model());
    assert_eq!(result.peak_ductilities, vec![vec![0.0], vec![0.0]]);
    assert!(
        std::panic::catch_unwind(|| ida_fragility_from_records(
            &suite,
            IdaIntensity::PeakGroundAcceleration,
            &[1.0],
            model()
        ))
        .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| ida_fragility_from_records(
            &suite,
            IdaIntensity::ScaleFactor,
            &[f64::NAN],
            model()
        ))
        .is_err()
    );
}

fn linear_controls() -> LinearSdofControls {
    LinearSdofControls {
        maximum_samples: 10_000,
        maximum_acceleration_zeros_per_interval: 128,
        maximum_transition_evaluations: 100_000,
        maximum_bisections_per_extremum: 128,
        extremum_time_tolerance_s: 1.0e-10,
    }
}

#[test]
fn continuous_spectrum_resolves_hidden_peaks_and_preserves_period_order() {
    // The one-second oscillator returns to rest after this whole-period
    // constant-force interval, hiding its displacement peak between samples.
    let motion = record(1.0, vec![3.0, 3.0]);
    let before = motion.clone();
    let periods = [1.0, 2.0, 0.5, 1.0];
    let report = motion
        .linear_response_spectrum(&periods, 0.0, periods.len(), linear_controls(), || false)
        .unwrap();
    assert_eq!(report.damping_ratio, 0.0);
    assert_eq!(report.duration_s, 1.0);
    for (&period, point) in periods.iter().zip(&report.points) {
        let omega = core::f64::consts::TAU / period;
        let expected_peak = 6.0 / (omega * omega);
        assert_eq!(point.period_s, period);
        assert!((point.response.peak_displacement_m - expected_peak).abs() < 5e-14);
        assert!((point.response.signed_peak_displacement_m + expected_peak).abs() < 5e-14);
        assert!((point.pseudo_spectral_velocity_m_s - 6.0 / omega).abs() < 3e-13);
        assert!((point.pseudo_spectral_acceleration_m_s2 - 6.0).abs() < 2e-12);
        assert_eq!(point.response.duration_s, 1.0);
        assert_eq!(point.response.intervals, 1);
        assert_eq!(
            point.pseudo_velocity_localization_allowance_m_s,
            omega * point.response.peak_localization_allowance_m,
        );
        assert_eq!(
            point.pseudo_acceleration_localization_allowance_m_s2,
            omega * omega * point.response.peak_localization_allowance_m,
        );
    }
    let first = report.points[0].response;
    assert!(first.final_displacement_m.abs() < 5e-14);
    assert!((first.peak_time_s - 0.5).abs() < 1e-10);
    assert_eq!(report.points[0], report.points[3]);
    assert_eq!(
        motion
            .linear_response(1.0, 0.0, linear_controls(), || false)
            .unwrap(),
        first,
    );
    assert_eq!(
        report.transition_evaluations,
        report
            .points
            .iter()
            .map(|p| p.response.transition_evaluations)
            .sum::<usize>(),
    );
    assert_eq!(motion, before);
}

#[test]
fn recorded_continuous_spectra_preserve_si_units_and_linear_interpolation() {
    let values_g = [0.0, 1.0, -0.5, 0.0];
    let motion =
        GroundMotionRecord::try_new(0.3, values_g.to_vec(), AccelerationUnit::StandardGravity)
            .unwrap();
    let centimetres = GroundMotionRecord::try_new(
        0.3,
        values_g.iter().map(|a| a * 980.665).collect(),
        AccelerationUnit::CentimetresPerSecondSquared,
    )
    .unwrap();
    // Add only points on the existing straight segments: the continuous
    // excitation is identical, despite the different record clock.
    let mut refined_samples = Vec::new();
    for pair in motion.acceleration_m_s2().windows(2) {
        for i in 0..8 {
            let weight = f64::from(i) / 8.0;
            refined_samples.push((1.0 - weight) * pair[0] + weight * pair[1]);
        }
    }
    refined_samples.push(*motion.acceleration_m_s2().last().unwrap());
    let refined = record(0.3 / 8.0, refined_samples);
    let periods = [0.4, 1.1, 2.0];
    let baseline = motion
        .linear_response_spectrum(&periods, 0.07, 3, linear_controls(), || false)
        .unwrap();
    for equivalent in [centimetres, refined] {
        let result = equivalent
            .linear_response_spectrum(&periods, 0.07, 3, linear_controls(), || false)
            .unwrap();
        assert!((result.duration_s - baseline.duration_s).abs() < 2e-16);
        for (actual, expected) in result.points.iter().zip(&baseline.points) {
            assert!(
                (actual.response.peak_displacement_m - expected.response.peak_displacement_m).abs()
                    < 2e-12
            );
            assert!(
                (actual.response.signed_peak_displacement_m
                    - expected.response.signed_peak_displacement_m)
                    .abs()
                    < 2e-12
            );
            assert!((actual.response.peak_time_s - expected.response.peak_time_s).abs() < 2e-9);
            assert!(
                (actual.pseudo_spectral_acceleration_m_s2
                    - expected.pseudo_spectral_acceleration_m_s2)
                    .abs()
                    < 5e-11
            );
        }
    }
}

#[test]
fn recorded_linear_response_ends_at_the_supplied_duration_without_decay_tail() {
    let motion = record(0.03, vec![0.0, -1.0]);
    let response = motion
        .linear_response(1.0, 0.05, linear_controls(), || false)
        .unwrap();
    assert_eq!(response.duration_s, 0.03);
    assert_eq!(response.peak_time_s, 0.03);
    assert_eq!(
        response.peak_displacement_m,
        response.final_displacement_m.abs()
    );
    assert!(response.final_velocity_m_s > 0.0);
    // The explicitly authored zero-acceleration tail is a different record.
    let extended = record(0.03, vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0]);
    let later = extended
        .linear_response(1.0, 0.05, linear_controls(), || false)
        .unwrap();
    assert!(later.peak_displacement_m > response.peak_displacement_m);
    assert!(later.peak_time_s > response.duration_s);
}

#[test]
fn recorded_spectrum_enforces_global_work_period_sample_and_cancellation_limits() {
    let motion = record(0.1, vec![0.0; 3]);
    let periods = [0.5, 1.0, 2.0];
    let mut controls = linear_controls();
    controls.maximum_transition_evaluations = 5;
    assert!(matches!(
        motion.linear_response_spectrum(&periods, 0.05, 3, controls, || false),
        Err(LinearSdofError::Budget(_)),
    ));
    controls.maximum_transition_evaluations = 6;
    let completed = motion
        .linear_response_spectrum(&periods, 0.05, 3, controls, || false)
        .unwrap();
    assert_eq!(completed.transition_evaluations, 6);
    assert!(
        completed
            .points
            .iter()
            .all(|p| p.response.transition_evaluations == 2)
    );
    assert!(matches!(
        motion.linear_response_spectrum(&[f64::NAN], 0.05, 0, controls, || false),
        Err(LinearSdofError::Budget(_)),
    ));
    controls.maximum_samples = 2;
    assert!(matches!(
        motion.linear_response_spectrum(&periods, 0.05, 3, controls, || false),
        Err(LinearSdofError::Budget(_)),
    ));
    for period in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE / 4.0] {
        assert!(
            motion
                .linear_response(period, 0.05, linear_controls(), || false)
                .is_err()
        );
    }
    assert!(
        motion
            .linear_response_spectrum(&[], 0.05, 3, linear_controls(), || false)
            .is_err()
    );
    assert!(
        motion
            .linear_response_spectrum(&periods, -0.1, 3, linear_controls(), || false)
            .is_err()
    );
    assert!(matches!(
        motion.linear_response(1.0, 0.05, linear_controls(), || true),
        Err(LinearSdofError::Cancelled),
    ));
    let active = record(0.01, vec![0.1; 50]);
    let before = active.clone();
    let mut polls = 0;
    assert!(matches!(
        active.linear_response_spectrum(&periods, 0.05, 3, linear_controls(), || {
            polls += 1;
            polls >= 80
        }),
        Err(LinearSdofError::Cancelled),
    ));
    assert_eq!(active, before);
}
