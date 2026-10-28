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
