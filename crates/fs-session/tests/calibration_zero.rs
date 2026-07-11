//! Zero-prediction calibration telemetry (bead gp3.21): rows the ratio
//! quantiles cannot see must be REPORTED, never silently dropped — and
//! no ratio is ever invented for them.

use fs_session::{CalibrationHealth, CalibrationPolicy, CalibrationReport, Estimate};

fn estimate_with(p50: f64, unmodeled: &[&str]) -> Estimate {
    Estimate {
        wall_p10_s: 0.0,
        wall_p50_s: p50,
        wall_p90_s: p50,
        mem_ask_bytes: None,
        energy_j: 0.0,
        unmodeled_ops: unmodeled.iter().map(ToString::to_string).collect(),
    }
}

/// ALL-ZERO: quantiles stay None (no invented ratios), but the JSON and
/// summary carry the counts and the raw actual-time distribution.
#[test]
fn all_zero_predictions_are_reported_not_hidden() {
    let calibration = CalibrationReport::new();
    calibration
        .record(&estimate_with(0.0, &[]), 3.0)
        .expect("true zero-cost row records");
    calibration
        .record(&estimate_with(0.0, &["fluid.solve"]), 9.0)
        .expect("unmodeled zero row records");
    assert!(
        calibration.ratio_quantiles().is_none(),
        "no ratio exists for zero predictions"
    );
    let zero = calibration.zero_prediction_summary();
    assert_eq!(zero.true_zero, 1, "fully modeled zero = true zero cost");
    assert_eq!(zero.unmodeled, 1, "unmodeled zero = coverage gap");
    let (q10, q50, q90) = zero.actual_quantiles_s.expect("actual distribution");
    assert!(q10 >= 3.0 && q50 >= 3.0 && q90 <= 9.0);
    let json = calibration.to_json();
    assert!(json.contains("\"zero_predictions\":{\"true_zero\":1,\"unmodeled\":1"));
    assert!(!json.contains("NaN") && !json.contains("inf"));
}

/// MIXED: positive rows keep their exact quantiles; zero rows appear
/// only in the zero section — the summary can no longer look healthier
/// than the evidence.
#[test]
fn mixed_zero_rows_do_not_distort_the_quantiles() {
    let calibration = CalibrationReport::new();
    for _ in 0..3 {
        calibration
            .record(&estimate_with(2.0, &[]), 4.0)
            .expect("positive row");
    }
    calibration
        .record(&estimate_with(0.0, &["a.op"]), 100.0)
        .expect("zero row");
    let (q10, q50, q90) = calibration.ratio_quantiles().expect("positive rows exist");
    assert_eq!(
        (q10, q50, q90),
        (2.0, 2.0, 2.0),
        "ratios from positive rows only"
    );
    let zero = calibration.zero_prediction_summary();
    assert_eq!((zero.true_zero, zero.unmodeled), (0, 1));
    assert_eq!(zero.actual_quantiles_s, Some((100.0, 100.0, 100.0)));
    // Governance: a quarter of the mass invisible is the default limit;
    // this report sits at exactly 1/4 (healthy), one more zero row tips it.
    let policy = CalibrationPolicy::default();
    assert_eq!(
        calibration.health(&policy).expect("valid policy"),
        CalibrationHealth::Healthy
    );
    calibration
        .record(&estimate_with(0.0, &[]), 1.0)
        .expect("second zero row");
    match calibration.health(&policy).expect("valid policy") {
        CalibrationHealth::Degraded {
            zero_fraction,
            limit,
        } => {
            assert!(zero_fraction > limit);
        }
        CalibrationHealth::Healthy => panic!("2/5 zero rows must degrade under a 0.25 limit"),
    }
    // An unusable threshold refuses instead of certifying.
    for bad in [f64::NAN, -0.1, 1.5] {
        assert!(
            calibration
                .health(&CalibrationPolicy {
                    max_zero_prediction_fraction: bad
                })
                .is_err(),
            "{bad} must refuse"
        );
    }
}

/// REPLAY: identical row sequences render byte-identical JSON, so the
/// telemetry is ledger-stable.
#[test]
fn replayed_rows_render_identical_json() {
    let build = || {
        let calibration = CalibrationReport::new();
        calibration
            .record(&estimate_with(1.5, &[]), 3.0)
            .expect("row");
        calibration
            .record(&estimate_with(0.0, &["x.y"]), 7.0)
            .expect("row");
        calibration.to_json()
    };
    assert_eq!(build(), build(), "replay is byte-identical");
}

/// NON-FINITE: refusal happens at record() before anything can poison
/// the JSON — including the zero-prediction actual distribution.
#[test]
fn non_finite_rows_are_refused_before_the_telemetry() {
    let calibration = CalibrationReport::new();
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(
            calibration.record(&estimate_with(0.0, &[]), bad).is_err(),
            "actual {bad} must refuse"
        );
        assert!(
            calibration.record(&estimate_with(bad, &[]), 1.0).is_err(),
            "prediction {bad} must refuse"
        );
    }
    let zero = calibration.zero_prediction_summary();
    assert_eq!((zero.true_zero, zero.unmodeled), (0, 0));
    assert!(zero.actual_quantiles_s.is_none());
    assert!(calibration.to_json().contains("\"rows\":[]"));
}
