//! End-to-end battery: BO finds the tilted-double-well optimum and stops on an
//! anytime-valid e-process certificate before exhausting its budget.

use fs_adaptbo_e2e::{objective, run_campaign};
use fs_evidence::{Color, StatisticalCertificate};

#[test]
fn bo_converges_and_stops_with_an_anytime_valid_certificate() {
    let report = run_campaign(30, 0.02, 0.05);
    // CONVERGENCE: the global minimum sits near x ≈ 3 (value ≈ −0.45).
    assert!(
        (report.best_x - 3.0).abs() < 0.3,
        "best_x {}",
        report.best_x
    );
    assert!(
        report.best_value < objective(1.0),
        "did not beat the other well"
    );
    assert!(report.best_value < -0.3, "best value {}", report.best_value);
    // ANYTIME-VALID STOP: the e-process rejected before the iteration cap.
    assert!(
        report.stopped_early,
        "did not stop early (iters {})",
        report.iterations
    );
    assert!(report.iterations < 30, "used the whole budget");
    // the e-value crossed the Ville threshold ln(1/α).
    assert!(
        report.log_e_value >= -(0.05_f64).ln() - 1e-9,
        "log-e {}",
        report.log_e_value
    );
    assert!(matches!(
        report.stop_evidence,
        StatisticalCertificate::EValue { e, alpha }
            if e >= 1.0 / alpha && alpha.to_bits() == 0.05_f64.to_bits()
    ));
    assert!(matches!(report.incumbent_color, Color::Estimated { .. }));
    assert!(matches!(report.surrogate_color, Color::Estimated { .. }));
    // A trace diagnostic is reported, not mislabeled as an optimum enclosure.
    assert!(report.ci_radius > 0.0 && report.ci_radius.is_finite());
    assert!(report.evaluations >= report.iterations + 3);
    println!(
        "{{\"campaign\":\"anytimebo\",\"best_x\":{:.3},\"best_value\":{:.4},\"iterations\":{},\
         \"evaluations\":{},\"stopped_early\":{},\"log_e_value\":{:.3},\"ci\":[{:.3},{:.3}]}}",
        report.best_x,
        report.best_value,
        report.iterations,
        report.evaluations,
        report.stopped_early,
        report.log_e_value,
        report.ci_center,
        report.ci_radius,
    );
}

#[test]
fn a_tiny_delta_never_declares_a_stall() {
    // with an impossibly small improvement threshold, no step counts as a stall,
    // so the e-process never rejects — the search runs to the cap.
    let report = run_campaign(8, 1e-12, 0.05);
    assert!(!report.stopped_early || report.iterations == 8);
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(30, 0.02, 0.05);
    let b = run_campaign(30, 0.02, 0.05);
    assert_eq!(a.best_x.to_bits(), b.best_x.to_bits());
    assert_eq!(a.iterations, b.iterations);
    assert_eq!(a.log_e_value.to_bits(), b.log_e_value.to_bits());
}

#[test]
fn malformed_or_unbounded_campaign_controls_are_refused() {
    assert!(std::panic::catch_unwind(|| run_campaign(65, 0.02, 0.05)).is_err());
    assert!(std::panic::catch_unwind(|| run_campaign(30, f64::NAN, 0.05)).is_err());
    assert!(std::panic::catch_unwind(|| run_campaign(30, -0.1, 0.05)).is_err());
    assert!(std::panic::catch_unwind(|| run_campaign(30, 0.02, 0.0)).is_err());
}
