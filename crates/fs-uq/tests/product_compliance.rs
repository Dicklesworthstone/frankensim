//! Fixed-event confidence assessment, refusal, and checkpoint interoperability.
//! Bead: frankensim-extreal-program-f85xj.6.7 (statistical execution prerequisite).

use fs_blake3::hash_domain;
use fs_evidence::Color;
use fs_uq::{
    ParameterUncertainty, PropagationMethod, UqComplianceError, UqExecution, UqPlan,
    UqStatus,
};

fn plan(samples: usize) -> UqPlan {
    UqPlan::new("temperature_margin", PropagationMethod::MonteCarlo, samples)
        .with_parameter(ParameterUncertainty::uniform("offset", -1.0, 1.0, "K"))
        .with_compliance_threshold(0.0)
}

#[test]
fn assessment_uses_all_indicators_and_matches_the_closed_form_boundary() {
    let mut execution = UqExecution::new(&plan(192)).unwrap();
    let mut ordinal = 0;
    execution.advance(192, || false, |_| {
        let value = [-1.0, 0.0, 1.0][ordinal % 3];
        ordinal += 1;
        Ok::<_, &str>(value)
    });
    let estimate = execution.assess_compliance(0.05, 0.15).unwrap().unwrap();
    let t = 192.0_f64;
    let v = 0.25 * t + 1.0;
    let radius = (v * (v.ln() + 2.0 * (1.0_f64 / 0.05).ln())).sqrt() / t;
    assert_eq!(estimate.n, 192);
    assert_eq!(estimate.mean, 2.0 / 3.0); // Equality to the threshold is compliant.
    assert!((estimate.lo - (2.0 / 3.0 - radius)).abs() < 1e-12);
    assert!((estimate.hi - (2.0 / 3.0 + radius)).abs() < 1e-12);
    assert!(estimate.converged);
    assert!(!execution.assess_compliance(0.05, 0.1).unwrap().unwrap().converged);
    assert_eq!(execution.evaluations_attempted(), 192);
}

#[test]
fn no_observations_means_no_invented_confidence_estimate() {
    let mut execution = UqExecution::new(&plan(10)).unwrap();
    assert_eq!(execution.assess_compliance(0.05, 0.1), Ok(None));
    execution.advance(10, || true, |_| -> Result<f64, &str> { panic!("cancelled") });
    assert_eq!(execution.report().status, UqStatus::Cancelled);
    assert_eq!(execution.assess_compliance(0.05, 0.1), Ok(None));
}

#[test]
fn invalid_confidence_parameters_refuse_without_modifying_retained_work() {
    let mut execution = UqExecution::new(&plan(10)).unwrap();
    execution.advance(4, || false, |x| Ok::<_, &str>(x[0]));
    let model = hash_domain("test-model", b"identity");
    let before = execution.checkpoint(model).unwrap();
    for alpha in [0.0, -0.1, 1.0, 2.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert_eq!(execution.assess_compliance(alpha, 0.1), Err(UqComplianceError::InvalidAlpha));
    }
    for width in [-0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(execution.assess_compliance(0.05, width), Err(UqComplianceError::InvalidHalfWidth));
    }
    assert_eq!(execution.checkpoint(model).unwrap(), before);
    assert_eq!(execution.evaluations_attempted(), 4);
}

#[test]
fn missing_threshold_and_failed_sample_prefixes_cannot_produce_confidence() {
    let mut undeclared = plan(10);
    undeclared.compliance_threshold = None;
    let execution = UqExecution::new(&undeclared).unwrap();
    assert_eq!(execution.assess_compliance(0.05, 0.1), Err(UqComplianceError::MissingThreshold));
    let mut execution = UqExecution::new(&plan(10)).unwrap();
    execution.advance(4, || false, |x| Ok::<_, &str>(x[0]));
    execution.advance(1, || false, |_| Err::<f64, _>("solver failed"));
    assert_eq!(execution.observations().len(), 4);
    assert_eq!(execution.assess_compliance(0.05, 0.1), Err(UqComplianceError::RefusedExecution));
}

#[test]
fn constant_results_and_interval_clipping_do_not_fake_statistical_precision() {
    for output in [-1.0, 1.0] {
        let mut execution = UqExecution::new(&plan(2)).unwrap();
        let report = execution.advance(2, || false, |_| Ok::<_, &str>(output));
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(report.sampling_error, 0.0);
        let estimate = execution.assess_compliance(0.05, 0.5).unwrap().unwrap();
        assert_eq!([estimate.lo, estimate.hi], [0.0, 1.0]);
        // The clipped interval's half-width equals 0.5, but the raw radius is
        // wider: clipping must not be mistaken for meeting the precision goal.
        assert!(!estimate.converged);
        assert!(matches!(report.evidence_color, Color::Estimated { .. }));
    }
}

#[test]
fn chunked_confidence_stop_leaves_budget_and_physical_authority_unchanged() {
    let mut execution = UqExecution::new(&plan(10_000)).unwrap();
    let estimate = loop {
        let report = execution.advance(64, || false, |x| Ok::<_, &str>(x[0]));
        let estimate = execution.assess_compliance(0.05, 0.05).unwrap().unwrap();
        assert!(matches!(report.evidence_color, Color::Estimated { .. }));
        if estimate.converged { break estimate; }
        assert_ne!(report.status, UqStatus::Complete, "target should fit the budget");
    };
    assert!(estimate.n < 10_000);
    assert_eq!(execution.report().status, UqStatus::BudgetTruncated);
    assert_eq!(execution.plan().budget_max_samples, 10_000);
    assert_eq!(execution.evaluations_attempted() as u64, estimate.n);
    assert!(estimate.hi - estimate.lo <= 0.1 + 1e-14);
}

#[test]
fn checkpoint_resume_preserves_confidence_and_assessment_is_read_only() {
    let plan = plan(512);
    let model = hash_domain("test-model", b"offset");
    let mut execution = UqExecution::new(&plan).unwrap();
    execution.advance(173, || false, |x| Ok::<_, &str>(x[0]));
    let bytes = execution.checkpoint(model).unwrap();
    let expected = execution.assess_compliance(0.05, 0.1).unwrap();
    assert_eq!(execution.checkpoint(model).unwrap(), bytes);
    let mut restored = UqExecution::restore(&plan, model, &bytes).unwrap();
    assert_eq!(restored.assess_compliance(0.05, 0.1).unwrap(), expected);
    restored.advance(339, || false, |x| Ok::<_, &str>(x[0]));
    execution.advance(339, || false, |x| Ok::<_, &str>(x[0]));
    assert_eq!(restored.assess_compliance(0.05, 0.1), execution.assess_compliance(0.05, 0.1));
    assert_eq!(restored.checkpoint(model), execution.checkpoint(model));
}

#[test]
fn tiny_admitted_alpha_is_finite_and_zero_target_never_reports_convergence() {
    let mut execution = UqExecution::new(&plan(16)).unwrap();
    execution.advance(16, || false, |x| Ok::<_, &str>(x[0]));
    for alpha in [f64::MIN_POSITIVE, 0.05, 0.5] {
        let estimate = execution.assess_compliance(alpha, 0.0).unwrap().unwrap();
        assert!(estimate.lo.is_finite() && estimate.hi.is_finite());
        assert!((0.0..=1.0).contains(&estimate.lo));
        assert!((0.0..=1.0).contains(&estimate.hi));
        assert!(estimate.lo <= estimate.mean && estimate.mean <= estimate.hi);
        assert!(!estimate.converged);
    }
}
