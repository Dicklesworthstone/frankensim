//! The compliance adapter reuses complete observations and the binary owner.
use fs_eproc::bernoulli::BernoulliMixtureCs;
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod,
    UqComplianceError, UqExecution, UqPlan, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("test-qoi", PropagationMethod::MonteCarlo, 2048)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("input", 0.0, 1.0, "1"))
        .with_compliance_threshold(0.5)
}

#[test]
fn assessment_matches_the_owner_and_does_not_change_work_statistics_or_checkpoint() {
    let p = plan();
    let model = fs_blake3::ContentHash([73; 32]);
    let mut execution = UqExecution::new(&p).unwrap();
    assert_eq!(execution.assess_bernoulli_compliance(0.05, 0.0).unwrap(), None);
    execution.advance(127, || false, |x| Ok::<_, &str>(x[0]));
    let before = execution.report();
    let checkpoint = execution.checkpoint(model).unwrap();
    let mut owner = BernoulliMixtureCs::new(0.05).unwrap();
    for &value in execution.observations() { owner.observe(value <= 0.5).unwrap(); }
    let expected = owner.interval().unwrap().unwrap();
    let observed = execution.assess_bernoulli_compliance(0.05, 0.0).unwrap().unwrap();
    assert_eq!((observed.mean, observed.lo, observed.hi, observed.n),
        (expected.mean, expected.lo, expected.hi, expected.n));
    assert!(!observed.converged);
    assert_eq!(execution.report(), before);
    assert_eq!(execution.checkpoint(model).unwrap(), checkpoint);
    let restored = UqExecution::restore(&p, model, &checkpoint).unwrap();
    assert_eq!(restored.assess_bernoulli_compliance(0.05, 0.0).unwrap(), Some(observed));
    let mut full = UqExecution::new(&p).unwrap();
    full.advance(2048, || false, |x| Ok::<_, &str>(x[0]));
    let mut resumed = restored;
    resumed.advance(2048, || false, |x| Ok::<_, &str>(x[0]));
    assert_eq!(resumed.assess_bernoulli_compliance(0.05, 0.0), full.assess_bernoulli_compliance(0.05, 0.0));
    assert_eq!(resumed.checkpoint(model), full.checkpoint(model));
}

#[test]
fn zero_failures_keep_nonzero_uncertainty_and_asymmetric_precision() {
    let mut execution = UqExecution::new(&plan()).unwrap();
    // Inclusive compliance at the exact threshold; no epsilon widening.
    execution.advance(1024, || false, |_| Ok::<_, &str>(0.5));
    let binary = execution.assess_bernoulli_compliance(0.05, 0.0).unwrap().unwrap();
    let generic = execution.assess_compliance(0.05, 0.0).unwrap().unwrap();
    assert_eq!(binary.mean, 1.0); assert_eq!(binary.hi, 1.0);
    assert!(binary.lo > 0.99 && binary.lo < 1.0);
    assert!(generic.lo < 0.99); assert!(!binary.converged);
    let width = binary.hi - binary.lo;
    assert!(!execution.assess_bernoulli_compliance(0.05, width*0.75).unwrap().unwrap().converged);
    assert!(execution.assess_bernoulli_compliance(0.05, width*1.01).unwrap().unwrap().converged);
}

#[test]
fn a_failed_sample_invalidates_inference_from_the_previously_successful_prefix() {
    for terminal in [Err("physical failure"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut execution = UqExecution::new(&plan()).unwrap();
        execution.advance(16, || false, |_| Ok::<_, &str>(0.0));
        assert!(execution.assess_bernoulli_compliance(0.05, 0.0).unwrap().is_some());
        assert_eq!(execution.advance(1, || false, |_| terminal).status, UqStatus::Refused);
        assert_eq!(execution.assess_bernoulli_compliance(0.05, 0.0), Err(UqComplianceError::RefusedExecution));
        assert_eq!(execution.evaluations_attempted(), 17);
    }
}

#[test]
fn inference_parameters_and_missing_events_refuse_without_mutation() {
    let execution = UqExecution::new(&plan()).unwrap();
    for alpha in [0.0, 1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert_eq!(execution.assess_bernoulli_compliance(alpha, 0.1), Err(UqComplianceError::InvalidAlpha));
    }
    for width in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(execution.assess_bernoulli_compliance(0.05, width), Err(UqComplianceError::InvalidHalfWidth));
    }
    let mut no_event = plan(); no_event.compliance_threshold = None;
    let no_event = UqExecution::new(&no_event).unwrap();
    assert_eq!(no_event.assess_bernoulli_compliance(0.05, 0.1), Err(UqComplianceError::MissingThreshold));
    assert_eq!(execution.evaluations_attempted(), 0);
}
