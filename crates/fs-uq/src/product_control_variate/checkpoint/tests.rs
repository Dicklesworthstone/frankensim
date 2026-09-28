use super::*;
use crate::{CorrelationModel, ParameterUncertainty, PropagationMethod, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::MonteCarlo, 19)
        .with_parameter(ParameterUncertainty::uniform("power", 2.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 0.5, "K"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(314.0)
}
fn model() -> ContentHash { ContentHash([31; 32]) }
fn response(x: &[f64]) -> Result<f64, &'static str> {
    Ok(x[1] + 3.0 * x[0] + 0.25 * (x[0] - 4.0).powi(2))
}
fn execution() -> (UqExecution, LinearControlVariate) {
    let run = UqExecution::new(&plan()).unwrap();
    let control = run.freeze_linear_control_variate(&[3.0, 1.0]).unwrap();
    (run, control)
}
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|x| x.to_bits()).collect() }

#[test]
fn every_prefix_resumes_the_original_control_and_only_the_missing_physics() {
    let (mut whole, original) = execution();
    whole.advance(usize::MAX, || false, response);
    let expected = whole.assess_linear_control_variate(&original).unwrap();
    for prefix in [0, 1, 7, 18, 19] {
        let (mut partial, control) = execution();
        partial.advance(prefix, || false, response);
        let raw_before = partial.checkpoint(model()).unwrap();
        let bytes = control.checkpoint(&partial, model()).unwrap();
        assert_eq!(partial.checkpoint(model()).unwrap(), raw_before);
        let (mut recovered, frozen) = LinearControlVariate::restore(&plan(), model(), &bytes).unwrap();
        assert_eq!(bits(frozen.gradient()), bits(original.gradient()));
        assert_eq!(bits(frozen.parameter_means()), bits(original.parameter_means()));
        assert_eq!(recovered.checkpoint(model()).unwrap(), raw_before);
        let mut calls = 0;
        recovered.advance(usize::MAX, || false, |x| { calls += 1; response(x) });
        assert_eq!(calls, 19 - prefix);
        assert_eq!(recovered.assess_linear_control_variate(&frozen).unwrap(), expected);
        assert_eq!(recovered.report(), whole.report());
        assert_eq!(frozen.checkpoint(&recovered, model()).unwrap(), original.checkpoint(&whole, model()).unwrap());
    }
}

#[test]
fn interrupted_sample_retries_the_same_ordinal_and_completed_recovery_is_terminal() {
    let (mut run, control) = execution();
    run.advance(5, || false, response);
    let mut interrupted_parameters = Vec::new();
    let report = run.advance_interruptible(1, || false, |x| {
        interrupted_parameters = bits(x);
        Ok::<_, &str>(None)
    });
    assert_eq!(report.status, UqStatus::Cancelled);
    let bytes = control.checkpoint(&run, model()).unwrap();
    let (mut recovered, frozen) = LinearControlVariate::restore(&plan(), model(), &bytes).unwrap();
    assert_eq!(recovered.evaluations_attempted(), 5);
    recovered.advance(1, || false, |x| {
        assert_eq!(bits(x), interrupted_parameters);
        response(x)
    });
    recovered.advance(usize::MAX, || false, response);
    let bytes = frozen.checkpoint(&recovered, model()).unwrap();
    let (mut complete, _) = LinearControlVariate::restore(&plan(), model(), &bytes).unwrap();
    let report = complete.advance(usize::MAX, || false, |_| -> Result<f64, &str> {
        panic!("complete recovery must not call a model")
    });
    assert_eq!(report.status, UqStatus::Complete);
}

#[test]
fn changed_plan_model_and_spliced_coefficient_headers_refuse() {
    let (mut run, control) = execution();
    run.advance(6, || false, response);
    let bytes = control.checkpoint(&run, model()).unwrap();
    let mut different = plan(); different.seed += 1;
    assert_eq!(LinearControlVariate::restore(&different, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    assert_eq!(LinearControlVariate::restore(&plan(), ContentHash([32; 32]), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    assert_eq!(control.checkpoint(&UqExecution::new(&different).unwrap(), model()).unwrap_err(), UqCheckpointError::IdentityMismatch);
    let mut changed = bytes.clone();
    changed[16..24].copy_from_slice(&4.0_f64.to_bits().to_le_bytes());
    assert_eq!(LinearControlVariate::restore(&plan(), model(), &changed).unwrap_err(), UqCheckpointError::IdentityMismatch);
    let other = UqExecution::new(&plan()).unwrap().freeze_linear_control_variate(&[4.0, 1.0]).unwrap();
    let other_bytes = other.checkpoint(&run, model()).unwrap();
    let mut spliced = bytes.clone(); spliced[32..].copy_from_slice(&other_bytes[32..]);
    assert_eq!(LinearControlVariate::restore(&plan(), model(), &spliced).unwrap_err(), UqCheckpointError::IdentityMismatch);
    // The plain checkpoint cannot be reinterpreted as already controlled, or vice versa.
    assert!(LinearControlVariate::restore(&plan(), model(), &run.checkpoint(model()).unwrap()).is_err());
    assert!(UqExecution::restore(&plan(), model(), &bytes).is_err());
}

#[test]
fn malformed_payloads_and_failed_runs_cannot_publish_a_valid_prefix() {
    let (mut run, control) = execution();
    run.advance(3, || false, response);
    let bytes = control.checkpoint(&run, model()).unwrap();
    for end in 0..bytes.len() {
        assert!(LinearControlVariate::restore(&plan(), model(), &bytes[..end]).is_err());
    }
    let mut extended = bytes.clone(); extended.push(0);
    assert!(LinearControlVariate::restore(&plan(), model(), &extended).is_err());
    for index in [0, 8, 16, 24, 32, bytes.len() - 1] {
        let mut changed = bytes.clone(); changed[index] ^= 1;
        assert!(LinearControlVariate::restore(&plan(), model(), &changed).is_err());
    }
    for coefficient in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut changed = bytes.clone();
        changed[16..24].copy_from_slice(&coefficient.to_bits().to_le_bytes());
        assert!(LinearControlVariate::restore(&plan(), model(), &changed).is_err());
    }
    let mut oversized = bytes;
    oversized.resize(4096, 0);
    assert!(LinearControlVariate::restore(&plan(), model(), &oversized).is_err());
    run.advance(1, || false, |_| Err::<f64, _>("physical sample refused"));
    assert_eq!(control.checkpoint(&run, model()).unwrap_err(), UqCheckpointError::NotResumable);
}

#[test]
fn joint_gaussian_and_signed_zero_bits_survive_without_new_draws() {
    let p = UqPlan::new("response", PropagationMethod::MonteCarlo, 8)
        .with_parameter(ParameterUncertainty::gaussian("a", 2.0, 1.0, "1"))
        .with_parameter(ParameterUncertainty::gaussian("b", 7.0, 3.0, "1"))
        .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 0.8], vec![0.8, 1.0]] });
    let mut run = UqExecution::new(&p).unwrap();
    let control = run.freeze_linear_control_variate(&[-0.0, 2.0]).unwrap();
    run.advance(3, || false, |_| Ok::<_, &str>(-0.0));
    let bytes = control.checkpoint(&run, model()).unwrap();
    let (recovered, frozen) = LinearControlVariate::restore(&p, model(), &bytes).unwrap();
    assert_eq!(bits(frozen.gradient()), bits(control.gradient()));
    assert_eq!(bits(recovered.observations()), bits(run.observations()));
    assert_eq!(run.assess_linear_control_variate(&control).unwrap(), recovered.assess_linear_control_variate(&frozen).unwrap());
    let mut changed = p; changed.correlation = CorrelationModel::Independent;
    assert!(LinearControlVariate::restore(&changed, model(), &bytes).is_err());
}
