use super::*;
use crate::{ParameterUncertainty, PropagationMethod};

fn plan() -> UqPlan {
    UqPlan::new("junction K", PropagationMethod::MonteCarlo, 31)
        .with_parameter(ParameterUncertainty::uniform("power", 2.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::uniform("inlet", 290.0, 310.0, "K"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(308.0)
}
fn matrix() -> Vec<Vec<f64>> { vec![vec![1.0, 0.75], vec![0.75, 1.0]] }
fn model() -> ContentHash { hash_domain("test-physical-copula-control", b"T=inlet+2*power") }
fn evaluate(x: &[f64]) -> Result<f64, &'static str> { Ok(x[1] + 2.0 * x[0]) }

#[test]
fn physical_adjoint_removes_affine_variance_not_latent_normal_variance() {
    let mut run = GaussianCopulaExecution::new(&plan(), &matrix()).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    assert_eq!(control.parameter_means(), [4.0, 300.0]);
    assert_eq!(run.monte_carlo().parameter_means(), [0.0, 0.0]);
    run.advance(31, || false, evaluate);
    let before = run.checkpoint(model()).unwrap();
    let raw = run.report();
    let estimate = run.assess_linear_control_variate(&control).unwrap().unwrap();
    assert_eq!(estimate.n, 31);
    assert!((estimate.mean - 308.0).abs() < 1e-11);
    assert!(estimate.raw_standard_error.unwrap() > 0.1);
    assert!(estimate.standard_error.unwrap() < 1e-11);
    assert!(estimate.variance_ratio.unwrap() < 1e-20);
    assert_eq!(run.checkpoint(model()).unwrap(), before);
    assert_eq!(run.report(), raw, "physical quantiles and pass indicators must stay raw");
}

#[test]
fn fixed_bad_control_increases_variance_and_is_not_silently_discarded() {
    let mut run = GaussianCopulaExecution::new(&plan(), &matrix()).unwrap();
    let bad = run.freeze_linear_control_variate(&[-2.0, -1.0]).unwrap();
    run.advance(31, || false, evaluate);
    let estimate = run.assess_linear_control_variate(&bad).unwrap().unwrap();
    assert!((estimate.variance_ratio.unwrap() - 4.0).abs() < 1e-11);
}

#[test]
fn every_prefix_restores_coefficients_raw_observations_and_physical_means() {
    let p = plan(); let m = matrix();
    let mut full = GaussianCopulaExecution::new(&p, &m).unwrap();
    let control = full.freeze_linear_control_variate(&[2.0, 0.75]).unwrap();
    full.advance(31, || false, evaluate);
    let expected = full.assess_linear_control_variate(&control).unwrap();
    for split in 0..=31 {
        let mut prefix = GaussianCopulaExecution::new(&p, &m).unwrap();
        let control = prefix.freeze_linear_control_variate(&[2.0, 0.75]).unwrap();
        prefix.advance(split, || false, evaluate);
        let bytes = control.checkpoint(&prefix, model()).unwrap();
        let (mut resumed, retained) = CopulaLinearControlVariate::restore(&p, &m, model(), &bytes).unwrap();
        assert_eq!(retained.gradient(), [2.0, 0.75]);
        assert_eq!(retained.parameter_means(), [4.0, 300.0]);
        resumed.advance(31, || false, evaluate);
        assert_eq!(resumed.assess_linear_control_variate(&retained).unwrap(), expected);
        assert_eq!(resumed.checkpoint(model()).unwrap(), full.checkpoint(model()).unwrap());
    }
}

#[test]
fn complete_law_and_coefficient_header_remain_bound() {
    let p = plan(); let m = matrix();
    let run = GaussianCopulaExecution::new(&p, &m).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    let bytes = control.checkpoint(&run, model()).unwrap();
    for change in 0..5 {
        let mut other = p.clone(); let mut matrix = m.clone();
        match change {
            0 => other.seed += 1,
            1 => other.budget_max_samples += 1,
            2 => other.parameters[0].unit = "kW".into(),
            3 => other.parameters[0].kind = crate::UncertaintyKind::AleatoryUniform { lo: 1.0, hi: 7.0 },
            _ => { matrix[0][1] = -0.75; matrix[1][0] = -0.75; }
        }
        let other_run = GaussianCopulaExecution::new(&other, &matrix).unwrap();
        assert_eq!(other_run.assess_linear_control_variate(&control), Err(UqControlError::PlanMismatch));
        assert!(control.checkpoint(&other_run, model()).is_err());
        assert!(CopulaLinearControlVariate::restore(&other, &matrix, model(), &bytes).is_err());
    }
    let mut changed = bytes.clone(); changed[16..24].copy_from_slice(&3.0_f64.to_le_bytes());
    assert!(CopulaLinearControlVariate::restore(&p, &m, model(), &changed).is_err());
    assert!(CopulaLinearControlVariate::restore(&p, &m, hash_domain("other", b"model"), &bytes).is_err());
    for length in [0, 8, 15, 16, bytes.len() - 1] {
        assert!(CopulaLinearControlVariate::restore(&p, &m, model(), &bytes[..length]).is_err());
    }
}

#[test]
fn no_late_freeze_no_refused_prefix_estimate_and_no_partial_cancellation() {
    let mut run = GaussianCopulaExecution::new(&plan(), &matrix()).unwrap();
    for gradient in [vec![], vec![f64::NAN, 1.0], vec![1.0, f64::INFINITY]] {
        assert!(matches!(run.freeze_linear_control_variate(&gradient), Err(UqControlError::InvalidGradient)));
    }
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    assert_eq!(run.assess_linear_control_variate(&control), Ok(None));
    run.advance(3, || false, evaluate);
    assert!(matches!(run.freeze_linear_control_variate(&[2.0, 1.0]), Err(UqControlError::AlreadySampled)));
    let before = control.checkpoint(&run, model()).unwrap();
    let mut polls = 0;
    run.assess_linear_control_variate_interruptible(&control, || { polls += 1; false }).unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert_eq!(run.assess_linear_control_variate_interruptible(&control, || {
            count += 1; count == stop
        }), Err(UqControlError::Cancelled));
        assert_eq!(control.checkpoint(&run, model()).unwrap(), before);
    }
    run.advance(1, || false, |_| Err::<f64, _>("actual model refusal"));
    assert_eq!(run.assess_linear_control_variate(&control), Err(UqControlError::RefusedExecution));
    assert!(matches!(run.freeze_linear_control_variate(&[2.0, 1.0]), Err(UqControlError::RefusedExecution)));
    assert!(control.checkpoint(&run, model()).is_err());
}

#[test]
fn interruption_retries_the_same_physical_sample_with_the_same_control() {
    let p = plan(); let m = matrix();
    let mut run = GaussianCopulaExecution::new(&p, &m).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    run.advance(2, || false, evaluate);
    let mut interrupted = Vec::new();
    run.advance_interruptible(1, || false, |x| { interrupted = x.to_vec(); Ok::<_, &str>(None) });
    let (mut resumed, control) = CopulaLinearControlVariate::restore(
        &p, &m, model(), &control.checkpoint(&run, model()).unwrap()).unwrap();
    resumed.advance(1, || false, |x| { assert_eq!(x, interrupted); evaluate(x) });
    assert_eq!(resumed.monte_carlo().evaluations_attempted(), 3);
    assert!((resumed.assess_linear_control_variate(&control).unwrap().unwrap().mean - 308.0).abs() < 1e-11);
}

#[test]
fn unrepresentable_control_refuses_without_poisoning_raw_execution() {
    let mut run = GaussianCopulaExecution::new(&plan(), &matrix()).unwrap();
    let control = run.freeze_linear_control_variate(&[0.0, f64::MAX]).unwrap();
    run.advance(31, || false, evaluate);
    let before = run.checkpoint(model()).unwrap();
    assert_eq!(run.assess_linear_control_variate(&control), Err(UqControlError::NumericalRange));
    assert_eq!(run.checkpoint(model()).unwrap(), before);
}
