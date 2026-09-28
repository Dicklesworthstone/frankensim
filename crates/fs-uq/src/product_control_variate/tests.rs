use super::*;
use crate::{CorrelationModel, ParameterUncertainty, PropagationMethod};
use fs_blake3::ContentHash;

fn plan(n: usize) -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::MonteCarlo, n)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("power", 2.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("inlet", 1.0, 0.25, "K"))
        .with_compliance_threshold(310.0)
}
fn affine(x: &[f64]) -> Result<f64, &'static str> { Ok(300.0 + 3.0*x[0] - 2.0*x[1]) }
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }

#[test]
fn affine_adjoint_removes_variation_but_never_changes_compliance_data() {
    let mut run = UqExecution::new(&plan(256)).unwrap();
    let control = run.freeze_linear_control_variate(&[3.0,-2.0]).unwrap();
    assert_eq!(control.parameter_means(), &[4.0,1.0]);
    assert_eq!(control.gradient(), &[3.0,-2.0]);
    run.advance(256, || false, affine);
    let raw = run.observations().to_vec();
    let probability = run.assess_bernoulli_compliance(0.05,0.1).unwrap();
    let before = run.checkpoint(ContentHash([7;32])).unwrap();
    let estimate = run.assess_linear_control_variate(&control).unwrap().unwrap();
    // Centering at zero would give 300 K, not the exact declared expectation.
    close(estimate.mean,310.0,1e-11);
    assert!((estimate.mean-300.0).abs()>9.0);
    assert!(estimate.raw_std_dev.unwrap()>1.0);
    assert!(estimate.standard_error.unwrap()<1e-11);
    assert!(estimate.variance_ratio.unwrap()<1e-20);
    assert_eq!(estimate.raw_mean,run.report().mean.unwrap());
    assert_eq!(estimate.raw_std_dev,run.report().std_dev);
    assert_eq!(run.observations(),raw);
    assert_eq!(before,run.checkpoint(ContentHash([7;32])).unwrap());
    let after = run.assess_bernoulli_compliance(0.05,0.1).unwrap().unwrap();
    let probability = probability.unwrap();
    assert_eq!((after.mean,after.lo,after.hi,after.n),
        (probability.mean,probability.lo,probability.hi,probability.n));
}

#[test]
fn nonlinear_residual_is_evaluated_not_replaced_by_the_taylor_model() {
    let mut run = UqExecution::new(&plan(128)).unwrap();
    let control = run.freeze_linear_control_variate(&[3.0,-2.0]).unwrap();
    let mut independent = Vec::new();
    run.advance(128,||false,|x| {
        let remainder = 0.5*(x[0]-4.0).powi(2);
        independent.push(310.0+remainder);
        Ok::<_, &'static str>(300.0+3.0*x[0]-2.0*x[1]+remainder)
    });
    let report=run.assess_linear_control_variate(&control).unwrap().unwrap();
    close(report.mean,independent.iter().sum::<f64>()/128.0,1e-11);
    assert!(report.mean>310.1,"must retain the nonlinear model remainder");
    assert!(report.std_dev.unwrap()>0.1);
    assert!(report.variance_ratio.unwrap()<0.1);
}

#[test]
fn a_harmful_frozen_gradient_reports_increased_variance_not_a_fallback() {
    let p=UqPlan::new("response",PropagationMethod::MonteCarlo,64)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("x",-1.0,1.0,"1"));
    let mut run=UqExecution::new(&p).unwrap();
    let control=run.freeze_linear_control_variate(&[-2.0]).unwrap();
    run.advance(64,||false,|x|Ok::<_,&str>(2.0*x[0]));
    let report=run.assess_linear_control_variate(&control).unwrap().unwrap();
    close(report.variance_ratio.unwrap(),4.0,1e-12);
    close(report.mean,2.0*report.raw_mean,1e-14);
}

#[test]
fn frozen_control_survives_raw_checkpoint_recovery_and_interrupted_assessment() {
    let p=plan(40);
    let mut split=UqExecution::new(&p).unwrap();
    let control=split.freeze_linear_control_variate(&[3.0,-2.0]).unwrap();
    split.advance(13,||false,affine);
    let bytes=split.checkpoint(ContentHash([3;32])).unwrap();
    let mut polls=0;
    assert_eq!(split.assess_linear_control_variate_interruptible(&control,|| {
        polls+=1; polls==8
    }),Err(UqControlError::Cancelled));
    assert_eq!(split.checkpoint(ContentHash([3;32])).unwrap(),bytes);
    let mut restored=UqExecution::restore(&p,ContentHash([3;32]),&bytes).unwrap();
    assert_eq!(split.assess_linear_control_variate(&control).unwrap(),
        restored.assess_linear_control_variate(&control).unwrap());
    restored.advance(40,||false,affine);
    let mut whole=UqExecution::new(&p).unwrap();
    whole.advance(40,||false,affine);
    assert_eq!(restored.assess_linear_control_variate(&control).unwrap(),
        whole.assess_linear_control_variate(&control).unwrap());
    assert_eq!(restored.observations(),whole.observations());
}

#[test]
fn joint_gaussian_uses_its_declared_marginal_means_and_original_sampler() {
    let p=UqPlan::new("response",PropagationMethod::MonteCarlo,128)
        .with_correlation(CorrelationModel::JointGaussian {matrix:vec![vec![1.0,0.8],vec![0.8,1.0]]})
        .with_parameter(ParameterUncertainty::gaussian("a",2.0,1.0,"1"))
        .with_parameter(ParameterUncertainty::gaussian("b",7.0,3.0,"1"));
    let mut run=UqExecution::new(&p).unwrap();
    let control=run.freeze_linear_control_variate(&[-1.0,2.0]).unwrap();
    run.advance(128,||false,|x|Ok::<_,&str>(4.0-x[0]+2.0*x[1]));
    let report=run.assess_linear_control_variate(&control).unwrap().unwrap();
    close(report.mean,16.0,1e-12);
    assert!(report.standard_error.unwrap()<1e-12);
}

#[test]
fn malformed_late_mismatched_and_failed_assessments_refuse() {
    let mut run=UqExecution::new(&plan(8)).unwrap();
    for coefficients in [vec![],vec![1.0],vec![1.0,2.0,3.0],vec![f64::NAN,1.0],vec![1.0,f64::INFINITY]] {
        assert!(matches!(run.freeze_linear_control_variate(&coefficients),Err(UqControlError::InvalidGradient)));
    }
    let control=run.freeze_linear_control_variate(&[3.0,-2.0]).unwrap();
    let mut changed=plan(8); changed.parameters[0].unit="kW".into();
    assert_eq!(UqExecution::new(&changed).unwrap().assess_linear_control_variate(&control),Err(UqControlError::PlanMismatch));
    changed=plan(8); changed.seed+=1;
    assert_eq!(UqExecution::new(&changed).unwrap().assess_linear_control_variate(&control),Err(UqControlError::PlanMismatch));
    run.advance(1,||false,affine);
    assert!(matches!(run.freeze_linear_control_variate(&[3.0,-2.0]),Err(UqControlError::AlreadySampled)));
    run.advance(1,||false,|_|Err::<f64,_>("failed physical solve"));
    assert_eq!(run.assess_linear_control_variate(&control),Err(UqControlError::RefusedExecution));
}

#[test]
fn small_prefix_and_zero_dispersion_are_not_confidence_claims() {
    let mut run=UqExecution::new(&plan(4)).unwrap();
    let control=run.freeze_linear_control_variate(&[0.0,0.0]).unwrap();
    assert_eq!(run.assess_linear_control_variate(&control).unwrap(),None);
    run.advance(1,||false,|_|Ok::<_,&str>(300.0));
    let one=run.assess_linear_control_variate(&control).unwrap().unwrap();
    assert_eq!(one.n,1); assert_eq!(one.standard_error,None); assert_eq!(one.variance_ratio,None);
    run.advance(3,||false,|_|Ok::<_,&str>(300.0));
    let all=run.assess_linear_control_variate(&control).unwrap().unwrap();
    assert_eq!(all.mean,all.raw_mean); assert_eq!(all.standard_error,Some(0.0));
    assert_eq!(all.variance_ratio,None);
}

#[test]
fn overflowing_control_is_not_silently_dropped_and_raw_data_remains_usable() {
    let p=UqPlan::new("response",PropagationMethod::MonteCarlo,8)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("x",-1e300,1e300,"1"));
    let mut run=UqExecution::new(&p).unwrap();
    let control=run.freeze_linear_control_variate(&[1e300]).unwrap();
    let zero=run.freeze_linear_control_variate(&[0.0]).unwrap();
    run.advance(8,||false,|_|Ok::<_,&str>(1.0));
    let before=run.observations().to_vec();
    assert_eq!(run.assess_linear_control_variate(&control),Err(UqControlError::NumericalRange));
    assert_eq!(run.observations(),before);
    assert_eq!(run.assess_linear_control_variate(&zero).unwrap().unwrap().mean,1.0);
}
