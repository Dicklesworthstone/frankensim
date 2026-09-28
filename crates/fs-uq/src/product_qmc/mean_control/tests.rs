use super::*;

fn config() -> QmcConfig { QmcConfig { replicates: 4, samples_per_replicate: 4 } }
fn plan() -> UqPlan {
    UqPlan::new("junction K", PropagationMethod::QuasiMonteCarlo, 16)
        .with_parameter(ParameterUncertainty::uniform("power", 2.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::uniform("inlet", 290.0, 310.0, "K"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(308.0)
}
fn matrix() -> Vec<Vec<f64>> { vec![vec![1.0, 0.75], vec![0.75, 1.0]] }
fn model() -> ContentHash { hash_domain("qmc-control-test", b"inlet+2*power") }
fn affine(x: &[f64]) -> Result<f64, &'static str> { Ok(x[1] + 2.0 * x[0]) }
fn close(a: f64, b: f64) { assert!((a - b).abs() < 1e-10, "{a:e} != {b:e}"); }

#[test]
fn affine_control_removes_variation_without_changing_raw_qmc_or_compliance() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    run.advance(16, || false, affine);
    let before = run.checkpoint(model()).unwrap();
    let report = run.report();
    let adjusted = run.assess_linear_control_variate(&control).unwrap();
    assert_eq!(adjusted.raw, report.estimate);
    assert_eq!(adjusted.raw_replicate_means, report.replicate_means);
    close(adjusted.controlled.as_ref().unwrap().mean, 308.0);
    assert!(adjusted.raw.as_ref().unwrap().standard_error.unwrap() > 1e-4);
    assert!(adjusted.controlled.as_ref().unwrap().standard_error.unwrap() < 1e-11);
    assert_eq!(run.report(), report);
    assert_eq!(run.checkpoint(model()).unwrap(), before);
}

#[test]
fn copula_controls_use_physical_units_and_means_not_latent_coordinates() {
    let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(), config()).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    assert_eq!(control.parameter_means(), [4.0, 300.0]);
    run.advance(16, || false, affine);
    let before = run.checkpoint(model()).unwrap();
    let adjusted = run.assess_linear_control_variate(&control).unwrap();
    close(adjusted.controlled.as_ref().unwrap().mean, 308.0);
    assert!(adjusted.variance_ratio.unwrap() < 1e-20);
    assert_eq!(adjusted.raw, run.report().estimate);
    assert_eq!(run.checkpoint(model()).unwrap(), before);
}

#[test]
fn nonlinear_standard_error_is_between_replicates_not_between_points() {
    let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(), config()).unwrap();
    let gradient = [2.0, 1.0];
    let control = run.freeze_linear_control_variate(&gradient).unwrap();
    let mut independently_adjusted = Vec::new();
    run.advance(16, || false, |x| {
        let y = x[1] + x[0] * x[0];
        independently_adjusted.push(y - 2.0 * (x[0] - 4.0) - (x[1] - 300.0));
        Ok::<_, &'static str>(y)
    });
    let estimate = run.assess_linear_control_variate(&control).unwrap();
    let means: Vec<f64> = independently_adjusted.chunks_exact(4)
        .map(|chunk| chunk.iter().sum::<f64>() / 4.0).collect();
    for (a, b) in estimate.controlled_replicate_means.iter().zip(&means) { close(*a, *b); }
    let mean = means.iter().sum::<f64>() / 4.0;
    let error = (means.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 12.0).sqrt();
    let wrong_point_error = (independently_adjusted.iter().map(|x| (x - mean).powi(2))
        .sum::<f64>() / (16.0 * 15.0)).sqrt();
    let actual = estimate.controlled.unwrap();
    close(actual.mean, mean);
    close(actual.standard_error.unwrap(), error);
    assert!((error - wrong_point_error).abs() > 1e-4);
}

#[test]
fn incomplete_nets_remain_paid_work_without_entering_estimates() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let control = run.freeze_linear_control_variate(&[1.0, 0.5]).unwrap();
    let empty = run.assess_linear_control_variate(&control).unwrap();
    assert_eq!(empty.samples_accepted, 0);
    assert!(empty.controlled.is_none());
    run.advance(3, || false, affine);
    let partial = run.assess_linear_control_variate(&control).unwrap();
    assert_eq!((partial.samples_accepted, partial.samples_in_estimate), (3, 0));
    assert!(partial.raw.is_none());
    run.advance(2, || false, affine);
    let one = run.assess_linear_control_variate(&control).unwrap();
    assert_eq!((one.samples_accepted, one.samples_in_estimate, one.completed_replicates), (5, 4, 1));
    assert_eq!(one.controlled.unwrap().standard_error, None);
    run.advance(3, || false, affine);
    assert!(run.assess_linear_control_variate(&control).unwrap().controlled.unwrap().standard_error.is_some());
}

#[test]
fn all_plain_and_copula_prefixes_restore_the_same_controls_and_raw_results() {
    let p = plan(); let m = matrix(); let cfg = config();
    let mut full = QmcExecution::new(&p, cfg).unwrap();
    let c = full.freeze_linear_control_variate(&[1.5, 0.75]).unwrap();
    full.advance(16, || false, affine);
    let mut full_copula = GaussianCopulaQmcExecution::new(&p, &m, cfg).unwrap();
    let cc = full_copula.freeze_linear_control_variate(&[1.5, 0.75]).unwrap();
    full_copula.advance(16, || false, affine);
    for split in 0..=16 {
        let mut part = QmcExecution::new(&p, cfg).unwrap();
        let cpart = part.freeze_linear_control_variate(&[1.5, 0.75]).unwrap();
        part.advance(split, || false, affine);
        let (mut resumed, retained) = QmcLinearControlVariate::restore(
            &p, cfg, model(), &cpart.checkpoint(&part, model()).unwrap()).unwrap();
        resumed.advance(16, || false, affine);
        assert_eq!(resumed.assess_linear_control_variate(&retained).unwrap(), full.assess_linear_control_variate(&c).unwrap());
        assert_eq!(resumed.checkpoint(model()).unwrap(), full.checkpoint(model()).unwrap());
        let mut part = GaussianCopulaQmcExecution::new(&p, &m, cfg).unwrap();
        let cpart = part.freeze_linear_control_variate(&[1.5, 0.75]).unwrap();
        part.advance(split, || false, affine);
        let (mut resumed, retained) = QmcLinearControlVariate::restore_copula(
            &p, &m, cfg, model(), &cpart.checkpoint_copula(&part, model()).unwrap()).unwrap();
        resumed.advance(16, || false, affine);
        assert_eq!(retained.gradient(), [1.5, 0.75]);
        assert_eq!(resumed.assess_linear_control_variate(&retained).unwrap(), full_copula.assess_linear_control_variate(&cc).unwrap());
        assert_eq!(resumed.checkpoint(model()).unwrap(), full_copula.checkpoint(model()).unwrap());
    }
}

#[test]
fn layout_dependence_physical_units_and_coefficients_cannot_be_substituted() {
    let p = plan(); let m = matrix(); let cfg = config();
    let run = GaussianCopulaQmcExecution::new(&p, &m, cfg).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 1.0]).unwrap();
    let saved = control.checkpoint_copula(&run, model()).unwrap();
    for change in 0..5 {
        let mut p = p.clone(); let mut m = m.clone(); let mut cfg = cfg;
        match change {
            0 => p.seed += 1,
            1 => cfg = QmcConfig { replicates: 2, samples_per_replicate: 8 },
            2 => p.parameters[0].unit = "kW".into(),
            3 => p.parameters[0].kind = crate::UncertaintyKind::AleatoryUniform { lo: 1.0, hi: 7.0 },
            _ => { m[0][1] = -0.75; m[1][0] = -0.75; }
        }
        let different = GaussianCopulaQmcExecution::new(&p, &m, cfg).unwrap();
        assert_eq!(different.assess_linear_control_variate(&control), Err(UqControlError::PlanMismatch));
        assert!(QmcLinearControlVariate::restore_copula(&p, &m, cfg, model(), &saved).is_err());
    }
    let plain = QmcExecution::new(&p, cfg).unwrap();
    assert_eq!(plain.assess_linear_control_variate(&control), Err(UqControlError::PlanMismatch));
    assert!(control.checkpoint(&plain, model()).is_err());
    assert!(QmcLinearControlVariate::restore(&p, cfg, model(), &saved).is_err());
    let mut changed = saved.clone(); changed[16..24].copy_from_slice(&3.0_f64.to_le_bytes());
    assert!(QmcLinearControlVariate::restore_copula(&p, &m, cfg, model(), &changed).is_err());
}

#[test]
fn unfavorable_coefficients_are_reported_not_dropped_for_a_better_result() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let control = run.freeze_linear_control_variate(&[-2.0, -1.0]).unwrap();
    run.advance(16, || false, affine);
    close(run.assess_linear_control_variate(&control).unwrap().variance_ratio.unwrap(), 4.0);
}

#[test]
fn cancelled_assessment_leaves_raw_execution_intact_and_failed_prefix_refuses() {
    let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(), config()).unwrap();
    assert!(matches!(run.freeze_linear_control_variate(&[f64::NAN, 1.0]), Err(UqControlError::InvalidGradient)));
    let control = run.freeze_linear_control_variate(&[1.5, 0.75]).unwrap();
    run.advance(8, || false, affine);
    assert!(matches!(run.freeze_linear_control_variate(&[2.0, 1.0]), Err(UqControlError::AlreadySampled)));
    let before = control.checkpoint_copula(&run, model()).unwrap();
    let mut polls = 0;
    run.assess_linear_control_variate_interruptible(&control, || { polls += 1; false }).unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert_eq!(run.assess_linear_control_variate_interruptible(&control, || {
            count += 1; count == stop
        }), Err(UqControlError::Cancelled));
        assert_eq!(control.checkpoint_copula(&run, model()).unwrap(), before);
    }
    run.advance(1, || false, |_| Err::<f64, _>("physical failure"));
    assert_eq!(run.assess_linear_control_variate(&control), Err(UqControlError::RefusedExecution));
    assert!(control.checkpoint_copula(&run, model()).is_err());
}

#[test]
fn extreme_constants_and_unused_coefficients_remain_finite() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let zero = run.freeze_linear_control_variate(&[0.0, 0.0]).unwrap();
    run.advance(16, || false, |_| Ok::<_, &str>(f64::MAX));
    let estimate = run.assess_linear_control_variate(&zero).unwrap();
    assert_eq!(estimate.raw, estimate.controlled);
    assert_eq!(estimate.controlled.unwrap().mean, f64::MAX);
    assert_eq!(estimate.variance_ratio, None);
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let huge = run.freeze_linear_control_variate(&[f64::MAX, f64::MAX]).unwrap();
    run.advance(16, || false, affine);
    let before = run.checkpoint(model()).unwrap();
    assert_eq!(run.assess_linear_control_variate(&huge), Err(UqControlError::NumericalRange));
    assert_eq!(run.checkpoint(model()).unwrap(), before);
}
