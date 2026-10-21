use fs_uq::{MultiQoiAssessmentError as Error, MultiQoiExecution, ParameterUncertainty,
    PropagationMethod, QoiDefinition, UqComplianceError, UqExecution, UqPlan, UqStatus};

fn plan(n: usize) -> UqPlan {
    UqPlan::new("paired assessment", PropagationMethod::MonteCarlo, n)
        .with_parameter(ParameterUncertainty::uniform("input", -1.0, 1.0, "1"))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())), "{a:e} != {b:e}");
}
fn variance(x: &[f64]) -> f64 {
    let center = x.iter().sum::<f64>() / x.len() as f64;
    x.iter().map(|v| (v-center).powi(2)).sum::<f64>() / (x.len()-1) as f64
}
fn definitions() -> Vec<QoiDefinition> {
    vec![QoiDefinition::new("x", "N").with_upper_limit(0.8),
         QoiDefinition::new("opposite", "Pa").with_upper_limit(0.5),
         QoiDefinition::new("unconstrained", "K")]
}
fn response(x: &[f64]) -> Result<Vec<f64>, &'static str> { Ok(vec![x[0], -x[0], 300.0+x[0]*x[0]]) }

#[test]
fn paired_linear_covariance_retains_signs_units_and_undefined_constant_correlations() {
    let definitions = [QoiDefinition::new("x", "N"), QoiDefinition::new("twice", "Pa"),
        QoiDefinition::new("negative", "K"), QoiDefinition::new("constant", "W")];
    let mut run = MultiQoiExecution::new(&plan(64), &definitions).unwrap();
    run.advance(64, || false, |x| Ok::<_, &str>(vec![x[0], 3.0+2.0*x[0], -x[0], 7.0]));
    let original = run.report();
    let result = run.assess_covariance(64*10, || false).unwrap().unwrap();
    let v = variance(run.observations(0).unwrap());
    assert_eq!(result.quantities, definitions);
    for (i, a) in [1.0, 2.0, -1.0, 0.0].into_iter().enumerate() {
        for (j, b) in [1.0, 2.0, -1.0, 0.0].into_iter().enumerate() {
            close(result.covariance[i][j].unwrap(), a*b*v, 1e-13);
            close(result.mean_covariance[i][j].unwrap(), a*b*v/64.0, 1e-13);
            assert_eq!(result.covariance[i][j], result.covariance[j][i]);
            if a*b == 0.0 { assert!(result.correlation[i][j].is_none()); }
            else { close(result.correlation[i][j].unwrap(), (a*b).signum(), 1e-13); }
        }
    }
    assert_eq!(run.report(), original);
    assert_eq!(run.assess_covariance(64*10-1, || false), Err(Error::Budget));
}

#[test]
fn centering_preserves_small_variations_on_large_offsets() {
    let mut run = MultiQoiExecution::new(&plan(128),
        &[QoiDefinition::new("a", "N"), QoiDefinition::new("b", "Pa")]).unwrap();
    run.advance(128, || false, |x| {
        let v = if x[0] < 0.0 { -2.0 } else { 2.0 };
        Ok::<_, &str>(vec![1e16+v, -1e16+3.0*v])
    });
    let changes: Vec<_> = run.observations(0).unwrap().iter().map(|x| x-1e16).collect();
    let v = variance(&changes);
    assert!(v > 0.0);
    let result = run.assess_covariance(384, || false).unwrap().unwrap();
    close(result.covariance[0][0].unwrap(), v, 1e-12);
    close(result.covariance[0][1].unwrap(), 3.0*v, 1e-12);
    close(result.covariance[1][1].unwrap(), 9.0*v, 1e-12);
    close(result.correlation[0][1].unwrap(), 1.0, 1e-14);
}

#[test]
fn extreme_output_units_keep_representable_cross_moments_and_mean_covariance() {
    let definitions = [QoiDefinition::new("large", "N"), QoiDefinition::new("small", "Pa")];
    let mut run = MultiQoiExecution::new(&plan(128), &definitions).unwrap();
    let mut inputs = Vec::new();
    let report = run.advance(128, || false, |x| {
        inputs.push(x[0]); Ok::<_, &str>(vec![x[0]*1e200, x[0]*1e-200])
    });
    assert_eq!(report.status, UqStatus::Complete);
    let result = run.assess_covariance(384, || false).unwrap().unwrap();
    assert!(result.covariance[0][0].is_none()); // overflow is not zero
    assert!(result.covariance[1][1].is_none()); // nonzero underflow is not zero
    close(result.covariance[0][1].unwrap(), variance(&inputs), 1e-12);
    close(result.correlation[0][1].unwrap(), 1.0, 1e-13);

    let mut run = MultiQoiExecution::new(&plan(128), &[QoiDefinition::new("wide", "N")]).unwrap();
    run.advance(128, || false, |x| Ok::<_, &str>(vec![if x[0]<0.0 { -2e154 } else { 2e154 }]));
    let result = run.assess_covariance(128, || false).unwrap().unwrap();
    assert!(result.covariance[0][0].is_none());
    assert!(result.mean_covariance[0][0].unwrap() > 1e306);
    assert!(result.mean_covariance[0][0].unwrap().is_finite());
}

#[test]
fn joint_confidence_matches_the_canonical_scalar_intersection_event() {
    let p = plan(128);
    let mut run = MultiQoiExecution::new(&p, &definitions()).unwrap();
    let mut scalar = UqExecution::new(&p.clone().with_compliance_threshold(0.0)).unwrap();
    assert!(run.assess_joint_compliance(0.05, 0.1, || false).unwrap().is_none());
    for chunk in [1, 7, 24, 96] {
        run.advance(chunk, || false, response);
        scalar.advance(chunk, || false, |x| {
            Ok::<_, &str>(if x[0] <= 0.8 && -x[0] <= 0.5 { 0.0 } else { 1.0 })
        });
        assert_eq!(run.assess_joint_compliance(0.05, 0.1, || false).unwrap(),
            scalar.assess_bernoulli_compliance(0.05, 0.1).unwrap());
    }
    let before = run.report();
    let result = run.assess_joint_compliance(0.05, 0.0, || false).unwrap().unwrap();
    close(result.mean, before.joint_compliance.unwrap().mean, 1e-15);
    assert!(!result.converged);
    assert_eq!(result.n, 128);
    assert_eq!(run.report().samples_evaluated, 128);
}

#[test]
fn no_observed_joint_failures_are_not_a_probability_one_certificate() {
    for constant in [0.0, 1.0] {
        let mut run = MultiQoiExecution::new(&plan(32),
            &[QoiDefinition::new("a", "N").with_upper_limit(0.0),
              QoiDefinition::new("b", "Pa").with_upper_limit(0.0)]).unwrap();
        run.advance(32, || false, |_| Ok::<_, &str>(vec![constant, 0.0]));
        let result = run.assess_joint_compliance(0.05, 0.0, || false).unwrap().unwrap();
        assert!(!result.converged);
        if constant == 0.0 { assert_eq!(result.hi, 1.0); assert!(result.lo < 1.0); }
        else { assert_eq!(result.lo, 0.0); assert!(result.hi > 0.0); }
    }
}

#[test]
fn covariance_and_confidence_cancel_without_mutation_at_every_checkpoint() {
    let mut run = MultiQoiExecution::new(&plan(16), &definitions()).unwrap();
    run.advance(16, || false, response);
    let before = run.report();
    let rows: Vec<_> = (0..16).map(|i| run.observation(i).unwrap()).collect();
    let mut covariance_polls = 0;
    let expected = run.assess_covariance(96, || { covariance_polls += 1; false }).unwrap();
    for stop in 0..covariance_polls {
        let mut polls = 0;
        assert_eq!(run.assess_covariance(96, || { polls += 1; polls > stop }), Err(Error::Cancelled));
    }
    assert_eq!(run.assess_covariance(96, || false).unwrap(), expected);
    let mut confidence_polls = 0;
    let expected = run.assess_joint_compliance(0.05, 0.1,
        || { confidence_polls += 1; false }).unwrap();
    for stop in 0..confidence_polls {
        let mut polls = 0;
        assert_eq!(run.assess_joint_compliance(0.05, 0.1,
            || { polls += 1; polls > stop }), Err(Error::Cancelled));
    }
    assert_eq!(run.assess_joint_compliance(0.05, 0.1, || false).unwrap(), expected);
    assert_eq!(run.report(), before);
    assert_eq!((0..16).map(|i| run.observation(i).unwrap()).collect::<Vec<_>>(), rows);
}

#[test]
fn empty_singleton_failed_and_invalid_assessments_do_not_manufacture_statistics() {
    let mut run = MultiQoiExecution::new(&plan(8), &definitions()).unwrap();
    assert!(run.assess_covariance(0, || false).unwrap().is_none());
    run.advance(1, || false, response);
    assert!(run.assess_covariance(0, || false).unwrap().is_none());
    for alpha in [0.0, 1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(run.assess_joint_compliance(alpha, 0.1, || false),
            Err(Error::Compliance(UqComplianceError::InvalidAlpha)));
    }
    for radius in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(run.assess_joint_compliance(0.05, radius, || false),
            Err(Error::Compliance(UqComplianceError::InvalidHalfWidth)));
    }
    run.advance(1, || false, |_| Err::<Vec<f64>, _>("physics refused"));
    assert_eq!(run.assess_covariance(usize::MAX, || false), Err(Error::RefusedExecution));
    assert_eq!(run.assess_joint_compliance(0.05, 0.1, || false), Err(Error::RefusedExecution));
    let run = MultiQoiExecution::new(&plan(8), &[QoiDefinition::new("free", "K")]).unwrap();
    assert_eq!(run.assess_joint_compliance(0.05, 0.1, || false),
        Err(Error::Compliance(UqComplianceError::MissingThreshold)));
}
