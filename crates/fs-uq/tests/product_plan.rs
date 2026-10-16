//! Tests for product uncertainty propagation and sampling execution.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`

use fs_evidence::Color;
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UqPlan, UqPropagator, UqStatus,
};

#[test]
fn test_uq_gaussian_sampling_and_compliance() {
    let plan = UqPlan::new("junction_maximum", PropagationMethod::MonteCarlo, 500)
        .with_parameter(ParameterUncertainty::gaussian(
            "ambient_temp",
            300.0,
            5.0,
            "K",
        ))
        .with_parameter(ParameterUncertainty::gaussian("die_power", 50.0, 2.0, "W"))
        .with_compliance_threshold(355.0)
        .with_correlation(CorrelationModel::Independent);

    // Linear thermal model: T_j = T_amb + R_th * P_die (with R_th = 0.8 K/W)
    let result = UqPropagator::run(&plan, |params| {
        let t_amb = params[0];
        let p_die = params[1];
        t_amb + 0.8 * p_die
    });

    assert_eq!(result.status, UqStatus::Complete);
    assert_eq!(result.samples_evaluated, 500);

    // Expected mean: 300 + 0.8 * 50 = 340.0 K
    let mean = result.mean.expect("mean computed");
    assert!(
        (mean - 340.0).abs() < 2.0,
        "mean should be near 340.0, got {mean}"
    );

    let p_comp = result
        .probability_of_compliance
        .expect("compliance probability computed");
    assert!(p_comp > 0.90, "P(T <= 355) should be > 0.90, got {p_comp}");

    assert!(matches!(result.evidence_color, Color::Estimated { .. }));
}

#[test]
fn test_uq_epistemic_interval_refuses_probability_measure() {
    let plan = UqPlan::new(
        "junction_maximum",
        PropagationMethod::EpistemicBounding,
        100,
    )
    .with_parameter(ParameterUncertainty::interval(
        "interface_conductance",
        1000.0,
        5000.0,
        "W/(m^2·K)",
    ))
    .with_compliance_threshold(350.0);

    let result = UqPropagator::run(&plan, |_| panic!("a band is not a sampling measure"));

    assert_eq!(result.status, UqStatus::Refused);
    assert_eq!(result.samples_evaluated, 0);
    // MUST NOT emit P(compliance) when input is an epistemic interval without a probability measure!
    assert!(
        result.probability_of_compliance.is_none(),
        "epistemic interval must refuse probability of compliance"
    );
}

#[test]
fn test_uq_unknown_correlation_refuses_multivariate_propagation() {
    let plan = UqPlan::new("junction_maximum", PropagationMethod::MonteCarlo, 100)
        .with_parameter(ParameterUncertainty::gaussian("p1", 10.0, 1.0, "W"))
        .with_parameter(ParameterUncertainty::gaussian("p2", 20.0, 2.0, "W"))
        .with_correlation(CorrelationModel::Unknown);

    let result = UqPropagator::run(&plan, |params| params[0] + params[1]);

    assert_eq!(result.status, UqStatus::Refused);
    assert!(result.rejection_reason.is_some());
}

#[test]
fn test_uq_determinism() {
    let plan1 = UqPlan::new("junction_maximum", PropagationMethod::MonteCarlo, 200)
        .with_parameter(ParameterUncertainty::gaussian(
            "ambient_temp",
            300.0,
            5.0,
            "K",
        ))
        .with_parameter(ParameterUncertainty::uniform(
            "fan_speed",
            2000.0,
            3000.0,
            "RPM",
        ))
        .with_correlation(CorrelationModel::Independent);

    let plan2 = plan1.clone();

    let res1 = UqPropagator::run(&plan1, |params| params[0] + 0.01 * params[1]);
    let res2 = UqPropagator::run(&plan2, |params| params[0] + 0.01 * params[1]);
    assert_eq!(res1.status, UqStatus::Complete);
    assert_eq!(res2.status, UqStatus::Complete);

    assert_eq!(
        res1.content_hash(),
        res2.content_hash(),
        "content hash bit-identical"
    );
}

fn gaussian_pair(rho: f64) -> UqPlan {
    UqPlan::new("sum", PropagationMethod::MonteCarlo, 20_000)
        .with_parameter(ParameterUncertainty::gaussian("a", 0.0, 1.0, "1"))
        .with_parameter(ParameterUncertainty::gaussian("b", 0.0, 1.0, "1"))
        .with_correlation(CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, rho], vec![rho, 1.0]],
        })
}

#[test]
fn g1_declared_gaussian_dependence_changes_the_actual_samples() {
    // Var(A+B) = 2+2*rho. Ignoring the matrix fails both nonzero-rho cases.
    for rho in [-0.75_f64, 0.0, 0.75] {
        let result = UqPropagator::run(&gaussian_pair(rho), |x| x[0] + x[1]);
        assert_eq!(result.status, UqStatus::Complete);
        let expected = (2.0 + 2.0 * rho).sqrt();
        assert!((result.std_dev.unwrap() / expected - 1.0).abs() < 0.03);
        assert!(result.mean.unwrap().abs() < 0.04);
        assert!(matches!(result.evidence_color, Color::Estimated { .. }));
    }
    // Singular PSD models are useful, and must retain their exact dependence.
    let result = UqPropagator::run(&gaussian_pair(-1.0), |x| x[0] + x[1]);
    assert_eq!(result.status, UqStatus::Complete);
    assert_eq!(result.interval_bounds, [0.0, 0.0]);
    assert_eq!(result.std_dev, Some(0.0));
}

#[test]
fn g0_invalid_or_unimplemented_plans_never_evaluate_the_model() {
    let base = gaussian_pair(0.5);
    let mut invalid = Vec::new();
    for method in [
        PropagationMethod::QuasiMonteCarlo,
        PropagationMethod::PolynomialChaos,
        PropagationMethod::MultilevelMonteCarlo,
        PropagationMethod::EpistemicBounding,
    ] {
        let mut plan = base.clone();
        plan.method = method;
        invalid.push(plan);
    }
    for count in [0, 1, 1_000_001] {
        let mut plan = base.clone();
        plan.budget_max_samples = count;
        invalid.push(plan);
    }
    invalid.push(base.clone().with_correlation(CorrelationModel::Unknown));
    invalid.push(base.clone().with_correlation(CorrelationModel::Correlated {
        matrix: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
    }));
    for matrix in [
        vec![vec![1.0]],
        vec![vec![1.0, f64::NAN], vec![f64::NAN, 1.0]],
        vec![vec![1.0, 0.5], vec![0.2, 1.0]],
        vec![vec![1.0, 1.1], vec![1.1, 1.0]],
    ] {
        invalid.push(
            base.clone()
                .with_correlation(CorrelationModel::JointGaussian { matrix }),
        );
    }
    let indefinite = base
        .clone()
        .with_parameter(ParameterUncertainty::gaussian("c", 0.0, 1.0, "1"))
        .with_correlation(CorrelationModel::JointGaussian {
            matrix: vec![
                vec![1.0, 0.9, 0.9],
                vec![0.9, 1.0, -0.9],
                vec![0.9, -0.9, 1.0],
            ],
        });
    invalid.push(indefinite);
    for parameter in [
        ParameterUncertainty::gaussian("a", 0.0, -1.0, "1"),
        ParameterUncertainty::gaussian("a", f64::INFINITY, 1.0, "1"),
        ParameterUncertainty::uniform("a", 2.0, 1.0, "1"),
        ParameterUncertainty::interval("a", 0.0, 1.0, "1"),
        ParameterUncertainty::uniform("a", 0.0, 1.0, "1"), // unspecified copula
    ] {
        let mut plan = base.clone();
        plan.parameters[0] = parameter;
        invalid.push(plan);
    }
    for plan in invalid {
        let result = UqPropagator::run(&plan, |_| panic!("invalid plan invoked model"));
        assert_eq!(result.status, UqStatus::Refused);
        assert_eq!(result.samples_evaluated, 0);
        assert!(result.probability_of_compliance.is_none());
        assert!(result.rejection_reason.is_some());
    }
    assert_eq!(
        UqPlan::new("q", PropagationMethod::MonteCarlo, 2).correlation,
        CorrelationModel::Unknown
    );
}

#[test]
fn g0_nonfinite_qoi_stops_without_publishing_statistics() {
    let calls = std::cell::Cell::new(0);
    let result = UqPropagator::run(&gaussian_pair(0.5), |_| {
        calls.set(calls.get() + 1);
        if calls.get() == 3 { f64::NAN } else { 1.0 }
    });
    assert_eq!(result.status, UqStatus::Refused);
    assert_eq!(calls.get(), 3);
    assert_eq!(result.samples_evaluated, 3);
    assert!(result.mean.is_none() && result.probability_of_compliance.is_none());
}

#[test]
fn g0_finite_constant_extremes_do_not_overflow_or_promote_authority() {
    let mut plan = gaussian_pair(0.0);
    plan.budget_max_samples = 100;
    for value in [0.0, f64::MAX, -f64::MAX] {
        let result = UqPropagator::run(&plan, |_| value);
        assert_eq!(result.status, UqStatus::Complete);
        assert_eq!(result.mean, Some(value));
        assert_eq!(result.std_dev, Some(0.0));
        assert!(matches!(result.evidence_color, Color::Estimated { .. }));
    }
}

#[test]
fn g0_result_identity_preserves_authority_and_small_range_changes() {
    let result = UqPropagator::run(&gaussian_pair(0.0), |x| x[0]);
    let mut changed = result.clone();
    changed.interval_bounds[0] += 1e-10;
    assert_ne!(result.content_hash(), changed.content_hash());
    changed = result.clone();
    changed.evidence_color = Color::Verified { lo: -1.0, hi: 1.0 };
    assert_ne!(result.content_hash(), changed.content_hash());
}

#[test]
fn g3_zero_uncertainty_collapses_without_inventing_dependence() {
    let plan = UqPlan::new("sum", PropagationMethod::MonteCarlo, 2)
        .with_parameter(ParameterUncertainty::gaussian("a", 3.0, 0.0, "1"))
        .with_parameter(ParameterUncertainty::uniform("b", 4.0, 4.0, "1"));
    let result = UqPropagator::run(&plan, |x| x[0] + x[1]);
    assert_eq!(result.status, UqStatus::Complete);
    assert_eq!(result.interval_bounds, [7.0, 7.0]);
    assert_eq!(result.sampling_error, 0.0);
}
