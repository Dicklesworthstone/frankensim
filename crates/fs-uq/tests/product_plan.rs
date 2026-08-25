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
        .with_parameter(ParameterUncertainty::gaussian("ambient_temp", 300.0, 5.0, "K"))
        .with_parameter(ParameterUncertainty::gaussian("die_power", 50.0, 2.0, "W"))
        .with_compliance_threshold(355.0);

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
    assert!((mean - 340.0).abs() < 2.0, "mean should be near 340.0, got {mean}");

    let p_comp = result.probability_of_compliance.expect("compliance probability computed");
    assert!(p_comp > 0.90, "P(T <= 355) should be > 0.90, got {p_comp}");

    assert!(matches!(result.evidence_color, Color::Verified { .. }));
}

#[test]
fn test_uq_epistemic_interval_refuses_probability_measure() {
    let plan = UqPlan::new("junction_maximum", PropagationMethod::EpistemicBounding, 100)
        .with_parameter(ParameterUncertainty::interval("interface_conductance", 1000.0, 5000.0, "W/(m^2·K)"))
        .with_compliance_threshold(350.0);

    let result = UqPropagator::run(&plan, |params| {
        let htc = params[0];
        300.0 + 10000.0 / htc
    });

    assert_eq!(result.status, UqStatus::Complete);
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
        .with_parameter(ParameterUncertainty::gaussian("ambient_temp", 300.0, 5.0, "K"))
        .with_parameter(ParameterUncertainty::uniform("fan_speed", 2000.0, 3000.0, "RPM"));

    let plan2 = plan1.clone();

    let res1 = UqPropagator::run(&plan1, |params| params[0] + 0.01 * params[1]);
    let res2 = UqPropagator::run(&plan2, |params| params[0] + 0.01 * params[1]);

    assert_eq!(res1.content_hash(), res2.content_hash(), "content hash bit-identical");
}
