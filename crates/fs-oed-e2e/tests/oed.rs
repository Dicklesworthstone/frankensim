//! End-to-end battery: sensors are placed on the decision-relevant candidates,
//! the campaign stops when the decision is robust, and the posterior remains
//! honestly model-form Estimated.

use fs_evidence::{Color, color_leaf_identity_reason};
use fs_oed_e2e::{
    Candidate, CandidateError, MAX_CAMPAIGN_CANDIDATES, MAX_CAMPAIGN_EVALUATIONS,
    MAX_CAMPAIGN_SENSORS, OedError, demo_candidates, run_campaign,
};

fn candidate(
    name: &str,
    truth: f64,
    prior_mean: f64,
    prior_var: f64,
    sensor_noise: f64,
    sensor_cost: f64,
) -> Candidate {
    Candidate::new(
        name,
        truth,
        prior_mean,
        prior_var,
        sensor_noise,
        sensor_cost,
    )
    .expect("test candidate must satisfy the checked constructor")
}

fn estimator(color: &Color) -> &str {
    match color {
        Color::Estimated { estimator, .. } => estimator,
        stronger => panic!("expected Estimated evidence, got {stronger:?}"),
    }
}

#[test]
fn sensors_target_the_decision_and_the_campaign_knows_when_to_stop() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let report = run_campaign(&candidates, 0.01, 12).expect("demo campaign succeeds");
    // sensors WERE placed, and on the decision-relevant contenders (A and/or B),
    // never on the clearly-dominated D.
    assert!(report.sensors_placed > 0, "no sensors placed");
    assert!(report.placements.iter().any(|n| n == "A" || n == "B"));
    assert!(
        !report.placements.contains(&"D".to_string()),
        "wasted a sensor on D"
    );
    // measurement sharpened the beliefs; EVPI fell.
    assert!(report.variance_reduction > 0.0);
    assert!(report.final_evpi < report.initial_evpi, "EVPI did not fall");
    // the campaign STOPPED because the decision became robust.
    assert!(report.decision_robust, "did not reach a robust decision");
    assert!(
        report.final_evpi <= 0.01 + 1e-9,
        "final EVPI {}",
        report.final_evpi
    );
    // A is the true best and is chosen.
    assert_eq!(report.chosen_design, "A");
    // Kalman variance remains model-form Estimated until independently
    // certified; EVPI is Estimated as well.
    assert!(matches!(report.variance_color, Color::Estimated { .. }));
    assert!(matches!(report.evpi_color, Color::Estimated { .. }));
    assert_eq!(report.evpi_trace.len(), report.sensors_placed + 1);
    assert_eq!(report.assimilation_colors.len(), report.sensors_placed);
    for color in &report.assimilation_colors {
        let identity = estimator(color);
        assert!(color_leaf_identity_reason(identity).is_none());
    }
    // a cost-optimal precision budget was allocated across candidates.
    assert_eq!(report.allocation.len(), 4);
    assert!(report.allocation.iter().all(|(_, t)| *t > 0.0));
    println!(
        "{{\"campaign\":\"sensorforge\",\"placements\":{:?},\"sensors\":{},\"var_reduction\":{:.3},\
         \"initial_evpi\":{:.4},\"final_evpi\":{:.4},\"robust\":{},\"chosen\":\"{}\"}}",
        report.placements,
        report.sensors_placed,
        report.variance_reduction,
        report.initial_evpi,
        report.final_evpi,
        report.decision_robust,
        report.chosen_design,
    );
}

#[test]
fn a_clear_winner_needs_no_sensors() {
    // A is far cheaper than B, well beyond any uncertainty. The initial STOP
    // check must run even when no placements are permitted.
    let clear = vec![
        candidate("A", 0.1, 0.1, 0.001, 0.001, 1.0),
        candidate("B", 2.0, 2.0, 0.001, 0.001, 1.0),
    ];
    let report = run_campaign(&clear, 0.01, 0).expect("clear campaign succeeds");
    assert_eq!(report.sensors_placed, 0);
    assert!(report.decision_robust);
    assert_eq!(report.chosen_design, "A");
    assert_eq!(report.evpi_trace, vec![report.initial_evpi]);
}

#[test]
fn the_campaign_is_deterministic() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let a = run_campaign(&candidates, 0.01, 12).expect("first campaign succeeds");
    let b = run_campaign(&candidates, 0.01, 12).expect("replay succeeds");
    assert_eq!(a, b);
}

#[test]
fn candidate_construction_is_fail_closed_and_canonicalizes_zero() {
    assert_eq!(
        Candidate::new(" A", 0.0, 0.0, 0.0, 1.0, 1.0),
        Err(CandidateError::InvalidName {
            reason: "surrounding-whitespace"
        })
    );
    assert!(matches!(
        Candidate::new("x".repeat(129), 0.0, 0.0, 0.0, 1.0, 1.0),
        Err(CandidateError::InvalidName { reason: "too-long" })
    ));
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            Candidate::new("A", invalid, 0.0, 0.0, 1.0, 1.0),
            Err(CandidateError::InvalidNumber { field: "truth", .. })
        ));
    }
    assert!(matches!(
        Candidate::new("A", 0.0, 0.0, -1.0, 1.0, 1.0),
        Err(CandidateError::InvalidNumber {
            field: "prior_var",
            ..
        })
    ));
    assert!(matches!(
        Candidate::new("A", 0.0, 0.0, 1.0, 0.0, 1.0),
        Err(CandidateError::InvalidNumber {
            field: "sensor_noise",
            ..
        })
    ));
    assert!(matches!(
        Candidate::new("A", 0.0, 0.0, 1.0, 1.0, f64::NAN),
        Err(CandidateError::InvalidNumber {
            field: "sensor_cost",
            ..
        })
    ));

    let zero = candidate("zero", -0.0, -0.0, -0.0, 1.0, 1.0);
    assert_eq!(zero.truth().to_bits(), 0.0_f64.to_bits());
    assert_eq!(zero.prior_mean().to_bits(), 0.0_f64.to_bits());
    assert_eq!(zero.prior_variance().to_bits(), 0.0_f64.to_bits());
}

#[test]
fn campaign_inputs_are_checked_before_work_starts() {
    assert_eq!(run_campaign(&[], 0.0, 0), Err(OedError::NoCandidates));
    let one = candidate("A", 0.0, 0.0, 1.0, 1.0, 1.0);
    for threshold in [f64::NAN, f64::INFINITY, -1.0] {
        assert_eq!(
            run_campaign(std::slice::from_ref(&one), threshold, 0),
            Err(OedError::InvalidThreshold)
        );
    }
    assert_eq!(
        run_campaign(&[one.clone(), one.clone()], 0.0, 0),
        Err(OedError::DuplicateCandidate {
            name: "A".to_string()
        })
    );
    assert_eq!(
        run_campaign(std::slice::from_ref(&one), 0.0, MAX_CAMPAIGN_SENSORS + 1),
        Err(OedError::TooManySensors {
            count: MAX_CAMPAIGN_SENSORS + 1,
            max: MAX_CAMPAIGN_SENSORS
        })
    );
    let too_many = vec![one.clone(); MAX_CAMPAIGN_CANDIDATES + 1];
    assert_eq!(
        run_campaign(&too_many, 0.0, 0),
        Err(OedError::TooManyCandidates {
            count: MAX_CAMPAIGN_CANDIDATES + 1,
            max: MAX_CAMPAIGN_CANDIDATES
        })
    );
    let candidates: Vec<Candidate> = (0..MAX_CAMPAIGN_CANDIDATES)
        .map(|index| candidate(&format!("C{index}"), 0.0, 0.0, 1.0, 1.0, 1.0))
        .collect();
    assert_eq!(
        run_campaign(&candidates, 0.0, MAX_CAMPAIGN_SENSORS),
        Err(OedError::WorkBudgetExceeded {
            candidates: MAX_CAMPAIGN_CANDIDATES,
            max_sensors: MAX_CAMPAIGN_SENSORS,
            evaluations: MAX_CAMPAIGN_CANDIDATES * MAX_CAMPAIGN_CANDIDATES * MAX_CAMPAIGN_SENSORS,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS
        })
    );
}

#[test]
fn zero_total_prior_variance_has_defined_reduction_and_unbounded_allocation() {
    let exact = vec![
        candidate("A", 0.0, 0.0, 0.0, 0.01, 1.0),
        candidate("B", 1.0, 1.0, 0.0, 0.01, 2.0),
    ];
    let report = run_campaign(&exact, 0.0, 0).expect("exact campaign succeeds");
    assert_eq!(report.prior_total_variance.to_bits(), 0.0_f64.to_bits());
    assert_eq!(report.posterior_total_variance.to_bits(), 0.0_f64.to_bits());
    assert_eq!(report.variance_reduction.to_bits(), 0.0_f64.to_bits());
    assert!(report.decision_robust);
    assert!(
        report
            .allocation
            .iter()
            .all(|(_, tolerance)| tolerance.is_infinite() && tolerance.is_sign_positive())
    );
}

#[test]
fn evidence_identities_bind_unmeasured_inputs_and_realized_updates() {
    let baseline = demo_candidates().expect("compiled demo candidates are valid");
    let report = run_campaign(&baseline, 0.01, 12).expect("baseline succeeds");
    assert_eq!(report.assimilation_colors.len(), report.placements.len());

    let mut changed = baseline.clone();
    let d = &baseline[3];
    changed[3] = candidate(
        d.name(),
        d.truth() + 0.01,
        d.prior_mean(),
        d.prior_variance(),
        d.sensor_noise(),
        d.sensor_cost(),
    );
    let changed_report = run_campaign(&changed, 0.01, 12).expect("changed campaign succeeds");
    assert_eq!(report.placements, changed_report.placements);
    assert_ne!(
        estimator(&report.variance_color),
        estimator(&changed_report.variance_color),
        "an unmeasured candidate's truth is still a semantic campaign input"
    );
    assert_ne!(
        estimator(&report.evpi_color),
        estimator(&changed_report.evpi_color)
    );
    assert_ne!(
        estimator(&report.variance_color),
        estimator(&report.evpi_color),
        "distinct quantities require domain-separated identities"
    );
    assert!(color_leaf_identity_reason(estimator(&report.variance_color)).is_none());
    assert!(color_leaf_identity_reason(estimator(&report.evpi_color)).is_none());
}

#[test]
fn a_zero_placement_cap_does_not_claim_an_ambiguous_decision_is_robust() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let report = run_campaign(&candidates, 0.0, 0).expect("bounded campaign succeeds");
    assert_eq!(report.sensors_placed, 0);
    assert!(!report.decision_robust);
    assert_eq!(report.initial_evpi.to_bits(), report.final_evpi.to_bits());
}
