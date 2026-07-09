//! End-to-end battery: sensors are placed on the decision-relevant candidates,
//! the campaign stops when the decision is robust, and the posterior is exact.

use fs_evidence::Color;
use fs_oed_e2e::{Candidate, demo_candidates, run_campaign};

#[test]
fn sensors_target_the_decision_and_the_campaign_knows_when_to_stop() {
    let report = run_campaign(&demo_candidates(), 0.01, 12);
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
    // the posterior variance is exact (Verified); EVPI is Estimated.
    assert!(matches!(report.variance_color, Color::Verified { .. }));
    assert!(matches!(report.evpi_color, Color::Estimated { .. }));
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
    // A is far cheaper than B, well beyond any uncertainty: place nothing, stop.
    let clear = vec![
        Candidate {
            name: "A".into(),
            truth: 0.1,
            prior_mean: 0.1,
            prior_var: 0.001,
            sensor_noise: 0.001,
            sensor_cost: 1.0,
        },
        Candidate {
            name: "B".into(),
            truth: 2.0,
            prior_mean: 2.0,
            prior_var: 0.001,
            sensor_noise: 0.001,
            sensor_cost: 1.0,
        },
    ];
    let report = run_campaign(&clear, 0.01, 8);
    assert_eq!(report.sensors_placed, 0);
    assert!(report.decision_robust);
    assert_eq!(report.chosen_design, "A");
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(&demo_candidates(), 0.01, 12);
    let b = run_campaign(&demo_candidates(), 0.01, 12);
    assert_eq!(a.placements, b.placements);
    assert_eq!(a.final_evpi.to_bits(), b.final_evpi.to_bits());
}
