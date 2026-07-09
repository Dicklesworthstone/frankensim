//! End-to-end battery: the makespan is an exact (Verified) tropical critical
//! path, the bottleneck + slack are named, and VoI recommends the next study or
//! STOP when the decision is already robust.

use fs_evidence::Color;
use fs_schedule_e2e::{Study, run_campaign};
use fs_voi::{Action, ActionKind, DesignEstimate, Uncertainty};

fn designs() -> Vec<DesignEstimate> {
    vec![
        DesignEstimate::new(
            "A",
            1.00,
            Uncertainty {
                numerical: 0.05,
                statistical: 0.05,
                model: 0.08,
            },
        ),
        DesignEstimate::new(
            "B",
            0.95,
            Uncertainty {
                numerical: 0.08,
                statistical: 0.06,
                model: 0.10,
            },
        ),
        DesignEstimate::new(
            "C",
            0.70,
            Uncertainty {
                numerical: 0.05,
                statistical: 0.05,
                model: 0.05,
            },
        ),
    ]
}

fn actions() -> Vec<Action> {
    vec![
        Action {
            name: "hifi-B".into(),
            kind: ActionKind::Simulate,
            target_design: "B".into(),
            reduction: 0.9,
            cost: 8.0,
        },
        Action {
            name: "sample-B".into(),
            kind: ActionKind::Sample,
            target_design: "B".into(),
            reduction: 0.7,
            cost: 4.0,
        },
        Action {
            name: "windtunnel-A".into(),
            kind: ActionKind::Test,
            target_design: "A".into(),
            reduction: 0.8,
            cost: 12.0,
        },
    ]
}

fn studies() -> Vec<Study> {
    vec![
        Study::new("surrogate-B", 2.0, vec![]),      // 0
        Study::new("hifi-B", 8.0, vec![0]),          // 1 (needs the surrogate first)
        Study::new("sample-scenarios", 4.0, vec![]), // 2
        Study::new("windtunnel-A", 12.0, vec![]),    // 3 (the long pole)
        Study::new("decide", 1.0, vec![1, 2, 3]),    // 4
    ]
}

#[test]
fn the_schedule_and_the_decision_are_both_certified() {
    let report = run_campaign(&studies(), &designs(), &actions(), 1e-6);
    // WHEN: exact tropical critical path (windtunnel-A -> decide = 13), Verified.
    assert!(
        (report.makespan - 13.0).abs() < 1e-9,
        "makespan {}",
        report.makespan
    );
    assert!(matches!(report.makespan_color, Color::Verified { .. }));
    assert_eq!(report.critical_path, vec![3, 4]);
    assert_eq!(report.bottleneck.as_deref(), Some("windtunnel-A"));
    assert!(report.slack_studies.contains(&"hifi-B".to_string()));
    // WHETHER: A leads, but B is a live contender — a real ranking ambiguity
    // (B carries more uncertainty, so the flip probability is substantial).
    assert_eq!(report.leading_design, "A");
    assert!(
        report.flip_risk > 0.05 && report.flip_risk < 1.0,
        "flip {}",
        report.flip_risk
    );
    assert!(report.evpi > 0.0);
    assert!(
        report.recommendation.starts_with("Act:"),
        "rec {}",
        report.recommendation
    );
    assert!(!report.should_stop);
    println!(
        "{{\"campaign\":\"campaign-schedule\",\"makespan\":{},\"bottleneck\":{:?},\"evpi\":{:.4},\
         \"leader\":\"{}\",\"flip_risk\":{:.3},\"rec\":\"{}\"}}",
        report.makespan,
        report.bottleneck,
        report.evpi,
        report.leading_design,
        report.flip_risk,
        report.recommendation,
    );
}

#[test]
fn a_robust_decision_recommends_stop() {
    // C is dominant and the others are far behind: no study is worth its cost.
    let clear = vec![
        DesignEstimate::new(
            "A",
            1.0,
            Uncertainty {
                numerical: 0.01,
                statistical: 0.01,
                model: 0.01,
            },
        ),
        DesignEstimate::new(
            "B",
            0.3,
            Uncertainty {
                numerical: 0.01,
                statistical: 0.01,
                model: 0.01,
            },
        ),
    ];
    let report = run_campaign(&studies(), &clear, &actions(), 0.5);
    assert!(report.should_stop, "rec {}", report.recommendation);
    assert!(report.recommendation.starts_with("Stop:"));
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(&studies(), &designs(), &actions(), 1e-6);
    let b = run_campaign(&studies(), &designs(), &actions(), 1e-6);
    assert_eq!(a.makespan.to_bits(), b.makespan.to_bits());
    assert_eq!(a.recommendation, b.recommendation);
}
