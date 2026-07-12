//! End-to-end battery: the makespan has a rigorous tropical enclosure, the
//! unique bottleneck + nominal slack are named, and VoI distinguishes action,
//! robust stop, and an inadequate action menu.

use fs_evidence::Color;
use fs_schedule_e2e::{
    MAX_CAMPAIGN_DECISION_ITEMS, MAX_CAMPAIGN_DECISION_WORK, ScheduleDisposition, ScheduleError,
    Study, run_campaign,
};
use fs_voi::{Action, ActionKind, DesignEstimate, Uncertainty};

// Candidate COSTS (fs-voi minimizes, so lower is better). A is cheapest but B
// is close and more uncertain — a live ranking ambiguity; C is clearly costlier.
fn designs() -> Vec<DesignEstimate> {
    vec![
        DesignEstimate::new(
            "A",
            0.60,
            Uncertainty {
                numerical: 0.05,
                statistical: 0.05,
                model: 0.08,
            },
        ),
        DesignEstimate::new(
            "B",
            0.65,
            Uncertainty {
                numerical: 0.08,
                statistical: 0.06,
                model: 0.10,
            },
        ),
        DesignEstimate::new(
            "C",
            0.90,
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
fn the_makespan_is_bounded_and_the_decision_stays_estimated() {
    let report =
        run_campaign(&studies(), &designs(), &actions(), 1e-6).expect("valid campaign schedule");
    // WHEN: tropical critical path (windtunnel-A -> decide = 13), enclosed.
    assert!(
        (report.makespan - 13.0).abs() < 1e-9,
        "makespan {}",
        report.makespan
    );
    let Color::Verified { lo, hi } = report.makespan_color else {
        panic!("makespan must carry enclosure evidence");
    };
    assert!(lo <= 13.0 && hi >= 13.0 && lo < hi);
    assert_eq!(report.critical_path, vec![3, 4]);
    assert!(report.critical_path_is_unique);
    assert_eq!(report.bottleneck_index, Some(3));
    assert_eq!(report.bottleneck.as_deref(), Some("windtunnel-A"));
    assert!(report.slack_studies.contains(&"hifi-B".to_string()));
    // WHETHER: A is the cheapest (leader), but B is a live contender — a real
    // ranking ambiguity, so the flip probability is substantial yet below 0.5
    // (the leader is still more likely to hold).
    assert_eq!(report.leading_design, "A");
    assert!(
        report.flip_risk > 0.05 && report.flip_risk < 0.5,
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
    assert_eq!(report.disposition, ScheduleDisposition::Act);
    assert!(report.recommendation_value_per_cost.is_some());
    assert!(matches!(
        report.recommendation_color,
        Color::Estimated { .. }
    ));
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
    // A is far cheaper than B and both are near-certain: no study is worth its
    // cost, so the campaign recommends STOP.
    let clear = vec![
        DesignEstimate::new(
            "A",
            0.30,
            Uncertainty {
                numerical: 0.01,
                statistical: 0.01,
                model: 0.01,
            },
        ),
        DesignEstimate::new(
            "B",
            1.00,
            Uncertainty {
                numerical: 0.01,
                statistical: 0.01,
                model: 0.01,
            },
        ),
    ];
    let report = run_campaign(&studies(), &clear, &actions(), 0.5).expect("valid robust campaign");
    assert!(report.should_stop, "rec {}", report.recommendation);
    assert_eq!(report.disposition, ScheduleDisposition::RobustStop);
    assert!(report.recommendation.starts_with("Stop:"));
    assert_eq!(report.leading_design, "A");
}

#[test]
fn a_fragile_decision_without_an_effective_action_expands_the_menu() {
    let report = run_campaign(&studies(), &designs(), &[], 1e-6)
        .expect("an empty action menu is a representable decision state");
    assert!(report.evpi > 1e-6);
    assert_eq!(report.disposition, ScheduleDisposition::NoEffectiveAction);
    assert!(!report.should_stop);
    assert!(report.recommendation.starts_with("Expand menu:"));
    assert_eq!(report.recommendation_value_per_cost, None);
}

#[test]
fn the_campaign_is_deterministic() {
    let a =
        run_campaign(&studies(), &designs(), &actions(), 1e-6).expect("first deterministic run");
    let b =
        run_campaign(&studies(), &designs(), &actions(), 1e-6).expect("second deterministic run");
    assert_eq!(a.makespan.to_bits(), b.makespan.to_bits());
    assert_eq!(a.recommendation, b.recommendation);
}

#[test]
fn malformed_campaign_inputs_refuse_before_verified_evidence() {
    assert_eq!(
        run_campaign(&[], &designs(), &actions(), 0.1),
        Err(ScheduleError::NoStudies)
    );
    assert_eq!(
        run_campaign(&studies(), &[], &actions(), 0.1),
        Err(ScheduleError::NoDesigns)
    );
    for latency in [f64::NAN, f64::INFINITY, -1.0] {
        let invalid = vec![Study::new("invalid", latency, vec![])];
        assert!(matches!(
            run_campaign(&invalid, &designs(), &actions(), 0.1),
            Err(ScheduleError::InvalidField {
                record: "study",
                field: "latency",
                ..
            })
        ));
    }
    let cyclic = vec![Study::new("a", 1.0, vec![1]), Study::new("b", 1.0, vec![0])];
    assert!(matches!(
        run_campaign(&cyclic, &designs(), &actions(), 0.1),
        Err(ScheduleError::Tropical(fs_tropical::TropicalError::Cyclic))
    ));
    assert!(matches!(
        run_campaign(&studies(), &designs(), &actions(), f64::NAN),
        Err(ScheduleError::InvalidField {
            field: "stop_threshold",
            ..
        })
    ));

    let mut invalid_designs = designs();
    invalid_designs[0].uncertainty.numerical = f64::MAX;
    invalid_designs[0].uncertainty.statistical = f64::MAX;
    assert!(matches!(
        run_campaign(&studies(), &invalid_designs, &actions(), 0.1),
        Err(ScheduleError::InvalidField {
            field: "uncertainty",
            ..
        })
    ));
    let mut invalid_actions = actions();
    invalid_actions[0].target_design = "missing".to_string();
    assert_eq!(
        run_campaign(&studies(), &designs(), &invalid_actions, 0.1),
        Err(ScheduleError::UnknownActionTarget { action: 0 })
    );
    let too_many_designs = vec![designs()[0].clone(); MAX_CAMPAIGN_DECISION_ITEMS + 1];
    assert!(matches!(
        run_campaign(&studies(), &too_many_designs, &actions(), 0.1),
        Err(ScheduleError::ResourceLimit {
            resource: "designs",
            ..
        })
    ));

    let mut invalid_name = studies();
    invalid_name[0].name = "bad\nname".to_string();
    assert!(matches!(
        run_campaign(&invalid_name, &designs(), &actions(), 0.1),
        Err(ScheduleError::InvalidField {
            record: "study",
            field: "name",
            ..
        })
    ));

    let over_work_designs = vec![designs()[0].clone(); 257];
    let over_work_actions = vec![actions()[0].clone(); 256];
    assert!(over_work_designs.len() * over_work_actions.len() > MAX_CAMPAIGN_DECISION_WORK);
    assert!(matches!(
        run_campaign(&studies(), &over_work_designs, &over_work_actions, 0.1),
        Err(ScheduleError::ResourceLimit {
            resource: "action-design evaluations",
            ..
        })
    ));
}

#[test]
fn exact_cartesian_work_boundary_is_admitted() {
    let designs: Vec<_> = (0..256)
        .map(|index| {
            DesignEstimate::new(
                format!("design-{index}"),
                f64::from(index),
                Uncertainty {
                    numerical: 0.01,
                    statistical: 0.01,
                    model: 0.01,
                },
            )
        })
        .collect();
    let actions: Vec<_> = (0..256)
        .map(|index| Action {
            name: format!("action-{index}"),
            kind: ActionKind::Simulate,
            target_design: "design-0".to_string(),
            reduction: 0.5,
            cost: f64::from(index) + 1.0,
        })
        .collect();
    assert_eq!(designs.len() * actions.len(), MAX_CAMPAIGN_DECISION_WORK);
    run_campaign(&studies(), &designs, &actions, 0.1).expect("exact work boundary");
}
