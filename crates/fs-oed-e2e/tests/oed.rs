//! End-to-end battery: sensors are placed on the decision-relevant candidates,
//! the campaign stops when the decision is robust, and the posterior remains
//! honestly model-form Estimated.

use fs_assimilate::AssimError;
use fs_evidence::{Color, color_leaf_identity_reason};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_oed_e2e::{
    Candidate, CandidateError, MAX_CAMPAIGN_CANDIDATES, MAX_CAMPAIGN_EVALUATIONS,
    MAX_CAMPAIGN_SENSORS, OedError, OedReport, demo_candidates, run_campaign,
};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x6f65_642d_6532_6501,
    kernel_id: 1,
    tile: 0,
    iteration: 0,
};

fn with_stream_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    mode: ExecMode,
    stream: StreamKey,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, stream, budget, mode);
        f(&cx)
    });
    let stats = pool.stats();
    assert!(
        stats.quiescent(),
        "Cx arena must be quiescent after scope: {}",
        stats.to_json()
    );
    result
}

fn with_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    mode: ExecMode,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    with_stream_cx(gate, budget, mode, TEST_STREAM, f)
}

fn campaign(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
) -> Result<fs_oed_e2e::OedReport, OedError> {
    with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| run_campaign(candidates, threshold, max_sensors, cx),
    )
}

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

fn assert_same_non_identity_report(left: &OedReport, right: &OedReport) {
    assert_eq!(left.placements(), right.placements());
    assert_eq!(left.sensors_placed(), right.sensors_placed());
    assert_eq!(
        left.prior_total_variance().to_bits(),
        right.prior_total_variance().to_bits()
    );
    assert_eq!(
        left.posterior_total_variance().to_bits(),
        right.posterior_total_variance().to_bits()
    );
    assert_eq!(
        left.variance_reduction().to_bits(),
        right.variance_reduction().to_bits()
    );
    assert_eq!(
        left.initial_evpi().to_bits(),
        right.initial_evpi().to_bits()
    );
    assert_eq!(left.final_evpi().to_bits(), right.final_evpi().to_bits());
    assert_eq!(left.decision_robust(), right.decision_robust());
    assert_eq!(left.chosen_design(), right.chosen_design());
    assert_eq!(left.allocation().len(), right.allocation().len());
    for ((left_name, left_tolerance), (right_name, right_tolerance)) in
        left.allocation().iter().zip(right.allocation())
    {
        assert_eq!(left_name, right_name);
        assert_eq!(left_tolerance.to_bits(), right_tolerance.to_bits());
    }
    assert_eq!(left.evpi_trace().len(), right.evpi_trace().len());
    for (left_evpi, right_evpi) in left.evpi_trace().iter().zip(right.evpi_trace()) {
        assert_eq!(left_evpi.to_bits(), right_evpi.to_bits());
    }
    assert_eq!(left.posteriors().len(), right.posteriors().len());
    for (left_posterior, right_posterior) in left.posteriors().iter().zip(right.posteriors()) {
        assert_eq!(left_posterior.name, right_posterior.name);
        assert_eq!(
            left_posterior.mean.to_bits(),
            right_posterior.mean.to_bits()
        );
        assert_eq!(
            left_posterior.variance.to_bits(),
            right_posterior.variance.to_bits()
        );
    }
    assert_eq!(
        left.assimilation_colors().len(),
        right.assimilation_colors().len()
    );
    for (left_color, right_color) in left
        .assimilation_colors()
        .iter()
        .zip(right.assimilation_colors())
    {
        match (left_color, right_color) {
            (
                Color::Estimated {
                    dispersion: left_dispersion,
                    ..
                },
                Color::Estimated {
                    dispersion: right_dispersion,
                    ..
                },
            ) => assert_eq!(left_dispersion.to_bits(), right_dispersion.to_bits()),
            other => panic!("expected matching Estimated assimilation colors, got {other:?}"),
        }
    }
    for (left_color, right_color) in [
        (left.variance_color(), right.variance_color()),
        (left.evpi_color(), right.evpi_color()),
    ] {
        match (left_color, right_color) {
            (
                Color::Estimated {
                    dispersion: left_dispersion,
                    ..
                },
                Color::Estimated {
                    dispersion: right_dispersion,
                    ..
                },
            ) => assert_eq!(left_dispersion.to_bits(), right_dispersion.to_bits()),
            other => panic!("expected matching Estimated report colors, got {other:?}"),
        }
    }
}

#[test]
fn sensors_target_the_decision_and_the_campaign_knows_when_to_stop() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let report = campaign(&candidates, 0.01, 12).expect("demo campaign succeeds");
    // sensors WERE placed, and on the decision-relevant contenders (A and/or B),
    // never on the clearly-dominated D.
    assert!(report.sensors_placed() > 0, "no sensors placed");
    assert!(report.placements().iter().any(|n| n == "A" || n == "B"));
    assert!(
        !report.placements().contains(&"D".to_string()),
        "wasted a sensor on D"
    );
    // measurement sharpened the beliefs; EVPI fell.
    assert!(report.variance_reduction() > 0.0);
    assert!(
        report.final_evpi() < report.initial_evpi(),
        "EVPI did not fall"
    );
    // the campaign STOPPED because the decision became robust.
    assert!(report.decision_robust(), "did not reach a robust decision");
    assert!(
        report.final_evpi() <= 0.01 + 1e-9,
        "final EVPI {}",
        report.final_evpi()
    );
    // A is the true best and is chosen.
    assert_eq!(report.chosen_design(), "A");
    // Kalman variance remains model-form Estimated until independently
    // certified; EVPI is Estimated as well.
    assert!(matches!(report.variance_color(), Color::Estimated { .. }));
    assert!(matches!(report.evpi_color(), Color::Estimated { .. }));
    assert_eq!(report.evpi_trace().len(), report.sensors_placed() + 1);
    assert_eq!(report.assimilation_colors().len(), report.sensors_placed());
    for color in report.assimilation_colors() {
        let identity = estimator(color);
        assert!(color_leaf_identity_reason(identity).is_none());
    }
    // a cost-optimal precision budget was allocated across candidates.
    assert_eq!(report.allocation().len(), 4);
    assert!(report.allocation().iter().all(|(_, t)| *t > 0.0));
    println!(
        "{{\"campaign\":\"sensorforge\",\"placements\":{:?},\"sensors\":{},\"var_reduction\":{:.3},\
         \"initial_evpi\":{:.4},\"final_evpi\":{:.4},\"robust\":{},\"chosen\":\"{}\"}}",
        report.placements(),
        report.sensors_placed(),
        report.variance_reduction(),
        report.initial_evpi(),
        report.final_evpi(),
        report.decision_robust(),
        report.chosen_design(),
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
    let report = campaign(&clear, 0.01, 0).expect("clear campaign succeeds");
    assert_eq!(report.sensors_placed(), 0);
    assert!(report.decision_robust());
    assert_eq!(report.chosen_design(), "A");
    assert_eq!(report.evpi_trace(), &[report.initial_evpi()]);
}

#[test]
fn the_campaign_is_deterministic() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let a = campaign(&candidates, 0.01, 12).expect("first campaign succeeds");
    let b = campaign(&candidates, 0.01, 12).expect("replay succeeds");
    assert_eq!(a, b);
}

#[test]
fn tied_final_means_use_candidate_identity_not_menu_order() {
    let candidates = vec![
        candidate("Zulu", 0.0, 0.0, 1.0, 1.0, 1.0),
        candidate("Alpha", 0.0, 0.0, 1.0, 1.0, 1.0),
    ];
    let forward = campaign(&candidates, 0.0, 0).expect("tied campaign succeeds");
    let reversed = campaign(&candidates.into_iter().rev().collect::<Vec<_>>(), 0.0, 0)
        .expect("permuted tied campaign succeeds");

    assert_eq!(forward.chosen_design(), "Alpha");
    assert_eq!(reversed.chosen_design(), "Alpha");
}

#[test]
fn tied_multicandidate_evpi_and_sensor_policy_ignore_menu_order() {
    let candidates = vec![
        candidate("Alpha", 0.0, 0.0, 9.0, 0.25, 1.0),
        candidate("Beta", 0.0, 0.0, 1.0, 0.25, 1.0),
        candidate("Gamma", 0.0, 0.0, 0.01, 0.25, 1.0),
    ];
    let forward = campaign(&candidates, 0.0, 1).expect("tied campaign succeeds");
    let reversed_candidates: Vec<Candidate> = candidates.into_iter().rev().collect();
    let reversed = campaign(&reversed_candidates, 0.0, 1).expect("permuted campaign succeeds");

    assert_eq!(
        forward.initial_evpi().to_bits(),
        reversed.initial_evpi().to_bits()
    );
    assert_eq!(
        forward.final_evpi().to_bits(),
        reversed.final_evpi().to_bits()
    );
    assert_eq!(forward.placements(), reversed.placements());
    assert_eq!(forward.chosen_design(), reversed.chosen_design());
}

#[test]
fn sensor_planning_uses_declared_noise_and_matches_the_kalman_variance() {
    // This order is adversarial to the former tie policy: the nearly useless
    // high-noise sensor is last, so the old fixed 85% reduction selected it.
    let candidates = vec![
        candidate("low-noise", 0.0, 0.0, 1.0, 0.01, 1.0),
        candidate("high-noise", 0.0, 0.0, 1.0, 100.0, 1.0),
    ];
    let report = campaign(&candidates, 0.0, 1).expect("noise-aware campaign succeeds");
    assert_eq!(
        report.placements().first().map(String::as_str),
        Some("low-noise")
    );
    assert_eq!(report.placements().len(), 1);
    let low = report
        .posteriors()
        .iter()
        .find(|posterior| posterior.name == "low-noise")
        .expect("low-noise posterior retained");
    let expected = 1.0 * 0.01 / (1.0 + 0.01);
    assert!(
        (low.variance - expected).abs() <= 8.0 * f64::EPSILON,
        "planned and realized scalar Kalman variance must agree: got {}, expected {expected}",
        low.variance
    );

    let reversed: Vec<Candidate> = candidates.into_iter().rev().collect();
    let replay = campaign(&reversed, 0.0, 1).expect("permuted campaign succeeds");
    assert_eq!(
        replay.placements().first().map(String::as_str),
        Some("low-noise"),
        "sensor informativeness, not menu order, must determine the action"
    );
    assert_eq!(replay.placements().len(), 1);
}

#[test]
fn sensor_value_has_the_correct_noise_limits() {
    let candidates = vec![
        candidate("informative", 0.0, 0.0, 1.0, 1.0e-12, 1.0),
        candidate("negligible", 0.0, 0.0, 1.0, 1.0e300, 1.0),
    ];
    let report = campaign(&candidates, 0.0, 1).expect("extreme finite noises are admitted");
    assert_eq!(
        report.placements().first().map(String::as_str),
        Some("informative")
    );
    assert_eq!(report.placements().len(), 1);
    let informative = report
        .posteriors()
        .iter()
        .find(|posterior| posterior.name == "informative")
        .expect("informative posterior retained");
    assert!(
        informative.variance <= 1.000_000_000_001e-12,
        "near-noiseless sensing should nearly collapse the declared variance"
    );
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
    assert_eq!(campaign(&[], 0.0, 0), Err(OedError::NoCandidates));
    let one = candidate("A", 0.0, 0.0, 1.0, 1.0, 1.0);
    for threshold in [f64::NAN, f64::INFINITY, -1.0] {
        assert_eq!(
            campaign(std::slice::from_ref(&one), threshold, 0),
            Err(OedError::InvalidThreshold)
        );
    }
    assert_eq!(
        campaign(&[one.clone(), one.clone()], 0.0, 0),
        Err(OedError::DuplicateCandidate {
            name: "A".to_string()
        })
    );
    assert_eq!(
        campaign(std::slice::from_ref(&one), 0.0, MAX_CAMPAIGN_SENSORS + 1),
        Err(OedError::TooManySensors {
            count: MAX_CAMPAIGN_SENSORS + 1,
            max: MAX_CAMPAIGN_SENSORS
        })
    );
    let too_many = vec![one.clone(); MAX_CAMPAIGN_CANDIDATES + 1];
    assert_eq!(
        campaign(&too_many, 0.0, 0),
        Err(OedError::TooManyCandidates {
            count: MAX_CAMPAIGN_CANDIDATES + 1,
            max: MAX_CAMPAIGN_CANDIDATES
        })
    );
    let candidates: Vec<Candidate> = (0..MAX_CAMPAIGN_CANDIDATES)
        .map(|index| candidate(&format!("C{index}"), 0.0, 0.0, 1.0, 1.0, 1.0))
        .collect();
    assert_eq!(
        campaign(&candidates, 0.0, MAX_CAMPAIGN_SENSORS),
        Err(OedError::WorkBudgetExceeded {
            candidates: MAX_CAMPAIGN_CANDIDATES,
            max_sensors: MAX_CAMPAIGN_SENSORS,
            evaluations: MAX_CAMPAIGN_CANDIDATES
                * MAX_CAMPAIGN_CANDIDATES
                * MAX_CAMPAIGN_SENSORS
                * 11,
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
    let report = campaign(&exact, 0.0, 0).expect("exact campaign succeeds");
    assert_eq!(report.prior_total_variance().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        report.posterior_total_variance().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(report.variance_reduction().to_bits(), 0.0_f64.to_bits());
    assert!(report.decision_robust());
    assert!(
        report
            .allocation()
            .iter()
            .all(|(_, tolerance)| tolerance.is_infinite() && tolerance.is_sign_positive())
    );
}

#[test]
fn mixed_zero_variance_priors_have_exact_realized_allocation_work() {
    let candidates = vec![
        candidate("A", 0.0, 0.0, 0.0, 1.0, 1.0),
        candidate("B", 1.0, 1.0, 1.0, 1.0, 1.0),
        candidate("C", 2.0, 2.0, 0.0, 1.0, 1.0),
    ];
    let report = campaign(&candidates, 1.0e9, 2)
        .expect("mixed positive-prior cardinality has an exact realized work ledger");

    assert_eq!(report.sensors_placed(), 0);
    assert_eq!(report.allocation().len(), 3);
    assert!(report.allocation()[0].1.is_infinite());
    assert!(report.allocation()[1].1.is_finite());
    assert!(report.allocation()[1].1 > 0.0);
    assert!(report.allocation()[2].1.is_infinite());
}

#[test]
fn evidence_identities_bind_unmeasured_inputs_and_realized_updates() {
    let baseline = demo_candidates().expect("compiled demo candidates are valid");
    let report = campaign(&baseline, 0.01, 12).expect("baseline succeeds");
    assert_eq!(
        report.assimilation_colors().len(),
        report.placements().len()
    );

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
    let changed_report = campaign(&changed, 0.01, 12).expect("changed campaign succeeds");
    assert_eq!(report.placements(), changed_report.placements());
    assert_ne!(
        estimator(report.variance_color()),
        estimator(changed_report.variance_color()),
        "an unmeasured candidate's truth is still a semantic campaign input"
    );
    assert_ne!(
        estimator(report.evpi_color()),
        estimator(changed_report.evpi_color())
    );
    assert_ne!(
        estimator(report.variance_color()),
        estimator(report.evpi_color()),
        "distinct quantities require domain-separated identities"
    );
    assert!(color_leaf_identity_reason(estimator(report.variance_color())).is_none());
    assert!(color_leaf_identity_reason(estimator(report.evpi_color())).is_none());
    assert!(estimator(report.variance_color()).starts_with("sensorforge-posterior-variance:v5:"));
    assert!(estimator(report.evpi_color()).starts_with("sensorforge-evpi:v5:"));
}

#[test]
fn a_zero_placement_cap_does_not_claim_an_ambiguous_decision_is_robust() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let report = campaign(&candidates, 0.0, 0).expect("bounded campaign succeeds");
    assert_eq!(report.sensors_placed(), 0);
    assert!(!report.decision_robust());
    assert_eq!(
        report.initial_evpi().to_bits(),
        report.final_evpi().to_bits()
    );
}

#[test]
fn a_completed_zero_value_action_round_has_exact_realized_work() {
    let candidates = vec![
        candidate("A", 0.0, 0.0, 1.0e-300, 1.0e300, 1.0),
        candidate("B", 0.0, 0.0, 1.0e-300, 1.0e300, 1.0),
    ];
    let report = campaign(&candidates, 0.0, 1)
        .expect("a no-positive-action STOP has a complete realized work ledger");
    assert!(report.initial_evpi() > 0.0);
    assert_eq!(report.sensors_placed(), 0);
    assert!(!report.decision_robust());
    assert_eq!(
        report.initial_evpi().to_bits(),
        report.final_evpi().to_bits()
    );
}

#[test]
fn g4_pre_cancel_is_bounded_and_never_publishes_a_partial_report() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let gate = CancelGate::new();
    gate.request();
    let cancelled = with_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        run_campaign(&candidates, 0.01, 12, cx)
    });
    match cancelled {
        Err(OedError::Cancelled {
            phase,
            completed_placements,
            completed_work_units,
            admitted_work_units,
        }) => {
            assert_eq!(phase, "campaign admission");
            assert_eq!(completed_placements, 0);
            assert_eq!(completed_work_units, candidates.len() as u128);
            assert_eq!(admitted_work_units, 2_545);
        }
        other => panic!("expected structured admission cancellation, got {other:?}"),
    }

    // The cancelled scratch run cannot leak state into a later replay.
    let recovered = campaign(&candidates, 0.01, 12).expect("healthy replay succeeds");
    let baseline = campaign(&candidates, 0.01, 12).expect("baseline succeeds");
    assert_eq!(recovered, baseline);
}

#[test]
fn g4_poll_quota_bounds_action_tile_latency_without_partial_publication() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let budget = Budget::new().with_poll_quota(7);
    let cancelled = with_cx(&CancelGate::new(), budget, ExecMode::Deterministic, |cx| {
        run_campaign(&candidates, 0.0, 12, cx)
    });
    match cancelled {
        Err(OedError::Cancelled {
            phase,
            completed_placements,
            completed_work_units,
            admitted_work_units,
        }) => {
            assert_eq!(phase, "action-value tile");
            assert_eq!(completed_placements, 0);
            // Five setup scans (20 units), dynamic action construction (4),
            // then one nine-node outcome-integrated action tile (45 units).
            assert_eq!(completed_work_units, 69);
            assert!(admitted_work_units > completed_work_units);
        }
        other => panic!("expected bounded action-tile cancellation, got {other:?}"),
    }
}

#[test]
fn g4_lower_assimilation_cancellation_preserves_both_progress_ledgers() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let mut nested = None;
    for quota in 1..=128 {
        let budget = Budget::new().with_poll_quota(quota);
        let result = with_cx(&CancelGate::new(), budget, ExecMode::Deterministic, |cx| {
            run_campaign(&candidates, 0.0, 1, cx)
        });
        if let Err(error @ OedError::AssimilationCancelled { .. }) = result {
            nested = Some(error);
            break;
        }
    }

    match nested.expect("a bounded quota reaches and cancels the lower assimilation") {
        OedError::AssimilationCancelled {
            candidate,
            completed_placements,
            completed_work_units,
            admitted_work_units,
            source,
        } => {
            let AssimError::Cancelled {
                phase,
                completed,
                planned,
            } = source.as_ref()
            else {
                panic!("expected nested cancellation source, got {source:?}");
            };
            assert!(candidates.iter().any(|entry| entry.name() == candidate));
            assert_eq!(completed_placements, 0);
            // The OED ledger includes the selected-action lookup and the
            // successfully constructed point-sensor observation, but no
            // posterior or placement has been committed.
            assert_eq!(completed_work_units, 209);
            assert!(admitted_work_units > completed_work_units);
            assert!(!phase.is_empty());
            assert!(completed < planned);
        }
        other => panic!("expected nested structured cancellation, got {other:?}"),
    }

    let recovered = campaign(&candidates, 0.0, 1).expect("healthy replay succeeds");
    let baseline = campaign(&candidates, 0.0, 1).expect("baseline succeeds");
    assert_eq!(recovered, baseline);
}

#[test]
fn g4_finalization_quota_sweep_covers_identity_and_publication_boundaries() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let baseline = campaign(&candidates, 0.0, 1).expect("one-placement baseline succeeds");
    let mut cancelled_phases = Vec::new();
    let mut completed = false;

    for quota in 1..=512 {
        let result = with_cx(
            &CancelGate::new(),
            Budget::new().with_poll_quota(quota),
            ExecMode::Deterministic,
            |cx| run_campaign(&candidates, 0.0, 1, cx),
        );
        match result {
            Err(OedError::Cancelled {
                phase,
                completed_work_units,
                admitted_work_units,
                ..
            }) => {
                assert_eq!(admitted_work_units, 301);
                assert!(completed_work_units < admitted_work_units);
                cancelled_phases.push(phase);
            }
            Err(OedError::AssimilationCancelled {
                completed_placements,
                completed_work_units,
                admitted_work_units,
                source,
                ..
            }) => {
                assert!(matches!(source.as_ref(), AssimError::Cancelled { .. }));
                assert_eq!(completed_placements, 0);
                assert_eq!(completed_work_units, 209);
                assert_eq!(admitted_work_units, 301);
            }
            Ok(report) => {
                assert_same_non_identity_report(&baseline, &report);
                completed = true;
                break;
            }
            other => panic!("unexpected quota-sweep outcome: {other:?}"),
        }
    }

    assert!(completed, "a sufficient finite quota must complete");
    assert!(cancelled_phases.contains(&"report identity placements"));
    assert!(cancelled_phases.contains(&"report identity allocation"));
    assert!(cancelled_phases.contains(&"report identity EVPI trace"));
    assert!(cancelled_phases.contains(&"report identity assimilation colors"));
    assert!(cancelled_phases.contains(&"report identity hash"));
    assert!(cancelled_phases.contains(&"report publication"));
    let recovered = campaign(&candidates, 0.0, 1).expect("healthy replay succeeds");
    assert_eq!(recovered, baseline);
}

#[test]
fn g5_execution_mode_and_complete_budget_are_bound_into_evidence() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let run = |budget, mode, stream| {
        with_stream_cx(&CancelGate::new(), budget, mode, stream, |cx| {
            // Zero placement isolates the OED manifest from lower-layer
            // assimilation estimator identities.
            run_campaign(&candidates, 1.0e9, 0, cx)
        })
        .expect("zero-placement campaign succeeds")
    };
    let deterministic = run(Budget::INFINITE, ExecMode::Deterministic, TEST_STREAM);
    let fast = run(Budget::INFINITE, ExecMode::Fast, TEST_STREAM);

    let mut deadline_budget = Budget::INFINITE;
    deadline_budget.deadline = Budget::ZERO.deadline;
    let mut poll_budget = Budget::INFINITE;
    poll_budget.poll_quota = 10_000;
    let mut cost_budget = Budget::INFINITE;
    cost_budget.cost_quota = Some(1_000_000);
    let mut priority_budget = Budget::INFINITE;
    priority_budget.priority = 7;
    let budget_variants = [
        run(deadline_budget, ExecMode::Deterministic, TEST_STREAM),
        run(poll_budget, ExecMode::Deterministic, TEST_STREAM),
        run(cost_budget, ExecMode::Deterministic, TEST_STREAM),
        run(priority_budget, ExecMode::Deterministic, TEST_STREAM),
    ];
    let stream_variants = [
        StreamKey {
            seed: TEST_STREAM.seed + 1,
            ..TEST_STREAM
        },
        StreamKey {
            kernel_id: TEST_STREAM.kernel_id + 1,
            ..TEST_STREAM
        },
        StreamKey {
            tile: TEST_STREAM.tile + 1,
            ..TEST_STREAM
        },
        StreamKey {
            iteration: TEST_STREAM.iteration + 1,
            ..TEST_STREAM
        },
    ]
    .map(|stream| run(Budget::INFINITE, ExecMode::Deterministic, stream));

    assert_same_non_identity_report(&deterministic, &fast);
    assert_ne!(
        estimator(deterministic.variance_color()),
        estimator(fast.variance_color())
    );
    assert_ne!(
        estimator(deterministic.evpi_color()),
        estimator(fast.evpi_color())
    );
    for variant in budget_variants {
        assert_same_non_identity_report(&deterministic, &variant);
        assert_ne!(
            estimator(deterministic.variance_color()),
            estimator(variant.variance_color())
        );
        assert_ne!(
            estimator(deterministic.evpi_color()),
            estimator(variant.evpi_color())
        );
    }
    for variant in stream_variants {
        assert_same_non_identity_report(&deterministic, &variant);
        assert_ne!(
            estimator(deterministic.variance_color()),
            estimator(variant.variance_color())
        );
        assert_ne!(
            estimator(deterministic.evpi_color()),
            estimator(variant.evpi_color())
        );
    }
}

#[test]
fn g5_admitted_work_shape_is_bound_when_realized_science_is_identical() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let no_capacity = campaign(&candidates, 1.0e9, 0).expect("zero-cap campaign succeeds");
    let unused_capacity =
        campaign(&candidates, 1.0e9, 12).expect("unused-capacity campaign succeeds");

    assert_same_non_identity_report(&no_capacity, &unused_capacity);
    assert_ne!(
        estimator(no_capacity.variance_color()),
        estimator(unused_capacity.variance_color())
    );
    assert_ne!(
        estimator(no_capacity.evpi_color()),
        estimator(unused_capacity.evpi_color())
    );
}

#[test]
fn g4_hostile_maximum_work_shape_is_admitted_and_the_next_shape_refuses() {
    let accepted: Vec<Candidate> = (0..15)
        .map(|index| candidate(&format!("C{index}"), 0.0, 0.0, 1.0, 1.0, 1.0))
        .collect();
    let gate = CancelGate::new();
    gate.request();
    let at_limit = with_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        run_campaign(&accepted, 0.0, MAX_CAMPAIGN_SENSORS, cx)
    });
    assert!(matches!(at_limit, Err(OedError::Cancelled { .. })));
    let bounded_later = with_cx(
        &CancelGate::new(),
        Budget::INFINITE.with_poll_quota(1),
        ExecMode::Deterministic,
        |cx| run_campaign(&accepted, 0.0, MAX_CAMPAIGN_SENSORS, cx),
    );
    assert!(matches!(
        bounded_later,
        Err(OedError::Cancelled {
            phase: "prior variance",
            completed_placements: 0,
            completed_work_units: 30,
            ..
        })
    ));

    let mut refused = accepted;
    refused.push(candidate("C15", 0.0, 0.0, 1.0, 1.0, 1.0));
    assert_eq!(
        with_cx(
            &CancelGate::new(),
            Budget::INFINITE,
            ExecMode::Deterministic,
            |cx| run_campaign(&refused, 0.0, MAX_CAMPAIGN_SENSORS, cx),
        ),
        Err(OedError::WorkBudgetExceeded {
            candidates: 16,
            max_sensors: MAX_CAMPAIGN_SENSORS,
            evaluations: 16 * 16 * MAX_CAMPAIGN_SENSORS * 11,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
        })
    );
}
