//! End-to-end battery: sensors are placed on the decision-relevant candidates,
//! the campaign stops when the decision is robust, and the posterior remains
//! honestly model-form Estimated.

use fs_assimilate::AssimError;
use fs_evidence::{Color, color_leaf_identity_reason};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_oed_e2e::{
    AcquisitionCost, AcquisitionCostSpec, Candidate, CandidateError, MAX_CAMPAIGN_CANDIDATES,
    MAX_CAMPAIGN_EVALUATIONS, MAX_CAMPAIGN_SENSORS, ObjectiveValue, OedError, OedReport,
    demo_candidates, run_campaign,
};
use fs_qty::parse::parse_qty;
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticQty, SemanticType, ValueForm};
use fs_qty::{Dims, QtyAny, QuantitySpec};

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
    let owned_clock = fs_exec::VirtualClock::new();
    let clock = &owned_clock;
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, stream, budget, mode).with_time_source(clock);
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
        |cx| {
            run_campaign(
                candidates,
                ObjectiveValue::dimensionless(threshold).expect("test threshold must be finite"),
                max_sensors,
                cx,
            )
        },
    )
}

fn objective(value: f64) -> ObjectiveValue {
    ObjectiveValue::dimensionless(value).expect("test objective must be finite")
}

fn semantic_objective(value: f64, kind: QuantityKind) -> ObjectiveValue {
    let semantic_type = SemanticType::new(kind, ValueForm::Static);
    ObjectiveValue::semantic(
        SemanticQty::new(QtyAny::new(value, kind.expected_dims()), semantic_type)
            .expect("test semantic objective must be admissible"),
    )
}

fn semantic_acquisition_cost(value: f64, kind: QuantityKind) -> AcquisitionCost {
    let semantic_type = SemanticType::new(kind, ValueForm::Static);
    AcquisitionCost::semantic(
        SemanticQty::new(QtyAny::new(value, kind.expected_dims()), semantic_type)
            .expect("test semantic cost must be admissible"),
    )
    .expect("test semantic cost must be positive")
}

fn candidate_with_acquisition_cost(
    name: &str,
    truth: ObjectiveValue,
    prior_mean: ObjectiveValue,
    prior_variance: QtyAny,
    sensor_noise_variance: QtyAny,
    acquisition_cost: AcquisitionCost,
    allocation_cost_weight: f64,
) -> Candidate {
    Candidate::with_acquisition_cost(
        name,
        truth,
        prior_mean,
        prior_variance,
        sensor_noise_variance,
        acquisition_cost,
        QtyAny::dimensionless(allocation_cost_weight),
    )
    .expect("typed acquisition-cost candidate must be admissible")
}

fn typed_candidate(
    name: &str,
    truth: ObjectiveValue,
    prior_mean: ObjectiveValue,
    prior_variance: QtyAny,
    sensor_noise_variance: QtyAny,
) -> Candidate {
    Candidate::new(
        name,
        truth,
        prior_mean,
        prior_variance,
        sensor_noise_variance,
        QtyAny::dimensionless(1.0),
    )
    .expect("typed test candidate must satisfy the checked constructor")
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
        objective(truth),
        objective(prior_mean),
        QtyAny::dimensionless(prior_var),
        QtyAny::dimensionless(sensor_noise),
        QtyAny::dimensionless(sensor_cost),
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
    // Retained bytes are determined by the realized science (names,
    // traces, summaries, identity strings); admitted/consumed byte
    // units deliberately vary with the admitted shape and the budget
    // encoding bound into the identity preimage, so identity-varying
    // fixtures compare only the retained ledger here.
    assert_eq!(left.retained_byte_units(), right.retained_byte_units());
    assert_eq!(
        left.prior_total_variance().value.to_bits(),
        right.prior_total_variance().value.to_bits()
    );
    assert_eq!(
        left.posterior_total_variance().value.to_bits(),
        right.posterior_total_variance().value.to_bits()
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
    for (left_evpi, right_evpi) in left.evpi_trace().zip(right.evpi_trace()) {
        assert_eq!(left_evpi.to_bits(), right_evpi.to_bits());
    }
    assert_eq!(left.posteriors().len(), right.posteriors().len());
    for (left_posterior, right_posterior) in left.posteriors().iter().zip(right.posteriors()) {
        assert_eq!(left_posterior.name(), right_posterior.name());
        assert_eq!(
            left_posterior.mean().to_bits(),
            right_posterior.mean().to_bits()
        );
        assert_eq!(
            left_posterior.variance().value.to_bits(),
            right_posterior.variance().value.to_bits()
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
    // sensors WERE placed, prioritizing the decision-relevant
    // contenders (A and/or B); every other alternative gets at most an
    // exclusion-certifying measurement.
    assert!(report.sensors_placed() > 0, "no sensors placed");
    assert!(report.placements().iter().any(|n| n == "A" || n == "B"));
    // Under the FULL opportunity-loss algebra (bead sj31i.5) the
    // campaign certifies GLOBAL robustness: after resolving the
    // contenders it spends single placements on C and finally D —
    // exactly the alternatives whose residual optimality probability
    // blocks the certificate — instead of declaring robustness while
    // ignoring them (the retired top-two surrogate's laundering).
    assert!(
        report.placements()[..2]
            .iter()
            .all(|n| n == "A" || n == "B"),
        "the first placements belong to the contenders: {:?}",
        report.placements()
    );
    assert!(
        report.placements().iter().filter(|n| *n == "D").count() <= 1,
        "D needs at most one exclusion-certifying measurement: {:?}",
        report.placements()
    );
    // measurement sharpened the beliefs; EVPI fell.
    assert!(report.variance_reduction() > 0.0);
    assert!(
        report.final_evpi().value() < report.initial_evpi().value(),
        "EVPI did not fall"
    );
    // the campaign STOPPED because the decision became robust.
    assert!(report.decision_robust(), "did not reach a robust decision");
    assert!(
        report.final_evpi().value() <= 0.01 + 1e-9,
        "final EVPI {}",
        report.final_evpi().value()
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
        report.initial_evpi().value(),
        report.final_evpi().value(),
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
    assert_eq!(
        report.evpi_trace().collect::<Vec<_>>(),
        [report.initial_evpi()]
    );
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
        .find(|posterior| posterior.name() == "low-noise")
        .expect("low-noise posterior retained");
    let expected = 1.0 * 0.01 / (1.0 + 0.01);
    assert!(
        (low.variance().value - expected).abs() <= 8.0 * f64::EPSILON,
        "planned and realized scalar Kalman variance must agree: got {}, expected {expected}",
        low.variance().value
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
        .find(|posterior| posterior.name() == "informative")
        .expect("informative posterior retained");
    assert!(
        informative.variance().value <= 1.000_000_000_001e-12,
        "near-noiseless sensing should nearly collapse the declared variance"
    );
}

#[test]
fn candidate_construction_is_fail_closed_and_canonicalizes_zero() {
    let valid = || {
        (
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(0.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(1.0),
        )
    };
    let (truth, mean, variance, noise, cost) = valid();
    assert_eq!(
        Candidate::new(" A", truth, mean, variance, noise, cost),
        Err(CandidateError::InvalidName {
            reason: "surrounding-whitespace"
        })
    );
    let (truth, mean, variance, noise, cost) = valid();
    assert!(matches!(
        Candidate::new("x".repeat(129), truth, mean, variance, noise, cost),
        Err(CandidateError::InvalidName { reason: "too-long" })
    ));
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            ObjectiveValue::dimensionless(invalid),
            Err(CandidateError::InvalidNumber {
                field: "objective value",
                requirement: "finite",
            })
        );
    }
    assert!(matches!(
        Candidate::new(
            "A",
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(-1.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(1.0),
        ),
        Err(CandidateError::InvalidNumber {
            field: "prior_variance",
            ..
        })
    ));
    assert!(matches!(
        Candidate::new(
            "A",
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(0.0),
            QtyAny::dimensionless(1.0),
        ),
        Err(CandidateError::InvalidNumber {
            field: "sensor_noise_variance",
            ..
        })
    ));
    assert!(matches!(
        Candidate::new(
            "A",
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(f64::NAN),
        ),
        Err(CandidateError::InvalidNumber {
            field: "sensor_cost",
            ..
        })
    ));

    let zero = candidate("zero", -0.0, -0.0, -0.0, 1.0, 1.0);
    assert_eq!(zero.truth().to_bits(), 0.0_f64.to_bits());
    assert_eq!(zero.prior_mean().to_bits(), 0.0_f64.to_bits());
    assert_eq!(zero.prior_variance().value.to_bits(), 0.0_f64.to_bits());
}

#[test]
fn candidate_refuses_wrong_variance_noise_and_cost_dimensions() {
    let length = Dims([1, 0, 0, 0, 0, 0]);
    let area = QtyAny::new(1.0, length)
        .powi(2)
        .expect("length squared is representable")
        .dims;
    let truth = ObjectiveValue::dimensional(QtyAny::new(1.0, length)).expect("finite length");
    let mean = ObjectiveValue::dimensional(QtyAny::new(1.0, length)).expect("finite length");

    assert!(matches!(
        Candidate::new(
            "A",
            truth,
            mean,
            QtyAny::new(0.1, length),
            QtyAny::new(0.1, area),
            QtyAny::dimensionless(1.0),
        ),
        Err(CandidateError::DimensionMismatch {
            field: "prior_variance",
            actual,
            expected,
        }) if actual == length && expected == area
    ));
    assert!(matches!(
        Candidate::new(
            "A",
            truth,
            mean,
            QtyAny::new(0.1, area),
            QtyAny::new(0.1, length),
            QtyAny::dimensionless(1.0),
        ),
        Err(CandidateError::DimensionMismatch {
            field: "sensor_noise_variance",
            actual,
            expected,
        }) if actual == length && expected == area
    ));
    assert!(matches!(
        Candidate::new(
            "A",
            truth,
            mean,
            QtyAny::new(0.1, area),
            QtyAny::new(0.1, area),
            QtyAny::new(1.0, length),
        ),
        Err(CandidateError::DimensionMismatch {
            field: "sensor_cost",
            actual,
            expected,
        }) if actual == length && expected == Dims::NONE
    ));
}

#[test]
fn acquisition_cost_roles_are_sealed_and_invalid_values_refuse() {
    let relative = AcquisitionCost::relative(1.0).expect("positive relative weight");
    let dimensionless_quantity =
        AcquisitionCost::dimensional(QtyAny::dimensionless(1.0)).expect("positive quantity");
    assert_eq!(relative.spec(), AcquisitionCostSpec::RelativeWeight);
    assert_eq!(
        dimensionless_quantity.spec(),
        AcquisitionCostSpec::Quantity(QuantitySpec::dimensional(Dims::NONE))
    );
    assert_ne!(relative.spec(), dimensionless_quantity.spec());

    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            AcquisitionCost::relative(invalid),
            Err(CandidateError::InvalidNumber {
                field: "acquisition_cost",
                ..
            })
        ));
        assert!(matches!(
            AcquisitionCost::dimensional(QtyAny::dimensionless(invalid)),
            Err(CandidateError::InvalidNumber {
                field: "acquisition_cost",
                ..
            })
        ));
    }

    let time_cost = AcquisitionCost::dimensional(parse_qty("1s").expect("time cost"))
        .expect("positive time cost");
    assert!(matches!(
        Candidate::with_acquisition_cost(
            "A",
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(1.0),
            time_cost,
            parse_qty("1s").expect("wrong allocation units"),
        ),
        Err(CandidateError::DimensionMismatch {
            field: "allocation_cost_weight",
            actual,
            expected,
        }) if actual == parse_qty("1s").expect("time dimensions").dims
            && expected == Dims::NONE
    ));
}

#[test]
fn campaign_refuses_acquisition_cost_role_and_semantic_aliases() {
    let relative = candidate("A", 0.0, 0.0, 1.0, 1.0, 1.0);
    let dimensionless_quantity = candidate_with_acquisition_cost(
        "B",
        objective(0.1),
        objective(0.1),
        QtyAny::dimensionless(1.0),
        QtyAny::dimensionless(1.0),
        AcquisitionCost::dimensional(QtyAny::dimensionless(1.0))
            .expect("dimensionless physical quantity"),
        1.0,
    );
    assert!(matches!(
        campaign(&[relative, dimensionless_quantity], 0.0, 0),
        Err(OedError::AcquisitionCostSchemaMismatch {
            candidate,
            actual: AcquisitionCostSpec::Quantity(QuantitySpec::Dimensional(Dims::NONE)),
            expected: AcquisitionCostSpec::RelativeWeight,
        }) if candidate == "B"
    ));

    let semantic_candidate = |name, kind| {
        candidate_with_acquisition_cost(
            name,
            objective(0.0),
            objective(0.0),
            QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(1.0),
            semantic_acquisition_cost(1.0, kind),
            1.0,
        )
    };
    let energy = semantic_candidate("energy", QuantityKind::Energy);
    let torque = semantic_candidate("torque", QuantityKind::Torque);
    assert_eq!(
        energy.acquisition_cost().quantity().dims,
        torque.acquisition_cost().quantity().dims
    );
    assert!(matches!(
        campaign(&[energy.clone(), torque], 0.0, 0),
        Err(OedError::AcquisitionCostSchemaMismatch { candidate, .. })
            if candidate == "torque"
    ));

    let time = candidate_with_acquisition_cost(
        "time",
        objective(0.0),
        objective(0.0),
        QtyAny::dimensionless(1.0),
        QtyAny::dimensionless(1.0),
        AcquisitionCost::dimensional(parse_qty("1s").expect("time cost"))
            .expect("positive time cost"),
        1.0,
    );
    assert!(matches!(
        campaign(&[energy.clone(), time], 0.0, 0),
        Err(OedError::AcquisitionCostSchemaMismatch {
            candidate,
            actual,
            expected,
        }) if candidate == "time"
            && actual.dims() == parse_qty("1s").expect("time dimensions").dims
            && expected == energy.acquisition_cost_spec()
    ));

    let dimension_only_energy = candidate_with_acquisition_cost(
        "dimension-only-energy",
        objective(0.0),
        objective(0.0),
        QtyAny::dimensionless(1.0),
        QtyAny::dimensionless(1.0),
        AcquisitionCost::dimensional(energy.acquisition_cost().quantity())
            .expect("positive dimension-only energy"),
        1.0,
    );
    assert!(matches!(
        campaign(&[energy, dimension_only_energy], 0.0, 0),
        Err(OedError::AcquisitionCostSchemaMismatch { candidate, .. })
            if candidate == "dimension-only-energy"
    ));
}

#[test]
fn acquisition_cost_units_normalize_and_role_changes_move_v9_identity() {
    let menu_with = |cost: AcquisitionCost| {
        vec![
            candidate_with_acquisition_cost(
                "A",
                objective(0.60),
                objective(0.60),
                QtyAny::dimensionless(0.10),
                QtyAny::dimensionless(0.01),
                cost,
                1.0,
            ),
            candidate_with_acquisition_cost(
                "B",
                objective(0.65),
                objective(0.65),
                QtyAny::dimensionless(0.12),
                QtyAny::dimensionless(0.01),
                cost,
                1.0,
            ),
        ]
    };
    let second =
        AcquisitionCost::dimensional(parse_qty("1s").expect("seconds")).expect("positive second");
    let thousand_ms = AcquisitionCost::dimensional(parse_qty("1000ms").expect("milliseconds"))
        .expect("positive milliseconds");
    assert_eq!(second, thousand_ms);
    let seconds_report = campaign(&menu_with(second), 0.0, 0).expect("seconds campaign");
    let milliseconds_report =
        campaign(&menu_with(thousand_ms), 0.0, 0).expect("milliseconds campaign");
    assert_eq!(seconds_report, milliseconds_report);

    let relative_report = campaign(
        &[
            candidate("A", 0.60, 0.60, 0.10, 0.01, 1.0),
            candidate("B", 0.65, 0.65, 0.12, 0.01, 1.0),
        ],
        0.0,
        0,
    )
    .expect("relative-weight campaign");
    let dimensionless_report = campaign(
        &menu_with(
            AcquisitionCost::dimensional(QtyAny::dimensionless(1.0))
                .expect("dimensionless declared quantity"),
        ),
        0.0,
        0,
    )
    .expect("dimensionless-quantity campaign");
    assert_same_non_identity_report(&relative_report, &dimensionless_report);
    assert_ne!(
        relative_report.acquisition_cost_spec(),
        dimensionless_report.acquisition_cost_spec()
    );
    assert_ne!(
        estimator(relative_report.variance_color()),
        estimator(dimensionless_report.variance_color())
    );
}

#[test]
fn physical_cost_derives_score_schema_and_never_becomes_allocator_weight() {
    let absolute = semantic_objective(300.0, QuantityKind::AbsoluteTemperature);
    let variance_dims = absolute
        .quantity()
        .powi(2)
        .expect("temperature variance dimensions")
        .dims;
    let energy_cost = semantic_acquisition_cost(2.0, QuantityKind::Energy);
    let temperature = candidate_with_acquisition_cost(
        "T",
        absolute,
        absolute,
        QtyAny::new(1.0, variance_dims),
        QtyAny::new(1.0, variance_dims),
        energy_cost,
        1.0,
    );
    let temperature_report = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                std::slice::from_ref(&temperature),
                temperature
                    .objective_spec()
                    .decision_value(0.0)
                    .expect("delta-temperature threshold"),
                0,
                cx,
            )
        },
    )
    .expect("typed physical-cost campaign");
    assert_eq!(
        temperature_report.acquisition_cost_spec(),
        AcquisitionCostSpec::Quantity(QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Energy,
            ValueForm::Static,
        )))
    );
    assert_eq!(
        temperature_report.score_spec().dims(),
        temperature
            .objective_spec()
            .decision_spec()
            .dims()
            .checked_minus(energy_cost.quantity().dims)
            .expect("temperature per energy dimensions")
    );
    assert!(temperature_report.score_spec().semantic_type().is_none());

    let make_menu = |first_cost: &str, first_weight: f64| {
        vec![
            candidate_with_acquisition_cost(
                "A",
                objective(0.0),
                objective(0.0),
                QtyAny::dimensionless(1.0),
                QtyAny::dimensionless(1.0),
                AcquisitionCost::dimensional(parse_qty(first_cost).expect("time cost"))
                    .expect("positive time cost"),
                first_weight,
            ),
            candidate_with_acquisition_cost(
                "B",
                objective(0.0),
                objective(0.0),
                QtyAny::dimensionless(4.0),
                QtyAny::dimensionless(1.0),
                AcquisitionCost::dimensional(parse_qty("1s").expect("time cost"))
                    .expect("positive time cost"),
                1.0,
            ),
        ]
    };
    let one_second = campaign(&make_menu("1s", 1.0), 0.0, 0).expect("one-second cost");
    let hundred_seconds = campaign(&make_menu("100s", 1.0), 0.0, 0).expect("hundred-second cost");
    assert_same_non_identity_report(&one_second, &hundred_seconds);

    let heavier_allocation =
        campaign(&make_menu("1s", 8.0), 0.0, 0).expect("heavier allocator coefficient");
    assert_ne!(one_second.allocation(), heavier_allocation.allocation());
}

#[test]
fn score_dimension_overflow_refuses_before_campaign_science() {
    let objective_dims = Dims([-64, 0, 0, 0, 0, 0]);
    let variance_dims = Dims([-128, 0, 0, 0, 0, 0]);
    let cost_dims = Dims([127, 0, 0, 0, 0, 0]);
    let candidate = candidate_with_acquisition_cost(
        "overflow",
        ObjectiveValue::dimensional(QtyAny::new(1.0, objective_dims)).expect("finite objective"),
        ObjectiveValue::dimensional(QtyAny::new(1.0, objective_dims)).expect("finite objective"),
        QtyAny::new(1.0, variance_dims),
        QtyAny::new(1.0, variance_dims),
        AcquisitionCost::dimensional(QtyAny::new(1.0, cost_dims)).expect("positive cost"),
        1.0,
    );
    let refusal = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                std::slice::from_ref(&candidate),
                ObjectiveValue::dimensional(QtyAny::new(0.0, objective_dims))
                    .expect("finite threshold"),
                0,
                cx,
            )
        },
    );
    assert_eq!(
        refusal,
        Err(OedError::ScoreDimensionOverflow {
            benefit: objective_dims,
            cost: cost_dims,
        })
    );
}

#[test]
fn campaign_refuses_semantic_aliases_and_affine_temperature_thresholds() {
    let pressure = semantic_objective(1.0, QuantityKind::Pressure);
    let stress = semantic_objective(1.0, QuantityKind::Stress);
    let dimension_only_pressure = ObjectiveValue::dimensional(pressure.quantity())
        .expect("finite dimension-only pressure dimensions");
    assert_eq!(pressure.quantity().dims, stress.quantity().dims);
    assert_ne!(pressure.spec(), stress.spec());
    assert_ne!(pressure.spec(), dimension_only_pressure.spec());
    let variance_dims = pressure
        .quantity()
        .powi(2)
        .expect("pressure squared is representable")
        .dims;
    let pressure_candidate = typed_candidate(
        "A",
        pressure,
        pressure,
        QtyAny::new(0.1, variance_dims),
        QtyAny::new(0.1, variance_dims),
    );
    let pressure_threshold = pressure
        .spec()
        .decision_value(0.01)
        .expect("finite pressure decision threshold");
    let stress_candidate = typed_candidate(
        "B",
        stress,
        stress,
        QtyAny::new(0.1, variance_dims),
        QtyAny::new(0.1, variance_dims),
    );
    let mixed = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                &[pressure_candidate.clone(), stress_candidate],
                pressure_threshold,
                0,
                cx,
            )
        },
    );
    assert!(matches!(
        mixed,
        Err(OedError::ObjectiveSchemaMismatch {
            candidate,
            actual,
            expected,
        }) if candidate == "B" && actual == stress.spec() && expected == pressure.spec()
    ));
    let dimension_only_candidate = typed_candidate(
        "B",
        dimension_only_pressure,
        dimension_only_pressure,
        QtyAny::new(0.1, variance_dims),
        QtyAny::new(0.1, variance_dims),
    );
    let untyped_alias = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                &[pressure_candidate.clone(), dimension_only_candidate],
                pressure_threshold,
                0,
                cx,
            )
        },
    );
    assert!(matches!(
        untyped_alias,
        Err(OedError::ObjectiveSchemaMismatch { .. })
    ));
    let wrong_threshold = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                std::slice::from_ref(&pressure_candidate),
                semantic_objective(0.01, QuantityKind::Stress),
                0,
                cx,
            )
        },
    );
    assert!(matches!(
        wrong_threshold,
        Err(OedError::ThresholdSchemaMismatch { .. })
    ));

    let menu_for = |kind| {
        let a = semantic_objective(0.60, kind);
        let b = semantic_objective(0.65, kind);
        let squared = a
            .quantity()
            .powi(2)
            .expect("semantic objective square")
            .dims;
        vec![
            typed_candidate(
                "A",
                a,
                a,
                QtyAny::new(0.10, squared),
                QtyAny::new(0.01, squared),
            ),
            typed_candidate(
                "B",
                b,
                b,
                QtyAny::new(0.12, squared),
                QtyAny::new(0.01, squared),
            ),
        ]
    };
    let pressure_menu = menu_for(QuantityKind::Pressure);
    let stress_menu = menu_for(QuantityKind::Stress);
    let pressure_report = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                &pressure_menu,
                pressure_menu[0]
                    .objective_spec()
                    .decision_value(0.0)
                    .expect("finite pressure decision threshold"),
                1,
                cx,
            )
        },
    )
    .expect("pressure campaign");
    let stress_report = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                &stress_menu,
                stress_menu[0]
                    .objective_spec()
                    .decision_value(0.0)
                    .expect("finite stress decision threshold"),
                1,
                cx,
            )
        },
    )
    .expect("stress campaign");
    assert_ne!(
        estimator(pressure_report.variance_color()),
        estimator(stress_report.variance_color())
    );
    assert_eq!(pressure_report.assimilation_colors().len(), 1);
    assert_eq!(stress_report.assimilation_colors().len(), 1);
    assert_ne!(
        estimator(&pressure_report.assimilation_colors()[0]),
        estimator(&stress_report.assimilation_colors()[0]),
        "nested assimilation identity must bind the objective schema token"
    );
    assert!(
        pressure_report
            .initial_evpi()
            .spec()
            .semantic_type()
            .is_none()
    );
    assert!(
        stress_report
            .initial_evpi()
            .spec()
            .semantic_type()
            .is_none()
    );

    let absolute = semantic_objective(300.0, QuantityKind::AbsoluteTemperature);
    let temperature_variance_dims = absolute
        .quantity()
        .powi(2)
        .expect("temperature squared is representable")
        .dims;
    let temperature = typed_candidate(
        "T",
        absolute,
        absolute,
        QtyAny::new(1.0, temperature_variance_dims),
        QtyAny::new(1.0, temperature_variance_dims),
    );
    let absolute_threshold = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                std::slice::from_ref(&temperature),
                semantic_objective(0.5, QuantityKind::AbsoluteTemperature),
                0,
                cx,
            )
        },
    );
    assert!(matches!(
        absolute_threshold,
        Err(OedError::ThresholdSchemaMismatch { actual, expected })
            if actual.semantic_type().is_some_and(|ty| ty.kind() == QuantityKind::AbsoluteTemperature)
                && expected.semantic_type().is_some_and(|ty| ty.kind() == QuantityKind::TemperatureDifference)
    ));
    let delta_threshold = with_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| {
            run_campaign(
                std::slice::from_ref(&temperature),
                temperature
                    .objective_spec()
                    .decision_value(0.5)
                    .expect("finite temperature-difference threshold"),
                0,
                cx,
            )
        },
    )
    .expect("delta-temperature threshold matches an absolute-temperature objective");
    assert_eq!(
        delta_threshold
            .decision_spec()
            .semantic_type()
            .expect("semantic decision schema")
            .kind(),
        QuantityKind::TemperatureDifference
    );
}

#[test]
fn semantic_measurements_do_not_launder_their_kind_or_form_into_decision_loss() {
    let composition = semantic_objective(
        0.4,
        QuantityKind::Composition(CompositionBasis::MassFraction),
    );
    let decision = composition
        .spec()
        .decision_value(0.05)
        .expect("finite composition decision threshold");

    assert_eq!(decision.quantity().dims, composition.quantity().dims);
    assert!(composition.spec().semantic_type().is_some());
    assert!(decision.spec().semantic_type().is_none());
}

#[test]
fn g3_metre_millimetre_rescaling_preserves_science_identity_and_output_units() {
    let metre_candidates = vec![
        typed_candidate(
            "A",
            ObjectiveValue::dimensional(parse_qty("1m").expect("metres")).expect("finite"),
            ObjectiveValue::dimensional(parse_qty("1m").expect("metres")).expect("finite"),
            parse_qty("0m^2").expect("square metres"),
            parse_qty("1m^2").expect("square metres"),
        ),
        typed_candidate(
            "B",
            ObjectiveValue::dimensional(parse_qty("2m").expect("metres")).expect("finite"),
            ObjectiveValue::dimensional(parse_qty("2m").expect("metres")).expect("finite"),
            parse_qty("0m^2").expect("square metres"),
            parse_qty("1m^2").expect("square metres"),
        ),
    ];
    let millimetre_candidates = vec![
        typed_candidate(
            "A",
            ObjectiveValue::dimensional(parse_qty("1000mm").expect("millimetres")).expect("finite"),
            ObjectiveValue::dimensional(parse_qty("1000mm").expect("millimetres")).expect("finite"),
            parse_qty("0mm^2").expect("square millimetres"),
            parse_qty("1000000mm^2").expect("square millimetres"),
        ),
        typed_candidate(
            "B",
            ObjectiveValue::dimensional(parse_qty("2000mm").expect("millimetres")).expect("finite"),
            ObjectiveValue::dimensional(parse_qty("2000mm").expect("millimetres")).expect("finite"),
            parse_qty("0mm^2").expect("square millimetres"),
            parse_qty("1000000mm^2").expect("square millimetres"),
        ),
    ];
    let run = |candidates: &[Candidate], threshold: &str| {
        with_cx(
            &CancelGate::new(),
            Budget::INFINITE,
            ExecMode::Deterministic,
            |cx| {
                run_campaign(
                    candidates,
                    ObjectiveValue::dimensional(parse_qty(threshold).expect("threshold units"))
                        .expect("finite threshold"),
                    0,
                    cx,
                )
            },
        )
        .expect("unit-rescaled campaign")
    };
    let metres = run(&metre_candidates, "0.1m");
    let millimetres = run(&millimetre_candidates, "100mm");
    assert_eq!(metres, millimetres);

    let length = Dims([1, 0, 0, 0, 0, 0]);
    let area = Dims([2, 0, 0, 0, 0, 0]);
    assert_eq!(metres.objective_spec().dims(), length);
    assert_eq!(metres.initial_evpi().quantity().dims, length);
    assert_eq!(metres.final_evpi().quantity().dims, length);
    assert_eq!(metres.prior_total_variance().dims, area);
    assert_eq!(metres.posterior_total_variance().dims, area);
    assert!(
        metres
            .posteriors()
            .iter()
            .all(|posterior| posterior.mean().quantity().dims == length
                && posterior.variance().dims == area)
    );
}

#[test]
fn g3_rescaling_schema_is_scale_independent_while_value_bits_remain_exact() {
    let metres = ObjectiveValue::dimensional(parse_qty("0.7m").expect("metres")).expect("finite");
    let millimetres =
        ObjectiveValue::dimensional(parse_qty("700mm").expect("millimetres")).expect("finite");

    assert_eq!(metres.spec(), millimetres.spec());
    assert_eq!(
        metres.spec().quantity_spec(),
        millimetres.spec().quantity_spec()
    );
    assert_eq!(metres.quantity().dims, millimetres.quantity().dims);
    assert_ne!(metres.to_bits(), millimetres.to_bits());
    assert_eq!(metres.to_bits(), 0x3fe6_6666_6666_6666);
    assert_eq!(millimetres.to_bits(), 0x3fe6_6666_6666_6667);
}

#[test]
fn campaign_inputs_are_checked_before_work_starts() {
    assert_eq!(campaign(&[], 0.0, 0), Err(OedError::NoCandidates));
    let one = candidate("A", 0.0, 0.0, 1.0, 1.0, 1.0);
    for threshold in [f64::NAN, f64::INFINITY] {
        assert!(ObjectiveValue::dimensionless(threshold).is_err());
    }
    assert_eq!(
        campaign(std::slice::from_ref(&one), -1.0, 0),
        Err(OedError::InvalidThreshold)
    );
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
    assert_eq!(
        report.prior_total_variance().value.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        report.posterior_total_variance().value.to_bits(),
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
        d.truth().value() + 0.01,
        d.prior_mean().value(),
        d.prior_variance().value,
        d.sensor_noise_variance().value,
        d.sensor_cost().value,
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
    assert!(estimator(report.variance_color()).starts_with("sensorforge-posterior-variance:v9:"));
    assert!(estimator(report.evpi_color()).starts_with("sensorforge-evpi:v9:"));
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
    assert!(report.initial_evpi().value() > 0.0);
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
        run_campaign(&candidates, objective(0.01), 12, cx)
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
        run_campaign(&candidates, objective(0.0), 12, cx)
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
            run_campaign(&candidates, objective(0.0), 1, cx)
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
            |cx| run_campaign(&candidates, objective(0.0), 1, cx),
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
            run_campaign(&candidates, objective(1.0e9), 0, cx)
        })
        .expect("zero-placement campaign succeeds")
    };
    let deterministic = run(Budget::INFINITE, ExecMode::Deterministic, TEST_STREAM);
    let fast = run(Budget::INFINITE, ExecMode::Fast, TEST_STREAM);

    let mut deadline_budget = Budget::INFINITE;
    // A zero deadline can never admit under the sj31i.6 accountant; a
    // future deadline still mutates the identity hash the same way.
    deadline_budget.deadline = Budget::with_deadline_at_ns(123_456).deadline;
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
        run_campaign(&accepted, objective(0.0), MAX_CAMPAIGN_SENSORS, cx)
    });
    assert!(matches!(at_limit, Err(OedError::Cancelled { .. })));
    let bounded_later = with_cx(
        &CancelGate::new(),
        Budget::INFINITE.with_poll_quota(1),
        ExecMode::Deterministic,
        |cx| run_campaign(&accepted, objective(0.0), MAX_CAMPAIGN_SENSORS, cx),
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
            |cx| { run_campaign(&refused, objective(0.0), MAX_CAMPAIGN_SENSORS, cx,) },
        ),
        Err(OedError::WorkBudgetExceeded {
            candidates: 16,
            max_sensors: MAX_CAMPAIGN_SENSORS,
            evaluations: 16 * 16 * MAX_CAMPAIGN_SENSORS * 11,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
        })
    );
}

/// G0/G5 (bead sj31i.62): every bounded seam charges the deterministic
/// byte ledger. The demo campaign's admitted envelope, consumed seam
/// total, and retained subset are exact pinned values, ordered
/// retained <= consumed <= admitted, and replay bit-identically.
#[test]
fn byte_ledger_is_admitted_charged_and_retained_exactly() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let report = campaign(&candidates, 0.01, 12).expect("demo campaign succeeds");
    assert!(report.retained_byte_units() <= report.consumed_byte_units());
    assert!(report.consumed_byte_units() <= report.admitted_byte_units());
    assert_eq!(report.admitted_byte_units(), 9_002_023);
    assert_eq!(report.consumed_byte_units(), 7_504_799);
    assert_eq!(report.retained_byte_units(), 782);
    let replay = campaign(&candidates, 0.01, 12).expect("replay succeeds");
    assert_eq!(replay.admitted_byte_units(), report.admitted_byte_units());
    assert_eq!(replay.consumed_byte_units(), report.consumed_byte_units());
    assert_eq!(replay.retained_byte_units(), report.retained_byte_units());
}

/// G4 (bead sj31i.62): sweep the poll quota from zero until the demo
/// campaign completes. Every deterministic boundary — admission, the
/// action-value tiles whose quadrature override views are constructed
/// and read inside them, placement commit, the canonical-menu and EVPI
/// refresh seams, finalization scans, both report identities, and
/// publication — either refuses typed with a consistent work ledger or
/// publishes the full deterministic report. No quota can publish a
/// partial report, and quota sufficiency is monotone.
#[test]
fn g4_poll_quota_sweep_covers_every_boundary_without_partial_publication() {
    let candidates = demo_candidates().expect("compiled demo candidates are valid");
    let baseline = campaign(&candidates, 0.01, 12).expect("unbounded baseline succeeds");
    let mut first_sufficient = None;
    for quota in 0u32..=4_096 {
        let budget = Budget::new().with_poll_quota(quota);
        let outcome = with_cx(&CancelGate::new(), budget, ExecMode::Deterministic, |cx| {
            run_campaign(&candidates, objective(0.01), 12, cx)
        });
        match outcome {
            Ok(report) => {
                assert_same_non_identity_report(&report, &baseline);
                first_sufficient = Some(quota);
                break;
            }
            Err(OedError::Cancelled {
                completed_work_units,
                admitted_work_units,
                ..
            }) => {
                assert!(
                    completed_work_units <= admitted_work_units,
                    "quota {quota}: the refused ledger exceeded its admitted plan"
                );
            }
            Err(OedError::AssimilationCancelled {
                completed_work_units,
                admitted_work_units,
                ..
            }) => {
                assert!(
                    completed_work_units <= admitted_work_units,
                    "quota {quota}: the refused lower-layer ledger exceeded its plan"
                );
            }
            Err(other) => panic!("quota {quota}: expected a typed cancellation, got {other:?}"),
        }
    }
    let sufficient = first_sufficient.expect("the demo campaign completes within the sweep");
    for margin in [1u32, 7, 64] {
        let budget = Budget::new().with_poll_quota(sufficient + margin);
        let replay = with_cx(&CancelGate::new(), budget, ExecMode::Deterministic, |cx| {
            run_campaign(&candidates, objective(0.01), 12, cx)
        })
        .expect("a strictly larger poll quota cannot regress to refusal");
        assert_same_non_identity_report(&replay, &baseline);
    }
}
