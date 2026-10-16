use super::*;
use crate::evaluate_robust_sampled_stress;

fn cases() -> Vec<RobustLoadCase> {
    vec![
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ]
}

fn setup(cases: &[RobustLoadCase], aggregate: RobustAggregate) -> MultiLoadProjectedOptimizer {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = geometry.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
    MultiLoadProjectedOptimizer::new(
        geometry, cases,
        OptimizeSettings {
            level: 3, iterations: 3, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default()
        }, aggregate, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 6, ..MultiLoadProjectedSettings::default() },
    ).unwrap()
}

fn constrained() -> MultiLoadProjectedOptimizer {
    setup(&cases(), RobustAggregate::WeightedSum)
        .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap()
}

fn accepted(progress: MultiLoadProjectedProgress) -> Box<MultiLoadProjectedStep> {
    match progress {
        MultiLoadProjectedProgress::Accepted(step) => step,
        other => panic!("must exercise an accepted constrained update: {other:?}"),
    }
}

#[test]
fn cached_stress_equals_an_independent_full_resolve_without_hidden_solves() {
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        let mut optimizer = setup(&cases(), aggregate);
        let before = optimizer.current();
        let solves = optimizer.solves_started();
        let limit = SampledStressLimit::new(1e8, 0.0).unwrap();
        optimizer = optimizer.with_sampled_stress_limit(limit).unwrap();
        let replay = evaluate_robust_sampled_stress(
            optimizer.geometry(), &cases(), optimizer.kernel.settings, aggregate,
        ).unwrap();
        assert_eq!(optimizer.baseline_stress(), Some(&replay));
        assert_eq!(optimizer.current_stress(), Some(&replay));
        assert_eq!(optimizer.current(), before);
        assert_eq!(optimizer.solves_started(), solves);
        assert_eq!(optimizer.stress_limit(), Some(limit));
        assert!(replay.case_sample_counts.iter().all(|&count| count > 0));
    }
}

#[test]
fn zero_weight_scenario_can_control_and_refuse_the_baseline() {
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap(),
    ];
    let optimizer = setup(&cases, RobustAggregate::WorstWeightedCase);
    let reference = evaluate_robust_sampled_stress(
        optimizer.geometry(), &cases, optimizer.kernel.settings, RobustAggregate::WorstWeightedCase,
    ).unwrap();
    assert_eq!(reference.worst_stress_case, 1);
    assert_eq!(optimizer.current().active_case, Some(0));
    let low = reference.case_sampled_max_von_mises[0];
    assert!(low > 0.0);
    assert!((reference.worst_sampled_von_mises - 4.0 * low).abs() < 1e-8 * low);
    let error = optimizer.with_sampled_stress_limit(SampledStressLimit::new(2.0 * low, 0.0).unwrap())
        .err().expect("zero-weight overload must refuse");
    assert!(error.to_string().contains("case 1"));
}

#[test]
fn accepted_endpoint_preserves_stress_and_replays_under_all_loads() {
    let mut optimizer = constrained();
    let baseline = optimizer.baseline_stress().unwrap().clone();
    let step = accepted(optimizer.advance_one().unwrap());
    let replay = evaluate_robust_sampled_stress(
        optimizer.geometry(), &cases(), optimizer.kernel.settings, RobustAggregate::WeightedSum,
    ).unwrap();
    assert_eq!(step.stress, Some(replay.clone()));
    assert_eq!(optimizer.current_stress(), Some(&replay));
    assert_eq!(optimizer.baseline_stress(), Some(&baseline));
    assert_eq!(step.state.snapshot, replay.snapshot);
    assert_eq!(step.state.case_compliances, replay.case_compliances);
    assert!(step.state.objective < step.previous.objective);
    assert!((step.state.volume - 0.6).abs() <= 1e-4);
    assert!(replay.worst_sampled_von_mises <= optimizer.stress_limit().unwrap().admitted_max());
}

#[test]
fn deliberately_tightened_candidate_gate_cannot_publish_compliance_only_success() {
    let mut optimizer = constrained();
    let before = optimizer.current();
    let before_stress = optimizer.current_stress().unwrap().clone();
    // A deliberate internal mutation isolates the candidate gate. Public API
    // forbids replacing a constraint or admitting this overstressed baseline.
    optimizer.stress_limit = Some(SampledStressLimit::new(1e-30, 0.0).unwrap());
    let progress = optimizer.advance_one().unwrap();
    let MultiLoadProjectedProgress::NoDescent(attempts) = progress else {
        panic!("an impossible stress limit cannot admit the nominal descent");
    };
    assert!(attempts.iter().any(|attempt| attempt.stress.is_some()));
    for attempt in attempts.iter().filter(|attempt| attempt.stress.is_some()) {
        assert!(attempt.refusal.as_deref().unwrap().contains("sampled stress limit exceeded"));
        assert!(attempt.stress.as_ref().unwrap().worst_sampled_von_mises > 1e-30);
    }
    assert_eq!(optimizer.current(), before);
    assert_eq!(optimizer.current_stress(), Some(&before_stress));
    assert_eq!(optimizer.next_iteration(), 0);
    assert!(optimizer.solves_started() > 2);
}

#[test]
fn late_stress_sampling_cancellation_keeps_old_evidence_and_retry_matches() {
    let mut optimizer = constrained();
    let before = optimizer.current();
    let before_stress = optimizer.current_stress().unwrap().clone();
    let before_ell = optimizer.ell.to_bits();
    let stopped = optimizer.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::StressCell { case: 1, cell: 1, .. }) {
            ControlFlow::Break("late stress")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(stopped, ControlFlow::Break("late stress")));
    assert_eq!(optimizer.solves_started(), 4);
    assert_eq!(optimizer.next_iteration(), 0);
    assert_eq!(optimizer.current(), before);
    assert_eq!(optimizer.current_stress(), Some(&before_stress));
    assert_eq!(optimizer.ell.to_bits(), before_ell);
    let resumed = accepted(optimizer.advance_one().unwrap());
    let mut clean = constrained();
    let uninterrupted = accepted(clean.advance_one().unwrap());
    assert_eq!(resumed.state, uninterrupted.state);
    assert_eq!(resumed.stress, uninterrupted.stress);
    assert_eq!(optimizer.geometry().nodes(), clean.geometry().nodes());
    assert_eq!(optimizer.solves_started(), clean.solves_started() + 2);
}

#[test]
fn no_limit_retains_the_original_path_without_stress_work_or_claims() {
    let mut optimizer = setup(&cases(), RobustAggregate::WeightedSum);
    assert_eq!(optimizer.stress_limit(), None);
    assert_eq!(optimizer.current_stress(), None);
    let progress = optimizer.advance_one_controlled(|stage| -> ControlFlow<()> {
        assert!(!matches!(stage, MultiLoadProjectedStage::StressCell { .. }));
        ControlFlow::Continue(())
    }).unwrap();
    let ControlFlow::Continue(progress) = progress else { panic!("not cancelled") };
    let step = accepted(progress);
    assert!(step.stress.is_none());
    assert!(step.attempts.iter().all(|attempt| attempt.stress.is_none()));
}

#[test]
fn absent_nonfinite_repeated_or_late_constraint_evidence_refuses() {
    let limit = SampledStressLimit::new(1e8, 0.0).unwrap();
    assert!(require_feasible(Some(limit), None).is_err());
    let optimizer = constrained();
    let mut bad = optimizer.current_stress().unwrap().clone();
    bad.worst_sampled_von_mises = f64::NAN;
    assert!(require_feasible(Some(limit), Some(&bad)).is_err());
    assert!(optimizer.with_sampled_stress_limit(limit).is_err());
    assert!(setup(&cases(), RobustAggregate::WeightedSum)
        .with_sampled_stress_limit(SampledStressLimit { max_von_mises: f64::NAN, absolute_tolerance: 0.0 })
        .is_err());
    let mut late = setup(&cases(), RobustAggregate::WeightedSum);
    accepted(late.advance_one().unwrap());
    assert!(late.with_sampled_stress_limit(limit).is_err());
}
