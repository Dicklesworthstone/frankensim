use fs_uq::{CorrelationModel, EpistemicConfig, EpistemicExecution, EpistemicStatus,
    ParameterUncertainty, PropagationMethod, UqPlan};

fn plan(budget: usize) -> UqPlan {
    UqPlan::new("response", PropagationMethod::EpistemicBounding, budget)
        .with_parameter(ParameterUncertainty::interval("x", -1.0, 1.0, "m"))
}
fn config() -> EpistemicConfig {
    EpistemicConfig { absolute_tolerance: 1e-3, maximum_cells: 512 }
}
fn peak(x: &[[f64; 2]]) -> Result<[f64; 2], &'static str> {
    let [lo, hi] = x[0];
    let square_min = if lo <= 0.0 && hi >= 0.0 { 0.0 }
        else { (lo * lo).min(hi * hi).next_down().max(0.0) };
    let square_max = (lo * lo).max(hi * hi).next_up();
    Ok([(1.0 - square_max).next_down(), (1.0 - square_min).next_up()])
}
fn contains(bounds: [f64; 2], value: f64) {
    assert!(bounds[0] <= value && value <= bounds[1], "{bounds:?} omits {value}");
}

#[test]
fn interior_worst_case_is_not_replaced_by_zero_corner_samples() {
    let mut run = EpistemicExecution::new(&plan(4096), config()).unwrap();
    let result = run.advance(usize::MAX, || false, peak);
    assert_eq!(result.status, EpistemicStatus::Complete);
    contains(result.minimum.unwrap(), 0.0);
    contains(result.maximum.unwrap(), 1.0);
    for interval in [result.minimum.unwrap(), result.maximum.unwrap()] {
        assert!(interval[1] - interval[0] <= config().absolute_tolerance);
    }
    for i in 0..1001 {
        let x = -1.0 + 2.0 * i as f64 / 1000.0;
        contains(result.output_range.unwrap(), 1.0 - x * x);
    }
    let witness = result.maximizing_witness.unwrap();
    assert_eq!(witness.parameters, vec![0.0]);
    contains(witness.value, 1.0);
    assert_eq!(run.plan(), &plan(4096));
}

#[test]
fn interval_heat_load_propagation_retains_physical_parameter_order() {
    let plan = UqPlan::new("junction temperature", PropagationMethod::EpistemicBounding, 16_384)
        .with_parameter(ParameterUncertainty::interval("ambient", 20.0, 25.0, "K"))
        .with_parameter(ParameterUncertainty::interval("power", 2.0, 4.0, "W"))
        .with_parameter(ParameterUncertainty::interval("resistance", 1.0, 2.0, "K/W"));
    let config = EpistemicConfig { absolute_tolerance: 0.05, maximum_cells: 4096 };
    let mut run = EpistemicExecution::new(&plan, config).unwrap();
    let result = run.advance(usize::MAX, || false, |b| {
        let lo = (b[1][0] * b[2][0]).next_down();
        let hi = (b[1][1] * b[2][1]).next_up();
        Ok::<_, &str>([(b[0][0] + lo).next_down(), (b[0][1] + hi).next_up()])
    });
    assert_eq!(result.status, EpistemicStatus::Complete);
    contains(result.minimum.unwrap(), 22.0);
    contains(result.maximum.unwrap(), 33.0);
    assert!(result.minimum.unwrap()[1] - 22.0 <= 0.05);
    assert!(33.0 - result.maximum.unwrap()[0] <= 0.05);
    let w = result.maximizing_witness.unwrap();
    contains(w.value, w.parameters[0] + w.parameters[1] * w.parameters[2]);
    for (x, [lo, hi]) in w.parameters.into_iter().zip([[20.0, 25.0], [2.0, 4.0], [1.0, 2.0]]) {
        assert!(lo <= x && x <= hi);
    }
}

#[test]
fn paid_child_enclosures_survive_cancel_interrupt_and_clone_without_losing_parent() {
    let mut full = EpistemicExecution::new(&plan(4096), config()).unwrap();
    let mut expected = Vec::new();
    let complete = full.advance(usize::MAX, || false, |b| { expected.push(b.to_vec()); peak(b) });
    assert_eq!(complete.status, EpistemicStatus::Complete);
    let mut split = EpistemicExecution::new(&plan(4096), config()).unwrap();
    let mut paid = Vec::new();
    let partial = split.advance(3, || false, |b| { paid.push(b.to_vec()); peak(b) });
    assert!(partial.pending_split);
    assert_eq!(partial.covered_cells, 1);
    contains(partial.output_range.unwrap(), 0.0);
    contains(partial.output_range.unwrap(), 1.0);
    let cancelled = split.advance(1, || true, |_| -> Result<[f64; 2], &str> { panic!("cancelled") });
    assert_eq!(cancelled.status, EpistemicStatus::Cancelled);
    let unfinished = split.advance_interruptible(1, || false, |b| {
        assert_eq!(b, expected[3].as_slice());
        Ok::<_, &str>(None)
    });
    assert_eq!(unfinished.evaluations_attempted, 3);
    let mut resumed = split.clone();
    for allowance in [1, 2, 3, 5, usize::MAX] {
        resumed.advance(allowance, || false, |b| { paid.push(b.to_vec()); peak(b) });
    }
    assert_eq!(paid, expected);
    assert_eq!(resumed.report(), complete);
    assert_eq!(split.report(), unfinished);
    assert_eq!(resumed.advance(1, || panic!("terminal"),
        |_| -> Result<[f64; 2], &str> { panic!("terminal") }), complete);
    assert_eq!(split.advance(0, || panic!("zero allowance"),
        |_| -> Result<[f64; 2], &str> { panic!("zero allowance") }), unfinished);
}

#[test]
fn every_call_boundary_preserves_range_and_matches_one_shot_execution() {
    let mut reference = EpistemicExecution::new(&plan(512), config()).unwrap();
    let expected = reference.advance(usize::MAX, || false, peak);
    let mut run = EpistemicExecution::new(&plan(512), config()).unwrap();
    let mut prior = [f64::NEG_INFINITY, f64::INFINITY];
    for _ in 0..512 {
        let report = run.advance(1, || false, peak);
        let range = report.output_range.unwrap();
        contains(range, 0.0); contains(range, 1.0);
        assert!(range[0] >= prior[0] && range[1] <= prior[1]);
        prior = range;
        if report.status != EpistemicStatus::Paused { break; }
    }
    assert_eq!(run.report(), expected);
}

#[test]
fn lifetime_and_cell_limits_publish_coverage_not_false_completion() {
    for budget in [1, 3, 4, 5] {
        let mut run = EpistemicExecution::new(&plan(budget), config()).unwrap();
        let report = run.advance(100, || false, peak);
        assert_eq!(report.status, EpistemicStatus::EvaluationLimit);
        assert_eq!(report.evaluations_attempted, budget);
        contains(report.output_range.unwrap(), 0.0); contains(report.output_range.unwrap(), 1.0);
        if budget == 1 { assert!(report.maximizing_witness.is_none()); }
    }
    let mut run = EpistemicExecution::new(&plan(100), EpistemicConfig {
        maximum_cells: 1, ..config()
    }).unwrap();
    let report = run.advance(100, || false, peak);
    assert_eq!(report.status, EpistemicStatus::CellLimit);
    assert_eq!(report.evaluations_attempted, 2);
}

#[test]
fn singleton_and_adjacent_float_domains_do_not_invent_resolution() {
    let mut fixed = plan(16);
    fixed.parameters[0] = ParameterUncertainty::interval("x", 2.0, 2.0, "m");
    let mut run = EpistemicExecution::new(&fixed, config()).unwrap();
    let report = run.advance(16, || false, |b| Ok::<_, &str>(b[0]));
    assert_eq!(report.status, EpistemicStatus::Complete);
    assert_eq!(report.evaluations_attempted, 1);
    assert_eq!(report.minimum, Some([2.0; 2]));
    let mut run = EpistemicExecution::new(&fixed, config()).unwrap();
    let report = run.advance(16, || false, |_| Ok::<_, &str>([1.0, 3.0]));
    assert_eq!(report.status, EpistemicStatus::ParameterResolution);
    fixed.parameters[0] = ParameterUncertainty::interval("x", 1.0, 1.0_f64.next_up(), "m");
    let mut run = EpistemicExecution::new(&fixed, EpistemicConfig {
        absolute_tolerance: 1e-30, ..config()
    }).unwrap();
    let report = run.advance(16, || false, |b| Ok::<_, &str>(b[0]));
    assert_eq!(report.status, EpistemicStatus::ParameterResolution);
    contains(report.maximum.unwrap(), 1.0_f64.next_up());
    fixed.parameters[0] = ParameterUncertainty::interval("x", -f64::MAX, f64::MAX, "m");
    let mut run = EpistemicExecution::new(&fixed, config()).unwrap();
    let report = run.advance(2, || false, |b| Ok::<_, &str>(b[0]));
    assert_eq!(report.maximizing_witness.unwrap().parameters, vec![0.0]);
    assert_eq!(report.output_range, Some([-f64::MAX, f64::MAX]));
}

#[test]
fn invalid_nested_enclosures_and_model_failures_are_terminal_without_claims() {
    for invalid in [[f64::NAN, 1.0], [0.0, f64::INFINITY], [1.0, -1.0], [2.0, 3.0]] {
        let mut run = EpistemicExecution::new(&plan(16), config()).unwrap();
        run.advance(1, || false, |_| Ok::<_, &str>([-1.0, 1.0]));
        let failed = run.advance(1, || false, |_| Ok::<_, &str>(invalid));
        assert_eq!(failed.status, EpistemicStatus::Refused);
        assert_eq!(failed.evaluations_attempted, 2);
        assert_eq!(failed.evaluations_accepted, 1);
        assert!(failed.output_range.is_none() && failed.maximum.is_none());
        assert!(failed.rejection_reason.is_some());
        assert_eq!(run.advance(1, || panic!("terminal"),
            |_| -> Result<[f64; 2], &str> { panic!("terminal") }), failed);
    }
    let mut run = EpistemicExecution::new(&plan(16), config()).unwrap();
    let failed = run.advance(1, || false, |_| Err::<[f64; 2], _>("solver refused"));
    assert!(failed.rejection_reason.unwrap().contains("solver refused"));
}

#[test]
fn probability_laws_and_malformed_domains_refuse_before_callbacks() {
    let mut p = plan(16);
    for n in [0, 65_537] {
        p.budget_max_samples = n;
        assert!(EpistemicExecution::new(&p, config()).is_err());
    }
    p = plan(16);
    for kind in [ParameterUncertainty::gaussian("x", 0.0, 1.0, "m"),
        ParameterUncertainty::interval("x", 2.0, 1.0, "m"),
        ParameterUncertainty::interval("x", f64::NEG_INFINITY, 1.0, "m"),
        ParameterUncertainty::interval("x", -1.0, 1.0, "")] {
        p.parameters[0] = kind;
        assert!(EpistemicExecution::new(&p, config()).is_err());
    }
    p = plan(16);
    p.parameters.push(p.parameters[0].clone());
    assert!(EpistemicExecution::new(&p, config()).is_err());
    p = plan(16);
    p.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0]] };
    assert!(EpistemicExecution::new(&p, config()).is_err());
    p = plan(16); p.method = PropagationMethod::MonteCarlo;
    assert!(EpistemicExecution::new(&p, config()).is_err());
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(EpistemicExecution::new(&plan(16), EpistemicConfig {
            absolute_tolerance: tolerance, ..config()
        }).is_err());
    }
}
