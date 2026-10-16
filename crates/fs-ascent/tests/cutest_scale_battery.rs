//! Bead 7tv.21.13: dimension-scaled CUTEst-class gradient-stack battery.
//!
//! The landed oracle tranche covers small two-variable fixtures. This battery
//! adds three scalable, smooth families from the standard optimization test
//! corpus at 16 and 64 variables. Every analytic gradient crosses an
//! independent `fsci_opt::check_grad` finite-difference gate before L-BFGS is
//! trusted, and a one-coordinate mutant proves that the gate is live. The
//! optimizer then has to reach the known zero objective inside an explicit
//! evaluation ceiling and reproduce every evaluated point bit for bit.
//!
//! These are deterministic G0/G5 conformance rows, not wall-clock or complete
//! CUTEst-distribution coverage. No external CUTEst runtime is used.

use core::cell::{Cell, RefCell};

use fs_ascent::{LbfgsReport, LbfgsState, StopReason, StopRule};
use fs_obs::ident::IdentityBuilder;
use fs_obs::{Emitter, EventKind, Severity};
use fsci_opt::check_grad;

const SUITE: &str = "fs-ascent/cutest-scale";
const FIXED_INPUT_SEED: u64 = 0;
const DIMENSIONS: [usize; 2] = [16, 64];
const MUTANT_FLOOR: f64 = 5e-3;
const MAX_ITERS: usize = 10_000;

type FnGrad = fn(&[f64]) -> (f64, Vec<f64>);
type PointFactory = fn(usize) -> Vec<f64>;

#[derive(Clone, Copy)]
struct Fixture {
    name: &'static str,
    fg: FnGrad,
    start: PointFactory,
    gradient_probe: PointFactory,
    gradient_tolerance: f64,
    objective_target: f64,
    evaluation_ceiling: usize,
}

struct RunRecord {
    state: LbfgsState,
    report: LbfgsReport,
    calls: usize,
    evaluated_points: Vec<Vec<u64>>,
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn emit_verdict(case: &str, pass: bool, detail: &str) {
    let mut emitter = Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed: FIXED_INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("failure row must reproduce from its detail");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("CUTEst-scale row must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "{case}: {detail}");
}

fn run_once(fixture: Fixture, dimension: usize) -> RunRecord {
    let start = (fixture.start)(dimension);
    assert_eq!(start.len(), dimension, "fixture start dimension drift");

    let calls = Cell::new(0usize);
    let evaluated_points = RefCell::new(Vec::new());
    let (state, report) = {
        let mut counted = |x: &[f64]| {
            calls.set(calls.get() + 1);
            evaluated_points.borrow_mut().push(bits(x));
            (fixture.fg)(x)
        };
        let mut state = LbfgsState::new(&start, 17, &mut counted);
        let stop = StopRule::Any(vec![
            StopRule::ObjectiveBelow(fixture.objective_target),
            StopRule::Budget(fixture.evaluation_ceiling),
        ]);
        let report = state.run(&mut counted, &stop, MAX_ITERS);
        (state, report)
    };

    RunRecord {
        state,
        report,
        calls: calls.get(),
        evaluated_points: evaluated_points.into_inner(),
    }
}

fn public_replay_equal(first: &RunRecord, repeat: &RunRecord) -> bool {
    first.evaluated_points == repeat.evaluated_points
        && bits(&first.state.x) == bits(&repeat.state.x)
        && first.state.f.to_bits() == repeat.state.f.to_bits()
        && bits(&first.state.g) == bits(&repeat.state.g)
        && bits(&first.state.history) == bits(&repeat.state.history)
        && first.state.iters == repeat.state.iters
        && first.state.evals == repeat.state.evals
        && first.calls == repeat.calls
        && first.report.reason == repeat.report.reason
        && first.report.grad_norm.to_bits() == repeat.report.grad_norm.to_bits()
        && first.report.f.to_bits() == repeat.report.f.to_bits()
        && first.report.iters == repeat.report.iters
        && first.report.evals == repeat.report.evals
}

fn stop_reason_name(reason: &StopReason) -> &'static str {
    match reason {
        StopReason::GradNorm => "gradient-norm",
        StopReason::ObjectiveBelow => "objective-below",
        StopReason::Budget => "budget",
        StopReason::Stall => "stall",
        StopReason::Composite => "composite",
        StopReason::IterationCap => "iteration-cap",
    }
}

fn trace_identity(fixture: Fixture, dimension: usize, run: &RunRecord) -> String {
    let mut trace_bytes = Vec::new();
    for point in &run.evaluated_points {
        assert_eq!(point.len(), dimension, "evaluated point dimension drift");
        for &coordinate_bits in point {
            trace_bytes.extend_from_slice(&coordinate_bits.to_le_bytes());
        }
    }
    let mut builder = IdentityBuilder::new("fs-ascent-cutest-scale-evaluation-trace-v1")
        .str("fixture", fixture.name)
        .str("engine", "fs-ascent/L-BFGS")
        .str("stop-reason", stop_reason_name(&run.report.reason))
        .str("fs-ascent-version", fs_ascent::VERSION)
        .u64("fixed-input-seed", FIXED_INPUT_SEED)
        .u64(
            "dimension",
            u64::try_from(dimension).expect("dimension fits u64"),
        )
        .u64("memory-pairs", 17)
        .u64(
            "maximum-iterations",
            u64::try_from(MAX_ITERS).expect("iteration cap fits u64"),
        )
        .u64(
            "evaluation-ceiling",
            u64::try_from(fixture.evaluation_ceiling).expect("evaluation ceiling fits u64"),
        )
        .f64_bits("objective-target", fixture.objective_target)
        .u64(
            "iterations",
            u64::try_from(run.report.iters).expect("iteration count fits u64"),
        )
        .u64(
            "evaluations",
            u64::try_from(run.calls).expect("evaluation count fits u64"),
        );
    for value in (fixture.start)(dimension) {
        builder = builder.f64_bits("initial-coordinate", value);
    }
    builder
        .u64(
            "evaluated-point-count",
            u64::try_from(run.evaluated_points.len()).expect("evaluated point count fits u64"),
        )
        .u64(
            "evaluated-point-width",
            u64::try_from(dimension).expect("point width fits u64"),
        )
        .bytes("evaluated-coordinate-bits-le", &trace_bytes)
        .finish()
        .hex()
}

fn exercise(fixture: Fixture) {
    for dimension in DIMENSIONS {
        let probe = (fixture.gradient_probe)(dimension);
        assert_eq!(probe.len(), dimension, "gradient probe dimension drift");
        let objective = |x: &[f64]| (fixture.fg)(x).0;
        let gradient = |x: &[f64]| (fixture.fg)(x).1;
        let gate_error = check_grad(&objective, &gradient, &probe)
            .expect("FrankenScipy finite-difference gradient gate must run");
        let corrupted_gradient = |x: &[f64]| {
            let mut g = (fixture.fg)(x).1;
            g[0] += 1e-2;
            g
        };
        let mutant_error = check_grad(&objective, &corrupted_gradient, &probe)
            .expect("FrankenScipy mutant gradient gate must run");

        let first = run_once(fixture, dimension);
        let repeat = run_once(fixture, dimension);
        let replay_equal = public_replay_equal(&first, &repeat);
        let trace_identity = trace_identity(fixture, dimension, &first);
        let finite = first.state.f.is_finite()
            && first.report.grad_norm.is_finite()
            && first.state.x.iter().all(|value| value.is_finite())
            && first.state.g.iter().all(|value| value.is_finite());
        let accounting_exact = first.calls == first.state.evals
            && first.calls == first.report.evals
            && first.evaluated_points.len() == first.calls;
        let pass = gate_error <= fixture.gradient_tolerance
            && mutant_error >= MUTANT_FLOOR
            && finite
            && first.state.f <= fixture.objective_target
            && first.report.reason == StopReason::ObjectiveBelow
            && first.calls > 1
            && first.calls <= fixture.evaluation_ceiling
            && accounting_exact
            && replay_equal;
        let case = format!("{}/n={dimension}", fixture.name);
        let detail = format!(
            "n={dimension}; finite-difference gradient error={gate_error:.6e} <= \
             {:.1e}; corrupted-gradient error={mutant_error:.6e} >= \
             {MUTANT_FLOOR:.1e}; f={:.6e} <= {:.1e}; grad_inf={:.6e}; \
             reason={:?}; iterations={}; evaluations={} in [2, {}]; \
             accounting_exact={accounting_exact}; full_evaluation_trace_replay={replay_equal}; \
             trace_identity={trace_identity}",
            fixture.gradient_tolerance,
            first.state.f,
            fixture.objective_target,
            first.report.grad_norm,
            first.report.reason,
            first.report.iters,
            first.calls,
            fixture.evaluation_ceiling,
        );
        emit_verdict(&case, pass, &detail);
    }
}

fn extended_rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
    assert_eq!(x.len() % 2, 0, "extended Rosenbrock requires even n");
    let mut f = 0.0;
    let mut g = vec![0.0; x.len()];
    for (x_block, g_block) in x.chunks_exact(2).zip(g.chunks_exact_mut(2)) {
        let a = x_block[0];
        let b = x_block[1];
        let valley = b - a * a;
        let center = 1.0 - a;
        f += 100.0 * valley * valley + center * center;
        g_block[0] = -400.0 * a * valley - 2.0 * center;
        g_block[1] = 200.0 * valley;
    }
    (f, g)
}

fn rosenbrock_start(dimension: usize) -> Vec<f64> {
    assert_eq!(dimension % 2, 0);
    (0..dimension)
        .map(|index| if index % 2 == 0 { -1.2 } else { 1.0 })
        .collect()
}

fn rosenbrock_probe(dimension: usize) -> Vec<f64> {
    assert_eq!(dimension % 2, 0);
    (0..dimension)
        .map(|index| if index % 2 == 0 { 0.8 } else { 0.7 })
        .collect()
}

fn extended_powell_singular(x: &[f64]) -> (f64, Vec<f64>) {
    assert_eq!(x.len() % 4, 0, "extended Powell requires n divisible by 4");
    let mut f = 0.0;
    let mut g = vec![0.0; x.len()];
    for (x_block, g_block) in x.chunks_exact(4).zip(g.chunks_exact_mut(4)) {
        let t1 = x_block[0] + 10.0 * x_block[1];
        let t2 = x_block[2] - x_block[3];
        let t3 = x_block[1] - 2.0 * x_block[2];
        let t4 = x_block[0] - x_block[3];
        let t3_sq = t3 * t3;
        let t4_sq = t4 * t4;
        let t3_cube = t3_sq * t3;
        let t4_cube = t4_sq * t4;
        f += t1 * t1 + 5.0 * t2 * t2 + t3_sq * t3_sq + 10.0 * t4_sq * t4_sq;
        g_block[0] = 2.0 * t1 + 40.0 * t4_cube;
        g_block[1] = 20.0 * t1 + 4.0 * t3_cube;
        g_block[2] = 10.0 * t2 - 8.0 * t3_cube;
        g_block[3] = -10.0 * t2 - 40.0 * t4_cube;
    }
    (f, g)
}

fn powell_start(dimension: usize) -> Vec<f64> {
    assert_eq!(dimension % 4, 0);
    const BLOCK: [f64; 4] = [3.0, -1.0, 0.0, 1.0];
    (0..dimension).map(|index| BLOCK[index % 4]).collect()
}

fn powell_probe(dimension: usize) -> Vec<f64> {
    assert_eq!(dimension % 4, 0);
    const BLOCK: [f64; 4] = [0.3, -0.1, 0.05, 0.1];
    (0..dimension).map(|index| BLOCK[index % 4]).collect()
}

fn variably_dimensioned(x: &[f64]) -> (f64, Vec<f64>) {
    let residuals: Vec<f64> = x.iter().map(|value| value - 1.0).collect();
    let sum_sq = residuals.iter().map(|value| value * value).sum::<f64>();
    let weighted = residuals
        .iter()
        .enumerate()
        .map(|(index, value)| (index + 1) as f64 * value)
        .sum::<f64>();
    let weighted_sq = weighted * weighted;
    let common = 2.0 * weighted + 4.0 * weighted_sq * weighted;
    let gradient = residuals
        .iter()
        .enumerate()
        .map(|(index, value)| 2.0 * value + (index + 1) as f64 * common)
        .collect();
    (sum_sq + weighted_sq + weighted_sq * weighted_sq, gradient)
}

fn variably_dimensioned_start(dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|index| 1.0 - (index + 1) as f64 / dimension as f64)
        .collect()
}

fn variably_dimensioned_probe(dimension: usize) -> Vec<f64> {
    (0..dimension)
        .map(|index| {
            let sign = if index % 2 == 0 { 1.0 } else { -1.0 };
            1.0 + sign * 0.01 / (index + 1) as f64
        })
        .collect()
}

#[test]
fn extended_rosenbrock_scales_with_live_gradient_and_replay_gates() {
    exercise(Fixture {
        name: "extended-rosenbrock",
        fg: extended_rosenbrock,
        start: rosenbrock_start,
        gradient_probe: rosenbrock_probe,
        gradient_tolerance: 5e-5,
        objective_target: 1e-12,
        evaluation_ceiling: 2_400,
    });
}

#[test]
fn extended_powell_scales_with_live_gradient_and_replay_gates() {
    exercise(Fixture {
        name: "extended-powell-singular",
        fg: extended_powell_singular,
        start: powell_start,
        gradient_probe: powell_probe,
        gradient_tolerance: 5e-5,
        objective_target: 1e-10,
        evaluation_ceiling: 6_000,
    });
}

#[test]
fn variably_dimensioned_scales_with_live_gradient_and_replay_gates() {
    exercise(Fixture {
        name: "variably-dimensioned",
        fg: variably_dimensioned,
        start: variably_dimensioned_start,
        gradient_probe: variably_dimensioned_probe,
        // Forward differencing has O(h * i^2) truncation on this weighted
        // family; the n=64 aggregate bound is still 20x below the mutant.
        gradient_tolerance: 5e-4,
        objective_target: 1e-12,
        evaluation_ceiling: 3_000,
    });
}
