//! G0/G3/G5 full-study replay for the standalone production L-BFGS path
//! (7tv.21.47).
//!
//! The fixture retains every objective/gradient callback, every public final
//! `LbfgsState` field, and every `LbfgsReport` field as exact IEEE-754 bits.
//! It proves one-shot, repeated, and checkpoint/resume executions produce the
//! same complete receipt. A seeded one-bit final-coordinate corruption must
//! fail while stale, after resealing against the retained reference, and under
//! an independent semantic oracle.
//!
//! This is one deterministic, same-ISA, fixed-input Rosenbrock study. It makes
//! no claim about arbitrary objectives, dimensions, memory lengths, stop rules,
//! cross-ISA equality, cancellation, private curvature-pair state,
//! authenticated ledger trust, or performance.

use fs_ascent::{LbfgsReport, LbfgsState, StopReason, StopRule};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-ascent/lbfgs-study-replay";
const CASE: &str = "fixed-rosenbrock-6d-full-trace";
const RED_CASE: &str = "seeded-final-coordinate-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-ascent-lbfgs-study-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-ascent-lbfgs-study-result-v1";

const DIMENSION: usize = 6;
const START: [f64; DIMENSION] = [-1.2; DIMENSION];
const MEMORY: usize = 10;
const GRADIENT_TOLERANCE: f64 = 1.0e-8;
const EVALUATION_BUDGET: usize = 5_000;
const MAX_ITERATIONS: usize = 2_000;
const RESUME_CUT: usize = 5;
const ORACLE_TOLERANCE_FACTOR: f64 = 4_096.0;

const MUTATION_SEED: u64 = 0x1B_F6_55_7D_21_47;
const MUTATION_KERNEL: u32 = 0xB647;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION >= 2);
const _: () = assert!(MEMORY > 0);
const _: () = assert!(RESUME_CUT < MAX_ITERATIONS);
const _: () = assert!(GRADIENT_TOLERANCE > 0.0);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    objective: u64,
    gradient: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StateBits {
    point: Vec<u64>,
    objective: u64,
    gradient: Vec<u64>,
    memory: usize,
    iterations: usize,
    evaluations: usize,
    history: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ReportBits {
    reason: StopReason,
    gradient_norm: u64,
    objective: u64,
    iterations: usize,
    evaluations: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    evaluations: Vec<EvaluationBits>,
    state: StateBits,
    report: ReportBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    record: StudyRecord,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
    SemanticInconsistency(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    coordinate: usize,
    mantissa_bit: u32,
    selector_draws: u64,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_error: AdmissionError,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixed fixture cardinality fits u64")
}

fn stop_rule() -> StopRule {
    StopRule::Any(vec![
        StopRule::GradNorm(GRADIENT_TOLERANCE),
        StopRule::Budget(EVALUATION_BUDGET),
    ])
}

fn stop_reason_name(reason: &StopReason) -> &'static str {
    match reason {
        StopReason::GradNorm => "gradient-norm",
        StopReason::ObjectiveBelow => "objective-below",
        StopReason::Budget => "evaluation-budget",
        StopReason::Stall => "stall",
        StopReason::Composite => "composite",
        StopReason::IterationCap => "iteration-cap",
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn values(bits: &[u64]) -> Vec<f64> {
    bits.iter().copied().map(f64::from_bits).collect()
}

fn gradient_inf_norm(gradient: &[f64]) -> f64 {
    gradient
        .iter()
        .map(|component| component.abs())
        .fold(0.0f64, f64::max)
}

/// Objective/gradient implementation passed to production L-BFGS. Its term
/// order matches the existing ASCENT Rosenbrock conformance fixture.
fn rosenbrock_under_test(point: &[f64]) -> (f64, Vec<f64>) {
    assert_eq!(point.len(), DIMENSION);
    let mut objective = 0.0f64;
    let mut gradient = vec![0.0f64; point.len()];
    for index in 0..point.len() - 1 {
        let a = 1.0 - point[index];
        let b = point[index + 1] - point[index] * point[index];
        objective += a.mul_add(a, 100.0 * b * b);
        gradient[index] += (-2.0f64).mul_add(a, -400.0 * point[index] * b);
        gradient[index + 1] += 200.0 * b;
    }
    (objective, gradient)
}

/// Algebraically expanded oracle used only to admit retained callbacks. It
/// deliberately does not share the producer's residual-square or fused form.
fn expanded_rosenbrock_oracle(point: &[f64]) -> Option<(f64, Vec<f64>)> {
    if point.len() != DIMENSION {
        return None;
    }
    let mut objective = 0.0f64;
    let mut gradient = vec![0.0f64; point.len()];
    for index in 0..point.len() - 1 {
        let x = point[index];
        let next = point[index + 1];
        let x2 = x * x;
        let x4 = x2 * x2;
        objective += (1.0 - 2.0 * x + x2) + (100.0 * next * next - 200.0 * next * x2 + 100.0 * x4);
        gradient[index] += (2.0 * x - 2.0) + (400.0 * x * x2 - 400.0 * x * next);
        gradient[index + 1] += 200.0 * next - 200.0 * x2;
    }
    Some((objective, gradient))
}

fn within_oracle_tolerance(produced: f64, reference: f64) -> bool {
    if !produced.is_finite() || !reference.is_finite() {
        return false;
    }
    let scale = produced.abs().max(reference.abs()).max(1.0);
    (produced - reference).abs() <= ORACLE_TOLERANCE_FACTOR * f64::EPSILON * scale
}

fn evaluation_bits(point: &[f64], objective: f64, gradient: &[f64]) -> EvaluationBits {
    EvaluationBits {
        point: bits(point),
        objective: objective.to_bits(),
        gradient: bits(gradient),
    }
}

fn state_bits(state: &LbfgsState) -> StateBits {
    StateBits {
        point: bits(&state.x),
        objective: state.f.to_bits(),
        gradient: bits(&state.g),
        memory: state.memory,
        iterations: state.iters,
        evaluations: state.evals,
        history: bits(&state.history),
    }
}

fn report_bits(report: &LbfgsReport) -> ReportBits {
    ReportBits {
        reason: report.reason.clone(),
        gradient_norm: report.grad_norm.to_bits(),
        objective: report.f.to_bits(),
        iterations: report.iters,
        evaluations: report.evals,
    }
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_ascent::LbfgsState")
        .str("objective", "rosenbrock-chain-6d-v1")
        .str(
            "producer-formula",
            "sum[(1-x_i)^2+100*(x_(i+1)-x_i^2)^2];fused-a-square",
        )
        .str(
            "oracle-formula",
            "expanded-polynomial-objective-and-gradient-v1",
        )
        .str(
            "oracle-bound",
            "4096*EPSILON*max(1,abs(produced),abs(reference))",
        )
        .f64_bits("oracle-tolerance-factor", ORACLE_TOLERANCE_FACTOR)
        .str("coordinate-units", "dimensionless")
        .str("objective-units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .u64("memory", usize_u64(MEMORY))
        .f64_bits("gradient-tolerance", GRADIENT_TOLERANCE)
        .u64("evaluation-budget", usize_u64(EVALUATION_BUDGET))
        .u64("maximum-iterations", usize_u64(MAX_ITERATIONS))
        .u64("resume-cut", usize_u64(RESUME_CUT))
        .str(
            "stop-rule",
            "Any([GradNorm(1e-8),Budget(5000)]);first-satisfied-child-attribution",
        )
        .str("line-search", "deterministic-strong-wolfe-c1=1e-4-c2=0.9")
        .str("checkpoint", "LbfgsState::clone-after-five-accepted-iterations")
        .str("optimizer-randomness", "none-fixed-input-no-algorithm-seed")
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-ascent-version", fs_ascent::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .str(
            "no-claims",
            "arbitrary-objectives;arbitrary-dimensions;arbitrary-memory;arbitrary-stop-rules;cross-ISA;Cx;cancellation;private-curvature-pairs;authenticated-ledger;performance",
        );
    for (coordinate, start) in START.into_iter().enumerate() {
        builder = builder
            .u64("start-coordinate-index", usize_u64(coordinate))
            .f64_bits("start-coordinate", start);
    }
    builder.finish()
}

#[allow(clippy::too_many_lines)] // Every retained callback and public field is identity-bearing.
fn result_identity(fixture: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .child("fixture-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .u64("callback-count", usize_u64(record.evaluations.len()));
    for (evaluation_index, evaluation) in record.evaluations.iter().enumerate() {
        builder = builder
            .u64("callback-index", usize_u64(evaluation_index))
            .u64("callback-point-length", usize_u64(evaluation.point.len()));
        for (coordinate, &coordinate_bits) in evaluation.point.iter().enumerate() {
            builder = builder
                .u64("callback-point-index", usize_u64(coordinate))
                .f64_bits("callback-point", f64::from_bits(coordinate_bits));
        }
        builder = builder
            .f64_bits("callback-objective", f64::from_bits(evaluation.objective))
            .u64(
                "callback-gradient-length",
                usize_u64(evaluation.gradient.len()),
            );
        for (coordinate, &gradient_bits) in evaluation.gradient.iter().enumerate() {
            builder = builder
                .u64("callback-gradient-index", usize_u64(coordinate))
                .f64_bits("callback-gradient", f64::from_bits(gradient_bits));
        }
    }

    let state = &record.state;
    builder = builder
        .u64("state-point-length", usize_u64(state.point.len()))
        .f64_bits("state-objective", f64::from_bits(state.objective))
        .u64("state-gradient-length", usize_u64(state.gradient.len()))
        .u64("state-memory", usize_u64(state.memory))
        .u64("state-iterations", usize_u64(state.iterations))
        .u64("state-evaluations", usize_u64(state.evaluations))
        .u64("state-history-length", usize_u64(state.history.len()));
    for (coordinate, &coordinate_bits) in state.point.iter().enumerate() {
        builder = builder
            .u64("state-point-index", usize_u64(coordinate))
            .f64_bits("state-point", f64::from_bits(coordinate_bits));
    }
    for (coordinate, &gradient_bits) in state.gradient.iter().enumerate() {
        builder = builder
            .u64("state-gradient-index", usize_u64(coordinate))
            .f64_bits("state-gradient", f64::from_bits(gradient_bits));
    }
    for (history_index, &objective_bits) in state.history.iter().enumerate() {
        builder = builder
            .u64("state-history-index", usize_u64(history_index))
            .f64_bits("state-history-objective", f64::from_bits(objective_bits));
    }

    let report = &record.report;
    builder
        .str("report-reason", stop_reason_name(&report.reason))
        .f64_bits("report-gradient-norm", f64::from_bits(report.gradient_norm))
        .f64_bits("report-objective", f64::from_bits(report.objective))
        .u64("report-iterations", usize_u64(report.iterations))
        .u64("report-evaluations", usize_u64(report.evaluations))
        .finish()
}

fn run_study(resume_at: Option<usize>) -> StudyRun {
    let mut evaluations = Vec::new();
    let mut callback = |point: &[f64]| {
        let (objective, gradient) = rosenbrock_under_test(point);
        evaluations.push(evaluation_bits(point, objective, &gradient));
        (objective, gradient)
    };
    let rule = stop_rule();
    let mut state = LbfgsState::new(&START, MEMORY, &mut callback);
    let report = match resume_at {
        Some(cut) => {
            let prefix = state.run(&mut callback, &rule, cut);
            assert_eq!(
                prefix.reason,
                StopReason::IterationCap,
                "resume cut must occur before the fixed fixture converges"
            );
            assert_eq!(state.iters, cut, "resume cut counts accepted iterations");
            let mut resumed = state.clone();
            let report = resumed.run(&mut callback, &rule, MAX_ITERATIONS - cut);
            state = resumed;
            report
        }
        None => state.run(&mut callback, &rule, MAX_ITERATIONS),
    };
    drop(callback);

    let record = StudyRecord {
        evaluations,
        state: state_bits(&state),
        report: report_bits(&report),
    };
    let fixture = fixture_identity();
    let result = result_identity(&fixture, &record);
    StudyRun {
        fixture,
        record,
        result,
    }
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let expected_fixture = fixture_identity();
    if run.fixture != expected_fixture {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.fixture.root(),
            computed: expected_fixture.root(),
        });
    }
    let computed_result = result_identity(&run.fixture, &run.record);
    if run.result != computed_result {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.result.root(),
            computed: computed_result.root(),
        });
    }
    Ok(())
}

fn history_is_callback_subsequence(history: &[u64], evaluations: &[EvaluationBits]) -> bool {
    let mut cursor = 0usize;
    for &objective in history {
        let Some(offset) = evaluations[cursor..]
            .iter()
            .position(|evaluation| evaluation.objective == objective)
        else {
            return false;
        };
        cursor += offset + 1;
    }
    true
}

#[allow(clippy::too_many_lines)] // The semantic gate checks each retained callback/public field.
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    let state = &record.state;
    let report = &record.report;
    if record.evaluations.is_empty() {
        return Some("empty-callback-trace".to_string());
    }
    if record.evaluations.len() != state.evaluations
        || record.evaluations.len() != report.evaluations
    {
        return Some(format!(
            "callback-accounting:trace={};state={};report={}",
            record.evaluations.len(),
            state.evaluations,
            report.evaluations
        ));
    }
    if state.memory != MEMORY {
        return Some(format!("state-memory:{}!=fixture-{MEMORY}", state.memory));
    }
    if state.iterations != report.iterations || state.evaluations != report.evaluations {
        return Some(format!(
            "state-report-accounting:iters={}/{};evals={}/{}",
            state.iterations, report.iterations, state.evaluations, report.evaluations
        ));
    }
    if state.objective != report.objective {
        return Some(format!(
            "state-report-objective:0x{:016x}!=0x{:016x}",
            state.objective, report.objective
        ));
    }
    if state.point.len() != DIMENSION || state.gradient.len() != DIMENSION {
        return Some(format!(
            "final-state-shape:point={};gradient={};expected={DIMENSION}",
            state.point.len(),
            state.gradient.len()
        ));
    }
    if state.history.len() != state.iterations + 1 {
        return Some(format!(
            "history-accounting:length={}!=iterations-plus-one={}",
            state.history.len(),
            state.iterations + 1
        ));
    }
    if state.iterations <= RESUME_CUT {
        return Some(format!(
            "fixture-converged-before-meaningful-resume:iterations={}",
            state.iterations
        ));
    }
    if report.reason != StopReason::GradNorm {
        return Some(format!(
            "unexpected-stop-reason:{}",
            stop_reason_name(&report.reason)
        ));
    }
    let final_gradient = values(&state.gradient);
    let recomputed_norm = gradient_inf_norm(&final_gradient);
    if recomputed_norm.to_bits() != report.gradient_norm {
        return Some(format!(
            "gradient-certificate:recomputed=0x{:016x};reported=0x{:016x}",
            recomputed_norm.to_bits(),
            report.gradient_norm
        ));
    }
    if recomputed_norm > GRADIENT_TOLERANCE {
        return Some(format!(
            "gradient-certificate-exceeds-threshold:0x{:016x}",
            recomputed_norm.to_bits()
        ));
    }
    if state.evaluations >= EVALUATION_BUDGET {
        return Some(format!(
            "fixed-fixture-reached-budget:{}",
            state.evaluations
        ));
    }

    let expected_start = bits(&START);
    if record.evaluations[0].point != expected_start {
        return Some("first-callback-is-not-the-bound-start-point".to_string());
    }
    if state.history.first().copied() != Some(record.evaluations[0].objective) {
        return Some("history-does-not-start-at-initial-callback".to_string());
    }
    if !history_is_callback_subsequence(&state.history, &record.evaluations) {
        return Some("accepted-history-is-not-an-ordered-callback-subsequence".to_string());
    }

    for (evaluation_index, evaluation) in record.evaluations.iter().enumerate() {
        if evaluation.point.len() != DIMENSION || evaluation.gradient.len() != DIMENSION {
            return Some(format!(
                "callback[{evaluation_index}]-shape:point={};gradient={}",
                evaluation.point.len(),
                evaluation.gradient.len()
            ));
        }
        let point = values(&evaluation.point);
        let recorded_objective = f64::from_bits(evaluation.objective);
        let recorded_gradient = values(&evaluation.gradient);
        if point.iter().any(|value| !value.is_finite())
            || !recorded_objective.is_finite()
            || recorded_gradient.iter().any(|value| !value.is_finite())
        {
            return Some(format!("callback[{evaluation_index}]-non-finite"));
        }

        let (producer_objective, producer_gradient) = rosenbrock_under_test(&point);
        if producer_objective.to_bits() != evaluation.objective
            || bits(&producer_gradient) != evaluation.gradient
        {
            return Some(format!(
                "callback[{evaluation_index}]-producer-replay-mismatch"
            ));
        }
        let Some((oracle_objective, oracle_gradient)) = expanded_rosenbrock_oracle(&point) else {
            return Some(format!("callback[{evaluation_index}]-oracle-shape-refusal"));
        };
        if !within_oracle_tolerance(recorded_objective, oracle_objective) {
            return Some(format!(
                "callback[{evaluation_index}]-objective-oracle-mismatch"
            ));
        }
        for coordinate in 0..DIMENSION {
            if !within_oracle_tolerance(recorded_gradient[coordinate], oracle_gradient[coordinate])
            {
                return Some(format!(
                    "callback[{evaluation_index}]-gradient-oracle-mismatch[{coordinate}]"
                ));
            }
        }
    }

    let last = record
        .evaluations
        .last()
        .expect("nonempty callback trace checked above");
    if state.point != last.point {
        return Some("final-state-point!=last-callback-point".to_string());
    }
    if state.objective != last.objective {
        return Some("final-state-objective!=last-callback-objective".to_string());
    }
    if state.gradient != last.gradient {
        return Some("final-state-gradient!=last-callback-gradient".to_string());
    }
    if state.history.last().copied() != Some(state.objective) {
        return Some("history-does-not-end-at-final-state-objective".to_string());
    }

    let final_point = values(&state.point);
    if final_point
        .iter()
        .any(|coordinate| (*coordinate - 1.0).abs() > 1.0e-5)
    {
        return Some("fixed-rosenbrock-solution-outside-coordinate-oracle".to_string());
    }
    if f64::from_bits(state.objective) > 1.0e-12 {
        return Some("fixed-rosenbrock-objective-outside-oracle".to_string());
    }
    None
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    match semantic_mismatch(&run.record) {
        Some(mismatch) => Err(AdmissionError::SemanticInconsistency(mismatch)),
        None => Ok(()),
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result == reference.result {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: reference.result.root(),
            found: run.result.root(),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, &run.record);
}

fn exact_state_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(&reference_bits) = reference.record.state.point.get(mutation.coordinate) else {
        return false;
    };
    let Some(&mutant_bits) = mutant.record.state.point.get(mutation.coordinate) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.state.point[mutation.coordinate] = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let coordinate = usize::try_from(selector.next_below(usize_u64(DIMENSION)))
        .expect("selected coordinate fits usize");
    let mantissa_bit = u32::try_from(selector.next_below(20)).expect("selected bit fits u32");
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.state.point[coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.state.point[coordinate] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed mutation must not match retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed state mutation must remain semantically invalid");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            coordinate,
            mantissa_bit,
            selector_draws,
            before,
            after,
        },
        stale_error,
        reference_error,
        semantic_error,
    }
}

fn green_receipt(run: &StudyRun) -> Event {
    let state = &run.record.state;
    let report = &run.record.report;
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "lbfgs-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"algorithm\":\"fs_ascent::LbfgsState\",",
                    "\"dimension\":{},\"memory\":{},\"resume_cut\":{},",
                    "\"callbacks\":{},\"iterations\":{},",
                    "\"stop_reason\":\"{}\",\"gradient_norm_bits\":\"0x{:016x}\",",
                    "\"objective_bits\":\"0x{:016x}\",",
                    "\"versions\":{{\"fs_ascent\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"arbitrary-objectives\",\"arbitrary-dimensions\",",
                    "\"cross-ISA\",\"cancellation\",\"private-curvature-pairs\",",
                    "\"authenticated-ledger\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                DIMENSION,
                state.memory,
                RESUME_CUT,
                run.record.evaluations.len(),
                state.iterations,
                stop_reason_name(&report.reason),
                report.gradient_norm,
                report.objective,
                fs_ascent::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
            ),
        },
        None,
    )
}

fn green_verdict(run: &StudyRun) -> Event {
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: format!(
                "fixture={}; result={}; callbacks={}; public-state-and-report=fully-bound; one-shot=repeat=resume-cut-{RESUME_CUT}",
                run.fixture.hex(),
                run.result.hex(),
                run.record.evaluations.len(),
            ),
            seed: 0,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=state.x[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result.hex(),
        corruption.run.result.hex(),
        corruption.mutation.seed,
        corruption.mutation.kernel,
        corruption.mutation.tile,
        corruption.mutation.selector_draws,
        corruption.mutation.coordinate,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
        },
        None,
    )
}

fn assert_mergeable(event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn assert_event_replays(first: &Event, second: &Event, label: &str) {
    assert_eq!(
        first.content_identity().canonical_bytes(),
        second.content_identity().canonical_bytes(),
        "{label} content identity must replay byte-for-byte"
    );
    for event in [first, second] {
        fs_obs::lint_failure_record(event).expect("study evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("study evidence is wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("fresh evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay and all refusal gates.
fn standalone_lbfgs_full_study_replays_and_seeded_failure_is_refused() {
    let original = run_study(None);
    let replay = run_study(None);
    let resumed = run_study(Some(RESUME_CUT));
    for run in [&original, &replay, &resumed] {
        assert_eq!(validate_payload(run), Ok(()));
        assert_eq!(validate_semantics(run), Ok(()));
        assert_eq!(admit_reference(run, &original), Ok(()));
    }
    assert_eq!(original.record, replay.record);
    assert_eq!(original.record, resumed.record);
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.fixture, resumed.fixture);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.result, resumed.result);
    assert_eq!(
        original.result.canonical_bytes(),
        resumed.result.canonical_bytes(),
        "checkpoint/resume must preserve the complete result frame"
    );

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    let resumed_receipt = green_receipt(&resumed);
    assert_event_replays(&first_receipt, &second_receipt, "green receipt");
    assert_event_replays(&first_receipt, &resumed_receipt, "resumed receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    let resumed_green = green_verdict(&resumed);
    assert_event_replays(&first_green, &second_green, "green verdict");
    assert_event_replays(&first_green, &resumed_green, "resumed verdict");
    for event in [&first_green, &second_green, &resumed_green] {
        assert_mergeable(event);
    }
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded corruption must replay exactly");
    assert!(
        exact_state_bit_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one retained final-coordinate bit"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed mutation must be internally self-consistent"
    );
    assert!(f64::from_bits(first.mutation.after).is_finite());
    assert!(matches!(
        &first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if *declared == original.result.root()
                && *computed == first.run.result.root()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if *expected == original.result.root()
                && *found == first.run.result.root()
    ));
    assert!(matches!(
        &first.semantic_error,
        AdmissionError::SemanticInconsistency(mismatch)
            if mismatch == "final-state-point!=last-callback-point"
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_replays(&first_red, &second_red, "red evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded L-BFGS corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!("state.x[{}]", first.mutation.coordinate)));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
