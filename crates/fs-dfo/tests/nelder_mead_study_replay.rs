//! G5 study-scale replay for the production Nelder–Mead path (7tv.21.29).
//!
//! The fixture reuses the fixed two-dimensional Rosenbrock leg from the
//! crate's retained golden. It records every objective callback and every
//! public tuple output bit, checks the initial simplex and the optimizer's
//! soft evaluation-budget accounting, and compares a complete canonical
//! result frame across independent runs. A disclosed `StreamKey` mutation
//! changes one returned-coordinate mantissa bit, is resealed into a valid
//! result frame, emits independently reproducible red fs-obs evidence, and is
//! refused by both the reference-identity gate and a test-local merge gate.
//!
//! This is one deterministic fixed-input study. It makes no optimizer-quality,
//! all-objective, all-dimension, all-configuration, seeded-optimizer,
//! strict-budget, cross-ISA, cancellation, checkpoint, authenticated-ledger,
//! or performance claim.

use fs_dfo::nelder_mead;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/nelder-mead-study-replay";
const CASE: &str = "fixed-rosenbrock-2d-full-trace";
const RED_CASE: &str = "seeded-returned-coordinate-corruption";

const DIMENSION: usize = 2;
const X0: [f64; DIMENSION] = [0.3, -0.2];
const SIMPLEX_SCALE: f64 = 0.2;
const MAX_EVALUATIONS: usize = 2_000;
const IMPOSSIBLE_TARGET: f64 = -1.0;

const REFLECTION: f64 = 1.0;
const EXPANSION: f64 = 2.0;
const CONTRACTION: f64 = 0.5;
const SHRINK: f64 = 0.5;

const MUTATION_SEED: u64 = 0xD0F0_7E1D_0000_0029;
const MUTATION_KERNEL: u32 = 0xD029;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION > 0);
const _: () = assert!(MAX_EVALUATIONS > DIMENSION + 1);
const _: () = assert!(IMPOSSIBLE_TARGET < 0.0);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    value: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    evaluations: Vec<EvaluationBits>,
    x_best: Vec<u64>,
    f_best: u64,
    evals: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    record: StudyRecord,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
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
    first_mismatch: String,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn rosenbrock(point: &[f64]) -> f64 {
    point
        .windows(2)
        .map(|window| {
            let a = 1.0 - window[0];
            let b = window[1] - window[0] * window[0];
            100.0f64.mul_add(b * b, a * a)
        })
        .sum()
}

fn evaluation_bits(point: &[f64], value: f64) -> EvaluationBits {
    EvaluationBits {
        point: point
            .iter()
            .map(|coordinate| coordinate.to_bits())
            .collect(),
        value: value.to_bits(),
    }
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-nelder-mead-study-fixture-v1")
        .str("algorithm", "fs_dfo::nelder_mead")
        .str("objective", "rosenbrock-sum-of-squares-v1")
        .str("coordinate-units", "dimensionless")
        .str("objective-units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .f64_bits("simplex-scale", SIMPLEX_SCALE)
        .u64("soft-max-evaluations", usize_u64(MAX_EVALUATIONS))
        .f64_bits("impossible-target", IMPOSSIBLE_TARGET)
        .f64_bits("reflection-alpha", REFLECTION)
        .f64_bits("expansion-gamma", EXPANSION)
        .f64_bits("contraction-rho", CONTRACTION)
        .f64_bits("shrink-sigma", SHRINK)
        .str("vertex-order", "f64-total-cmp-then-lowest-index")
        .str(
            "budget-semantics",
            "loop-entry-soft-cap;one-in-flight-nelder-mead-transition-may-finish",
        )
        .str("optimizer-randomness", "none-fixed-input")
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "fs-rand-stream-position-domain",
            fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
        )
        .str(
            "no-claims",
            "quality;all-objectives;all-dimensions;all-configurations;seeded-optimizer;strict-budget;cross-ISA;Cx;checkpoint;authenticated-ledger;performance",
        );
    for (coordinate, value) in X0.into_iter().enumerate() {
        builder = builder
            .u64("initial-coordinate-index", usize_u64(coordinate))
            .f64_bits("initial-coordinate", value);
    }
    builder.finish()
}

fn result_identity(fixture: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-nelder-mead-study-result-v1")
        .child("fixture", fixture)
        .u64("reported-evals", usize_u64(record.evals))
        .u64(
            "evaluation-trace-length",
            usize_u64(record.evaluations.len()),
        )
        .u64("returned-point-length", usize_u64(record.x_best.len()))
        .f64_bits("returned-objective", f64::from_bits(record.f_best));
    for (evaluation_index, evaluation) in record.evaluations.iter().enumerate() {
        builder = builder
            .u64("evaluation-index", usize_u64(evaluation_index))
            .u64("evaluation-point-length", usize_u64(evaluation.point.len()));
        for (coordinate, &value) in evaluation.point.iter().enumerate() {
            builder = builder
                .u64("evaluation-coordinate-index", usize_u64(coordinate))
                .f64_bits("evaluation-coordinate", f64::from_bits(value));
        }
        builder = builder.f64_bits("evaluation-objective", f64::from_bits(evaluation.value));
    }
    for (coordinate, &value) in record.x_best.iter().enumerate() {
        builder = builder
            .u64("returned-coordinate-index", usize_u64(coordinate))
            .f64_bits("returned-coordinate", f64::from_bits(value));
    }
    builder.finish()
}

fn run_study() -> StudyRun {
    let mut evaluations = Vec::with_capacity(MAX_EVALUATIONS + DIMENSION + 1);
    let (x_best, f_best, evals) = {
        let mut objective = |point: &[f64]| {
            let value = rosenbrock(point);
            evaluations.push(evaluation_bits(point, value));
            value
        };
        nelder_mead(
            &mut objective,
            &X0,
            SIMPLEX_SCALE,
            MAX_EVALUATIONS,
            IMPOSSIBLE_TARGET,
        )
    };
    let record = StudyRecord {
        evaluations,
        x_best: x_best
            .iter()
            .map(|coordinate| coordinate.to_bits())
            .collect(),
        f_best: f_best.to_bits(),
        evals,
    };
    let fixture = fixture_identity();
    let result = result_identity(&fixture, &record);
    StudyRun {
        fixture,
        record,
        result,
    }
}

fn expected_initial_simplex() -> Vec<Vec<u64>> {
    let mut simplex = vec![X0.to_vec()];
    for coordinate in 0..DIMENSION {
        let mut vertex = X0.to_vec();
        vertex[coordinate] += SIMPLEX_SCALE;
        simplex.push(vertex);
    }
    simplex
        .into_iter()
        .map(|vertex| vertex.into_iter().map(f64::to_bits).collect())
        .collect()
}

#[allow(clippy::too_many_lines)] // Full callback/output accounting is the receipt.
fn accounting_mismatch(record: &StudyRecord) -> Option<String> {
    if record.evals != record.evaluations.len() {
        return Some(format!(
            "reported-evals:{}!=closure-count:{}",
            record.evals,
            record.evaluations.len()
        ));
    }
    let maximum_soft_overshoot = MAX_EVALUATIONS + DIMENSION + 1;
    if !(MAX_EVALUATIONS..=maximum_soft_overshoot).contains(&record.evals) {
        return Some(format!(
            "soft-budget-accounting:{} not in {MAX_EVALUATIONS}..={maximum_soft_overshoot}",
            record.evals
        ));
    }
    if record.x_best.len() != DIMENSION {
        return Some(format!(
            "returned-point-dimension:{}!=expected-{DIMENSION}",
            record.x_best.len()
        ));
    }

    let initial_simplex = expected_initial_simplex();
    if record.evaluations.len() < initial_simplex.len() {
        return Some(format!(
            "trace-too-short-for-simplex:{}<{}",
            record.evaluations.len(),
            initial_simplex.len()
        ));
    }
    for (index, expected) in initial_simplex.iter().enumerate() {
        if record.evaluations[index].point != *expected {
            return Some(format!(
                "initial-simplex[{index}]:actual={:016x?};expected={expected:016x?}",
                record.evaluations[index].point
            ));
        }
    }

    let mut minimum = f64::INFINITY;
    for (index, evaluation) in record.evaluations.iter().enumerate() {
        if evaluation.point.len() != DIMENSION {
            return Some(format!(
                "evaluation[{index}]-dimension:{}!=expected-{DIMENSION}",
                evaluation.point.len()
            ));
        }
        let point: Vec<f64> = evaluation
            .point
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        if point.iter().any(|coordinate| !coordinate.is_finite()) {
            return Some(format!(
                "evaluation[{index}]-non-finite-point:{:016x?}",
                evaluation.point
            ));
        }
        let value = f64::from_bits(evaluation.value);
        if !value.is_finite() || value <= IMPOSSIBLE_TARGET {
            return Some(format!(
                "evaluation[{index}]-invalid-objective:0x{:016x}",
                evaluation.value
            ));
        }
        let recomputed = rosenbrock(&point).to_bits();
        if recomputed != evaluation.value {
            return Some(format!(
                "evaluation[{index}]-objective:recomputed=0x{recomputed:016x};recorded=0x{:016x}",
                evaluation.value
            ));
        }
        minimum = minimum.min(value);
    }

    let returned_point: Vec<f64> = record.x_best.iter().copied().map(f64::from_bits).collect();
    if returned_point
        .iter()
        .any(|coordinate| !coordinate.is_finite())
    {
        return Some(format!("returned-point-non-finite:{:016x?}", record.x_best));
    }
    let returned_recomputed = rosenbrock(&returned_point).to_bits();
    if returned_recomputed != record.f_best {
        return Some(format!(
            "returned-objective:recomputed=0x{returned_recomputed:016x};reported=0x{:016x}",
            record.f_best
        ));
    }
    if minimum.to_bits() != record.f_best {
        return Some(format!(
            "global-trace-minimum=0x{:016x};returned=0x{:016x}",
            minimum.to_bits(),
            record.f_best
        ));
    }
    if !record
        .evaluations
        .iter()
        .any(|evaluation| evaluation.point == record.x_best && evaluation.value == record.f_best)
    {
        return Some("returned-best-is-not-an-evaluated-point".to_string());
    }
    None
}

fn first_record_mismatch(left: &StudyRecord, right: &StudyRecord) -> Option<String> {
    if left.evaluations.len() != right.evaluations.len() {
        return Some(format!(
            "evaluations.length:{}!={}",
            left.evaluations.len(),
            right.evaluations.len()
        ));
    }
    for (index, (a, b)) in left.evaluations.iter().zip(&right.evaluations).enumerate() {
        if a != b {
            return Some(format!("evaluations[{index}]:left={a:?};right={b:?}"));
        }
    }
    if left.x_best.len() != right.x_best.len() {
        return Some(format!(
            "x_best.length:{}!={}",
            left.x_best.len(),
            right.x_best.len()
        ));
    }
    if let Some((coordinate, (a, b))) = left
        .x_best
        .iter()
        .zip(&right.x_best)
        .enumerate()
        .find(|(_, (a, b))| a != b)
    {
        return Some(format!("x_best[{coordinate}]:0x{a:016x}!=0x{b:016x}"));
    }
    if left.f_best != right.f_best {
        return Some(format!(
            "f_best:0x{:016x}!=0x{:016x}",
            left.f_best, right.f_best
        ));
    }
    (left.evals != right.evals).then(|| format!("evals:{}!={}", left.evals, right.evals))
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed = result_identity(&run.fixture, &run.record);
    if computed.canonical_bytes() == run.result.canonical_bytes() {
        Ok(())
    } else {
        Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.result.root(),
            computed: computed.root(),
        })
    }
}

fn admit_against(run: &StudyRun, reference: &ReplayIdentity) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.canonical_bytes() {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: reference.root(),
            found: run.result.root(),
        })
    }
}

fn exact_returned_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(&reference_bits) = reference.record.x_best.get(mutation.coordinate) else {
        return false;
    };
    let Some(&mutant_bits) = mutant.record.x_best.get(mutation.coordinate) else {
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
    expected.x_best[mutation.coordinate] = mutation.after;
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
        .expect("coordinate index fits usize");
    let mantissa_bit = u32::try_from(selector.next_below(20)).expect("mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.x_best[coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.x_best[coordinate] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed result mutation must refuse");
    run.result = result_identity(&run.fixture, &run.record);
    let reference_error = admit_against(&run, &reference.result)
        .expect_err("resealed mutation must not match the retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &run.record)
        .expect("seeded mutation changes the returned record");
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
        first_mismatch,
    }
}

fn corruption_detail(reference: &StudyRun, corruption: &SeededCorruption) -> String {
    format!(
        "fixture={}; reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=x_best[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}",
        reference.fixture.hex(),
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
        corruption.first_mismatch,
    )
}

fn failure_event(detail: &str, mutation: Mutation) -> Event {
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail: detail.to_string(),
            seed: mutation.seed,
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

fn emit_green_receipt(run: &StudyRun) {
    let overshoot = run.record.evals - MAX_EVALUATIONS;
    let mut emitter = Emitter::new(SUITE, CASE);
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "nelder-mead-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::nelder_mead\",",
                    "\"objective\":\"rosenbrock-2d\",\"units\":\"dimensionless\",",
                    "\"soft_max_evaluations\":{},\"actual_evaluations\":{},",
                    "\"soft_budget_overshoot\":{},\"trace_length\":{},",
                    "\"returned_point_length\":{},\"returned_objective_bits\":\"0x{:016x}\",",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"optimizer-quality\",\"all-objectives\",",
                    "\"all-dimensions\",\"all-configurations\",\"seeded-optimizer\",",
                    "\"strict-budget\",\"cross-ISA\",\"cancellation\",",
                    "\"checkpointing\",\"authenticated-ledger\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                MAX_EVALUATIONS,
                run.record.evals,
                overshoot,
                run.record.evaluations.len(),
                run.record.x_best.len(),
                run.record.f_best,
                fs_dfo::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
            ),
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("Nelder-Mead receipt must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("Nelder-Mead receipt must use the fs-obs wire schema");
    let receipt = event.content_identity_receipt();
    event
        .admit_content_identity(&receipt)
        .expect("fresh receipt identity must admit exactly");
    println!("{line}");
}

fn emit_green_verdict(run: &StudyRun) -> Event {
    let detail = format!(
        "fixture={}; result={}; evaluations={}; overshoot={}; complete_trace=bit-exact; public_tuple=fully-bound",
        run.fixture.hex(),
        run.result.hex(),
        run.record.evals,
        run.record.evals - MAX_EVALUATIONS,
    );
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail,
            seed: 0,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("Nelder-Mead verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("Nelder-Mead verdict must use the fs-obs wire schema");
    println!("{line}");
    event
}

fn exercise_seeded_corruption(original: &StudyRun, replay: &StudyRun) {
    let first = seeded_corruption(original);
    let second = seeded_corruption(replay);
    assert_eq!(first, second, "seeded red state must replay exactly");
    assert!(
        exact_returned_bit_delta(original, &first.run, first.mutation),
        "the corruption must change exactly one returned bit"
    );
    assert!(
        exact_returned_bit_delta(replay, &second.run, second.mutation),
        "the replay corruption must change exactly one returned bit"
    );
    assert!(first.mutation.coordinate < DIMENSION);
    assert!((0..20).contains(&first.mutation.mantissa_bit));
    assert!(f64::from_bits(first.mutation.after).is_finite());
    assert!(matches!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == original.result.root() && computed == first.run.result.root()
    ));
    assert!(matches!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result.root() && found == first.run.result.root()
    ));
    assert!(
        first
            .first_mismatch
            .starts_with(&format!("x_best[{}]", first.mutation.coordinate))
    );

    let first_detail = corruption_detail(original, &first);
    let second_detail = corruption_detail(replay, &second);
    assert_eq!(first_detail, second_detail);
    let first_event = failure_event(&first_detail, first.mutation);
    let second_event = failure_event(&second_detail, second.mutation);
    for event in [&first_event, &second_event] {
        fs_obs::lint_failure_record(event)
            .expect("seeded Nelder-Mead corruption must retain replay inputs");
        fs_obs::validate_line(&event.to_jsonl())
            .expect("seeded Nelder-Mead corruption must remain wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("red event identity must admit its exact content");
    }
    assert_eq!(first_event, second_event);
    assert_eq!(
        first_event.content_identity().canonical_bytes(),
        second_event.content_identity().canonical_bytes()
    );
    println!("{}", first_event.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_event))
        .expect_err("the merge gate must refuse seeded returned-bit corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!("x_best[{}]", first.mutation.coordinate)));
    assert!(message.contains("ReferenceIdentityMismatch"));
}

#[test]
fn nelder_mead_full_study_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();
    let original_accounting = accounting_mismatch(&original.record);
    let replay_accounting = accounting_mismatch(&replay.record);
    assert_eq!(original_accounting, None, "original accounting failed");
    assert_eq!(replay_accounting, None, "replay accounting failed");
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));

    let mismatch = first_record_mismatch(&original.record, &replay.record);
    assert_eq!(mismatch, None, "full study replay drifted");
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.result, replay.result);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete result frames must replay byte-for-byte"
    );

    emit_green_receipt(&original);
    let green = emit_green_verdict(&original);
    assert_mergeable(&green);
    exercise_seeded_corruption(&original, &replay);
}
