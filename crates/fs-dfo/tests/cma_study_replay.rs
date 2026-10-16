//! G0/G3/G5 full-study replay for the standalone full-covariance CMA-ES
//! surface (7tv.21.41).
//!
//! The fixture retains every callback and every public `CmaReport` field as
//! exact IEEE-754 bits, binds the complete frame under domain-separated
//! BLAKE3, and independently checks callback accounting and the earliest
//! `total_cmp` minimum. A seeded one-bit returned-coordinate corruption is
//! refused both before and after resealing.
//!
//! This is one deterministic, same-ISA, fixed-input study. It makes no claim
//! about optimizer quality, arbitrary objectives/dimensions/parameters/seeds,
//! validation or fallibility of the legacy `cmaes` API, cross-ISA equality,
//! cancellation, checkpoint/resume, authenticated ledger trust, or
//! performance. It also does not retain or certify the optimizer's internal
//! covariance/eigensystem/evolution paths, RNG counter, candidate ordering,
//! or unexposed stop reason.

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::{CmaParams, CmaReport, cmaes};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/cma-study-replay";
const CASE: &str = "fixed-quadratic-4d-full-trace";
const RED_CASE: &str = "seeded-returned-coordinate-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-dfo-cma-study-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-dfo-cma-study-result-v1";
const FIXTURE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.cma-study-fixture.v1";
const RESULT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.cma-study-result.v1";
const EVENT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.cma-study-event.v1";

const DIMENSION: usize = 4;
const START: [f64; DIMENSION] = [1.5, -1.25, 2.0, 0.5];
const CENTER: [f64; DIMENSION] = [0.25, -0.5, 0.75, -1.0];
const WEIGHTS: [f64; DIMENSION] = [1.0, 2.0, 4.0, 8.0];
const OBJECTIVE_OFFSET: f64 = 0.125;

const INPUT_SEED: u64 = 0xC0A5_7A11_4D00_0041;
const LAMBDA: usize = 8;
const EXPECTED_GENERATIONS: usize = 8;
const MAX_EVALUATIONS: usize = 1 + EXPECTED_GENERATIONS * LAMBDA;
const SIGMA0: f64 = 0.4;
const IMPOSSIBLE_TARGET: f64 = -1.0;
const EIGEN_INTERVAL: usize = 1;
const ORACLE_TOLERANCE_FACTOR: f64 = 64.0;

const MUTATION_SEED: u64 = 0xC0A5_FA11_0000_0041;
const MUTATION_KERNEL: u32 = 0xC041;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION > 0);
const _: () = assert!(LAMBDA >= 4);
const _: () = assert!(MAX_EVALUATIONS == 65);
const _: () = assert!(OBJECTIVE_OFFSET > IMPOSSIBLE_TARGET);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    value: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ReportBits {
    x_best: Vec<u64>,
    f_best: u64,
    evals: usize,
    generations: usize,
    converged: bool,
    sigma: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    evaluations: Vec<EvaluationBits>,
    report: ReportBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    fixture_digest: ContentHash,
    record: StudyRecord,
    result: ReplayIdentity,
    result_digest: ContentHash,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    ReferenceIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
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

fn digest_bytes(digest: ContentHash) -> [u8; 32] {
    *digest.as_bytes()
}

fn params() -> CmaParams {
    CmaParams {
        lambda: LAMBDA,
        sigma0: SIGMA0,
        max_evals: MAX_EVALUATIONS,
        f_target: IMPOSSIBLE_TARGET,
        eigen_interval: EIGEN_INTERVAL,
    }
}

/// Objective implementation passed to production CMA. The explicit `mul_add`
/// chain fixes term order and avoids an ambient reassociation assumption.
fn objective_under_test(point: &[f64]) -> f64 {
    assert_eq!(point.len(), DIMENSION);
    let mut value = OBJECTIVE_OFFSET;
    for ((&coordinate, &center), &weight) in point.iter().zip(&CENTER).zip(&WEIGHTS) {
        let residual = coordinate - center;
        value = weight.mul_add(residual * residual, value);
    }
    value
}

/// Algebraically expanded objective oracle used only by retained-record
/// admission. This deliberately does not share the producer's residual-square
/// loop or fused operations.
fn reference_objective(point: &[f64]) -> Option<f64> {
    let [x0, x1, x2, x3] = point else {
        return None;
    };
    let q0 = *x0 * *x0 - 0.5 * *x0;
    let q1 = 2.0 * *x1 * *x1 + 2.0 * *x1;
    let q2 = 4.0 * *x2 * *x2 - 6.0 * *x2;
    let q3 = 8.0 * *x3 * *x3 + 16.0 * *x3;
    Some((q0 + q1) + (q2 + q3) + 10.9375)
}

/// Forward-error allowance for the expanded oracle, scaled by the absolute
/// terms in that independent expression rather than by the possibly
/// cancellation-reduced result alone.
fn reference_tolerance(point: &[f64], produced: f64, reference: f64) -> Option<f64> {
    let [x0, x1, x2, x3] = point else {
        return None;
    };
    let expanded_scale = 10.9375
        + (*x0 * *x0).abs()
        + (0.5 * *x0).abs()
        + (2.0 * *x1 * *x1).abs()
        + (2.0 * *x1).abs()
        + (4.0 * *x2 * *x2).abs()
        + (6.0 * *x2).abs()
        + (8.0 * *x3 * *x3).abs()
        + (16.0 * *x3).abs();
    Some(
        ORACLE_TOLERANCE_FACTOR
            * f64::EPSILON
            * expanded_scale
                .max(produced.abs())
                .max(reference.abs())
                .max(1.0),
    )
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

fn report_bits(report: CmaReport) -> ReportBits {
    ReportBits {
        x_best: report.x_best.into_iter().map(f64::to_bits).collect(),
        f_best: report.f_best.to_bits(),
        evals: report.evals,
        generations: report.generations,
        converged: report.converged,
        sigma: report.sigma.to_bits(),
    }
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_dfo::cmaes/full-covariance")
        .str("objective", "weighted-shifted-quadratic-mul-add-v1")
        .str(
            "objective-producer-formula",
            "offset;for-i:value=weight[i].mul_add((x[i]-center[i])^2,value)",
        )
        .str(
            "objective-oracle-formula",
            "expanded-diagonal-quadratic-grouped-(q0+q1)+(q2+q3)+10.9375",
        )
        .str(
            "objective-oracle-bound",
            "64*EPSILON*max(1,abs(produced),abs(oracle),sum-abs-expanded-terms)",
        )
        .f64_bits("objective-oracle-tolerance-factor", ORACLE_TOLERANCE_FACTOR)
        .str("coordinate-units", "dimensionless")
        .str("objective-units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .u64("seed", INPUT_SEED)
        .u64("lambda", usize_u64(LAMBDA))
        .f64_bits("sigma0", SIGMA0)
        .u64("max-evaluations", usize_u64(MAX_EVALUATIONS))
        .f64_bits("impossible-target", IMPOSSIBLE_TARGET)
        .u64("eigen-interval", usize_u64(EIGEN_INTERVAL))
        .str("ranking", "f64-total-cmp-then-lowest-candidate-index")
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-la-version", fs_la::VERSION)
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
        .str("fixture-digest-domain", FIXTURE_DIGEST_DOMAIN)
        .str("result-digest-domain", RESULT_DIGEST_DOMAIN)
        .str("event-digest-domain", EVENT_DIGEST_DOMAIN)
        .str(
            "no-claims",
            "quality;independent-CMA-oracle;internal-CMA-state;arbitrary-objectives;arbitrary-dimensions;arbitrary-params;arbitrary-seeds;legacy-api-validation;cross-ISA;Cx;checkpoint;authenticated-ledger;performance",
        );
    for (index, ((start, center), weight)) in START.into_iter().zip(CENTER).zip(WEIGHTS).enumerate()
    {
        builder = builder
            .u64("coordinate-index", usize_u64(index))
            .f64_bits("start-coordinate", start)
            .f64_bits("objective-center", center)
            .f64_bits("objective-weight", weight);
    }
    builder
        .f64_bits("objective-offset", OBJECTIVE_OFFSET)
        .finish()
}

fn fixture_digest(fixture: &ReplayIdentity) -> ContentHash {
    hash_domain(FIXTURE_DIGEST_DOMAIN, fixture.canonical_bytes())
}

fn result_identity(
    fixture: &ReplayIdentity,
    strong_fixture: ContentHash,
    record: &StudyRecord,
) -> ReplayIdentity {
    let report = &record.report;
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .child("fixture-compatibility-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .bytes("fixture-blake3", strong_fixture.as_bytes())
        .u64("evaluation-count", usize_u64(record.evaluations.len()))
        .u64("reported-evals", usize_u64(report.evals))
        .u64("reported-generations", usize_u64(report.generations))
        .flag("reported-converged", report.converged)
        .f64_bits("reported-f-best", f64::from_bits(report.f_best))
        .f64_bits("reported-sigma", f64::from_bits(report.sigma))
        .u64("returned-point-length", usize_u64(report.x_best.len()));
    for (evaluation_index, evaluation) in record.evaluations.iter().enumerate() {
        builder = builder
            .u64("evaluation-index", usize_u64(evaluation_index))
            .u64("evaluation-point-length", usize_u64(evaluation.point.len()));
        for (coordinate, &bits) in evaluation.point.iter().enumerate() {
            builder = builder
                .u64("evaluation-coordinate-index", usize_u64(coordinate))
                .f64_bits("evaluation-coordinate", f64::from_bits(bits));
        }
        builder = builder.f64_bits("evaluation-objective", f64::from_bits(evaluation.value));
    }
    for (coordinate, &bits) in report.x_best.iter().enumerate() {
        builder = builder
            .u64("returned-coordinate-index", usize_u64(coordinate))
            .f64_bits("returned-coordinate", f64::from_bits(bits));
    }
    builder.finish()
}

fn result_digest(result: &ReplayIdentity) -> ContentHash {
    hash_domain(RESULT_DIGEST_DOMAIN, result.canonical_bytes())
}

fn event_digest(event: &Event) -> ContentHash {
    let identity = event.content_identity();
    hash_domain(EVENT_DIGEST_DOMAIN, identity.canonical_bytes())
}

fn run_study() -> StudyRun {
    let mut evaluations = Vec::with_capacity(MAX_EVALUATIONS);
    let report = {
        let mut objective = |point: &[f64]| {
            let value = objective_under_test(point);
            evaluations.push(evaluation_bits(point, value));
            value
        };
        cmaes(&mut objective, &START, &params(), INPUT_SEED)
    };
    let record = StudyRecord {
        evaluations,
        report: report_bits(report),
    };
    let fixture = fixture_identity();
    let fixture_digest = fixture_digest(&fixture);
    let result = result_identity(&fixture, fixture_digest, &record);
    let result_digest = result_digest(&result);
    StudyRun {
        fixture,
        fixture_digest,
        record,
        result,
        result_digest,
    }
}

#[allow(clippy::too_many_lines)] // Every retained callback/report field is admitted here.
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    let report = &record.report;
    if record.evaluations.len() != report.evals {
        return Some(format!(
            "reported-evals:{}!=closure-count:{}",
            report.evals,
            record.evaluations.len()
        ));
    }
    let expected_evals = report
        .generations
        .checked_mul(LAMBDA)
        .and_then(|population_evals| population_evals.checked_add(1));
    if expected_evals != Some(report.evals) {
        return Some(format!(
            "generation-accounting:1+{}*{LAMBDA}!={}",
            report.generations, report.evals
        ));
    }
    if report.evals != MAX_EVALUATIONS || report.generations != EXPECTED_GENERATIONS {
        return Some(format!(
            "fixed-budget-shape:evals={};generations={}",
            report.evals, report.generations
        ));
    }
    let sigma = f64::from_bits(report.sigma);
    if !sigma.is_finite() || sigma <= 0.0 {
        return Some(format!("invalid-final-sigma:0x{:016x}", report.sigma));
    }
    if report.x_best.len() != DIMENSION {
        return Some(format!(
            "returned-point-dimension:{}!=expected-{DIMENSION}",
            report.x_best.len()
        ));
    }
    if record.evaluations.is_empty() {
        return Some("empty-callback-trace".to_string());
    }
    let expected_start: Vec<u64> = START.into_iter().map(f64::to_bits).collect();
    if record.evaluations[0].point != expected_start {
        return Some("first-callback-is-not-the-bound-start-point".to_string());
    }

    let mut earliest_best = 0usize;
    let mut reached_target = false;
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
        if !value.is_finite() {
            return Some(format!(
                "evaluation[{index}]-invalid-objective:0x{:016x}",
                evaluation.value
            ));
        }
        reached_target |= value <= IMPOSSIBLE_TARGET;
        let reproduced = objective_under_test(&point);
        if reproduced.to_bits() != evaluation.value {
            return Some(format!(
                "evaluation[{index}]-producer-replay:recomputed=0x{:016x};recorded=0x{:016x}",
                reproduced.to_bits(),
                evaluation.value
            ));
        }
        let Some(recomputed) = reference_objective(&point) else {
            return Some(format!("evaluation[{index}]-oracle-shape-refusal"));
        };
        let Some(tolerance) = reference_tolerance(&point, value, recomputed) else {
            return Some(format!("evaluation[{index}]-oracle-bound-shape-refusal"));
        };
        let error = (recomputed - value).abs();
        if !recomputed.is_finite() || !tolerance.is_finite() || error > tolerance {
            return Some(format!(
                "evaluation[{index}]-oracle:error=0x{:016x};tolerance=0x{:016x};oracle=0x{:016x};recorded=0x{:016x}",
                error.to_bits(),
                tolerance.to_bits(),
                recomputed.to_bits(),
                evaluation.value
            ));
        }
        if index > 0
            && value
                .total_cmp(&f64::from_bits(record.evaluations[earliest_best].value))
                .is_lt()
        {
            earliest_best = index;
        }
    }
    if report.converged != reached_target {
        return Some(format!(
            "convergence-witness:reported={};trace-witness={reached_target}",
            report.converged
        ));
    }
    if reached_target {
        return Some("impossible-target-was-reached".to_string());
    }

    let expected_best = &record.evaluations[earliest_best];
    if report.f_best != expected_best.value {
        return Some(format!(
            "returned-objective:0x{:016x}!=earliest-trace-minimum[{}]=0x{:016x}",
            report.f_best, earliest_best, expected_best.value
        ));
    }
    if report.x_best != expected_best.point {
        return Some(format!(
            "returned-point!=earliest-trace-minimum[{earliest_best}]"
        ));
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed_fixture = fixture_digest(&run.fixture);
    if computed_fixture != run.fixture_digest {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: digest_bytes(run.fixture_digest),
            computed: digest_bytes(computed_fixture),
        });
    }
    let computed_result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    let computed_digest = result_digest(&computed_result);
    if computed_result.canonical_bytes() != run.result.canonical_bytes()
        || computed_digest != run.result_digest
    {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: digest_bytes(run.result_digest),
            computed: digest_bytes(computed_digest),
        });
    }
    Ok(())
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    match semantic_mismatch(&run.record) {
        Some(mismatch) => Err(AdmissionError::SemanticInconsistency(mismatch)),
        None => Ok(()),
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.result.canonical_bytes()
        && run.result_digest == reference.result_digest
    {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: digest_bytes(reference.result_digest),
            found: digest_bytes(run.result_digest),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    run.result_digest = result_digest(&run.result);
}

fn exact_returned_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(&reference_bits) = reference.record.report.x_best.get(mutation.coordinate) else {
        return false;
    };
    let Some(&mutant_bits) = mutant.record.report.x_best.get(mutation.coordinate) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference.fixture_digest != mutant.fixture_digest
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.report.x_best[mutation.coordinate] = mutation.after;
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
    let before = run.record.report.x_best[coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.report.x_best[coordinate] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed mutation must not match retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed returned-point mutation must remain semantically invalid");
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
    let report = &run.record.report;
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "cma-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"fixture_blake3\":\"{}\",",
                    "\"result_identity\":\"{}\",\"result_blake3\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::cmaes\",\"seed\":{},",
                    "\"dimension\":{},\"lambda\":{},\"max_evaluations\":{},",
                    "\"actual_evaluations\":{},\"generations\":{},",
                    "\"converged\":{},\"returned_objective_bits\":\"0x{:016x}\",",
                    "\"final_sigma_bits\":\"0x{:016x}\",",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_la\":\"{}\",",
                    "\"fs_math\":\"{}\",\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"optimizer-quality\",\"arbitrary-fixtures\",",
                    "\"legacy-api-validation\",\"cross-ISA\",\"cancellation\",",
                    "\"internal-CMA-state\",\"checkpointing\",",
                    "\"authenticated-ledger\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.fixture_digest.to_hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                INPUT_SEED,
                DIMENSION,
                LAMBDA,
                MAX_EVALUATIONS,
                report.evals,
                report.generations,
                report.converged,
                report.f_best,
                report.sigma,
                fs_dfo::VERSION,
                fs_la::VERSION,
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
                "fixture={}; result={}; blake3={}; callbacks={}; report=fully-bound",
                run.fixture.hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                run.record.evaluations.len(),
            ),
            seed: INPUT_SEED,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.x_best[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result_digest.to_hex(),
        corruption.run.result_digest.to_hex(),
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

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay plus all three refusal gates.
fn standalone_cma_full_study_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(validate_semantics(&original), Ok(()));
    assert_eq!(validate_semantics(&replay), Ok(()));
    assert_eq!(admit_reference(&original, &replay), Ok(()));
    assert_eq!(admit_reference(&replay, &original), Ok(()));
    assert_eq!(original.record, replay.record);
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.fixture_digest, replay.fixture_digest);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.result_digest, replay.result_digest);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete result frames must replay byte-for-byte"
    );

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    assert_eq!(
        first_receipt.content_identity().canonical_bytes(),
        second_receipt.content_identity().canonical_bytes(),
        "green CMA receipt content must replay byte-for-byte"
    );
    assert_eq!(event_digest(&first_receipt), event_digest(&second_receipt));
    for event in [&first_receipt, &second_receipt] {
        fs_obs::lint_failure_record(event).expect("CMA receipt must retain replay inputs");
        fs_obs::validate_line(&event.to_jsonl())
            .expect("CMA receipt must use the fs-obs wire schema");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("fresh CMA receipt content identity must admit");
    }
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    assert_eq!(
        first_green.content_identity().canonical_bytes(),
        second_green.content_identity().canonical_bytes(),
        "green CMA verdict content must replay byte-for-byte"
    );
    assert_eq!(event_digest(&first_green), event_digest(&second_green));
    for event in [&first_green, &second_green] {
        fs_obs::lint_failure_record(event).expect("green CMA verdict retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("green CMA verdict is wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("green CMA verdict content identity admits exactly");
        assert_mergeable(event);
    }
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded corruption must replay exactly");
    assert!(
        exact_returned_bit_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one retained returned-coordinate bit"
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
            if declared == original.result_digest.as_bytes()
                && computed == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result_digest.as_bytes()
                && found == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.semantic_error,
        AdmissionError::SemanticInconsistency(mismatch)
            if mismatch.starts_with("returned-point!=earliest-trace-minimum")
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_eq!(
        first_red.content_identity().canonical_bytes(),
        second_red.content_identity().canonical_bytes(),
        "red CMA evidence content must replay byte-for-byte"
    );
    assert_eq!(event_digest(&first_red), event_digest(&second_red));
    for event in [&first_red, &second_red] {
        fs_obs::lint_failure_record(event).expect("red CMA evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("red CMA evidence is wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("red CMA evidence content identity admits exactly");
    }
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded CMA corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!("report.x_best[{}]", first.mutation.coordinate)));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
