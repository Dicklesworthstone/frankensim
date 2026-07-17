//! G5 study-scale replay for the production BIPOP-CMA path (7tv.21.22).
//!
//! The fixture captures every objective call and binds that trace together
//! with every public `BipopReport`/`CmaReport` field.  A disclosed seeded
//! mutation changes one returned coordinate bit. The unsealed edit is refused
//! as a stale payload, the self-consistently resealed edit is refused against
//! the retained reference identity, and the resulting red fs-obs evidence is
//! independently reproducible. This is one finite deterministic study, not an
//! optimizer-quality or performance claim.

use fs_dfo::{BipopReport, bipop_cmaes};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use std::fmt::Write as _;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/bipop-study-replay";
const CASE: &str = "shifted-rastrigin-4d-full-public-state";
const RED_CASE: &str = "seeded-returned-coordinate-corruption";

const INPUT_SEED: u64 = 0xDF0A_2100_0000_0001;
const CORRUPTION_SEED: u64 = 0xDF0A_F11E_0000_0001;
const DIMENSION: usize = 4;
const X0: [f64; DIMENSION] = [2.75, -1.25, 3.5, -2.0];
const SHIFT: [f64; DIMENSION] = [0.25, -0.5, 1.0, -1.5];
const SIGMA0: f64 = 1.25;
const TOTAL_BUDGET: usize = 6_000;
// Shifted Rastrigin is non-negative, so this target keeps every restart in
// the non-converged evidence path without depending on solution quality.
const F_TARGET: f64 = -1.0;
const BASE_LAMBDA: usize = 8;

// These are the logical stream coordinates and restart rule used by
// `fs_dfo::cma`; recording them makes the private implementation choice
// explicit in the fixture identity.  A change also changes the captured trace.
const CMA_STREAM_KERNEL: u64 = 0xD1F0;
const CMA_SAMPLE_TILE: u64 = 0;
const CMA_RESTART_TILE: u64 = 1;
const RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

#[derive(Debug, Clone)]
struct Evaluation {
    x: Vec<f64>,
    value: f64,
}

#[derive(Debug, Clone)]
struct StudyRun {
    fixture: ReplayIdentity,
    report: BipopReport,
    evaluations: Vec<Evaluation>,
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
    coordinate: usize,
    mantissa_bit: u32,
    before: u64,
    after: u64,
}

#[derive(Debug)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    mismatch: String,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn shifted_rastrigin(x: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        DIMENSION,
        "fixture dimension is part of the contract"
    );
    x.iter()
        .zip(SHIFT)
        .map(|(&value, shift)| {
            let z = value - shift;
            10.0 + z.mul_add(z, -10.0 * fs_math::det::cos(std::f64::consts::TAU * z))
        })
        .sum()
}

fn fixture_identity(seed: u64) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-bipop-study-fixture-v1")
        .str("algorithm", "fs_dfo::bipop_cmaes")
        .str("objective", "shifted-rastrigin")
        .str("units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .f64_bits("sigma0", SIGMA0)
        .u64("total-evaluation-budget", usize_u64(TOTAL_BUDGET))
        .f64_bits("f-target", F_TARGET)
        .u64("input-seed", seed)
        .u64("base-lambda", usize_u64(BASE_LAMBDA))
        .str("base-lambda-rule", "4+floor(3*ln(dimension))")
        .str("per-restart-budget-rule", "min(lambda*250,remaining)")
        .str("large-restart-rule", "large-budget-used<=small-budget-used")
        .str("large-population-rule", "base-lambda*2^large-runs")
        .str("stagnation-rule", "tol-x-or-120-generation-tol-f")
        .u64("cma-stream-kernel", CMA_STREAM_KERNEL)
        .u64("sample-stream-tile", CMA_SAMPLE_TILE)
        .u64("restart-stream-tile", CMA_RESTART_TILE)
        .u64("restart-seed-stride", RESTART_SEED_STRIDE)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "fs-rand-stream-position-domain",
            fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
        )
        .str(
            "capabilities",
            "safe-rust;strict-fs-math;keyed-fs-rand;canonical-fs-obs",
        )
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-la-version", fs_la::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str("fs-obs-version", fs_obs::VERSION);
    for (coordinate, (&x0, &shift)) in X0.iter().zip(&SHIFT).enumerate() {
        builder = builder
            .u64("coordinate-index", usize_u64(coordinate))
            .f64_bits("initial-coordinate", x0)
            .f64_bits("objective-shift", shift);
    }
    builder.finish()
}

fn result_identity(
    fixture: &ReplayIdentity,
    report: &BipopReport,
    evaluations: &[Evaluation],
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-bipop-study-result-v1")
        .child("fixture", fixture)
        .u64("total-evals", usize_u64(report.total_evals))
        .u64("schedule-length", usize_u64(report.schedule.len()))
        .u64("best-x-length", usize_u64(report.best.x_best.len()))
        .f64_bits("best-f", report.best.f_best)
        .u64("best-run-evals", usize_u64(report.best.evals))
        .u64("best-run-generations", usize_u64(report.best.generations))
        .flag("best-run-converged", report.best.converged)
        .f64_bits("best-run-sigma", report.best.sigma)
        .u64("evaluation-trace-length", usize_u64(evaluations.len()));
    for (restart, &lambda) in report.schedule.iter().enumerate() {
        builder = builder
            .u64("restart-index", usize_u64(restart))
            .u64("restart-lambda", usize_u64(lambda));
    }
    for (coordinate, &value) in report.best.x_best.iter().enumerate() {
        builder = builder
            .u64("best-coordinate-index", usize_u64(coordinate))
            .f64_bits("best-coordinate", value);
    }
    for (evaluation_index, evaluation) in evaluations.iter().enumerate() {
        builder = builder
            .u64("evaluation-index", usize_u64(evaluation_index))
            .u64("evaluation-dimension", usize_u64(evaluation.x.len()));
        for (coordinate, &value) in evaluation.x.iter().enumerate() {
            builder = builder
                .u64("evaluation-coordinate-index", usize_u64(coordinate))
                .f64_bits("evaluation-coordinate", value);
        }
        builder = builder.f64_bits("evaluation-objective", evaluation.value);
    }
    builder.finish()
}

fn run_study(seed: u64) -> StudyRun {
    let mut evaluations = Vec::with_capacity(TOTAL_BUDGET);
    let report = {
        let mut objective = |x: &[f64]| {
            let value = shifted_rastrigin(x);
            evaluations.push(Evaluation {
                x: x.to_vec(),
                value,
            });
            value
        };
        bipop_cmaes(&mut objective, &X0, SIGMA0, TOTAL_BUDGET, F_TARGET, seed)
    };
    let fixture = fixture_identity(seed);
    let result = result_identity(&fixture, &report, &evaluations);
    StudyRun {
        fixture,
        report,
        evaluations,
        result,
    }
}

fn same_point_bits(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(a, b)| a.to_bits() == b.to_bits())
}

#[allow(clippy::too_many_lines)] // Complete trace and public-report accounting is the oracle.
fn accounting_mismatch(run: &StudyRun) -> Option<String> {
    if run.report.total_evals != run.evaluations.len() {
        return Some(format!(
            "reported-total-evals:{}!=closure-count:{}",
            run.report.total_evals,
            run.evaluations.len()
        ));
    }
    if !(1..=TOTAL_BUDGET).contains(&run.report.total_evals) {
        return Some(format!(
            "total-evals:{} not in 1..={TOTAL_BUDGET}",
            run.report.total_evals
        ));
    }
    if run.report.schedule.len() < 2 {
        return Some(format!(
            "schedule-too-short:{};fixture-must-exercise-a-restart",
            run.report.schedule.len()
        ));
    }
    if run.report.schedule.len() > run.report.total_evals {
        return Some(format!(
            "schedule-length:{}>total-evals:{}",
            run.report.schedule.len(),
            run.report.total_evals
        ));
    }
    if run.report.schedule[..2] != [BASE_LAMBDA, BASE_LAMBDA] {
        return Some(format!(
            "first-two-populations:{:?}!=[{BASE_LAMBDA},{BASE_LAMBDA}]",
            &run.report.schedule[..2]
        ));
    }
    for (restart, &lambda) in run.report.schedule.iter().enumerate() {
        if lambda < BASE_LAMBDA
            || lambda % BASE_LAMBDA != 0
            || !(lambda / BASE_LAMBDA).is_power_of_two()
        {
            return Some(format!(
                "schedule[{restart}]={lambda};expected-base-times-power-of-two"
            ));
        }
    }

    let best = &run.report.best;
    if best.x_best.len() != DIMENSION {
        return Some(format!(
            "best-point-dimension:{}!=expected-{DIMENSION}",
            best.x_best.len()
        ));
    }
    if best.x_best.iter().any(|value| !value.is_finite()) {
        return Some(format!(
            "best-point-non-finite:{:016x?}",
            best.x_best
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        ));
    }
    if !best.f_best.is_finite() || best.f_best <= F_TARGET {
        return Some(format!(
            "best-objective-invalid:0x{:016x};target=0x{:016x}",
            best.f_best.to_bits(),
            F_TARGET.to_bits()
        ));
    }
    if !(1..=run.report.total_evals).contains(&best.evals) {
        return Some(format!(
            "best-run-evals:{} not in 1..={}",
            best.evals, run.report.total_evals
        ));
    }
    let generation_accounting_matches = run.report.schedule.iter().any(|&lambda| {
        best.generations
            .checked_mul(lambda)
            .and_then(|population_evals| population_evals.checked_add(1))
            == Some(best.evals)
    });
    if !generation_accounting_matches {
        return Some(format!(
            "best-run-accounting:evals-{}!=1+generations-{}*any-scheduled-population",
            best.evals, best.generations
        ));
    }
    if best.converged {
        return Some("impossible-target-was-reported-converged".to_string());
    }
    if !best.sigma.is_finite() || best.sigma <= 0.0 {
        return Some(format!(
            "best-run-invalid-sigma:0x{:016x}",
            best.sigma.to_bits()
        ));
    }

    let mut trace_minimum: Option<f64> = None;
    for (evaluation_index, evaluation) in run.evaluations.iter().enumerate() {
        if evaluation.x.len() != DIMENSION {
            return Some(format!(
                "trace[{evaluation_index}]-dimension:{}!=expected-{DIMENSION}",
                evaluation.x.len()
            ));
        }
        if evaluation.x.iter().any(|value| !value.is_finite()) {
            return Some(format!(
                "trace[{evaluation_index}]-non-finite-point:{:016x?}",
                evaluation
                    .x
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            ));
        }
        if !evaluation.value.is_finite() {
            return Some(format!(
                "trace[{evaluation_index}]-non-finite-objective:0x{:016x}",
                evaluation.value.to_bits()
            ));
        }
        let recomputed = shifted_rastrigin(&evaluation.x).to_bits();
        if evaluation.value.to_bits() != recomputed {
            return Some(format!(
                "trace[{evaluation_index}]-objective:recorded=0x{:016x};recomputed=0x{recomputed:016x}",
                evaluation.value.to_bits()
            ));
        }
        trace_minimum = Some(match trace_minimum {
            Some(current) if current <= evaluation.value => current,
            _ => evaluation.value,
        });
    }
    let trace_minimum = trace_minimum.expect("positive total-eval accounting makes trace nonempty");
    if trace_minimum.to_bits() != best.f_best.to_bits() {
        return Some(format!(
            "complete-trace-minimum=0x{:016x};reported-best=0x{:016x}",
            trace_minimum.to_bits(),
            best.f_best.to_bits()
        ));
    }
    let recomputed_best = shifted_rastrigin(&best.x_best).to_bits();
    if recomputed_best != best.f_best.to_bits() {
        return Some(format!(
            "best-point-objective:recomputed=0x{recomputed_best:016x};reported=0x{:016x}",
            best.f_best.to_bits()
        ));
    }
    if !run.evaluations.iter().any(|evaluation| {
        evaluation.value.to_bits() == best.f_best.to_bits()
            && same_point_bits(&evaluation.x, &best.x_best)
    }) {
        return Some("reported-best-is-not-an-evaluated-point".to_string());
    }
    None
}

#[allow(clippy::too_many_lines)] // Exhaustive field-by-field public-state audit.
fn first_public_mismatch(left: &StudyRun, right: &StudyRun) -> Option<String> {
    if left.fixture.canonical_bytes() != right.fixture.canonical_bytes() {
        return Some("fixture-identity".to_string());
    }
    if left.report.total_evals != right.report.total_evals {
        return Some(format!(
            "total-evals:{}!={}",
            left.report.total_evals, right.report.total_evals
        ));
    }
    if left.report.schedule.len() != right.report.schedule.len() {
        return Some(format!(
            "schedule-length:{}!={}",
            left.report.schedule.len(),
            right.report.schedule.len()
        ));
    }
    for (restart, (&a, &b)) in left
        .report
        .schedule
        .iter()
        .zip(&right.report.schedule)
        .enumerate()
    {
        if a != b {
            return Some(format!("schedule[{restart}]:{a}!={b}"));
        }
    }

    let a = &left.report.best;
    let b = &right.report.best;
    if a.x_best.len() != b.x_best.len() {
        return Some(format!(
            "best.x-length:{}!={}",
            a.x_best.len(),
            b.x_best.len()
        ));
    }
    for (coordinate, (&x, &y)) in a.x_best.iter().zip(&b.x_best).enumerate() {
        if x.to_bits() != y.to_bits() {
            return Some(format!(
                "best.x[{coordinate}]:0x{:016x}!=0x{:016x}",
                x.to_bits(),
                y.to_bits()
            ));
        }
    }
    if a.f_best.to_bits() != b.f_best.to_bits() {
        return Some(format!(
            "best.f:0x{:016x}!=0x{:016x}",
            a.f_best.to_bits(),
            b.f_best.to_bits()
        ));
    }
    if a.evals != b.evals {
        return Some(format!("best.evals:{}!={}", a.evals, b.evals));
    }
    if a.generations != b.generations {
        return Some(format!(
            "best.generations:{}!={}",
            a.generations, b.generations
        ));
    }
    if a.converged != b.converged {
        return Some(format!("best.converged:{}!={}", a.converged, b.converged));
    }
    if a.sigma.to_bits() != b.sigma.to_bits() {
        return Some(format!(
            "best.sigma:0x{:016x}!=0x{:016x}",
            a.sigma.to_bits(),
            b.sigma.to_bits()
        ));
    }

    if left.evaluations.len() != right.evaluations.len() {
        return Some(format!(
            "trace-length:{}!={}",
            left.evaluations.len(),
            right.evaluations.len()
        ));
    }
    for (evaluation_index, (a, b)) in left.evaluations.iter().zip(&right.evaluations).enumerate() {
        if a.x.len() != b.x.len() {
            return Some(format!(
                "trace[{evaluation_index}].x-length:{}!={}",
                a.x.len(),
                b.x.len()
            ));
        }
        for (coordinate, (&x, &y)) in a.x.iter().zip(&b.x).enumerate() {
            if x.to_bits() != y.to_bits() {
                return Some(format!(
                    "trace[{evaluation_index}].x[{coordinate}]:0x{:016x}!=0x{:016x}",
                    x.to_bits(),
                    y.to_bits()
                ));
            }
        }
        if a.value.to_bits() != b.value.to_bits() {
            return Some(format!(
                "trace[{evaluation_index}].f:0x{:016x}!=0x{:016x}",
                a.value.to_bits(),
                b.value.to_bits()
            ));
        }
    }
    if left.result.canonical_bytes() != right.result.canonical_bytes() {
        return Some("result-identity".to_string());
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed = result_identity(&run.fixture, &run.report, &run.evaluations);
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
    let Some(reference_coordinate) = reference.report.best.x_best.get(mutation.coordinate) else {
        return false;
    };
    let Some(mutant_coordinate) = mutant.report.best.x_best.get(mutation.coordinate) else {
        return false;
    };
    if reference_coordinate.to_bits() != mutation.before
        || mutant_coordinate.to_bits() != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }

    let mut expected = reference.clone();
    expected.report.best.x_best[mutation.coordinate] = f64::from_bits(mutation.after);
    expected.result = result_identity(&expected.fixture, &expected.report, &expected.evaluations);
    first_public_mismatch(&expected, mutant).is_none()
}

fn schedule_json(schedule: &[usize]) -> String {
    let mut json = String::from("[");
    for (index, lambda) in schedule.iter().enumerate() {
        if index > 0 {
            json.push(',');
        }
        write!(&mut json, "{lambda}").expect("String writes are infallible");
    }
    json.push(']');
    json
}

fn emit_green_receipt(run: &StudyRun) {
    let mut emitter = Emitter::new(SUITE, CASE);
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "bipop-cma-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::bipop_cmaes\",\"objective\":\"shifted-rastrigin\",",
                    "\"units\":\"dimensionless\",\"input_seed\":{},\"dimension\":{},",
                    "\"total_budget\":{},\"total_evals\":{},\"schedule\":{},",
                    "\"best\":{{\"x_len\":{},\"f_bits\":\"0x{:016x}\",",
                    "\"evals\":{},\"generations\":{},\"converged\":{},",
                    "\"sigma_bits\":\"0x{:016x}\"}},\"trace_len\":{},",
                    "\"stream_semantics_version\":{},\"versions\":{{",
                    "\"fs_dfo\":\"{}\",\"fs_la\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_rand\":\"{}\",\"fs_obs\":\"{}\"}},",
                    "\"no_claims\":[\"optimizer-quality\",\"all-objectives\",",
                    "\"all-dimensions\",\"all-budgets\",\"all-seeds\",",
                    "\"cross-ISA-equality\",\"cancellation\",\"checkpointing\",",
                    "\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                INPUT_SEED,
                DIMENSION,
                TOTAL_BUDGET,
                run.report.total_evals,
                schedule_json(&run.report.schedule),
                run.report.best.x_best.len(),
                run.report.best.f_best.to_bits(),
                run.report.best.evals,
                run.report.best.generations,
                run.report.best.converged,
                run.report.best.sigma.to_bits(),
                run.evaluations.len(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                fs_dfo::VERSION,
                fs_la::VERSION,
                fs_math::VERSION,
                fs_rand::VERSION,
                fs_obs::VERSION,
            ),
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("BIPOP study receipt must use the fs-obs wire schema");
    let receipt = event.content_identity_receipt();
    event
        .admit_content_identity(&receipt)
        .expect("fresh retained event identity must admit exactly");
    println!("{line}");
}

fn emit_green_verdict(run: &StudyRun) -> Event {
    let detail = format!(
        "fixture={}; result={}; total_evals={}; restarts={}; trace=bit-exact; public_report=fully-bound",
        run.fixture.hex(),
        run.result.hex(),
        run.report.total_evals,
        run.report.schedule.len()
    );
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail,
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("BIPOP study verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("BIPOP study verdict must use the fs-obs wire schema");
    println!("{line}");
    event
}

fn failure_event(detail: &str, corruption_seed: u64) -> Event {
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail: detail.to_string(),
            seed: corruption_seed,
        },
        None,
    )
}

fn assert_mergeable(run: &StudyRun, reference: &ReplayIdentity, event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    if let Err(error) = admit_against(run, reference) {
        panic!("merge gate refused {case}: {error:?}; {detail}");
    }
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn seeded_corruption(canonical: &StudyRun, seed: u64) -> SeededCorruption {
    let coordinate =
        usize::try_from(seed % usize_u64(DIMENSION)).expect("corruption coordinate fits usize");
    let mantissa_bit = u32::try_from((seed >> 32) & 0x1f).expect("corruption bit fits u32");

    let mut run = canonical.clone();
    let before = run.report.best.x_best[coordinate].to_bits();
    let after = before ^ (1_u64 << mantissa_bit);
    run.report.best.x_best[coordinate] = f64::from_bits(after);
    assert!(run.report.best.x_best[coordinate].is_finite());
    let stale_error = validate_payload(&run).expect_err("unsealed result mutation must refuse");
    run.result = result_identity(&run.fixture, &run.report, &run.evaluations);
    validate_payload(&run).expect("resealed mutation must be internally self-consistent");
    let reference_error = admit_against(&run, &canonical.result)
        .expect_err("resealed semantic mutation must not match the retained reference");

    let mismatch = first_public_mismatch(canonical, &run)
        .expect("the disclosed mutation must change public replay state");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed,
            coordinate,
            mantissa_bit,
            before,
            after,
        },
        stale_error,
        reference_error,
        mismatch,
    }
}

fn corruption_detail(canonical: &StudyRun, corruption: &SeededCorruption) -> String {
    format!(
        "input_seed=0x{INPUT_SEED:016x}; corruption_seed=0x{:016x}; fixture={}; coordinate={}; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}; canonical={}; corrupted={}",
        corruption.mutation.seed,
        canonical.fixture.hex(),
        corruption.mutation.coordinate,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.mismatch,
        canonical.result.hex(),
        corruption.run.result.hex()
    )
}

fn exercise_disclosed_corruption(canonical: &StudyRun, replay: &StudyRun) {
    let first_corruption = seeded_corruption(canonical, CORRUPTION_SEED);
    let replay_corruption = seeded_corruption(replay, CORRUPTION_SEED);
    assert_eq!(
        (
            first_corruption.mutation.coordinate,
            first_corruption.mutation.mantissa_bit
        ),
        (1, 30)
    );
    assert!(
        first_corruption.mismatch.starts_with("best.x[1]"),
        "unexpected mismatch: {}",
        first_corruption.mismatch
    );
    assert_eq!(first_corruption.mutation, replay_corruption.mutation);
    assert_eq!(first_corruption.stale_error, replay_corruption.stale_error);
    assert_eq!(
        first_corruption.reference_error,
        replay_corruption.reference_error
    );
    assert_eq!(first_corruption.mismatch, replay_corruption.mismatch);
    assert!(exact_returned_bit_delta(
        canonical,
        &first_corruption.run,
        first_corruption.mutation
    ));
    assert!(exact_returned_bit_delta(
        replay,
        &replay_corruption.run,
        replay_corruption.mutation
    ));
    assert!(matches!(
        first_corruption.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == canonical.result.root()
                && computed == first_corruption.run.result.root()
    ));
    assert!(matches!(
        first_corruption.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == canonical.result.root()
                && found == first_corruption.run.result.root()
    ));
    assert_eq!(validate_payload(&first_corruption.run), Ok(()));
    assert!(matches!(
        admit_against(&first_corruption.run, &canonical.result),
        Err(AdmissionError::ReferenceIdentityMismatch { expected, found })
            if expected == canonical.result.root()
                && found == first_corruption.run.result.root()
    ));
    assert_ne!(
        first_corruption.mutation.before,
        first_corruption.mutation.after
    );
    assert!(f64::from_bits(first_corruption.mutation.after).is_finite());
    assert_ne!(canonical.result.root(), first_corruption.run.result.root());
    assert_ne!(replay.result.root(), replay_corruption.run.result.root());
    assert_eq!(
        first_public_mismatch(&first_corruption.run, &replay_corruption.run),
        None,
        "the corruption seed must independently reproduce the complete red state"
    );
    assert_eq!(
        first_corruption.run.result.canonical_bytes(),
        replay_corruption.run.result.canonical_bytes()
    );

    let first_detail = corruption_detail(canonical, &first_corruption);
    let replay_detail = corruption_detail(replay, &replay_corruption);
    assert_eq!(first_detail, replay_detail);
    let first_event = failure_event(&first_detail, first_corruption.mutation.seed);
    let replay_event = failure_event(&replay_detail, replay_corruption.mutation.seed);
    for event in [&first_event, &replay_event] {
        fs_obs::lint_failure_record(event)
            .expect("disclosed BIPOP corruption must retain its replay seed and detail");
        fs_obs::validate_line(&event.to_jsonl())
            .expect("disclosed BIPOP corruption must remain wire-valid");
    }
    assert_eq!(
        first_event, replay_event,
        "independent seeded red evidence construction replays"
    );
    assert_eq!(
        first_event.content_identity().canonical_bytes(),
        replay_event.content_identity().canonical_bytes()
    );
    let retained = first_event.content_identity_receipt();
    first_event
        .admit_content_identity(&retained)
        .expect("red evidence identity must admit exactly");
    println!("{}", first_event.to_jsonl());

    let panic = catch_unwind(|| {
        assert_mergeable(&first_corruption.run, &canonical.result, &first_event);
    })
    .expect_err("the merge gate must reject the disclosed returned-bit corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{CORRUPTION_SEED:016x}")));
    assert!(message.contains("best.x[1]"));
    assert!(message.contains("ReferenceIdentityMismatch"));
}

#[test]
fn bipop_full_study_replays_and_seeded_failure_is_refused() {
    let first = run_study(INPUT_SEED);
    let replay = run_study(INPUT_SEED);
    let first_accounting = accounting_mismatch(&first);
    let replay_accounting = accounting_mismatch(&replay);
    assert_eq!(first_accounting, None, "original accounting failed");
    assert_eq!(replay_accounting, None, "replay accounting failed");
    assert_eq!(validate_payload(&first), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(admit_against(&first, &first.result), Ok(()));
    assert_eq!(admit_against(&replay, &first.result), Ok(()));

    let mismatch = first_public_mismatch(&first, &replay);
    assert_eq!(
        mismatch, None,
        "same-seed study must replay every public bit"
    );
    assert_eq!(first.fixture.root(), replay.fixture.root());
    assert_eq!(first.result.root(), replay.result.root());
    assert_eq!(
        first.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "the retained complete result frame must replay byte-for-byte"
    );

    emit_green_receipt(&first);
    let green_verdict = emit_green_verdict(&first);
    assert_mergeable(&first, &first.result, &green_verdict);
    exercise_disclosed_corruption(&first, &replay);
}
