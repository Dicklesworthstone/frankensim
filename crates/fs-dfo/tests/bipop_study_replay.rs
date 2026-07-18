//! G5 study-scale replay for the production BIPOP-CMA path (7tv.21.22).
//!
//! The fixture captures every objective call and binds that trace together
//! with the complete ordered public restart ledger, its named best restart,
//! and every legacy `BipopReport`/nested `CmaReport` projection. A test-local
//! algebraic oracle independently checks objective semantics. A separate
//! trace-driven CMA shadow reconstructs every restart boundary, stream sample,
//! budget decision, terminal reason, and report field before linking the first
//! stable global minimum to the exact earliest winning restart. A disclosed
//! seeded mutation changes one returned coordinate bit. The unsealed edit is
//! refused as a stale payload, the self-consistently resealed edit is refused
//! both against the retained reference and by semantic admission, and the
//! resulting red fs-obs evidence is independently reproducible. This is one
//! finite source-snapshot-bound, same-target deterministic study, not an optimizer-quality,
//! refreshed cross-ISA, or performance claim.

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::{
    BIPOP_ADMISSION_SCHEMA_VERSION, BIPOP_RESTART_SCHEMA_VERSION, BipopAdmission, BipopLane,
    BipopReport, BipopRestartRecord, CmaReport, CmaStopReason, admit_bipop, bipop_cmaes,
};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::{StreamCheckpoint, StreamKey};
use std::fmt::Write as _;
use std::panic::catch_unwind;

#[path = "support/bipop_study_replay_goldens.rs"]
mod goldens;

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
const PER_RESTART_GENERATIONS: usize = 250;
const LARGE_RUN_CAP_TRIGGER: u32 = 8;
const CMA_EIGEN_INTERVAL: usize = 1;
const CMA_SPREAD_RELATIVE_LIMIT: f64 = 1e-12;
const CMA_MEANINGFUL_IMPROVEMENT_RELATIVE: f64 = 1e-12;
const CMA_STAGNATION_GENERATIONS: usize = 120;
const OBJECTIVE_ORACLE_ROUNDOFF_SCALE: f64 = 64.0;
const SEMANTIC_ORACLE_VERSION: &str =
    "bipop-algebraic-objective-trace-driven-cma-ledger-accounting-v4";
const STRONG_IDENTITY_DOMAIN: &str = "frankensim.fs-dfo.bipop-study.replay-identity.v1";
const STRONG_EVENT_IDENTITY_DOMAIN: &str =
    "frankensim.fs-dfo.bipop-study.event-content-identity.v1";
const SOURCE_FILE_IDENTITY_DOMAIN: &str = "frankensim.fs-dfo.bipop-study.source-file-identity.v1";
const SOURCE_SNAPSHOT_FILES: &[(&str, &[u8])] = &[
    ("crates/fs-dfo/Cargo.toml", include_bytes!("../Cargo.toml")),
    (
        "crates/fs-dfo/tests/bipop_study_replay.rs",
        include_bytes!("bipop_study_replay.rs"),
    ),
    ("crates/fs-dfo/src/lib.rs", include_bytes!("../src/lib.rs")),
    ("crates/fs-dfo/src/cma.rs", include_bytes!("../src/cma.rs")),
    (
        "crates/fs-la/src/lib.rs",
        include_bytes!("../../fs-la/src/lib.rs"),
    ),
    (
        "crates/fs-la/src/eigen.rs",
        include_bytes!("../../fs-la/src/eigen.rs"),
    ),
    (
        "crates/fs-rand/src/lib.rs",
        include_bytes!("../../fs-rand/src/lib.rs"),
    ),
    (
        "crates/fs-rand/src/philox.rs",
        include_bytes!("../../fs-rand/src/philox.rs"),
    ),
    (
        "crates/fs-math/src/lib.rs",
        include_bytes!("../../fs-math/src/lib.rs"),
    ),
    (
        "crates/fs-math/src/det.rs",
        include_bytes!("../../fs-math/src/det.rs"),
    ),
    (
        "crates/fs-math/src/dd.rs",
        include_bytes!("../../fs-math/src/dd.rs"),
    ),
    (
        "crates/fs-math/src/eft.rs",
        include_bytes!("../../fs-math/src/eft.rs"),
    ),
    (
        "crates/fs-math/src/payne.rs",
        include_bytes!("../../fs-math/src/payne.rs"),
    ),
    (
        "crates/fs-obs/src/ident.rs",
        include_bytes!("../../fs-obs/src/ident.rs"),
    ),
    (
        "crates/fs-obs/src/lib.rs",
        include_bytes!("../../fs-obs/src/lib.rs"),
    ),
    (
        "crates/fs-blake3/src/lib.rs",
        include_bytes!("../../fs-blake3/src/lib.rs"),
    ),
];

// These are the logical stream coordinates and restart rule used by
// `fs_dfo::cma`; recording them makes the private implementation choice
// explicit in the fixture identity.  A change also changes the captured trace.
const CMA_STREAM_KERNEL: u32 = 0xD1F0;
const CMA_SAMPLE_TILE: u32 = 0;
const CMA_RESTART_TILE: u32 = 1;
const RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

#[derive(Debug, Clone)]
struct Evaluation {
    x: Vec<f64>,
    value: f64,
}

#[derive(Debug, Clone)]
struct StudyRun {
    input_seed: u64,
    admission: BipopAdmission,
    fixture: ReplayIdentity,
    report: BipopReport,
    evaluations: Vec<Evaluation>,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    AdmissionReceiptMismatch { field: &'static str },
    FixtureIdentityMismatch { declared: u64, computed: u64 },
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
}

#[derive(Debug, Clone, Copy)]
struct RestartSlice {
    start: usize,
    end: usize,
    lambda: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StreamWitness {
    checkpoint: StreamCheckpoint,
    exact_blocks: u128,
}

#[derive(Debug, Clone)]
struct ReconstructedRestart {
    slice: RestartSlice,
    report: CmaReport,
    stop_reason: CmaStopReason,
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

fn strong_identity_hash(identity: &ReplayIdentity) -> ContentHash {
    hash_domain(STRONG_IDENTITY_DOMAIN, identity.canonical_bytes())
}

fn strong_event_identity_hash(event: &Event) -> ContentHash {
    hash_domain(
        STRONG_EVENT_IDENTITY_DOMAIN,
        event.content_identity().canonical_bytes(),
    )
}

fn source_file_identity_hash(path: &str, bytes: &[u8]) -> ContentHash {
    let domain = format!("{SOURCE_FILE_IDENTITY_DOMAIN}:{path}");
    hash_domain(&domain, bytes)
}

fn lane_name(lane: BipopLane) -> &'static str {
    match lane {
        BipopLane::Large => "large",
        BipopLane::Small => "small",
    }
}

fn stop_reason_name(reason: CmaStopReason) -> &'static str {
    match reason {
        CmaStopReason::TargetReached => "target-reached",
        CmaStopReason::BudgetExhausted => "budget-exhausted",
        CmaStopReason::Stagnated => "stagnated",
    }
}

fn shifted_rastrigin_callback(x: &[f64]) -> f64 {
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

fn shifted_rastrigin_oracle(x: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        DIMENSION,
        "oracle dimension is part of the fixture contract"
    );
    x.iter()
        .zip(SHIFT)
        .map(|(&value, shift)| {
            let z = value - shift;
            let periodic_penalty = 10.0 * (1.0 - fs_math::det::cos(std::f64::consts::TAU * z));
            z * z + periodic_penalty
        })
        .sum()
}

fn objective_oracle_tolerance(recorded: f64, oracle: f64) -> f64 {
    let scale = recorded.abs().max(oracle.abs()).max(1.0);
    OBJECTIVE_ORACLE_ROUNDOFF_SCALE * f64::EPSILON * (DIMENSION as f64) * scale
}

fn expected_base_lambda() -> usize {
    4 + (3.0 * fs_math::det::ln(DIMENSION as f64)).floor() as usize
}

fn admission_mismatch(left: &BipopAdmission, right: &BipopAdmission) -> Option<&'static str> {
    if left.schema_version() != right.schema_version() {
        return Some("schema-version");
    }
    if left.stream_semantics_version() != right.stream_semantics_version() {
        return Some("stream-semantics-version");
    }
    match (left.jacobi_admission(), right.jacobi_admission()) {
        (Some(left_jacobi), Some(right_jacobi)) => {
            if left_jacobi.schema_version() != right_jacobi.schema_version() {
                return Some("jacobi.schema-version");
            }
            if left_jacobi.dimension() != right_jacobi.dimension() {
                return Some("jacobi.dimension");
            }
            if left_jacobi.matrix_entries() != right_jacobi.matrix_entries() {
                return Some("jacobi.matrix-entries");
            }
            if left_jacobi.aggregate_work_elements() != right_jacobi.aggregate_work_elements() {
                return Some("jacobi.aggregate-work-elements");
            }
            if left_jacobi.work_element_cap() != right_jacobi.work_element_cap() {
                return Some("jacobi.work-element-cap");
            }
        }
        (None, None) => {}
        (Some(_), None) | (None, Some(_)) => return Some("jacobi.presence"),
    }
    if left.dimension() != right.dimension() {
        return Some("dimension");
    }
    if left.total_budget() != right.total_budget() {
        return Some("total-budget");
    }
    if left.base_lambda() != right.base_lambda() {
        return Some("base-lambda");
    }
    if left.max_large_lambda() != right.max_large_lambda() {
        return Some("max-large-lambda");
    }
    if left.max_local_budget() != right.max_local_budget() {
        return Some("max-local-budget");
    }
    if left.max_restart_ordinal() != right.max_restart_ordinal() {
        return Some("max-restart-ordinal");
    }
    if left.max_matrix_entries() != right.max_matrix_entries() {
        return Some("max-matrix-entries");
    }
    if left.max_population_entries() != right.max_population_entries() {
        return Some("max-population-entries");
    }
    if left.max_restart_stream_blocks() != right.max_restart_stream_blocks() {
        return Some("max-restart-stream-blocks");
    }
    if left.max_cma_stream_blocks() != right.max_cma_stream_blocks() {
        return Some("max-cma-stream-blocks");
    }
    (left != right).then_some("unclassified-private-field")
}

fn fixture_admission_kat_mismatch(admission: &BipopAdmission) -> Option<String> {
    let checks = [
        ("schema-version", u128::from(admission.schema_version()), 3),
        (
            "stream-semantics-version",
            u128::from(admission.stream_semantics_version()),
            1,
        ),
        ("dimension", admission.dimension() as u128, 4),
        ("total-budget", admission.total_budget() as u128, 6_000),
        ("base-lambda", admission.base_lambda() as u128, 8),
        (
            "max-large-lambda",
            admission.max_large_lambda() as u128,
            2_048,
        ),
        (
            "max-local-budget",
            admission.max_local_budget() as u128,
            512_000,
        ),
        (
            "max-restart-ordinal",
            u128::from(admission.max_restart_ordinal()),
            5_999,
        ),
        (
            "max-matrix-entries",
            admission.max_matrix_entries() as u128,
            16,
        ),
        (
            "max-population-entries",
            admission.max_population_entries() as u128,
            8_192,
        ),
        (
            "max-restart-stream-blocks",
            admission.max_restart_stream_blocks(),
            47_992,
        ),
        (
            "max-cma-stream-blocks",
            admission.max_cma_stream_blocks(),
            47_872,
        ),
    ];
    for (field, found, expected) in checks {
        if found != expected {
            return Some(format!("{field}:{found}!={expected}"));
        }
    }
    let Some(jacobi) = admission.jacobi_admission() else {
        return Some("jacobi:missing".to_string());
    };
    let jacobi_checks = [
        ("schema-version", u128::from(jacobi.schema_version()), 1),
        ("dimension", jacobi.dimension() as u128, 4),
        ("matrix-entries", jacobi.matrix_entries() as u128, 16),
        (
            "aggregate-work-elements",
            jacobi.aggregate_work_elements() as u128,
            76,
        ),
        (
            "work-element-cap",
            jacobi.work_element_cap() as u128,
            67_108_864,
        ),
    ];
    for (field, found, expected) in jacobi_checks {
        if found != expected {
            return Some(format!("jacobi.{field}:{found}!={expected}"));
        }
    }
    None
}

fn fixture_identity(seed: u64, admission: &BipopAdmission) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-bipop-study-fixture-v4")
        .str("suite", SUITE)
        .str("case", CASE)
        .str("red-case", RED_CASE)
        .str("semantic-oracle-version", SEMANTIC_ORACLE_VERSION)
        .str("strong-identity-domain", STRONG_IDENTITY_DOMAIN)
        .str("strong-event-identity-domain", STRONG_EVENT_IDENTITY_DOMAIN)
        .str("source-file-identity-domain", SOURCE_FILE_IDENTITY_DOMAIN)
        .str("algorithm", "fs_dfo::bipop_cmaes")
        .str("objective", "shifted-rastrigin")
        .str(
            "objective-callback-form",
            "sum(10+fma(z,z,-10*cos(tau*z)));z=x-shift",
        )
        .str(
            "objective-oracle-form",
            "sum(z*z+10*(1-cos(tau*z)));z=x-shift",
        )
        .str(
            "objective-oracle-independence",
            "algebraic-regrouping-only;shared-fs-math-cosine",
        )
        .str(
            "objective-oracle-tolerance-rule",
            "roundoff-scale*epsilon*dimension*max(abs(recorded),abs(oracle),1)",
        )
        .str("units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .f64_bits("sigma0", SIGMA0)
        .u64("total-evaluation-budget", usize_u64(TOTAL_BUDGET))
        .f64_bits("f-target", F_TARGET)
        .u64("input-seed", seed)
        .u64("corruption-seed", CORRUPTION_SEED)
        .u64("base-lambda", usize_u64(expected_base_lambda()))
        .str("base-lambda-rule", "4+floor(3*ln(dimension))")
        .u64(
            "per-restart-generations",
            usize_u64(PER_RESTART_GENERATIONS),
        )
        .str(
            "per-restart-budget-rule",
            "min(lambda*per-restart-generations,remaining)",
        )
        .str("large-restart-rule", "large-budget-used<=small-budget-used")
        .str("large-population-rule", "base-lambda*2^large-runs")
        .str(
            "cma-generation-budget-rule",
            "one-initial-evaluation-then-only-complete-lambda-sized-generations",
        )
        .str(
            "cma-candidate-ranking-rule",
            "objective-total-cmp-then-lowest-candidate-index",
        )
        .str(
            "cma-best-update-rule",
            "strictly-lower-objective-only;first-stable-tie",
        )
        .str(
            "bipop-best-update-rule",
            "f64-total-cmp-strictly-less;earliest-restart-wins-ties",
        )
        .u64(
            "bipop-restart-record-schema-version",
            u64::from(BIPOP_RESTART_SCHEMA_VERSION),
        )
        .u64(
            "bipop-admission-schema-version",
            u64::from(BIPOP_ADMISSION_SCHEMA_VERSION),
        )
        .u64(
            "admission-receipt-schema-version",
            u64::from(admission.schema_version()),
        )
        .u64(
            "admission-stream-semantics-version",
            u64::from(admission.stream_semantics_version()),
        )
        .u64(
            "supported-fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64(
            "supported-fs-rand-stream-checkpoint-version",
            u64::from(fs_rand::STREAM_CHECKPOINT_VERSION),
        )
        .str(
            "fs-rand-stream-checkpoint-identity-domain",
            fs_rand::STREAM_CHECKPOINT_IDENTITY_DOMAIN,
        )
        .u64("admission-dimension", usize_u64(admission.dimension()))
        .u64(
            "admission-total-budget",
            usize_u64(admission.total_budget()),
        )
        .u64("admission-base-lambda", usize_u64(admission.base_lambda()))
        .u64(
            "admission-max-large-lambda",
            usize_u64(admission.max_large_lambda()),
        )
        .u64(
            "admission-max-local-budget",
            usize_u64(admission.max_local_budget()),
        )
        .u64(
            "admission-max-restart-ordinal",
            admission.max_restart_ordinal(),
        )
        .u64(
            "admission-max-matrix-entries",
            usize_u64(admission.max_matrix_entries()),
        )
        .u64(
            "admission-max-population-entries",
            usize_u64(admission.max_population_entries()),
        )
        .u64(
            "admission-max-restart-stream-blocks-low",
            admission.max_restart_stream_blocks() as u64,
        )
        .u64(
            "admission-max-restart-stream-blocks-high",
            (admission.max_restart_stream_blocks() >> 64) as u64,
        )
        .u64(
            "admission-max-cma-stream-blocks-low",
            admission.max_cma_stream_blocks() as u64,
        )
        .u64(
            "admission-max-cma-stream-blocks-high",
            (admission.max_cma_stream_blocks() >> 64) as u64,
        )
        .flag(
            "jacobi-admission-present",
            admission.jacobi_admission().is_some(),
        )
        .u64(
            "supported-jacobi-admission-schema-version",
            u64::from(fs_la::eigen::JACOBI_EIGH_ADMISSION_SCHEMA_VERSION),
        );
    for &(path, bytes) in SOURCE_SNAPSHOT_FILES {
        builder = builder.str("source-file-path", path).bytes(
            "source-file-blake3",
            source_file_identity_hash(path, bytes).as_bytes(),
        );
    }
    if let Some(jacobi) = admission.jacobi_admission() {
        builder = builder
            .u64(
                "jacobi-admission-schema-version",
                u64::from(jacobi.schema_version()),
            )
            .u64("jacobi-admission-dimension", usize_u64(jacobi.dimension()))
            .u64(
                "jacobi-admission-matrix-entries",
                usize_u64(jacobi.matrix_entries()),
            )
            .u64(
                "jacobi-admission-aggregate-work-elements",
                usize_u64(jacobi.aggregate_work_elements()),
            )
            .u64(
                "jacobi-admission-work-element-cap",
                usize_u64(jacobi.work_element_cap()),
            );
    }
    builder = builder
        .str(
            "bipop-restart-record-fields",
            "schema-version;ordinal;lane;lambda;allocated-budget;seed;start;half-open-trace-interval;terminal-reason;complete-cma-report",
        )
        .str(
            "bipop-legacy-projections",
            "schedule=ordered-record-lambdas;total-evals=terminal-trace-offset;best=named-best-record-report",
        )
        .u64("large-run-cap-trigger", u64::from(LARGE_RUN_CAP_TRIGGER))
        .u64("cma-eigen-interval", usize_u64(CMA_EIGEN_INTERVAL))
        .str(
            "stagnation-rule",
            "spread<relative-limit*sigma0 OR generations-since-meaningful-improvement>limit",
        )
        .f64_bits("cma-spread-relative-limit", CMA_SPREAD_RELATIVE_LIMIT)
        .f64_bits(
            "cma-meaningful-improvement-relative",
            CMA_MEANINGFUL_IMPROVEMENT_RELATIVE,
        )
        .str(
            "cma-meaningful-improvement-rule",
            "previous-best-minus-candidate>relative*(1+abs(previous-best))",
        )
        .u64(
            "cma-stagnation-generations",
            usize_u64(CMA_STAGNATION_GENERATIONS),
        )
        .f64_bits(
            "objective-oracle-roundoff-scale",
            OBJECTIVE_ORACLE_ROUNDOFF_SCALE,
        )
        .u64("cma-stream-kernel", u64::from(CMA_STREAM_KERNEL))
        .u64("sample-stream-tile", u64::from(CMA_SAMPLE_TILE))
        .u64("restart-stream-tile", u64::from(CMA_RESTART_TILE))
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
        .str(
            "determinism-scope",
            "same-source-snapshot-same-target-test-binary",
        )
        .str("compiler-fingerprint", "not-bound-no-claim")
        .str("target-architecture", std::env::consts::ARCH)
        .str("target-operating-system", std::env::consts::OS)
        .u64("target-pointer-width-bits", u64::from(usize::BITS))
        .str(
            "target-endianness",
            if cfg!(target_endian = "little") {
                "little"
            } else {
                "big"
            },
        )
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
    let mut builder = IdentityBuilder::new("fs-dfo-bipop-study-result-v4")
        .child("fixture", fixture)
        // IdentityBuilder::child retains the legacy 64-bit child root. Carry
        // the complete child preimage as well so the domain-separated result
        // BLAKE3 is transitively bound to every fixture byte rather than only
        // to that compact compatibility projection.
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .bytes(
            "fixture-blake3",
            strong_identity_hash(fixture).as_bytes(),
        )
        .u64("total-evals", usize_u64(report.total_evals))
        .u64("total-budget", usize_u64(report.total_budget()))
        .u64("schedule-length", usize_u64(report.schedule.len()))
        .u64("restart-record-count", usize_u64(report.records().len()))
        .u64("best-restart", usize_u64(report.best_restart()))
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
    for (record_index, record) in report.records().iter().enumerate() {
        builder = builder
            .u64("record-index", usize_u64(record_index))
            .u64("record-schema-version", u64::from(record.schema_version()))
            .u64("record-ordinal", record.ordinal())
            .str("record-lane", lane_name(record.lane()))
            .u64("record-lambda", usize_u64(record.lambda()))
            .u64(
                "record-allocated-budget",
                usize_u64(record.allocated_budget()),
            )
            .u64("record-seed", record.seed())
            .u64("record-start-length", usize_u64(record.start().len()))
            .u64("record-trace-start", usize_u64(record.trace_start()))
            .u64("record-trace-end", usize_u64(record.trace_end()))
            .str(
                "record-terminal-reason",
                stop_reason_name(record.stop_reason()),
            )
            .u64(
                "record-report-x-length",
                usize_u64(record.report().x_best.len()),
            )
            .f64_bits("record-report-best-f", record.report().f_best)
            .u64("record-report-evals", usize_u64(record.report().evals))
            .u64(
                "record-report-generations",
                usize_u64(record.report().generations),
            )
            .flag("record-report-converged", record.report().converged)
            .f64_bits("record-report-sigma", record.report().sigma);
        for (coordinate, &value) in record.start().iter().enumerate() {
            builder = builder
                .u64("record-start-coordinate-index", usize_u64(coordinate))
                .f64_bits("record-start-coordinate", value);
        }
        for (coordinate, &value) in record.report().x_best.iter().enumerate() {
            builder = builder
                .u64("record-report-coordinate-index", usize_u64(coordinate))
                .f64_bits("record-report-coordinate", value);
        }
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

fn execute_study_payload(seed: u64, total_budget: usize) -> (BipopReport, Vec<Evaluation>) {
    let mut evaluations = Vec::with_capacity(total_budget);
    let report = {
        let mut objective = |x: &[f64]| {
            let value = shifted_rastrigin_callback(x);
            evaluations.push(Evaluation {
                x: x.to_vec(),
                value,
            });
            value
        };
        bipop_cmaes(&mut objective, &X0, SIGMA0, total_budget, F_TARGET, seed)
    };
    (report, evaluations)
}

fn run_study(seed: u64) -> StudyRun {
    let admission = admit_bipop(&X0, SIGMA0, TOTAL_BUDGET, Some(F_TARGET), seed)
        .expect("the retained study fixture admits before callbacks");
    let (report, evaluations) = execute_study_payload(seed, TOTAL_BUDGET);
    let fixture = fixture_identity(seed, &admission);
    let result = result_identity(&fixture, &report, &evaluations);
    StudyRun {
        input_seed: seed,
        admission,
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

fn expected_restart_starts(
    input_seed: u64,
    restart_count: usize,
) -> (Vec<Vec<f64>>, StreamCheckpoint) {
    let key = StreamKey {
        seed: input_seed,
        kernel: CMA_STREAM_KERNEL,
        tile: CMA_RESTART_TILE,
    };
    let mut stream = key.stream();
    let starts = (0..restart_count)
        .map(|restart| {
            if restart == 0 {
                X0.to_vec()
            } else {
                X0.iter()
                    .map(|&value| SIGMA0.mul_add(stream.next_normal(), value))
                    .collect()
            }
        })
        .collect();
    (starts, stream.checkpoint())
}

fn exact_cma_stream_witness(
    checkpoint: StreamCheckpoint,
    expected_key: StreamKey,
    admitted_stream_semantics_version: u32,
    dimension: usize,
    lambda: usize,
    generations: usize,
) -> Result<StreamWitness, String> {
    if checkpoint.checkpoint_version != fs_rand::STREAM_CHECKPOINT_VERSION {
        return Err(format!(
            "CMA-shadow-checkpoint-version-{}!=supported-{}",
            checkpoint.checkpoint_version,
            fs_rand::STREAM_CHECKPOINT_VERSION,
        ));
    }
    if checkpoint.stream_semantics_version != admitted_stream_semantics_version {
        return Err(format!(
            "CMA-shadow-stream-semantics-{}!=admitted-{admitted_stream_semantics_version}",
            checkpoint.stream_semantics_version
        ));
    }
    if checkpoint.key != expected_key {
        return Err("CMA-shadow-stream-key-mismatch".to_string());
    }
    let expected = 2_u128
        .checked_mul(dimension as u128)
        .and_then(|blocks| blocks.checked_mul(lambda as u128))
        .and_then(|blocks| blocks.checked_mul(generations as u128))
        .ok_or_else(|| "CMA-shadow-counter-formula-overflow".to_string())?;
    if u128::from(checkpoint.index) != expected {
        return Err(format!(
            "CMA-shadow-counter-index-{}!=formula-{expected}",
            checkpoint.index
        ));
    }
    Ok(StreamWitness {
        checkpoint,
        exact_blocks: expected,
    })
}

fn report_mismatch(prefix: &str, actual: &CmaReport, expected: &CmaReport) -> Option<String> {
    if actual.x_best.len() != expected.x_best.len() {
        return Some(format!(
            "{prefix}.x-length:{}!={}",
            actual.x_best.len(),
            expected.x_best.len()
        ));
    }
    for (coordinate, (&actual, &expected)) in actual.x_best.iter().zip(&expected.x_best).enumerate()
    {
        if actual.to_bits() != expected.to_bits() {
            return Some(format!(
                "{prefix}.x[{coordinate}]:0x{:016x}!=0x{:016x}",
                actual.to_bits(),
                expected.to_bits()
            ));
        }
    }
    if actual.f_best.to_bits() != expected.f_best.to_bits() {
        return Some(format!(
            "{prefix}.f:0x{:016x}!=0x{:016x}",
            actual.f_best.to_bits(),
            expected.f_best.to_bits()
        ));
    }
    if actual.evals != expected.evals {
        return Some(format!(
            "{prefix}.evals:{}!={}",
            actual.evals, expected.evals
        ));
    }
    if actual.generations != expected.generations {
        return Some(format!(
            "{prefix}.generations:{}!={}",
            actual.generations, expected.generations
        ));
    }
    if actual.converged != expected.converged {
        return Some(format!(
            "{prefix}.converged:{}!={}",
            actual.converged, expected.converged
        ));
    }
    if actual.sigma.to_bits() != expected.sigma.to_bits() {
        return Some(format!(
            "{prefix}.sigma:0x{:016x}!=0x{:016x}",
            actual.sigma.to_bits(),
            expected.sigma.to_bits()
        ));
    }
    None
}

/// Reconstruct one CMA run from its declared root inputs and callback trace.
///
/// This deliberately does not call `cmaes`, `bipop_cmaes`, or the production
/// ledger validator. It shares the declared deterministic math/eigensolver
/// primitives, but independently advances the public algorithm equations and
/// requires every generated candidate bit to occur at the exact trace slot.
/// The fixture therefore claims a trace-driven shadow, not implementation-total
/// independence from the production numerical substrate.
#[allow(clippy::too_many_lines)]
fn reconstruct_cma_from_trace(
    evaluations: &[Evaluation],
    start: &[f64],
    lambda: usize,
    allocated_budget: usize,
    seed: u64,
    jacobi_admission: fs_la::eigen::JacobiEighAdmission,
    admitted_stream_semantics_version: u32,
) -> Result<(CmaReport, CmaStopReason, StreamWitness), String> {
    let Some(initial) = evaluations.first() else {
        return Err("empty-callback-trace".to_string());
    };
    if !same_point_bits(&initial.x, start) {
        return Err("initial-callback-point-does-not-match-restart-start".to_string());
    }
    if start.len() != DIMENSION {
        return Err(format!(
            "restart-dimension:{}!=expected-{DIMENSION}",
            start.len()
        ));
    }
    if allocated_budget == 0 {
        return Err("zero-local-budget-cannot-admit-initial-callback".to_string());
    }

    let n = start.len();
    let recomputed_jacobi = fs_la::eigen::admit_jacobi_eigh(n)
        .map_err(|error| format!("shadow-jacobi-admission-refused:{error}"))?;
    if recomputed_jacobi != jacobi_admission {
        return Err("shadow-jacobi-admission-receipt-mismatch".to_string());
    }
    let lambda = lambda.max(4);
    let mu = lambda / 2;
    let raw: Vec<f64> = (0..mu)
        .map(|index| {
            fs_math::det::ln(f64::midpoint(lambda as f64, 1.0))
                - fs_math::det::ln(index as f64 + 1.0)
        })
        .collect();
    let weight_sum: f64 = raw.iter().sum();
    let weights: Vec<f64> = raw.iter().map(|weight| weight / weight_sum).collect();
    let mu_eff = 1.0 / weights.iter().map(|weight| weight * weight).sum::<f64>();
    let n_f64 = n as f64;
    let c_c = (4.0 + mu_eff / n_f64) / (n_f64 + 4.0 + 2.0 * mu_eff / n_f64);
    let c_s = (mu_eff + 2.0) / (n_f64 + mu_eff + 5.0);
    let c_1 = 2.0 / ((n_f64 + 1.3) * (n_f64 + 1.3) + mu_eff);
    let c_mu = (1.0 - c_1)
        .min(2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((n_f64 + 2.0) * (n_f64 + 2.0) + mu_eff));
    let damping =
        1.0 + 2.0 * (fs_math::det::sqrt((mu_eff - 1.0) / (n_f64 + 1.0)) - 1.0).max(0.0) + c_s;
    let expected_normal_norm =
        fs_math::det::sqrt(n_f64) * (1.0 - 1.0 / (4.0 * n_f64) + 1.0 / (21.0 * n_f64 * n_f64));

    let mut mean = start.to_vec();
    let mut sigma = SIGMA0;
    let mut covariance = vec![0.0_f64; n * n];
    for coordinate in 0..n {
        covariance[coordinate * n + coordinate] = 1.0;
    }
    let mut covariance_path = vec![0.0_f64; n];
    let mut sigma_path = vec![0.0_f64; n];
    let mut eigenvectors = covariance.clone();
    let mut eigenvalue_roots = vec![1.0_f64; n];
    let stream_key = StreamKey {
        seed,
        kernel: CMA_STREAM_KERNEL,
        tile: CMA_SAMPLE_TILE,
    };
    let mut stream = stream_key.stream();

    let mut x_best = mean.clone();
    let mut f_best = initial.value;
    let mut evals = 1usize;
    let mut generations = 0usize;
    if f_best <= F_TARGET {
        if evaluations.len() != 1 {
            return Err(format!(
                "target-reached-at-initial-callback-but-trace-has-{}-entries",
                evaluations.len()
            ));
        }
        let stream_witness = exact_cma_stream_witness(
            stream.checkpoint(),
            stream_key,
            admitted_stream_semantics_version,
            n,
            lambda,
            generations,
        )?;
        return Ok((
            CmaReport {
                x_best,
                f_best,
                evals,
                generations,
                converged: true,
                sigma,
            },
            CmaStopReason::TargetReached,
            stream_witness,
        ));
    }
    let mut generations_since_improvement = 0usize;
    let mut normal_samples = vec![vec![0.0_f64; n]; lambda];
    let mut transformed_samples = vec![vec![0.0_f64; n]; lambda];
    let mut fitness = vec![0.0_f64; lambda];

    while evals + lambda <= allocated_budget {
        generations += 1;
        if generations % CMA_EIGEN_INTERVAL.max(1) == 1 || CMA_EIGEN_INTERVAL <= 1 {
            for row in 0..n {
                for column in row + 1..n {
                    let average =
                        f64::midpoint(covariance[row * n + column], covariance[column * n + row]);
                    covariance[row * n + column] = average;
                    covariance[column * n + row] = average;
                }
            }
            let (eigenvalues, next_eigenvectors) = fs_la::eigen::jacobi_eigh(&covariance, n);
            let maximum = eigenvalues
                .last()
                .copied()
                .unwrap_or(1.0)
                .max(f64::MIN_POSITIVE);
            for (index, &eigenvalue) in eigenvalues.iter().enumerate() {
                eigenvalue_roots[index] = fs_math::det::sqrt(eigenvalue.max(1e-14 * maximum));
            }
            eigenvectors.copy_from_slice(&next_eigenvectors);
        }

        for (candidate, normal) in normal_samples.iter_mut().enumerate() {
            for coordinate in normal.iter_mut() {
                *coordinate = stream.next_normal();
            }
            let transformed = &mut transformed_samples[candidate];
            for row in 0..n {
                let mut accumulator = 0.0_f64;
                for column in 0..n {
                    accumulator = (eigenvectors[row * n + column] * eigenvalue_roots[column])
                        .mul_add(normal[column], accumulator);
                }
                transformed[row] = accumulator;
            }
            let point: Vec<f64> = mean
                .iter()
                .zip(transformed.iter())
                .map(|(location, displacement)| sigma.mul_add(*displacement, *location))
                .collect();
            let Some(observed) = evaluations.get(evals) else {
                return Err(format!(
                    "shadow-required-callback-{evals}-for-generation-{generations}"
                ));
            };
            if !same_point_bits(&point, &observed.x) {
                return Err(format!(
                    "generation-{generations}-candidate-{candidate}-stream-point-mismatch"
                ));
            }
            fitness[candidate] = observed.value;
            evals += 1;
            if fitness[candidate] < f_best {
                if f_best - fitness[candidate]
                    > CMA_MEANINGFUL_IMPROVEMENT_RELATIVE * (1.0 + f_best.abs())
                {
                    generations_since_improvement = 0;
                }
                f_best = fitness[candidate];
                x_best = point;
            }
        }
        generations_since_improvement += 1;
        if f_best <= F_TARGET {
            if evals != evaluations.len() {
                return Err(format!(
                    "target-terminal-offset-{evals}!=trace-length-{}",
                    evaluations.len()
                ));
            }
            let stream_witness = exact_cma_stream_witness(
                stream.checkpoint(),
                stream_key,
                admitted_stream_semantics_version,
                n,
                lambda,
                generations,
            )?;
            return Ok((
                CmaReport {
                    x_best,
                    f_best,
                    evals,
                    generations,
                    converged: true,
                    sigma,
                },
                CmaStopReason::TargetReached,
                stream_witness,
            ));
        }

        let mut order: Vec<usize> = (0..lambda).collect();
        order.sort_by(|&left, &right| {
            fitness[left]
                .total_cmp(&fitness[right])
                .then(left.cmp(&right))
        });
        let mut weighted_step = vec![0.0_f64; n];
        for (weight, &candidate) in weights.iter().zip(&order) {
            for coordinate in 0..n {
                weighted_step[coordinate] = weight.mul_add(
                    transformed_samples[candidate][coordinate],
                    weighted_step[coordinate],
                );
            }
        }
        for coordinate in 0..n {
            mean[coordinate] = sigma.mul_add(weighted_step[coordinate], mean[coordinate]);
        }

        let mut inverse_root_step = vec![0.0_f64; n];
        for row in 0..n {
            let mut accumulator = 0.0_f64;
            for column in 0..n {
                accumulator =
                    eigenvectors[column * n + row].mul_add(weighted_step[column], accumulator);
            }
            inverse_root_step[row] = accumulator / eigenvalue_roots[row];
        }
        let mut whitened_step = vec![0.0_f64; n];
        for row in 0..n {
            let mut accumulator = 0.0_f64;
            for column in 0..n {
                accumulator =
                    eigenvectors[row * n + column].mul_add(inverse_root_step[column], accumulator);
            }
            whitened_step[row] = accumulator;
        }
        let sigma_path_scale = fs_math::det::sqrt(c_s * (2.0 - c_s) * mu_eff);
        for coordinate in 0..n {
            sigma_path[coordinate] = (1.0 - c_s).mul_add(
                sigma_path[coordinate],
                sigma_path_scale * whitened_step[coordinate],
            );
        }
        let sigma_path_norm =
            fs_math::det::sqrt(sigma_path.iter().map(|value| value * value).sum::<f64>());
        sigma *=
            fs_math::det::exp((c_s / damping) * (sigma_path_norm / expected_normal_norm - 1.0));
        let spread = sigma
            * eigenvalue_roots
                .iter()
                .fold(0.0_f64, |maximum, &root| maximum.max(root));
        if spread < CMA_SPREAD_RELATIVE_LIMIT * SIGMA0
            || generations_since_improvement > CMA_STAGNATION_GENERATIONS
        {
            if evals != evaluations.len() {
                return Err(format!(
                    "stagnation-terminal-offset-{evals}!=trace-length-{}",
                    evaluations.len()
                ));
            }
            let stream_witness = exact_cma_stream_witness(
                stream.checkpoint(),
                stream_key,
                admitted_stream_semantics_version,
                n,
                lambda,
                generations,
            )?;
            return Ok((
                CmaReport {
                    x_best,
                    f_best,
                    evals,
                    generations,
                    converged: false,
                    sigma,
                },
                CmaStopReason::Stagnated,
                stream_witness,
            ));
        }

        let covariance_path_active = sigma_path_norm
            / fs_math::det::sqrt(
                1.0 - fs_math::det::powi(
                    1.0 - c_s,
                    2 * i32::try_from(generations.min(100_000))
                        .expect("fixture generation count fits i32"),
                ),
            )
            < (1.4 + 2.0 / (n_f64 + 1.0)) * expected_normal_norm;
        let covariance_path_scale = fs_math::det::sqrt(c_c * (2.0 - c_c) * mu_eff);
        for coordinate in 0..n {
            let innovation = if covariance_path_active {
                covariance_path_scale * weighted_step[coordinate]
            } else {
                0.0
            };
            covariance_path[coordinate] =
                (1.0 - c_c).mul_add(covariance_path[coordinate], innovation);
        }
        let inactive_path_correction = if covariance_path_active {
            0.0
        } else {
            c_c * (2.0 - c_c)
        };
        for row in 0..n {
            for column in 0..n {
                let mut rank_mu = 0.0_f64;
                for (weight, &candidate) in weights.iter().zip(&order) {
                    rank_mu = (weight * transformed_samples[candidate][row])
                        .mul_add(transformed_samples[candidate][column], rank_mu);
                }
                let rank_one = covariance_path[row] * covariance_path[column];
                covariance[row * n + column] = (1.0 - c_1 - c_mu).mul_add(
                    covariance[row * n + column],
                    c_1.mul_add(
                        rank_one + inactive_path_correction * covariance[row * n + column],
                        c_mu * rank_mu,
                    ),
                );
            }
        }
    }

    if evals != evaluations.len() {
        return Err(format!(
            "budget-terminal-offset-{evals}!=trace-length-{}",
            evaluations.len()
        ));
    }
    let stream_witness = exact_cma_stream_witness(
        stream.checkpoint(),
        stream_key,
        admitted_stream_semantics_version,
        n,
        lambda,
        generations,
    )?;
    Ok((
        CmaReport {
            x_best,
            f_best,
            evals,
            generations,
            converged: false,
            sigma,
        },
        CmaStopReason::BudgetExhausted,
        stream_witness,
    ))
}

#[allow(clippy::too_many_lines)] // Every public restart field is checked in schema order.
fn reconstruct_restart_ledger(run: &StudyRun) -> Result<Vec<ReconstructedRestart>, String> {
    let records = run.report.records();
    if records.is_empty() {
        return Err("restart-ledger-is-empty".to_string());
    }
    if run.report.schedule.len() != records.len() {
        return Err(format!(
            "legacy-schedule-length:{}!=record-count:{}",
            run.report.schedule.len(),
            records.len()
        ));
    }
    let (expected_starts, restart_checkpoint) =
        expected_restart_starts(run.input_seed, records.len());
    let expected_restart_key = StreamKey {
        seed: run.input_seed,
        kernel: CMA_STREAM_KERNEL,
        tile: CMA_RESTART_TILE,
    };
    if restart_checkpoint.checkpoint_version != fs_rand::STREAM_CHECKPOINT_VERSION {
        return Err(format!(
            "restart-checkpoint-version-{}!=supported-{}",
            restart_checkpoint.checkpoint_version,
            fs_rand::STREAM_CHECKPOINT_VERSION,
        ));
    }
    if restart_checkpoint.stream_semantics_version != run.admission.stream_semantics_version() {
        return Err(format!(
            "restart-stream-semantics-{}!=admitted-{}",
            restart_checkpoint.stream_semantics_version,
            run.admission.stream_semantics_version()
        ));
    }
    if restart_checkpoint.key != expected_restart_key {
        return Err("restart-stream-key-mismatch".to_string());
    }
    let expected_restart_stream_blocks = 2_u128
        .checked_mul(DIMENSION as u128)
        .and_then(|blocks| blocks.checked_mul((records.len() - 1) as u128))
        .ok_or_else(|| "restart-stream-counter-formula-overflow".to_string())?;
    if u128::from(restart_checkpoint.index) != expected_restart_stream_blocks {
        return Err(format!(
            "restart-stream-index-{}!=formula-{expected_restart_stream_blocks}",
            restart_checkpoint.index
        ));
    }
    if expected_restart_stream_blocks > run.admission.max_restart_stream_blocks() {
        return Err(format!(
            "restart-stream-blocks-{expected_restart_stream_blocks}>admitted-{}",
            run.admission.max_restart_stream_blocks()
        ));
    }
    let jacobi_admission = run
        .admission
        .jacobi_admission()
        .ok_or_else(|| "study-admission-omits-reachable-jacobi-authority".to_string())?;
    if run.admission.dimension() != DIMENSION
        || jacobi_admission.dimension() != DIMENSION
        || jacobi_admission.matrix_entries() != DIMENSION * DIMENSION
    {
        return Err("study-admission-dimension-authority-mismatch".to_string());
    }
    let mut boundaries = Vec::with_capacity(expected_starts.len());
    for (restart, expected) in expected_starts.iter().enumerate() {
        let matches: Vec<usize> = run
            .evaluations
            .iter()
            .enumerate()
            .filter_map(|(index, evaluation)| {
                same_point_bits(&evaluation.x, expected).then_some(index)
            })
            .collect();
        if matches.len() != 1 {
            return Err(format!(
                "restart[{restart}]-start-occurrences:{}!=1",
                matches.len()
            ));
        }
        boundaries.push(matches[0]);
    }
    if boundaries.first().copied() != Some(0) {
        return Err(format!(
            "first-restart-boundary:{:?}!=0",
            boundaries.first()
        ));
    }
    if boundaries.windows(2).any(|window| window[0] >= window[1]) {
        return Err(format!(
            "restart-boundaries-not-strictly-increasing:{boundaries:?}"
        ));
    }

    let mut reconstructed = Vec::with_capacity(boundaries.len());
    let mut large_runs = 0u32;
    let mut small_budget_used = 0usize;
    let mut large_budget_used = 0usize;
    let base_lambda = expected_base_lambda();
    for (restart, (&start, record)) in boundaries.iter().zip(records).enumerate() {
        let end = boundaries
            .get(restart + 1)
            .copied()
            .unwrap_or(run.evaluations.len());
        let run_large = large_budget_used <= small_budget_used;
        let expected_lane = if run_large {
            BipopLane::Large
        } else {
            BipopLane::Small
        };
        let expected_lambda = if run_large {
            let multiplier = 1usize.checked_shl(large_runs).ok_or_else(|| {
                format!("restart[{restart}]-large-population-shift-overflow:{large_runs}")
            })?;
            base_lambda.checked_mul(multiplier).ok_or_else(|| {
                format!("restart[{restart}]-large-population-multiplication-overflow")
            })?
        } else {
            base_lambda
        };
        let expected_ordinal =
            u64::try_from(restart).map_err(|_| format!("restart[{restart}]-ordinal-overflow"))?;
        let seed_delta = expected_ordinal
            .checked_mul(RESTART_SEED_STRIDE)
            .ok_or_else(|| format!("restart[{restart}]-seed-stride-multiplication-overflow"))?;
        let expected_seed = run
            .input_seed
            .checked_add(seed_delta)
            .ok_or_else(|| format!("restart[{restart}]-seed-range-overflow"))?;
        if record.schema_version() != BIPOP_RESTART_SCHEMA_VERSION {
            return Err(format!(
                "restart[{restart}]-schema-version:{}!={BIPOP_RESTART_SCHEMA_VERSION}",
                record.schema_version()
            ));
        }
        if record.ordinal() != expected_ordinal {
            return Err(format!(
                "restart[{restart}]-ordinal:{}!={expected_ordinal}",
                record.ordinal()
            ));
        }
        if record.lane() != expected_lane {
            return Err(format!(
                "restart[{restart}]-lane:{}!={}",
                lane_name(record.lane()),
                lane_name(expected_lane)
            ));
        }
        if record.lambda() != expected_lambda {
            return Err(format!(
                "restart[{restart}]-lambda:{}!={expected_lambda}",
                record.lambda()
            ));
        }
        if run.report.schedule[restart] != expected_lambda {
            return Err(format!(
                "legacy-schedule[{restart}]={}!=reconstructed-{expected_lambda}",
                run.report.schedule[restart]
            ));
        }
        if record.seed() != expected_seed {
            return Err(format!(
                "restart[{restart}]-seed:0x{:016x}!=0x{expected_seed:016x}",
                record.seed()
            ));
        }
        if !same_point_bits(record.start(), &expected_starts[restart]) {
            return Err(format!(
                "restart[{restart}]-start-does-not-match-declared-restart-stream"
            ));
        }
        if start >= end {
            return Err(format!(
                "restart[{restart}]-empty-or-reversed-slice:{start}..{end}"
            ));
        }
        if record.trace_start() != start {
            return Err(format!(
                "restart[{restart}]-trace-start:{}!={start}",
                record.trace_start()
            ));
        }
        if record.trace_end() != end {
            return Err(format!(
                "restart[{restart}]-trace-end:{}!={end}",
                record.trace_end()
            ));
        }
        let run_evals = end - start;
        let remaining = TOTAL_BUDGET.checked_sub(start).ok_or_else(|| {
            format!("restart[{restart}]-start-{start}-exceeds-budget-{TOTAL_BUDGET}")
        })?;
        let nominal_budget = expected_lambda
            .checked_mul(PER_RESTART_GENERATIONS)
            .ok_or_else(|| format!("restart[{restart}]-nominal-budget-overflow"))?;
        let admitted_budget = nominal_budget.min(remaining);
        if record.allocated_budget() != admitted_budget {
            return Err(format!(
                "restart[{restart}]-allocated-budget:{}!={admitted_budget}",
                record.allocated_budget()
            ));
        }
        if run_evals > admitted_budget {
            return Err(format!(
                "restart[{restart}]-evals-{run_evals}>admitted-budget-{admitted_budget}"
            ));
        }
        if (run_evals - 1) % expected_lambda != 0 {
            return Err(format!(
                "restart[{restart}]-evals-{run_evals}-not-1-plus-whole-generations-of-{expected_lambda}"
            ));
        }
        if admitted_budget >= 1 + expected_lambda && run_evals < 1 + expected_lambda {
            return Err(format!(
                "restart[{restart}]-omitted-admissible-first-generation"
            ));
        }
        let slice = RestartSlice {
            start,
            end,
            lambda: expected_lambda,
        };
        let (expected_report, expected_stop_reason, cma_stream_witness) =
            reconstruct_cma_from_trace(
                &run.evaluations[start..end],
                &expected_starts[restart],
                expected_lambda,
                admitted_budget,
                expected_seed,
                jacobi_admission,
                run.admission.stream_semantics_version(),
            )
            .map_err(|mismatch| format!("restart[{restart}]-{mismatch}"))?;
        if cma_stream_witness.checkpoint.stream_semantics_version
            != run.admission.stream_semantics_version()
        {
            return Err(format!(
                "restart[{restart}]-CMA-stream-semantics-witness-mismatch"
            ));
        }
        if cma_stream_witness.exact_blocks > run.admission.max_cma_stream_blocks() {
            return Err(format!(
                "restart[{restart}]-CMA-stream-blocks-{}>admitted-{}",
                cma_stream_witness.exact_blocks,
                run.admission.max_cma_stream_blocks()
            ));
        }
        if let Some(mismatch) = report_mismatch(
            &format!("restart[{restart}].report"),
            record.report(),
            &expected_report,
        ) {
            return Err(mismatch);
        }
        if record.stop_reason() != expected_stop_reason {
            return Err(format!(
                "restart[{restart}]-terminal-reason:{}!={}",
                stop_reason_name(record.stop_reason()),
                stop_reason_name(expected_stop_reason)
            ));
        }
        reconstructed.push(ReconstructedRestart {
            slice,
            report: expected_report,
            stop_reason: expected_stop_reason,
        });

        if run_large {
            large_budget_used = large_budget_used
                .checked_add(run_evals)
                .ok_or_else(|| format!("restart[{restart}]-large-budget-overflow"))?;
            large_runs += 1;
        } else {
            small_budget_used = small_budget_used
                .checked_add(run_evals)
                .ok_or_else(|| format!("restart[{restart}]-small-budget-overflow"))?;
        }
        let has_next = restart + 1 < boundaries.len();
        if has_next
            && (end >= TOTAL_BUDGET
                || large_runs > LARGE_RUN_CAP_TRIGGER
                || record.report().converged)
        {
            return Err(format!(
                "restart[{restart}]-schedule-continued-after-terminal-condition"
            ));
        }
    }
    let terminal_converged = reconstructed
        .last()
        .is_some_and(|restart| restart.stop_reason == CmaStopReason::TargetReached);
    if run.evaluations.len() < TOTAL_BUDGET
        && large_runs <= LARGE_RUN_CAP_TRIGGER
        && !terminal_converged
    {
        return Err(format!(
            "schedule-ended-before-budget-or-large-run-cap:evals={};large-runs={large_runs}",
            run.evaluations.len()
        ));
    }
    Ok(reconstructed)
}

#[allow(clippy::too_many_lines)] // Complete trace and public-report accounting is the oracle.
fn accounting_mismatch(run: &StudyRun) -> Option<String> {
    let expected_admission =
        match admit_bipop(&X0, SIGMA0, TOTAL_BUDGET, Some(F_TARGET), run.input_seed) {
            Ok(admission) => admission,
            Err(error) => return Some(format!("fixture-admission-refused:{error}")),
        };
    if let Some(field) = admission_mismatch(&run.admission, &expected_admission) {
        return Some(format!("admission-receipt.{field}"));
    }
    if let Some(mismatch) = fixture_admission_kat_mismatch(&run.admission) {
        return Some(format!("admission-KAT.{mismatch}"));
    }
    if run.report.total_budget() != run.admission.total_budget() {
        return Some(format!(
            "report-total-budget:{}!=admission-budget:{}",
            run.report.total_budget(),
            run.admission.total_budget()
        ));
    }
    if run.report.total_budget() != TOTAL_BUDGET {
        return Some(format!(
            "retained-total-budget:{}!=fixture-budget:{TOTAL_BUDGET}",
            run.report.total_budget()
        ));
    }
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
    if run.report.records().len() < 2 {
        return Some(format!(
            "restart-ledger-too-short:{};fixture-must-exercise-a-restart",
            run.report.records().len()
        ));
    }
    if run.report.records().len() > run.report.total_evals {
        return Some(format!(
            "restart-record-count:{}>total-evals:{}",
            run.report.records().len(),
            run.report.total_evals
        ));
    }
    let mut first_trace_best: Option<(usize, &Evaluation)> = None;
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
        let oracle = shifted_rastrigin_oracle(&evaluation.x);
        let tolerance = objective_oracle_tolerance(evaluation.value, oracle);
        if !oracle.is_finite() || (evaluation.value - oracle).abs() > tolerance {
            return Some(format!(
                "trace[{evaluation_index}]-objective:recorded=0x{:016x};oracle=0x{:016x};tolerance=0x{:016x}",
                evaluation.value.to_bits(),
                oracle.to_bits(),
                tolerance.to_bits()
            ));
        }
        if first_trace_best
            .is_none_or(|(_, current)| evaluation.value.total_cmp(&current.value).is_lt())
        {
            first_trace_best = Some((evaluation_index, evaluation));
        }
    }
    let reconstructed = match reconstruct_restart_ledger(run) {
        Ok(reconstructed) => reconstructed,
        Err(mismatch) => return Some(mismatch),
    };

    let mut expected_best_restart = 0usize;
    for restart in 1..reconstructed.len() {
        if reconstructed[restart]
            .report
            .f_best
            .total_cmp(&reconstructed[expected_best_restart].report.f_best)
            .is_lt()
        {
            expected_best_restart = restart;
        }
    }
    if run.report.best_restart() != expected_best_restart {
        return Some(format!(
            "best-restart:{}!=earliest-total-cmp-winner-{expected_best_restart}",
            run.report.best_restart()
        ));
    }
    let Some(named_best) = run.report.best_record() else {
        return Some(format!(
            "best-record-missing-at-index-{}",
            run.report.best_restart()
        ));
    };
    if let Some(mismatch) = report_mismatch(
        "named-best-record",
        named_best.report(),
        &reconstructed[expected_best_restart].report,
    ) {
        return Some(mismatch);
    }
    if let Some(mismatch) = report_mismatch(
        "legacy-best",
        &run.report.best,
        &reconstructed[expected_best_restart].report,
    ) {
        return Some(mismatch);
    }

    let best = &run.report.best;
    let (first_trace_best_index, first_trace_best) =
        first_trace_best.expect("positive total-eval accounting makes trace nonempty");
    if first_trace_best.value.to_bits() != best.f_best.to_bits() {
        return Some(format!(
            "complete-trace-minimum=0x{:016x};reported-best=0x{:016x}",
            first_trace_best.value.to_bits(),
            best.f_best.to_bits()
        ));
    }
    let best_oracle = shifted_rastrigin_oracle(&best.x_best);
    let best_tolerance = objective_oracle_tolerance(best.f_best, best_oracle);
    if !best_oracle.is_finite() || (best.f_best - best_oracle).abs() > best_tolerance {
        return Some(format!(
            "best-point-objective:oracle=0x{:016x};reported=0x{:016x};tolerance=0x{:016x}",
            best_oracle.to_bits(),
            best.f_best.to_bits(),
            best_tolerance.to_bits()
        ));
    }
    if !same_point_bits(&first_trace_best.x, &best.x_best) {
        return Some("reported-best-is-not-the-first-stable-trace-minimum".to_string());
    }
    let Some((winning_restart, winning)) = reconstructed.iter().enumerate().find(|(_, restart)| {
        restart.slice.start <= first_trace_best_index && first_trace_best_index < restart.slice.end
    }) else {
        return Some(format!(
            "trace-minimum-index-{first_trace_best_index}-is-outside-restart-slices"
        ));
    };
    if winning_restart != expected_best_restart {
        return Some(format!(
            "trace-minimum-restart-{winning_restart}!=named-best-restart-{expected_best_restart}"
        ));
    }
    let winning_run_evals = winning.slice.end - winning.slice.start;
    let winning_generations = (winning_run_evals - 1) / winning.slice.lambda;
    if best.evals != winning_run_evals || best.generations != winning_generations {
        return Some(format!(
            "best-run-accounting:reported-evals-{}-generations-{}!=restart-{winning_restart}-evals-{winning_run_evals}-generations-{winning_generations}",
            best.evals, best.generations
        ));
    }
    if let Err(error) = run.report.validate_ledger() {
        return Some(format!(
            "production-ledger-validator-refused:restart={:?};invariant={}",
            error.restart(),
            error.invariant()
        ));
    }
    None
}

fn restart_record_mismatch(
    restart: usize,
    left: &BipopRestartRecord,
    right: &BipopRestartRecord,
) -> Option<String> {
    if left.schema_version() != right.schema_version() {
        return Some(format!(
            "record[{restart}].schema-version:{}!={}",
            left.schema_version(),
            right.schema_version()
        ));
    }
    if left.ordinal() != right.ordinal() {
        return Some(format!(
            "record[{restart}].ordinal:{}!={}",
            left.ordinal(),
            right.ordinal()
        ));
    }
    if left.lane() != right.lane() {
        return Some(format!(
            "record[{restart}].lane:{}!={}",
            lane_name(left.lane()),
            lane_name(right.lane())
        ));
    }
    if left.lambda() != right.lambda() {
        return Some(format!(
            "record[{restart}].lambda:{}!={}",
            left.lambda(),
            right.lambda()
        ));
    }
    if left.allocated_budget() != right.allocated_budget() {
        return Some(format!(
            "record[{restart}].allocated-budget:{}!={}",
            left.allocated_budget(),
            right.allocated_budget()
        ));
    }
    if left.seed() != right.seed() {
        return Some(format!(
            "record[{restart}].seed:0x{:016x}!=0x{:016x}",
            left.seed(),
            right.seed()
        ));
    }
    if !same_point_bits(left.start(), right.start()) {
        return Some(format!("record[{restart}].start"));
    }
    if left.trace_start() != right.trace_start() {
        return Some(format!(
            "record[{restart}].trace-start:{}!={}",
            left.trace_start(),
            right.trace_start()
        ));
    }
    if left.trace_end() != right.trace_end() {
        return Some(format!(
            "record[{restart}].trace-end:{}!={}",
            left.trace_end(),
            right.trace_end()
        ));
    }
    if left.stop_reason() != right.stop_reason() {
        return Some(format!(
            "record[{restart}].terminal-reason:{}!={}",
            stop_reason_name(left.stop_reason()),
            stop_reason_name(right.stop_reason())
        ));
    }
    report_mismatch(
        &format!("record[{restart}].report"),
        left.report(),
        right.report(),
    )
}

#[allow(clippy::too_many_lines)] // Exhaustive field-by-field public-state audit.
fn first_public_mismatch(left: &StudyRun, right: &StudyRun) -> Option<String> {
    if left.input_seed != right.input_seed {
        return Some(format!(
            "input-seed:0x{:016x}!=0x{:016x}",
            left.input_seed, right.input_seed
        ));
    }
    if let Some(field) = admission_mismatch(&left.admission, &right.admission) {
        return Some(format!("admission-receipt.{field}"));
    }
    if left.fixture.canonical_bytes() != right.fixture.canonical_bytes() {
        return Some("fixture-identity".to_string());
    }
    if left.report.total_evals != right.report.total_evals {
        return Some(format!(
            "total-evals:{}!={}",
            left.report.total_evals, right.report.total_evals
        ));
    }
    if left.report.total_budget() != right.report.total_budget() {
        return Some(format!(
            "total-budget:{}!={}",
            left.report.total_budget(),
            right.report.total_budget()
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
    if left.report.records().len() != right.report.records().len() {
        return Some(format!(
            "record-count:{}!={}",
            left.report.records().len(),
            right.report.records().len()
        ));
    }
    if left.report.best_restart() != right.report.best_restart() {
        return Some(format!(
            "best-restart:{}!={}",
            left.report.best_restart(),
            right.report.best_restart()
        ));
    }
    for (restart, (left, right)) in left
        .report
        .records()
        .iter()
        .zip(right.report.records())
        .enumerate()
    {
        if let Some(mismatch) = restart_record_mismatch(restart, left, right) {
            return Some(mismatch);
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
    let expected_admission = admit_bipop(&X0, SIGMA0, TOTAL_BUDGET, Some(F_TARGET), run.input_seed)
        .expect("the fixed study input domain is admitted");
    if let Some(field) = admission_mismatch(&run.admission, &expected_admission) {
        return Err(AdmissionError::AdmissionReceiptMismatch { field });
    }
    let expected_fixture = fixture_identity(run.input_seed, &expected_admission);
    if run.fixture.canonical_bytes() != expected_fixture.canonical_bytes() {
        return Err(AdmissionError::FixtureIdentityMismatch {
            declared: run.fixture.root(),
            computed: expected_fixture.root(),
        });
    }
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

fn realized_stream_block_witnesses(run: &StudyRun) -> (u128, u128) {
    let restart_blocks = 2_u128
        .checked_mul(DIMENSION as u128)
        .and_then(|blocks| blocks.checked_mul(run.report.records().len().saturating_sub(1) as u128))
        .expect("fixed study restart witness fits u128");
    let max_cma_blocks = run
        .report
        .records()
        .iter()
        .map(|record| {
            2_u128
                .checked_mul(DIMENSION as u128)
                .and_then(|blocks| blocks.checked_mul(record.lambda() as u128))
                .and_then(|blocks| blocks.checked_mul(record.report().generations as u128))
                .expect("fixed study CMA witness fits u128")
        })
        .max()
        .unwrap_or(0);
    (restart_blocks, max_cma_blocks)
}

#[allow(clippy::too_many_lines)] // The emitted object is one field-complete receipt schema.
fn emit_green_receipt(run: &StudyRun) -> Event {
    let admission = &run.admission;
    let jacobi = admission
        .jacobi_admission()
        .expect("the retained study reaches the admitted Jacobi dependency");
    let (_, restart_checkpoint) =
        expected_restart_starts(run.input_seed, run.report.records().len());
    let (realized_restart_stream_blocks, max_realized_cma_stream_blocks) =
        realized_stream_block_witnesses(run);
    let mut emitter = Emitter::new(SUITE, CASE);
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "bipop-cma-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{fixture}\",",
                    "\"result_identity\":\"{result}\",",
                    "\"fixture_blake3\":\"{fixture_blake3}\",",
                    "\"result_blake3\":\"{result_blake3}\",",
                    "\"strong_identity_domain\":\"{strong_identity_domain}\",",
                    "\"strong_event_identity_domain\":\"{strong_event_identity_domain}\",",
                    "\"source_file_identity_domain\":\"{source_file_identity_domain}\",",
                    "\"source_file_count\":{source_file_count},",
                    "\"algorithm\":\"fs_dfo::bipop_cmaes\",\"objective\":\"shifted-rastrigin\",",
                    "\"semantic_oracle\":\"{semantic_oracle}\",\"units\":\"dimensionless\",",
                    "\"input_seed\":{input_seed},\"corruption_seed\":{corruption_seed},",
                    "\"dimension\":{dimension},\"total_budget\":{total_budget},",
                    "\"total_evals\":{total_evals},\"schedule\":{schedule},",
                    "\"admission\":{{\"schema_version\":{admission_schema},",
                    "\"supported_schema_version\":{supported_admission_schema},",
                    "\"stream_semantics_version\":{stream_semantics},",
                    "\"supported_stream_semantics_version\":{supported_stream_semantics},",
                    "\"stream_position_domain\":\"{stream_position_domain}\",",
                    "\"checkpoint_version\":{checkpoint_version},",
                    "\"supported_checkpoint_version\":{supported_checkpoint_version},",
                    "\"checkpoint_identity_domain\":\"{checkpoint_identity_domain}\",",
                    "\"dimension\":{admission_dimension},",
                    "\"total_budget\":{admission_budget},",
                    "\"base_lambda\":{base_lambda},",
                    "\"max_large_lambda\":{max_large_lambda},",
                    "\"max_local_budget\":{max_local_budget},",
                    "\"max_restart_ordinal\":{max_restart_ordinal},",
                    "\"max_matrix_entries\":{max_matrix_entries},",
                    "\"max_population_entries\":{max_population_entries},",
                    "\"max_restart_stream_blocks\":\"{max_restart_stream_blocks}\",",
                    "\"max_cma_stream_blocks\":\"{max_cma_stream_blocks}\",",
                    "\"realized_restart_stream_blocks\":\"{realized_restart_stream_blocks}\",",
                    "\"max_realized_cma_stream_blocks\":\"{max_realized_cma_stream_blocks}\",",
                    "\"jacobi_present\":{jacobi_present},\"jacobi\":{{",
                    "\"schema_version\":{jacobi_schema},",
                    "\"supported_schema_version\":{supported_jacobi_schema},",
                    "\"dimension\":{jacobi_dimension},",
                    "\"matrix_entries\":{jacobi_matrix_entries},",
                    "\"aggregate_work_elements\":{jacobi_work},",
                    "\"work_element_cap\":{jacobi_cap}}}}},",
                    "\"restart_schema_version\":{restart_schema},",
                    "\"record_count\":{record_count},\"best_restart\":{best_restart},",
                    "\"ledger\":\"complete-ordered-public-records\",",
                    "\"best\":{{\"x_len\":{best_x_len},",
                    "\"f_bits\":\"0x{best_f_bits:016x}\",",
                    "\"evals\":{best_evals},\"generations\":{best_generations},",
                    "\"converged\":{best_converged},",
                    "\"sigma_bits\":\"0x{best_sigma_bits:016x}\"}},",
                    "\"trace_len\":{trace_len},",
                    "\"target\":{{\"architecture\":\"{target_architecture}\",",
                    "\"operating_system\":\"{target_operating_system}\",",
                    "\"pointer_width_bits\":{target_pointer_width_bits},",
                    "\"endianness\":\"{target_endianness}\"}},",
                    "\"determinism_scope\":\"same-source-snapshot-same-target-test-binary\",",
                    "\"compiler_fingerprint\":\"not-bound-no-claim\",\"versions\":{{",
                    "\"fs_dfo\":\"{fs_dfo_version}\",\"fs_la\":\"{fs_la_version}\",",
                    "\"fs_math\":\"{fs_math_version}\",",
                    "\"fs_rand\":\"{fs_rand_version}\",\"fs_obs\":\"{fs_obs_version}\"}},",
                    "\"no_claims\":[\"optimizer-quality\",\"all-objectives\",",
                    "\"all-dimensions\",\"all-budgets\",\"all-seeds\",",
                    "\"cross-ISA-equality\",\"cancellation\",\"checkpointing\",",
                    "\"compiler-build-identity\",\"signed-authentication\",\"performance\"]}}"
                ),
                fixture = run.fixture.hex(),
                result = run.result.hex(),
                fixture_blake3 = strong_identity_hash(&run.fixture),
                result_blake3 = strong_identity_hash(&run.result),
                strong_identity_domain = STRONG_IDENTITY_DOMAIN,
                strong_event_identity_domain = STRONG_EVENT_IDENTITY_DOMAIN,
                source_file_identity_domain = SOURCE_FILE_IDENTITY_DOMAIN,
                source_file_count = SOURCE_SNAPSHOT_FILES.len(),
                semantic_oracle = SEMANTIC_ORACLE_VERSION,
                input_seed = run.input_seed,
                corruption_seed = CORRUPTION_SEED,
                dimension = DIMENSION,
                total_budget = run.report.total_budget(),
                total_evals = run.report.total_evals,
                schedule = schedule_json(&run.report.schedule),
                admission_schema = admission.schema_version(),
                supported_admission_schema = BIPOP_ADMISSION_SCHEMA_VERSION,
                stream_semantics = admission.stream_semantics_version(),
                supported_stream_semantics = fs_rand::STREAM_SEMANTICS_VERSION,
                stream_position_domain = fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
                checkpoint_version = restart_checkpoint.checkpoint_version,
                supported_checkpoint_version = fs_rand::STREAM_CHECKPOINT_VERSION,
                checkpoint_identity_domain = fs_rand::STREAM_CHECKPOINT_IDENTITY_DOMAIN,
                admission_dimension = admission.dimension(),
                admission_budget = admission.total_budget(),
                base_lambda = admission.base_lambda(),
                max_large_lambda = admission.max_large_lambda(),
                max_local_budget = admission.max_local_budget(),
                max_restart_ordinal = admission.max_restart_ordinal(),
                max_matrix_entries = admission.max_matrix_entries(),
                max_population_entries = admission.max_population_entries(),
                max_restart_stream_blocks = admission.max_restart_stream_blocks(),
                max_cma_stream_blocks = admission.max_cma_stream_blocks(),
                realized_restart_stream_blocks = realized_restart_stream_blocks,
                max_realized_cma_stream_blocks = max_realized_cma_stream_blocks,
                jacobi_present = admission.jacobi_admission().is_some(),
                jacobi_schema = jacobi.schema_version(),
                supported_jacobi_schema = fs_la::eigen::JACOBI_EIGH_ADMISSION_SCHEMA_VERSION,
                jacobi_dimension = jacobi.dimension(),
                jacobi_matrix_entries = jacobi.matrix_entries(),
                jacobi_work = jacobi.aggregate_work_elements(),
                jacobi_cap = jacobi.work_element_cap(),
                restart_schema = BIPOP_RESTART_SCHEMA_VERSION,
                record_count = run.report.records().len(),
                best_restart = run.report.best_restart(),
                best_x_len = run.report.best.x_best.len(),
                best_f_bits = run.report.best.f_best.to_bits(),
                best_evals = run.report.best.evals,
                best_generations = run.report.best.generations,
                best_converged = run.report.best.converged,
                best_sigma_bits = run.report.best.sigma.to_bits(),
                trace_len = run.evaluations.len(),
                target_architecture = std::env::consts::ARCH,
                target_operating_system = std::env::consts::OS,
                target_pointer_width_bits = usize::BITS,
                target_endianness = if cfg!(target_endian = "little") {
                    "little"
                } else {
                    "big"
                },
                fs_dfo_version = fs_dfo::VERSION,
                fs_la_version = fs_la::VERSION,
                fs_math_version = fs_math::VERSION,
                fs_rand_version = fs_rand::VERSION,
                fs_obs_version = fs_obs::VERSION,
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
    event
}

fn assert_green_receipt_json_mutation_is_detected(canonical: &Event) {
    let retained = canonical.content_identity_receipt();
    let mutations = [
        (
            "\"admission\":{\"schema_version\":3,",
            "\"admission\":{\"schema_version\":4,",
        ),
        (
            "\"stream_semantics_version\":1,",
            "\"stream_semantics_version\":2,",
        ),
        ("\"checkpoint_version\":1,", "\"checkpoint_version\":2,"),
        (
            "\"checkpoint_identity_domain\":\"org.frankensim.fs-rand.stream-checkpoint.v1\"",
            "\"checkpoint_identity_domain\":\"org.frankensim.fs-rand.stream-checkpoint.v2\"",
        ),
        (
            "\"dimension\":4,\"total_budget\":6000,\"base_lambda\":8,",
            "\"dimension\":4,\"total_budget\":6001,\"base_lambda\":8,",
        ),
        ("\"jacobi_present\":true", "\"jacobi_present\":false"),
        (
            "\"jacobi\":{\"schema_version\":1,",
            "\"jacobi\":{\"schema_version\":2,",
        ),
        (
            "\"max_cma_stream_blocks\":\"47872\"",
            "\"max_cma_stream_blocks\":\"47873\"",
        ),
        (
            "\"realized_restart_stream_blocks\":\"",
            "\"realized_restart_stream_blocks\":\"0",
        ),
        (
            "\"max_realized_cma_stream_blocks\":\"",
            "\"max_realized_cma_stream_blocks\":\"0",
        ),
        (
            "\"strong_identity_domain\":\"frankensim.fs-dfo.bipop-study.replay-identity.v1\"",
            "\"strong_identity_domain\":\"frankensim.fs-dfo.bipop-study.replay-identity.v2\"",
        ),
        (
            "\"strong_event_identity_domain\":\"frankensim.fs-dfo.bipop-study.event-content-identity.v1\"",
            "\"strong_event_identity_domain\":\"frankensim.fs-dfo.bipop-study.event-content-identity.v2\"",
        ),
        (
            "\"source_file_identity_domain\":\"frankensim.fs-dfo.bipop-study.source-file-identity.v1\"",
            "\"source_file_identity_domain\":\"frankensim.fs-dfo.bipop-study.source-file-identity.v2\"",
        ),
        (
            "\"target\":{\"architecture\":\"",
            "\"target\":{\"architecture\":\"mutated-",
        ),
    ];
    for (needle, replacement) in mutations {
        let mut mutant = canonical.clone();
        let EventKind::Custom { json, .. } = &mut mutant.kind else {
            panic!("green BIPOP receipt must be a Custom event");
        };
        assert!(
            json.contains(needle),
            "receipt is missing mutation field {needle}"
        );
        *json = json.replacen(needle, replacement, 1);
        let line = mutant.to_jsonl();
        fs_obs::validate_line(&line)
            .expect("mutated receipt remains structurally valid fs-obs wire JSON");
        mutant
            .admit_content_identity(&retained)
            .expect_err("an admission-envelope mutation must stale the retained event identity");
    }
}

fn green_verdict_detail(run: &StudyRun) -> String {
    let jacobi = run
        .admission
        .jacobi_admission()
        .expect("the retained study reaches the Jacobi dependency");
    let (realized_restart_blocks, max_realized_cma_blocks) = realized_stream_block_witnesses(run);
    format!(
        "fixture={}; result={}; fixture_blake3={}; result_blake3={}; admission_schema={}; stream_semantics={}; checkpoint_version={}; jacobi_present={}; jacobi_schema={}; restart_blocks={}/{}; max_CMA_blocks={}/{}; total_budget={}; total_evals={}; records={}; best_restart={}; target={}-{}-{}-bit-{}-endian; source_files={}; trace=bit-exact; ordered_ledger=fully-bound; legacy_projections=exact; determinism=same-source-snapshot-same-target-test-binary; compiler-fingerprint=no-claim; semantic_oracle={SEMANTIC_ORACLE_VERSION}",
        run.fixture.hex(),
        run.result.hex(),
        strong_identity_hash(&run.fixture),
        strong_identity_hash(&run.result),
        run.admission.schema_version(),
        run.admission.stream_semantics_version(),
        fs_rand::STREAM_CHECKPOINT_VERSION,
        run.admission.jacobi_admission().is_some(),
        jacobi.schema_version(),
        realized_restart_blocks,
        run.admission.max_restart_stream_blocks(),
        max_realized_cma_blocks,
        run.admission.max_cma_stream_blocks(),
        run.report.total_budget(),
        run.report.total_evals,
        run.report.records().len(),
        run.report.best_restart(),
        std::env::consts::ARCH,
        std::env::consts::OS,
        usize::BITS,
        if cfg!(target_endian = "little") {
            "little"
        } else {
            "big"
        },
        SOURCE_SNAPSHOT_FILES.len(),
    )
}

fn emit_green_verdict(run: &StudyRun) -> Event {
    let detail = green_verdict_detail(run);
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail,
            seed: run.input_seed,
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
        suite,
        case,
        pass,
        detail,
        seed,
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    if let Err(error) = admit_against(run, reference) {
        panic!("merge gate refused {case}: {error:?}; {detail}");
    }
    if let Some(mismatch) = accounting_mismatch(run) {
        panic!("merge gate refused {case}: semantic mismatch {mismatch}; {detail}");
    }
    assert_eq!(
        event.session.as_str(),
        SUITE,
        "merge gate refused an event from the wrong session"
    );
    let expected_scope = format!("{CASE}/verdict");
    assert_eq!(
        event.scope.as_str(),
        expected_scope.as_str(),
        "merge gate refused an event from the wrong scope"
    );
    assert_eq!(
        event.seq, 0,
        "merge gate requires the canonical verdict slot"
    );
    assert_eq!(
        event.severity,
        Severity::Info,
        "merge gate requires an informational green verdict"
    );
    assert!(
        event.wall_ns.is_none(),
        "merge gate requires deterministic evidence without a wall-clock envelope"
    );
    assert_eq!(suite.as_str(), SUITE, "merge gate refused the wrong suite");
    assert_eq!(case.as_str(), CASE, "merge gate refused the wrong case");
    assert_eq!(
        *seed, run.input_seed,
        "merge gate refused a verdict with the wrong causal input seed"
    );
    let expected_detail = green_verdict_detail(run);
    assert_eq!(
        detail.as_str(),
        expected_detail.as_str(),
        "merge gate refused a verdict that does not name the admitted run"
    );
    assert!(
        *pass,
        "merge gate requires a passing verdict for {case}: {detail}"
    );
}

fn assert_verdict_mutation_refused(
    run: &StudyRun,
    reference: &ReplayIdentity,
    canonical: &Event,
    expected_message: &str,
    mutate: impl FnOnce(&mut Event),
) {
    let mut mutant = canonical.clone();
    mutate(&mut mutant);
    let panic = catch_unwind(|| {
        assert_mergeable(run, reference, &mutant);
    })
    .expect_err("a mutated verdict envelope must fail closed");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("verdict-refusal panic carries text");
    assert!(message.contains(expected_message), "{message}");
}

fn assert_verdict_envelope_mutations_refused(
    run: &StudyRun,
    reference: &ReplayIdentity,
    canonical: &Event,
) {
    assert_verdict_mutation_refused(run, reference, canonical, "wrong session", |event| {
        event.session.push_str("/wrong");
    });
    assert_verdict_mutation_refused(run, reference, canonical, "wrong scope", |event| {
        event.scope.push_str("/wrong");
    });
    assert_verdict_mutation_refused(
        run,
        reference,
        canonical,
        "canonical verdict slot",
        |event| event.seq = event.seq.wrapping_add(1),
    );
    assert_verdict_mutation_refused(
        run,
        reference,
        canonical,
        "informational green verdict",
        |event| event.severity = Severity::Warn,
    );
    assert_verdict_mutation_refused(
        run,
        reference,
        canonical,
        "without a wall-clock envelope",
        |event| event.wall_ns = Some(1),
    );
    assert_verdict_mutation_refused(run, reference, canonical, "wrong suite", |event| {
        let EventKind::ConformanceCase { suite, .. } = &mut event.kind else {
            unreachable!("canonical verdict is a ConformanceCase");
        };
        suite.push_str("/wrong");
    });
    assert_verdict_mutation_refused(run, reference, canonical, "wrong case", |event| {
        let EventKind::ConformanceCase { case, .. } = &mut event.kind else {
            unreachable!("canonical verdict is a ConformanceCase");
        };
        case.push_str("/wrong");
    });
    assert_verdict_mutation_refused(run, reference, canonical, "passing verdict", |event| {
        let EventKind::ConformanceCase { pass, .. } = &mut event.kind else {
            unreachable!("canonical verdict is a ConformanceCase");
        };
        *pass = false;
    });
    assert_verdict_mutation_refused(
        run,
        reference,
        canonical,
        "does not name the admitted run",
        |event| {
            let EventKind::ConformanceCase { detail, .. } = &mut event.kind else {
                unreachable!("canonical verdict is a ConformanceCase");
            };
            detail.push_str("; wrong");
        },
    );
    assert_verdict_mutation_refused(
        run,
        reference,
        canonical,
        "wrong causal input seed",
        |event| {
            let EventKind::ConformanceCase { seed, .. } = &mut event.kind else {
                unreachable!("canonical verdict is a ConformanceCase");
            };
            *seed ^= 1;
        },
    );
    assert_verdict_mutation_refused(run, reference, canonical, "only ConformanceCase", |event| {
        event.kind = EventKind::Custom {
            name: "wrong-kind".to_string(),
            json: "{}".to_string(),
        };
    });
}

fn assert_resealed_semantic_refusal(
    mut mutant: StudyRun,
    event: &Event,
    expected_mismatch_fragment: &str,
) {
    let stale = validate_payload(&mutant)
        .expect_err("an unsealed semantic mutation must fail payload admission");
    assert!(
        matches!(stale, AdmissionError::PayloadIdentityMismatch { .. }),
        "semantic mutation must retain a valid fixture but stale the result: {stale:?}"
    );
    mutant.result = result_identity(&mutant.fixture, &mutant.report, &mutant.evaluations);
    validate_payload(&mutant).expect("resealed semantic mutant must be identity-consistent");
    let panic = catch_unwind(|| {
        assert_mergeable(&mutant, &mutant.result, event);
    })
    .expect_err("resealed semantic mutant must fail admission against its own identity");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("semantic-refusal panic carries text");
    assert!(message.contains("semantic mismatch"), "{message}");
    assert!(message.contains(expected_mismatch_fragment), "{message}");
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
        "input_seed=0x{:016x}; corruption_seed=0x{:016x}; fixture={}; fixture_blake3={}; coordinate={}; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}; canonical={}; canonical_blake3={}; corrupted={}; corrupted_blake3={}",
        canonical.input_seed,
        corruption.mutation.seed,
        canonical.fixture.hex(),
        strong_identity_hash(&canonical.fixture),
        corruption.mutation.coordinate,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.mismatch,
        canonical.result.hex(),
        strong_identity_hash(&canonical.result),
        corruption.run.result.hex(),
        strong_identity_hash(&corruption.run.result),
    )
}

fn exercise_disclosed_corruption(canonical: &StudyRun, replay: &StudyRun) -> ContentHash {
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

    let semantic_panic = catch_unwind(|| {
        assert_mergeable(
            &first_corruption.run,
            &first_corruption.run.result,
            &first_event,
        );
    })
    .expect_err("a resealed mutant must fail semantic admission even against itself");
    let semantic_message = semantic_panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| semantic_panic.downcast_ref::<&str>().copied())
        .expect("semantic merge-gate panic carries text");
    assert!(semantic_message.contains(RED_CASE));
    assert!(semantic_message.contains("semantic mismatch"));

    strong_event_identity_hash(&first_event)
}

#[test]
fn admission_comparison_checks_outer_fields_when_jacobi_is_unreachable() {
    let one = admit_bipop(&[0.0], 0.5, 1, None, INPUT_SEED).expect("budget one admits");
    let two = admit_bipop(&[0.0], 0.5, 2, None, INPUT_SEED).expect("budget two admits");
    assert_eq!(one.jacobi_admission(), None);
    assert_eq!(two.jacobi_admission(), None);
    assert_eq!(
        admission_mismatch(&one, &two),
        Some("total-budget"),
        "absence of a nested capability cannot skip outer receipt comparison"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One test exercises the complete green/red merge contract.
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

    let green_receipt = emit_green_receipt(&first);
    assert_green_receipt_json_mutation_is_detected(&green_receipt);
    let green_verdict = emit_green_verdict(&replay);
    assert_mergeable(&replay, &first.result, &green_verdict);
    assert_verdict_envelope_mutations_refused(&replay, &first.result, &green_verdict);

    let mut objective_mutant = first.clone();
    objective_mutant.evaluations[0].value += 0.25;
    assert_resealed_semantic_refusal(objective_mutant, &green_verdict, "trace[0]-objective");

    let mut schedule_mutant = first.clone();
    schedule_mutant.report.schedule[0] = schedule_mutant.report.schedule[0]
        .checked_mul(2)
        .expect("fixture schedule mutation fits usize");
    assert_resealed_semantic_refusal(schedule_mutant, &green_verdict, "schedule[0]");

    let mut total_mutant = first.clone();
    total_mutant.report.total_evals = total_mutant
        .report
        .total_evals
        .checked_add(1)
        .expect("fixture total-evaluation mutation fits usize");
    assert_resealed_semantic_refusal(total_mutant, &green_verdict, "reported-total-evals");

    let mut admission_budget_mutant = first.clone();
    admission_budget_mutant.admission =
        admit_bipop(&X0, SIGMA0, TOTAL_BUDGET + 1, Some(F_TARGET), INPUT_SEED)
            .expect("adjacent-budget admission remains structurally valid");
    assert_eq!(
        validate_payload(&admission_budget_mutant),
        Err(AdmissionError::AdmissionReceiptMismatch {
            field: "total-budget"
        })
    );
    admission_budget_mutant.fixture = fixture_identity(
        admission_budget_mutant.input_seed,
        &admission_budget_mutant.admission,
    );
    admission_budget_mutant.result = result_identity(
        &admission_budget_mutant.fixture,
        &admission_budget_mutant.report,
        &admission_budget_mutant.evaluations,
    );
    assert_eq!(
        validate_payload(&admission_budget_mutant),
        Err(AdmissionError::AdmissionReceiptMismatch {
            field: "total-budget"
        }),
        "self-consistently resealing the wrong public receipt cannot bypass input admission"
    );

    let mut jacobi_receipt_mutant = first.clone();
    jacobi_receipt_mutant.admission = admit_bipop(
        &[0.0; DIMENSION + 1],
        SIGMA0,
        TOTAL_BUDGET,
        Some(F_TARGET),
        INPUT_SEED,
    )
    .expect("adjacent-dimension admission remains structurally valid");
    assert_eq!(
        validate_payload(&jacobi_receipt_mutant),
        Err(AdmissionError::AdmissionReceiptMismatch {
            field: "jacobi.dimension"
        })
    );

    let mut jacobi_presence_mutant = first.clone();
    jacobi_presence_mutant.admission = admit_bipop(
        &X0,
        SIGMA0,
        expected_base_lambda(),
        Some(F_TARGET),
        INPUT_SEED,
    )
    .expect("callback-only BIPOP fixture admits without Jacobi authority");
    assert_eq!(jacobi_presence_mutant.admission.jacobi_admission(), None);
    assert_eq!(
        validate_payload(&jacobi_presence_mutant),
        Err(AdmissionError::AdmissionReceiptMismatch {
            field: "jacobi.presence"
        })
    );
    let canonical_fixture_hash = strong_identity_hash(&first.fixture);
    jacobi_presence_mutant.fixture = fixture_identity(
        jacobi_presence_mutant.input_seed,
        &jacobi_presence_mutant.admission,
    );
    jacobi_presence_mutant.result = result_identity(
        &jacobi_presence_mutant.fixture,
        &jacobi_presence_mutant.report,
        &jacobi_presence_mutant.evaluations,
    );
    assert_ne!(
        strong_identity_hash(&jacobi_presence_mutant.fixture),
        canonical_fixture_hash,
        "Jacobi authority presence is semantic fixture provenance"
    );
    assert_eq!(
        validate_payload(&jacobi_presence_mutant),
        Err(AdmissionError::AdmissionReceiptMismatch {
            field: "jacobi.presence"
        }),
        "resealing the absent authority cannot bypass fixed-input admission"
    );

    let (alternate_budget_report, alternate_budget_evaluations) =
        execute_study_payload(INPUT_SEED, TOTAL_BUDGET + 1);
    let mut report_budget_mutant = first.clone();
    report_budget_mutant.report = alternate_budget_report;
    report_budget_mutant.evaluations = alternate_budget_evaluations;
    assert_resealed_semantic_refusal(report_budget_mutant, &green_verdict, "report-total-budget");

    // Public records are intentionally immutable, so mutate the complete
    // ledger by substituting a separately generated report from a different
    // causal seed while retaining the canonical fixture and callback trace.
    // Identity admission must first see a stale payload; after correct
    // resealing, independent semantics must reject the substituted ledger at
    // its first incompatible record/boundary provenance fact.
    let alternate_ledger = run_study(INPUT_SEED ^ 0x0100_0000_0000_0000);
    let mut ledger_mutant = first.clone();
    ledger_mutant.report = alternate_ledger.report;
    assert_resealed_semantic_refusal(ledger_mutant, &green_verdict, "restart[");

    let mut causal_seed_mutant = first.clone();
    causal_seed_mutant.input_seed ^= 1;
    assert!(matches!(
        validate_payload(&causal_seed_mutant),
        Err(AdmissionError::FixtureIdentityMismatch { declared, computed })
            if declared == first.fixture.root()
                && computed == fixture_identity(
                    causal_seed_mutant.input_seed,
                    &causal_seed_mutant.admission,
                )
                .root()
    ));
    causal_seed_mutant.fixture =
        fixture_identity(causal_seed_mutant.input_seed, &causal_seed_mutant.admission);
    assert_resealed_semantic_refusal(causal_seed_mutant, &green_verdict, "restart[");

    let observed = [
        strong_identity_hash(&first.fixture).to_string(),
        strong_identity_hash(&first.result).to_string(),
        strong_event_identity_hash(&green_receipt).to_string(),
        strong_event_identity_hash(&green_verdict).to_string(),
        exercise_disclosed_corruption(&first, &replay).to_string(),
    ];
    println!(
        "BIPOP_STUDY_FREEZE fixture_blake3={} result_blake3={} green_receipt_event_blake3={} green_verdict_event_blake3={} red_event_blake3={}",
        observed[0], observed[1], observed[2], observed[3], observed[4],
    );
    if let Some(expected) = goldens::EXPECTED_FOR_TARGET {
        let labels = [
            "fixture_blake3",
            "result_blake3",
            "green_receipt_event_blake3",
            "green_verdict_event_blake3",
            "red_event_blake3",
        ];
        for ((label, observed), expected) in labels.iter().zip(&observed).zip(expected) {
            assert_eq!(
                observed, expected,
                "target-specific BIPOP study golden moved at {label}"
            );
        }
    } else {
        println!(
            "BIPOP_STUDY_UNFROZEN target={}-{}-{}-bit-{}-endian; emitted hashes are sentinels, not regression goldens",
            std::env::consts::ARCH,
            std::env::consts::OS,
            usize::BITS,
            if cfg!(target_endian = "little") {
                "little"
            } else {
                "big"
            },
        );
    }
}
