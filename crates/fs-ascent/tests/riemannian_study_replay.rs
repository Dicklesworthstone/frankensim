//! G5 study-scale replay and seeded-failure self-test for Riemannian L-BFGS.
//!
//! The production engine minimizes a seeded Rayleigh-quotient fixture on the
//! sphere.  The retained receipt binds the complete public optimizer state,
//! complete public report, engine configuration, objective matrix, derived
//! start point, and logical RNG coordinates.  An independently replayed run
//! must produce the same receipt byte for byte.  A deterministic red mutation
//! flips one mantissa bit in the returned decision, remains a valid evidence
//! payload after being resealed, and is refused by the test-local merge gate.
//!
//! This is one sphere fixture.  It does not claim all manifolds, persisted
//! checkpoints, cancellation recovery, cross-ISA equality, or performance.

use fs_ascent::{RiemannianLbfgs, RiemannianReport, StopReason, StopRule};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};
use fs_opt::Manifold;
use fs_rand::StreamKey;

const SUITE: &str = "fs-ascent/riemannian-study-replay";
// The same known-convergent logical fixture coordinates used by the original
// sphere battery, retained here instead of silently inventing another corpus.
const INPUT_SEED: u64 = 41;
const MUTATION_SEED: u64 = 0x5245_445F_4249_5401;
const RNG_KERNEL: u32 = 0xA5C3;
const MATRIX_TILE: u32 = 20;
const START_TILE: u32 = 21;
const DIMENSION: usize = 12;
const MEMORY: usize = 10;
const GRADIENT_TOLERANCE: f64 = 1e-9;
const EVALUATION_BUDGET: usize = 3_000;
const MAX_ITERATIONS: usize = 1_000;

#[derive(Debug)]
struct RunRecord {
    state: RiemannianLbfgs,
    report: RiemannianReport,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    input_seed: u64,
    matrix_bits: Vec<u64>,
    start_bits: Vec<u64>,
    final_x_bits: Vec<u64>,
    final_f_bits: u64,
    final_gradient_bits: Vec<u64>,
    history_bits: Vec<u64>,
    state_iterations: usize,
    state_evaluations: usize,
    stop_reason: &'static str,
    report_gradient_norm_bits: u64,
    report_f_bits: u64,
    report_iterations: usize,
    report_evaluations: usize,
    report_worst_violation_bits: u64,
}

impl ReceiptPayload {
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-riemannian-study-receipt-v1")
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("fs-opt-version", fs_opt::VERSION)
            .str("fs-rand-version", fs_rand::VERSION)
            .str("engine", "RiemannianLbfgs")
            .str("manifold", "Sphere")
            .u64("ambient-dimension", DIMENSION as u64)
            .u64("memory", MEMORY as u64)
            .f64_bits("gradient-tolerance", GRADIENT_TOLERANCE)
            .u64("evaluation-budget", EVALUATION_BUDGET as u64)
            .u64("maximum-iterations", MAX_ITERATIONS as u64)
            .u64("input-seed", self.input_seed)
            .u64("rng-kernel", u64::from(RNG_KERNEL))
            .u64("matrix-tile", u64::from(MATRIX_TILE))
            .u64("start-tile", u64::from(START_TILE))
            .u64("matrix-values", self.matrix_bits.len() as u64);
        for &value_bits in &self.matrix_bits {
            builder = builder.u64("matrix-value-bits", value_bits);
        }
        builder = builder.u64("start-values", self.start_bits.len() as u64);
        for &value_bits in &self.start_bits {
            builder = builder.u64("start-value-bits", value_bits);
        }
        builder = builder.u64("final-values", self.final_x_bits.len() as u64);
        for &value_bits in &self.final_x_bits {
            builder = builder.u64("final-value-bits", value_bits);
        }
        builder = builder.u64("final-objective-bits", self.final_f_bits).u64(
            "final-gradient-values",
            self.final_gradient_bits.len() as u64,
        );
        for &value_bits in &self.final_gradient_bits {
            builder = builder.u64("final-gradient-value-bits", value_bits);
        }
        builder = builder.u64("history-values", self.history_bits.len() as u64);
        for &value_bits in &self.history_bits {
            builder = builder.u64("history-value-bits", value_bits);
        }
        builder
            .u64("state-iterations", self.state_iterations as u64)
            .u64("state-evaluations", self.state_evaluations as u64)
            .str("stop-reason", self.stop_reason)
            .u64("report-gradient-norm-bits", self.report_gradient_norm_bits)
            .u64("report-objective-bits", self.report_f_bits)
            .u64("report-iterations", self.report_iterations as u64)
            .u64("report-evaluations", self.report_evaluations as u64)
            .u64(
                "report-worst-violation-bits",
                self.report_worst_violation_bits,
            )
            .finish()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RetainedReceipt {
    payload: ReceiptPayload,
    declared_identity: ReplayIdentity,
}

impl RetainedReceipt {
    fn new(payload: ReceiptPayload) -> Self {
        let declared_identity = payload.identity();
        Self {
            payload,
            declared_identity,
        }
    }

    fn reseal(&mut self) {
        self.declared_identity = self.payload.identity();
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MergeRefusal {
    UnsupportedIdentityVersion,
    PayloadIdentityMismatch,
    ReferenceIdentityMismatch,
}

fn admit_receipt(
    reference: &ReplayIdentity,
    candidate: &RetainedReceipt,
) -> Result<(), MergeRefusal> {
    check_version(candidate.declared_identity.version())
        .map_err(|_| MergeRefusal::UnsupportedIdentityVersion)?;
    if candidate.payload.identity() != candidate.declared_identity {
        return Err(MergeRefusal::PayloadIdentityMismatch);
    }
    if &candidate.declared_identity != reference {
        return Err(MergeRefusal::ReferenceIdentityMismatch);
    }
    Ok(())
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
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

fn stream_vector(seed: u64, tile: u32, length: usize) -> Vec<f64> {
    let mut stream = StreamKey {
        seed,
        kernel: RNG_KERNEL,
        tile,
    }
    .stream();
    (0..length)
        .map(|_| 2.0f64.mul_add(stream.next_f64(), -1.0))
        .collect()
}

fn seeded_fixture(seed: u64) -> (Vec<f64>, Vec<f64>) {
    let raw = stream_vector(seed, MATRIX_TILE, DIMENSION * DIMENSION);
    let mut matrix = vec![0.0; DIMENSION * DIMENSION];
    for row in 0..DIMENSION {
        for column in 0..DIMENSION {
            matrix[row * DIMENSION + column] =
                raw[row * DIMENSION + column] + raw[column * DIMENSION + row];
        }
        matrix[row * DIMENSION + row] += 2.0;
    }

    let mut start = stream_vector(seed, START_TILE, DIMENSION);
    let norm = fs_math::det::sqrt(start.iter().map(|value| value * value).sum());
    assert!(
        norm.is_finite() && norm > 0.0,
        "seeded start must normalize"
    );
    for value in &mut start {
        *value /= norm;
    }
    (matrix, start)
}

fn rayleigh(matrix: &[f64], x: &[f64]) -> (f64, Vec<f64>) {
    assert_eq!(matrix.len(), DIMENSION * DIMENSION);
    assert_eq!(x.len(), DIMENSION);
    let mut ax = vec![0.0; DIMENSION];
    for row in 0..DIMENSION {
        for column in 0..DIMENSION {
            ax[row] = matrix[row * DIMENSION + column].mul_add(x[column], ax[row]);
        }
    }
    let objective = x
        .iter()
        .zip(&ax)
        .map(|(coordinate, product)| coordinate * product)
        .sum();
    let gradient = ax.into_iter().map(|value| 2.0 * value).collect();
    (objective, gradient)
}

fn run_once(matrix: &[f64], start: &[f64]) -> RunRecord {
    let mut objective = |x: &[f64]| rayleigh(matrix, x);
    let ambient = u32::try_from(DIMENSION).expect("fixture dimension must fit the manifold API");
    let mut state =
        RiemannianLbfgs::new(Manifold::Sphere { ambient }, start, MEMORY, &mut objective);
    let stop = StopRule::Any(vec![
        StopRule::GradNorm(GRADIENT_TOLERANCE),
        StopRule::Budget(EVALUATION_BUDGET),
    ]);
    let report = state.run(&mut objective, &stop, MAX_ITERATIONS);
    RunRecord { state, report }
}

fn receipt(matrix: &[f64], start: &[f64], run: &RunRecord) -> RetainedReceipt {
    let ambient = u32::try_from(DIMENSION).expect("fixture dimension must fit the manifold API");
    assert!(matches!(
        &run.state.manifold,
        Manifold::Sphere { ambient: actual } if *actual == ambient
    ));
    RetainedReceipt::new(ReceiptPayload {
        input_seed: INPUT_SEED,
        matrix_bits: bits(matrix),
        start_bits: bits(start),
        final_x_bits: bits(&run.state.x),
        final_f_bits: run.state.f.to_bits(),
        final_gradient_bits: bits(&run.state.g),
        history_bits: bits(&run.state.history),
        state_iterations: run.state.iters,
        state_evaluations: run.state.evals,
        stop_reason: stop_reason_name(&run.report.reason),
        report_gradient_norm_bits: run.report.grad_norm.to_bits(),
        report_f_bits: run.report.f.to_bits(),
        report_iterations: run.report.iters,
        report_evaluations: run.report.evals,
        report_worst_violation_bits: run.report.worst_violation.to_bits(),
    })
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, usize, u64) {
    let mut mutant = receipt.clone();
    let coordinate = (MUTATION_SEED as usize) % mutant.payload.final_x_bits.len();
    let mantissa_bit = (MUTATION_SEED >> 8) % 52;
    let mask = 1_u64 << mantissa_bit;
    mutant.payload.final_x_bits[coordinate] ^= mask;
    assert!(
        f64::from_bits(mutant.payload.final_x_bits[coordinate]).is_finite(),
        "mantissa-only mutation must remain a finite wire-valid decision"
    );
    mutant.reseal();
    (mutant, coordinate, mask)
}

fn emit_receipt(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    coordinate: usize,
    mask: u64,
) {
    let json = format!(
        "{{\"input_seed\":{INPUT_SEED},\"mutation_seed\":{MUTATION_SEED},\
         \"reference_identity\":\"{}\",\"mutant_identity\":\"{}\",\
         \"mutated_coordinate\":{coordinate},\"mantissa_mask\":\"{mask:#018x}\",\
         \"merge_refusal\":\"reference-identity-mismatch\"}}",
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
    );
    let mut emitter = Emitter::new(SUITE, "sphere-rayleigh-full-study");
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "riemannian-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("Riemannian study receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "sphere-rayleigh-full-study".to_string(),
            pass: true,
            detail: format!(
                "input seed {INPUT_SEED:#018x} replayed every public state/report bit; mutation seed {MUTATION_SEED:#018x} flipped coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and the merge gate refused it",
                mutant.declared_identity.hex(),
            ),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("Riemannian seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("Riemannian seeded-failure verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
fn riemannian_sphere_study_replays_and_rejects_seeded_red_mutation() {
    let (matrix, start) = seeded_fixture(INPUT_SEED);
    let (eigenvalues, _) = fs_la::eigen::jacobi_eigh(&matrix, DIMENSION);
    let expected_minimum = eigenvalues.iter().copied().fold(f64::INFINITY, f64::min);

    let reference_run = run_once(&matrix, &start);
    assert_eq!(reference_run.report.reason, StopReason::GradNorm);
    assert!(
        (reference_run.report.f - expected_minimum).abs() < 1e-7,
        "Rayleigh minimum {} differs from Jacobi oracle {expected_minimum}",
        reference_run.report.f,
    );
    assert!(
        reference_run.report.worst_violation < 1e-14,
        "sphere retraction violated the manifold: {}",
        reference_run.report.worst_violation,
    );
    assert_eq!(
        reference_run.state.history.len(),
        reference_run.state.iters + 1,
        "every completed state transition must retain its objective"
    );
    assert_eq!(reference_run.report.iters, reference_run.state.iters);
    assert_eq!(reference_run.report.evals, reference_run.state.evals);
    assert_eq!(
        reference_run.report.f.to_bits(),
        reference_run.state.f.to_bits()
    );

    let reference = receipt(&matrix, &start, &reference_run);
    admit_receipt(&reference.declared_identity, &reference)
        .expect("the internally consistent reference receipt must admit");

    let replay_run = run_once(&matrix, &start);
    let replay = receipt(&matrix, &start, &replay_run);
    assert_eq!(
        replay, reference,
        "same logical seed must replay the complete public study receipt"
    );

    let (mutant, coordinate, mask) = mutate_returned_decision(&reference);
    let (mutant_repeat, repeat_coordinate, repeat_mask) = mutate_returned_decision(&reference);
    assert_eq!((coordinate, mask), (repeat_coordinate, repeat_mask));
    assert_eq!(
        mutant, mutant_repeat,
        "the seeded red mutation and its evidence identity must be stable"
    );
    assert_ne!(
        mutant.declared_identity, reference.declared_identity,
        "one returned-decision bit must move the retained evidence identity"
    );
    let mut stale_identity_mutant = mutant.clone();
    stale_identity_mutant.declared_identity = reference.declared_identity.clone();
    assert_eq!(
        admit_receipt(&reference.declared_identity, &stale_identity_mutant),
        Err(MergeRefusal::PayloadIdentityMismatch),
        "changing result bits without resealing must fail the payload-integrity check"
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, &mutant),
        Err(MergeRefusal::ReferenceIdentityMismatch),
        "a self-consistent but semantically different result must fail closed"
    );

    emit_receipt(&reference, &mutant, coordinate, mask);
}
