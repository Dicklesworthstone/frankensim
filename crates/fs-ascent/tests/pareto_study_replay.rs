//! G5 full-study replay and seeded-failure self-test for Pareto tracing.
//!
//! The production weighted-sum continuation traces the known convex quadratic
//! front, then production epsilon-constraint continuation traces the known
//! concave Fonseca-Fleming front. The retained receipt binds both schedules,
//! starts, every objective callback input/value/gradient, and every public
//! `ParetoPoint` decision/objective/gradient/KKT field. An independent repeat
//! must reproduce the receipt byte for byte. A deterministic red mutation
//! flips one finite returned-decision mantissa bit and is refused both before
//! and after self-consistent resealing.
//!
//! This is the two-objective tracing family only. It does not claim
//! tri-objective behavior, the full WFG transformation stack, cancellation,
//! checkpointing, cross-ISA equality, persistence, or performance.

use core::cell::RefCell;

use fs_ascent::{ParetoPoint, epsilon_constraint_sweep, weighted_sum_sweep};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-ascent/pareto-study-replay";
const INPUT_SEED: u64 = 0;
const MUTATION_SEED: u64 = 0x5041_5245_544F_5244;
const EPSILON_TOLERANCE: f64 = 1e-7;
const WEIGHTED_START: [f64; 3] = [0.5, 0.5, 0.5];
const EPSILON_START: [f64; 2] = [0.0, 0.0];

#[derive(Clone, Debug, PartialEq, Eq)]
struct ObjectiveCall {
    phase: &'static str,
    objective: &'static str,
    point_bits: Vec<u64>,
    value_bits: u64,
    gradient_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct PointPayload {
    x_bits: Vec<u64>,
    objective_bits: [u64; 2],
    kkt_bits: Option<[u64; 4]>,
    gradient_norm_bits: u64,
}

impl From<&ParetoPoint> for PointPayload {
    fn from(point: &ParetoPoint) -> Self {
        Self {
            x_bits: bits(&point.x),
            objective_bits: [point.f[0].to_bits(), point.f[1].to_bits()],
            kkt_bits: point.kkt.as_ref().map(|kkt| {
                [
                    kkt.stationarity.to_bits(),
                    kkt.feasibility.to_bits(),
                    kkt.dual_feasibility.to_bits(),
                    kkt.complementarity.to_bits(),
                ]
            }),
            gradient_norm_bits: point.grad_norm.to_bits(),
        }
    }
}

#[derive(Debug)]
struct RunRecord {
    objective_calls: Vec<ObjectiveCall>,
    weighted_points: Vec<PointPayload>,
    epsilon_points: Vec<PointPayload>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    input_seed: u64,
    weights_bits: Vec<u64>,
    weighted_start_bits: Vec<u64>,
    epsilon_bits: Vec<u64>,
    epsilon_start_bits: Vec<u64>,
    objective_calls: Vec<ObjectiveCall>,
    weighted_points: Vec<PointPayload>,
    epsilon_points: Vec<PointPayload>,
}

impl ReceiptPayload {
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-pareto-study-receipt-v1")
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("family", "two-objective-pareto-tracing")
            .str("weighted-engine", "weighted_sum_sweep/L-BFGS")
            .str(
                "epsilon-engine",
                "epsilon_constraint_sweep/augmented-lagrangian",
            )
            .u64("input-seed", self.input_seed)
            .f64_bits("epsilon-tolerance", EPSILON_TOLERANCE)
            .u64("weights", self.weights_bits.len() as u64);
        for &value_bits in &self.weights_bits {
            builder = builder.u64("weight-bits", value_bits);
        }
        builder = builder.u64(
            "weighted-start-values",
            self.weighted_start_bits.len() as u64,
        );
        for &value_bits in &self.weighted_start_bits {
            builder = builder.u64("weighted-start-bits", value_bits);
        }
        builder = builder.u64("epsilons", self.epsilon_bits.len() as u64);
        for &value_bits in &self.epsilon_bits {
            builder = builder.u64("epsilon-bits", value_bits);
        }
        builder = builder.u64("epsilon-start-values", self.epsilon_start_bits.len() as u64);
        for &value_bits in &self.epsilon_start_bits {
            builder = builder.u64("epsilon-start-bits", value_bits);
        }

        builder = builder.u64("objective-calls", self.objective_calls.len() as u64);
        for call in &self.objective_calls {
            builder = builder
                .str("call-phase", call.phase)
                .str("call-objective", call.objective)
                .u64("call-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("call-point-bits", value_bits);
            }
            builder = builder
                .u64("call-value-bits", call.value_bits)
                .u64("call-gradient-values", call.gradient_bits.len() as u64);
            for &value_bits in &call.gradient_bits {
                builder = builder.u64("call-gradient-bits", value_bits);
            }
        }

        builder = append_points(builder, "weighted", &self.weighted_points);
        append_points(builder, "epsilon", &self.epsilon_points).finish()
    }
}

fn append_points(
    mut builder: IdentityBuilder,
    path: &'static str,
    points: &[PointPayload],
) -> IdentityBuilder {
    builder = builder
        .str("point-path", path)
        .u64("points", points.len() as u64);
    for point in points {
        builder = builder.u64("point-x-values", point.x_bits.len() as u64);
        for &value_bits in &point.x_bits {
            builder = builder.u64("point-x-bits", value_bits);
        }
        for &value_bits in &point.objective_bits {
            builder = builder.u64("point-objective-bits", value_bits);
        }
        builder = builder.flag("point-has-kkt", point.kkt_bits.is_some());
        if let Some(kkt_bits) = point.kkt_bits {
            for value_bits in kkt_bits {
                builder = builder.u64("point-kkt-bits", value_bits);
            }
        }
        builder = builder.u64("point-gradient-norm-bits", point.gradient_norm_bits);
    }
    builder
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

fn record_call(
    calls: &RefCell<Vec<ObjectiveCall>>,
    phase: &'static str,
    objective: &'static str,
    point: &[f64],
    value: f64,
    gradient: &[f64],
) {
    calls.borrow_mut().push(ObjectiveCall {
        phase,
        objective,
        point_bits: bits(point),
        value_bits: value.to_bits(),
        gradient_bits: bits(gradient),
    });
}

fn weights() -> Vec<f64> {
    (1..10).map(|index| f64::from(index) / 10.0).collect()
}

fn epsilons() -> Vec<f64> {
    (0..8)
        .map(|index| 0.1f64.mul_add(f64::from(index), 0.15))
        .collect()
}

fn run_once() -> RunRecord {
    let calls = RefCell::new(Vec::new());
    let weighted_points = {
        let f1 = |x: &[f64]| {
            let value: f64 = x.iter().map(|coordinate| coordinate * coordinate).sum();
            let gradient: Vec<f64> = x.iter().map(|coordinate| 2.0 * coordinate).collect();
            record_call(&calls, "weighted", "f1", x, value, &gradient);
            (value, gradient)
        };
        let f2 = |x: &[f64]| {
            let value: f64 = x
                .iter()
                .map(|coordinate| (coordinate - 1.0) * (coordinate - 1.0))
                .sum();
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate - 1.0))
                .collect();
            record_call(&calls, "weighted", "f2", x, value, &gradient);
            (value, gradient)
        };
        weighted_sum_sweep(&f1, &f2, &weights(), &WEIGHTED_START)
            .iter()
            .map(PointPayload::from)
            .collect()
    };

    let epsilon_points = {
        let center = 1.0 / fs_math::det::sqrt(2.0);
        let f1 = |x: &[f64]| {
            let squared: f64 = x
                .iter()
                .map(|coordinate| (coordinate - center) * (coordinate - center))
                .sum();
            let exponential = fs_math::det::exp(-squared);
            let value = 1.0 - exponential;
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate - center) * exponential)
                .collect();
            record_call(&calls, "epsilon", "f1", x, value, &gradient);
            (value, gradient)
        };
        let f2 = |x: &[f64]| {
            let squared: f64 = x
                .iter()
                .map(|coordinate| (coordinate + center) * (coordinate + center))
                .sum();
            let exponential = fs_math::det::exp(-squared);
            let value = 1.0 - exponential;
            let gradient: Vec<f64> = x
                .iter()
                .map(|coordinate| 2.0 * (coordinate + center) * exponential)
                .collect();
            record_call(&calls, "epsilon", "f2", x, value, &gradient);
            (value, gradient)
        };
        epsilon_constraint_sweep(&f1, &f2, &epsilons(), &EPSILON_START, EPSILON_TOLERANCE)
            .iter()
            .map(PointPayload::from)
            .collect()
    };

    RunRecord {
        objective_calls: calls.into_inner(),
        weighted_points,
        epsilon_points,
    }
}

fn receipt(run: &RunRecord) -> RetainedReceipt {
    RetainedReceipt::new(ReceiptPayload {
        input_seed: INPUT_SEED,
        weights_bits: bits(&weights()),
        weighted_start_bits: bits(&WEIGHTED_START),
        epsilon_bits: bits(&epsilons()),
        epsilon_start_bits: bits(&EPSILON_START),
        objective_calls: run.objective_calls.clone(),
        weighted_points: run.weighted_points.clone(),
        epsilon_points: run.epsilon_points.clone(),
    })
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, usize, usize, u64) {
    let mut mutant = receipt.clone();
    let point = (MUTATION_SEED as usize) % mutant.payload.epsilon_points.len();
    let coordinate =
        ((MUTATION_SEED >> 8) as usize) % mutant.payload.epsilon_points[point].x_bits.len();
    let mask = 1_u64 << ((MUTATION_SEED >> 16) % 52);
    mutant.payload.epsilon_points[point].x_bits[coordinate] ^= mask;
    assert!(
        f64::from_bits(mutant.payload.epsilon_points[point].x_bits[coordinate]).is_finite(),
        "mantissa-only mutation must remain a finite wire-valid decision"
    );
    mutant.reseal();
    (mutant, point, coordinate, mask)
}

fn emit_receipt(
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    point: usize,
    coordinate: usize,
    mask: u64,
) {
    let json = format!(
        "{{\"input_seed\":{INPUT_SEED},\"mutation_seed\":{MUTATION_SEED},\
         \"reference_identity\":\"{}\",\"mutant_identity\":\"{}\",\
         \"mutated_path\":\"epsilon\",\"mutated_point\":{point},\
         \"mutated_coordinate\":{coordinate},\"mantissa_mask\":\"{mask:#018x}\",\
         \"merge_refusal\":\"reference-identity-mismatch\"}}",
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
    );
    let mut emitter = Emitter::new(SUITE, "two-objective-tracing");
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "pareto-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("Pareto study receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "two-objective-tracing".to_string(),
            pass: true,
            detail: format!(
                "fixed input seed {INPUT_SEED} replayed weighted and epsilon callback/result receipts; mutation seed {MUTATION_SEED:#018x} flipped epsilon point {point} coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and both merge gates refused it",
                mutant.declared_identity.hex(),
            ),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("Pareto seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line).expect("Pareto verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

fn assert_quality(run: &RunRecord) {
    let schedule = weights();
    assert_eq!(run.weighted_points.len(), schedule.len());
    let mut worst_closed_form = 0.0f64;
    for (point, weight) in run.weighted_points.iter().zip(schedule) {
        let f1 = f64::from_bits(point.objective_bits[0]);
        let f2 = f64::from_bits(point.objective_bits[1]);
        let expected_f1 = 3.0 * (1.0 - weight) * (1.0 - weight);
        let expected_f2 = 3.0 * weight * weight;
        worst_closed_form = worst_closed_form
            .max((f1 - expected_f1).abs())
            .max((f2 - expected_f2).abs());
        assert!(point.kkt_bits.is_none());
        assert!(f64::from_bits(point.gradient_norm_bits) < 1e-9);
    }
    assert!(worst_closed_form < 1e-7);

    assert_eq!(run.epsilon_points.len(), epsilons().len());
    let mut lowest_f1 = f64::INFINITY;
    let mut highest_f1 = f64::NEG_INFINITY;
    for point in &run.epsilon_points {
        let x0 = f64::from_bits(point.x_bits[0]);
        let x1 = f64::from_bits(point.x_bits[1]);
        assert!((x0 - x1).abs() < 1e-4);
        let kkt = point
            .kkt_bits
            .expect("epsilon path must retain KKT evidence");
        assert!(
            kkt.into_iter()
                .map(f64::from_bits)
                .all(|residual| residual.is_finite() && residual < 1e-5)
        );
        let f1 = f64::from_bits(point.objective_bits[0]);
        lowest_f1 = lowest_f1.min(f1);
        highest_f1 = highest_f1.max(f1);
    }
    assert!(highest_f1 - lowest_f1 > 0.6);
    for (phase, objective) in [
        ("weighted", "f1"),
        ("weighted", "f2"),
        ("epsilon", "f1"),
        ("epsilon", "f2"),
    ] {
        assert!(
            run.objective_calls
                .iter()
                .any(|call| call.phase == phase && call.objective == objective),
            "missing callback trace for {phase}/{objective}"
        );
    }
}

#[test]
fn pareto_tracing_replays_and_rejects_seeded_red_mutation() {
    let reference_run = run_once();
    assert_quality(&reference_run);
    let reference = receipt(&reference_run);
    admit_receipt(&reference.declared_identity, &reference)
        .expect("the internally consistent reference receipt must admit");

    let replay = receipt(&run_once());
    assert_eq!(
        replay, reference,
        "complete Pareto callback and result receipts did not replay"
    );

    let (mutant, point, coordinate, mask) = mutate_returned_decision(&reference);
    let (mutant_repeat, repeat_point, repeat_coordinate, repeat_mask) =
        mutate_returned_decision(&reference);
    assert_eq!(
        (point, coordinate, mask),
        (repeat_point, repeat_coordinate, repeat_mask)
    );
    assert_eq!(mutant, mutant_repeat, "seeded mutation was not stable");
    assert_ne!(mutant.declared_identity, reference.declared_identity);
    let mut stale_identity_mutant = mutant.clone();
    stale_identity_mutant.declared_identity = reference.declared_identity.clone();
    assert_eq!(
        admit_receipt(&reference.declared_identity, &stale_identity_mutant),
        Err(MergeRefusal::PayloadIdentityMismatch)
    );
    assert_eq!(
        admit_receipt(&reference.declared_identity, &mutant),
        Err(MergeRefusal::ReferenceIdentityMismatch)
    );

    emit_receipt(&reference, &mutant, point, coordinate, mask);
}
