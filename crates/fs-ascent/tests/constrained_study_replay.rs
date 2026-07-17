//! G5 replay and seeded-failure self-tests for all constrained ASCENT engines.
//!
//! Augmented Lagrangian, log-barrier interior point, and active-set SQP solve
//! the same analytic equality-plus-active-inequality fixture. Per engine, a
//! retained receipt binds the exact ordered objective, constraint, and
//! Jacobian-transpose callback trace together with every public report and KKT
//! field. Same-input repeats must reproduce the receipt byte for byte. An
//! engine-keyed deterministic red mutation flips one finite returned-decision
//! mantissa bit and is refused both before and after self-consistent resealing.
//!
//! This is one small dense fixture. It does not claim all constrained problems,
//! large-scale sparse behavior, cancellation, checkpointing, cross-ISA
//! equality, ledger persistence, or performance.

use core::cell::RefCell;

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{KktResidual, augmented_lagrangian, interior_point, sqp};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity, check_version};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-ascent/constrained-study-replay";
const INPUT_SEED: u64 = 0;
const MUTATION_SEED: u64 = 0x434F_4E53_5452_5244;
const START: [f64; 2] = [0.0, 0.0];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Engine {
    AugmentedLagrangian,
    InteriorPoint,
    Sqp,
}

impl Engine {
    const ALL: [Self; 3] = [Self::AugmentedLagrangian, Self::InteriorPoint, Self::Sqp];

    const fn name(self) -> &'static str {
        match self {
            Self::AugmentedLagrangian => "augmented-lagrangian",
            Self::InteriorPoint => "interior-point",
            Self::Sqp => "active-set-sqp",
        }
    }

    const fn tag(self) -> u64 {
        match self {
            Self::AugmentedLagrangian => 0x414C,
            Self::InteriorPoint => 0x4950,
            Self::Sqp => 0x5351,
        }
    }

    const fn tolerance(self) -> f64 {
        match self {
            Self::InteriorPoint => 1e-6,
            Self::AugmentedLagrangian | Self::Sqp => 1e-7,
        }
    }

    const fn iteration_cap(self) -> usize {
        match self {
            Self::AugmentedLagrangian => 40,
            Self::InteriorPoint | Self::Sqp => 60,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CallbackCall {
    kind: &'static str,
    point_bits: Vec<u64>,
    weight_bits: Vec<u64>,
    scalar_bits: Option<u64>,
    output_bits: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReportPayload {
    x_bits: Vec<u64>,
    f_bits: u64,
    kkt_bits: [u64; 4],
    lambda_bits: Vec<u64>,
    nu_bits: Vec<u64>,
    iterations: usize,
    evaluations: usize,
    converged: bool,
}

impl ReportPayload {
    #[allow(clippy::too_many_arguments)] // Mirrors the closed public report surface exactly.
    fn from_parts(
        x: &[f64],
        f: f64,
        kkt: &KktResidual,
        lambda: &[f64],
        nu: &[f64],
        iterations: usize,
        evaluations: usize,
        converged: bool,
    ) -> Self {
        Self {
            x_bits: bits(x),
            f_bits: f.to_bits(),
            kkt_bits: [
                kkt.stationarity.to_bits(),
                kkt.feasibility.to_bits(),
                kkt.dual_feasibility.to_bits(),
                kkt.complementarity.to_bits(),
            ],
            lambda_bits: bits(lambda),
            nu_bits: bits(nu),
            iterations,
            evaluations,
            converged,
        }
    }
}

#[derive(Debug)]
struct RunRecord {
    callbacks: Vec<CallbackCall>,
    report: ReportPayload,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReceiptPayload {
    engine: Engine,
    input_seed: u64,
    start_bits: Vec<u64>,
    callbacks: Vec<CallbackCall>,
    report: ReportPayload,
}

impl ReceiptPayload {
    fn identity(&self) -> ReplayIdentity {
        let mut builder = IdentityBuilder::new("fs-ascent-constrained-study-receipt-v1")
            .str("fs-ascent-version", fs_ascent::VERSION)
            .str("engine", self.engine.name())
            .u64("engine-tag", self.engine.tag())
            .u64("input-seed", self.input_seed)
            .f64_bits("kkt-tolerance", self.engine.tolerance())
            .u64("iteration-cap", self.engine.iteration_cap() as u64)
            .u64("start-values", self.start_bits.len() as u64);
        for &value_bits in &self.start_bits {
            builder = builder.u64("start-value-bits", value_bits);
        }

        builder = builder.u64("callback-calls", self.callbacks.len() as u64);
        for call in &self.callbacks {
            builder = builder
                .str("callback-kind", call.kind)
                .u64("callback-point-values", call.point_bits.len() as u64);
            for &value_bits in &call.point_bits {
                builder = builder.u64("callback-point-bits", value_bits);
            }
            builder = builder.u64("callback-weight-values", call.weight_bits.len() as u64);
            for &value_bits in &call.weight_bits {
                builder = builder.u64("callback-weight-bits", value_bits);
            }
            builder = builder.flag("callback-has-scalar", call.scalar_bits.is_some());
            if let Some(scalar_bits) = call.scalar_bits {
                builder = builder.u64("callback-scalar-bits", scalar_bits);
            }
            builder = builder.u64("callback-output-values", call.output_bits.len() as u64);
            for &value_bits in &call.output_bits {
                builder = builder.u64("callback-output-bits", value_bits);
            }
        }

        builder = builder.u64("report-x-values", self.report.x_bits.len() as u64);
        for &value_bits in &self.report.x_bits {
            builder = builder.u64("report-x-bits", value_bits);
        }
        builder = builder
            .u64("report-objective-bits", self.report.f_bits)
            .u64("report-kkt-values", self.report.kkt_bits.len() as u64);
        for &value_bits in &self.report.kkt_bits {
            builder = builder.u64("report-kkt-bits", value_bits);
        }
        builder = builder.u64("report-lambda-values", self.report.lambda_bits.len() as u64);
        for &value_bits in &self.report.lambda_bits {
            builder = builder.u64("report-lambda-bits", value_bits);
        }
        builder = builder.u64("report-nu-values", self.report.nu_bits.len() as u64);
        for &value_bits in &self.report.nu_bits {
            builder = builder.u64("report-nu-bits", value_bits);
        }
        builder
            .u64("report-iterations", self.report.iterations as u64)
            .u64("report-evaluations", self.report.evaluations as u64)
            .flag("report-converged", self.report.converged)
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

fn push_call(
    calls: &RefCell<Vec<CallbackCall>>,
    kind: &'static str,
    point: &[f64],
    weights: &[f64],
    scalar: Option<f64>,
    output: &[f64],
) {
    calls.borrow_mut().push(CallbackCall {
        kind,
        point_bits: bits(point),
        weight_bits: bits(weights),
        scalar_bits: scalar.map(f64::to_bits),
        output_bits: bits(output),
    });
}

fn run_once(engine: Engine) -> RunRecord {
    let calls = RefCell::new(Vec::new());
    let report = {
        let mut objective = |x: &[f64]| {
            let dx = x[0] - 2.0;
            let dy = x[1] - 1.0;
            let value = dx * dx + dy * dy;
            let gradient = vec![2.0 * dx, 2.0 * dy];
            push_call(&calls, "objective-gradient", x, &[], Some(value), &gradient);
            (value, gradient)
        };
        let equality = |x: &[f64]| {
            let value = vec![x[0] + x[1] - 2.0];
            push_call(&calls, "equality", x, &[], None, &value);
            value
        };
        let equality_jt = |x: &[f64], weights: &[f64]| {
            let value = vec![weights[0], weights[0]];
            push_call(&calls, "equality-jt", x, weights, None, &value);
            value
        };
        let inequality = |x: &[f64]| {
            let value = vec![x[0] - 1.2];
            push_call(&calls, "inequality", x, &[], None, &value);
            value
        };
        let inequality_jt = |x: &[f64], weights: &[f64]| {
            let value = vec![weights[0], 0.0];
            push_call(&calls, "inequality-jt", x, weights, None, &value);
            value
        };
        let mut problem = ConstrainedProblem {
            fg: &mut objective,
            ce: &equality,
            ce_jt: &equality_jt,
            ci: &inequality,
            ci_jt: &inequality_jt,
        };
        match engine {
            Engine::AugmentedLagrangian => {
                let report = augmented_lagrangian(
                    &mut problem,
                    &START,
                    engine.tolerance(),
                    engine.iteration_cap(),
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.outer_iters,
                    report.evals,
                    report.converged,
                )
            }
            Engine::InteriorPoint => {
                let report = interior_point(
                    &mut problem,
                    &START,
                    engine.tolerance(),
                    engine.iteration_cap(),
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.outer_iters,
                    report.evals,
                    report.converged,
                )
            }
            Engine::Sqp => {
                let report = sqp(
                    &mut problem,
                    &START,
                    engine.tolerance(),
                    engine.iteration_cap(),
                );
                ReportPayload::from_parts(
                    &report.x,
                    report.f,
                    &report.kkt,
                    &report.lambda,
                    &report.nu,
                    report.iters,
                    report.evals,
                    report.converged,
                )
            }
        }
    };
    RunRecord {
        callbacks: calls.into_inner(),
        report,
    }
}

fn receipt(engine: Engine, run: &RunRecord) -> RetainedReceipt {
    RetainedReceipt::new(ReceiptPayload {
        engine,
        input_seed: INPUT_SEED,
        start_bits: bits(&START),
        callbacks: run.callbacks.clone(),
        report: run.report.clone(),
    })
}

fn mutate_returned_decision(receipt: &RetainedReceipt) -> (RetainedReceipt, u64, usize, u64) {
    let mutation_seed = MUTATION_SEED ^ receipt.payload.engine.tag();
    let mut mutant = receipt.clone();
    let coordinate = (mutation_seed as usize) % mutant.payload.report.x_bits.len();
    let mask = 1_u64 << ((mutation_seed >> 8) % 52);
    mutant.payload.report.x_bits[coordinate] ^= mask;
    assert!(
        f64::from_bits(mutant.payload.report.x_bits[coordinate]).is_finite(),
        "mantissa-only mutation must remain a finite wire-valid decision"
    );
    mutant.reseal();
    (mutant, mutation_seed, coordinate, mask)
}

fn emit_receipt(
    engine: Engine,
    reference: &RetainedReceipt,
    mutant: &RetainedReceipt,
    mutation_seed: u64,
    coordinate: usize,
    mask: u64,
) {
    let json = format!(
        "{{\"engine\":\"{}\",\"input_seed\":{INPUT_SEED},\
         \"mutation_seed\":{mutation_seed},\"reference_identity\":\"{}\",\
         \"mutant_identity\":\"{}\",\"mutated_coordinate\":{coordinate},\
         \"mantissa_mask\":\"{mask:#018x}\",\
         \"merge_refusal\":\"reference-identity-mismatch\"}}",
        engine.name(),
        reference.declared_identity.hex(),
        mutant.declared_identity.hex(),
    );
    let mut emitter = Emitter::new(SUITE, engine.name());
    let receipt_event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "constrained-study-replay-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt_event.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("constrained study receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: engine.name().to_string(),
            pass: true,
            detail: format!(
                "fixed input seed {INPUT_SEED} replayed the complete {} callback/report receipt; mutation seed {mutation_seed:#018x} flipped coordinate {coordinate} mask {mask:#018x}, produced stable identity {}, and both merge gates refused it",
                engine.name(),
                mutant.declared_identity.hex(),
            ),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict)
        .expect("constrained seeded-failure verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("constrained verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

fn assert_quality(engine: Engine, run: &RunRecord) {
    assert!(run.report.converged, "{} did not converge", engine.name());
    let x: Vec<f64> = run
        .report
        .x_bits
        .iter()
        .map(|&value_bits| f64::from_bits(value_bits))
        .collect();
    assert!(
        (x[0] - 1.2).abs() < 1e-4 && (x[1] - 0.8).abs() < 1e-4,
        "{} missed analytic optimum: {x:?}",
        engine.name(),
    );
    assert!(
        run.report
            .kkt_bits
            .iter()
            .map(|&value_bits| f64::from_bits(value_bits))
            .all(|residual| residual.is_finite() && residual < engine.tolerance()),
        "{} returned a non-certifying KKT receipt",
        engine.name(),
    );
    assert!(
        run.report
            .nu_bits
            .first()
            .is_some_and(|&value_bits| f64::from_bits(value_bits) > 0.0),
        "{} must retain a positive active multiplier",
        engine.name(),
    );
    assert!(run.report.iterations > 0 && run.report.evaluations > 0);
    for required_kind in [
        "objective-gradient",
        "equality",
        "equality-jt",
        "inequality",
        "inequality-jt",
    ] {
        assert!(
            run.callbacks.iter().any(|call| call.kind == required_kind),
            "{} did not exercise callback {required_kind}",
            engine.name(),
        );
    }
}

#[test]
fn constrained_families_replay_and_reject_seeded_red_mutations() {
    for engine in Engine::ALL {
        let reference_run = run_once(engine);
        assert_quality(engine, &reference_run);
        let reference = receipt(engine, &reference_run);
        admit_receipt(&reference.declared_identity, &reference)
            .expect("the internally consistent reference receipt must admit");

        let replay = receipt(engine, &run_once(engine));
        assert_eq!(
            replay,
            reference,
            "{} callback trace and report did not replay",
            engine.name(),
        );

        let (mutant, mutation_seed, coordinate, mask) = mutate_returned_decision(&reference);
        let (mutant_repeat, repeat_seed, repeat_coordinate, repeat_mask) =
            mutate_returned_decision(&reference);
        assert_eq!(
            (mutation_seed, coordinate, mask),
            (repeat_seed, repeat_coordinate, repeat_mask)
        );
        assert_eq!(
            mutant,
            mutant_repeat,
            "{} red mutation and identity were unstable",
            engine.name(),
        );
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

        emit_receipt(engine, &reference, &mutant, mutation_seed, coordinate, mask);
    }
}
