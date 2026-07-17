//! G4 conformance for Problem-IR Study cancellation boundaries (7tv.21.21).
//!
//! Cancellation is orchestration state, not an optimizer stop reason. The
//! Study observes `fs_exec::Cx` only at replay-safe iteration boundaries,
//! returns a complete pause receipt without mutation, and resumes under a new
//! context. This battery injects repeated pre-requested gates after genuine
//! progress and requires bit-identical final state versus uninterrupted and
//! independently repeated storm schedules.

use fs_ascent::{
    STUDY_CANCELLATION_BOUNDARY_VERSION, StopReason, StopRule, Study, StudyPauseReceipt,
    StudyReport, StudyRunProgress,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder, Sense};
use fs_qty::Dims;

const SUITE: &str = "fs-ascent/runner-cancel";
const STORM_SEED: u64 = 0x5354_5544_59ca_1101;
const STORM_KERNEL_ID: u64 = 0x4153_4345_4e54_ca11;
const STORM_TILE: u64 = 0;
const FINAL_CONTEXT_ITERATION: u64 = 99;
const FD_H: f64 = 1e-6;
const LEARNING_RATE: f64 = 0.2;
const GRADIENT_TOLERANCE: f64 = 1e-9;
const FINAL_STEP_CAP: usize = 200;
const STORM_SEGMENTS: [usize; 4] = [1, 2, 3, 4];
const EXPECTED_PAUSE_STEPS: [usize; 5] = [0, 1, 3, 6, 10];
const D0: Dims = Dims([0, 0, 0, 0, 0, 0]);

fn with_cx<R>(cancelled: bool, iteration: u64, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: STORM_SEED,
                kernel_id: STORM_KERNEL_ID,
                tile: STORM_TILE,
                iteration,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    });
    assert!(
        pool.stats().quiescent(),
        "Study cancellation context must release its arena scope"
    );
    result
}

fn affine(builder: &mut ProblemBuilder, x: NodeId, slope: f64, offset: f64) -> NodeId {
    let slope = builder.konst(slope, D0).expect("finite slope");
    let scaled = builder.mul(slope, x).expect("scalar product");
    let offset = builder.konst(offset, D0).expect("finite offset");
    builder.add(scaled, offset).expect("scalar affine sum")
}

/// One-dimensional `(x - 1)^2`, started at `x = -3`.
fn quadratic_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let variable = builder
        .var("x", Manifold::Rn { dim: 1 }, D0)
        .expect("one-dimensional variable");
    let variable_ref = builder.var_ref(variable).expect("variable reference");
    let x = builder
        .component(variable_ref, 0)
        .expect("scalar component");
    let shifted = affine(&mut builder, x, 1.0, -1.0);
    let objective = builder.mul(shifted, shifted).expect("quadratic objective");
    builder
        .objective(objective, Sense::Minimize, 1.0)
        .expect("unit-weight objective");
    builder.finish()
}

fn stop_reason_name(reason: &StopReason) -> &'static str {
    match reason {
        StopReason::GradNorm => "grad-norm",
        StopReason::ObjectiveBelow => "objective-below",
        StopReason::Budget => "budget",
        StopReason::Stall => "stall",
        StopReason::Composite => "composite",
        StopReason::IterationCap => "iteration-cap",
    }
}

fn fixture_identity(problem: &Problem) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-cancellation-fixture-v1")
        .str("fs-ascent-version", fs_ascent::VERSION)
        .str("fs-exec-version", fs_exec::VERSION)
        .str(
            "problem-semantic-id",
            &problem
                .admit()
                .expect("fixture admits")
                .semantic_id()
                .to_hex(),
        )
        .f64_bits("initial-x", -3.0)
        .f64_bits("finite-difference-h", FD_H)
        .f64_bits("learning-rate", LEARNING_RATE)
        .f64_bits("gradient-tolerance", GRADIENT_TOLERANCE)
        .u64(
            "final-step-cap",
            u64::try_from(FINAL_STEP_CAP).expect("step cap fits u64"),
        )
        .u64("storm-seed", STORM_SEED)
        .u64("storm-kernel-id", STORM_KERNEL_ID)
        .u64("storm-tile", STORM_TILE)
        .u64("entry-cancellation-iteration", 0)
        .u64("final-context-iteration", FINAL_CONTEXT_ITERATION)
        .flag("deterministic-mode", true)
        .flag("infinite-context-budget", true)
        .u64(
            "cancellation-boundary-version",
            u64::from(STUDY_CANCELLATION_BOUNDARY_VERSION),
        );
    for (index, &steps) in STORM_SEGMENTS.iter().enumerate() {
        builder = builder.u64(
            "clean-segment-steps",
            u64::try_from(steps).expect("segment steps fit u64"),
        );
        builder = builder.u64(
            "segment-and-cancellation-iteration",
            u64::try_from(index).expect("storm index fits u64") + 1,
        );
    }
    builder.finish()
}

fn public_state_identity(study: &Study) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-public-cancel-state-v1")
        .str("problem-semantic-id", &study.problem_id().to_hex())
        .u64(
            "steps",
            u64::try_from(study.steps).expect("step count fits u64"),
        )
        .u64(
            "evaluations",
            u64::try_from(study.evals).expect("evaluation count fits u64"),
        )
        .u64(
            "point-values",
            u64::try_from(study.x.len()).expect("point length fits u64"),
        )
        .u64(
            "history-values",
            u64::try_from(study.history.len()).expect("history length fits u64"),
        );
    for &value in &study.x {
        builder = builder.f64_bits("point", value);
    }
    for &value in &study.history {
        builder = builder.f64_bits("objective-history", value);
    }
    builder.finish()
}

fn pause_identity(receipt: &StudyPauseReceipt) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-pause-receipt-v1")
        .u64("boundary-version", u64::from(receipt.boundary_version))
        .str("problem-semantic-id", &receipt.problem_id.to_hex())
        .u64(
            "steps",
            u64::try_from(receipt.steps).expect("step count fits u64"),
        )
        .u64(
            "evaluations",
            u64::try_from(receipt.evals).expect("evaluation count fits u64"),
        )
        .u64("finite-difference-h-bits", receipt.fd_h_bits)
        .u64("learning-rate-bits", receipt.learning_rate_bits)
        .u64(
            "point-values",
            u64::try_from(receipt.point_bits.len()).expect("point length fits u64"),
        )
        .u64(
            "history-values",
            u64::try_from(receipt.history_bits.len()).expect("history length fits u64"),
        )
        .flag(
            "current-objective-present",
            receipt.current_objective_bits.is_some(),
        )
        .flag(
            "current-gradient-present",
            receipt.current_gradient_norm_bits.is_some(),
        );
    if let Some(bits) = receipt.current_objective_bits {
        builder = builder.u64("current-objective-bits", bits);
    }
    if let Some(bits) = receipt.current_gradient_norm_bits {
        builder = builder.u64("current-gradient-bits", bits);
    }
    for &bits in &receipt.point_bits {
        builder = builder.u64("point-bits", bits);
    }
    for &bits in &receipt.history_bits {
        builder = builder.u64("objective-history-bits", bits);
    }
    builder.finish()
}

fn final_identity(fixture: &ReplayIdentity, study: &Study, report: &StudyReport) -> ReplayIdentity {
    IdentityBuilder::new("fs-ascent-study-cancel-final-v1")
        .child("fixture", fixture)
        .child("state", &public_state_identity(study))
        .str("stop-reason", stop_reason_name(&report.reason))
        .u64(
            "report-evaluations",
            u64::try_from(report.evals).expect("report evaluations fit u64"),
        )
        .f64_bits("report-objective", report.f)
        .f64_bits("report-gradient-norm", report.grad_norm)
        .finish()
}

fn storm_result_identity(
    fixture: &ReplayIdentity,
    final_state: &ReplayIdentity,
    pauses: &[StudyPauseReceipt],
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-cancel-storm-result-v1")
        .child("fixture", fixture)
        .child("final-state", final_state)
        .u64(
            "pause-receipts",
            u64::try_from(pauses.len()).expect("pause count fits u64"),
        );
    for receipt in pauses {
        builder = builder.child("pause", &pause_identity(receipt));
    }
    builder.finish()
}

fn assert_pause_matches(study: &Study, receipt: &StudyPauseReceipt) {
    assert_eq!(
        receipt.boundary_version,
        STUDY_CANCELLATION_BOUNDARY_VERSION
    );
    assert_eq!(receipt.problem_id, study.problem_id());
    assert_eq!(receipt.steps, study.steps);
    assert_eq!(receipt.evals, study.evals);
    assert_eq!(
        receipt.point_bits,
        study
            .x
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        receipt.history_bits,
        study
            .history
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
}

fn observe_cancel(
    study: &mut Study,
    problem: &Problem,
    rule: &StopRule,
    storm_index: u64,
) -> StudyPauseReceipt {
    let before = public_state_identity(study);
    let first = with_cx(true, storm_index, |cx| {
        study
            .try_run_cancellable(problem, rule, FINAL_STEP_CAP, cx)
            .expect("matching problem")
    });
    let StudyRunProgress::Paused(first) = first else {
        panic!("a pre-requested Study context must pause")
    };
    let after_first = public_state_identity(study);
    let second = with_cx(true, storm_index, |cx| {
        study
            .try_run_cancellable(problem, rule, FINAL_STEP_CAP, cx)
            .expect("matching problem")
    });
    let StudyRunProgress::Paused(second) = second else {
        panic!("repeated cancellation observation must remain paused")
    };
    let after_second = public_state_identity(study);
    assert_eq!(
        before.canonical_bytes(),
        after_first.canonical_bytes(),
        "first cancellation observation must not mutate public state"
    );
    assert_eq!(
        before.canonical_bytes(),
        after_second.canonical_bytes(),
        "repeated cancellation observation must not mutate public state"
    );
    assert_eq!(
        first, second,
        "private cache state must also stay unchanged"
    );
    assert_pause_matches(study, &first);
    first
}

fn run_storm(problem: &Problem) -> (Study, StudyReport, Vec<StudyPauseReceipt>) {
    let rule = StopRule::GradNorm(GRADIENT_TOLERANCE);
    let mut study = Study::new(problem, &[-3.0], FD_H, LEARNING_RATE);
    let mut pauses = Vec::with_capacity(STORM_SEGMENTS.len() + 1);

    // Entry cancellation has precedence and cannot even populate the objective
    // cache. The same checkpoint then proceeds normally under a fresh context.
    let entry_pause = observe_cancel(&mut study, problem, &rule, 0);
    assert_eq!(entry_pause.steps, 0);
    assert_eq!(entry_pause.evals, 0);
    assert!(entry_pause.current_objective_bits.is_none());
    pauses.push(entry_pause);

    for (index, &steps) in STORM_SEGMENTS.iter().enumerate() {
        let iteration = u64::try_from(index).expect("storm index fits u64") + 1;
        let progress = with_cx(false, iteration, |cx| {
            study
                .try_run_cancellable(problem, &rule, steps, cx)
                .expect("matching problem")
        });
        let StudyRunProgress::Stopped(segment) = progress else {
            panic!("a clean Study segment must not pause")
        };
        assert_eq!(
            segment.reason,
            StopReason::IterationCap,
            "the finite clean prefix must land before convergence"
        );
        let receipt = observe_cancel(&mut study, problem, &rule, iteration);
        assert!(
            receipt.steps > 0,
            "storm must follow genuine study progress"
        );
        pauses.push(receipt);
    }
    assert_eq!(
        pauses
            .iter()
            .map(|receipt| receipt.steps)
            .collect::<Vec<_>>(),
        EXPECTED_PAUSE_STEPS.to_vec(),
        "the G4 schedule must observe the intended completed boundaries"
    );

    let progress = with_cx(false, FINAL_CONTEXT_ITERATION, |cx| {
        study
            .try_run_cancellable(problem, &rule, FINAL_STEP_CAP, cx)
            .expect("matching problem")
    });
    let StudyRunProgress::Stopped(report) = progress else {
        panic!("fresh context must resume the checkpoint to an optimizer stop")
    };
    (study, report, pauses)
}

fn emit_receipt(result: &ReplayIdentity, reference: &ReplayIdentity, pauses: &[StudyPauseReceipt]) {
    let pause_roots = pauses
        .iter()
        .map(|receipt| format!("\"{}\"", pause_identity(receipt).hex()))
        .collect::<Vec<_>>()
        .join(",");
    let mut emitter = Emitter::new(SUITE, "g4-study-cancel-storm");
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "study-cancellation-storm-receipt".to_string(),
            json: format!(
                "{{\"identity\":\"{}\",\"reference_identity\":\"{}\",\
                 \"input_seed\":{STORM_SEED},\"pause_receipts\":[{pause_roots}],\
                 \"boundary_version\":{STUDY_CANCELLATION_BOUNDARY_VERSION}}}",
                result.hex(),
                reference.hex(),
            ),
        },
        None,
    );
    let receipt_line = receipt.to_jsonl();
    fs_obs::validate_line(&receipt_line)
        .expect("Study cancellation receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::StormAssertion {
            name: "study-pause-resume-bit-replay".to_string(),
            pass: true,
            seed: STORM_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict).expect("Study storm verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line)
        .expect("Study storm verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
fn g4_study_cancellation_storm_resumes_bit_exactly() {
    let problem = quadratic_problem();
    let fixture = fixture_identity(&problem);
    let rule = StopRule::GradNorm(GRADIENT_TOLERANCE);

    let mut reference_study = Study::new(&problem, &[-3.0], FD_H, LEARNING_RATE);
    let reference_report = reference_study.run(&problem, &rule, FINAL_STEP_CAP);
    let reference_identity = final_identity(&fixture, &reference_study, &reference_report);

    let (storm_study, storm_report, pauses) = run_storm(&problem);
    let storm_final_identity = final_identity(&fixture, &storm_study, &storm_report);
    let storm_identity = storm_result_identity(&fixture, &storm_final_identity, &pauses);
    let (repeat_study, repeat_report, repeat_pauses) = run_storm(&problem);
    let repeat_final_identity = final_identity(&fixture, &repeat_study, &repeat_report);
    let repeat_identity = storm_result_identity(&fixture, &repeat_final_identity, &repeat_pauses);

    assert_eq!(
        storm_final_identity.canonical_bytes(),
        reference_identity.canonical_bytes(),
        "repeated cancellation boundaries must be invisible to the final Study trajectory"
    );
    assert_eq!(
        storm_identity.canonical_bytes(),
        repeat_identity.canonical_bytes(),
        "the complete cancellation schedule must replay bit for bit"
    );
    assert_eq!(
        pauses, repeat_pauses,
        "every pause receipt must replay bit for bit"
    );
    assert_eq!(storm_report.reason, reference_report.reason);
    assert_ne!(
        storm_report.reason,
        StopReason::IterationCap,
        "the final resumed run must reach the optimizer's own stop rule"
    );

    emit_receipt(&storm_identity, &reference_identity, &pauses);
}
