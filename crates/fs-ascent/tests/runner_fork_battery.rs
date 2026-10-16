//! G5/P9 conformance for gradient-study world forks (bead 7tv.21.15).
//!
//! A `Study` checkpoint is semantically bound to the admitted `Problem` whose
//! objective and projected-gradient caches it contains. Steering is therefore
//! explicit: `fork_for` copies the branch point and driver configuration into a
//! fresh child, clears every problem-dependent cache/counter, records both
//! problem identities, and leaves the parent untouched. This battery proves
//! opposite objective reweightings diverge, each sibling replays bitwise, and
//! accidental continuation or coordinate reinterpretation refuses before
//! mutating public state. It does not claim persisted ledger replay, `Cx`
//! cancellation recovery, worker-count replay, cross-ISA equality, or runtime
//! performance.

use fs_ascent::{StopReason, StopRule, Study, StudyError, StudyForkReceipt, StudyReport};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder, Sense};
use fs_qty::Dims;

const SUITE: &str = "fs-ascent/runner-fork";
const INPUT_SEED: u64 = 0;
const FD_H: f64 = 1e-6;
const LEARNING_RATE: f64 = 0.2;
const GRADIENT_TOLERANCE: f64 = 1e-9;
const MAX_BRANCH_STEPS: usize = 200;
const D0: Dims = Dims([0, 0, 0, 0, 0, 0]);

fn affine(builder: &mut ProblemBuilder, x: NodeId, slope: f64, offset: f64) -> NodeId {
    let slope = builder.konst(slope, D0).expect("finite slope");
    let scaled = builder.mul(slope, x).expect("scalar product");
    let offset = builder.konst(offset, D0).expect("finite offset");
    builder.add(scaled, offset).expect("scalar affine sum")
}

/// Two objectives on one shared coordinate:
///
/// `left_weight * (x + 1)^2 + right_weight * (x - 1)^2`.
///
/// Equal weights minimize at zero; 0.95/0.05 minimizes at -0.9 and the
/// reversed weights minimize at +0.9.
fn weighted_problem(variable_name: &str, left_weight: f64, right_weight: f64) -> Problem {
    let mut builder = ProblemBuilder::new();
    let variable = builder
        .var(variable_name, Manifold::Rn { dim: 1 }, D0)
        .expect("one-dimensional variable");
    let variable_ref = builder.var_ref(variable).expect("variable reference");
    let x = builder
        .component(variable_ref, 0)
        .expect("scalar component");
    let left_offset = affine(&mut builder, x, 1.0, 1.0);
    let right_offset = affine(&mut builder, x, 1.0, -1.0);
    let left = builder.mul(left_offset, left_offset).expect("left square");
    let right = builder
        .mul(right_offset, right_offset)
        .expect("right square");
    builder
        .objective(left, Sense::Minimize, left_weight)
        .expect("left objective");
    builder
        .objective(right, Sense::Minimize, right_weight)
        .expect("right objective");
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

fn state_identity(study: &Study) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-public-state-v1")
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

fn branch_identity(study: &Study, report: &StudyReport) -> ReplayIdentity {
    IdentityBuilder::new("fs-ascent-study-fork-branch-v1")
        .child("state", &state_identity(study))
        .str("stop-reason", stop_reason_name(&report.reason))
        .u64(
            "report-evaluations",
            u64::try_from(report.evals).expect("report evaluation count fits u64"),
        )
        .f64_bits("report-objective", report.f)
        .f64_bits("report-gradient-norm", report.grad_norm)
        .finish()
}

fn receipt_identity(receipt: &StudyForkReceipt) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-study-fork-receipt-v1")
        .str(
            "parent-problem-semantic-id",
            &receipt.parent_problem_id.to_hex(),
        )
        .str(
            "child-problem-semantic-id",
            &receipt.child_problem_id.to_hex(),
        )
        .u64(
            "parent-steps",
            u64::try_from(receipt.parent_steps).expect("parent step count fits u64"),
        )
        .u64(
            "parent-evaluations",
            u64::try_from(receipt.parent_evals).expect("parent evaluation count fits u64"),
        )
        .u64(
            "parent-history-values",
            u64::try_from(receipt.parent_history_len).expect("parent history length fits u64"),
        )
        .u64("finite-difference-step-bits", receipt.fd_h_bits)
        .u64("learning-rate-bits", receipt.learning_rate_bits)
        .u64(
            "branch-point-values",
            u64::try_from(receipt.branch_point_bits.len()).expect("branch point length fits u64"),
        );
    for &bits in &receipt.branch_point_bits {
        builder = builder.u64("branch-point-bits", bits);
    }
    builder.finish()
}

fn emit_receipt(
    parent: &ReplayIdentity,
    left_receipt: &ReplayIdentity,
    right_receipt: &ReplayIdentity,
    left_branch: &ReplayIdentity,
    right_branch: &ReplayIdentity,
) {
    let json = format!(
        "{{\"input_seed\":{INPUT_SEED},\"parent_state\":\"{}\",\
         \"left_fork\":\"{}\",\"right_fork\":\"{}\",\
         \"left_branch\":\"{}\",\"right_branch\":\"{}\"}}",
        parent.hex(),
        left_receipt.hex(),
        right_receipt.hex(),
        left_branch.hex(),
        right_branch.hex(),
    );
    let mut emitter = Emitter::new(SUITE, "gradient-study-world-fork");
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "gradient-study-world-fork-receipt".to_string(),
            json,
        },
        None,
    );
    let receipt_line = receipt.to_jsonl();
    fs_obs::validate_line(&receipt_line).expect("fork receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: "gradient-study-world-fork".to_string(),
            pass: true,
            detail: "semantic binding refused stale-cache continuation; opposite reweight forks left the parent immutable, cleared branch-local state, diverged, and replayed bitwise"
                .to_string(),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict).expect("fork verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line).expect("fork verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
#[allow(clippy::too_many_lines)] // one end-to-end world-fork receipt and its refusal mutants
fn gradient_study_world_forks_are_independent_and_replayable() {
    let parent_problem = weighted_problem("x", 0.5, 0.5);
    let left_problem = weighted_problem("x", 0.95, 0.05);
    let right_problem = weighted_problem("x", 0.05, 0.95);
    let renamed_problem = weighted_problem("renamed-x", 0.95, 0.05);
    let rule = StopRule::GradNorm(GRADIENT_TOLERANCE);

    let mut parent =
        Study::try_new(&parent_problem, &[0.6], FD_H, LEARNING_RATE).expect("parent study admits");
    let trunk_report = parent
        .try_run(&parent_problem, &rule, 2)
        .expect("parent problem matches its checkpoint");
    assert_eq!(trunk_report.reason, StopReason::IterationCap);
    assert_eq!(parent.steps, 2, "fixture must hold a live cached objective");
    assert_eq!(parent.history.len(), 2);
    let parent_at_fork = state_identity(&parent);

    // The pre-fix hazard: a clone carried the parent's current_f cache into a
    // differently weighted problem. Semantic binding now refuses before any
    // public state moves and points callers to the explicit fork operation.
    let mut accidental_continuation = parent.clone();
    let accidental_before = state_identity(&accidental_continuation);
    let mismatch = accidental_continuation
        .try_run(&left_problem, &rule, 1)
        .expect_err("a checkpoint cannot silently change semantic problems");
    assert!(matches!(
        mismatch,
        StudyError::ProblemMismatch {
            bound,
            supplied
        } if bound == parent.problem_id() && supplied == left_problem.admit().expect("left admits").semantic_id()
    ));
    assert_eq!(
        state_identity(&accidental_continuation).canonical_bytes(),
        accidental_before.canonical_bytes(),
        "typed mismatch refusal is mutation-free"
    );
    assert!(matches!(
        parent
            .fork_for(&parent_problem)
            .expect_err("an unchanged problem cannot reset branch-local accounting"),
        StudyError::ForkProblemUnchanged { problem_id }
            if problem_id == parent.problem_id()
    ));

    let (mut left, left_fork) = parent
        .fork_for(&left_problem)
        .expect("left reweight preserves the variable schema");
    let (mut right, right_fork) = parent
        .fork_for(&right_problem)
        .expect("right reweight preserves the variable schema");
    assert_eq!(
        state_identity(&parent).canonical_bytes(),
        parent_at_fork.canonical_bytes(),
        "world-forking borrows and cannot mutate the parent"
    );
    for child in [&left, &right] {
        assert!(
            child
                .x
                .iter()
                .map(|value| value.to_bits())
                .eq(parent.x.iter().map(|value| value.to_bits())),
            "child starts at the exact branch point"
        );
        assert!(
            child.history.is_empty(),
            "old objective values are invalid after reweighting"
        );
        assert_eq!(child.steps, 0, "steps are branch-local");
        assert_eq!(child.evals, 0, "evaluations are branch-local");
    }
    assert_eq!(
        left_fork.branch_point_bits,
        parent
            .x
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(left_fork.parent_problem_id, parent.problem_id());
    assert_eq!(right_fork.parent_problem_id, parent.problem_id());
    assert_ne!(left_fork.child_problem_id, right_fork.child_problem_id);

    let left_report = left
        .try_run(&left_problem, &rule, MAX_BRANCH_STEPS)
        .expect("left branch remains bound");
    let right_report = right
        .try_run(&right_problem, &rule, MAX_BRANCH_STEPS)
        .expect("right branch remains bound");
    assert_eq!(left_report.reason, StopReason::GradNorm);
    assert_eq!(right_report.reason, StopReason::GradNorm);
    assert!(left.x[0] < -0.899_999, "left weight must steer toward -0.9");
    assert!(
        right.x[0] > 0.899_999,
        "right weight must steer toward +0.9"
    );

    let left_branch = branch_identity(&left, &left_report);
    let right_branch = branch_identity(&right, &right_report);
    assert_ne!(
        left_branch.canonical_bytes(),
        right_branch.canonical_bytes(),
        "opposite steering must produce genuinely different worlds"
    );

    let (mut left_repeat, left_fork_repeat) = parent
        .fork_for(&left_problem)
        .expect("repeat left fork admits");
    let (mut right_repeat, right_fork_repeat) = parent
        .fork_for(&right_problem)
        .expect("repeat right fork admits");
    let left_repeat_report = left_repeat
        .try_run(&left_problem, &rule, MAX_BRANCH_STEPS)
        .expect("repeat left branch runs");
    let right_repeat_report = right_repeat
        .try_run(&right_problem, &rule, MAX_BRANCH_STEPS)
        .expect("repeat right branch runs");
    assert_eq!(
        left_fork, left_fork_repeat,
        "left fork receipt replays exactly"
    );
    assert_eq!(
        right_fork, right_fork_repeat,
        "right fork receipt replays exactly"
    );
    assert_eq!(
        branch_identity(&left_repeat, &left_repeat_report).canonical_bytes(),
        left_branch.canonical_bytes(),
        "left branch replays bitwise"
    );
    assert_eq!(
        branch_identity(&right_repeat, &right_repeat_report).canonical_bytes(),
        right_branch.canonical_bytes(),
        "right branch replays bitwise"
    );

    let schema_error = parent
        .fork_for(&renamed_problem)
        .expect_err("a fork cannot reinterpret packed coordinates");
    assert!(matches!(
        schema_error,
        StudyError::ForkVariableSchemaMismatch {
            first_mismatch: 0,
            parent_variables: 1,
            child_variables: 1,
        }
    ));
    assert_eq!(
        state_identity(&parent).canonical_bytes(),
        parent_at_fork.canonical_bytes(),
        "all successful and refused forks leave the parent immutable"
    );

    emit_receipt(
        &parent_at_fork,
        &receipt_identity(&left_fork),
        &receipt_identity(&right_fork),
        &left_branch,
        &right_branch,
    );
}
