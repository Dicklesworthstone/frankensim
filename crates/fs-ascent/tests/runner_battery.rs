//! Bead ijil (d): the Problem-IR study runner — manifold-product
//! variable packing, budget threading through the stop algebra,
//! bitwise-resumable studies, and constrained routing through the
//! packed adapters. Aggregate verdicts use the canonical fs-obs schema;
//! the G5 case additionally retains the admitted problem identity and
//! exact public replay state.

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{Packing, StopReason, StopRule, Study, StudyReport, augmented_lagrangian};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_opt::{ConstraintKind, EvalLimit, Manifold, NodeId, ProblemBuilder, Sense};
use fs_qty::Dims;
use std::fmt::Write as _;
use std::num::NonZeroU64;

const SUITE: &str = "fs-ascent/runner";
const FIXED_INPUT_SEED: u64 = 0;
const REPLAY_FD_H: f64 = 1e-6;
const REPLAY_LR: f64 = 0.2;
const REPLAY_GRAD_TOL: f64 = 1e-9;
const REPLAY_STEPS: usize = 300;
const REPLAY_CUTS: [usize; 3] = [1, 17, 64];

fn emit_verdict(emitter: &mut Emitter, name: &str, pass: bool, details: &str, seed: u64) {
    let event = emitter.emit(
        if pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: name.to_string(),
            pass,
            detail: details.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("runner verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("runner verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "{name}: {details}");
}

fn verdict(name: &str, pass: bool, details: &str) {
    let mut emitter = Emitter::new(SUITE, name);
    emit_verdict(&mut emitter, name, pass, details, FIXED_INPUT_SEED);
}

fn replay_verdict(name: &str, pass: bool, details: &str, receipt_json: String) {
    let mut emitter = Emitter::new(SUITE, name);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "problem-ir-study-replay".to_string(),
            json: receipt_json,
        },
        None,
    );
    let line = receipt.to_jsonl();
    fs_obs::validate_line(&line).expect("study replay receipt must use the fs-obs wire schema");
    println!("{line}");
    emit_verdict(&mut emitter, name, pass, details, FIXED_INPUT_SEED);
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

fn study_state_identity(
    fixture: &ReplayIdentity,
    study: &Study,
    report: &StudyReport,
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-ascent-problem-ir-study-state-v1")
        .child("fixture", fixture)
        .u64(
            "steps",
            u64::try_from(study.steps).expect("study steps fit u64"),
        )
        .u64(
            "study-evals",
            u64::try_from(study.evals).expect("study evals fit u64"),
        )
        .u64(
            "report-evals",
            u64::try_from(report.evals).expect("report evals fit u64"),
        )
        .str("stop-reason", stop_reason_name(&report.reason))
        .f64_bits("report-objective", report.f)
        .f64_bits("report-gradient-norm", report.grad_norm)
        .u64(
            "point-values",
            u64::try_from(study.x.len()).expect("study point length fits u64"),
        )
        .u64(
            "history-values",
            u64::try_from(study.history.len()).expect("study history length fits u64"),
        );
    for &value in &study.x {
        builder = builder.f64_bits("point", value);
    }
    for &value in &study.history {
        builder = builder.f64_bits("objective-history", value);
    }
    builder.finish()
}

fn expected_study_evals(study: &Study) -> usize {
    1 + 2 * study.x.len() * study.history.len() + study.steps
}

fn bit_vector_json(values: &[f64]) -> String {
    let mut json = String::from("[");
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            json.push(',');
        }
        write!(&mut json, "\"0x{:016x}\"", value.to_bits()).expect("String writes are infallible");
    }
    json.push(']');
    json
}

fn identity_vector_json(identities: &[ReplayIdentity]) -> String {
    let mut json = String::from("[");
    for (index, identity) in identities.iter().enumerate() {
        if index > 0 {
            json.push(',');
        }
        write!(&mut json, "\"{}\"", identity.hex()).expect("String writes are infallible");
    }
    json.push(']');
    json
}

fn usize_vector_json(values: &[usize]) -> String {
    let mut json = String::from("[");
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            json.push(',');
        }
        write!(&mut json, "{value}").expect("String writes are infallible");
    }
    json.push(']');
    json
}

const D0: Dims = Dims([0, 0, 0, 0, 0, 0]);

/// (a·x + b) as a scalar node from a component.
fn affine(b: &mut ProblemBuilder, x: NodeId, a: f64, off: f64) -> NodeId {
    let ca = b.konst(a, D0).expect("finite konst");
    let m = b.mul(ca, x).expect("scalar mul");
    let co = b.konst(off, D0).expect("finite konst");
    b.add(m, co).expect("scalar add")
}

/// A product-manifold problem: v on the unit sphere S², z in R²;
/// minimize 0.6·v₀ + 0.8·v₁ + (z₀−1)² + (z₁+2)² — linear on the
/// sphere (optimum v = −(0.6, 0.8, 0)) plus a shifted bowl.
fn product_problem() -> (fs_opt::Problem, Vec<f64>) {
    product_problem_budgeted(EvalLimit::Unlimited)
}

/// [`product_problem`] with a P4 budget attached at BUILD time — the
/// sealed `Problem` no longer exposes a mutable budget field.
fn product_problem_budgeted(eval_limit: EvalLimit) -> (fs_opt::Problem, Vec<f64>) {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("v", Manifold::Sphere { ambient: 3 }, D0)
        .expect("var v");
    let z = b.var("z", Manifold::Rn { dim: 2 }, D0).expect("var z");
    let vref = b.var_ref(v).expect("v node");
    let zref = b.var_ref(z).expect("z node");
    let v0 = b.component(vref, 0).expect("v0");
    let v1 = b.component(vref, 1).expect("v1");
    let z0 = b.component(zref, 0).expect("z0");
    let z1 = b.component(zref, 1).expect("z1");
    let lin0 = affine(&mut b, v0, 0.6, 0.0);
    let lin1 = affine(&mut b, v1, 0.8, 0.0);
    let lin = b.add(lin0, lin1).expect("lin");
    let dz0 = affine(&mut b, z0, 1.0, -1.0);
    let dz1 = affine(&mut b, z1, 1.0, 2.0);
    let q0 = b.mul(dz0, dz0).expect("q0");
    let q1 = b.mul(dz1, dz1).expect("q1");
    let bowl = b.add(q0, q1).expect("bowl");
    let total = b.add(lin, bowl).expect("total");
    b.objective(total, Sense::Minimize, 1.0).expect("objective");
    b.set_eval_limit(eval_limit);
    let problem = b.finish();
    let x0 = vec![1.0, 0.0, 0.0, 0.0, 0.0];
    (problem, x0)
}

#[test]
fn runner_product_manifold_packing() {
    let (problem, x0) = product_problem();
    let mut study = Study::new(&problem, &x0, 1e-6, 0.2);
    let rep = study.run(&problem, &StopRule::GradNorm(1e-6), 4000);
    // Sphere block stays unit; optimum v = -(0.6, 0.8, 0), z = (1, -2).
    let vn = fs_math::det::sqrt(study.x[..3].iter().map(|v| v * v).sum::<f64>());
    let ok_v = (study.x[0] + 0.6).abs() < 1e-3 && (study.x[1] + 0.8).abs() < 1e-3;
    let ok_z = (study.x[3] - 1.0).abs() < 1e-3 && (study.x[4] + 2.0).abs() < 1e-3;
    verdict(
        "ijil-runner-product",
        (vn - 1.0).abs() < 1e-12 && ok_v && ok_z && rep.f < -1.0 + 1e-3,
        &format!(
            "sphere x R2 study: |v| = {vn:.12} (unit along the whole path), v = ({:.4},{:.4},{:.4}), z = ({:.4},{:.4}), f = {:.6}",
            study.x[0], study.x[1], study.x[2], study.x[3], study.x[4], rep.f
        ),
    );
}

#[test]
fn runner_budget_threads_into_stop_algebra() {
    let (problem, x0) = product_problem_budgeted(EvalLimit::Limited(
        NonZeroU64::new(50).expect("fixture evaluation limit is nonzero"),
    ));
    let mut study = Study::new(&problem, &x0, 1e-6, 0.2);
    let rep = study.run(&problem, &StopRule::GradNorm(1e-12), 4000);
    verdict(
        "ijil-runner-budget",
        rep.reason == StopReason::Budget && rep.evals >= 50 && rep.evals < 80,
        &format!(
            "budget 50 evals: stopped with {:?} at {} evals (the problem's own P4 budget, not the caller's rule)",
            rep.reason, rep.evals
        ),
    );
}

#[test]
fn runner_resume_is_bitwise() {
    let (problem, x0) = product_problem();
    let admission = problem.admit().expect("runner fixture must admit");
    let problem_semantic_id = admission.semantic_id().to_hex();
    let mut fixture_builder = IdentityBuilder::new("fs-ascent-problem-ir-study-fixture-v1")
        .str("problem-semantic-id", &problem_semantic_id)
        .u64(
            "admission-schema-version",
            u64::from(admission.schema_version()),
        )
        .str("engine", "Study/projected-gradient-v1")
        .f64_bits("fd-h", REPLAY_FD_H)
        .f64_bits("learning-rate", REPLAY_LR)
        .f64_bits("gradient-stop", REPLAY_GRAD_TOL)
        .u64(
            "maximum-steps",
            u64::try_from(REPLAY_STEPS).expect("replay steps fit u64"),
        )
        .u64("input-seed", FIXED_INPUT_SEED)
        .str("fs-ascent-version", fs_ascent::VERSION);
    for &value in &x0 {
        fixture_builder = fixture_builder.f64_bits("initial-point", value);
    }
    let fixture_identity = fixture_builder.finish();
    let rule = StopRule::GradNorm(REPLAY_GRAD_TOL);

    let mut straight = Study::new(&problem, &x0, REPLAY_FD_H, REPLAY_LR);
    let straight_report = straight.run(&problem, &rule, REPLAY_STEPS);
    let straight_identity = study_state_identity(&fixture_identity, &straight, &straight_report);

    let mut repeat = Study::new(&problem, &x0, REPLAY_FD_H, REPLAY_LR);
    let repeat_report = repeat.run(&problem, &rule, REPLAY_STEPS);
    let repeat_identity = study_state_identity(&fixture_identity, &repeat, &repeat_report);
    let repeat_exact = straight_identity.canonical_bytes() == repeat_identity.canonical_bytes();

    let sphere_norm = fs_math::det::sqrt(straight.x[..3].iter().map(|v| v * v).sum::<f64>());
    let expected_evals = expected_study_evals(&straight);
    let quality_pass = straight_report.reason == StopReason::GradNorm
        && straight_report.grad_norm <= REPLAY_GRAD_TOL
        && straight_report.f < -1.0 + 1e-8
        && straight_report.evals == straight.evals
        && straight.history.len() == straight.steps + 1
        && straight.evals == expected_evals
        && (sphere_norm - 1.0).abs() <= 1e-12;
    let mut pass = repeat_exact && quality_pass;
    let mut first_failure = if !quality_pass {
        Some("straight-quality".to_string())
    } else if !repeat_exact {
        Some("independent-repeat".to_string())
    } else {
        None
    };
    let mut cut_identities = Vec::with_capacity(REPLAY_CUTS.len());
    let mut checkpoint_receipts = Vec::with_capacity(REPLAY_CUTS.len());
    for cut in REPLAY_CUTS {
        let mut first = Study::new(&problem, &x0, REPLAY_FD_H, REPLAY_LR);
        let first_report = first.run(&problem, &rule, cut);
        let checkpoint_identity = study_state_identity(&fixture_identity, &first, &first_report);
        let genuine_split = first_report.reason == StopReason::IterationCap
            && first.steps == cut
            && first.history.len() == first.steps
            && first_report.evals == first.evals
            && first.evals == expected_study_evals(&first);
        let mut resumed = first.clone(); // clone = checkpoint
        let resumed_report = resumed.run(&problem, &rule, REPLAY_STEPS - cut);
        let resumed_identity = study_state_identity(&fixture_identity, &resumed, &resumed_report);
        let exact = resumed_identity.canonical_bytes() == straight_identity.canonical_bytes();
        let checkpoint_resume_identity =
            IdentityBuilder::new("fs-ascent-problem-ir-study-checkpoint-resume-v1")
                .child("fixture", &fixture_identity)
                .u64(
                    "cut-steps",
                    u64::try_from(cut).expect("replay cut fits u64"),
                )
                .child("checkpoint-state", &checkpoint_identity)
                .child("resumed-state", &resumed_identity)
                .finish();
        if (!genuine_split || !exact) && first_failure.is_none() {
            first_failure = Some(format!(
                "cut-{cut}: genuine_split={genuine_split}, exact={exact}"
            ));
        }
        pass &= genuine_split && exact;
        checkpoint_receipts.push(format!(
            "{{\"cut\":{cut},\"checkpoint_state\":\"{}\",\
             \"checkpoint_resume_identity\":\"{}\",\
             \"reason\":\"{}\",\"steps\":{},\"history_values\":{},\
             \"study_evals\":{},\"report_evals\":{},\"genuine_split\":{genuine_split}}}",
            checkpoint_identity.hex(),
            checkpoint_resume_identity.hex(),
            stop_reason_name(&first_report.reason),
            first.steps,
            first.history.len(),
            first.evals,
            first_report.evals,
        ));
        cut_identities.push(checkpoint_resume_identity);
    }

    let first_failure_json = first_failure
        .as_ref()
        .map_or_else(|| "null".to_string(), |failure| format!("\"{failure}\""));
    let receipt_json = format!(
        "{{\"problem_semantic_id\":\"{problem_semantic_id}\",\
         \"fixture_identity\":\"{}\",\"input_seed\":{FIXED_INPUT_SEED},\
         \"cuts\":{},\"reference_state\":\"{}\",\
         \"repeat_state\":\"{}\",\"cut_states\":{},\
         \"checkpoints\":[{}],\
         \"initial_point_bits\":{},\"final_point_bits\":{},\
         \"objective_history_bits\":{},\"steps\":{},\"study_evals\":{},\
         \"expected_evals\":{},\"report_evals\":{},\
         \"report_objective_bits\":\"0x{:016x}\",\
         \"report_gradient_norm_bits\":\"0x{:016x}\",\"stop_reason\":\"{}\",\
         \"first_failure\":{first_failure_json},\"pass\":{pass}}}",
        fixture_identity.hex(),
        usize_vector_json(&REPLAY_CUTS),
        straight_identity.hex(),
        repeat_identity.hex(),
        identity_vector_json(&cut_identities),
        checkpoint_receipts.join(","),
        bit_vector_json(&x0),
        bit_vector_json(&straight.x),
        bit_vector_json(&straight.history),
        straight.steps,
        straight.evals,
        expected_evals,
        straight_report.evals,
        straight_report.f.to_bits(),
        straight_report.grad_norm.to_bits(),
        stop_reason_name(&straight_report.reason),
    );
    replay_verdict(
        "ijil-runner-resume",
        pass,
        &format!(
            "problem={problem_semantic_id}; fixture={}; final={}; repeat={}; cuts={}; \
             reason={}; steps={}; evals={}; history={}; f={:.17e}; grad_norm={:.17e}; \
             expected_evals={expected_evals}; sphere_norm={sphere_norm:.17e}; \
             first_failure={first_failure:?}",
            fixture_identity.hex(),
            straight_identity.hex(),
            repeat_identity.hex(),
            identity_vector_json(&cut_identities),
            stop_reason_name(&straight_report.reason),
            straight.steps,
            straight.evals,
            straight.history.len(),
            straight_report.f,
            straight_report.grad_norm,
        ),
        receipt_json,
    );
}

#[test]
fn runner_constraints_route_to_al() {
    // The landed fixture expressed AS DATA: minimize (x−2)² + (y−1)²
    // s.t. x + y = 2 (EqZero), x − 1.2 ≤ 0 (LeZero), solved through
    // the packed adapters.
    let mut b = ProblemBuilder::new();
    let xy = b.var("xy", Manifold::Rn { dim: 2 }, D0).expect("var xy");
    let xyref = b.var_ref(xy).expect("xy node");
    let x0c = b.component(xyref, 0).expect("x");
    let x1c = b.component(xyref, 1).expect("y");
    let dx = affine(&mut b, x0c, 1.0, -2.0);
    let dy = affine(&mut b, x1c, 1.0, -1.0);
    let qx = b.mul(dx, dx).expect("qx");
    let qy = b.mul(dy, dy).expect("qy");
    let obj = b.add(qx, qy).expect("obj");
    b.objective(obj, Sense::Minimize, 1.0).expect("objective");
    let sum = b.add(x0c, x1c).expect("sum");
    let ce_node = affine(&mut b, sum, 1.0, -2.0);
    b.constraint(ce_node, ConstraintKind::EqZero, "sum-to-2")
        .expect("eq constraint");
    let ci_node = affine(&mut b, x0c, 1.0, -1.2);
    b.constraint(ci_node, ConstraintKind::LeZero, "x-cap")
        .expect("ineq constraint");
    let problem = b.finish();

    let packing = Packing::new(&problem);
    let (ce, ce_jt, ci, ci_jt) = Study::constraint_adapters(&problem, &packing, 1e-6);
    let mut fg = |x: &[f64]| -> (f64, Vec<f64>) {
        (
            (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
            vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)],
        )
    };
    let mut p = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };
    let rep = augmented_lagrangian(&mut p, &[0.0, 0.0], 1e-6, 40);
    verdict(
        "ijil-runner-constraints",
        rep.converged && (rep.x[0] - 1.2).abs() < 1e-4 && (rep.x[1] - 0.8).abs() < 1e-4,
        &format!(
            "IR-declared constraints through packed adapters: x = ({:.5}, {:.5}), kkt = ({:.1e},{:.1e},{:.1e})",
            rep.x[0], rep.x[1], rep.kkt.stationarity, rep.kkt.feasibility, rep.kkt.complementarity
        ),
    );
}
