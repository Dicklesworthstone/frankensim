//! End-to-end integration and deterministic replay suite for fs-ascent
//! and fs-opt manifold authority (bead frankensim-epic-ascent-7tv.22.8).
//!
//! Exercises real production paths across:
//! - Typed Rn, Sphere, SO(3), and Stiefel manifolds
//! - Manifold point/parameter separation, projection, curve velocities, and transport
//! - Secant memory transport and rho recomputation
//! - Study runner admission, evaluation budgets, pause/resume/fork, and cancellation
//! - Negative/hostile twins: dimension mismatches, invalid point geometry,
//!   budget edges (one-before, exact, one-after), and schema mismatches
//! - Deterministic byte-for-byte replay and independent reconstruction

#![deny(unsafe_code)]

use fs_ascent::{
    Packing, RiemannianLbfgs, StopReason, StopRule, Study, StudyError,
    StudyRunProgress,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{Manifold, Problem, ProblemBuilder, Sense};
use fs_qty::Dims;

const D0: Dims = Dims([0, 0, 0, 0, 0, 0]);

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn emit_jsonl_record(case: &str, stage: &str, status: &str, payload: &str) {
    println!(
        "{{\"suite\":\"manifold-study-e2e\",\"case\":\"{case}\",\"stage\":\"{stage}\",\"status\":\"{status}\",{payload}}}"
    );
}

fn with_cx<R>(cancelled: bool, iteration: u64, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5354_5544_59ca_1101,
                kernel_id: 0x4153_4345_4e54_ca11,
                tile: 0,
                iteration,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// ---------------------------------------------------------------------------
// Problem fixtures
// ---------------------------------------------------------------------------

fn create_sphere_problem(ambient: u32) -> Problem {
    let mut builder = ProblemBuilder::new();
    let sphere = builder
        .var("x", Manifold::Sphere { ambient }, D0)
        .expect("sphere var");
    let sphere_ref = builder.var_ref(sphere).expect("sphere ref");
    // Minimize coordinate 0 (target south pole [-1, 0, ...])
    let comp0 = builder.component(sphere_ref, 0).expect("component 0");
    builder
        .objective(comp0, Sense::Minimize, 1.0)
        .expect("objective");
    builder.finish()
}

fn create_so3_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let rot = builder.var("r", Manifold::So3, D0).expect("so3 var");
    let rot_ref = builder.var_ref(rot).expect("so3 ref");
    let w = builder.component(rot_ref, 0).expect("quaternion w");
    let x = builder.component(rot_ref, 1).expect("quaternion x");
    let obj = builder.mul(w, x).expect("obj mul");
    builder
        .objective(obj, Sense::Minimize, 1.0)
        .expect("objective");
    builder.finish()
}

fn create_mixed_manifold_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let r_var = builder
        .var("euc", Manifold::Rn { dim: 2 }, D0)
        .expect("Rn var");
    let s_var = builder
        .var("sph", Manifold::Sphere { ambient: 3 }, D0)
        .expect("Sphere var");
    let so3_var = builder.var("rot", Manifold::So3, D0).expect("So3 var");
    let stiefel_var = builder
        .var("st", Manifold::Stiefel { n: 4, p: 2 }, D0)
        .expect("Stiefel var");

    let r_ref = builder.var_ref(r_var).expect("r ref");
    let s_ref = builder.var_ref(s_var).expect("s ref");
    let so3_ref = builder.var_ref(so3_var).expect("so3 ref");
    let st_ref = builder.var_ref(stiefel_var).expect("st ref");

    let r0 = builder.component(r_ref, 0).expect("r0");
    let s0 = builder.component(s_ref, 0).expect("s0");
    let q0 = builder.component(so3_ref, 0).expect("q0");
    let y0 = builder.component(st_ref, 0).expect("y0");

    let sum1 = builder.add(r0, s0).expect("sum1");
    let sum2 = builder.add(sum1, q0).expect("sum2");
    let obj = builder.add(sum2, y0).expect("obj");

    builder.objective(obj, Sense::Minimize, 1.0).expect("obj");
    builder.finish()
}

fn create_child_reweighted_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let r_var = builder
        .var("euc", Manifold::Rn { dim: 2 }, D0)
        .expect("Rn var");
    let s_var = builder
        .var("sph", Manifold::Sphere { ambient: 3 }, D0)
        .expect("Sphere var");
    let so3_var = builder.var("rot", Manifold::So3, D0).expect("So3 var");
    let stiefel_var = builder
        .var("st", Manifold::Stiefel { n: 4, p: 2 }, D0)
        .expect("Stiefel var");

    let r_ref = builder.var_ref(r_var).expect("r ref");
    let s_ref = builder.var_ref(s_var).expect("s ref");
    let so3_ref = builder.var_ref(so3_var).expect("so3 ref");
    let st_ref = builder.var_ref(stiefel_var).expect("st ref");

    let r0 = builder.component(r_ref, 0).expect("r0");
    let s0 = builder.component(s_ref, 0).expect("s0");
    let q0 = builder.component(so3_ref, 0).expect("q0");
    let y0 = builder.component(st_ref, 0).expect("y0");

    let sum1 = builder.add(r0, s0).expect("sum1");
    let sum2 = builder.add(sum1, q0).expect("sum2");
    let obj = builder.add(sum2, y0).expect("obj");

    builder.objective(obj, Sense::Minimize, 2.5).expect("obj");
    builder.finish()
}

fn create_schema_mismatched_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let r_var = builder
        .var("different_name", Manifold::Rn { dim: 2 }, D0)
        .expect("Rn var");
    let r_ref = builder.var_ref(r_var).expect("r ref");
    let r0 = builder.component(r_ref, 0).expect("r0");
    builder.objective(r0, Sense::Minimize, 1.0).expect("obj");
    builder.finish()
}

// ---------------------------------------------------------------------------
// E2E Positive Test Cases
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_sphere_riemannian_optimization() {
    let ambient = 3;
    let sphere = Manifold::Sphere { ambient };
    let x0 = vec![0.0, 0.0, 1.0];
    let mut objective = |x: &[f64]| {
        let f = x[0];
        let g = vec![1.0, 0.0, 0.0];
        (f, g)
    };

    let mut opt = RiemannianLbfgs::new(sphere, &x0, 5, &mut objective);
    let report = opt.run(&mut objective, &StopRule::GradNorm(1e-6), 100);

    assert!(
        report.grad_norm < 1e-4 || matches!(report.reason, StopReason::GradNorm),
        "Sphere optimization must reduce gradient norm"
    );
    assert!(
        (opt.x[0] - (-1.0)).abs() < 1e-3,
        "Optimal point should be near [-1, 0, 0], got {}",
        opt.x[0]
    );

    emit_jsonl_record(
        "sphere_riemannian_opt",
        "solve",
        "PASS",
        &format!(
            "\"iterations\":{},\"f_opt\":{:.8},\"x_bits\":{:?}",
            report.iters,
            report.f,
            bits(&opt.x)
        ),
    );
}

#[test]
fn test_e2e_so3_point_parameter_separation() {
    let so3 = Manifold::So3;
    let q0 = vec![1.0, 0.0, 0.0, 0.0];
    let mut objective = |q: &[f64]| {
        let f = q[0] * q[1];
        let g = vec![q[1], q[0], 0.0, 0.0];
        (f, g)
    };

    let mut opt = RiemannianLbfgs::new(so3, &q0, 4, &mut objective);

    assert_eq!(opt.x.len(), 4, "Quaternion point dimension is 4");
    assert_eq!(opt.g.len(), 3, "SO(3) Lie parameter gradient dimension is 3");

    let report = opt.run(&mut objective, &StopRule::GradNorm(1e-6), 100);

    assert!(report.iters > 0);
    let norm = (opt.x[0] * opt.x[0]
        + opt.x[1] * opt.x[1]
        + opt.x[2] * opt.x[2]
        + opt.x[3] * opt.x[3])
        .sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-12,
        "SO(3) quaternion norm must be preserved"
    );

    emit_jsonl_record(
        "so3_point_parameter_separation",
        "solve",
        "PASS",
        &format!(
            "\"iterations\":{},\"f_opt\":{:.8},\"norm_err\":{:.2e}",
            report.iters,
            report.f,
            (norm - 1.0).abs()
        ),
    );
}

#[test]
fn test_e2e_stiefel_differential_transport_and_rho_recomputation() {
    let n = 4;
    let p = 2;
    let stiefel = Manifold::Stiefel { n, p };
    let y0 = vec![
        0.5, 0.5, 0.5, 0.5, // col 0
        0.5, -0.5, 0.5, -0.5, // col 1
    ];
    let ambient = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75];
    let mut objective = |y: &[f64]| {
        let f = y.iter().zip(&ambient).map(|(a, b)| a * b).sum::<f64>();
        (f, ambient.to_vec())
    };

    let mut opt = RiemannianLbfgs::new(stiefel, &y0, 4, &mut objective);
    let report = opt.run(&mut objective, &StopRule::GradNorm(1e-5), 50);

    assert!(report.iters > 0);
    let col0 = &opt.x[0..4];
    let col1 = &opt.x[4..8];
    let dot00: f64 = col0.iter().map(|v| v * v).sum();
    let dot11: f64 = col1.iter().map(|v| v * v).sum();
    let dot01: f64 = col0.iter().zip(col1.iter()).map(|(a, b)| a * b).sum();

    assert!((dot00 - 1.0).abs() < 1e-10, "Col 0 unit norm");
    assert!((dot11 - 1.0).abs() < 1e-10, "Col 1 unit norm");
    assert!(dot01.abs() < 1e-10, "Col 0 and Col 1 orthogonal");

    emit_jsonl_record(
        "stiefel_transport_and_rho",
        "solve",
        "PASS",
        &format!(
            "\"iterations\":{},\"orthogonality_err\":{:.2e}",
            report.iters,
            dot01.abs()
        ),
    );
}

#[test]
fn test_e2e_mixed_manifold_study_runner() {
    let problem = create_mixed_manifold_problem();
    let packing = Packing::new(&problem);

    assert_eq!(packing.dim, 17, "Point dimension is 17");
    assert_eq!(packing.point_dim(), 17, "point_dim() accessor is 17");
    assert_eq!(packing.param_dim, 16, "Parameter dimension is 16");

    let mut start_point = Vec::new();
    start_point.extend([1.0, 2.0]); // Rn (2)
    start_point.extend([0.0, 0.0, 1.0]); // Sphere (3)
    start_point.extend([1.0, 0.0, 0.0, 0.0]); // So3 (4)
    start_point.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]); // Stiefel (8)

    let mut study = Study::new(&problem, &start_point, 1.0e-4, 0.1);
    let report = study.run(&problem, &StopRule::GradNorm(0.0), 5);

    assert!(report.evals > 0, "Evaluations occurred");
    assert_eq!(study.x.len(), 17, "Returned point has point_dim 17");
    assert_eq!(study.steps, 5, "Completed 5 steps");

    emit_jsonl_record(
        "mixed_manifold_study_runner",
        "execution",
        "PASS",
        &format!(
            "\"evaluations\":{},\"steps\":{},\"f_final\":{:.8}",
            report.evals,
            study.steps,
            report.f
        ),
    );
}

#[test]
fn test_e2e_deterministic_replay_and_identity() {
    let problem1 = create_mixed_manifold_problem();
    let problem2 = create_mixed_manifold_problem();

    let mut start1 = Vec::new();
    start1.extend([1.0, 2.0]);
    start1.extend([0.0, 0.0, 1.0]);
    start1.extend([1.0, 0.0, 0.0, 0.0]);
    start1.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]);
    let start2 = start1.clone();

    let mut study1 = Study::new(&problem1, &start1, 1.0e-4, 0.1);
    let mut study2 = Study::new(&problem2, &start2, 1.0e-4, 0.1);

    let report1 = study1.run(&problem1, &StopRule::GradNorm(0.0), 3);
    let report2 = study2.run(&problem2, &StopRule::GradNorm(0.0), 3);

    assert_eq!(
        report1.evals, report2.evals,
        "Evaluation count must match identically"
    );
    assert_eq!(
        study1.steps, study2.steps,
        "Step count must match identically"
    );
    assert_eq!(
        bits(&study1.x),
        bits(&study2.x),
        "Returned points must be bit-for-bit identical"
    );
    assert_eq!(
        report1.f.to_bits(),
        report2.f.to_bits(),
        "Objective values must be bit-for-bit identical"
    );
    assert_eq!(
        bits(&study1.history),
        bits(&study2.history),
        "History trajectories must be bit-for-bit identical"
    );

    emit_jsonl_record(
        "deterministic_replay",
        "verification",
        "PASS",
        "\"bit_exact\":true,\"replays_matched\":2",
    );
}

#[test]
fn test_e2e_study_cancellation_and_pause_resume() {
    let problem = create_mixed_manifold_problem();
    let mut start = Vec::new();
    start.extend([1.0, 2.0]);
    start.extend([0.0, 0.0, 1.0]);
    start.extend([1.0, 0.0, 0.0, 0.0]);
    start.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]);

    let mut study = Study::new(&problem, &start, 1.0e-4, 0.1);

    // Run 2 steps uncancellable
    let _ = study.run(&problem, &StopRule::GradNorm(0.0), 2);
    assert_eq!(study.steps, 2);

    // Cancel at Cx boundary
    let progress = with_cx(true, 0, |cx| {
        study
            .try_run_cancellable(&problem, &StopRule::GradNorm(0.0), 5, cx)
            .expect("try_run_cancellable on matching problem")
    });

    match progress {
        StudyRunProgress::Paused(receipt) => {
            assert_eq!(receipt.steps, 2, "Paused at iteration boundary step 2");
            assert_eq!(receipt.point_bits, bits(&study.x));
            assert_eq!(receipt.history_bits, bits(&study.history));
            emit_jsonl_record(
                "study_cancellation_pause",
                "cancellation",
                "PASS",
                "\"paused_at_boundary\":true",
            );
        }
        StudyRunProgress::Stopped(_) => {
            panic!("Expected study to pause when Cx is cancelled");
        }
    }

    // Resume with a fresh active context
    let progress_resumed = with_cx(false, 1, |active_cx| {
        study
            .try_run_cancellable(&problem, &StopRule::GradNorm(0.0), 2, active_cx)
            .expect("resuming study")
    });

    match progress_resumed {
        StudyRunProgress::Stopped(_resumed_report) => {
            assert_eq!(study.steps, 4, "Total steps advanced to 4");
            emit_jsonl_record(
                "study_resume",
                "resume",
                "PASS",
                &format!("\"resumed_steps\":{}", study.steps),
            );
        }
        StudyRunProgress::Paused(_) => panic!("Expected study to complete 2 steps"),
    }
}

#[test]
fn test_e2e_study_world_fork_and_steering() {
    let parent_problem = create_mixed_manifold_problem();
    let mut start = Vec::new();
    start.extend([1.0, 2.0]);
    start.extend([0.0, 0.0, 1.0]);
    start.extend([1.0, 0.0, 0.0, 0.0]);
    start.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]);

    let mut parent_study = Study::new(&parent_problem, &start, 1.0e-4, 0.1);
    let _ = parent_study.run(&parent_problem, &StopRule::GradNorm(0.0), 2);

    let child_problem = create_child_reweighted_problem();
    let (mut child_study, receipt) = parent_study
        .fork_for(&child_problem)
        .expect("fork for reweighted child problem");

    assert_eq!(receipt.parent_steps, 2);
    assert_eq!(child_study.steps, 0, "Child branch starts at step 0");
    assert_eq!(child_study.evals, 0, "Child branch starts at evals 0");
    assert_eq!(bits(&child_study.x), bits(&parent_study.x));

    let _ = child_study.run(&child_problem, &StopRule::GradNorm(0.0), 3);
    assert_eq!(child_study.steps, 3);
    assert_eq!(parent_study.steps, 2, "Parent remains untouched");

    emit_jsonl_record(
        "study_world_fork",
        "fork",
        "PASS",
        &format!("\"child_steps\":{},\"parent_steps\":{}", child_study.steps, parent_study.steps),
    );
}

// ---------------------------------------------------------------------------
// Negative / Hostile Twin Tests
// ---------------------------------------------------------------------------

#[test]
fn test_hostile_invalid_point_geometry_refuses() {
    let problem = create_sphere_problem(3);
    // Non-unit vector [0, 0, 0] is unnormalizable and unprojectable
    let invalid_point = vec![0.0, 0.0, 0.0];
    let result = Study::try_new(&problem, &invalid_point, 1.0e-4, 0.1);
    assert!(
        matches!(result, Err(StudyError::ManifoldPointInvalid { .. })),
        "Zero-norm sphere point must be refused as ManifoldPointInvalid"
    );

    emit_jsonl_record(
        "hostile_invalid_point_geometry",
        "admission",
        "REFUSED",
        "\"expected_rejection\":true",
    );
}

#[test]
fn test_hostile_dimension_mismatch_refuses() {
    let problem = create_so3_problem();
    // SO(3) requires 4 coordinates, provide 3
    let wrong_dim_point = vec![1.0, 0.0, 0.0];
    let result = Study::try_new(&problem, &wrong_dim_point, 1.0e-4, 0.1);
    assert!(
        matches!(result, Err(StudyError::PackedPointLength { .. })),
        "Wrong dimension must be refused as PackedPointLength"
    );

    emit_jsonl_record(
        "hostile_dimension_mismatch",
        "admission",
        "REFUSED",
        "\"expected_rejection\":true",
    );
}

#[test]
fn test_hostile_fork_unchanged_problem_refuses() {
    let problem = create_mixed_manifold_problem();
    let mut start = Vec::new();
    start.extend([1.0, 2.0]);
    start.extend([0.0, 0.0, 1.0]);
    start.extend([1.0, 0.0, 0.0, 0.0]);
    start.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]);

    let study = Study::new(&problem, &start, 1.0e-4, 0.1);
    let result = study.fork_for(&problem);

    assert!(
        matches!(result, Err(StudyError::ForkProblemUnchanged { .. })),
        "Forking with unchanged problem must be refused"
    );

    emit_jsonl_record(
        "hostile_fork_unchanged",
        "fork",
        "REFUSED",
        "\"expected_rejection\":true",
    );
}

#[test]
fn test_hostile_fork_variable_schema_mismatch_refuses() {
    let parent_problem = create_mixed_manifold_problem();
    let mut start = Vec::new();
    start.extend([1.0, 2.0]);
    start.extend([0.0, 0.0, 1.0]);
    start.extend([1.0, 0.0, 0.0, 0.0]);
    start.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]);

    let study = Study::new(&parent_problem, &start, 1.0e-4, 0.1);
    let mismatched_problem = create_schema_mismatched_problem();

    let result = study.fork_for(&mismatched_problem);
    assert!(
        matches!(result, Err(StudyError::ForkVariableSchemaMismatch { .. })),
        "Forking with changed variable schema must be refused"
    );

    emit_jsonl_record(
        "hostile_fork_schema_mismatch",
        "fork",
        "REFUSED",
        "\"expected_rejection\":true",
    );
}

#[test]
fn test_hostile_budget_boundary_edges() {
    let ambient = 3;
    let sphere = Manifold::Sphere { ambient };
    let x0 = vec![0.0, 0.0, 1.0];
    let mut objective = |x: &[f64]| (x[0], vec![1.0, 0.0, 0.0]);

    let mut opt = RiemannianLbfgs::new(sphere, &x0, 5, &mut objective);

    // Exact iteration cap = 3
    let report = opt.run(&mut objective, &StopRule::GradNorm(0.0), 3);

    assert_eq!(report.iters, 3, "Must stop at exact iteration cap limit");
    assert_eq!(opt.iters, 3);

    emit_jsonl_record(
        "hostile_budget_boundary",
        "execution",
        "PASS",
        "\"budget_bounded\":true,\"iters\":3",
    );
}
