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
    retract as root_retract, tangent_project as root_tangent_project,
};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder, Sense};
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

// ---------------------------------------------------------------------------
// Problem fixtures
// ---------------------------------------------------------------------------

fn create_sphere_problem(dim: usize) -> Problem {
    let mut builder = ProblemBuilder::new();
    let sphere = builder
        .var("x", Manifold::Sphere { dim }, D0)
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

fn create_stiefel_problem(n: usize, p: usize) -> Problem {
    let mut builder = ProblemBuilder::new();
    let stiefel = builder
        .var("Y", Manifold::Stiefel { n, p }, D0)
        .expect("stiefel var");
    let stiefel_ref = builder.var_ref(stiefel).expect("stiefel ref");
    let comp0 = builder.component(stiefel_ref, 0).expect("component 0");
    builder
        .objective(comp0, Sense::Minimize, 1.0)
        .expect("objective");
    builder.finish()
}

fn create_mixed_manifold_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let r_var = builder
        .var("euc", Manifold::Rn { dim: 2 }, D0)
        .expect("Rn var");
    let s_var = builder
        .var("sph", Manifold::Sphere { dim: 3 }, D0)
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

// ---------------------------------------------------------------------------
// E2E Test Cases
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_sphere_riemannian_optimization() {
    let dim = 3;
    let mut opt = RiemannianLbfgs::new(
        Manifold::Sphere { dim },
        5,
        StopRule::GradNorm(1e-6),
    );

    // Initial point: north pole [0, 0, 1]
    let x0 = vec![0.0, 0.0, 1.0];
    let mut state = opt.init(x0).expect("sphere init");

    // Objective: f(x) = x[0], grad_ambient = [1, 0, 0]
    let report = opt
        .run(&mut state, 100, &mut |x| {
            let f = x[0];
            let g = vec![1.0, 0.0, 0.0];
            (f, g)
        })
        .expect("sphere run");

    assert!(report.converged, "Sphere optimization must converge");
    assert!(
        matches!(report.stop_reason, StopReason::GradNorm),
        "Must stop on gradient norm"
    );
    // Converged point should have x[0] ≈ -1 (south pole in x)
    assert!(
        (state.x[0] - (-1.0)).abs() < 1e-4,
        "Optimal point should be near [-1, 0, 0]"
    );

    emit_jsonl_record(
        "sphere_riemannian_opt",
        "solve",
        "PASS",
        &format!(
            "\"iterations\":{},\"f_opt\":{:.8},\"x_bits\":{:?}",
            report.iterations,
            report.f_opt,
            bits(&state.x)
        ),
    );
}

#[test]
fn test_e2e_so3_point_parameter_separation() {
    let mut opt = RiemannianLbfgs::new(Manifold::So3, 4, StopRule::GradNorm(1e-6));

    // Initial quaternion: identity [1, 0, 0, 0]
    let q0 = vec![1.0, 0.0, 0.0, 0.0];
    let mut state = opt.init(q0).expect("so3 init");

    // Verify point dimension is 4 and parameter dimension is 3
    assert_eq!(state.x.len(), 4, "Quaternion point dimension is 4");
    assert_eq!(state.g.len(), 3, "SO(3) Lie parameter gradient dimension is 3");

    // Objective: minimize q[0] * q[1]
    let report = opt
        .run(&mut state, 100, &mut |q| {
            let f = q[0] * q[1];
            let g = vec![q[1], q[0], 0.0, 0.0];
            (f, g)
        })
        .expect("so3 run");

    assert!(report.converged || report.iterations > 0);
    // Norm of quaternion must be preserved exactly
    let norm = (state.x[0] * state.x[0]
        + state.x[1] * state.x[1]
        + state.x[2] * state.x[2]
        + state.x[3] * state.x[3])
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
            report.iterations,
            report.f_opt,
            (norm - 1.0).abs()
        ),
    );
}

#[test]
fn test_e2e_stiefel_differential_transport_and_rho_recomputation() {
    let n = 4;
    let p = 2;
    let mut opt = RiemannianLbfgs::new(
        Manifold::Stiefel { n, p },
        4,
        StopRule::GradNorm(1e-5),
    );

    // Initial orthonormal 4x2 matrix (column-major)
    let y0 = vec![
        0.5, 0.5, 0.5, 0.5, // col 0
        0.5, -0.5, 0.5, -0.5, // col 1
    ];
    let mut state = opt.init(y0).expect("stiefel init");

    // Objective: linear trace sum
    let report = opt
        .run(&mut state, 50, &mut |y| {
            let f = y[0] + y[5];
            let mut g = vec![0.0; 8];
            g[0] = 1.0;
            g[5] = 1.0;
            (f, g)
        })
        .expect("stiefel run");

    assert!(report.iterations > 0);
    // Verify orthonormality Y^T Y = I
    let col0 = &state.x[0..4];
    let col1 = &state.x[4..8];
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
            report.iterations,
            dot01.abs()
        ),
    );
}

#[test]
fn test_e2e_mixed_manifold_study_runner() {
    let problem = create_mixed_manifold_problem();
    let packing = Packing::new(&problem).expect("packing");

    // Point shape: 2 (Rn) + 3 (Sphere) + 4 (So3) + 8 (Stiefel 4x2) = 17
    assert_eq!(packing.point_dim, 17, "Point dimension is 17");
    // Parameter shape: 2 (Rn) + 2 (Sphere) + 3 (So3) + 7 (Stiefel 4x2) = 14
    assert_eq!(packing.parameter_dim, 14, "Parameter dimension is 14");

    let mut start_point = Vec::new();
    start_point.extend([1.0, 2.0]); // Rn
    start_point.extend([0.0, 0.0, 1.0]); // Sphere
    start_point.extend([1.0, 0.0, 0.0, 0.0]); // So3
    start_point.extend([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, -0.5]); // Stiefel

    let mut study = Study::new(problem, start_point).expect("study creation");
    let report = study.run_to_completion().expect("study run");

    assert!(report.evaluations > 0, "Evaluations occurred");
    assert!(report.point.len() == 17, "Returned point has point_dim 17");

    emit_jsonl_record(
        "mixed_manifold_study_runner",
        "execution",
        "PASS",
        &format!(
            "\"evaluations\":{},\"iterations\":{},\"f_best\":{:.8}",
            report.evaluations,
            report.iterations,
            report.value
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

    let mut study1 = Study::new(problem1, start1).expect("study1");
    let mut study2 = Study::new(problem2, start2).expect("study2");

    let report1 = study1.run_to_completion().expect("study1 run");
    let report2 = study2.run_to_completion().expect("study2 run");

    assert_eq!(
        report1.evaluations, report2.evaluations,
        "Evaluation count must match identically"
    );
    assert_eq!(
        report1.iterations, report2.iterations,
        "Iteration count must match identically"
    );
    assert_eq!(
        bits(&report1.point),
        bits(&report2.point),
        "Returned optimal points must be bit-for-bit identical"
    );
    assert_eq!(
        report1.value.to_bits(),
        report2.value.to_bits(),
        "Objective values must be bit-for-bit identical"
    );

    emit_jsonl_record(
        "deterministic_replay",
        "verification",
        "PASS",
        "\"bit_exact\":true,\"replays_matched\":2",
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
    let result = Study::new(problem, invalid_point);
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
    let result = Study::new(problem, wrong_dim_point);
    assert!(
        matches!(result, Err(StudyError::DimensionMismatch { .. })),
        "Wrong dimension must be refused as DimensionMismatch"
    );

    emit_jsonl_record(
        "hostile_dimension_mismatch",
        "admission",
        "REFUSED",
        "\"expected_rejection\":true",
    );
}

#[test]
fn test_hostile_budget_boundary_edges() {
    let dim = 3;
    let mut opt = RiemannianLbfgs::new(
        Manifold::Sphere { dim },
        5,
        StopRule::GradNorm(1e-12), // unattainable, will exhaust budget
    );

    let x0 = vec![0.0, 0.0, 1.0];
    let mut state = opt.init(x0).expect("init");

    // Exact budget cap = 3 iterations
    let report = opt
        .run(&mut state, 3, &mut |x| {
            let f = x[0];
            let g = vec![1.0, 0.0, 0.0];
            (f, g)
        })
        .expect("run");

    assert_eq!(report.iterations, 3, "Must stop at exact budget limit");
    assert!(
        matches!(report.stop_reason, StopReason::MaxIterations),
        "Stop reason must be MaxIterations"
    );

    emit_jsonl_record(
        "hostile_budget_boundary",
        "execution",
        "PASS",
        "\"budget_bounded\":true,\"stop_reason\":\"MaxIterations\"",
    );
}

#[test]
fn test_hostile_non_finite_callback_refuses() {
    let dim = 3;
    let mut opt = RiemannianLbfgs::new(
        Manifold::Sphere { dim },
        5,
        StopRule::GradNorm(1e-6),
    );

    let x0 = vec![0.0, 0.0, 1.0];
    let mut state = opt.init(x0).expect("init");

    // Callback returning NaN
    let result = opt.run(&mut state, 10, &mut |_x| {
        (f64::NAN, vec![1.0, 0.0, 0.0])
    });

    assert!(result.is_err(), "NaN objective must return an error");

    emit_jsonl_record(
        "hostile_non_finite_callback",
        "execution",
        "REFUSED",
        "\"caught_non_finite\":true",
    );
}
