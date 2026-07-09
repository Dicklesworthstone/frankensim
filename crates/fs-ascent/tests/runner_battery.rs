//! Bead ijil (d): the Problem-IR study runner — manifold-product
//! variable packing, budget threading through the stop algebra,
//! bitwise-resumable studies, and constrained routing through the
//! packed adapters.

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{Packing, StopReason, StopRule, Study, augmented_lagrangian};
use fs_opt::{ConstraintKind, Manifold, NodeId, ProblemBuilder, Sense};
use fs_qty::Dims;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

const D0: Dims = Dims([0, 0, 0, 0, 0]);

/// (a·x + b) as a scalar node from a component.
fn affine(b: &mut ProblemBuilder, x: NodeId, a: f64, off: f64) -> NodeId {
    let ca = b.konst(a, D0);
    let m = b.mul(ca, x).expect("scalar mul");
    let co = b.konst(off, D0);
    b.add(m, co).expect("scalar add")
}

/// A product-manifold problem: v on the unit sphere S², z in R²;
/// minimize 0.6·v₀ + 0.8·v₁ + (z₀−1)² + (z₁+2)² — linear on the
/// sphere (optimum v = −(0.6, 0.8, 0)) plus a shifted bowl.
fn product_problem() -> (fs_opt::Problem, Vec<f64>) {
    let mut b = ProblemBuilder::new();
    let v = b.var("v", Manifold::Sphere { ambient: 3 }, D0);
    let z = b.var("z", Manifold::Rn { dim: 2 }, D0);
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
    let (mut problem, x0) = product_problem();
    problem.budget = fs_opt::EvalBudget { max_evals: 50 };
    let mut study = Study::new(&problem, &x0, 1e-6, 0.2);
    let rep = study.run(&problem, &StopRule::GradNorm(1e-12), 4000);
    verdict(
        "ijil-runner-budget",
        rep.reason == StopReason::Composite && rep.evals >= 50 && rep.evals < 80,
        &format!(
            "budget 50 evals: stopped with {:?} at {} evals (the problem's own P4 budget, not the caller's rule)",
            rep.reason, rep.evals
        ),
    );
}

#[test]
fn runner_resume_is_bitwise() {
    let (problem, x0) = product_problem();
    let mut straight = Study::new(&problem, &x0, 1e-6, 0.2);
    straight.run(&problem, &StopRule::GradNorm(1e-9), 300);
    for cut in [1usize, 17, 111] {
        let mut a = Study::new(&problem, &x0, 1e-6, 0.2);
        a.run(&problem, &StopRule::GradNorm(1e-9), cut);
        let mut b = a.clone(); // clone = checkpoint
        b.run(&problem, &StopRule::GradNorm(1e-9), 300 - cut);
        let bitwise = b
            .x
            .iter()
            .zip(&straight.x)
            .all(|(u, v)| u.to_bits() == v.to_bits());
        assert!(bitwise, "resume not bitwise at cut {cut}");
        assert_eq!(
            b.evals, straight.evals,
            "eval accounting differs at cut {cut}"
        );
    }
    verdict(
        "ijil-runner-resume",
        true,
        "cuts 1/17/111: checkpointed studies bitwise == straight run (evals included)",
    );
}

#[test]
fn runner_constraints_route_to_al() {
    // The landed fixture expressed AS DATA: minimize (x−2)² + (y−1)²
    // s.t. x + y = 2 (EqZero), x − 1.2 ≤ 0 (LeZero), solved through
    // the packed adapters.
    let mut b = ProblemBuilder::new();
    let xy = b.var("xy", Manifold::Rn { dim: 2 }, D0);
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
            rep.x[0],
            rep.x[1],
            rep.kkt.stationarity,
            rep.kkt.feasibility,
            rep.kkt.complementarity
        ),
    );
}
