//! End-to-end constrained IR optimization without finite-difference gradients.
//!
//! cargo run -p fs-ascent --example reverse_ir_optimize
//!
//! This closed, globally smooth quadratic/affine fixture supplies fallible IR
//! evaluations to the existing infallible SQP callback API via expect. It is
//! not a general panic-free adapter for arbitrary domain-limited expressions.
//! The objective evaluation budget is explicitly unlimited; SQP has an outer
//! iteration cap. ReverseProblem itself does not own a solve-loop budget.
use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{SqpReport, sqp};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, EvalLimit, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

fn solve() -> SqpReport {
    let mut b = ProblemBuilder::new();
    let v = b.var("design", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    let y = b.component(r, 1).unwrap();
    let one = b.konst(1.0, Dims::NONE).unwrap();
    let two = b.konst(2.0, Dims::NONE).unwrap();
    let upper = b.konst(1.2, Dims::NONE).unwrap();
    let dx = b.sub(x, two).unwrap(); let dy = b.sub(y, one).unwrap();
    let dx2 = b.powi(dx, 2).unwrap(); let dy2 = b.powi(dy, 2).unwrap();
    let f = b.add(dx2, dy2).unwrap();
    let sum = b.add(x, y).unwrap();
    let equality = b.sub(sum, two).unwrap();
    let inequality = b.sub(x, upper).unwrap();
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.constraint(equality, ConstraintKind::EqZero, "x+y=2").unwrap();
    b.constraint(inequality, ConstraintKind::LeZero, "x<=1.2").unwrap();
    b.set_eval_limit(EvalLimit::Unlimited);
    let problem = b.finish();
    let oracle = ReverseProblem::new(&problem, ReverseLimits { max_nodes: 1 << 18, max_scalar_slots: 1 << 22 }).unwrap();

    let mut fg = |point: &[f64]| {
        let bindings = vec![point.to_vec()];
        let tape = oracle.evaluate(&bindings).expect("finite quadratic fixture");
        (tape.objective_value(), tape.objective_gradient().unwrap().remove(0))
    };
    let ce = |point: &[f64]| {
        let bindings = vec![point.to_vec()];
        vec![oracle.evaluate(&bindings).unwrap().constraint_values()[0]]
    };
    let ci = |point: &[f64]| {
        let bindings = vec![point.to_vec()];
        vec![oracle.evaluate(&bindings).unwrap().constraint_values()[1]]
    };
    let ce_jt = |point: &[f64], weights: &[f64]| {
        let bindings = vec![point.to_vec()];
        oracle.evaluate(&bindings).unwrap().constraint_pullback(&[weights[0], 0.0]).unwrap().remove(0)
    };
    let ci_jt = |point: &[f64], weights: &[f64]| {
        let bindings = vec![point.to_vec()];
        oracle.evaluate(&bindings).unwrap().constraint_pullback(&[0.0, weights[0]]).unwrap().remove(0)
    };
    let mut callbacks = ConstrainedProblem { fg: &mut fg, ce: &ce, ce_jt: &ce_jt, ci: &ci, ci_jt: &ci_jt };
    sqp(&mut callbacks, &[0.0, 0.0], 1e-7, 60)
}

fn main() {
    let report = solve();
    assert!(report.converged && report.kkt.within_tolerance(1e-7));
    assert!((report.x[0] - 1.2).abs() < 1e-6 && (report.x[1] - 0.8).abs() < 1e-6);
    println!("x={:?}, objective={}, converged={}, KKT={:?}", report.x, report.f, report.converged, report.kkt);
}

#[cfg(test)]
mod tests {
    #[test]
    fn compiled_ir_supplies_sqp_objective_and_constraint_derivatives() {
        let result = super::solve();
        assert!(result.converged);
        assert!(result.kkt.within_tolerance(1e-7));
        assert!((result.x[0] - 1.2).abs() < 1e-6);
        assert!((result.x[1] - 0.8).abs() < 1e-6);
    }
}
