//! Exact-IR Hessians plus declared Jacobi preconditioning for a stiff quadratic.
//! This is synthetic dimensionless algebra, not a calibrated physical model.
use fs_ascent::{ReverseNewtonStop, ReverseNewtonStudy, StopReason, StopRule};
use fs_ascent::trust::KrylovOptions;
use fs_qty::Dims;
use fs_opt::{ Manifold, ProblemBuilder, ReverseProblem, Sense, reverse::ReverseLimits};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let weights: Vec<f64> = (0..16).map(|i| (2.0f64).powi(i)).collect();
    let mut builder = ProblemBuilder::new();
    let x = builder.var("design", Manifold::Rn { dim: 16 }, Dims::NONE)?;
    let values = builder.var_ref(x)?;
    let mut objective = builder.konst(0.0, Dims::NONE)?;
    for (i, weight) in weights.iter().enumerate() {
        let x_i = builder.component(values, i as u32)?;
        let square = builder.powi(x_i, 2)?;
        let stiffness = builder.konst(0.5 * weight, Dims::NONE)?;
        let energy = builder.mul(stiffness, square)?;
        objective = builder.add(objective, energy)?;
    }
    builder.objective(objective, Sense::Minimize, 1.0)?;
    let problem = builder.finish();
    let oracle = ReverseProblem::new(&problem, ReverseLimits {
        max_nodes: 1000, max_scalar_slots: 10_000,
    })?;
    let mut study = ReverseNewtonStudy::new(&oracle, &[0.125; 16], None)?;
    // The local cap controls how accurately to solve one Newton subproblem.
    // The separate cumulative limit still needs a final model HVP as well.
    study.set_krylov_options(KrylovOptions::new(1e-8, 1)?);
    let mut inverse = |_: &[f64], r: &[f64]| {
        r.iter().zip(&weights).map(|(r, diagonal)| r / diagonal).collect()
    };
    let report = study.run_preconditioned(&mut inverse, &StopRule::GradNorm(1e-10),
        1, 2, 1, None)?;
    println!("stop={:?} f={:.6e} grad_inf={:.6e} objective_attempts={} hessian_products={} inverse_actions={}",
        report.stop, report.solution.f, report.solution.grad_norm,
        report.solution.evals, report.solution.hv_evals, report.preconditioner_evals);
    if report.stop != ReverseNewtonStop::Stopped(StopReason::GradNorm) {
        return Err("the declared diagonal inverse did not converge within the exact allowance".into());
    }
    Ok(())
}
