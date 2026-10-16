//! cargo run -p fs-ascent --example reverse_study
//! Solve a 256-coordinate IR objective with three budgeted reverse evaluations.
use fs_ascent::{ReverseStudy, StopRule};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{EvalLimit, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut builder = ProblemBuilder::new();
    let variable = builder.var("design", Manifold::Rn { dim: 256 }, Dims::NONE)?;
    let reference = builder.var_ref(variable)?;
    let objective = builder.norm_sq(reference)?;
    builder.objective(objective, Sense::Minimize, 1.0)?;
    builder.set_eval_limit(EvalLimit::Limited(NonZeroU64::new(3).expect("positive constant")));
    let problem = builder.finish();
    let oracle = ReverseProblem::new(&problem, ReverseLimits {
        max_nodes: 1024, max_scalar_slots: 16_384,
    })?;
    let mut study = ReverseStudy::new(&oracle, &[1.0; 256], 7, None)?;
    let report = study.run(&StopRule::GradNorm(0.0), 100, None)?;
    println!("stop={:?}, objective={}, gradient={}, evaluations={}, accepted_steps={}",
        report.reason, report.f, report.grad_norm, report.evals, report.iters);
    assert_eq!(report.evals, 3);
    assert_eq!(study.optimizer().x, vec![0.0; 256]);
    Ok(())
}
