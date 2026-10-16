//! Nonconvex, dimensionless quartic chain solved without callbacks or a Hessian.
//! Illustrative algebra, not an identified material or mesh-based physics solve.
use fs_ascent::{ReverseNewtonStop, ReverseNewtonStudy, StopReason, StopRule};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{EvalLimit, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

fn problem(n: u32) -> Problem {
    let mut b=ProblemBuilder::new();let v=b.var("chain",Manifold::Rn{dim:n},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let one=b.konst(1.0,Dims::NONE).unwrap();
    let quarter=b.konst(0.25,Dims::NONE).unwrap();let half=b.konst(0.5,Dims::NONE).unwrap();
    for i in 0..n {
        let x=b.component(r,i).unwrap();let xx=b.powi(x,2).unwrap();let a=b.sub(xx,one).unwrap();
        let square=b.powi(a,2).unwrap();let potential=b.mul(quarter,square).unwrap();
        b.objective(potential,Sense::Minimize,1.0).unwrap();
        if i>0 {
            let previous=b.component(r,i-1).unwrap();let d=b.sub(x,previous).unwrap();
            let square=b.powi(d,2).unwrap();let coupling=b.mul(half,square).unwrap();
            b.objective(coupling,Sense::Minimize,1.0).unwrap();
        }
    }
    b.set_eval_limit(EvalLimit::Limited(NonZeroU64::new(100).unwrap()));b.finish()
}
fn run() -> Result<(),Box<dyn std::error::Error>> {
    let p=problem(64);let oracle=ReverseProblem::new(&p,ReverseLimits{max_nodes:10_000,max_scalar_slots:100_000})?;
    let mut study=ReverseNewtonStudy::new(&oracle,&vec![0.25;64],None)?;
    let report=study.run(&StopRule::GradNorm(1e-9),80,1000,None)?;
    println!("stop={:?} f={:.12e} grad_inf={:.6e} evaluations={} hessian_products={} negative_curvature_hits={}",
        report.stop,report.solution.f,report.solution.grad_norm,report.solution.evals,
        report.solution.hv_evals,report.solution.negative_curvature_hits);
    println!("first_coordinate={:.12} last_coordinate={:.12}",report.solution.x[0],report.solution.x[63]);
    if report.stop!=ReverseNewtonStop::Stopped(StopReason::GradNorm) { return Err("study stopped before gradient convergence".into()); }
    Ok(())
}
fn main() {
    if let Err(error)=run(){eprintln!("reverse Newton study: {error}");std::process::exit(1);}
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn true_negative_curvature_escapes_to_the_known_chain_minimum() {
        let p=problem(64);let oracle=ReverseProblem::new(&p,ReverseLimits{max_nodes:10_000,max_scalar_slots:100_000}).unwrap();
        let mut study=ReverseNewtonStudy::new(&oracle,&vec![0.25;64],None).unwrap();
        let report=study.run(&StopRule::GradNorm(1e-9),80,1000,None).unwrap();
        assert_eq!(report.stop,ReverseNewtonStop::Stopped(StopReason::GradNorm));
        assert!(report.solution.negative_curvature_hits>0);
        assert!(report.solution.x.iter().all(|x|(x-1.0).abs()<1e-8));assert!(report.solution.f<1e-16);
    }
}
