//! Synthetic normalized design: stage-dependent stiffness and an active inequality.
//! cargo run --release -p fs-ascent --example preconditioned_continuation
//! This diagonal fixture demonstrates the API; it is not a PDE speed benchmark.
use fs_ascent::newton_continuation::{NewtonContinuation, NewtonStage};
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions};
use fs_ascent::trust::KrylovOptions;
use fs_opt::reverse::ReverseLimits;
use fs_opt::{AugmentedLagrangian, ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let n=12;
    let diagonals:Vec<Vec<f64>>=(0..3).map(|stage|(0..n).map(|i|(4.0+stage as f64).powi(i as i32)).collect()).collect();
    let mut problems=Vec::new();
    for (stage,d) in diagonals.iter().enumerate() {
        let mut b=ProblemBuilder::new();let var=b.var("design",Manifold::Rn {dim:n as u32},Dims::NONE)?;
        let x=b.var_ref(var)?;let mut total=b.konst(0.0,Dims::NONE)?;
        for (i,&d) in d.iter().enumerate() {
            let xi=b.component(x,i as u32)?;
            let t=b.konst(if i==0 {2.0} else {0.005*(stage+1) as f64},Dims::NONE)?;
            let r=b.sub(xi,t)?;let square=b.powi(r,2)?;let weight=b.konst(0.5*d,Dims::NONE)?;
            let term=b.mul(weight,square)?;total=b.add(total,term)?;
        }
        b.objective(total,Sense::Minimize,1.0)?;
        let first=b.component(x,0)?;let one=b.konst(1.0,Dims::NONE)?;let c=b.sub(first,one)?;
        b.constraint(c,ConstraintKind::LeZero,"first coordinate allowance")?;
        problems.push(b.finish());
    }
    let limits=ReverseLimits {max_nodes:1024,max_scalar_slots:8192};
    let oracles=problems.iter().map(|p|ReverseProblem::new(p,limits)).collect::<Result<Vec<_>,_>>()?;
    let stages:Vec<_>=oracles.iter().map(|oracle|NewtonStage {oracle,tolerance:1e-6,
        options:AugmentedNewtonOptions::default(),krylov:KrylovOptions::default()}).collect();
    let mut study=NewtonContinuation::new(&stages,&vec![0.0;n],None)?;
    let mut inverse=|stage:usize,_:&[f64],r:&[f64],merit:&AugmentedLagrangian<'_>,c:&[f64]| {
        let active=merit.multipliers()[0]+merit.penalty()*c[0]>0.0;
        r.iter().zip(&diagonals[stage]).enumerate().map(|(i,(&r,&d))|
            r/(d+if i==0 && active {merit.penalty()} else {0.0})).collect()
    };
    let report=study.run_preconditioned(&mut inverse,500,AugmentedNewtonBudget {
        objective_evals:2000,hessian_products:20_000,multiplier_updates:300},5000,None)?;
    println!("stop={:?}; stage={}; lifetime_work={:?}",report.stop,report.stage,report.work);
    if let Some(solution)=report.solution {
        println!("original_objective={:.12e}; kkt={:?}; x={:?}",solution.objective,solution.kkt,solution.x);
    }
    Ok(())
}
