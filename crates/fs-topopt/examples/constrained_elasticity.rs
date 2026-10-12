//! Per-element stiffness sizing under two independent dead-load cases.
//! Run: cargo run --release -p fs-topopt --features newton-design --example constrained_elasticity
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStudy};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::DensityElasticity;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget, ComplianceNewtonModel};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let (mesh,positions)=fs_feec::kuhn_cube(2);
    let elasticity=DensityElasticity::new(&mesh,&positions,1.0,0.29,&|p|p[0]==0.0);
    let n=elasticity.cells();let mut vertical=vec![0.0;elasticity.n()];let mut lateral=vertical.clone();
    for (i,p) in positions.iter().enumerate() {
        if p[0]==1.0 {vertical[3*i+2]=-0.2;lateral[3*i+1]=0.13*(1.0+p[2]);}
    }
    let model=ComplianceNewtonModel::new(&elasticity,"elastic compliance",
        &[LoadCase {force:&vertical,weight:0.7},LoadCase {force:&lateral,weight:1.3}],
        ComplianceSolveOptions {tolerance:1e-11,max_iterations:2000},
        ComplianceNewtonBudget {solve_reservations:20_000,linear_iterations:2_000_000})?;
    let pool=TilePool::for_host(1,37);let gate=CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| -> Result<(),Box<dyn std::error::Error>> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:37,kernel_id:82,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let initial=vec![0.0;n];let baseline=model.value_gradient(&initial,Some(&cx))?.value;
        let mut b=ProblemBuilder::new();let var=b.var("element log stiffness",Manifold::Rn {dim:n as u32},Dims::NONE)?;
        let r=b.var_ref(var)?;let node=b.pde_residual("elastic compliance",var,true,COMPLIANCE_DIMS)?;
        let scale=b.konst(baseline,COMPLIANCE_DIMS)?;let normalized=b.div(node,scale)?;
        let mut amount=b.konst(0.0,Dims::NONE)?;let mut regularizer=b.konst(0.0,Dims::NONE)?;
        for e in 0..n {let x=b.component(r,e as u32)?;let m=b.exp(x)?;amount=b.add(amount,m)?;
            let square=b.mul(x,x)?;regularizer=b.add(regularizer,square)?;}
        let count=b.konst(n as f64,Dims::NONE)?;let mean=b.div(amount,count)?;
        let ridge=b.konst(0.001/n as f64,Dims::NONE)?;let regularizer=b.mul(ridge,regularizer)?;
        let objective=b.add(normalized,regularizer)?;b.objective(objective,Sense::Minimize,1.0)?;
        let allowance=b.konst(1.5,Dims::NONE)?;let residual=b.sub(mean,allowance)?;
        b.constraint(residual,ConstraintKind::LeZero,"mean stiffness multiplier <= 1.5")?;
        let problem=b.finish();let oracle=ReverseProblem::new_with_physics(&problem,
            ReverseLimits {max_nodes:4096,max_scalar_slots:65536},&[PhysicsBinding {node,model:&model}],Some(&cx))?;
        let mut study=AugmentedNewtonStudy::new(&oracle,&initial,AugmentedNewtonOptions::default(),Some(&cx))?;
        let report=study.run(1e-5,150,AugmentedNewtonBudget {
            objective_evals:300,hessian_products:1000,multiplier_updates:50},Some(&cx))?;
        let final_compliance=model.value_gradient(&report.x,Some(&cx))?.value;
        println!("fixed-mesh numerical stiffness sizing; not a topology/stress/continuum certificate");
        println!("parameters={n} independent_loads=2 baseline_compliance_J={baseline:.9e}");
        println!("stop={:?} objective={:.9e} compliance_J={final_compliance:.9e}",report.stop,report.objective);
        println!("original_constraints={:?} candidate_multipliers={:?} kkt={:?}",report.constraints,report.multipliers,report.kkt);
        println!("physics_work={:?} optimizer_work={:?}",model.work(),study.inner().snapshot());
        Ok(())
    })
}
