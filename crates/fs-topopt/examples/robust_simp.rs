//! Eroded compliance / dilated volume with 48 independent raw-density variables.
//! cargo run --release -p fs-topopt --features newton-design --example robust_simp
//! Synthetic fixed-mesh linear elasticity; no manufacturing or continuum certificate.
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStudy};
use fs_ascent::trust::KrylovOptions;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsModel};
use fs_opt::reverse::physics::box_map::LogisticBox;
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::SimpNewtonModel;
use fs_topopt::{DensityElasticity, DensityFilter, RobustPipeline, SimpParams};

#[path = "robust_simp/continuation.rs"]
mod continuation;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--continuation"] { return continuation::run(); }
    if !args.is_empty() {
        return Err(std::io::Error::other("usage: robust_simp [--continuation]").into());
    }
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions {p[0] *= p[0];}
    let elasticity = DensityElasticity::new(&mesh,&positions,1.0,0.29,&|p|p[0]==0.0);
    let n = elasticity.cells();
    let volumes:Vec<_> = fs_feec::element_geometry(&mesh,&positions).vol_signed.iter().map(|v|v.abs()).collect();
    let mut a = vec![0.0;elasticity.n()]; let mut b = a.clone();
    for (i,p) in positions.iter().enumerate() {
        if p[0]==1.0 {a[3*i+2]=-0.2;b[3*i+1]=0.13*(1.0+p[2]);}
    }
    let loads = [LoadCase {force:&a,weight:0.7},LoadCase {force:&b,weight:1.3}];
    let pipeline = RobustPipeline {filter:DensityFilter::new(&mesh,&positions,0.15),
        params:SimpParams {e_min:0.01,penal:3.0,beta:2.0,eta:0.5},eta_offset:0.1};
    let physics = SimpNewtonModel::new_robust(&pipeline,&elasticity,&loads,&volumes,
        ComplianceSolveOptions {tolerance:1e-11,max_iterations:4000},
        ComplianceNewtonBudget {solve_reservations:100_000,linear_iterations:10_000_000})?;
    let raw_c = physics.compliance("eroded-compliance")?;
    let raw_v = physics.volume("dilated-material")?;
    let chart = LogisticBox::new(&vec![0.05;n],&vec![0.95;n],Dims::NONE)?;
    let c = chart.bind(&raw_c)?; let v = chart.bind(&raw_v)?;
    let pool = TilePool::for_host(1,74); let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| -> Result<(),Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate,arena,StreamKey {seed:74,kernel_id:95,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let start = chart.encode(&vec![0.5;n],Some(&cx))?;
        let baseline = c.value_gradient(&start,Some(&cx))?.value;
        let mut builder = ProblemBuilder::new();
        let x = builder.var("bounded density",Manifold::Rn {dim:n as u32},Dims::NONE)?;
        let cn = builder.pde_residual("eroded-compliance",x,true,COMPLIANCE_DIMS)?;
        let vn = builder.pde_residual("dilated-material",x,true,Dims::NONE)?;
        let scale = builder.konst(baseline,COMPLIANCE_DIMS)?; let objective = builder.div(cn,scale)?;
        builder.objective(objective,Sense::Minimize,1.0)?;
        // This allowance is dilated physical material, NOT nominal material.
        let cap = 0.6; let bound = builder.konst(cap,Dims::NONE)?; let residual = builder.sub(vn,bound)?;
        builder.constraint(residual,ConstraintKind::LeZero,"dilated volume fraction")?;
        let problem = builder.finish();
        let oracle = ReverseProblem::new_with_physics(&problem,ReverseLimits {max_nodes:256,max_scalar_slots:8192},
            &[PhysicsBinding {node:cn,model:&c},PhysicsBinding {node:vn,model:&v}],None)?;
        let mut study = AugmentedNewtonStudy::new(&oracle,&start,AugmentedNewtonOptions::default(),Some(&cx))?;
        study.set_krylov_options(KrylovOptions::new(0.05,8)?);
        let report = study.run(1e-5,300,AugmentedNewtonBudget {
            objective_evals:2500,hessian_products:10_000,multiplier_updates:100},Some(&cx))?;
        let raw = chart.decode(&report.x,Some(&cx))?;
        // Independent first-order, reassembled physical check, not a cache read.
        let mut reference = DensityElasticity::new(&mesh,&positions,1.0,0.29,&|p|p[0]==0.0);
        let mut fresh_c = 0.0;
        for load in &loads {
            fresh_c += load.weight*pipeline.eroded_compliance_and_gradient(&mut reference,&raw,load.force).0;
        }
        let fields = pipeline.three_fields(&raw);
        let total:f64 = volumes.iter().sum();
        let volume = |rho:&[f64]| rho.iter().zip(&volumes).map(|(r,v)|r*v).sum::<f64>()/total;
        let (ve,vn,vd) = (volume(&fields.eroded),volume(&fields.nominal),volume(&fields.dilated));
        if !fresh_c.is_finite() || !vd.is_finite()
            || (fresh_c/baseline-report.objective).abs()>1e-7*(1.0+report.objective.abs())
            || (vd-cap-report.constraints[0]).abs()>1e-10 {
            return Err(std::io::Error::other("fresh robust physics disagreed with Newton report").into());
        }
        println!("stop={:?}; chart_kkt={:?}",report.stop,report.kkt);
        println!("normalized_eroded_compliance={:.12}; volumes(eroded,nominal,dilated)=({ve:.9},{vn:.9},{vd:.9})",fresh_c/baseline);
        println!("dilated_allowance={cap}; thresholds={:?}; physics_work={:?}",physics.projection_thresholds(),physics.work());
        println!("raw_densities={raw:?}");
        println!("A numerical candidate, not a binary/feature-size certificate. A budget or stall stop is not convergence.");
        Ok(())
    })
}
