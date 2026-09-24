//! Raw 3-D implicit-domain density optimization, with no geometric rebuilding.
//! cargo run -p fs-topopt --features cutfem-marquee --release --example elastic_sdf3 -- 5 250000
//! Arguments: accepted-update cap, cumulative Krylov-iteration cap. Prints a
//! CSV history followed by accepted cell densities; interrupted work exits
//! nonzero after printing its accepted prefix. All units are dimensionless.
//! Prefix --continuation to run (p,beta)=(1,1),(2,2),(3,8), restoring
//! volume and checking compliance/volume gradients at every stage baseline.
//! Prefix --adaptive-continuation for the same schedule with real DWR-guided
//! octree background refinement and freshly checked transfers between stages.
//! The curved design domain is z < 0.7 + 0.1*x*(1-x) inside the unit box.
//! This field is implicit, not an exact distance function. The two dead body
//! loads are independent. No continuum, manufacturing or optimality claim.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{ContinuationTermination,GradientCheckOptions,MultiLoadOcOptions,MultiLoadOcTermination,SimpParams,SolveBudget,SolveControl};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria};
use fs_topopt::sdf3::continuation::controlled_gradient_checked_sdf3_continuation;
#[path = "elastic_sdf3/adaptive.rs"]
mod adaptive;

struct CurvedCantilever;
impl CutSdf3 for CurvedCantilever {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.7-0.1*p[0]*(1.0-p[0])}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);let one=Interval::new(1.0,1.0);
        Interval::new(lo[2],hi[2])-Interval::new(0.7,0.7)-Interval::new(0.1,0.1)*x*(one-x)
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        match a {
            HeightAxis::X=>Interval::new(-0.1,-0.1)*(Interval::new(1.0,1.0)-Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])),
            HeightAxis::Y=>Interval::new(0.0,0.0),HeightAxis::Z=>Interval::new(1.0,1.0),
        }
    }
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let mut args:Vec<String>=std::env::args().skip(1).collect();
    let continuation=args.first().is_some_and(|arg|arg=="--continuation");
    let adapt=args.first().is_some_and(|arg|arg=="--adaptive-continuation");
    if continuation || adapt {args.remove(0);}
    if args.len()>2 {return Err("usage: elastic_sdf3 [--continuation|--adaptive-continuation] [UPDATES_PER_STAGE [TOTAL_KRYLOV_ITERATIONS]]".into());}
    let updates=args.first().map_or(Ok(5),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250_000),|s|s.parse::<usize>())?;
    if adapt {return adaptive::run(updates,iterations);}
    let mut poll=|_|ControlFlow::Continue(());
    let mut quadrature=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll)?;
    let op=CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3])?,[4,2,2],&CurvedCantilever,
        &IsotropicElastic::new(1.0,0.3,1.0)?,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut quadrature)?;
    let y=op.body_load(&|_|[0.0,-1.0,0.0],||ControlFlow::Continue(()))?;
    let z=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(()))?;
    let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());
    let rho=vec![0.5;study.cells()];
    let geometry_work=quadrature.work();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget {total_iterations:iterations,..Default::default()},&mut poll);
    let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
    let options=MultiLoadOcOptions {max_iterations:updates,..Default::default()};
    let (report,finished)=if continuation {
        let schedule=[(1.0,1.0),(2.0,2.0),(3.0,8.0)]
            .map(|(penal,beta)|SimpParams {penal,beta,..Default::default()});
        let stages=controlled_gradient_checked_sdf3_continuation(&mut study,&loads,&rho,&schedule,
            options,GradientCheckOptions::default(),&mut control);
        println!("stage,penal,beta,iteration,compliance,volume_fraction,max_change");
        for stage in &stages.stages {
            for row in &stage.history {println!("{},{},{},{},{:.17e},{:.17e},{:.17e}",
                stage.stage,stage.params.penal,stage.params.beta,row.iteration,row.compliance,row.volume_fraction,row.max_change);}
            let check=stage.gradient_check.as_ref().expect("gradient-gated stage");
            for probe in &check.probes {eprintln!("stage={} direction={:?} gradient_relative_error={:.9e} volume_gradient_relative_error={:.9e} passed={}",
                stage.stage,probe.direction,probe.compliance_relative_error,probe.volume_relative_error,check.passed());}
            eprintln!("stage={} incoming_volume={:.9e} restoration_scale={:.9e} stage_stop={:?}",
                stage.stage,stage.incoming_volume_fraction,stage.restoration_scale,stage.termination);
        }
        eprintln!("continuation_stop={:?}; stopped_stage={:?}; evaluation_stop={:?}; rejected_gradient_check={:?}; cumulative_linear_iterations={}; cross_model_descent_claimed=false",
            stages.termination,stages.stopped_stage,stages.evaluation_stop,stages.rejected_gradient_check,stages.work.linear_iterations);
        (stages.last,stages.termination==ContinuationTermination::ScheduleComplete)
    } else {
        let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,options,&mut control);
        println!("iteration,compliance,volume_fraction,max_change");
        for row in &report.history {println!("{},{:.17e},{:.17e},{:.17e}",row.iteration,row.compliance,row.volume_fraction,row.max_change);}
        let finished=!matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure);
        (Some(report),finished)
    };
    let Some(report)=report else {return Err("study stopped before any checked equilibrium baseline; no design is exported".into());};
    eprintln!("termination={:?}; cut_cells={}; geometry_points={}; linear_iterations={}; geometry_rebuilt=false; optimality_certified=false",
        report.termination,study.cells(),geometry_work.points,report.work.linear_iterations);
    if !report.history.is_empty() {
        println!("cell_x,cell_y,cell_z,raw_density,projected_density");
        for ((key,raw),physical) in study.operator().cell_keys().iter().zip(&report.rho).zip(&report.projected_rho) {
            println!("{},{},{},{raw:.17e},{physical:.17e}",key[0],key[1],key[2]);
        }
    }
    if !finished {
        return Err("study interrupted or continuation gradient gate failed; exported fields describe the retained prefix".into());
    }
    Ok(())
}
