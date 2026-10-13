//! Normalize each independent eroded load by a fixed baseline reference, then
//! minimize the smooth maximum under the actual dilated material-volume cap.
//! --minimum-material CAP_Y CAP_Z instead minimizes dilated material under
//! explicit hard eroded-compliance requirements; no cap is baseline-inferred.
//! --refine-material CAP_Y CAP_Z adds cap-wise enrichment and one selective refit.
//! --adapt-material CAP_Y CAP_Z MAX_REFITS repeats that process in bounded segments.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams, SolveBudget, SolveControl};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::robust::envelope::{LoadEnvelope3,
    optimize::{EnvelopeDesignOptions3, EnvelopeDesignStudy3}};
use fs_ascent::projected_al::ProjectedAlStop;
#[path = "worst_load_sdf3/limits.rs"]
mod limits;
#[path = "worst_load_sdf3/refinement.rs"]
mod refinement;
#[path = "worst_load_sdf3/campaign.rs"]
mod campaign;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 { p[2]-0.73 }
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)
    }
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z {1.0} else {0.0}; Interval::new(d,d)
    }
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg == "--minimum-material") { return limits::run(&args[1..]); }
    if args.first().is_some_and(|arg| arg == "--refine-material") { return refinement::run(&args[1..]); }
    if args.first().is_some_and(|arg| arg == "--adapt-material") { return campaign::run(&args[1..]); }
    if args.len()>2 { return Err("usage: worst_load_sdf3 [ACCEPTED_UPDATES [TOTAL_KRYLOV_ITERATIONS]]".into()); }
    let updates=args.first().map_or(Ok(200),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250_000),|s|s.parse::<usize>())?;
    let mut poll=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll)?;
    let op=CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3])?,[3,2,2],&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0)?,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut geometry)?;
    let y=op.body_load(&|p|[0.0,-1.0-p[1],0.0],||ControlFlow::Continue(()))?;
    let z=op.body_load(&|p|[0.0,0.0,-1.0-p[0]],||ControlFlow::Continue(()))?;
    let keys=op.cell_keys();let nodes=op.nodes().to_vec();
    let params=SimpParams {penal:2.0,beta:2.0,e_min:1e-3,..Default::default()};
    let mut study=CutDensityStudy3::new(op,0.15,params);
    let raw=vec![0.3;study.cells()];
    let options=EnvelopeDesignOptions3::default();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget {total_iterations:iterations,..Default::default()},&mut poll);
    // Explicit normalization policy: each case's own initial eroded compliance.
    // These references are measured once, charged, printed and NEVER updated.
    let baseline=study.evaluate_robust(&raw,&[LoadCase {force:&y,weight:0.5},LoadCase {force:&z,weight:0.5}],
        options.robustness,&mut control)?;
    let policy=LoadEnvelope3::new(baseline.eroded.objective.case_compliances,0.05,2)?;
    let forces=[&y[..],&z[..]];
    let mut design=EnvelopeDesignStudy3::new(&mut study,&forces,&raw,policy.clone(),options,&mut control)?;
    let outcome=design.run(updates);
    eprintln!("units=dimensionless; seeds=none; params={params:?}; options={options:?}; references={:?}; temperature={}; algebraic_smoothing_bias_bound={:.17e}; solve_budget={iterations}; geometry_work={:?}",
        policy.references(),policy.temperature(),policy.bias_bound(),geometry.work());
    eprintln!("gradient_passed={}; optimizer_outcome={outcome:?}; solve_work={:?}; optimizer_work={:?}",
        design.gradient_check().numerical.passed(),design.work(),design.optimizer_work());
    println!("iteration,envelope,actual_maximum,worst_case,dilated_volume,material_feasible");
    for row in design.history() {
        println!("{},{:.17e},{:.17e},{},{:.17e},{}",row.iteration,row.envelope,row.maximum,
            row.worst_case,row.dilated_volume,row.feasible);
    }
    let feasible=design.best_feasible().is_some();
    let point=design.best_feasible().unwrap_or_else(||design.accepted());
    let e=&point.evaluation;
    eprintln!("exported_best_feasible={feasible}; reference_values_are_hard_constraints=false; geometry_rebuilt=false; continuum_bound=false; global_optimum=false; durable_resume=false");
    println!("load_case,eroded_compliance,fixed_reference,normalized_compliance");
    for (i,(&c,&r)) in e.physical.eroded.objective.case_compliances.iter().zip(policy.references()).enumerate() {
        println!("{i},{c:.17e},{r:.17e},{:.17e}",e.envelope.ratios[i]);
    }
    println!("cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i,k) in keys.iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}",k[0],k[1],k[2],point.rho[i],
            e.physical.eroded.projected_rho[i],e.physical.nominal_rho[i],e.physical.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,eroded_uy_x,eroded_uy_y,eroded_uy_z,eroded_uz_x,eroded_uz_y,eroded_uz_z");
    for (i,p) in nodes.iter().enumerate() {
        let uy=&e.physical.eroded.objective.displacements[0];let uz=&e.physical.eroded.objective.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],p[1],p[2],uy[3*i],uy[3*i+1],uy[3*i+2],uz[3*i],uz[3*i+1],uz[3*i+2]);
    }
    let terminal=outcome?;
    if !feasible {return Err("no accepted design meets the dilated material cap; exported point is infeasible".into());}
    if !matches!(terminal.stop,ProjectedAlStop::Converged|ProjectedAlStop::IterationLimit) {
        return Err("optimizer stopped before its requested allowance; exported fields retain the best feasible prefix".into());
    }
    Ok(())
}
