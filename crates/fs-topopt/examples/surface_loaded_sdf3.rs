//! Optimize and selectively refine a surface-loaded implicit 3-D cantilever.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example surface_loaded_sdf3 -- 3 250000 2
//! Arguments: accepted updates per round, cumulative Krylov iterations,
//! optional rounds (1..=4, default 1). The previous two-argument run is retained.
//! Dimensionless unit box; solid x < 0.7+0.1*y*(1-y); clamp x=0. Two independent
//! dead loads: inward pressure 1+z and uniform -z surface traction, weights 1/2.
//! Surface laws are reintegrated on each new grid, never replaced by nodal forces.
//! Enrichment inherits physical stiffness; the next optimizer inherits only raw
//! densities, restores its own volume, and solves a fresh baseline. Numerical
//! DWR differences are not continuum bounds; no cross-grid descent is claimed.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::{Octree3,Octant3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl,MultiLoadOcOptions,MultiLoadOcTermination};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria,inherit_raw_densities3};
use fs_topopt::sdf3_goal::{GoalReferenceLoad3,GoalRefinementOptions3,GoalPreconditioner3};
struct CurvedEnd;
impl CutSdf3 for CurvedEnd {
    fn value(&self,p:[f64;3])->f64{p[0]-0.7-0.1*p[1]*(1.0-p[1])}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{
        let y=Interval::new(lo[1],hi[1]);
        Interval::new(lo[0],hi[0])-Interval::new(0.7,0.7)
            -Interval::new(0.1,0.1)*y*(Interval::new(1.0,1.0)-y)
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval{
        match a {
            HeightAxis::X=>Interval::new(1.0,1.0),
            HeightAxis::Y=>Interval::new(-0.1,-0.1)*(Interval::new(1.0,1.0)-Interval::new(2.0,2.0)*Interval::new(lo[1],hi[1])),
            HeightAxis::Z=>Interval::new(0.0,0.0),
        }
    }
}
fn build(tree:&Octree3,geometry:&mut QuadratureControl3<'_>)->Result<AdaptiveElasticity3,Box<dyn std::error::Error>>{
    Ok(AdaptiveElasticity3::build_with_surface(HexCell::try_new([0.0;3],[1.0;3])?,tree,&CurvedEnd,
        &IsotropicElastic::new(1.0,0.3,1.0)?,&|p|p[0]==0.0,ElasticityOptions3::default(),Default::default(),geometry)?)
}
fn main()->Result<(),Box<dyn std::error::Error>>{
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>3{return Err("usage: surface_loaded_sdf3 [UPDATES [TOTAL_KRYLOV_ITERATIONS [ROUNDS]]]".into());}
    let updates=args.first().map_or(Ok(3),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let rounds=args.get(2).map_or(Ok(1),|s|s.parse::<usize>())?;
    if !(1..=4).contains(&rounds){return Err("rounds must lie in 1..=4".into());}
    let mut gp=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let mut sp=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:iterations,..Default::default()},&mut sp);
    let mut tree=Octree3::uniform(1,4,4096)?;
    let mut previous:Option<(Vec<Octant3>,Vec<f64>)>=None;
    let pressure_law=|p:[f64;3]|1.0+p[2];
    let shear_law=|_:[f64;3],_:[f64;3]|[0.0,0.0,-1.0];
    let laws=[GoalReferenceLoad3{load:ReferenceLoad3::pressure(&pressure_law),weight:0.5},
        GoalReferenceLoad3{load:ReferenceLoad3::traction(&shear_law),weight:0.5}];
    println!("round,iteration,compliance,pressure_compliance,shear_compliance,volume_fraction");
    for round in 0..rounds {
        let op=build(&tree,&mut geometry)?;
        let pressure=op.pressure_load(&pressure_law,||ControlFlow::Continue(()))?;
        let shear=op.surface_load(&shear_law,||ControlFlow::Continue(()))?;
        eprintln!("round={round}; reference_area={:.17e}; pressure_resultant={:?}; pressure_moment={:?}; shear_resultant={:?}; surface_is_reference=true",
            pressure.area,pressure.resultant,pressure.moment,shear.resultant);
        let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
        let raw=match &previous {
            Some((leaves,rho))=>inherit_raw_densities3(leaves,rho,study.operator().elasticity().leaves(),&mut control)?,
            None=>vec![0.5;study.cells()],
        };
        let rho=study.feasible_start(&raw,0.5,1e-8,&mut control)?;
        let report=controlled_sdf3_optimality_criteria(&mut study,
            &[LoadCase{force:&pressure.rhs,weight:0.5},LoadCase{force:&shear.rhs,weight:0.5}],&rho,
            MultiLoadOcOptions{max_iterations:updates,..Default::default()},&mut control);
        for r in &report.history{println!("{round},{},{:.17e},{:.17e},{:.17e},{:.17e}",r.iteration,r.compliance,r.case_compliances[0],r.case_compliances[1],r.volume_fraction);}
        eprintln!("round={round}; termination={:?}; cumulative_krylov={}; geometry_fixed_within_optimization=true; continuum_certified=false",report.termination,report.work.linear_iterations);
        if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure){
            return Err(format!("stopped with accepted prefix: {:?}",report.evaluation_stop).into());
        }
        if round+1==rounds{break;}
        let probe_tree=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(()))?;
        let estimate=study.estimate_reference_compliance_enrichment(build(&probe_tree,&mut geometry)?,&laws,&report.displacements,
            GoalRefinementOptions3{preconditioner:GoalPreconditioner3::Jacobi{max_contributions:100_000_000},..Default::default()},&mut control)?;
        let marks=estimate.mark(0.5,2,||ControlFlow::Continue(()))?;
        eprintln!("round={round}; two_grid_change={:.17e}; marking_fraction={:.9}; target_met={}; marks={:?}; cumulative_krylov={}; cross_grid_descent_claimed=false",
            estimate.correction,marks.achieved_fraction,marks.target_met,marks.marked,estimate.work.linear_iterations);
        if marks.marked.is_empty(){eprintln!("no numerical refinement signal; not an accuracy certificate");break;}
        previous=Some((study.operator().elasticity().leaves().to_vec(),report.rho));
        tree=tree.refined(&marks.marked,||ControlFlow::Continue(()))?;
    }
    Ok(())
}
