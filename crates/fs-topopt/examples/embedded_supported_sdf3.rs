//! Real implicit-surface supports, with no background-box clamp.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example embedded_supported_sdf3 -- 3 250000 2
//! Arguments: updates per grid, cumulative Krylov budget, rounds (1..=3).
//! `--motion` instead solves an affine prescribed-motion experiment and prints
//! full nodal fields. Its linear iteration cap is 10000; it does not optimize.
//!
//! Dimensionless slab 0.17<x<0.83 in the unit box. The LEFT implicit face
//! supplies the support; right-face compression/shear are independent loads.
//! Homogeneous supports are used for optimization. Nonzero motion requires
//! recomputing the density-dependent lifting and is not a fixed-load OC problem.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::{Octant3, Octree3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{MultiLoadOcOptions, MultiLoadOcTermination, SimpParams, SolveBudget, SolveControl};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3, controlled_sdf3_optimality_criteria, inherit_raw_densities3};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalReferenceLoad3, GoalRefinementOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],axis:HeightAxis)->Interval {
        if axis==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(0.17,0.17)-Interval::new(0.83,0.83)}
        else {Interval::new(0.0,0.0)}
    }
}
fn build(tree:&Octree3,nu:f64,control:&mut QuadratureControl3<'_>) -> Result<AdaptiveElasticity3,Box<dyn std::error::Error>> {
    Ok(AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,tree,&Slab,
        &IsotropicElastic::new(1.0,nu,1.0)?,&|_|false,&|_,n|n[0]<0.0,
        ElasticityOptions3::default(),Default::default(),Default::default(),control)?)
}
fn compression(_: [f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[-1.0,0.0,0.0]}else{[0.0;3]}}
fn shear(_: [f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.0,0.0,-1.0]}else{[0.0;3]}}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let mut gp=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let mut tree=Octree3::uniform(1,5,4096)?;
    if args.len()==1 && args[0]=="--motion" {
        let op=build(&tree,0.0,&mut geometry)?;
        let f=op.surface_load(&|_,n|if n[0]>0.0{[1.0,0.0,0.0]}else{[0.0;3]},||ControlFlow::Continue(()))?.rhs;
        let lifting=op.prescribed_displacement_load(&|_,_|[0.125,0.0,0.0],||ControlFlow::Continue(()))?;
        let rhs:Vec<_>=f.iter().zip(lifting).map(|(f,g)|f+g).collect();
        let solved=op.solve_controlled(&rhs,1e-10,10000,16,|_|ControlFlow::Continue(()))?;
        let physical=op.physical_displacements(solved.coefficients())?;
        let external_work:f64=f.iter().zip(solved.coefficients()).map(|(f,u)|f*u).sum();
        println!("x,y,z,ux,uy,uz");
        for (p,u) in op.physical_nodes().iter().zip(physical.chunks_exact(3)) {
            println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",p[0],p[1],p[2],u[0],u[1],u[2]);
        }
        eprintln!("external_work={external_work:.17e}; augmented_rhs_work={:.17e}; iterations={}; true_relative_residual={:.9e}; expected_ux=x-0.17+0.125; actuator_work_not_claimed=true",
            solved.compliance(),solved.iterations(),solved.residual_claim().euclidean().expect("true residual"));
        return Ok(());
    }
    if args.len()>3 {return Err("usage: embedded_supported_sdf3 [UPDATES [TOTAL_KRYLOV [ROUNDS]]] or --motion".into());}
    let updates=args.first().map_or(Ok(3),|s|s.parse::<usize>())?;
    let budget=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let rounds=args.get(2).map_or(Ok(1),|s|s.parse::<usize>())?;
    if !(1..=3).contains(&rounds) {return Err("rounds must lie in 1..=3".into());}
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut poll);
    let mut previous:Option<(Vec<Octant3>,Vec<f64>)>=None;
    println!("round,active_cells,iteration,compliance,volume_fraction");
    for round in 0..rounds {
        let op=build(&tree,0.3,&mut geometry)?;
        let x=op.surface_load(&compression,||ControlFlow::Continue(()))?.rhs;
        let z=op.surface_load(&shear,||ControlFlow::Continue(()))?.rhs;
        eprintln!("round={round}; embedded_support_area={:.9e}; box_clamped_nodes={}",op.embedded_dirichlet_area().unwrap_or(0.0),op.fixed().iter().filter(|b|**b).count());
        let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
        let raw=match &previous {
            None=>vec![0.5;study.cells()],
            Some((leaves,rho))=>inherit_raw_densities3(leaves,rho,study.operator().elasticity().leaves(),&mut control)?,
        };
        let rho=study.feasible_start(&raw,0.5,1e-8,&mut control)?;
        let report=controlled_sdf3_optimality_criteria(&mut study,&[LoadCase{force:&x,weight:0.5},LoadCase{force:&z,weight:0.5}],
            &rho,MultiLoadOcOptions{max_iterations:updates,..Default::default()},&mut control);
        for row in &report.history {println!("{round},{},{},{:.17e},{:.17e}",study.cells(),row.iteration,row.compliance,row.volume_fraction);}
        eprintln!("round={round}; termination={:?}; cumulative_krylov={}; continuum_certified=false",report.termination,report.work.linear_iterations);
        if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure) {
            return Err(format!("stopped with accepted prefix: {:?}",report.evaluation_stop).into());
        }
        if round+1==rounds {break;}
        let probe_tree=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(()))?;
        let probe=build(&probe_tree,0.3,&mut geometry)?;
        let goal=study.estimate_reference_compliance_enrichment(probe,&[
            GoalReferenceLoad3{load:ReferenceLoad3::traction(&compression),weight:0.5},
            GoalReferenceLoad3{load:ReferenceLoad3::traction(&shear),weight:0.5}],&report.displacements,
            GoalRefinementOptions3{preconditioner:GoalPreconditioner3::Jacobi{max_contributions:100_000_000},..Default::default()},&mut control)?;
        let mark=goal.mark(0.5,2,||ControlFlow::Continue(()))?;
        eprintln!("round={round}; two_grid_difference={:.9e}; marking_fraction={:.6}; target_met={}; no_cross_grid_descent_claim=true",goal.correction,mark.achieved_fraction,mark.target_met);
        if mark.marked.is_empty() {break;}
        previous=Some((study.operator().elasticity().leaves().to_vec(),report.rho));
        tree=tree.refined(&mark.marked,||ControlFlow::Continue(()))?;
    }
    Ok(())
}
