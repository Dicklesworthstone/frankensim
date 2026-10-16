//! Run: cargo run -p fs-topopt --features cutfem-marquee --release
//!      --example adaptive_elastic_sdf3 -- 2 3 250000 multilevel multilevel 2
//! Arguments: rounds (1..=4), updates/round, total Krylov iterations,
//! enrichment mode, optimization mode, initial uniform level (1..=3).
//! Defaults remain two-level enrichment, Jacobi optimization, initial level 1.
//! Dimensionless slab z<0.73, independent y/z body loads, volume cap 0.5.
//! One geometric hierarchy per grid; current-density factors shared by loads.
//! Between rounds, enriched DWR marks at most two cells on the ORIGINAL tree.
//! Numerical setup products, setup applies and outer iterations are distinct.
//! No cross-grid descent, continuum-bound or wall-clock scalability claim.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::octree3::{Octree3,Octant3};
use fs_cutfem::elastic3::{ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{AdaptiveSolveOptions3,AdaptiveMultilevelOptions3,AdaptiveSolveSpace3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveControl,SolveBudget,MultiLoadOcOptions,MultiLoadOcTermination};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria,inherit_raw_densities3};
use fs_topopt::sdf3_goal::{GoalBodyLoad3,GoalRefinementOptions3,GoalPreconditioner3};
use fs_solver::op::two_level::TwoLevelBudget;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],axis:HeightAxis)->Interval {
        let d=if axis==HeightAxis::Z {1.0}else{0.0};Interval::new(d,d)
    }
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>6{return Err("usage: adaptive_elastic_sdf3 [ROUNDS [UPDATES [TOTAL_KRYLOV [ENRICHMENT_MODE [OPTIMIZATION_MODE [INITIAL_LEVEL]]]]]]".into());}
    let rounds=args.first().map_or(Ok(2),|s|s.parse::<usize>())?;
    let updates=args.get(1).map_or(Ok(3),|s|s.parse::<usize>())?;
    let budget=args.get(2).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let initial=args.get(5).map_or(Ok(1),|s|s.parse::<u8>())?;
    if !(1..=4).contains(&rounds)||!(1..=3).contains(&initial){return Err("rounds must be 1..=4 and initial level 1..=3".into());}
    let solver=args.get(3).map(String::as_str).unwrap_or("two-level");
    let preconditioner=match solver {
        "identity"=>GoalPreconditioner3::Identity,
        "jacobi"=>GoalPreconditioner3::Jacobi{max_contributions:100_000_000},
        "two-level"=>GoalPreconditioner3::TwoLevel{budget:TwoLevelBudget::default(),max_diagonal_contributions:100_000_000},
        "multilevel"=>GoalPreconditioner3::Multilevel{options:AdaptiveMultilevelOptions3::default()},
        _=>return Err("enrichment mode must be identity, jacobi, two-level, or multilevel".into()),
    };
    let optimization=args.get(4).map(String::as_str).unwrap_or("jacobi");
    if !matches!(optimization,"jacobi"|"two-level"|"multilevel") {return Err("optimization mode must be jacobi, two-level, or multilevel".into());}
    let leaf_cap=if initial==1{2048}else{8192};
    let mut tree=Octree3::uniform(initial,5,leaf_cap)?;
    let domain=HexCell::try_new([0.0;3],[1.0;3])?;
    let material=IsotropicElastic::new(1.0,0.3,1.0)?;
    let mut geometry_poll=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut geometry_poll)?;
    // Coarse geometries supply interpolation, never substituted stiffnesses.
    // Reuse their geometry across rounds and bind transfers to each new grid.
    let correction_levels:Vec<u8>=if solver=="multilevel"||optimization=="multilevel" {
        (0..initial).rev().collect()
    } else if optimization=="two-level" {vec![0]} else {Vec::new()};
    let mut corrections=Vec::with_capacity(correction_levels.len());
    for level in correction_levels {
        corrections.push(AdaptiveElasticity3::build(domain,&Octree3::uniform(level,5,leaf_cap)?,&Slab,&material,
            &|p|p[0]==0.0,ElasticityOptions3::default(),&mut geometry)?);
    }
    let correction_refs:Vec<_>=corrections.iter().collect();
    let mut solve_poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut solve_poll);
    let mut previous:Option<(Vec<Octant3>,Vec<f64>)>=None;
    let body_y=|_:[f64;3]|[0.0,-1.0,0.0];let body_z=|_:[f64;3]|[0.0,0.0,-1.0];
    println!("round,background_cells,active_cells,master_nodes,iteration,compliance,volume_fraction");
    for round in 0..rounds {
        let op=AdaptiveElasticity3::build(domain,&tree,&Slab,&material,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut geometry)?;
        let y=op.body_load(&body_y,||ControlFlow::Continue(()))?;let z=op.body_load(&body_z,||ControlFlow::Continue(()))?;
        let op=match optimization {
            "multilevel"=>AdaptiveSolveSpace3::multilevel(op,&correction_refs,AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(()))?,
            "two-level"=>AdaptiveSolveSpace3::two_level(op,corrections.last().ok_or("missing correction geometry")?,AdaptiveSolveOptions3::default(),||ControlFlow::Continue(()))?,
            _=>AdaptiveSolveSpace3::jacobi(op,100_000_000),
        };
        eprintln!("round={round} optimization_levels={:?}",op.level_sizes());
        let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());
        let raw=match &previous {
            Some((leaves,rho))=>inherit_raw_densities3(leaves,rho,study.operator().elasticity().leaves(),&mut control)?,
            None=>vec![0.5;study.cells()],
        };
        let rho=study.feasible_start(&raw,0.5,1e-8,&mut control)?;
        let loads=[LoadCase{force:&y,weight:0.3},LoadCase{force:&z,weight:0.7}];
        let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,
            MultiLoadOcOptions{max_iterations:updates,..Default::default()},&mut control);
        for row in &report.history {println!("{round},{},{},{},{},{:.17e},{:.17e}",
            tree.leaves().len(),study.cells(),study.operator().elasticity().nodes().len(),row.iteration,row.compliance,row.volume_fraction);}
        eprintln!("round={round} termination={:?} optimization_preconditioner={optimization} cumulative_krylov={} cumulative_setup_applications={} cumulative_galerkin_products={} continuum_certified=false",
            report.termination,report.work.linear_iterations,report.work.preconditioner_operator_applications,report.work.preconditioner_galerkin_products);
        if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure) {
            return Err(format!("stopped with accepted prefix: {:?}",report.evaluation_stop).into());
        }
        if round+1==rounds {break;}
        let probe_tree=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(()))?;
        let probe=AdaptiveElasticity3::build(domain,&probe_tree,&Slab,&material,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut geometry)?;
        let extra=if solver=="multilevel" {correction_refs.as_slice()}else{&[]};
        let goal=study.estimate_compliance_enrichment_with_coarse_levels(probe,extra,
            &[GoalBodyLoad3{density:&body_y,weight:0.3},GoalBodyLoad3{density:&body_z,weight:0.7}],
            &report.displacements,GoalRefinementOptions3{preconditioner,..Default::default()},&mut control)?;
        let marked=goal.mark(0.5,2,||ControlFlow::Continue(()))?;
        let dwr:f64=goal.cases.iter().zip(&goal.weights).map(|(g,w)|w*g.dwr).sum();
        let consistency:f64=goal.cases.iter().zip(&goal.weights).map(|(g,w)|w*g.coarse_space).sum();
        eprintln!("round={round} hierarchical_dwr={dwr:.9e} coarse_space={consistency:.9e} two_grid_change={:.9e} marking_fraction={:.6} target_met={} cumulative_krylov={} preconditioner={solver} cumulative_setup_applications={} cumulative_galerkin_products={}",
            goal.correction,marked.achieved_fraction,marked.target_met,goal.work.linear_iterations,goal.work.preconditioner_operator_applications,goal.work.preconditioner_galerkin_products);
        if marked.marked.is_empty(){eprintln!("no numerical refinement signal; this is not an accuracy certificate");break;}
        previous=Some((study.operator().elasticity().leaves().to_vec(),report.rho));
        tree=tree.refined(&marked.marked,||ControlFlow::Continue(()))?;
    }
    Ok(())
}
