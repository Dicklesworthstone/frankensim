//! Run: cargo run -p fs-topopt --features cutfem-marquee --release
//!      --example adaptive_elastic_sdf3 -- 2 3 250000
//! Arguments: refinement rounds (1..=4), updates/round, total Krylov iterations.
//! Dimensionless slab z<0.73, two independent body loads, projected volume cap 0.5.
//! Marks the two highest discrete-energy-density cells between rounds. This is
//! an energy heuristic, NOT DWR. Histories have fresh baselines: neither cross-
//! mesh compliance descent nor mesh-converged accuracy is asserted. Only raw
//! densities transfer; every refined mesh re-integrates and re-solves physics.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::octree3::{Octree3,Octant3};
use fs_cutfem::elastic3::{ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveControl,SolveBudget,MultiLoadOcOptions,MultiLoadOcTermination};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria,inherit_raw_densities3};
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
    if args.len()>3{return Err("usage: adaptive_elastic_sdf3 [ROUNDS [UPDATES [TOTAL_KRYLOV_ITERATIONS]]]".into());}
    let rounds=args.first().map_or(Ok(2),|s|s.parse::<usize>())?;
    let updates=args.get(1).map_or(Ok(3),|s|s.parse::<usize>())?;
    let budget=args.get(2).map_or(Ok(250000),|s|s.parse::<usize>())?;
    if !(1..=4).contains(&rounds){return Err("rounds must lie in 1..=4".into());}
    let mut tree=Octree3::uniform(1,5,2048)?;
    let domain=HexCell::try_new([0.0;3],[1.0;3])?;
    let material=IsotropicElastic::new(1.0,0.3,1.0)?;
    let mut geometry_poll=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut geometry_poll)?;
    let mut solve_poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut solve_poll);
    let mut previous:Option<(Vec<Octant3>,Vec<f64>)>=None;
    println!("round,background_cells,active_cells,master_nodes,iteration,compliance,volume_fraction");
    for round in 0..rounds {
        let op=AdaptiveElasticity3::build(domain,&tree,&Slab,&material,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut geometry)?;
        let y=op.body_load(&|_|[0.0,-1.0,0.0],||ControlFlow::Continue(()))?;
        let z=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(()))?;
        let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());
        let raw=match &previous {
            Some((leaves,rho))=>inherit_raw_densities3(leaves,rho,study.operator().leaves(),&mut control)?,
            None=>vec![0.5;study.cells()],
        };
        let rho=study.feasible_start(&raw,0.5,1e-8,&mut control)?;
        let loads=[LoadCase{force:&y,weight:0.3},LoadCase{force:&z,weight:0.7}];
        let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,
            MultiLoadOcOptions{max_iterations:updates,..Default::default()},&mut control);
        for row in &report.history {println!("{round},{},{},{},{},{:.17e},{:.17e}",
            tree.leaves().len(),study.cells(),study.operator().nodes().len(),row.iteration,row.compliance,row.volume_fraction);}
        eprintln!("round={round} termination={:?} cumulative_krylov={} DWR=false continuum_certified=false",report.termination,report.work.linear_iterations);
        if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure) {
            return Err(format!("stopped with accepted prefix: {:?}",report.evaluation_stop).into());
        }
        if round+1==rounds {break;}
        let volumes=study.operator().volumes();let mut indicators=vec![0.0;study.cells()];
        for (load,u) in loads.iter().zip(&report.displacements) {
            let energies=study.operator().scale_quadratic_forms(u)?;
            for i in 0..indicators.len(){indicators[i]+=load.weight*study.operator().scales()[i]*energies[i]/volumes[i];}
        }
        if !indicators.iter().all(|v|v.is_finite()){return Err("nonfinite refinement indicator".into());}
        let mut candidates:Vec<_>=study.operator().leaves().iter().copied().zip(indicators).filter(|(c,_)|c.level()<5).collect();
        candidates.sort_by(|a,b|b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        let marks:Vec<_>=candidates.iter().take(2).map(|(c,_)|*c).collect();
        if marks.is_empty(){return Err("refinement level budget exhausted".into());}
        previous=Some((study.operator().leaves().to_vec(),report.rho));
        tree=tree.refined(&marks,||ControlFlow::Continue(()))?;
    }
    Ok(())
}
