//! Minimum material under sampled stress limits with prescribed actuator motion.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example motion_stress_sdf3 -- 40 250000 0.1
//! Arguments: accepted updates, cumulative Krylov allowance, stress cap.
//! Dimensionless slab, E=1, nu=0, axial strain +/-0.02 plus separate body loads.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::stress::{StressMeasure3,StressMotions3};
use fs_topopt::sdf3::design::{StressDesignOptions3,StressDesignStudy3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],axis:HeightAxis)->Interval {
        if axis==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}
        else {Interval::new(0.0,0.0)}
    }
}
fn extension(p:[f64;3],_:[f64;3])->[f64;3] {[0.02*(p[0]-0.17),0.0,0.0]}
fn compression(p:[f64;3],n:[f64;3])->[f64;3] {extension(p,n).map(|v|-v)}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>3 {return Err("usage: motion_stress_sdf3 [UPDATES [KRYLOV [STRESS_CAP]]]".into());}
    let steps=args.first().map_or(Ok(40),|s|s.parse::<usize>())?;
    let budget=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let cap=args.get(2).map_or(Ok(0.1),|s|s.parse::<f64>())?;
    if !cap.is_finite()||cap<=0.0 {return Err("stress cap must be finite and positive".into());}
    let mut gp=|_|ControlFlow::Continue(());let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
        &Octree3::uniform(1,4,4096)?,&Slab,&IsotropicElastic::new(1.0,0.0,1.0)?,&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut geometry)?;
    let body=op.body_load(&|_|[0.0,0.0,-0.001],||ControlFlow::Continue(()))?;
    let opposite:Vec<f64>=body.iter().map(|v|-v).collect();
    let loads=[LoadCase{force:&body,weight:1.0},LoadCase{force:&opposite,weight:1.0}];
    let motions:&StressMotions3<'_>=&[Some(&extension),Some(&compression)];
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,
        SimpParams{beta:2.0,penal:3.0,..Default::default()});
    let rho=vec![0.65;study.cells()];
    let options=StressDesignOptions3{stress_limit:cap,density_floor:0.1,..Default::default()};
    let mut p=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut p);
    let mut session=StressDesignStudy3::new_with_prescribed_motion(&mut study,&loads,motions,&rho,options,
        StressMeasure3::SampledPeakBound,&mut control)?;
    let result=session.run(steps);
    println!("iteration,volume_fraction,stress_aggregate,sampled_relaxed_max,sampled_physical_max,violation,feasible");
    for r in session.history() {println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",r.iteration,
        r.volume_fraction,r.stress_aggregate,r.sampled_relaxed_max,r.sampled_physical_max,r.constraint_violation,r.feasible);}
    let selected=session.best_feasible().unwrap_or(session.accepted());
    println!("case,external_force_work,sampled_relaxed_max,sampled_physical_max");
    for i in 0..loads.len() {println!("{i},{:.17e},{:.17e},{:.17e}",selected.case_compliances[i],selected.case_relaxed_max[i],selected.case_physical_max[i]);}
    println!("cell,raw_density,projected_density");
    for (i,(r,p)) in selected.rho.iter().zip(&selected.projected_rho).enumerate() {println!("{i},{r:.17e},{p:.17e}");}
    eprintln!("selected_latest={}; selected_volume={:.17e}; stress_cap={cap:.17e}; work={:?}; continuum_certified=false; actuator_force_constrained=false; actuator_energy_claimed=false",
        selected.rho==session.accepted().rho,selected.volume_fraction,session.work());
    let report=result?;eprintln!("latest_optimizer_stop={:?}; latest_optimizer_kkt={:?}",report.stop,report.kkt);
    if session.best_feasible().is_none() {return Err("no accepted design satisfies the numerical stress cap".into());}
    Ok(())
}
