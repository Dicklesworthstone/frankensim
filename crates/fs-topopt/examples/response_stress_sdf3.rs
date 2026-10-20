//! Minimum material with simultaneous actuator windows and sampled stress caps.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example response_stress_sdf3 -- 80 500000 0.0012 0.008 0.1
//! Arguments: accepted updates, total Krylov allowance, min/max axial reaction,
//! stress cap. All units are dimensionless reference units (E=1, nu=0).
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::limits::*;
use fs_topopt::sdf3::response::limits::strength::ResponseStressOptions3;
use fs_topopt::sdf3::response::limits::strength::design::{ResponseStressDesignOptions3,ResponseStressDesignStudy3};
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
fn motion(p:[f64;3],_:[f64;3])->[f64;3] {[0.02*(p[0]-0.17),0.0,0.0]}
fn opposite(p:[f64;3],n:[f64;3])->[f64;3] {motion(p,n).map(|v|-v)}
fn right(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[1.0,0.0,0.0]}else{[0.0;3]}}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>5 {return Err("usage: response_stress_sdf3 [UPDATES [KRYLOV [MIN_REACTION [MAX_REACTION [STRESS_CAP]]]]]".into());}
    let steps=args.first().map_or(Ok(80),|s|s.parse::<usize>())?;
    let budget=args.get(1).map_or(Ok(500000),|s|s.parse::<usize>())?;
    let lo=args.get(2).map_or(Ok(0.0012),|s|s.parse::<f64>())?;
    let hi=args.get(3).map_or(Ok(0.008),|s|s.parse::<f64>())?;
    let cap=args.get(4).map_or(Ok(0.1),|s|s.parse::<f64>())?;
    if ![lo,hi,cap].iter().all(|v|v.is_finite())||lo<=0.0||hi<=lo||cap<=0.0 {
        return Err("require finite 0 < MIN_REACTION < MAX_REACTION and positive STRESS_CAP".into());
    }
    let mut gp=|_|ControlFlow::Continue(());let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
        &Octree3::uniform(1,4,4096)?,&Slab,&IsotropicElastic::new(1.0,0.0,1.0)?,&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut geometry)?;
    let force=op.body_load(&|_|[0.0,0.0,-0.001],||ControlFlow::Continue(()))?;
    let neg:Vec<_>=force.iter().map(|v|-v).collect();
    let q=op.body_load(&|_|[1.0,0.0,0.0],||ControlFlow::Continue(()))?;
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,
        SimpParams{beta:2.0,penal:3.0,..Default::default()});
    let tension=[ResponseLimit3{response:BoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtLeast,bound:lo,scale:0.002},
        ResponseLimit3{response:BoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtMost,bound:hi,scale:0.002},
        ResponseLimit3{response:BoundedResponse3::Displacement(&q),sense:ResponseLimitSense3::AtMost,bound:0.006,scale:0.005}];
    let compression=[ResponseLimit3{sense:ResponseLimitSense3::AtMost,bound:-lo,..tension[0]},
        ResponseLimit3{sense:ResponseLimitSense3::AtLeast,bound:-hi,..tension[1]},
        ResponseLimit3{sense:ResponseLimitSense3::AtLeast,bound:-0.006,..tension[2]}];
    let cases=[ResponseLimitCase3{force:&force,prescribed:Some(&motion),limits:&tension},
        ResponseLimitCase3{force:&neg,prescribed:Some(&opposite),limits:&compression}];
    let options=ResponseStressDesignOptions3{limits:ResponseStressOptions3{stress_limit:cap,..Default::default()},..Default::default()};
    let rho=vec![0.65;study.cells()];let mut p=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut p);
    let mut session=ResponseStressDesignStudy3::new(&mut study,&cases,&rho,options,&mut control)?;
    let initial=session.accepted().response.volume_fraction;let outcome=session.run(steps);
    println!("iteration,volume_fraction,joint_constraint,worst_response_residual,stress_aggregate,physical_max,feasible");
    for r in session.history() {println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",r.iteration,
        r.volume_fraction,r.constraint,r.worst_response_residual,r.stress_aggregate,r.sampled_physical_max,r.feasible);}
    let e=session.best_feasible().unwrap_or(session.accepted());
    println!("case,limit,sense,value,bound,scale,residual");
    for (i,c) in cases.iter().enumerate() {for (j,l) in c.limits.iter().enumerate() {
        println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e}",l.sense,e.response.measurements[i][j],l.bound,l.scale,e.response.residuals[i][j]);
    }}
    println!("case,relaxed_max,physical_max,stress_cap,external_force_work");
    for i in 0..cases.len() {println!("{i},{:.17e},{:.17e},{cap:.17e},{:.17e}",
        e.stress.case_relaxed_max[i],e.stress.case_physical_max[i],e.stress.case_compliances[i]);}
    println!("cell,raw_density,projected_density");
    for (i,(r,v)) in e.response.rho.iter().zip(&e.response.projected_rho).enumerate() {println!("{i},{r:.17e},{v:.17e}");}
    eprintln!("initial_volume={initial:.17e}; retained_volume={:.17e}; joint_constraint={:.17e}; worst_original_residual={:.17e}; selected_latest={}; feasible={}; work={:?}; continuum_certified=false; actuator_energy_claimed=false",
        e.response.volume_fraction,e.constraint,e.worst_residual,e.response.rho==session.point(),e.feasible(options.optimizer.tolerance),session.work());
    let report=outcome?;eprintln!("latest_al_stop={:?}; latest_al_kkt={:?}; selected_stationarity_not_claimed=true",report.stop,report.kkt);
    let _retained=session.into_best_feasible()?;
    Ok(())
}
