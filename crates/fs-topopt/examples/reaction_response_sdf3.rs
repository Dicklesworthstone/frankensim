//! Fit an imposed-motion actuator reaction and a displacement observation.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example reaction_response_sdf3 -- 40 250000 0.0012
//! Arguments: accepted updates, total Krylov budget, target axial reaction.
//! Defaults are declared design targets, NOT measured/calibrated specimen data.
//! All dimensions are reference/dimensionless: E=1, nu=0, strain=0.02, slab
//! thickness=0.66, unit cross-section, projected volume cap=0.5. Reaction sign
//! is the traction exerted on the solid by its right embedded actuator.
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
use fs_topopt::sdf3::response::{ReactionTarget3,ResponseCase3,ResponseTarget3,ProjectedResponseStudy3,ProjectedResponseOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval{
        if a==HeightAxis::X{Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}else{Interval::new(0.0,0.0)}
    }
}
fn motion(p:[f64;3],_:[f64;3])->[f64;3]{[0.02*(p[0]-0.17),0.0,0.0]}
fn right(_:[f64;3],n:[f64;3])->[f64;3]{if n[0]>0.0{[1.0,0.0,0.0]}else{[0.0;3]}}
fn main()->Result<(),Box<dyn std::error::Error>>{
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>3{return Err("usage: reaction_response_sdf3 [UPDATES [KRYLOV [TARGET_REACTION]]]".into());}
    let steps=args.first().map_or(Ok(40),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let target=args.get(2).map_or(Ok(0.0012),|s|s.parse::<f64>())?;
    if !target.is_finite(){return Err("reaction target must be finite".into());}
    let mut gp=|_|ControlFlow::Continue(());let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
        &Octree3::uniform(1,4,4096)?,&Slab,&IsotropicElastic::new(1.0,0.0,1.0)?,&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut geometry)?;
    let force=op.body_load(&|_|[0.0;3],||ControlFlow::Continue(()))?;
    let q=op.body_load(&|_|[1.0,0.0,0.0],||ControlFlow::Continue(()))?;
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
    let displacements=[ResponseTarget3{q:&q,target:0.02*0.66*0.66/2.0,scale:0.005,weight:0.1}];
    let reactions=[ReactionTarget3{mode:&right,target,scale:0.002,weight:1.0}];
    let family:[&[ReactionTarget3<'_>];1]=[&reactions];
    let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&displacements}];
    let options=ProjectedResponseOptions3::default();let initial=vec![0.5;study.cells()];
    let mut sp=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget{total_iterations:iterations,..Default::default()},&mut sp);
    let mut session=ProjectedResponseStudy3::new_with_reactions(&mut study,&cases,&family,&initial,options,&mut control)?;
    let outcome=session.run(steps);
    println!("iteration,objective,volume_fraction,constraint_violation");
    for row in session.history(){println!("{},{:.17e},{:.17e},{:.17e}",row.iteration,row.objective,row.volume_fraction,row.constraint_violation);}
    println!("reaction,target_reaction,displacement_integral");
    println!("{:.17e},{target:.17e},{:.17e}",session.accepted().reaction_responses[0][0],session.accepted().responses[0][0]);
    println!("cell,raw_density,projected_density");
    for (i,(raw,physical)) in session.point().iter().zip(&session.accepted().projected_rho).enumerate(){println!("{i},{raw:.17e},{physical:.17e}");}
    eprintln!("targets_source=declared_design_objective; continuum_certified=false; unique_recovery_claimed=false; actuator_energy_claimed=false; work={:?}",session.work());
    let report=outcome?;eprintln!("stop={:?}; numerical_kkt={:?}; volume_violation={:.9e}",report.stop,report.kkt,session.constraint_violation());
    if session.constraint_violation()>options.optimizer.tolerance{return Err("retained design still violates volume cap".into());}
    Ok(())
}
