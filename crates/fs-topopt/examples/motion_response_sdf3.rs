//! Fit volume/surface displacement responses under prescribed embedded motion.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example motion_response_sdf3 -- 40 250000
//! Optional four target values follow the two budgets, in case/observation order.
//! Without targets a disclosed synthetic forward design supplies them, using
//! the SAME work budget. Synthetic fitting is not experimental calibration or
//! unique density recovery. The geometry, loads and length units are dimensionless.
//! Append --estimate to compare the accepted objective with an enriched grid
//! under the remaining solve budget and mark up to two ORIGINAL-grid cells.
//! This reports a two-grid discrepancy, not a reoptimized or certified design.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveControl,SolveBudget};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{ResponseCase3,ResponseTarget3,ResponseDesignStudy3,ResponseDesignOptions3};
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3,ReferenceResponseTarget3,ResponseRefinementOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval{
        if a==HeightAxis::X{Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(0.17,0.17)-Interval::new(0.83,0.83)}else{Interval::new(0.0,0.0)}
    }
}
fn motion(p:[f64;3],n:[f64;3])->[f64;3]{if n[0]>0.0{[0.02,0.0,0.0]}else{[0.0,0.0,0.007*p[1]]}}
fn other(p:[f64;3],n:[f64;3])->[f64;3]{if n[0]>0.0{[0.0,0.01,0.005*p[2]]}else{[-0.01,0.0,0.0]}}
fn main()->Result<(),Box<dyn std::error::Error>>{
    let mut args:Vec<_>=std::env::args().skip(1).collect();
    let estimate=args.last().is_some_and(|s|s=="--estimate");
    if estimate {let _=args.pop();}
    if args.len()>2&&args.len()!=6{return Err("usage: motion_response_sdf3 [STEPS [TOTAL_KRYLOV [T00 T01 T10 T11]]]".into());}
    let steps=args.first().map_or(Ok(40),|s|s.parse::<usize>())?;
    let budget=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let supplied=if args.len()==6{
        let v=args[2..].iter().map(|s|s.parse::<f64>()).collect::<Result<Vec<_>,_>>()?;
        if v.iter().any(|v|!v.is_finite()){return Err("targets must be finite".into());}Some(v)
    }else{None};
    let mut gp=|_|ControlFlow::Continue(());let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
        &Octree3::uniform(1,4,4096)?,&Slab,&IsotropicElastic::new(1.0,0.3,1.0)?,&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut geometry)?;
    let force=op.body_load(&|_|[0.002,0.0,-0.003],||ControlFlow::Continue(()))?;
    let other_force:Vec<_>=force.iter().map(|v|-0.5*v).collect();
    let q=op.body_load(&|p|[0.0,0.0,1.0+p[0]],||ControlFlow::Continue(()))?;
    let r=op.body_load(&|p|[1.0+p[1],0.0,0.0],||ControlFlow::Continue(()))?;
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
    let mut poll=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut poll);
    let template=[ResponseTarget3{q:&q,target:0.0,scale:0.02,weight:0.7},ResponseTarget3{q:&r,target:0.0,scale:0.02,weight:0.3}];
    let values=match supplied{
        Some(v)=>{eprintln!("targets_source=caller; units=dimensionless_reference_integrals");v}
        None=>{
            eprintln!("targets_source=synthetic_forward_design; experimental_calibration=false; unique_recovery_claimed=false");
            let reference=(0..study.cells()).map(|i|0.30+0.05*(i%5)as f64).collect::<Vec<_>>();
            study.evaluate_responses(&reference,&[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&template},
                ResponseCase3{force:&other_force,prescribed:Some(&other),targets:&template}],Default::default(),&mut control)?
                .responses.into_iter().flatten().collect()
        }
    };
    let t0=[ResponseTarget3{target:values[0],..template[0]},ResponseTarget3{target:values[1],..template[1]}];
    let t1=[ResponseTarget3{target:values[2],..template[0]},ResponseTarget3{target:values[3],..template[1]}];
    let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&t0},ResponseCase3{force:&other_force,prescribed:Some(&other),targets:&t1}];
    let initial=vec![0.5;study.cells()];let options=ResponseDesignOptions3::default();
    let mut design=ResponseDesignStudy3::new(&mut study,&cases,&initial,options,&mut control)?;
    let outcome=design.run(steps);
    println!("iteration,objective,volume_fraction,constraint_violation");
    for row in design.history(){println!("{},{:.17e},{:.17e},{:.17e}",row.iteration,row.objective,row.volume_fraction,row.constraint_violation);}
    println!("case,observation,response,target");
    for (i,case) in design.accepted().responses.iter().enumerate(){for (j,value) in case.iter().enumerate(){println!("{i},{j},{value:.17e},{:.17e}",values[2*i+j]);}}
    println!("cell,raw_density,projected_density");
    for (i,(r,p)) in design.point().iter().zip(&design.accepted().projected_rho).enumerate(){println!("{i},{r:.17e},{p:.17e}");}
    eprintln!("evaluations={}; linear_iterations={}; setup_applications={}; galerkin_products={}; feasible={}; continuum_certified=false",
        design.evaluations(),design.work().linear_iterations,design.work().preconditioner_operator_applications,
        design.work().preconditioner_galerkin_products,design.constraint_violation()<=options.tolerance);
    let report=outcome?;eprintln!("stop={:?}; numerical_kkt={:?}; converged={}",report.stop,report.solution.kkt,report.solution.converged);
    if design.constraint_violation()>options.tolerance{return Err("accepted SQP point is still infeasible; no feasible design claimed".into());}
    if estimate {
        let accepted=design.accepted().clone();
        drop(design);
        let fine=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
            &Octree3::uniform(2,4,4096)?,&Slab,&IsotropicElastic::new(1.0,0.3,1.0)?,&|_|false,&|_,_|true,
            ElasticityOptions3::default(),Default::default(),Default::default(),&mut geometry)?;
        let f=|_:[f64;3]|[0.002,0.0,-0.003];
        let f_other=|p|f(p).map(|v|-0.5*v);
        let observe_q=|p:[f64;3]|[0.0,0.0,1.0+p[0]];
        let observe_r=|p:[f64;3]|[1.0+p[1],0.0,0.0];
        let targets0=[
            ReferenceResponseTarget3{observation:ReferenceLoad3::body(&observe_q),target:values[0],scale:0.02,weight:0.7},
            ReferenceResponseTarget3{observation:ReferenceLoad3::body(&observe_r),target:values[1],scale:0.02,weight:0.3},
        ];
        let targets1=[ReferenceResponseTarget3{target:values[2],..targets0[0]},ReferenceResponseTarget3{target:values[3],..targets0[1]}];
        let experiments=[
            ReferenceResponseCase3{load:ReferenceLoad3::body(&f),prescribed:Some(&motion),targets:&targets0},
            ReferenceResponseCase3{load:ReferenceLoad3::body(&f_other),prescribed:Some(&other),targets:&targets1},
        ];
        let evidence=study.estimate_response_enrichment(AdaptiveSolveSpace3::jacobi(fine,100_000_000),
            &accepted,&experiments,ResponseRefinementOptions3{response:options.response,..Default::default()},&mut control)?;
        let marked=evidence.mark(0.5,2,||ControlFlow::Continue(()))?;
        eprintln!("coarse_objective={:.17e}; enriched_objective={:.17e}; two_grid_change={:.17e}; identity_defect={:.3e}; marking_fraction={:.6}; target_met={}; marked={:?}; linear_iterations={}; continuum_certified=false",
            evidence.coarse_objective,evidence.fine_objective,evidence.correction(),evidence.identity_relative_defect,
            marked.achieved_fraction,marked.target_met,marked.marked,evidence.work.linear_iterations);
    }
    Ok(())
}
