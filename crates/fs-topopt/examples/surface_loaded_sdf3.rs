//! Optimize a pressure/shear-loaded implicit 3-D cantilever without remeshing.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example surface_loaded_sdf3 -- 3 250000
//! Arguments: maximum accepted updates and cumulative Krylov iterations.
//! Dimensionless unit box; solid x < 0.7+0.1*y*(1-y); clamp x=0. Two independent
//! dead loads: inward pressure 1+z and uniform -z surface traction, weights 1/2.
//! No body force, nodal point-load substitution, follower-pressure or DWR claim.
//! Prints accepted history to stdout; complete load totals and stop reason to
//! stderr. A resource/numerical stop prints the accepted prefix and exits nonzero.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl,MultiLoadOcOptions,MultiLoadOcTermination};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria};
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
fn main()->Result<(),Box<dyn std::error::Error>>{
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>2{return Err("usage: surface_loaded_sdf3 [UPDATES [TOTAL_KRYLOV_ITERATIONS]]".into());}
    let updates=args.first().map_or(Ok(3),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250000),|s|s.parse::<usize>())?;
    let mut gp=|_|ControlFlow::Continue(());
    let mut geometry=QuadratureControl3::new(QuadratureOptions3::default(),&mut gp)?;
    let op=AdaptiveElasticity3::build_with_surface(HexCell::try_new([0.0;3],[1.0;3])?,
        &Octree3::uniform(1,4,4096)?,&CurvedEnd,&IsotropicElastic::new(1.0,0.3,1.0)?,
        &|p|p[0]==0.0,ElasticityOptions3::default(),Default::default(),&mut geometry)?;
    let pressure=op.pressure_load(&|p|1.0+p[2],||ControlFlow::Continue(()))?;
    let shear=op.surface_load(&|_,_|[0.0,0.0,-1.0],||ControlFlow::Continue(()))?;
    eprintln!("reference_area={:.17e}; pressure_resultant={:?}; pressure_moment={:?}; shear_resultant={:?}; surface_is_reference=true",
        pressure.area,pressure.resultant,pressure.moment,shear.resultant);
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
    let rho=vec![0.5;study.cells()];let mut sp=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget{total_iterations:iterations,..Default::default()},&mut sp);
    let report=controlled_sdf3_optimality_criteria(&mut study,
        &[LoadCase{force:&pressure.rhs,weight:0.5},LoadCase{force:&shear.rhs,weight:0.5}],&rho,
        MultiLoadOcOptions{max_iterations:updates,..Default::default()},&mut control);
    println!("iteration,compliance,pressure_compliance,shear_compliance,volume_fraction");
    for r in &report.history{println!("{},{:.17e},{:.17e},{:.17e},{:.17e}",r.iteration,r.compliance,r.case_compliances[0],r.case_compliances[1],r.volume_fraction);}
    eprintln!("termination={:?}; iterations={}; geometry_rebuilt=false; continuum_certified=false",report.termination,report.work.linear_iterations);
    if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure){
        return Err(format!("stopped with accepted prefix: {:?}",report.evaluation_stop).into());
    }
    Ok(())
}
