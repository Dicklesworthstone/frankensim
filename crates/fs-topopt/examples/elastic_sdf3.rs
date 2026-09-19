//! Raw 3-D implicit-domain density optimization, with no geometric rebuilding.
//! cargo run -p fs-topopt --features cutfem-marquee --release --example elastic_sdf3 -- 5 250000
//! Arguments: accepted-update cap, cumulative Krylov-iteration cap. Prints a
//! CSV history followed by accepted cell densities; interrupted work exits
//! nonzero after printing its accepted prefix. All units are dimensionless.
//! The curved design domain is z < 0.7 + 0.1*x*(1-x) inside the unit box.
//! This field is implicit, not an exact distance function. The two dead body
//! loads are independent. No continuum, manufacturing or optimality claim.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{MultiLoadOcOptions,MultiLoadOcTermination,SimpParams,SolveBudget,SolveControl};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria};

struct CurvedCantilever;
impl CutSdf3 for CurvedCantilever {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.7-0.1*p[0]*(1.0-p[0])}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);let one=Interval::new(1.0,1.0);
        Interval::new(lo[2],hi[2])-Interval::new(0.7,0.7)-Interval::new(0.1,0.1)*x*(one-x)
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        match a {
            HeightAxis::X=>Interval::new(-0.1,-0.1)*(Interval::new(1.0,1.0)-Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])),
            HeightAxis::Y=>Interval::new(0.0,0.0),HeightAxis::Z=>Interval::new(1.0,1.0),
        }
    }
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<String>=std::env::args().skip(1).collect();
    if args.len()>2 {return Err("usage: elastic_sdf3 [UPDATES [TOTAL_KRYLOV_ITERATIONS]]".into());}
    let updates=args.first().map_or(Ok(5),|s|s.parse::<usize>())?;
    let iterations=args.get(1).map_or(Ok(250_000),|s|s.parse::<usize>())?;
    let mut poll=|_|ControlFlow::Continue(());
    let mut quadrature=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll)?;
    let op=CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3])?,[4,2,2],&CurvedCantilever,
        &IsotropicElastic::new(1.0,0.3,1.0)?,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut quadrature)?;
    let y=op.body_load(&|_|[0.0,-1.0,0.0],||ControlFlow::Continue(()))?;
    let z=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(()))?;
    let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());
    let rho=vec![0.5;study.cells()];
    let geometry_work=quadrature.work();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget {total_iterations:iterations,..Default::default()},&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,
        &[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}],&rho,
        MultiLoadOcOptions {max_iterations:updates,..Default::default()},&mut control);
    println!("iteration,compliance,volume_fraction,max_change");
    for row in &report.history {println!("{},{:.17e},{:.17e},{:.17e}",row.iteration,row.compliance,row.volume_fraction,row.max_change);}
    eprintln!("termination={:?}; cut_cells={}; geometry_points={}; linear_iterations={}; geometry_rebuilt=false; optimality_certified=false",
        report.termination,study.cells(),geometry_work.points,report.work.linear_iterations);
    if !report.history.is_empty() {
        println!("cell_x,cell_y,cell_z,raw_density,projected_density");
        for ((key,raw),physical) in study.operator().cell_keys().iter().zip(&report.rho).zip(&report.projected_rho) {
            println!("{},{},{},{raw:.17e},{physical:.17e}",key[0],key[1],key[2]);
        }
    }
    if matches!(report.termination,MultiLoadOcTermination::Cancelled|MultiLoadOcTermination::LinearBudget|MultiLoadOcTermination::NumericalFailure) {
        return Err(format!("study interrupted: {:?}",report.evaluation_stop).into());
    }
    Ok(())
}
