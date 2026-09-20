//! G0/G1/G3/G4: real surface loading, constraints, solves and density derivatives.
use std::cell::Cell;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityError3,ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3,QuadratureError3,QuadratureWork3};
use fs_cutfem::quad3::surface::surface_cell_rules3;
use fs_cutfem::octree3::Octree3;
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{p[0]-self.0}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{Interval::new(lo[0],hi[0])-Interval::new(self.0,self.0)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval{
        let d=if a==HeightAxis::X{1.0}else{0.0};Interval::new(d,d)
    }
}
fn domain()->HexCell{HexCell::try_new([0.0;3],[1.0;3]).unwrap()}
fn material(nu:f64)->IsotropicElastic{IsotropicElastic::new(1.0,nu,1.0).unwrap()}
fn quad()->QuadratureOptions3{QuadratureOptions3{depth:1,..Default::default()}}
fn cart(surface:bool,nu:f64)->CutElasticity3{
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(quad(),&mut p).unwrap();
    if surface {CutElasticity3::build_with_surface(domain(),[2;3],&Slab(0.73),&material(nu),&|p|p[0]==0.0,
        Default::default(),Default::default(),&mut q).unwrap()}
    else {CutElasticity3::build(domain(),[2;3],&Slab(0.73),&material(nu),&|p|p[0]==0.0,Default::default(),&mut q).unwrap()}
}
fn adaptive()->AdaptiveElasticity3{
    let t=Octree3::uniform(1,4,4096).unwrap();
    let mark=*t.leaves().iter().find(|c|c.index()==[1,0,0]).unwrap();
    let t=t.refined(&[mark],||ControlFlow::Continue(())).unwrap();
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(quad(),&mut p).unwrap();
    AdaptiveElasticity3::build_with_surface(domain(),&t,&Slab(0.73),&material(0.0),&|p|p[0]==0.0,
        Default::default(),Default::default(),&mut q).unwrap()
}
fn near(a:f64,b:f64,t:f64){assert!((a-b).abs()<=t,"{a} vs {b}");}
#[test]
fn g1_pressure_sign_and_analytic_extension_on_cartesian_and_hanging_cuts(){
    let op=cart(true,0.0);let load=op.pressure_load(&|_|1.0,||ControlFlow::Continue(())).unwrap();
    near(load.area,1.0,1e-12);near(load.resultant[0],-1.0,1e-12);
    near(load.moment[1],-0.5,1e-12);near(load.moment[2],0.5,1e-12);
    let u=op.solve_controlled(&load.rhs,1e-10,10000,32,|_|ControlFlow::Continue(())).unwrap();
    near(u.compliance(),0.73,1e-8);
    for (x,u) in op.nodes().iter().zip(u.coefficients().chunks_exact(3)){near(u[0],-x[0],1e-8);near(u[1],0.0,1e-8);near(u[2],0.0,1e-8);}
    let op=adaptive();assert!(op.physical_nodes().len()>op.nodes().len());
    let load=op.pressure_load(&|_|1.0,||ControlFlow::Continue(())).unwrap();
    let u=op.solve_controlled(&load.rhs,1e-10,10000,32,|_|ControlFlow::Continue(())).unwrap();
    let physical=op.physical_displacements(u.coefficients()).unwrap();
    near(u.compliance(),0.73,1e-8);
    for (x,u) in op.physical_nodes().iter().zip(physical.chunks_exact(3)){near(u[0],-x[0],1e-8);near(u[1],0.0,1e-8);near(u[2],0.0,1e-8);}
}
#[test]
fn g0_nonuniform_traction_virtual_work_matches_direct_physical_trace_integration(){
    let op=adaptive();let law=|p:[f64;3],_:[f64;3]|[p[1],p[2]*p[2],-1.0-p[1]*p[2]];
    let load=op.surface_load(&law,||ControlFlow::Continue(())).unwrap();
    let u:Vec<f64>=(0..op.n()).map(|i|((i*13)%23) as f64/23.0).collect();
    let physical=op.physical_displacements(&u).unwrap();let mut integrated=0.0;
    let mut cp=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(quad(),&mut cp).unwrap();
    for nodes in op.cell_nodes(){
        let lo=op.physical_nodes()[nodes[0]];let hi=op.physical_nodes()[nodes[7]];
        let rules=surface_cell_rules3(&Slab(0.73),HexCell::try_new(lo,hi).unwrap(),Default::default(),&mut q).unwrap();
        for node in rules.points(){
            let t:[f64;3]=std::array::from_fn(|i|(node.position[i]-lo[i])/(hi[i]-lo[i]));
            let f=law(node.position,node.normal);
            for a in 0..8{
                let basis=(0..3).map(|i|if a&(1<<i)==0{1.0-t[i]}else{t[i]}).product::<f64>();
                for c in 0..3{integrated+=node.weight*basis*physical[3*nodes[a]+c]*f[c];}
            }
        }
    }
    near(load.rhs.iter().zip(&u).map(|(f,u)|f*u).sum(),integrated,1e-12);
    near(load.resultant[0],0.5,1e-12);near(load.resultant[1],1.0/3.0,1e-12);near(load.resultant[2],-1.25,1e-12);
    for (i,f) in load.rhs.iter().enumerate(){if op.fixed()[i/3]{assert_eq!(*f,0.0);}}
}
#[test]
fn g3_dead_surface_load_reuses_geometry_and_retains_existing_density_pullback(){
    let mut op=cart(true,0.3);let load=op.surface_load(&|_,_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let scales=vec![0.6;op.cells()];op.set_scales(&scales).unwrap();
    let solve=|op:&CutElasticity3|op.solve_controlled(&load.rhs,1e-11,10000,32,|_|ControlFlow::Continue(())).unwrap();
    let base=solve(&op);let gradient=op.scale_quadratic_forms(base.coefficients()).unwrap();
    for i in 0..op.cells(){
        let mut s=scales.clone();s[i]+=1e-4;op.set_scales(&s).unwrap();let plus=solve(&op).compliance();
        s[i]-=2e-4;op.set_scales(&s).unwrap();let minus=solve(&op).compliance();let fd=(plus-minus)/2e-4;
        assert!((fd+gradient[i]).abs()<2e-4*fd.abs());
    }
    let same=op.surface_load(&|_,_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();assert_eq!(same.rhs,load.rhs);
}
#[test]
fn g3_opt_in_surface_does_not_change_stiffness_body_loads_or_geometry(){
    let a=cart(false,0.3);let b=cart(true,0.3);assert_eq!(a.nodes(),b.nodes());assert_eq!(a.volumes(),b.volumes());
    let x:Vec<f64>=(0..a.n()).map(|i|(i%7) as f64).collect();let mut ax=vec![0.0;a.n()];let mut bx=ax.clone();a.apply(&x,&mut ax);b.apply(&x,&mut bx);assert_eq!(ax,bx);
    assert_eq!(a.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap(),b.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap());
    let called=Cell::new(false);
    assert!(matches!(a.surface_load(&|_,_|{called.set(true);[0.0;3]},||ControlFlow::Continue(())),Err(ElasticityError3::Invalid(_))));assert!(!called.get());
}
#[test]
fn g4_bad_laws_and_mid_load_cancellation_never_publish_partial_vectors(){
    let op=adaptive();let before=op.scales().to_vec();let calls=Cell::new(0);
    assert!(matches!(op.surface_load(&|_,_|{calls.set(calls.get()+1);[1.0;3]},||if calls.get()>=3{ControlFlow::Break(())}else{ControlFlow::Continue(())}),Err(ElasticityError3::Cancelled)));
    assert_eq!(calls.get(),3);assert_eq!(op.scales(),before);
    assert!(matches!(op.pressure_load(&|_|f64::NAN,||ControlFlow::Continue(())),Err(ElasticityError3::Invalid(_))));
    let a=op.pressure_load(&|_|1.0,||ControlFlow::Continue(())).unwrap();let b=op.pressure_load(&|_|1.0,||ControlFlow::Continue(())).unwrap();assert_eq!(a.rhs,b.rhs);
    let zero=op.surface_load(&|_,_|[0.0;3],||ControlFlow::Continue(())).unwrap();assert!(zero.rhs.iter().all(|f|*f==0.0));near(zero.area,1.0,1e-12);
}
#[test]
fn g4_surface_preparation_can_stop_after_bulk_assembly_without_returning_operator(){
    let mut cp=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(quad(),&mut cp).unwrap();
    CutElasticity3::build(domain(),[2;3],&Slab(0.73),&material(0.0),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap();let bulk=q.work();
    let mut cp=|w:QuadratureWork3|if w.boxes>bulk.boxes{ControlFlow::Break(())}else{ControlFlow::Continue(())};let mut q=QuadratureControl3::new(quad(),&mut cp).unwrap();
    assert!(matches!(CutElasticity3::build_with_surface(domain(),[2;3],&Slab(0.73),&material(0.0),&|p|p[0]==0.0,
        Default::default(),Default::default(),&mut q),Err(ElasticityError3::Quadrature(QuadratureError3::Cancelled))));
    assert!(q.work().points>=bulk.points);
}
