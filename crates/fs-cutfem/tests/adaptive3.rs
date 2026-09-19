use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityOptions3,ElasticityError3};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-self.0}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(self.0,self.0)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z {1.0}else{0.0};Interval::new(d,d)
    }
}
fn tree(refine:bool)->Octree3 {
    let t=Octree3::uniform(1,4,1000).unwrap();
    if !refine{return t;}
    let mark=*t.leaves().iter().find(|c|c.index()==[0,0,1]).unwrap();
    t.refined(&[mark],||ControlFlow::Continue(())).unwrap()
}
fn build(t:&Octree3,cut:f64,nu:f64)->AdaptiveElasticity3 {
    let mut poll=|_|ControlFlow::Continue(());
    let mut c=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),t,&Slab(cut),
        &IsotropicElastic::new(1.0,nu,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut c).unwrap()
}
fn solve(op:&AdaptiveElasticity3,rhs:&[f64])->fs_cutfem::elastic3::ElasticitySolution3 {
    op.solve_controlled(rhs,1e-10,20000,32,|_|ControlFlow::Continue(())).unwrap()
}
#[test]
fn affine_extension_reconstructs_hanging_displacements() {
    let op=build(&tree(true),1.1,0.0);assert!(op.physical_nodes().len()>op.nodes().len());
    let mut rhs=vec![0.0;op.n()];
    for (i,p) in op.nodes().iter().enumerate(){if p[0]==1.0 {
        let a=if p[1]==0.0||p[1]==1.0 {0.5}else{1.0};
        let b=if p[2]==0.0||p[2]==1.0 {0.5}else{1.0};rhs[3*i]=a*b/4.0;
    }}
    let s=solve(&op,&rhs);let full=op.physical_displacements(s.coefficients()).unwrap();
    for (p,u) in op.physical_nodes().iter().zip(full.chunks_exact(3)) {
        assert!((p[0]-u[0]).abs()<1e-8);assert!(u[1].abs()<1e-8&&u[2].abs()<1e-8);
    }
    assert!((s.compliance()-1.0).abs()<1e-8);
}
#[test]
fn uniform_tree_agrees_with_existing_cartesian_operator() {
    let a=build(&tree(false),0.73,0.3);
    let mut poll=|_|ControlFlow::Continue(());
    let mut c=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    let b=CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),[2;3],&Slab(0.73),
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut c).unwrap();
    assert_eq!(a.nodes(),b.nodes());
    let x:Vec<f64>=(0..a.n()).map(|i|(i%11)as f64/11.0).collect();
    let mut ax=vec![0.0;a.n()];let mut bx=ax.clone();a.apply(&x,&mut ax);b.apply(&x,&mut bx);
    for (a,b) in ax.iter().zip(&bx){assert!((a-b).abs()<1e-12);}
}
#[test]
fn adaptive_ghost_density_derivatives_and_energy_identity() {
    let mut op=build(&tree(true),0.73,0.3);
    let base=vec![0.6;op.cells()];op.set_scales(&base).unwrap();
    let rhs=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let s=solve(&op,&rhs);let energies=op.scale_quadratic_forms(s.coefficients()).unwrap();
    let value:f64=energies.iter().zip(&base).map(|(e,s)|e*s).sum();
    assert!((value-s.compliance()).abs()<1e-8*s.compliance());
    for i in 0..op.cells(){
        let mut v=base.clone();v[i]+=1e-4;op.set_scales(&v).unwrap();let plus=solve(&op,&rhs).compliance();
        v[i]-=2e-4;op.set_scales(&v).unwrap();let minus=solve(&op,&rhs).compliance();
        let fd=(plus-minus)/2e-4;assert!((fd+energies[i]).abs()/fd.abs()<5e-5,"cell {i}");
    }
}
#[test]
fn mixed_level_operator_remains_symmetric_and_replays() {
    let op=build(&tree(true),0.73,0.3);
    let x:Vec<f64>=(0..op.n()).map(|i|(i%11)as f64/11.0).collect();
    let y:Vec<f64>=(0..op.n()).map(|i|(i%7)as f64/7.0).collect();
    let mut ax=vec![0.0;op.n()];let mut ay=ax.clone();op.apply(&x,&mut ax);op.apply_transpose(&y,&mut ay);
    let lhs:f64=ax.iter().zip(&y).map(|(a,b)|a*b).sum();let rhs:f64=ay.iter().zip(&x).map(|(a,b)|a*b).sum();
    assert!((lhs-rhs).abs()<1e-12);
    let rhs=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let a=solve(&op,&rhs);let b=solve(&op,&rhs);assert_eq!(a.coefficients(),b.coefficients());
    assert_eq!(a.compliance().to_bits(),b.compliance().to_bits());
}
#[test]
fn cancellation_and_exhaustion_never_return_a_partial_field() {
    let op=build(&tree(true),0.73,0.3);
    let rhs=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    assert!(matches!(op.solve_controlled(&rhs,1e-10,1000,1,|i|if i>=2 {ControlFlow::Break(())}else{ControlFlow::Continue(())}),Err(ElasticityError3::Cancelled)));
    assert!(matches!(op.solve_controlled(&rhs,1e-10,1,1,|_|ControlFlow::Continue(())),Err(ElasticityError3::NotConverged{iterations:1,..})));
    assert!(op.physical_displacements(&[0.0]).is_err());
    assert!(op.body_load(&|_|[f64::NAN;3],||ControlFlow::Continue(())).is_err());
}
