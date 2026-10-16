//! G0/G1/G3/G4: real raw-surface supports, nonzero motion and density lifting.
use std::cell::Cell;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityError3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::dirichlet::EmbeddedDirichletOptions3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3, QuadratureError3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { (p[0]-self.0)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        let x=Interval::new(lo[0],hi[0]);
        (x-Interval::new(self.0,self.0))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo:[f64;3], hi:[f64;3], a:HeightAxis)->Interval {
        if a==HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])
            -Interval::new(self.0,self.0)-Interval::new(0.83,0.83) }
        else {Interval::new(0.0,0.0)}
    }
}
fn domain()->HexCell {HexCell::try_new([0.0;3],[1.0;3]).unwrap()}
fn tree(mixed:bool)->Octree3 {
    let t=Octree3::uniform(1,4,4096).unwrap();
    if mixed {t.refined(&[*t.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap()}else{t}
}
fn adaptive(a:f64,nu:f64,mixed:bool)->AdaptiveElasticity3 {
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(domain(),&tree(mixed),&Slab(a),
        &IsotropicElastic::new(1.0,nu,1.0).unwrap(),&|_|false,&|_,n|n[0]<0.0,
        ElasticityOptions3::default(),EmbeddedDirichletOptions3::default(),Default::default(),&mut q).unwrap()
}
fn traction(_: [f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[1.0,0.0,0.0]} else {[0.0;3]}}
fn translation(_: [f64;3],_:[f64;3])->[f64;3] {[0.125,0.0,0.0]}
fn add(a:&[f64],b:&[f64])->Vec<f64>{a.iter().zip(b).map(|(a,b)|a+b).collect()}

#[test]
fn g1_cartesian_embedded_support_and_nonzero_motion_reproduce_axial_solution() {
    for a in [0.17,0.4999] {
        let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
        let op=CutElasticity3::build_with_embedded_dirichlet(domain(),[2;3],&Slab(a),
            &IsotropicElastic::new(1.0,0.0,1.0).unwrap(),&|_|false,&|_,n|n[0]<0.0,
            Default::default(),Default::default(),Default::default(),&mut q).unwrap();
        assert!(!op.fixed().iter().any(|v|*v));
        assert!((op.embedded_dirichlet_area().unwrap()-1.0).abs()<1e-8);
        let f=op.surface_load(&traction,||ControlFlow::Continue(())).unwrap().rhs;
        let g=op.prescribed_displacement_load(&translation,||ControlFlow::Continue(())).unwrap();
        let solved=op.solve_controlled(&add(&f,&g),1e-10,10_000,16,|_|ControlFlow::Continue(())).unwrap();
        for (p,u) in op.nodes().iter().zip(solved.coefficients().chunks_exact(3)) {
            assert!((u[0]-(p[0]-a+0.125)).abs()<2e-7,"a={a} p={p:?} u={u:?}");
            assert!(u[1].abs()<2e-7&&u[2].abs()<2e-7);
        }
        assert!(solved.residual_claim().euclidean().unwrap()<1e-10);
    }
}
#[test]
fn g1_mixed_octree_uses_the_same_boundary_stiffness_and_lifting_transpose() {
    let op=adaptive(0.17,0.0,true);
    assert!(op.physical_nodes().len()>op.nodes().len());
    assert!(!op.fixed().iter().any(|v|*v));
    let f=op.surface_load(&traction,||ControlFlow::Continue(())).unwrap().rhs;
    let g=op.prescribed_displacement_load(&translation,||ControlFlow::Continue(())).unwrap();
    let solved=op.solve_controlled(&add(&f,&g),1e-10,10_000,16,|_|ControlFlow::Continue(())).unwrap();
    let physical=op.physical_displacements(solved.coefficients()).unwrap();
    for (p,u) in op.physical_nodes().iter().zip(physical.chunks_exact(3)) {
        assert!((u[0]-(p[0]-0.17+0.125)).abs()<2e-7);
        assert!(u[1].abs()<2e-7&&u[2].abs()<2e-7);
    }
    let zeros=op.prescribed_displacement_load(&|_,_|[0.0;3],||ControlFlow::Continue(())).unwrap();
    assert!(zeros.iter().all(|v|*v==0.0));
}
#[test]
fn g3_nonzero_motion_load_derivatives_are_required_and_match_full_resolves() {
    let mut op=adaptive(0.17,0.3,true);
    let prescribed=|p:[f64;3],_:[f64;3]|[0.03+0.02*p[1],-0.01*p[2],0.04];
    let s:Vec<f64>=(0..op.cells()).map(|i|0.35+0.07*(i%7) as f64).collect();op.set_scales(&s).unwrap();
    let f=op.surface_load(&traction,||ControlFlow::Continue(())).unwrap().rhs;
    let work=|op:&AdaptiveElasticity3| {
        let g=op.prescribed_displacement_load(&prescribed,||ControlFlow::Continue(())).unwrap();
        op.solve_controlled(&add(&f,&g),1e-11,30_000,16,|_|ControlFlow::Continue(())).unwrap()
    };
    let base=work(&op);let stiffness=op.scale_quadratic_forms(base.coefficients()).unwrap();
    let lifting=op.prescribed_displacement_scale_work(&prescribed,base.coefficients(),||ControlFlow::Continue(())).unwrap();
    let mut omission_detected=false;
    for i in 0..s.len() {
        let h=1e-4;let mut p=s.clone();p[i]+=h;op.set_scales(&p).unwrap();let plus=work(&op).compliance();
        p[i]-=2.0*h;op.set_scales(&p).unwrap();let minus=work(&op).compliance();
        let fd=(plus-minus)/(2.0*h);let expected=2.0*lifting[i]-stiffness[i];
        assert!((fd-expected).abs()<5e-4*fd.abs().max(1e-8),"cell={i}: {fd} vs {expected}");
        omission_detected|=(fd+stiffness[i]).abs()>1e-2*fd.abs().max(1e-8);
    }
    assert!(omission_detected);op.set_scales(&s).unwrap();assert_eq!(work(&op).coefficients(),base.coefficients());
}
#[test]
fn g0_original_support_admission_remains_and_empty_or_invalid_patches_refuse() {
    let material=IsotropicElastic::new(1.0,0.3,1.0).unwrap();
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    assert!(matches!(CutElasticity3::build(domain(),[2;3],&Slab(0.17),&material,&|_|false,Default::default(),&mut q),Err(ElasticityError3::Invalid(_))));
    let calls=Cell::new(0);
    let bad=CutElasticity3::build_with_embedded_dirichlet(domain(),[2;3],&Slab(0.17),&material,&|_|false,
        &|_,_|{calls.set(calls.get()+1);true},Default::default(),EmbeddedDirichletOptions3{beta:f64::NAN},Default::default(),&mut q);
    assert!(matches!(bad,Err(ElasticityError3::Invalid(_))));assert_eq!(calls.get(),0);
    assert!(matches!(CutElasticity3::build_with_embedded_dirichlet(domain(),[2;3],&Slab(0.17),&material,&|_|false,&|_,_|false,
        Default::default(),Default::default(),Default::default(),&mut q),Err(ElasticityError3::Invalid(_))));
}
#[test]
fn g4_cancelled_boundary_build_and_invalid_lifting_never_publish_partial_results() {
    let stop=Cell::new(false);let mut poll=|_|if stop.get(){ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let result=CutElasticity3::build_with_embedded_dirichlet(domain(),[2;3],&Slab(0.17),&IsotropicElastic::new(1.0,0.3,1.0).unwrap(),
        &|_|false,&|_,_|{stop.set(true);true},Default::default(),Default::default(),Default::default(),&mut q);
    assert!(matches!(result,Err(ElasticityError3::Quadrature(QuadratureError3::Cancelled))));
    let op=adaptive(0.17,0.3,true);let scales=op.scales().to_vec();
    assert!(matches!(op.prescribed_displacement_load(&translation,||ControlFlow::Break(())),Err(ElasticityError3::Cancelled)));
    assert!(op.prescribed_displacement_load(&|_,_|[f64::NAN;3],||ControlFlow::Continue(())).is_err());
    assert!(op.prescribed_displacement_scale_work(&translation,&[1.0],||ControlFlow::Continue(())).is_err());
    assert_eq!(op.scales(),scales);
}
#[test]
fn g0_transfer_rejects_changed_boundary_method_or_penalty() {
    let coarse=adaptive(0.17,0.3,false);let fine=adaptive(0.17,0.3,true);
    assert!(AdaptiveTransfer3::new(&coarse,&fine,100_000,||ControlFlow::Continue(())).is_ok());
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let other=AdaptiveElasticity3::build_with_embedded_dirichlet(domain(),&tree(true),&Slab(0.17),&IsotropicElastic::new(1.0,0.3,1.0).unwrap(),
        &|_|false,&|_,n|n[0]<0.0,Default::default(),EmbeddedDirichletOptions3{beta:64.0},Default::default(),&mut q).unwrap();
    assert!(AdaptiveTransfer3::new(&coarse,&other,100_000,||ControlFlow::Continue(())).is_err());
}
