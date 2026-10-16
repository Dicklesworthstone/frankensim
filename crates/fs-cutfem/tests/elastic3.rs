use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityError3,ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;

struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 { p[2]-self.0 }
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        Interval::new(lo[2],hi[2])-Interval::new(self.0,self.0)
    }
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],axis:HeightAxis)->Interval {
        let d=if axis==HeightAxis::Z {1.0} else {0.0}; Interval::new(d,d)
    }
}
fn build(cut:f64,poisson:f64)->CutElasticity3 {
    let mut callback=|_|ControlFlow::Continue(());
    let mut control=QuadratureControl3::new(QuadratureOptions3 {depth:1,..Default::default()},&mut callback).unwrap();
    CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),[2;3],&Slab(cut),
        &IsotropicElastic::new(1.0,poisson,1.0).unwrap(),&|p|p[0]==0.0,
        ElasticityOptions3::default(),&mut control).unwrap()
}
fn load(op:&CutElasticity3)->Vec<f64> {
    op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap()
}
fn solve(op:&CutElasticity3,rhs:&[f64])->fs_cutfem::elastic3::ElasticitySolution3 {
    op.solve_controlled(rhs,1e-10,10_000,32,|_|ControlFlow::Continue(())).unwrap()
}
#[test]
fn full_domain_affine_extension_matches_the_analytic_3d_solution() {
    let op=build(1.1,0.0);
    let mut rhs=vec![0.0;op.n()];
    for (id,p) in op.nodes().iter().enumerate() {
        if p[0]==1.0 {
            let a=if p[1]==0.0||p[1]==1.0 {0.5} else {1.0};
            let b=if p[2]==0.0||p[2]==1.0 {0.5} else {1.0};
            rhs[3*id]=a*b/4.0;
        }
    }
    let solution=solve(&op,&rhs);
    for (p,u) in op.nodes().iter().zip(solution.coefficients().chunks_exact(3)) {
        assert!((u[0]-p[0]).abs()<1e-8);
        assert!(u[1].abs()<1e-8&&u[2].abs()<1e-8);
    }
    assert!((solution.compliance()-1.0).abs()<1e-8);
    assert!(solution.residual_claim().euclidean().unwrap()<1e-10);
}
#[test]
fn cut_domain_support_mask_symmetry_and_independent_opposite_loads() {
    let op=build(0.73,0.3);
    let volumes=op.volumes();
    assert!((volumes.iter().sum::<f64>()-0.73).abs()<1e-8);
    assert!(op.volume_bounds().lo()<=0.73&&op.volume_bounds().hi()>=0.73);
    let rhs=load(&op);let a=solve(&op,&rhs);
    let neg:Vec<f64>=rhs.iter().map(|f|-f).collect();let b=solve(&op,&neg);
    assert!(a.compliance()>0.0);
    assert!((a.compliance()-b.compliance()).abs()<1e-10*a.compliance());
    for (i,(&u,&v)) in a.coefficients().iter().zip(b.coefficients()).enumerate() {
        assert!((u+v).abs()<1e-9);
        if op.fixed()[i/3] {assert_eq!(u.to_bits(),0.0_f64.to_bits());}
    }
    let x:Vec<f64>=(0..op.n()).map(|i|(i%7)as f64/7.0).collect();
    let y:Vec<f64>=(0..op.n()).map(|i|(i%11)as f64/11.0).collect();
    let mut ax=vec![0.0;op.n()];let mut ay=ax.clone();op.apply(&x,&mut ax);op.apply_transpose(&y,&mut ay);
    let lhs:f64=ax.iter().zip(&y).map(|(a,b)|a*b).sum();
    let rhs:f64=ay.iter().zip(&x).map(|(a,b)|a*b).sum();
    assert!((lhs-rhs).abs()<1e-12*lhs.abs().max(1.0));
}
#[test]
fn scale_pullback_includes_ghost_terms_and_matches_resolve_differences() {
    let mut op=build(0.73,0.3);let base=vec![0.6;op.cells()];op.set_scales(&base).unwrap();
    let rhs=load(&op);let solution=solve(&op,&rhs);
    let energy=op.scale_quadratic_forms(solution.coefficients()).unwrap();
    let total:f64=energy.iter().zip(&base).map(|(e,s)|e*s).sum();
    assert!((total-solution.compliance()).abs()<1e-8*solution.compliance());
    for cell in 0..op.cells() {
        let h=1e-4;let mut scales=base.clone();scales[cell]+=h;op.set_scales(&scales).unwrap();
        let plus=solve(&op,&rhs).compliance();scales[cell]-=2.0*h;op.set_scales(&scales).unwrap();
        let minus=solve(&op,&rhs).compliance();let fd=(plus-minus)/(2.0*h);
        assert!((fd+energy[cell]).abs()/fd.abs().max(1e-12)<2e-5,"cell {cell}: {fd} vs {}",-energy[cell]);
    }
    op.set_scales(&base).unwrap();
    let before=op.scales().to_vec();let mut bad=base;bad[0]=f64::NAN;
    assert!(op.set_scales(&bad).is_err());assert_eq!(op.scales(),before);
}
#[test]
fn thin_cut_cells_solve_with_explicit_residual_verification() {
    for cut in [0.5001,0.51,0.73] {
        let op=build(cut,0.3);let rhs=load(&op);let solution=solve(&op,&rhs);
        assert!(solution.compliance().is_finite()&&solution.compliance()>0.0);
        assert!(solution.residual_claim().euclidean().unwrap()<1e-10);
    }
}
#[test]
fn cancellation_and_iteration_exhaustion_return_no_usable_field() {
    let op=build(0.73,0.3);let rhs=load(&op);
    assert!(matches!(op.solve_controlled(&rhs,1e-10,10_000,1,|i|if i>=2 {ControlFlow::Break(())}else{ControlFlow::Continue(())}),Err(ElasticityError3::Cancelled)));
    assert!(matches!(op.solve_controlled(&rhs,1e-10,1,1,|_|ControlFlow::Continue(())),Err(ElasticityError3::NotConverged{iterations:1,..})));
    assert!(matches!(op.body_load(&|_|[0.0;3],||ControlFlow::Break(())),Err(ElasticityError3::Cancelled)));
    let a=solve(&op,&rhs);let b=solve(&op,&rhs);assert_eq!(a.coefficients(),b.coefficients());
    assert_eq!(a.compliance().to_bits(),b.compliance().to_bits());
    assert_eq!(a.iterations(),b.iterations());
}
#[test]
fn malformed_geometry_material_and_shapes_refuse_before_solver_work() {
    let mut callback=|_|ControlFlow::Continue(());
    let mut control=QuadratureControl3::new(QuadratureOptions3::default(),&mut callback).unwrap();
    let material=IsotropicElastic::new(1.0,0.3,1.0).unwrap();
    assert!(CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),[0,2,2],&Slab(0.73),&material,&|_|true,ElasticityOptions3::default(),&mut control).is_err());
    assert_eq!(control.work().boxes,0);
    let op=build(0.73,0.3);
    assert!(op.scale_quadratic_forms(&[1.0]).is_err());
    assert!(op.solve_controlled(&[1.0],1e-10,100,1,|_|ControlFlow::Continue(())).is_err());
    assert!(op.body_load(&|_|[f64::NAN;3],||ControlFlow::Continue(())).is_err());
}
