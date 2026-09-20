//! G0/G1/G3/G4: mixed physical load laws, virtual work and stale geometry refusal.
use std::cell::Cell;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::surface::{ReferenceLoad3, SurfaceForce3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2]-0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73) }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z {1.0} else {0.0}; Interval::new(d,d)
    }
}
fn build(mixed: bool, surface: bool) -> AdaptiveElasticity3 {
    let t = Octree3::uniform(1,4,4096).unwrap();
    let t = if mixed {t.refined(&[*t.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap()} else {t};
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 {depth:1,..Default::default()},&mut poll).unwrap();
    let domain = HexCell::try_new([0.0;3],[1.0;3]).unwrap();
    let material = IsotropicElastic::new(1.0,0.3,1.0).unwrap();
    if surface { AdaptiveElasticity3::build_with_surface(domain,&t,&Slab,&material,&|p|p[0]==0.0,
        ElasticityOptions3::default(),Default::default(),&mut q).unwrap() }
    else { AdaptiveElasticity3::build(domain,&t,&Slab,&material,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap() }
}
fn body(_: [f64;3]) -> [f64;3] { [0.0,0.0,0.2] }
fn pressure(p: [f64;3]) -> f64 { 1.0+0.2*p[0] }
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }
#[test]
fn g1_mixed_load_has_analytic_trace_work_on_uniform_and_hanging_grids() {
    for mixed in [false,true] {
        let op=build(mixed,true);
        let law=ReferenceLoad3 {body:Some(&body),surface:Some(SurfaceForce3::Pressure(&pressure))};
        let rhs=op.reference_load(law,||ControlFlow::Continue(())).unwrap();
        let w:Vec<f64>=op.nodes().iter().flat_map(|p|[0.0,0.0,p[0]]).collect();
        let exact=0.2*0.73/2.0-0.5-0.2/3.0;
        assert!((dot(&w,&rhs)-exact).abs()<1e-9);
        let zero=vec![0.0;op.n()];
        let cells=op.cell_reference_residuals(&zero,&w,law,||ControlFlow::Continue(())).unwrap();
        assert!((cells.iter().map(|r|r.load).sum::<f64>()-exact).abs()<1e-9);
        assert!(cells.iter().all(|r|r.bulk==0.0&&r.ghost==0.0));
    }
}
#[test]
fn g0_local_residual_equals_actual_nodal_rhs_minus_stiffness_work() {
    let mut op=build(true,true);
    let scales:Vec<_>=(0..op.cells()).map(|i|0.25+0.1*(i%6) as f64).collect();op.set_scales(&scales).unwrap();
    let law=ReferenceLoad3 {body:Some(&body),surface:Some(SurfaceForce3::Pressure(&pressure))};
    let u:Vec<_>=(0..op.n()).map(|i|if op.fixed()[i/3]{0.0}else{(i%7) as f64/7.0}).collect();
    let w:Vec<_>=(0..op.n()).map(|i|if op.fixed()[i/3]{0.0}else{(i%11) as f64/11.0-0.5}).collect();
    let rhs=op.reference_load(law,||ControlFlow::Continue(())).unwrap();let mut au=vec![0.0;op.n()];op.apply(&u,&mut au);
    let local=op.cell_reference_residuals(&u,&w,law,||ControlFlow::Continue(())).unwrap();
    let expected=dot(&w,&rhs)-dot(&w,&au);
    assert!((local.iter().map(|r|r.residual()).sum::<f64>()-expected).abs()<1e-11*expected.abs().max(1.0));
}
#[test]
fn g3_body_only_and_pressure_only_reuse_existing_nodal_assembly() {
    let op=build(true,true);
    assert_eq!(op.reference_load(ReferenceLoad3::body(&body),||ControlFlow::Continue(())).unwrap(),op.body_load(&body,||ControlFlow::Continue(())).unwrap());
    assert_eq!(op.reference_load(ReferenceLoad3::pressure(&pressure),||ControlFlow::Continue(())).unwrap(),op.pressure_load(&pressure,||ControlFlow::Continue(())).unwrap().rhs);
    let f=|p:[f64;3],n:[f64;3]|n.map(|v|-pressure(p)*v);
    assert_eq!(op.reference_load(ReferenceLoad3::traction(&f),||ControlFlow::Continue(())).unwrap(),op.pressure_load(&pressure,||ControlFlow::Continue(())).unwrap().rhs);
    assert_eq!(op.reference_load(ReferenceLoad3::default(),||ControlFlow::Continue(())).unwrap(),vec![0.0;op.n()]);
}
#[test]
fn g4_missing_surface_nonfinite_force_and_mid_callback_cancel_do_not_return_partial_loads() {
    let legacy=build(false,false);let calls=Cell::new(0);
    let counted=|p|{calls.set(calls.get()+1);body(p)};
    let law=ReferenceLoad3 {body:Some(&counted),surface:Some(SurfaceForce3::Pressure(&pressure))};
    assert!(matches!(legacy.reference_load(law,||ControlFlow::Continue(())),Err(ElasticityError3::Invalid(_))));assert_eq!(calls.get(),0);
    let zero=vec![0.0;legacy.n()];
    assert!(legacy.cell_reference_residuals(&zero,&zero,law,||ControlFlow::Continue(())).is_err());assert_eq!(calls.get(),0);
    let op=build(true,true);let before=op.scales().to_vec();
    let p=|x|{calls.set(calls.get()+1);pressure(x)};
    let result=op.reference_load(ReferenceLoad3::pressure(&p),||if calls.get()>2{ControlFlow::Break(())}else{ControlFlow::Continue(())});
    assert!(matches!(result,Err(ElasticityError3::Cancelled)));assert_eq!(op.scales(),before);
    let bad=|_|f64::NAN;
    assert!(op.reference_load(ReferenceLoad3::pressure(&bad),||ControlFlow::Continue(())).is_err());
    let zero=vec![0.0;op.n()];
    assert!(op.cell_reference_residuals(&zero,&zero,ReferenceLoad3::pressure(&bad),||ControlFlow::Continue(())).is_err());
}
