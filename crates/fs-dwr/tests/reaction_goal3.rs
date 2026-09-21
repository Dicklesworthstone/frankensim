//! Actual cut-cell physics with independent reaction adjoints and Nitsche offsets.
use std::{cell::Cell, ops::ControlFlow};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_dwr::elasticity3::{GoalError3, GoalFields3, estimate_motion_goal3};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
use fs_dwr::elasticity3::reaction::{AffineMotionGoal3, assemble_affine_motion_goal3, estimate_affine_motion_goal3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        let x=Interval::new(lo[0],hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64;3], hi: [f64;3], a: HeightAxis) -> Interval {
        if a==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}
        else {Interval::new(0.0,0.0)}
    }
}
fn build(t: &Octree3) -> AdaptiveElasticity3 {
    let mut p=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),t,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn motion(p:[f64;3],n:[f64;3])->[f64;3] {
    if n[0]>0.0 {[0.02,0.003*p[1],0.0]} else {[0.0,0.0,0.004*p[2]]}
}
fn mode(p:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[1.0+0.2*p[1],0.0,0.0]}else{[0.0;3]}}
fn body(_:[f64;3])->[f64;3] {[0.003,0.001,-0.002]}
fn observe(p:[f64;3])->[f64;3] {[0.0,0.0,0.2+p[0]]}
fn solve(op:&AdaptiveElasticity3,b:&[f64])->Vec<f64> {
    op.solve_controlled(b,1e-11,30_000,32,|_|ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
fn pair(mixed:bool)->(AdaptiveElasticity3,AdaptiveElasticity3) {
    let t=Octree3::uniform(1,4,4096).unwrap();
    let marks=if mixed {vec![*t.leaves().iter().next().unwrap()]}else{t.leaves().iter().copied().collect()};
    let mut c=build(&t);let mut f=build(&t.refined(&marks,||ControlFlow::Continue(())).unwrap());
    c.set_scales(&(0..c.cells()).map(|i|0.3+0.1*(i%5)as f64).collect::<Vec<_>>()).unwrap();
    let s=AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap().inherited_scales();
    f.set_scales(&s).unwrap();(c,f)
}
#[test]
fn affine_forms_match_actual_reactions_and_complete_two_grid_difference() {
    for mixed in [false,true] {
        let(c,f)=pair(mixed);let transfer=AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap();
        let load=EquilibriumLoad3{external:ReferenceLoad3::body(&body),prescribed:Some(&motion)};
        let goal=AffineMotionGoal3{displacement:ReferenceLoad3::body(&observe),reaction_mode:Some(&mode)};
        let qc=assemble_affine_motion_goal3(&c,goal,Some(&motion),||ControlFlow::Continue(())).unwrap();
        let qf=assemble_affine_motion_goal3(&f,goal,Some(&motion),||ControlFlow::Continue(())).unwrap();
        let uc=solve(&c,&load.assemble(&c,||ControlFlow::Continue(())).unwrap());
        let uf=solve(&f,&load.assemble(&f,||ControlFlow::Continue(())).unwrap());
        let zc=solve(&c,&qc.gradient);let zf=solve(&f,&qf.gradient);
        let fields=GoalFields3{coarse_primal:&uc,fine_primal:&uf,coarse_adjoint:&zc,fine_adjoint:&zf};
        let report=estimate_affine_motion_goal3(&transfer,load,goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
        for (op,u,value) in [(&c,&uc,report.coarse_value),(&f,&uf,report.fine_value)] {
            let reaction=op.embedded_reaction(u,Some(&motion),&mode,||ControlFlow::Continue(())).unwrap().value;
            let q=op.reference_load(goal.displacement,||ControlFlow::Continue(())).unwrap();
            let actual=reaction+q.iter().zip(u).map(|(q,u)|q*u).sum::<f64>();
            assert!((actual-value).abs()<1e-9*actual.abs().max(1e-8));
        }
        let delta=report.fine_value-report.coarse_value;
        assert!((report.correction()-delta).abs()<1e-9);
        if !mixed {
            assert!(report.offset_transfer.abs()>1e-3);
            assert!((report.linearized.correction()-delta).abs()>1e-3,"omitted offset must fail");
        }
        assert!(!report.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
        let again=estimate_affine_motion_goal3(&transfer,load,goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
        assert_eq!(report.correction().to_bits(),again.correction().to_bits());
        let wrong=GoalFields3{coarse_adjoint:&uc,..fields};
        assert!(matches!(estimate_affine_motion_goal3(&transfer,load,goal,wrong,Default::default(),||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{field:"coarse-adjoint",..})));
    }
}
#[test]
fn displacement_only_form_keeps_original_goal_decomposition() {
    let(c,f)=pair(false);let transfer=AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap();
    let load=EquilibriumLoad3{external:ReferenceLoad3::body(&body),prescribed:Some(&motion)};
    let goal=AffineMotionGoal3{displacement:ReferenceLoad3::body(&observe),reaction_mode:None};
    let uc=solve(&c,&load.assemble(&c,||ControlFlow::Continue(())).unwrap());let uf=solve(&f,&load.assemble(&f,||ControlFlow::Continue(())).unwrap());
    let qc=assemble_affine_motion_goal3(&c,goal,load.prescribed,||ControlFlow::Continue(())).unwrap();
    let qf=assemble_affine_motion_goal3(&f,goal,load.prescribed,||ControlFlow::Continue(())).unwrap();
    let zc=solve(&c,&qc.gradient);let zf=solve(&f,&qf.gradient);
    let fields=GoalFields3{coarse_primal:&uc,fine_primal:&uf,coarse_adjoint:&zc,fine_adjoint:&zf};
    let a=estimate_affine_motion_goal3(&transfer,load,goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    let b=estimate_motion_goal3(&transfer,load.external,load.prescribed,goal.displacement,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    assert_eq!(a.offset_transfer,0.0);assert_eq!(a.linearized.correction().to_bits(),b.correction().to_bits());
}
#[test]
fn interrupted_reaction_assembly_and_changed_material_return_no_estimate() {
    let(c,mut f)=pair(false);let before=c.scales().to_vec();let calls=Cell::new(0);
    let h=|p,n|{calls.set(calls.get()+1);mode(p,n)};
    let goal=AffineMotionGoal3{displacement:ReferenceLoad3::default(),reaction_mode:Some(&h)};
    assert!(matches!(assemble_affine_motion_goal3(&c,goal,Some(&motion),||if calls.get()>0{ControlFlow::Break(())}else{ControlFlow::Continue(())}),Err(GoalError3::Physics(fs_cutfem::elastic3::ElasticityError3::Cancelled))));
    assert_eq!(c.scales(),before);
    f.set_scales(&vec![0.9;f.cells()]).unwrap();let transfer=AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap();
    let fields=GoalFields3{coarse_primal:&[],fine_primal:&[],coarse_adjoint:&[],fine_adjoint:&[]};
    assert!(matches!(estimate_affine_motion_goal3(&transfer,EquilibriumLoad3{external:ReferenceLoad3::default(),prescribed:None},goal,fields,Default::default(),||ControlFlow::Continue(())),Err(GoalError3::Invalid(_))));
}
