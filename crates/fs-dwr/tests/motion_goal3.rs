//! G0/G3/G4: imposed motion belongs in primal residuals, not adjoint loads.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_dwr::elasticity3::{estimate_motion_goal3, estimate_reference_goal3, GoalError3, GoalFields3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64;3], hi: [f64;3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X {
            Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)
        } else { Interval::new(0.0,0.0) }
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0;3],[1.0;3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(), &|_| false, &|_,n| n[0] < 0.0,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap()
}
fn force(_: [f64;3], n: [f64;3]) -> [f64;3] {
    if n[0] > 0.0 { [0.0,0.0,-0.02] } else { [0.0;3] }
}
fn observation(p: [f64;3], n: [f64;3]) -> [f64;3] {
    if n[0] > 0.0 { [0.0, 0.25*p[2], 1.0+p[1]] } else { [0.0;3] }
}
fn motion(p: [f64;3], _: [f64;3]) -> [f64;3] { [0.004*p[1], 0.007*p[2], 0.003] }
fn solve(op: &AdaptiveElasticity3, b: &[f64]) -> Vec<f64> {
    op.solve_controlled(b, 1e-10, 20_000, 32, |_| ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }

#[test]
fn g0_motion_localization_matches_the_actual_operator_on_hanging_nodes() {
    let t = Octree3::uniform(1,4,4096).unwrap();
    let t = t.refined(&[*t.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap();
    let mut op = build(&t);
    op.set_scales(&(0..op.cells()).map(|i|0.2+0.1*(i%7) as f64).collect::<Vec<_>>()).unwrap();
    let u: Vec<_> = (0..op.n()).map(|i| (i%7) as f64/11.0).collect();
    let w: Vec<_> = (0..op.n()).map(|i| (i%11) as f64/13.0-0.2).collect();
    let law = ReferenceLoad3::traction(&force);
    let rhs = op.reference_load_with_motion(law, Some(&motion), || ControlFlow::Continue(())).unwrap();
    let mut au = vec![0.0;op.n()]; op.apply(&u,&mut au);
    let expected = dot(&w,&rhs)-dot(&w,&au);
    let cells = op.cell_residuals_with_motion(&u,&w,law,Some(&motion),||ControlFlow::Continue(())).unwrap();
    let actual: f64 = cells.iter().map(|c|c.residual()).sum();
    assert!((actual-expected).abs() < 1e-11*expected.abs().max(1.0));
    let omitted: f64 = op.cell_reference_residuals(&u,&w,law,||ControlFlow::Continue(())).unwrap()
        .iter().map(|c|c.residual()).sum();
    assert!((omitted-actual).abs() > 1e-5);
}

#[test]
fn g3_motion_two_grid_identity_uses_independent_homogeneous_observation_adjoints() {
    let t = Octree3::uniform(1,4,4096).unwrap();
    let c = build(&t);
    let mut f = build(&t.refined(&t.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(())).unwrap());
    let scales = AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap().inherited_scales();
    f.set_scales(&scales).unwrap();
    let transfer = AdaptiveTransfer3::new(&c,&f,2_000_000,||ControlFlow::Continue(())).unwrap();
    let law = ReferenceLoad3::traction(&force); let goal = ReferenceLoad3::traction(&observation);
    let uc = solve(&c,&c.reference_load_with_motion(law,Some(&motion),||ControlFlow::Continue(())).unwrap());
    let uf = solve(&f,&f.reference_load_with_motion(law,Some(&motion),||ControlFlow::Continue(())).unwrap());
    let zc = solve(&c,&c.reference_load(goal,||ControlFlow::Continue(())).unwrap());
    let zf = solve(&f,&f.reference_load(goal,||ControlFlow::Continue(())).unwrap());
    let fields = GoalFields3 { coarse_primal:&uc, fine_primal:&uf, coarse_adjoint:&zc, fine_adjoint:&zf };
    let a = estimate_motion_goal3(&transfer,law,Some(&motion),goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    assert!(a.identity_relative_defect < 1e-8);
    assert!((a.correction()-(a.fine_value-a.coarse_value)).abs() < 1e-8*a.fine_value.abs().max(a.coarse_value.abs()));
    assert!(!a.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
    assert!(matches!(estimate_reference_goal3(&transfer,law,goal,fields,Default::default(),||ControlFlow::Continue(())),
        Err(GoalError3::FieldResidual { field:"coarse-primal", .. })));
    let b = estimate_motion_goal3(&transfer,law,Some(&motion),goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    assert_eq!(a.correction().to_bits(),b.correction().to_bits());
}

#[test]
fn g4_zero_motion_path_is_compatible_and_interrupted_lifting_is_not_returned() {
    let op = build(&Octree3::uniform(1,4,4096).unwrap());
    let law = ReferenceLoad3::traction(&force);
    assert_eq!(op.reference_load(law,||ControlFlow::Continue(())).unwrap(),
        op.reference_load_with_motion(law,None,||ControlFlow::Continue(())).unwrap());
    let u=vec![0.0;op.n()]; let w=vec![1.0;op.n()];
    let a=op.cell_reference_residuals(&u,&w,law,||ControlFlow::Continue(())).unwrap();
    let b=op.cell_residuals_with_motion(&u,&w,law,None,||ControlFlow::Continue(())).unwrap();
    for (a,b) in a.iter().zip(b) { assert_eq!(a.residual().to_bits(),b.residual().to_bits()); }
    let scales=op.scales().to_vec();let calls=std::cell::Cell::new(0);
    let g=|p,n|{calls.set(calls.get()+1);motion(p,n)};
    let result=op.reference_load_with_motion(law,Some(&g),||if calls.get()>0 {ControlFlow::Break(())}else{ControlFlow::Continue(())});
    assert!(matches!(result,Err(fs_cutfem::elastic3::ElasticityError3::Cancelled)));
    assert_eq!(op.scales(),scales);
}
