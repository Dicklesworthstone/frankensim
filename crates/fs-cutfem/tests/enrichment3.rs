use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn pair() -> (AdaptiveElasticity3, AdaptiveElasticity3) {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mark = *tree.leaves().iter().find(|c| c.index() == [0, 0, 1]).unwrap();
    let coarse = tree.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    let fine = coarse.refined(&coarse.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap();
    (build(&coarse), build(&fine))
}
fn polynomial(p: &[f64; 3]) -> [f64; 3] {
    [p[0] * (1.0 + p[1]) * (1.0 + p[2]), p[0] * p[1], p[0] * p[2]]
}
#[test]
fn nested_mixed_level_transfer_preserves_trilinear_fields_and_physical_constraints() {
    let (coarse, fine) = pair();
    let transfer = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap();
    let u: Vec<f64> = coarse.nodes().iter().flat_map(polynomial).collect();
    let fine_u = transfer.prolongate(&u, || ControlFlow::Continue(())).unwrap();
    let physical = fine.physical_displacements(&fine_u).unwrap();
    for (p, u) in fine.physical_nodes().iter().zip(physical.chunks_exact(3)) {
        for (a, b) in polynomial(p).iter().zip(u) { assert!((a - b).abs() < 1e-12); }
    }
    assert_eq!(transfer.parents().len(), fine.cells());
    assert_eq!(transfer.prolongate(&u, || ControlFlow::Continue(())).unwrap(), fine_u);
}
#[test]
fn transfer_inherits_the_physical_stiffness_field_not_a_refiltered_design() {
    let (mut coarse, fine) = pair();
    let scales: Vec<_> = (0..coarse.cells()).map(|i| 0.2 + 0.03 * i as f64).collect();
    coarse.set_scales(&scales).unwrap();
    let map = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap();
    for (&value, &parent) in map.inherited_scales().iter().zip(map.parents()) {
        assert_eq!(value.to_bits(), scales[parent].to_bits());
    }
}
#[test]
fn localized_weak_residual_matches_actual_reduced_operator_including_ghosts() {
    let (mut op, _) = pair();
    op.set_scales(&vec![0.6; op.cells()]).unwrap();
    let body = |p: [f64; 3]| [p[1], -1.0, p[0] - 0.5];
    let u: Vec<_> = (0..op.n()).map(|i| if op.fixed()[i / 3] { 0.0 } else { (i % 13) as f64 / 13.0 }).collect();
    let w: Vec<_> = (0..op.n()).map(|i| if op.fixed()[i / 3] { 0.0 } else { (i % 7) as f64 / 7.0 }).collect();
    let rhs = op.body_load(&body, || ControlFlow::Continue(())).unwrap();
    let mut au = vec![0.0; op.n()]; op.apply(&u, &mut au);
    let expected: f64 = w.iter().zip(&rhs).zip(au).map(|((w, b), a)| w * (b - a)).sum();
    let terms = op.cell_residuals(&u, &w, &body, || ControlFlow::Continue(())).unwrap();
    let actual: f64 = terms.iter().map(|r| r.residual()).sum();
    assert!((expected - actual).abs() < 1e-11 * expected.abs().max(1.0));
    assert!(terms.iter().map(|r| r.ghost.abs()).sum::<f64>() > 1e-5);
    let no_ghost: f64 = terms.iter().map(|r| r.load - r.bulk).sum();
    assert!((no_ghost - expected).abs() > 1e-5);
}
#[test]
fn interrupted_transfer_and_localization_publish_no_partial_result() {
    let (coarse, fine) = pair();
    assert!(matches!(AdaptiveTransfer3::new(&coarse, &fine, 0, || ControlFlow::Continue(())), Err(ElasticityError3::Invalid(_))));
    assert!(matches!(AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Break(())), Err(ElasticityError3::Cancelled)));
    let map = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap();
    let u: Vec<f64> = coarse.nodes().iter().flat_map(polynomial).collect();
    let mut calls = 0;
    assert!(matches!(map.prolongate(&u, || { calls += 1; if calls > 3 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) } }), Err(ElasticityError3::Cancelled)));
    assert!(map.prolongate(&[0.0], || ControlFlow::Continue(())).is_err());
    assert!(coarse.cell_residuals(&u, &u, &|_| [f64::NAN; 3], || ControlFlow::Continue(())).is_err());
    assert!(matches!(coarse.cell_residuals(&u, &u, &|_| [1.0; 3], || ControlFlow::Break(())), Err(ElasticityError3::Cancelled)));
}
#[test]
fn a_stale_zero_field_is_not_an_equilibrium_and_coarsening_is_not_enrichment() {
    let (coarse, fine) = pair();
    let load = coarse.body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(())).unwrap();
    assert!((coarse.field_residual(&vec![0.0; coarse.n()], &load, || ControlFlow::Continue(())).unwrap() - 1.0).abs() < 1e-12);
    assert!(AdaptiveTransfer3::new(&fine, &coarse, 100_000, || ControlFlow::Continue(())).is_err());
}
