//! Native density loads in actual two-grid forces, goal transfer and markers.
use std::{cell::Cell, collections::BTreeMap, ops::ControlFlow};
use fs_blake3::ContentHash;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticMaterialField3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_dwr::elasticity3::{GoalError3, GoalFields3};
use fs_dwr::elasticity3::reaction::{AffineGoalForm3, MaterialBodyLoads3, estimate_affine_assembled_goal3};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2]-0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2])-Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    // Synthetic equal-stiffness states with independent source density.
    let field = |cell: HexCell| {
        let density = if cell.lo()[1] < 0.5 { 2.0 } else { 8.0 };
        let c = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1000.0 } else { 0.0 }));
        Ok(ElasticAcousticState::try_new_mandel(c, density,
            ContentHash([density as u8; 32]), 1.0).unwrap())
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(QuadratureOptions3 { depth: 0, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        ElasticMaterialField3::ByCell(&field), &|p| p[0] == 0.0,
        Default::default(), &mut quadrature).unwrap()
}
fn pair() -> (AdaptiveElasticity3, AdaptiveElasticity3) {
    let tree = Octree3::uniform(1, 3, 512).unwrap();
    let marked = *tree.leaves().iter().find(|cell| cell.index() == [1, 1, 1]).unwrap();
    let refined = tree.refined(&[marked], || ControlFlow::Continue(())).unwrap();
    let (mut coarse, mut fine) = (build(&tree), build(&refined));
    let scales: Vec<_> = (0..coarse.cells()).map(|i| 0.45+0.04*i as f64).collect();
    coarse.set_scales(&scales).unwrap();
    let inherited = AdaptiveTransfer3::new(&coarse, &fine, 100_000,
        || ControlFlow::Continue(())).unwrap().inherited_scales();
    fine.set_scales(&inherited).unwrap();
    (coarse, fine)
}
fn acceleration(p: [f64; 3]) -> [f64; 3] { [0.2*p[0].powi(8), 0.0, -1.0+0.3*p[1].powi(8)] }
fn dead(p: [f64; 3]) -> [f64; 3] { [0.1+0.03*p[2], 0.0, -0.1] }
fn solve(op: &AdaptiveElasticity3, rhs: &[f64]) -> Vec<f64> {
    op.solve_controlled(rhs, 1e-11, 30_000, 32, |_| ControlFlow::Continue(()))
        .unwrap().coefficients().to_vec()
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

#[test]
fn native_material_work_localizes_and_reconstructs_total_force_compliance() {
    let (coarse, fine) = pair();
    let transfer = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap();
    let bc = coarse.prepare_material_body_load(&acceleration, 100_000, || ControlFlow::Continue(())).unwrap();
    let bf = fine.prepare_material_body_load(&acceleration, 100_000, || ControlFlow::Continue(())).unwrap();
    let mut rc: Vec<_> = (0..coarse.cells()).map(|i| 0.25+0.09*i as f64).collect();
    rc[0] = 0.0; rc[1] = 1.0;
    let rf: Vec<_> = transfer.parents().iter().map(|&p| rc[p]).collect();
    let body = MaterialBodyLoads3 { coarse: &bc, fine: &bf, coarse_occupancy: &rc, fine_occupancy: &rf };
    let gc = bc.load(&rc, || ControlFlow::Continue(())).unwrap();
    let gf = bf.load(&rf, || ControlFlow::Continue(())).unwrap();
    let load = EquilibriumLoad3 { external: ReferenceLoad3::body(&dead), prescribed: None };
    let add = |op, g: &[f64]| load.assemble(op, || ControlFlow::Continue(())).unwrap()
        .iter().zip(g).map(|(f,g)| f+g).collect::<Vec<_>>();
    let qc = AffineGoalForm3 { gradient: add(&coarse, &gc), offset: 0.0 };
    let qf = AffineGoalForm3 { gradient: add(&fine, &gf), offset: 0.0 };
    let uc = solve(&coarse, &qc.gradient); let uf = solve(&fine, &qf.gradient);
    let fields = GoalFields3::compliance(&uc, &uf);
    let report = estimate_affine_assembled_goal3(&transfer, load, Some(body), &qc, &qf,
        fields, Default::default(), || ControlFlow::Continue(())).unwrap();
    let delta = report.fine_value-report.coarse_value;
    assert!((report.correction()-delta).abs() < 1e-8*delta.abs().max(1.0));
    assert!(report.identity_relative_defect < 1e-8);
    assert!(report.linearized.field_residuals.iter().all(|r| *r <= 1e-8));

    let v = transfer.prolongate(&uc, || ControlFlow::Continue(())).unwrap();
    let w: Vec<_> = uf.iter().zip(&v).map(|(u,v)| u-v).collect();
    let rows = fine.cell_reference_residuals(&v, &w, load.external, || ControlFlow::Continue(())).unwrap();
    let work = bf.occupancy_pullback(&w, || ControlFlow::Continue(())).unwrap();
    let mut expected = BTreeMap::<_, (f64, f64)>::new();
    let mut body_signal = 0.0;
    for (i, row) in rows.iter().enumerate() {
        let extra = rf[i]*work[i]; body_signal += extra.abs();
        let term = (row.load+extra)-row.bulk-row.ghost;
        let entry = expected.entry(coarse.leaves()[transfer.parents()[i]]).or_default();
        entry.0 += term; entry.1 += term.abs();
    }
    assert!(body_signal > 1e-10);
    for (leaf, (signed, mass)) in expected {
        let actual = report.linearized.cells[&leaf];
        assert!((actual.signed_dwr-signed).abs() < 1e-12);
        assert!((actual.marking_mass-mass).abs() < 1e-12);
    }
    // The varying acceleration intentionally gives different native cut-rule
    // goal integrals; omitting this load-dependent goal transfer loses identity.
    let body_transfer = dot(&gf, &v)-dot(&gc, &uc);
    assert!(body_transfer.abs() > 1e-12);
    assert!((report.linearized.goal_transfer-(dot(&qf.gradient, &v)-dot(&qc.gradient, &uc))).abs() < 1e-12);
    assert!(matches!(estimate_affine_assembled_goal3(&transfer, load, None, &qc, &qf,
        fields, Default::default(), || ControlFlow::Continue(())), Err(GoalError3::FieldResidual { field: "coarse-primal", .. })));

    // Affine offsets remain separate even with the material force enabled.
    let ac = AffineGoalForm3 { offset: 0.03, ..qc };
    let af = AffineGoalForm3 { offset: -0.04, ..qf };
    let shifted = estimate_affine_assembled_goal3(&transfer, load, Some(body), &ac, &af,
        fields, Default::default(), || ControlFlow::Continue(())).unwrap();
    assert_eq!(shifted.offset_transfer, -0.04-0.03);
    assert!((shifted.correction()-(shifted.fine_value-shifted.coarse_value)).abs() < 1e-8);
}

#[test]
fn noninherited_occupancy_refuses_before_reference_callbacks_or_field_use() {
    let (coarse, fine) = pair();
    let transfer = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap();
    let bc = coarse.prepare_material_body_load(&acceleration, 100_000, || ControlFlow::Continue(())).unwrap();
    let bf = fine.prepare_material_body_load(&acceleration, 100_000, || ControlFlow::Continue(())).unwrap();
    let rc = vec![0.7; coarse.cells()]; let mut rf = vec![0.7; fine.cells()]; rf[0] = 0.6;
    let calls = Cell::new(0usize);
    let law = |p| { calls.set(calls.get()+1); dead(p) };
    let load = EquilibriumLoad3 { external: ReferenceLoad3::body(&law), prescribed: None };
    let qc = AffineGoalForm3 { gradient: vec![0.0; 3*coarse.nodes().len()], offset: 0.0 };
    let qf = AffineGoalForm3 { gradient: vec![0.0; 3*fine.nodes().len()], offset: 0.0 };
    let body = MaterialBodyLoads3 { coarse: &bc, fine: &bf, coarse_occupancy: &rc, fine_occupancy: &rf };
    let result = estimate_affine_assembled_goal3(&transfer, load, Some(body), &qc, &qf,
        GoalFields3::compliance(&[], &[]), Default::default(), || ControlFlow::Continue(()));
    assert!(matches!(result, Err(GoalError3::Invalid("enriched physical body-load occupancy is not parent-inherited"))));
    assert_eq!(calls.get(), 0);
}
