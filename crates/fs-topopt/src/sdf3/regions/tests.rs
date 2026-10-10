use super::*;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{MultiLoadOcOptions, MultiLoadOcTermination, SimpParams, SolveBudget, SolveProgress};
use crate::pipeline::LoadCase;
use crate::sdf3::controlled_sdf3_optimality_criteria;

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
fn fixture(radius: f64, beta: f64) -> (CutDensityStudy3, Vec<f64>, Vec<f64>) {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let op = CutElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), [2; 3], &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap();
    let y = op.body_load(&|_| [0.0, -1.0, 0.0], || ControlFlow::Continue(())).unwrap();
    let z = op.body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(op, radius, SimpParams { penal: 3.0, beta, e_min: 0.05, ..Default::default() }), y, z)
}
fn mask(n: usize) -> Vec<PhysicalRegion3> {
    let mut regions = vec![PhysicalRegion3::Design; n];
    regions[0] = PhysicalRegion3::Solid;
    regions[n - 1] = PhysicalRegion3::Void;
    regions
}

#[test]
fn physical_zero_and_one_survive_projection_and_have_full_chain_gradients() {
    for beta in [1.0, 8.0] {
        let (study, y, z) = fixture(0.15, beta);
        let n = study.cells();
        let mut study = study.with_physical_regions(mask(n)).unwrap();
        let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
        let rho: Vec<_> = (0..n).map(|i| 0.4 + i as f64 * 0.02).collect();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let base = study.evaluate(&rho, &loads, &mut control).unwrap();
        assert_eq!(base.projected_rho[0], 1.0);
        assert_eq!(base.projected_rho[n - 1], 0.0);
        assert_eq!(study.operator().scales()[0], 1.0);
        assert_eq!(study.operator().scales()[n - 1], study.params().e_min);
        assert!(base.volume_gradient[0] > 0.0, "the raw filter control still affects neighboring free cells");
        for i in 0..n {
            let mut plus = rho.clone(); let mut minus = rho.clone();
            let h = 1e-4;
            plus[i] += h; minus[i] -= h;
            let a = study.evaluate(&plus, &loads, &mut control).unwrap();
            let b = study.evaluate(&minus, &loads, &mut control).unwrap();
            assert_eq!((a.projected_rho[0], a.projected_rho[n - 1]), (1.0, 0.0));
            let fd = (a.objective.compliance - b.objective.compliance) / (2.0 * h);
            let vd = (a.volume_fraction - b.volume_fraction) / (2.0 * h);
            assert!((fd - base.objective.gradient[i]).abs() <= 5e-4 * fd.abs().max(1e-8), "compliance beta={beta} cell={i}");
            assert!((vd - base.volume_gradient[i]).abs() <= 5e-5 * vd.abs().max(1e-8), "volume beta={beta} cell={i}");
        }
    }
}

#[test]
fn zero_filter_inactive_controls_do_not_break_actual_oc_descent() {
    let (study, y, z) = fixture(0.0, 1.0);
    let n = study.cells();
    let mut study = study.with_physical_regions(mask(n)).unwrap();
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let rho = study.feasible_start(&vec![0.5; n], 0.55, 1e-8, &mut control).unwrap();
    let report = controlled_sdf3_optimality_criteria(&mut study, &loads, &rho,
        MultiLoadOcOptions { volume_fraction: 0.55, max_iterations: 3, change_tolerance: 0.0, ..Default::default() }, &mut control);
    assert!(report.history.len() > 1, "require a real update, not a no-op: {report:?}");
    assert!(report.history.last().unwrap().compliance < report.history[0].compliance);
    assert_eq!((report.projected_rho[0], report.projected_rho[n - 1]), (1.0, 0.0));
    assert_eq!((report.rho[0], report.rho[n - 1]), (rho[0], rho[n - 1]));
    assert!(report.history.iter().all(|h| h.volume_fraction <= 0.55000001));
    let replay = study.evaluate(&report.rho, &loads, &mut control).unwrap();
    assert_eq!(report.displacements, replay.objective.displacements);
    assert_eq!(report.projected_rho, replay.projected_rho);
}

#[test]
fn prescribed_solids_are_charged_and_cannot_be_removed_to_restore_feasibility() {
    let (study, _, z) = fixture(0.15, 2.0);
    let n = study.cells();
    let mut study = study.with_physical_regions(vec![PhysicalRegion3::Solid; n]).unwrap();
    assert!((study.prescribed_solid_fraction() - 1.0).abs() < 1e-14);
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.feasible_start(&vec![0.5; n], 0.5, 1e-8, &mut c),
        Err(EvaluationStop::Breakdown { stage: "sdf3-volume-restoration" })));
    let result = controlled_sdf3_optimality_criteria(&mut study, &[LoadCase { force: &z, weight: 1.0 }],
        &vec![0.5; n], MultiLoadOcOptions { volume_fraction: 1.0, max_iterations: 2, ..Default::default() }, &mut c);
    assert_eq!(result.termination, MultiLoadOcTermination::NoAcceptableStep);
    assert_eq!(result.history.len(), 1);
    assert!(result.projected_rho.iter().all(|r| *r == 1.0));
}

#[test]
fn interruption_before_evaluation_publication_preserves_scales_and_regions() {
    let (study, _, z) = fixture(0.15, 2.0);
    let n = study.cells(); let regions = mask(n);
    let mut study = study.with_physical_regions(regions.clone()).unwrap();
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&vec![0.5; n], &loads, &mut c).unwrap();
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-evaluation-publish" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.evaluate(&vec![0.6; n], &loads, &mut c), Err(EvaluationStop::Cancelled)));
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.physical_regions(), Some(regions.as_slice()));
}

#[test]
fn region_shape_parent_order_and_cancelled_transfer_are_checked() {
    let (study, _, _) = fixture(0.0, 1.0);
    assert!(study.with_physical_regions(vec![]).is_err());
    let source = [PhysicalRegion3::Solid, PhysicalRegion3::Void, PhysicalRegion3::Design];
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert_eq!(inherit_physical_regions(&source, &[1, 0, 0, 2, 1], &mut c).unwrap(),
        [PhysicalRegion3::Void, PhysicalRegion3::Solid, PhysicalRegion3::Solid, PhysicalRegion3::Design, PhysicalRegion3::Void]);
    assert!(inherit_physical_regions(&source, &[3], &mut c).is_err());
    let mut poll = |_| ControlFlow::Break(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(inherit_physical_regions(&source, &[0, 1], &mut c), Err(EvaluationStop::Cancelled)));
    assert_eq!(source, [PhysicalRegion3::Solid, PhysicalRegion3::Void, PhysicalRegion3::Design]);
}
