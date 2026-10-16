//! Real cut-geometry evaluations and accepted optimizer/recovery states.
use super::*;
use crate::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
use crate::sdf3::PhysicalRegion3;
use crate::{SolveBudget, SolveProgress};
use fs_ascent::projected_al::ProjectedAlOptions;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;

const PEAK: StressMeasure3 = StressMeasure3::SampledPeakBound;
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
fn fixture(level: u8, hanging: bool) -> (CutDensityStudy3<AdaptiveSolveSpace3>, Vec<f64>, Vec<f64>) {
    let mut tree = Octree3::uniform(level, 4, 1000).unwrap();
    if hanging {
        let mark = *tree.leaves().iter().find(|c| c.index() == [0, 0, 1]).unwrap();
        tree = tree.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    }
    let mut poll = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let op = AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree, &Slab, &IsotropicElastic::new(3.0, 0.27, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut geometry).unwrap();
    let y = op.body_load(&|p| [0.0, -(0.4 + p[0]), 0.0], || ControlFlow::Continue(())).unwrap();
    let z = op.body_load(&|p| [0.1 * p[1], 0.0, -1.0], || ControlFlow::Continue(())).unwrap();
    let study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15,
        SimpParams { penal: 3.0, beta: 2.0, e_min: 1e-3, ..Default::default() });
    (study, y, z)
}
fn bound(e: &StressEvaluation3, p: f64) {
    let maximum = e.sampled_relaxed_max.max(e.sampled_physical_max);
    assert!(e.aggregate >= maximum, "every returned peak aggregate must bound its numerical samples");
    assert!(e.aggregate <= maximum * (2.0 * e.point_count as f64).powf(1.0 / p) * (1.0 + 1e-12));
    assert!(e.case_relaxed_max.iter().chain(&e.case_physical_max).all(|m| *m <= e.aggregate));
}

#[test]
fn peak_constraint_cannot_hide_a_hot_zero_weight_case_and_legacy_is_unchanged() {
    let (mut study, y, z) = fixture(1, false);
    let hot: Vec<_> = z.iter().map(|f| 20.0 * f).collect();
    let rho = vec![0.6; study.cells()];
    let options = StressOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let loads = [LoadCase { force: &y, weight: 1.0 }, LoadCase { force: &hot, weight: 0.0 }];
    let original = study.operator().elasticity().scales().to_vec();
    let legacy = study.evaluate_stress(&rho, &loads, options, &mut control).unwrap();
    let explicit = study.evaluate_stress_with_measure(&rho, &loads, options,
        StressMeasure3::NormalizedAverage, &mut control).unwrap();
    assert_eq!(legacy.aggregate, explicit.aggregate);
    assert_eq!(legacy.gradient, explicit.gradient);
    assert_eq!(legacy.adjoints, explicit.adjoints);
    assert_eq!(legacy.displacements, explicit.displacements);
    let peak = study.evaluate_stress_with_measure(&rho, &loads, options, PEAK, &mut control).unwrap();
    bound(&peak, options.aggregation_power);
    assert!(peak.case_physical_max[1] > legacy.aggregate);
    assert!(peak.aggregate > 5.0 * legacy.aggregate);
    assert!(peak.adjoints[1].iter().any(|z| z.abs() > 0.0));
    assert_eq!(peak.displacements, legacy.displacements);
    for weight in [1e-250, 0.3, 100.0] {
        let loads = [LoadCase { force: &y, weight: 1.0 }, LoadCase { force: &hot, weight }];
        let reweighted = study.evaluate_stress_with_measure(&rho, &loads, options, PEAK, &mut control).unwrap();
        assert_eq!(peak.aggregate.to_bits(), reweighted.aggregate.to_bits());
        assert_eq!(peak.gradient, reweighted.gradient);
    }
    assert_eq!(study.operator().elasticity().scales(), original);
}

#[test]
fn peak_gradient_includes_physical_relaxation_hanging_constraints_and_fixed_regions() {
    let (study, y, z) = fixture(1, true);
    let n = study.cells();
    let mut labels = vec![PhysicalRegion3::Design; n];
    labels[0] = PhysicalRegion3::Solid;
    labels[n - 1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(labels).unwrap();
    assert!(study.operator().elasticity().physical_nodes().len() > study.operator().elasticity().nodes().len());
    let loads = [LoadCase { force: &y, weight: 0.0 }, LoadCase { force: &z, weight: 1.0 }];
    let options = StressOptions3 { relaxation_power: 1.4, aggregation_power: 16.0, ..Default::default() };
    let rho: Vec<_> = (0..n).map(|i| 0.45 + 0.012 * i as f64).collect();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let base = study.evaluate_stress_with_measure(&rho, &loads, options, PEAK, &mut control).unwrap();
    bound(&base, options.aggregation_power);
    assert_eq!(base.projected_rho[0], 1.0);
    assert_eq!(base.projected_rho[n - 1], 0.0);
    for index in 0..n {
        let h = 2e-5;
        let mut plus = rho.clone(); let mut minus = rho.clone();
        plus[index] += h; minus[index] -= h;
        let a = study.evaluate_stress_with_measure(&plus, &loads, options, PEAK, &mut control).unwrap();
        let b = study.evaluate_stress_with_measure(&minus, &loads, options, PEAK, &mut control).unwrap();
        let fd = (a.aggregate - b.aggregate) / (2.0 * h);
        let scale = fd.abs().max(base.aggregate * 1e-7);
        assert!((fd - base.gradient[index]).abs() / scale < 5e-4,
            "cell={index}, fd={fd}, adjoint={}", base.gradient[index]);
    }
}

#[test]
fn physical_ersatz_stress_is_not_erased_in_prescribed_voids() {
    let (study, _, z) = fixture(1, false);
    let n = study.cells();
    let mut study = study.with_physical_regions(vec![PhysicalRegion3::Void; n]).unwrap();
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = StressOptions3 { aggregation_power: 64.0, ..Default::default() };
    let e = study.evaluate_stress_with_measure(&vec![0.6; n], &loads, options, PEAK, &mut control).unwrap();
    assert_eq!(e.sampled_relaxed_max, 0.0);
    assert!(e.sampled_physical_max > 0.0);
    bound(&e, options.aggregation_power);
    assert!(e.gradient.iter().all(|g| *g == 0.0));
    assert!(e.projected_rho.iter().all(|r| *r == 0.0));
}

#[test]
fn peak_caps_and_cancellation_preserve_incoming_scales_without_partial_results() {
    let (mut study, y, z) = fixture(1, false);
    let rho = vec![0.6; study.cells()];
    let loads = [LoadCase { force: &y, weight: 1.0 }, LoadCase { force: &z, weight: 0.0 }];
    let original = study.operator().elasticity().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_stress_with_measure(&rho, &loads, Default::default(), PEAK, &mut control).unwrap();
    let options = StressOptions3 { max_points: e.point_count - 1, ..Default::default() };
    assert!(matches!(study.evaluate_stress_with_measure(&rho, &loads, options, PEAK, &mut control),
        Err(StressError3::PointBudget)));
    assert_eq!(study.operator().elasticity().scales(), original);
    for stop in ["sdf3-stress-peak-scale", "sdf3-stress-adjoint", "sdf3-stress-publish"] {
        let mut poll = |p: SolveProgress| if p.stage == stop { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        assert!(matches!(study.evaluate_stress_with_measure(&rho, &loads, Default::default(), PEAK, &mut control),
            Err(StressError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(study.operator().elasticity().scales(), original);
    }
}

#[test]
fn peak_design_accepts_real_updates_and_restores_the_same_constraint_and_fields() {
    let (mut study, _, z) = fixture(0, false);
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = study.evaluate_stress_with_measure(&[0.55], &loads, Default::default(), PEAK, &mut control).unwrap().aggregate;
    let options = StressDesignOptions3 {
        stress_limit: cap, density_floor: 0.25,
        optimizer: ProjectedAlOptions { tolerance: 2e-6, max_evaluations: 2000, ..Default::default() },
        ..Default::default()
    };
    let mut session = StressDesignStudy3::new_with_measure(&mut study, &loads, &[0.85], options, PEAK, &mut control).unwrap();
    session.run(2).unwrap();
    assert!(session.history().len() > 1, "must retain real accepted material updates");
    assert_eq!(session.stress_measure(), PEAK);
    let checkpoint = session.checkpoint();
    session.run(2).unwrap();
    let expected = session.accepted().clone();
    let history = session.history().to_vec();
    let best = session.best_feasible().unwrap().clone();
    assert!(best.volume_fraction < session.history()[0].volume_fraction);
    bound(&best, options.stress.aggregation_power);
    assert!(best.sampled_physical_max <= cap * (1.0 + options.optimizer.tolerance));
    drop(session);
    let (mut other, _, _) = fixture(0, false);
    let mut poll = |_| ControlFlow::Continue(());
    let mut resumed_control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut restored = StressDesignStudy3::restore_with_measure(&mut other, &loads, checkpoint.clone(),
        options, PEAK, &mut resumed_control).unwrap();
    restored.run(2).unwrap();
    assert_eq!(restored.accepted().rho, expected.rho);
    assert_eq!(restored.accepted().displacements, expected.displacements);
    assert_eq!(restored.accepted().gradient, expected.gradient);
    assert_eq!(restored.history(), history);
    drop(restored);
    let original = other.operator().elasticity().scales().to_vec();
    assert!(StressDesignStudy3::restore(&mut other, &loads, checkpoint, options, &mut resumed_control).is_err(),
        "a peak checkpoint cannot silently restore as a normalized-average problem");
    assert_eq!(other.operator().elasticity().scales(), original);
}
