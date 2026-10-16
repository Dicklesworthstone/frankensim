//! Real Cartesian/adaptive CutFEM regressions for the three-field design map.
use super::*;
use fs_topopt::sdf3::PhysicalRegion3;
use fs_topopt::sdf3::robust::{RobustOptions3, controlled_sdf3_robust_gradient_check};

mod optimize;

#[test]
fn zero_offset_matches_nominal_fields_gradients_and_independent_loads() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho: Vec<f64> = (0..study.cells()).map(|i| 0.35 + 0.03*i as f64).collect();
    let loads = [LoadCase { force: &y, weight: 2.0 }, LoadCase { force: &z, weight: 3.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let nominal = study.evaluate(&rho, &loads, &mut control).unwrap();
    let scales = study.operator().scales().to_vec();
    let before = control.work();
    let actual = study.evaluate_robust(&rho, &loads, RobustOptions3 { threshold_offset: 0.0 }, &mut control).unwrap();
    assert_eq!(actual.eroded.objective.compliance.to_bits(), nominal.objective.compliance.to_bits());
    assert_eq!(actual.eroded.objective.gradient, nominal.objective.gradient);
    assert_eq!(actual.eroded.objective.displacements, nominal.objective.displacements);
    assert_eq!(actual.eroded.projected_rho, nominal.projected_rho);
    assert_eq!(actual.nominal_rho, actual.dilated_rho);
    assert_eq!(actual.nominal_rho, nominal.projected_rho);
    assert_eq!(actual.dilated_volume_gradient, nominal.volume_gradient);
    assert_eq!(actual.nominal_volume_gradient, nominal.volume_gradient);
    assert_eq!(actual.work.linear_solves - before.linear_solves, loads.len() + 5);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn shifted_realizations_match_separate_physical_models_not_nominal_sensitivities() {
    let params = SimpParams { beta: 5.0, ..Default::default() };
    let (mut study, y, z) = fixture(params);
    let rho: Vec<f64> = (0..study.cells()).map(|i| 0.4 + 0.02*i as f64).collect();
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let before = study.operator().scales().to_vec();
    let options = RobustOptions3::default();
    let actual = study.evaluate_robust(&rho, &loads, options, &mut control).unwrap();
    let (mut eroded, _, _) = fixture(SimpParams { eta: params.eta + options.threshold_offset, ..params });
    let reference = eroded.evaluate(&rho, &loads, &mut control).unwrap();
    let (dilated, _, _) = fixture(SimpParams { eta: params.eta - options.threshold_offset, ..params });
    let (volume, gradient) = dilated.volume_and_gradient(&rho, &mut control).unwrap();
    assert_eq!(actual.eroded.objective.compliance.to_bits(), reference.objective.compliance.to_bits());
    assert_eq!(actual.eroded.objective.gradient, reference.objective.gradient);
    assert_eq!(actual.eroded.projected_rho, reference.projected_rho);
    assert_eq!(actual.dilated_volume_fraction.to_bits(), volume.to_bits());
    assert_eq!(actual.dilated_volume_gradient, gradient);
    assert!(actual.eroded.volume_fraction < actual.nominal_volume_fraction);
    assert!(actual.nominal_volume_fraction < actual.dilated_volume_fraction);
    for ((e, n), d) in actual.eroded.projected_rho.iter().zip(&actual.nominal_rho).zip(&actual.dilated_rho) {
        assert!(e <= n && n <= d);
    }
    assert_eq!(study.operator().scales(), before);
    assert_eq!(study.params().eta, params.eta);
}

#[test]
fn full_eroded_compliance_and_dilated_volume_chain_matches_coordinate_differences() {
    for (penal, beta) in [(1.0, 1.0), (3.0, 7.0)] {
        let (mut study, y, z) = fixture(SimpParams { penal, beta, ..Default::default() });
        let rho: Vec<f64> = (0..study.cells()).map(|i| 0.37 + 0.025*i as f64).collect();
        let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let options = RobustOptions3::default();
        let exact = study.evaluate_robust(&rho, &loads, options, &mut control).unwrap();
        for i in 0..rho.len() {
            let mut plus = rho.clone(); let mut minus = rho.clone();
            let h = 1e-4; plus[i] += h; minus[i] -= h;
            let a = study.evaluate_robust(&plus, &loads, options, &mut control).unwrap();
            let b = study.evaluate_robust(&minus, &loads, options, &mut control).unwrap();
            let dc = (a.eroded.objective.compliance - b.eroded.objective.compliance)/(2.0*h);
            let dv = (a.dilated_volume_fraction - b.dilated_volume_fraction)/(2.0*h);
            assert!((dc - exact.eroded.objective.gradient[i]).abs()/dc.abs() < 3e-4,
                "penal={penal} beta={beta} cell={i}");
            assert!((dv - exact.dilated_volume_gradient[i]).abs()/dv.abs() < 2e-5);
        }
    }
}

#[test]
fn protected_physical_regions_survive_all_realizations_without_freezing_filter_controls() {
    let (study, y, z) = fixture(SimpParams { beta: 3.0, ..Default::default() });
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let rho = vec![0.5; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let actual = study.evaluate_robust(&rho, &loads, RobustOptions3::default(), &mut control).unwrap();
    for field in [&actual.eroded.projected_rho, &actual.nominal_rho, &actual.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
    assert!(actual.dilated_volume_gradient[0] > 0.0, "protected raw controls still influence adjacent free cells");
    let check = controlled_sdf3_robust_gradient_check(&mut study, &loads, &rho,
        RobustOptions3::default(), Default::default(), &mut control).unwrap();
    assert!(check.numerical.passed(), "{check:?}");
}

#[test]
fn cancelled_publication_restores_incoming_scales_without_refunding_the_solve() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho = vec![0.5; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut proceed = |_| ControlFlow::Continue(());
    let mut initial = SolveControl::new(SolveBudget::default(), &mut proceed);
    study.evaluate(&rho, &loads, &mut initial).unwrap();
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-robust-publish" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.evaluate_robust(&rho, &loads, RobustOptions3::default(), &mut control),
        Err(EvaluationStop::Cancelled)));
    assert_eq!(study.operator().scales(), scales);
    assert!(control.work().linear_iterations > 0);
    assert_eq!(control.work().linear_solves, loads.len() + 5);
}

#[test]
fn adaptive_hanging_node_geometry_uses_the_same_robust_gradient_map() {
    use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
    use fs_cutfem::octree3::Octree3;
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree.leaves().iter().find(|c| c.index() == [0, 0, 1]).unwrap();
    let tree = tree.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut qc = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let op = AdaptiveElasticity3::build(HexCell::try_new([0.0;3], [1.0;3]).unwrap(), &tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut qc).unwrap();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let force = op.body_load(&|p| [0.0, -p[0], -1.0], || ControlFlow::Continue(())).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.15, SimpParams { beta: 4.0, ..Default::default() });
    let rho: Vec<f64> = (0..study.cells()).map(|i| 0.42 + 0.01*(i%7) as f64).collect();
    let loads = [LoadCase { force: &force, weight: 1.0 }];
    let before = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let check = controlled_sdf3_robust_gradient_check(&mut study, &loads, &rho,
        RobustOptions3::default(), Default::default(), &mut control).unwrap();
    assert!(check.numerical.passed(), "{check:?}");
    assert_eq!(study.operator().scales(), before);
}
