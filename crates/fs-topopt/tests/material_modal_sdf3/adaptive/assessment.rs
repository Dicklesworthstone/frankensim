//! The two grids assess one physical design and retain their incoming state.
use super::*;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_topopt::sdf3::AdaptiveSdf3Elasticity;
use fs_topopt::sdf3::modal::design::Sdf3EigenfrequencyGoal3;
use fs_topopt::sdf3::modal::refinement::{
    Sdf3EigenfrequencyAdaptiveError3, Sdf3EigenfrequencyRefinementOptions3,
};

#[test]
fn enrichment_inherits_physical_fields_and_rolls_back_on_every_exit() {
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree.leaves().iter().find(|leaf| leaf.index() == [1, 0, 0]).unwrap();
    let fine_tree = tree.refined(&[mark], poll).unwrap();
    let mut operator = build(&tree, 3.0, &mut poll).unwrap();
    let source_scales = vec![0.37; operator.cells()];
    operator.set_scales(&source_scales).unwrap();
    let mass = operator.prepare_material_mass(100_000, poll).unwrap();
    let mut regions = vec![PhysicalRegion3::Design; operator.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let mut study = CutDensityStudy3::new(operator, 0.25, params())
        .with_physical_regions(regions.clone()).unwrap();
    let rho: Vec<_> = (0..study.cells()).map(|i| 0.45 + 0.03 * (i % 5) as f64).collect();
    let mut fine = build(&fine_tree, 3.0, &mut poll).unwrap();
    let fine_scales = vec![0.61; fine.cells()];
    fine.set_scales(&fine_scales).unwrap();
    let parents = AdaptiveTransfer3::new(
        study.operator().adaptive(), fine.adaptive(), 100_000, poll,
    ).unwrap().parents().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = study.evaluate_eigenfrequency(&mass, &rho, frequency(), &mut control).unwrap();
    let options = design_settings(baseline.objective.aggregate);
    let goal = Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit: 1.0 };
    let report = study.estimate_eigenfrequency_enrichment(
        &mass, &mut fine, &rho, goal, options,
        Sdf3EigenfrequencyRefinementOptions3::default(), &mut control,
    ).unwrap();
    assert_eq!(report.coarse.physical_rho, baseline.projected_rho);
    assert_eq!(report.fine.physical_rho.len(), fine.cells());
    for (child, &parent) in parents.iter().enumerate() {
        assert_eq!(report.fine.physical_rho[child], report.coarse.physical_rho[parent]);
        assert_eq!(report.fine.scales[child], report.coarse.scales[parent]);
    }
    assert_eq!(report.coarse.physical_rho[0], 1.0);
    close(
        report.fine.physical_mass_kg,
        dot(fine.reference_cell_masses_kg().unwrap(), &report.fine.physical_rho),
        1e-12,
    );
    assert!(report.subspace_change.is_finite() && (0.0..=1.0 + 1e-8).contains(&report.subspace_change));
    assert!(report.eigenvalue_change.is_finite());
    assert!(report.marking_mass.values().all(|v| v.is_finite() && *v >= 0.0));
    assert!(report.mark(0.5, 4, poll).is_ok());
    assert_eq!(study.operator().scales(), source_scales);
    assert_eq!(fine.scales(), fine_scales);

    // A real refit must rebuild the filter. That is a DIFFERENT design from
    // the physical parent-cell copy used to assess discretization change.
    let rebuilt = build(&fine_tree, 3.0, &mut poll).unwrap();
    let rebuilt_mass = rebuilt.prepare_material_mass(100_000, poll).unwrap();
    let mut refiltered = CutDensityStudy3::new(rebuilt, 0.25, params())
        .with_physical_regions(parents.iter().map(|&i| regions[i]).collect()).unwrap();
    let inherited: Vec<_> = parents.iter().map(|&i| rho[i]).collect();
    let changed = refiltered.evaluate_eigenfrequency(
        &rebuilt_mass, &inherited, frequency(), &mut control,
    ).unwrap();
    assert!(changed.projected_rho.iter().zip(&report.fine.physical_rho)
        .any(|(a, b)| (a - b).abs() > 1e-5));

    let mut foreign = build(&fine_tree, 4.0, &mut poll).unwrap();
    let foreign_scales = vec![0.83; foreign.cells()];
    foreign.set_scales(&foreign_scales).unwrap();
    assert!(matches!(
        study.estimate_eigenfrequency_enrichment(
            &mass, &mut foreign, &rho, goal, options,
            Sdf3EigenfrequencyRefinementOptions3::default(), &mut control,
        ),
        Err(Sdf3EigenfrequencyAdaptiveError3::Refinement(GoalRefinementError3::Physics(_))),
    ));
    assert_eq!(study.operator().scales(), source_scales);
    assert_eq!(foreign.scales(), foreign_scales);

    let mut stop = |progress: SolveProgress| {
        if progress.stage == "sdf3-eigen-refinement-fields" { ControlFlow::Break(()) } else { poll() }
    };
    let mut cancelled = SolveControl::new(SolveBudget::default(), &mut stop);
    assert!(matches!(
        study.estimate_eigenfrequency_enrichment(
            &mass, &mut fine, &rho, goal, options,
            Sdf3EigenfrequencyRefinementOptions3::default(), &mut cancelled,
        ),
        Err(Sdf3EigenfrequencyAdaptiveError3::Refinement(
            GoalRefinementError3::Evaluation(EvaluationStop::Cancelled),
        )),
    ));
    assert_eq!(study.operator().scales(), source_scales);
    assert_eq!(fine.scales(), fine_scales);
}
