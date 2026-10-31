//! Physical mass can change across grids even when geometric volume agrees.
use std::{cell::Cell, ops::ControlFlow};

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::ElasticMaterialField3;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::limits::ResponseLimitSense3;
use fs_topopt::sdf3::response::limits::design::{
    MinimumMaterialObjective3, MinimumMaterialOptions3,
};
use fs_topopt::sdf3::response::limits::refinement::adaptive::{
    AdaptiveMinimumMaterialOptions3, AdaptiveMinimumMaterialStop3, AdaptiveMinimumMaterialStudy3,
};
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3, ReferenceResponseLimitCase3,
    ResponseLimitRefinementOptions3,
};
use fs_topopt::{SimpParams, SolveBudget, SolveControl};

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;

// The odd height perturbation has equal and opposite integration errors on
// the two material regions. Their different densities prevent kg cancellation.
struct CurvedTop;
impl CutSdf3 for CurvedTop {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.45 - 2.0 * (p[1] - 0.5).powi(7)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let y = Interval::new(lo[1], hi[1]) - Interval::new(0.5, 0.5);
        let y2 = y * y;
        Interval::new(lo[2], hi[2])
            - Interval::new(0.45, 0.45)
            - Interval::new(2.0, 2.0) * y2 * y2 * y2 * y
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        match axis {
            HeightAxis::Z => Interval::new(1.0, 1.0),
            HeightAxis::X => Interval::new(0.0, 0.0),
            HeightAxis::Y => {
                let y = Interval::new(lo[1], hi[1]) - Interval::new(0.5, 0.5);
                let y2 = y * y;
                Interval::new(-14.0, -14.0) * y2 * y2 * y2
            }
        }
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn tree(level: u8) -> Octree3 {
    Octree3::uniform(level, 4, 4096).unwrap()
}
fn build(tree: &Octree3, native: bool) -> AdaptiveSolveSpace3 {
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            root_relative_tolerance: 1e-12,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let domain = HexCell::try_new([0.0; 3], [1.0; 3]).unwrap();
    let op = if native {
        let field = |cell: HexCell| {
            let id = if cell.lo()[1] < 0.5 { 2 } else { 8 };
            let c =
                std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
            Ok(
                ElasticAcousticState::try_new_mandel(c, f64::from(id), ContentHash([id; 32]), 0.1)
                    .unwrap(),
            )
        };
        AdaptiveElasticity3::build_with_material_field(
            domain,
            tree,
            &CurvedTop,
            ElasticMaterialField3::ByCell(&field),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut q,
        )
        .unwrap()
    } else {
        AdaptiveElasticity3::build(
            domain,
            tree,
            &CurvedTop,
            &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut q,
        )
        .unwrap()
    };
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}
fn study(tree: &Octree3, native: bool) -> Study {
    CutDensityStudy3::new(
        build(tree, native),
        0.1,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.02,
            eta: 0.5,
        },
    )
}
fn limits() -> [ReferenceResponseLimit3<'static>; 1] {
    [ReferenceResponseLimit3 {
        response: ReferenceBoundedResponse3::Displacement(ReferenceLoad3::default()),
        sense: ResponseLimitSense3::AtMost,
        bound: 1.0,
        scale: 1.0,
    }]
}
fn cases<'a>(limits: &'a [ReferenceResponseLimit3<'a>]) -> [ReferenceResponseLimitCase3<'a>; 1] {
    [ReferenceResponseLimitCase3 {
        load: ReferenceLoad3::default(),
        prescribed: None,
        limits,
    }]
}
fn mass_options() -> MinimumMaterialOptions3 {
    MinimumMaterialOptions3 {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
        ..Default::default()
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn equal_geometric_volume_does_not_hide_native_material_mass_error() {
    let mut s = study(&tree(1), true);
    let fine = build(&tree(2), true);
    close(
        s.operator().elasticity().volumes().iter().sum(),
        0.45,
        1e-11,
    );
    close(fine.elasticity().volumes().iter().sum(), 0.45, 1e-11);
    let rows = limits();
    let cases = cases(&rows);
    let rho = vec![0.7; s.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let fit = s
        .fit_reference_minimum_material(&cases, &rho, mass_options(), 0, &mut control)
        .unwrap();
    let before = s.operator().elasticity().scales().to_vec();
    let q = fit.accepted.projected_rho[0];
    let exact = q * (2.0 * (0.225 - 1.0 / 1024.0) + 8.0 * (0.225 + 1.0 / 1024.0));
    let report = s
        .estimate_response_limit_enrichment(
            fine,
            &fit.accepted,
            &cases,
            ResponseLimitRefinementOptions3 {
                objective: mass_options().objective,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    let coarse = report.coarse_mass_kg.unwrap();
    let fine = report.fine_mass_kg.unwrap();
    close(coarse, fit.accepted.mass_kg.unwrap(), 1e-12);
    assert!(fine - coarse > 1e-5);
    assert!((fine - exact).abs() < 0.25 * (coarse - exact).abs());
    close(
        report.coarse_volume_fraction,
        report.fine_volume_fraction,
        1e-11,
    );
    close(report.max_normalized_change(), 0.0, 1e-13);
    close(report.coarse_objective, coarse / 2.0, 1e-13);
    close(report.fine_objective, fine / 2.0, 1e-13);
    close(
        report.mass_cell_corrections_kg.values().sum(),
        fine - coarse,
        1e-11,
    );
    assert!(report.mass_cell_corrections_kg.values().any(|x| *x > 0.0));
    assert!(report.mass_cell_corrections_kg.values().any(|x| *x < 0.0));
    for (leaf, correction) in &report.mass_cell_corrections_kg {
        close(report.marking_mass[leaf], correction.abs() / 2.0, 1e-12);
    }
    assert!(!report.mark(0.5, 4, poll).unwrap().marked.is_empty());
    assert_eq!(s.operator().elasticity().scales(), before);
}

#[test]
fn mass_controls_stopping_and_survives_selective_refit_while_volume_remains_default() {
    for mass in [false, true] {
        let background = tree(1);
        let mut s = study(&background, true);
        let rows = limits();
        let cases = cases(&rows);
        let rho = vec![0.7; s.cells()];
        let fit = if mass {
            mass_options()
        } else {
            MinimumMaterialOptions3::default()
        };
        let options = AdaptiveMinimumMaterialOptions3 {
            fit,
            steps_per_grid: 0,
            max_refinements: 1,
            theta: 0.25,
            max_marks: 1,
            mass_tolerance: 1e-10,
            ..Default::default()
        };
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let mut campaign = AdaptiveMinimumMaterialStudy3::new(
            &mut s,
            background,
            &cases,
            &rho,
            options,
            &mut control,
        )
        .unwrap();
        let stop = campaign
            .run(0, |t, _| Ok::<_, ()>(build(t, true)), &mut control)
            .unwrap();
        let assessed = campaign.assessment().unwrap().unwrap();
        assert!(assessed.fine_feasible);
        assert!(assessed.response_change <= options.response_tolerance);
        assert!(assessed.volume_change <= options.volume_tolerance);
        if !mass {
            assert_eq!(stop, AdaptiveMinimumMaterialStop3::TwoGridLimitsSatisfied);
            assert_eq!(assessed.mass_change_kg, None);
            let probe = campaign.probe().unwrap();
            assert_eq!(probe.objective, MinimumMaterialObjective3::VolumeFraction);
            assert!(probe.mass_cell_corrections_kg.is_empty());
            close(probe.coarse_objective, probe.coarse_volume_fraction, 1e-13);
            continue;
        }
        assert_eq!(stop, AdaptiveMinimumMaterialStop3::Paused);
        assert!(assessed.mass_change_kg.unwrap() > 1e-5);
        assert!(assessed.objective_change > options.mass_tolerance);
        assert_eq!(campaign.probe().unwrap().objective, fit.objective);
        let stop = campaign
            .run(1, |t, _| Ok::<_, ()>(build(t, true)), &mut control)
            .unwrap();
        assert_eq!(stop, AdaptiveMinimumMaterialStop3::RefinementLimit);
        assert_eq!(campaign.rounds().len(), 1);
        let round = &campaign.rounds()[0];
        assert!(round.cells[1] > round.cells[0]);
        for (objective, kg) in round.objectives.into_iter().zip(round.masses_kg.unwrap()) {
            close(objective, kg / 2.0, 1e-13);
        }
        assert_eq!(campaign.probe().unwrap().objective, fit.objective);
        close(
            campaign.probe().unwrap().coarse_mass_kg.unwrap(),
            campaign.accepted().mass_kg.unwrap(),
            1e-12,
        );
        for row in &campaign.fit().history {
            close(row.objective, row.mass_kg.unwrap() / 2.0, 1e-13);
        }
    }
}

#[test]
fn changed_objective_or_missing_native_inventory_is_refused_before_reference_callbacks() {
    let mut s = study(&tree(1), true);
    let rows = limits();
    let base = cases(&rows);
    let rho = vec![0.7; s.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let fit = s
        .fit_reference_minimum_material(&base, &rho, mass_options(), 0, &mut control)
        .unwrap();
    let scales = s.operator().elasticity().scales().to_vec();
    let calls = Cell::new(0);
    let zero = |_| {
        calls.set(calls.get() + 1);
        [0.0; 3]
    };
    let counted = [ReferenceResponseLimitCase3 {
        load: ReferenceLoad3::body(&zero),
        prescribed: None,
        limits: &rows,
    }];
    for objective in [
        MinimumMaterialObjective3::VolumeFraction,
        MinimumMaterialObjective3::MassKg { scale_kg: 0.0 },
    ] {
        assert!(
            s.estimate_response_limit_enrichment(
                build(&tree(2), true),
                &fit.accepted,
                &counted,
                ResponseLimitRefinementOptions3 {
                    objective,
                    ..Default::default()
                },
                &mut control
            )
            .is_err()
        );
        assert_eq!(calls.get(), 0);
    }
    assert!(
        s.refit_reference_minimum_material(
            study(&tree(2), true),
            &fit.accepted,
            &counted,
            MinimumMaterialOptions3::default(),
            0,
            2_000_000,
            &mut control
        )
        .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert!(
        s.refit_reference_minimum_material(
            study(&tree(2), false),
            &fit.accepted,
            &counted,
            mass_options(),
            0,
            2_000_000,
            &mut control
        )
        .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert_eq!(s.operator().elasticity().scales(), scales);
}
