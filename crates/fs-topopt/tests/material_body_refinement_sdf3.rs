//! G1/G3: native self-weight survives goal estimation and selective material refits.
use std::{cell::Cell, ops::ControlFlow};

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::ElasticMaterialField3;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_topopt::sdf3::response::limits::ResponseLimitSense3;
use fs_topopt::sdf3::response::limits::design::{
    MinimumMaterialObjective3, MinimumMaterialOptions3,
};
use fs_topopt::sdf3::response::limits::refinement::adaptive::{
    AdaptiveMinimumMaterialOptions3, AdaptiveMinimumMaterialStop3, AdaptiveMinimumMaterialStudy3,
    MinimumMaterialBuild3,
};
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3, ReferenceResponseLimitCase3,
    ResponseLimitRefinementOptions3,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::sdf3_goal::GoalRefinementError3;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        (p[0] - 0.17) * (p[0] - 0.83)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X {
            Interval::new(2.0, 2.0) * Interval::new(lo[0], hi[0]) - Interval::new(1.0, 1.0)
        } else {
            Interval::new(0.0, 0.0)
        }
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn acceleration(p: [f64; 3]) -> [f64; 3] {
    [0.004 * (1.0 + p[1]), 0.001, -0.003]
}
fn dead(_: [f64; 3]) -> [f64; 3] {
    [0.0002, 0.0, 0.0]
}
fn observe(p: [f64; 3]) -> [f64; 3] {
    [1.0 + p[1], 0.5, -0.25]
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [0.003 * (p[0] - 0.17), 0.0, 0.0]
}
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    if n[0] > 0.0 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0; 3]
    }
}
fn tree(level: u8) -> Octree3 {
    Octree3::uniform(level, 4, 4096).unwrap()
}
fn build(tree: &Octree3) -> AdaptiveSolveSpace3 {
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            root_relative_tolerance: 1e-12,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    let field = |cell: HexCell| {
        let id = if cell.lo()[1] < 0.5 { 2 } else { 8 };
        let c = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
        Ok(
            ElasticAcousticState::try_new_mandel(c, f64::from(id), ContentHash([id; 32]), 0.25)
                .unwrap(),
        )
    };
    let op = AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        tree,
        &Slab,
        ElasticMaterialField3::ByCell(&field),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut q,
    )
    .unwrap();
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}
fn study(tree: &Octree3) -> Study {
    let s = CutDensityStudy3::new(
        build(tree),
        0.15,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.03,
            eta: 0.5,
        },
    );
    let mut regions = vec![PhysicalRegion3::Design; s.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = PhysicalRegion3::Void;
    s.with_physical_regions(regions).unwrap()
}
fn raw(s: &Study) -> Vec<f64> {
    (0..s.cells()).map(|i| 0.33 + 0.055 * i as f64).collect()
}
fn limits() -> [ReferenceResponseLimit3<'static>; 2] {
    [
        ReferenceResponseLimit3 {
            response: ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
            sense: ResponseLimitSense3::AtMost,
            bound: 10.0,
            scale: 1.0,
        },
        ReferenceResponseLimit3 {
            response: ReferenceBoundedResponse3::Reaction(&right),
            sense: ResponseLimitSense3::AtMost,
            bound: 10.0,
            scale: 1.0,
        },
    ]
}
fn cases<'a>(limits: &'a [ReferenceResponseLimit3<'a>]) -> [ReferenceResponseLimitCase3<'a>; 1] {
    [ReferenceResponseLimitCase3 {
        load: ReferenceLoad3::body(&dead),
        prescribed: Some(&motion),
        limits,
    }]
}
fn fit_options() -> MinimumMaterialOptions3 {
    MinimumMaterialOptions3 {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
        ..Default::default()
    }
}
fn probe_options() -> ResponseLimitRefinementOptions3 {
    ResponseLimitRefinementOptions3 {
        objective: fit_options().objective,
        ..Default::default()
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * a.abs().max(b.abs()).max(1e-6),
        "{a} != {b}"
    );
}
// Independently prepare the actual native load on an independently built owner.
// Removing self-weight must fail this force-balance check by a visible margin.
fn audit(op: &mut AdaptiveSolveSpace3, rho: &[f64], u: &[f64], p: SimpParams) -> [f64; 2] {
    let scales: Vec<_> = rho
        .iter()
        .map(|&r| p.e_min + (1.0 - p.e_min) * fs_math::det::pow(r, p.penal))
        .collect();
    op.set_scales(&scales).unwrap();
    let native = op.elasticity();
    let body = native
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let force = body.load(rho, poll).unwrap();
    let dead = native
        .reference_load_with_motion(ReferenceLoad3::body(&dead), Some(&motion), poll)
        .unwrap();
    let rhs: Vec<_> = dead.iter().zip(force).map(|(a, b)| a + b).collect();
    assert!(native.field_residual(u, &rhs, poll).unwrap() < 1e-8);
    assert!(native.field_residual(u, &dead, poll).unwrap() > 1e-6);
    let q = native
        .reference_load(ReferenceLoad3::body(&observe), poll)
        .unwrap();
    [
        q.iter().zip(u).map(|(q, u)| q * u).sum(),
        native
            .embedded_reaction(u, Some(&motion), &right, poll)
            .unwrap()
            .value,
    ]
}

#[test]
fn self_weight_enrichment_reconstructs_affine_goals_and_native_force_balance() {
    let coarse_tree = tree(1);
    let mut s = study(&coarse_tree)
        .with_material_acceleration(acceleration, 100_000, poll)
        .unwrap();
    let rows = limits();
    let cases = cases(&rows);
    let rho = raw(&s);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let accepted = s
        .fit_reference_minimum_material(&cases, &rho, fit_options(), 0, &mut control)
        .unwrap()
        .accepted;
    let scales = s.operator().elasticity().scales().to_vec();
    let fingerprint = s.material_body_load_fingerprint();
    let report = s
        .estimate_response_limit_enrichment(
            build(&tree(2)),
            &accepted,
            &cases,
            probe_options(),
            &mut control,
        )
        .unwrap();
    assert_eq!(s.operator().elasticity().scales(), scales);
    assert_eq!(s.material_body_load_fingerprint(), fingerprint);
    let mut fine = build(&tree(2));
    let parents = AdaptiveTransfer3::new(
        s.operator().elasticity(),
        fine.elasticity(),
        2_000_000,
        poll,
    )
    .unwrap()
    .parents()
    .to_vec();
    assert_eq!(
        report.fine_projected_rho,
        parents
            .iter()
            .map(|&i| accepted.projected_rho[i])
            .collect::<Vec<_>>()
    );
    assert!(report.fine_projected_rho.contains(&0.0));
    assert!(report.fine_projected_rho.contains(&1.0));
    let measured = audit(
        &mut fine,
        &report.fine_projected_rho,
        &report.cases[0].fine_displacement,
        s.params(),
    );
    for (j, change) in report.cases[0].limits.iter().enumerate() {
        close(change.goal.coarse_value, accepted.measurements[0][j]);
        close(change.goal.fine_value, measured[j]);
        close(
            change.goal.correction(),
            measured[j] - accepted.measurements[0][j],
        );
        close(
            change.correction(),
            change.fine_residual - change.coarse_residual,
        );
        close(
            change
                .goal
                .linearized
                .cells
                .values()
                .map(|c| c.signed_dwr)
                .sum(),
            change.goal.linearized.dwr,
        );
        assert!(change.goal.identity_relative_defect < 1e-8);
    }
    assert!(report.max_normalized_change() > 1e-10);
    assert!(report.marking_mass.values().sum::<f64>() > 0.0);
    let marked = report.mark(0.5, 2, poll).unwrap();
    assert!(!marked.marked.is_empty());
    assert!(
        marked
            .marked
            .iter()
            .all(|leaf| s.operator().elasticity().leaves().contains(leaf))
    );
    let kg: f64 = fine
        .elasticity()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&report.fine_projected_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(report.fine_mass_kg.unwrap(), kg);
    close(report.fine_objective, kg / 2.0);
}

#[test]
fn adaptive_mass_design_reprepares_self_weight_on_probe_and_selective_refit() {
    let background = tree(1);
    let mut s = study(&background)
        .with_material_acceleration(acceleration, 100_000, poll)
        .unwrap();
    let initial_cells = s.cells();
    let initial_fingerprint = s.material_body_load_fingerprint();
    let rows = limits();
    let cases = cases(&rows);
    let rho = raw(&s);
    let options = AdaptiveMinimumMaterialOptions3 {
        fit: fit_options(),
        steps_per_grid: 0,
        max_refinements: 1,
        max_marks: 1,
        response_tolerance: 0.0,
        volume_tolerance: 1.0,
        mass_tolerance: 1.0,
        ..Default::default()
    };
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let mut campaign =
        AdaptiveMinimumMaterialStudy3::new(&mut s, background, &cases, &rho, options, &mut control)
            .unwrap();
    let mut builds = Vec::new();
    let stop = campaign
        .run(
            0,
            |tree, kind| {
                builds.push(kind);
                Ok::<_, ()>(build(tree))
            },
            &mut control,
        )
        .unwrap();
    assert_eq!(stop, AdaptiveMinimumMaterialStop3::Paused);
    assert_eq!(builds, [MinimumMaterialBuild3::Probe]);
    assert!(campaign.assessment().unwrap().unwrap().response_change > 1e-10);
    let stop = campaign
        .run(
            1,
            |tree, kind| {
                builds.push(kind);
                Ok::<_, ()>(build(tree))
            },
            &mut control,
        )
        .unwrap();
    assert_eq!(stop, AdaptiveMinimumMaterialStop3::RefinementLimit);
    assert_eq!(
        builds,
        [
            MinimumMaterialBuild3::Probe,
            MinimumMaterialBuild3::Refit,
            MinimumMaterialBuild3::Probe
        ]
    );
    assert_eq!(campaign.rounds().len(), 1);
    assert!(campaign.study().cells() > initial_cells);
    assert!(campaign.study().has_material_body_load());
    assert_ne!(
        campaign.study().material_body_load_fingerprint(),
        initial_fingerprint
    );
    let selected = campaign.accepted();
    assert!(selected.feasible(options.fit.optimizer.tolerance));
    let mut audit_op = build(campaign.tree());
    let measured = audit(
        &mut audit_op,
        &selected.projected_rho,
        &selected.displacements[0],
        campaign.study().params(),
    );
    for (actual, expected) in selected.measurements[0].iter().zip(measured) {
        close(*actual, expected);
    }
    let labels = campaign.study().physical_regions().unwrap();
    assert!(labels.contains(&PhysicalRegion3::Solid));
    assert!(labels.contains(&PhysicalRegion3::Void));
    for (label, &r) in labels.iter().zip(&selected.projected_rho) {
        match label {
            PhysicalRegion3::Solid => assert_eq!(r, 1.0),
            PhysicalRegion3::Void => assert_eq!(r, 0.0),
            _ => {}
        }
    }
    for (value, kg) in campaign.rounds()[0]
        .objectives
        .into_iter()
        .zip(campaign.rounds()[0].masses_kg.unwrap())
    {
        close(value, kg / 2.0);
    }
    close(
        campaign.probe().unwrap().coarse_mass_kg.unwrap(),
        selected.mass_kg.unwrap(),
    );
    assert_eq!(campaign.probe().unwrap().objective, options.fit.objective);
}

#[test]
fn fixed_stencil_refusal_and_cancelled_recipe_probe_do_not_solve_or_change_source() {
    let background = tree(1);
    let fixed = study(&background);
    let body = fixed
        .operator()
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let mut fixed = fixed.with_material_body_load(body).unwrap();
    let rows = limits();
    let cases = cases(&rows);
    let rho = raw(&fixed);
    let before = fixed.operator().elasticity().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let work = control.work();
    assert!(
        AdaptiveMinimumMaterialStudy3::new(
            &mut fixed,
            background,
            &cases,
            &rho,
            AdaptiveMinimumMaterialOptions3 {
                fit: fit_options(),
                steps_per_grid: 0,
                ..Default::default()
            },
            &mut control
        )
        .is_err()
    );
    assert_eq!(control.work(), work);
    assert_eq!(fixed.operator().elasticity().scales(), before);

    let mut s = study(&tree(1))
        .with_material_acceleration(acceleration, 100_000, poll)
        .unwrap();
    let accepted = s
        .fit_reference_minimum_material(&cases, &rho, fit_options(), 0, &mut control)
        .unwrap()
        .accepted;
    let scales = s.operator().elasticity().scales().to_vec();
    let fingerprint = s.material_body_load_fingerprint();
    let stages = Cell::new(0);
    let mut cancel = |p: SolveProgress| {
        if p.stage == "response-limit-goal-start" {
            stages.set(stages.get() + 1);
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut stopped = SolveControl::new(SolveBudget::default(), &mut cancel);
    let work = stopped.work();
    assert!(matches!(
        s.estimate_response_limit_enrichment(
            build(&tree(2)),
            &accepted,
            &cases,
            probe_options(),
            &mut stopped
        ),
        Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))
    ));
    assert_eq!(stages.get(), 1);
    assert_eq!(stopped.work(), work);
    assert_eq!(s.operator().elasticity().scales(), scales);
    assert_eq!(s.material_body_load_fingerprint(), fingerprint);
}
