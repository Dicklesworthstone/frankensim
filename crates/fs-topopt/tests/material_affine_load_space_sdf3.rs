//! G1/G3: an affine uncertainty ball includes native material self-weight.
use std::{cell::RefCell, ops::ControlFlow};

#[path = "material_affine_load_space_sdf3/response_refinement.rs"]
mod response_refinement;

#[path = "material_affine_load_space_sdf3/stress.rs"]
mod stress;

use fs_ascent::quadratic_ball::QuadraticBallBranch;
use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::sdf3::response::limits::{BoundedResponse3, ResponseLimit3, ResponseLimitSense3};
use fs_topopt::sdf3::robust::RobustOptions3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::campaign::{
    AffineLoadSpaceBuild3, AffineLoadSpaceCampaign3, AffineLoadSpaceCampaignOptions3,
    AffineLoadSpaceCampaignStop3, AffineLoadSpaceMarking3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::design::{
    AffineLoadSpaceDesign3, AffineLoadSpaceDesignGoal3, AffineLoadSpaceDesignOptions3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::refinement::AffineLoadSpaceRefinementOptions3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::refit::AffineLoadSpaceRefitOptions3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::{
    AffineResponseLimits3, AffineResponseOptions3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::{
    AffineLoadSpace3, AffineLoadSpaceEvaluation3,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::sdf3_goal::GoalRefinementError3;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn acceleration(p: [f64; 3]) -> [f64; 3] {
    [0.0, -0.003 * (1.0 + p[0]), -0.004]
}
fn native_material(cell: HexCell) -> Result<ElasticAcousticState, ElasticityError3> {
    let id = if cell.lo()[0] < 0.5 { 2 } else { 8 };
    let tensor = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
    Ok(
        ElasticAcousticState::try_new_mandel(tensor, f64::from(id), ContentHash([id; 32]), 0.25)
            .unwrap(),
    )
}
fn operator() -> CutElasticity3 {
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
    CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Top,
        ElasticMaterialField3::ByCell(&native_material),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut q,
    )
    .unwrap()
}
fn fixture() -> (CutDensityStudy3, Vec<f64>, Vec<f64>, Vec<f64>) {
    fixture_on_operator(operator())
}
fn fixture_on_operator(op: CutElasticity3) -> (CutDensityStudy3, Vec<f64>, Vec<f64>, Vec<f64>) {
    let first = op
        .body_load(&|p| [0.0, -0.001 * (1.0 + p[1]), 0.0], poll)
        .unwrap();
    let second = op
        .body_load(
            &|p| [0.0, -0.0006 * (1.0 + p[1]), -0.0008 * (1.0 + p[0])],
            poll,
        )
        .unwrap();
    let dead = op.body_load(&|_| [0.0002, -0.0001, -0.0003], poll).unwrap();
    let s = CutDensityStudy3::new(
        op,
        0.15,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.03,
            eta: 0.5,
        },
    )
    .with_material_acceleration(acceleration, 100_000, poll)
    .unwrap();
    let mut regions = vec![PhysicalRegion3::Design; s.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = PhysicalRegion3::Void;
    (
        s.with_physical_regions(regions).unwrap(),
        dead,
        first,
        second,
    )
}
fn raw(s: &CutDensityStudy3) -> Vec<f64> {
    (0..s.cells()).map(|i| 0.36 + 0.03 * i as f64).collect()
}
fn policy() -> AffineLoadSpace3 {
    AffineLoadSpace3::new(2, 1e-3).unwrap()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-12),
        "{a} != {b}"
    );
}
fn close_vector(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let error: f64 = a.iter().zip(b).map(|(a, b)| (a - b) * (a - b)).sum();
    assert!(
        error <= 1e-16 * dot(a, a).max(dot(b, b)).max(1e-30),
        "vector discrepancy {error}"
    );
}
fn audit(e: &AffineLoadSpaceEvaluation3, dead: &[f64], basis: &[&[f64]]) {
    let mut native = operator();
    native.set_scales(&e.physical.eroded_scales).unwrap();
    let body = native
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let self_weight = body.load(&e.physical.eroded_rho, poll).unwrap();
    let offset: Vec<_> = dead.iter().zip(&self_weight).map(|(a, b)| a + b).collect();
    close_vector(&offset, &e.offset_force);
    let force: Vec<_> = (0..native.n())
        .map(|d| {
            offset[d]
                + basis
                    .iter()
                    .zip(&e.worst_coefficients)
                    .map(|(f, q)| f[d] * q)
                    .sum::<f64>()
        })
        .collect();
    close_vector(&force, &e.worst_force);
    let u: Vec<_> = (0..native.n())
        .map(|d| {
            e.offset_displacement[d]
                + e.basis_displacements
                    .iter()
                    .zip(&e.worst_coefficients)
                    .map(|(u, q)| u[d] * q)
                    .sum::<f64>()
        })
        .collect();
    close_vector(&u, &e.worst_displacement);
    let mut ku = vec![0.0; native.n()];
    native.apply(&e.worst_displacement, &mut ku);
    let residual: f64 = ku.iter().zip(&force).map(|(a, b)| (a - b) * (a - b)).sum();
    assert!(residual.sqrt() <= 1e-8 * dot(&force, &force).sqrt());
    close(
        dot(&force, &e.worst_displacement),
        e.worst_compliance,
        1e-10,
    );
    close(
        dot(&e.worst_coefficients, &e.worst_coefficients),
        1.0,
        1e-10,
    );
    assert!(e.witness_relative_defect < 1e-8);
    assert!(e.reciprocity_defect < 1e-8);
}

#[test]
fn worst_total_field_obeys_native_equilibrium_and_dominates_dense_authored_directions() {
    let (mut s, dead, first, second) = fixture();
    let rho = raw(&s);
    let basis = [&first[..], &second[..]];
    let scales = s.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let e = s
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    assert_eq!(e.oracle.branch, QuadraticBallBranch::Regular);
    assert!(e.oracle.gradient_admissible());
    assert!(e.gradient.is_some());
    audit(&e, &dead, &basis);
    close(e.maximum, e.worst_compliance / policy().reference(), 1e-9);
    // Independent direction sampling uses the actual solved offset/basis fields,
    // not the oracle's eigensystem or its stationarity equation.
    let mut sampled = f64::NEG_INFINITY;
    for radius in [0.0, 0.5, 1.0] {
        for i in 0..720 {
            let theta = std::f64::consts::TAU * i as f64 / 720.0;
            let a = radius * theta.cos();
            let b = radius * theta.sin();
            let work: f64 = (0..dead.len())
                .map(|d| {
                    let f = e.offset_force[d] + a * first[d] + b * second[d];
                    let u = e.offset_displacement[d]
                        + a * e.basis_displacements[0][d]
                        + b * e.basis_displacements[1][d];
                    f * u
                })
                .sum();
            let value = work / policy().reference();
            let quadratic = e.offset
                + 2.0 * (a * e.linear[0] + b * e.linear[1])
                + a * a * e.gram[0]
                + a * b * (e.gram[1] + e.gram[2])
                + b * b * e.gram[3];
            close(value, quadratic, 1e-9);
            assert!(value <= e.maximum + 1e-9 * e.maximum.abs());
            sampled = sampled.max(value);
        }
    }
    assert!(sampled >= e.maximum * (1.0 - 2e-4));
    assert_eq!(s.operator().scales(), scales);
    for field in [
        &e.physical.eroded_rho,
        &e.physical.nominal_rho,
        &e.physical.dilated_rho,
    ] {
        assert_eq!(field[0], 1.0);
        assert_eq!(*field.last().unwrap(), 0.0);
    }
    assert_eq!(e.physical.dilated_mass_kg, None);
    assert_eq!(e.physical.dilated_mass_gradient, None);
}

#[test]
fn regular_branch_gradient_includes_both_self_weight_work_terms_through_filter_and_projection() {
    let (mut s, dead, first, second) = fixture();
    let rho = raw(&s);
    let basis = [&first[..], &second[..]];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let e = s
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let gradient = e.gradient.as_ref().unwrap();
    let directions = [
        vec![1.0; rho.len()],
        (0..rho.len())
            .map(|i| ((i % 3) as f64 - 1.0) * 0.5)
            .collect(),
    ];
    for direction in directions {
        let h = 2e-5;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = s
            .evaluate_robust_affine_load_space(
                &plus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                &mut control,
            )
            .unwrap();
        let b = s
            .evaluate_robust_affine_load_space(
                &minus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                &mut control,
            )
            .unwrap();
        assert_eq!(a.oracle.branch, QuadraticBallBranch::Regular);
        assert_eq!(b.oracle.branch, QuadraticBallBranch::Regular);
        close(
            (a.maximum - b.maximum) / (2.0 * h),
            dot(gradient, &direction),
            5e-4,
        );
        close(
            (a.physical.dilated_volume_fraction - b.physical.dilated_volume_fraction) / (2.0 * h),
            dot(&e.physical.dilated_volume_gradient, &direction),
            5e-5,
        );
    }
}

#[test]
fn hard_case_retains_a_physical_witness_but_refuses_a_smooth_optimizer_gradient() {
    let (mut s, _, first, second) = fixture();
    let rho = raw(&s);
    let basis = [&first[..], &second[..]];
    let zero = vec![0.0; first.len()];
    let scales = s.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let regular = s
        .evaluate_robust_affine_load_space(
            &rho,
            &zero,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    // An authored dead load cancels gravity at this one design only. The body
    // load is still present and still depends on the physical design.
    let dead: Vec<_> = regular.offset_force.iter().map(|f| -f).collect();
    let hard = s
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    assert_eq!(hard.oracle.branch, QuadraticBallBranch::HardCase);
    assert!(!hard.oracle.gradient_admissible());
    assert!(hard.gradient.is_none());
    assert!(hard.maximum > 0.0);
    assert!(hard.offset_force.iter().all(|f| *f == 0.0));
    audit(&hard, &dead, &basis);
    assert!(
        AffineLoadSpaceDesign3::new(
            &mut s,
            &dead,
            &basis,
            &rho,
            policy(),
            AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit: 0.9 },
            AffineLoadSpaceDesignOptions3::default(),
            &mut control
        )
        .is_err()
    );
    assert_eq!(s.operator().scales(), scales);
}

#[test]
fn minimum_native_mass_under_the_affine_compliance_cap_rechecks_the_actual_endpoint() {
    let (mut s, dead, first, second) = fixture();
    let rho = vec![0.7; s.cells()];
    let basis = [&first[..], &second[..]];
    let scales = s.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let initial = s
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let cap = 2.0 * initial.worst_compliance;
    let cap_policy = AffineLoadSpace3::new(2, cap).unwrap();
    let options = AffineLoadSpaceDesignOptions3::default();
    let mut session = AffineLoadSpaceDesign3::new(
        &mut s,
        &dead,
        &basis,
        &rho,
        cap_policy,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
            volume_scale: 1.0,
        },
        options,
        &mut control,
    )
    .unwrap();
    assert!(session.gradient_check().passed());
    let initial_mass = session
        .accepted()
        .evaluation
        .physical
        .dilated_mass_kg
        .unwrap();
    let terminal = session.run(6).unwrap();
    assert!(session.history().len() > 1);
    let accepted = &session.accepted().evaluation;
    close(
        terminal.objective,
        accepted.physical.dilated_mass_kg.unwrap() / 2.0,
        1e-10,
    );
    close(terminal.constraint, accepted.maximum - 1.0, 1e-10);
    let best = session.best_feasible().unwrap();
    assert!(best.evaluation.physical.dilated_mass_kg.unwrap() < initial_mass);
    assert!(best.evaluation.maximum <= 1.0);
    assert!(best.evaluation.worst_compliance <= cap);
    let inventory: f64 = session
        .study()
        .operator()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&best.evaluation.physical.dilated_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(
        best.evaluation.physical.dilated_mass_kg.unwrap(),
        inventory,
        1e-12,
    );
    assert!(best.evaluation.physical.dilated_mass_gradient.is_some());
    let best_rho = best.rho.clone();
    let best_evaluation = best.evaluation.clone();
    let history_len = session.history().len();
    let optimizer_work = session.optimizer_work();
    let rechecked = session.recheck_best().unwrap().unwrap();
    assert_eq!(session.history().len(), history_len);
    assert_eq!(session.optimizer_work(), optimizer_work);
    assert_eq!(rechecked.rho, best_rho);
    close(rechecked.evaluation.maximum, best_evaluation.maximum, 1e-9);
    close_vector(
        &rechecked.evaluation.worst_displacement,
        &best_evaluation.worst_displacement,
    );
    audit(&rechecked.evaluation, &dead, &basis);
    drop(session);
    assert_eq!(s.operator().scales(), scales);
}

fn adaptive_operator(level: u8) -> AdaptiveSolveSpace3 {
    let tree = Octree3::uniform(level, 4, 4096).unwrap();
    adaptive_on_tree(&tree)
}
fn adaptive_on_tree(tree: &Octree3) -> AdaptiveSolveSpace3 {
    adaptive_on_tree_controlled(tree, &mut poll).unwrap()
}
fn adaptive_on_tree_controlled(
    tree: &Octree3,
    checkpoint: &mut dyn FnMut() -> ControlFlow<()>,
) -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
    let mut callback = |_| checkpoint();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            root_relative_tolerance: 1e-12,
            ..Default::default()
        },
        &mut callback,
    )
    .map_err(ElasticityError3::from)?;
    let native = AdaptiveElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        tree,
        &Top,
        ElasticMaterialField3::ByCell(&native_material),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut q,
    )?;
    Ok(AdaptiveSolveSpace3::jacobi(native, 100_000_000))
}
fn adaptive_material_study(tree: &Octree3) -> CutDensityStudy3<AdaptiveSolveSpace3> {
    let study = CutDensityStudy3::new(
        adaptive_on_tree(tree),
        0.15,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.03,
            eta: 0.5,
        },
    )
    .with_material_acceleration(acceleration, 100_000, poll)
    .unwrap();
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = PhysicalRegion3::Void;
    study.with_physical_regions(regions).unwrap()
}
fn design_first(p: [f64; 3]) -> [f64; 3] {
    [0.0, -0.001 * (1.0 + p[1]), 0.0]
}
fn design_second(p: [f64; 3]) -> [f64; 3] {
    [0.0, -0.0006 * (1.0 + p[1]), -0.0008 * (1.0 + p[0])]
}
fn design_laws() -> [ReferenceLoad3<'static>; 2] {
    [
        ReferenceLoad3::body(&design_first),
        ReferenceLoad3::body(&design_second),
    ]
}
fn mass_goal() -> AffineLoadSpaceDesignGoal3 {
    AffineLoadSpaceDesignGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
        volume_scale: 1.0,
    }
}
fn refinement_dead(_: [f64; 3]) -> [f64; 3] {
    [0.0002, -0.0001, -0.0003]
}
fn refinement_first(_: [f64; 3]) -> [f64; 3] {
    [0.0, -0.001, 0.0]
}
fn refinement_second(p: [f64; 3]) -> [f64; 3] {
    // P2 has zero moments against the coarse Q1 functions on each y-half.
    // Its genuinely new fine-grid loading must participate in remaximization.
    let t = if p[1] < 0.5 {
        4.0 * p[1] - 1.0
    } else {
        4.0 * p[1] - 3.0
    };
    [0.0, 0.0, 0.0002 + 0.08 * (3.0 * t * t - 1.0) / 2.0]
}
fn combine(offset: &[f64], columns: &[Vec<f64>], coefficients: &[f64]) -> Vec<f64> {
    (0..offset.len())
        .map(|d| {
            offset[d]
                + columns
                    .iter()
                    .zip(coefficients)
                    .map(|(column, a)| column[d] * a)
                    .sum::<f64>()
        })
        .collect()
}

#[test]
fn affine_refinement_remaximizes_native_loads_and_decomposes_regular_and_hard_witnesses() {
    for self_weight in [true, false] {
        let mut source = CutDensityStudy3::new(
            adaptive_operator(1),
            0.15,
            SimpParams {
                beta: 2.0,
                penal: 3.0,
                e_min: 0.03,
                eta: 0.5,
            },
        )
        .with_material_acceleration(
            move |p| {
                if self_weight {
                    acceleration(p)
                } else {
                    [0.0; 3]
                }
            },
            100_000,
            poll,
        )
        .unwrap();
        let mut regions = vec![PhysicalRegion3::Design; source.cells()];
        regions[0] = PhysicalRegion3::Solid;
        let last = regions.len() - 1;
        regions[last] = PhysicalRegion3::Void;
        source = source.with_physical_regions(regions).unwrap();
        let rho: Vec<_> = (0..source.cells()).map(|i| 0.5 + 0.03 * i as f64).collect();
        let dead = if self_weight {
            ReferenceLoad3::body(&refinement_dead)
        } else {
            ReferenceLoad3::default()
        };
        let laws = [
            ReferenceLoad3::body(&refinement_first),
            ReferenceLoad3::body(&refinement_second),
        ];
        let scales = source.operator().elasticity().scales().to_vec();
        let fingerprint = source.material_body_load_fingerprint();
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let before = control.work();
        let report = source
            .estimate_affine_load_space_enrichment(
                adaptive_operator(2),
                &rho,
                dead,
                &laws,
                RobustOptions3::default(),
                &policy(),
                AffineLoadSpaceRefinementOptions3::default(),
                &mut control,
            )
            .unwrap();
        // One forward filter and at most m+1 primal solves per grid. Frozen
        // directions use the retained basis fields, with no extra adjoints.
        assert!(report.work.linear_solves - before.linear_solves <= 1 + 2 * (laws.len() + 1));
        assert_eq!(source.operator().elasticity().scales(), scales);
        assert_eq!(source.material_body_load_fingerprint(), fingerprint);
        if !self_weight {
            assert_eq!(report.coarse.oracle.branch, QuadraticBallBranch::HardCase);
            assert_eq!(report.fine.oracle.branch, QuadraticBallBranch::HardCase);
            assert!(!report.coarse.oracle.gradient_admissible());
            assert!(!report.fine.oracle.gradient_admissible());
        }

        let mut fine = adaptive_operator(2);
        let parents = AdaptiveTransfer3::new(
            source.operator().elasticity(),
            fine.elasticity(),
            2_000_000,
            poll,
        )
        .unwrap()
        .parents()
        .to_vec();
        assert_eq!(
            report.fine_eroded_rho,
            parents
                .iter()
                .map(|&i| report.coarse_eroded_rho[i])
                .collect::<Vec<_>>()
        );
        assert_eq!(
            report.fine_eroded_scales,
            parents
                .iter()
                .map(|&i| report.coarse_eroded_scales[i])
                .collect::<Vec<_>>()
        );
        assert!(report.fine_eroded_rho.contains(&0.0));
        assert!(report.fine_eroded_rho.contains(&1.0));
        fine.set_scales(&report.fine_eroded_scales).unwrap();
        let native = fine.elasticity();
        let body = native
            .prepare_material_body_load(
                &|p| {
                    if self_weight {
                        acceleration(p)
                    } else {
                        [0.0; 3]
                    }
                },
                100_000,
                poll,
            )
            .unwrap();
        let body_force = body.load(&report.fine_eroded_rho, poll).unwrap();
        let fixed = native.reference_load(dead, poll).unwrap();
        let fine_offset: Vec<_> = fixed.iter().zip(body_force).map(|(f, b)| f + b).collect();
        close_vector(&fine_offset, &report.fine.offset_force);
        let fine_basis: Vec<_> = laws
            .iter()
            .map(|&law| native.reference_load(law, poll).unwrap())
            .collect();
        let fine_force = combine(&fine_offset, &fine_basis, &report.fine.worst_coefficients);
        close_vector(&fine_force, &report.fine.worst_force);
        assert!(
            native
                .field_residual(&report.fine.worst_displacement, &fine_force, poll)
                .unwrap()
                < 1e-8
        );
        if self_weight {
            let omitted = combine(&fixed, &fine_basis, &report.fine.worst_coefficients);
            assert!(
                native
                    .field_residual(&report.fine.worst_displacement, &omitted, poll)
                    .unwrap()
                    > 1e-6
            );
        }

        // A=Jc(ac), B=Jf(ac), C=Jc(af), D=Jf(af), all in physical units.
        let a = dot(
            &report.coarse.worst_force,
            &report.coarse.worst_displacement,
        );
        let old_fine_force = combine(&fine_offset, &fine_basis, &report.coarse.worst_coefficients);
        let old_fine_u = combine(
            &report.fine.offset_displacement,
            &report.fine.basis_displacements,
            &report.coarse.worst_coefficients,
        );
        assert!(
            native
                .field_residual(&old_fine_u, &old_fine_force, poll)
                .unwrap()
                < 1e-8
        );
        let b = dot(&old_fine_force, &old_fine_u);
        let coarse_basis: Vec<_> = laws
            .iter()
            .map(|&law| {
                source
                    .operator()
                    .elasticity()
                    .reference_load(law, poll)
                    .unwrap()
            })
            .collect();
        let new_coarse_force = combine(
            &report.coarse.offset_force,
            &coarse_basis,
            &report.fine.worst_coefficients,
        );
        let new_coarse_u = combine(
            &report.coarse.offset_displacement,
            &report.coarse.basis_displacements,
            &report.fine.worst_coefficients,
        );
        let c = dot(&new_coarse_force, &new_coarse_u);
        let d = dot(&fine_force, &report.fine.worst_displacement);
        let reference = policy().reference();
        close(a, report.coarse.worst_compliance, 1e-9);
        close(d, report.fine.worst_compliance, 1e-9);
        // This deliberately unresolved load is more adverse on the fine grid
        // than the old direction. Checking only that old direction cannot pass.
        assert!((d - b) / reference > 1e-8);
        assert!(report.fine.maximum > b / reference + 1e-8);
        close(report.fine_switch_correction, (d - b) / reference, 1e-8);
        close(report.coarse_switch_correction, (a - c) / reference, 1e-8);
        for (witness, coarse, fine) in
            [(&report.coarse_witness, a, b), (&report.fine_witness, c, d)]
        {
            close(witness.coarse_value, coarse, 1e-8);
            close(witness.fine_value, fine, 1e-8);
            close(witness.correction(), fine - coarse, 1e-8);
            close(
                witness.cells.values().map(|cell| cell.signed_dwr).sum(),
                witness.dwr,
                1e-8,
            );
            assert!(witness.identity_relative_defect < 1e-8);
        }
        let symmetric = 0.5
            * (report.coarse_witness.correction() + report.fine_witness.correction())
            / reference
            + 0.5 * (report.fine_switch_correction - report.coarse_switch_correction);
        close(report.correction(), symmetric, 1e-10);
        close(report.correction(), (d - a) / reference, 1e-8);
        close(report.normalized_witness_change, (d - a) / reference, 1e-10);
        close(
            report.normalized_maximum_change,
            report.fine.maximum - report.coarse.maximum,
            1e-10,
        );
        close(report.normalized_maximum_change, report.correction(), 1e-8);
        assert!(
            report.normalized_maximum_change.abs() <= report.normalized_model_change * (1.0 + 1e-8)
        );
        assert!(report.identity_relative_defect < 1e-8);
        assert!(report.marking_mass.values().sum::<f64>() > 0.0);
        for (&leaf, &mass) in &report.marking_mass {
            assert!(source.operator().elasticity().leaves().contains(&leaf));
            let expected = 0.5
                * (report.coarse_witness.cells[&leaf].marking_mass
                    + report.fine_witness.cells[&leaf].marking_mass)
                / reference;
            close(mass, expected, 1e-10);
        }
    }
}

#[test]
fn explicit_affine_refit_reoptimizes_native_mass_on_its_own_grid_and_rechecks_the_full_cap() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let source = adaptive_material_study(&tree);
    let rho = vec![0.7; source.cells()];
    let fine_tree = tree
        .refined(&[source.operator().elasticity().leaves()[2]], poll)
        .unwrap();
    let candidate = adaptive_on_tree(&fine_tree);
    let parents = AdaptiveTransfer3::new(
        source.operator().elasticity(),
        candidate.elasticity(),
        2_000_000,
        poll,
    )
    .unwrap()
    .parents()
    .to_vec();
    let inherited: Vec<_> = parents.iter().map(|&p| rho[p]).collect();
    let labels: Vec<_> = parents
        .iter()
        .map(|&p| source.physical_regions().unwrap()[p])
        .collect();
    let mut direct = adaptive_material_study(&fine_tree)
        .with_physical_regions(labels.clone())
        .unwrap();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let laws = design_laws();
    let direct_dead = direct
        .operator()
        .elasticity()
        .reference_load(dead, poll)
        .unwrap();
    let direct_forces: Vec<_> = laws
        .iter()
        .map(|&law| {
            direct
                .operator()
                .elasticity()
                .reference_load(law, poll)
                .unwrap()
        })
        .collect();
    let basis: Vec<_> = direct_forces.iter().map(Vec::as_slice).collect();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    // Author the fixed cap from an independently refiltered fine baseline.
    // No assumption that increasing material restores self-weight feasibility.
    let baseline = direct
        .evaluate_robust_affine_load_space(
            &inherited,
            &direct_dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let cap = AffineLoadSpace3::new(laws.len(), 1.8 * baseline.worst_compliance).unwrap();
    let initial = source
        .solve_reference_affine_load_space(
            &rho,
            dead,
            &laws,
            &cap,
            mass_goal(),
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 0,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    let old_scales = initial.study.operator().elasticity().scales().to_vec();
    let old_fingerprint = initial.study.material_body_load_fingerprint();
    let old_point = initial.accepted.as_ref().unwrap().clone();
    let result = initial
        .study
        .refit_reference_affine_load_space(
            candidate,
            &rho,
            dead,
            &laws,
            &cap,
            mass_goal(),
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 6,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    assert!(result.completed_and_rechecked(), "{:?}", result.outcome);
    assert_eq!(result.policy, cap);
    assert_eq!(result.goal, mass_goal());
    assert_eq!(result.inherited_rho, inherited);
    assert_eq!(result.study.physical_regions().unwrap(), labels);
    assert!(result.study.cells() > initial.study.cells());
    assert!(result.study.can_refine_material_body_load());
    assert_ne!(
        result.study.material_body_load_fingerprint(),
        old_fingerprint
    );
    assert!(result.gradient_check.as_ref().unwrap().passed());
    let first = result.initial.as_ref().unwrap();
    let accepted = result.accepted.as_ref().unwrap();
    let best = result.best.as_ref().unwrap();
    let checked = result.recheck.as_ref().unwrap().as_ref().unwrap();
    assert_eq!(first.rho, inherited);
    assert_eq!(result.history[0].iteration, 0);
    assert!(result.history.len() > 1);
    close(
        first.evaluation.worst_compliance,
        baseline.worst_compliance,
        1e-9,
    );
    close(
        result.history.last().unwrap().objective,
        accepted.evaluation.physical.dilated_mass_kg.unwrap() / 2.0,
        1e-10,
    );
    assert!(
        best.evaluation.physical.dilated_mass_kg.unwrap()
            < first.evaluation.physical.dilated_mass_kg.unwrap()
    );
    assert_eq!(checked.rho, best.rho);
    assert!(checked.evaluation.worst_compliance <= cap.reference());
    assert!(checked.evaluation.maximum <= 1.0);
    close_vector(
        &checked.evaluation.worst_displacement,
        &best.evaluation.worst_displacement,
    );

    let fresh = direct
        .evaluate_robust_affine_load_space(
            &checked.rho,
            &direct_dead,
            &basis,
            RobustOptions3::default(),
            &cap,
            &mut control,
        )
        .unwrap();
    close_vector(&fresh.gram, &checked.evaluation.gram);
    close_vector(&fresh.linear, &checked.evaluation.linear);
    close_vector(&fresh.worst_force, &checked.evaluation.worst_force);
    close_vector(
        &fresh.worst_displacement,
        &checked.evaluation.worst_displacement,
    );
    close(fresh.maximum, checked.evaluation.maximum, 1e-9);
    let kg: f64 = direct
        .operator()
        .elasticity()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&fresh.physical.dilated_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(
        kg,
        checked.evaluation.physical.dilated_mass_kg.unwrap(),
        1e-12,
    );
    let mut native = adaptive_on_tree(&fine_tree);
    native.set_scales(&fresh.physical.eroded_scales).unwrap();
    let body = native
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let self_weight = body.load(&fresh.physical.eroded_rho, poll).unwrap();
    let offset: Vec<_> = direct_dead
        .iter()
        .zip(self_weight)
        .map(|(f, b)| f + b)
        .collect();
    let force = combine(&offset, &direct_forces, &fresh.worst_coefficients);
    assert!(
        native
            .elasticity()
            .field_residual(&checked.evaluation.worst_displacement, &force, poll)
            .unwrap()
            < 1e-8
    );
    assert_eq!(initial.study.operator().elasticity().scales(), old_scales);
    assert_eq!(
        initial.study.material_body_load_fingerprint(),
        old_fingerprint
    );
    assert_eq!(initial.accepted.as_ref().unwrap().rho, old_point.rho);
    assert_eq!(
        initial
            .accepted
            .as_ref()
            .unwrap()
            .evaluation
            .worst_displacement,
        old_point.evaluation.worst_displacement
    );
}

#[test]
fn affine_mass_campaign_reoptimizes_promotes_and_probes_its_actual_new_design() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let old_cells = study.cells();
    let old_fingerprint = study.material_body_load_fingerprint();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let laws = design_laws();
    // A deliberately loose physical cap isolates material optimization from
    // feasibility restoration. Each grid still has to pass its fresh cap check.
    let cap = AffineLoadSpace3::new(laws.len(), 1.0).unwrap();
    // The unused volume normalization must not restrict a kilogram objective.
    let goal = AffineLoadSpaceDesignGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
        volume_scale: 0.0,
    };
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let initial = study
        .solve_reference_affine_load_space(
            &rho,
            dead,
            &laws,
            &cap,
            goal,
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 1,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    let builds = RefCell::new(Vec::new());
    let mut campaign = AffineLoadSpaceCampaign3::new(
        initial,
        tree,
        dead,
        &laws,
        AffineLoadSpaceCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 1,
            change_tolerance: 0.0,
            volume_tolerance: 1e-8,
            objective_tolerance: 1e-8,
            marking: AffineLoadSpaceMarking3::FrozenWitness {
                theta: 0.5,
                max_marks: 1,
            },
            ..Default::default()
        },
        &mut control,
        |tree, kind, checkpoint| {
            builds.borrow_mut().push(kind);
            adaptive_on_tree_controlled(tree, checkpoint)
        },
    )
    .unwrap();
    assert_eq!(
        campaign.run(0).unwrap(),
        AffineLoadSpaceCampaignStop3::Paused
    );
    assert_eq!(campaign.grids().len(), 1);
    assert_eq!(campaign.probes().len(), 1);
    assert_eq!(campaign.candidate_attempts(), 0);
    assert_eq!(builds.borrow().as_slice(), &[AffineLoadSpaceBuild3::Probe]);
    let work = campaign.work();
    assert_eq!(
        campaign.run(0).unwrap(),
        AffineLoadSpaceCampaignStop3::Paused
    );
    assert_eq!(campaign.work(), work);
    assert_eq!(campaign.probes().len(), 1);
    assert_eq!(
        campaign.run(1).unwrap(),
        AffineLoadSpaceCampaignStop3::RefinementLimit
    );
    assert_eq!(
        builds.borrow().as_slice(),
        &[
            AffineLoadSpaceBuild3::Probe,
            AffineLoadSpaceBuild3::Refit,
            AffineLoadSpaceBuild3::Probe,
        ]
    );
    assert_eq!(campaign.grids().len(), 2);
    assert_eq!(campaign.probes().len(), 2);
    assert_eq!(campaign.probes()[1].source_grid, 1);
    assert_eq!(campaign.candidate_attempts(), 1);
    assert!(campaign.pending().is_none());
    assert!(campaign.rejected().is_none());
    let current = campaign.current();
    assert!(current.completed_and_rechecked());
    assert_eq!(current.policy, cap);
    assert_eq!(current.goal, goal);
    assert!(current.study.cells() > old_cells);
    assert!(current.study.can_refine_material_body_load());
    assert_ne!(
        current.study.material_body_load_fingerprint(),
        old_fingerprint
    );
    assert!(current.gradient_check.as_ref().unwrap().passed());
    assert!(current.history.len() > 1);
    let start = current.initial.as_ref().unwrap();
    let checked = current.recheck.as_ref().unwrap().as_ref().unwrap();
    assert_eq!(checked.rho, current.best.as_ref().unwrap().rho);
    assert!(
        checked.evaluation.physical.dilated_mass_kg.unwrap()
            < start.evaluation.physical.dilated_mass_kg.unwrap()
    );
    assert!(checked.evaluation.maximum <= 1.0);
    assert!(checked.evaluation.worst_compliance <= cap.reference());
    let probe = campaign.current_probe().unwrap();
    assert_eq!(probe.rho, checked.rho);
    close(probe.coarse.maximum, checked.evaluation.maximum, 1e-9);
    let inventory = probe.inventory.as_ref().unwrap();
    assert_eq!(
        inventory.objective,
        MinimumMaterialObjective3::MassKg { scale_kg: 2.0 }
    );
    assert_eq!(inventory.volume_scale, 0.0);
    close(
        inventory.coarse_mass_kg.unwrap(),
        checked.evaluation.physical.dilated_mass_kg.unwrap(),
        1e-12,
    );
    close(
        inventory.coarse_objective,
        inventory.coarse_mass_kg.unwrap() / 2.0,
        1e-12,
    );
    close(
        inventory.coarse_volume_fraction,
        checked.evaluation.physical.dilated_volume_fraction,
        1e-12,
    );

    // Reconstruct B and the authored loads on the installed native owner.
    // This audits the promoted field without another optimization or PDE solve.
    let mut native = adaptive_on_tree(campaign.tree());
    native
        .set_scales(&checked.evaluation.physical.eroded_scales)
        .unwrap();
    let body = native
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let weight = body
        .load(&checked.evaluation.physical.eroded_rho, poll)
        .unwrap();
    let offset: Vec<_> = native
        .elasticity()
        .reference_load(dead, poll)
        .unwrap()
        .iter()
        .zip(weight)
        .map(|(f, b)| f + b)
        .collect();
    let columns: Vec<_> = laws
        .iter()
        .map(|&law| native.elasticity().reference_load(law, poll).unwrap())
        .collect();
    let force = combine(&offset, &columns, &checked.evaluation.worst_coefficients);
    close_vector(&force, &checked.evaluation.worst_force);
    assert!(
        native
            .elasticity()
            .field_residual(&checked.evaluation.worst_displacement, &force, poll)
            .unwrap()
            < 1e-8
    );
    let kg: f64 = native
        .elasticity()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&checked.evaluation.physical.dilated_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(kg, inventory.coarse_mass_kg.unwrap(), 1e-12);
}

fn switching_dead(p: [f64; 3]) -> [f64; 3] {
    let ordinary = refinement_dead(p);
    let unresolved = refinement_second(p);
    // The fixed offset has a component along the new fine-grid mode. Both
    // endpoints remain regular; rejection must come from the physical cap.
    std::array::from_fn(|i| ordinary[i] + 2.0 * unresolved[i])
}

fn reference_baseline(
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    rho: &[f64],
    dead: ReferenceLoad3<'_>,
    laws: &[ReferenceLoad3<'_>],
    control: &mut SolveControl<'_>,
) -> AffineLoadSpaceEvaluation3 {
    let offset = study
        .operator()
        .elasticity()
        .reference_load(dead, poll)
        .unwrap();
    let columns: Vec<_> = laws
        .iter()
        .map(|&law| {
            study
                .operator()
                .elasticity()
                .reference_load(law, poll)
                .unwrap()
        })
        .collect();
    let basis: Vec<_> = columns.iter().map(Vec::as_slice).collect();
    study
        .evaluate_robust_affine_load_space(
            rho,
            &offset,
            &basis,
            RobustOptions3::default(),
            &policy(),
            control,
        )
        .unwrap()
}

#[test]
fn affine_mass_campaign_retains_the_source_when_a_regular_fine_candidate_exceeds_the_cap() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut source = adaptive_material_study(&tree);
    let rho = vec![0.7; source.cells()];
    let fine_tree = tree
        .refined(source.operator().elasticity().leaves(), poll)
        .unwrap();
    let fine = adaptive_on_tree(&fine_tree);
    let parents = AdaptiveTransfer3::new(
        source.operator().elasticity(),
        fine.elasticity(),
        2_000_000,
        poll,
    )
    .unwrap()
    .parents()
    .to_vec();
    let inherited: Vec<_> = parents.iter().map(|&p| rho[p]).collect();
    let labels: Vec<_> = parents
        .iter()
        .map(|&p| source.physical_regions().unwrap()[p])
        .collect();
    let mut direct = adaptive_material_study(&fine_tree)
        .with_physical_regions(labels)
        .unwrap();
    let dead = ReferenceLoad3::body(&switching_dead);
    let laws = [
        ReferenceLoad3::body(&refinement_first),
        ReferenceLoad3::body(&refinement_second),
    ];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    // Evaluate the refiltered raw candidate, not the inherited-physical probe:
    // it is this distinct full affine model that a candidate must reoptimize.
    let coarse = reference_baseline(&mut source, &rho, dead, &laws, &mut control);
    let fine = reference_baseline(&mut direct, &inherited, dead, &laws, &mut control);
    assert!(coarse.oracle.gradient_admissible());
    assert!(fine.oracle.gradient_admissible());
    assert!(fine.worst_compliance > coarse.worst_compliance);
    let cap = AffineLoadSpace3::new(
        laws.len(),
        0.5 * (coarse.worst_compliance + fine.worst_compliance),
    )
    .unwrap();
    let initial = source
        .solve_reference_affine_load_space(
            &rho,
            dead,
            &laws,
            &cap,
            mass_goal(),
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 0,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    let source_point = initial.recheck.as_ref().unwrap().as_ref().unwrap().clone();
    let source_scales = initial.study.operator().elasticity().scales().to_vec();
    let source_fingerprint = initial.study.material_body_load_fingerprint();
    let source_regions = initial.study.physical_regions().unwrap().to_vec();
    let builds = RefCell::new(Vec::new());
    let mut campaign = AffineLoadSpaceCampaign3::new(
        initial,
        tree.clone(),
        dead,
        &laws,
        AffineLoadSpaceCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 1,
            change_tolerance: 0.0,
            volume_tolerance: 1e-8,
            objective_tolerance: 1e-8,
            marking: AffineLoadSpaceMarking3::UniformActive,
            ..Default::default()
        },
        &mut control,
        |tree, kind, checkpoint| {
            builds.borrow_mut().push(kind);
            adaptive_on_tree_controlled(tree, checkpoint)
        },
    )
    .unwrap();
    assert_eq!(
        campaign.run(1).unwrap(),
        AffineLoadSpaceCampaignStop3::CandidateStopped
    );
    assert_eq!(
        builds.borrow().as_slice(),
        &[AffineLoadSpaceBuild3::Probe, AffineLoadSpaceBuild3::Refit]
    );
    assert_eq!(campaign.grids().len(), 1);
    assert_eq!(campaign.probes().len(), 1);
    assert_eq!(campaign.candidate_attempts(), 1);
    assert!(campaign.pending().is_none());
    let (rejected_tree, rejected) = campaign.rejected().unwrap();
    assert_eq!(rejected_tree.leaves(), fine_tree.leaves());
    assert!(!rejected.completed_and_rechecked());
    assert!(rejected.outcome.is_ok());
    assert!(rejected.gradient_check.as_ref().unwrap().passed());
    assert_eq!(rejected.history.len(), 1);
    assert_eq!(rejected.inherited_rho, inherited);
    let rejected_point = rejected.initial.as_ref().unwrap();
    assert!(rejected_point.evaluation.oracle.gradient_admissible());
    close(
        rejected_point.evaluation.worst_compliance,
        fine.worst_compliance,
        1e-9,
    );
    assert!(rejected_point.evaluation.worst_compliance > cap.reference());
    assert!(rejected_point.evaluation.maximum > 1.0);
    assert!(rejected.best.is_none());
    assert!(rejected.recheck.as_ref().unwrap().is_err());
    let current = campaign.current();
    assert!(current.completed_and_rechecked());
    assert_eq!(campaign.tree().leaves(), tree.leaves());
    assert_eq!(
        current.study.operator().elasticity().scales(),
        source_scales
    );
    assert_eq!(
        current.study.material_body_load_fingerprint(),
        source_fingerprint
    );
    assert_eq!(current.study.physical_regions().unwrap(), source_regions);
    let still_checked = current.recheck.as_ref().unwrap().as_ref().unwrap();
    assert_eq!(still_checked.rho, source_point.rho);
    assert_eq!(
        still_checked.evaluation.worst_force,
        source_point.evaluation.worst_force
    );
    assert_eq!(
        still_checked.evaluation.worst_displacement,
        source_point.evaluation.worst_displacement
    );
    assert!(still_checked.evaluation.maximum <= 1.0);
    assert!(still_checked.evaluation.worst_compliance <= cap.reference());
}

// Fix one genuine linear displacement observation whose affine response vector
// is orthogonal to the compliance adversary at this baseline. The observation
// stays fixed through every subsequent raw-density perturbation.
fn transverse_observation(e: &AffineLoadSpaceEvaluation3) -> Vec<f64> {
    let u = &e.basis_displacements;
    let g00 = dot(&u[0], &u[0]);
    let g01 = dot(&u[0], &u[1]);
    let g11 = dot(&u[1], &u[1]);
    let determinant = g00 * g11 - g01 * g01;
    assert!(determinant > 1e-8 * g00 * g11);
    let target = [
        -1e-4 * e.worst_coefficients[1],
        1e-4 * e.worst_coefficients[0],
    ];
    let x = (g11 * target[0] - g01 * target[1]) / determinant;
    let y = (g00 * target[1] - g01 * target[0]) / determinant;
    u[0].iter().zip(&u[1]).map(|(a, b)| x * a + y * b).collect()
}

#[test]
fn affine_response_ball_has_its_own_adversaries_and_factored_native_body_adjoints() {
    let (mut study, dead, first, second) = fixture();
    let rho = raw(&study);
    let basis = [&first[..], &second[..]];
    let scales = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = study
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let ordinary_solves = control.work().linear_solves;
    let q = transverse_observation(&baseline);
    let center = dot(&q, &baseline.offset_displacement);
    let response_basis: Vec<_> = baseline
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let radius = dot(&response_basis, &response_basis).sqrt();
    let limits = [
        ResponseLimit3 {
            response: BoundedResponse3::Displacement(&q),
            sense: ResponseLimitSense3::AtMost,
            bound: center + 1.1 * radius,
            scale: radius,
        },
        ResponseLimit3 {
            response: BoundedResponse3::Displacement(&q),
            sense: ResponseLimitSense3::AtLeast,
            bound: center - 2.0 * radius,
            scale: radius,
        },
    ];
    let options = AffineResponseOptions3 {
        sharpness: 20.0,
        ..Default::default()
    };
    let pair = study
        .evaluate_robust_affine_response_limits(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            AffineResponseLimits3 {
                limits: &limits,
                options,
            },
            &mut control,
        )
        .unwrap();
    let response = pair.responses.as_ref().unwrap();
    assert_eq!(response.adjoints.len(), 2);
    assert!(response.physically_feasible());
    for (i, witness) in response.witnesses.iter().enumerate() {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        close(witness.baseline, center, 1e-10);
        close_vector(&witness.basis_responses, &response_basis);
        close(witness.radius, radius, 1e-10);
        close(witness.worst_value, center + sign * radius, 1e-10);
        close(
            witness.residual,
            sign * (witness.worst_value - witness.bound) / witness.scale,
            1e-10,
        );
        close(
            dot(&witness.coefficients, &witness.coefficients),
            1.0,
            1e-10,
        );
        assert!(dot(&witness.coefficients, &pair.worst_coefficients).abs() < 1e-7);
        let field = combine(
            &pair.offset_displacement,
            &pair.basis_displacements,
            &witness.coefficients,
        );
        close(dot(&q, &field), witness.worst_value, 1e-10);
        for j in 0..360 {
            let angle = std::f64::consts::TAU * j as f64 / 360.0;
            let value = center + response_basis[0] * angle.cos() + response_basis[1] * angle.sin();
            assert!(sign * (value - witness.worst_value) <= 1e-10 * radius);
        }
    }
    // Four original inequalities share only three affine-field adjoints. Exact
    // duplication adds log(2)/beta to this unnormalized upper aggregate, while
    // preserving its derivative. It must not duplicate the primal family.
    let duplicated = [limits[0], limits[1], limits[0], limits[1]];
    let family = AffineResponseLimits3 {
        limits: &duplicated,
        options,
    };
    let before = control.work().linear_solves;
    let factored = study
        .evaluate_robust_affine_response_limits(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            family,
            &mut control,
        )
        .unwrap();
    assert!(control.work().linear_solves - before <= ordinary_solves + 4);
    let full = factored.responses.as_ref().unwrap();
    assert_eq!(full.adjoints.len(), 3);
    close(
        full.constraint - response.constraint,
        2.0_f64.ln() / options.sharpness,
        1e-10,
    );
    let gradient = full.gradient.as_ref().unwrap();
    for (&a, &b) in gradient.iter().zip(response.gradient.as_ref().unwrap()) {
        close(a, b, 3e-6);
    }
    // Freezing exactly this baseline's body force preserves the current fields
    // but removes B^T z from their derivative. The two gradients must differ.
    let mut frozen = CutDensityStudy3::new(operator(), 0.15, study.params())
        .with_physical_regions(study.physical_regions().unwrap().to_vec())
        .unwrap();
    let frozen = frozen
        .evaluate_robust_affine_response_limits(
            &rho,
            &baseline.offset_force,
            &basis,
            RobustOptions3::default(),
            &policy(),
            family,
            &mut control,
        )
        .unwrap();
    let frozen = frozen.responses.as_ref().unwrap();
    close(full.constraint, frozen.constraint, 1e-9);
    let difference = gradient
        .iter()
        .zip(frozen.gradient.as_ref().unwrap())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let gradient_scale = gradient.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
    assert!(difference > 1e-4 * gradient_scale);
    for direction in [
        vec![1.0; rho.len()],
        (0..rho.len()).map(|i| (i % 3) as f64 * 0.5 - 0.5).collect(),
    ] {
        let h = 2e-5;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = study
            .evaluate_robust_affine_response_limits(
                &plus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                family,
                &mut control,
            )
            .unwrap();
        let b = study
            .evaluate_robust_affine_response_limits(
                &minus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                family,
                &mut control,
            )
            .unwrap();
        close(
            (a.responses.unwrap().constraint - b.responses.unwrap().constraint) / (2.0 * h),
            dot(gradient, &direction),
            5e-4,
        );
    }
    assert_eq!(study.operator().scales(), scales);
}

fn embedded_response_operator() -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            root_relative_tolerance: 1e-12,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    CutElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Top,
        ElasticMaterialField3::ByCell(&native_material),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut quadrature,
    )
    .unwrap()
}
fn weighted_support_mode(p: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [0.0, 1.0 + p[0], 0.0]
}

#[test]
fn affine_reaction_limits_include_direct_material_and_body_terms_and_roll_back_stops() {
    let (mut study, dead, first, second) = fixture_on_operator(embedded_response_operator());
    let rho = raw(&study);
    let basis = [&first[..], &second[..]];
    let scales = study.operator().scales().to_vec();
    let fingerprint = study.material_body_load_fingerprint();
    let limits = [
        ResponseLimit3 {
            response: BoundedResponse3::Reaction(&weighted_support_mode),
            sense: ResponseLimitSense3::AtMost,
            bound: 0.01,
            scale: 0.01,
        },
        ResponseLimit3 {
            response: BoundedResponse3::Reaction(&weighted_support_mode),
            sense: ResponseLimitSense3::AtLeast,
            bound: -0.02,
            scale: 0.01,
        },
    ];
    let family = AffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let evaluated = study
        .evaluate_robust_affine_response_limits(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            family,
            &mut control,
        )
        .unwrap();
    let response = evaluated.responses.as_ref().unwrap();
    let mut native = embedded_response_operator();
    native
        .set_scales(&evaluated.physical.eroded_scales)
        .unwrap();
    for witness in &response.witnesses {
        let field = combine(
            &evaluated.offset_displacement,
            &evaluated.basis_displacements,
            &witness.coefficients,
        );
        let direct = native
            .embedded_reaction(&field, None, &weighted_support_mode, poll)
            .unwrap();
        close(direct.value, witness.worst_value, 1e-9);
        assert!(direct.scale_gradient.iter().any(|v| v.abs() > 1e-10));
    }
    for direction in [
        vec![1.0; rho.len()],
        (0..rho.len()).map(|i| (i % 3) as f64 * 0.5 - 0.5).collect(),
    ] {
        let h = 2e-5;
        let plus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let minus: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
        let a = study
            .evaluate_robust_affine_response_limits(
                &plus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                family,
                &mut control,
            )
            .unwrap();
        let b = study
            .evaluate_robust_affine_response_limits(
                &minus,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                family,
                &mut control,
            )
            .unwrap();
        close(
            (a.responses.unwrap().constraint - b.responses.unwrap().constraint) / (2.0 * h),
            dot(response.gradient.as_ref().unwrap(), &direction),
            5e-4,
        );
    }
    for stage in [
        "sdf3-affine-response-adjoint",
        "sdf3-affine-response-publish",
    ] {
        let mut stop = |p: SolveProgress| {
            if p.stage == stage && (stage.ends_with("publish") || p.solve_iterations > 0) {
                ControlFlow::Break(())
            } else {
                poll()
            }
        };
        let mut stopped = SolveControl::new(SolveBudget::default(), &mut stop);
        assert!(matches!(
            study.evaluate_robust_affine_response_limits(
                &rho,
                &dead,
                &basis,
                RobustOptions3::default(),
                &policy(),
                family,
                &mut stopped,
            ),
            Err(EvaluationStop::Cancelled)
        ));
        assert!(stopped.work().linear_iterations > 0);
        assert_eq!(study.operator().scales(), scales);
        assert_eq!(study.material_body_load_fingerprint(), fingerprint);
    }
}

#[test]
fn affine_mass_design_keeps_robust_response_and_compliance_caps_in_its_fresh_recheck() {
    let (mut study, dead, first, second) = fixture();
    let rho = vec![0.7; study.cells()];
    let basis = [&first[..], &second[..]];
    let q: Vec<_> = first.iter().map(|v| -v).collect();
    let scales = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = study
        .evaluate_robust_affine_load_space(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            &mut control,
        )
        .unwrap();
    let center = dot(&q, &baseline.offset_displacement);
    let projected: Vec<_> = baseline
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let lower = center - dot(&projected, &projected).sqrt();
    assert!(lower < 0.0);
    let cap = AffineLoadSpace3::new(2, 2.0 * baseline.worst_compliance).unwrap();
    let mut limits = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&q),
        sense: ResponseLimitSense3::AtLeast,
        bound: 0.5 * lower,
        scale: lower.abs(),
    }];
    {
        // Compliance already passes, but this original response bound does not.
        // It must not enter the feasible mass archive merely because it is stiff.
        let mut rejected = AffineLoadSpaceDesign3::new_with_response_limits(
            &mut study,
            &dead,
            &basis,
            &rho,
            cap.clone(),
            mass_goal(),
            AffineLoadSpaceDesignOptions3::default(),
            AffineResponseLimits3 {
                limits: &limits,
                options: AffineResponseOptions3::default(),
            },
            &mut control,
        )
        .unwrap();
        assert!(rejected.accepted().evaluation.maximum <= 1.0);
        assert!(
            rejected
                .accepted()
                .evaluation
                .responses
                .as_ref()
                .unwrap()
                .worst_residual
                > 0.0
        );
        assert!(rejected.best_feasible().is_none());
        assert!(rejected.recheck_best().unwrap().is_none());
    }
    limits[0].bound = 2.0 * lower;
    let mut design = AffineLoadSpaceDesign3::new_with_response_limits(
        &mut study,
        &dead,
        &basis,
        &rho,
        cap.clone(),
        mass_goal(),
        AffineLoadSpaceDesignOptions3::default(),
        AffineResponseLimits3 {
            limits: &limits,
            options: AffineResponseOptions3::default(),
        },
        &mut control,
    )
    .unwrap();
    assert!(design.gradient_check().passed());
    let start_mass = design
        .accepted()
        .evaluation
        .physical
        .dilated_mass_kg
        .unwrap();
    design.run(6).unwrap();
    assert!(design.history().len() > 1);
    let best = design.best_feasible().unwrap();
    assert!(best.evaluation.physical.dilated_mass_kg.unwrap() < start_mass);
    assert!(
        best.evaluation
            .responses
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(best.evaluation.maximum <= 1.0);
    let best_rho = best.rho.clone();
    let count = design.history().len();
    let work = design.optimizer_work();
    let checked = design.recheck_best().unwrap().unwrap();
    assert_eq!(checked.rho, best_rho);
    assert_eq!(design.history().len(), count);
    assert_eq!(design.optimizer_work(), work);
    let response = checked.evaluation.responses.as_ref().unwrap();
    assert!(response.physically_feasible());
    assert!(checked.evaluation.worst_compliance <= cap.reference());
    let actual_center = dot(&q, &checked.evaluation.offset_displacement);
    let actual_basis: Vec<_> = checked
        .evaluation
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let actual_lower = actual_center - dot(&actual_basis, &actual_basis).sqrt();
    close(actual_lower, response.witnesses[0].worst_value, 1e-9);
    assert!(actual_lower >= limits[0].bound);
    let inventory: f64 = design
        .study()
        .operator()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&checked.evaluation.physical.dilated_rho)
        .map(|(m, r)| m * r)
        .sum();
    close(
        checked.evaluation.physical.dilated_mass_kg.unwrap(),
        inventory,
        1e-12,
    );
    drop(design);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn a_zero_affine_response_radius_retains_values_without_inventing_a_smooth_design_gradient() {
    let (mut study, dead, first, second) = fixture();
    let rho = raw(&study);
    let basis = [&first[..], &second[..]];
    let zero = vec![0.0; study.operator().n()];
    let limits = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&zero),
        sense: ResponseLimitSense3::AtMost,
        bound: 1.0,
        scale: 1.0,
    }];
    let family = AffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    let scales = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let evaluated = study
        .evaluate_robust_affine_response_limits(
            &rho,
            &dead,
            &basis,
            RobustOptions3::default(),
            &policy(),
            family,
            &mut control,
        )
        .unwrap();
    assert!(evaluated.gradient.is_some());
    let response = evaluated.responses.as_ref().unwrap();
    assert_eq!(response.witnesses[0].radius, 0.0);
    assert_eq!(response.witnesses[0].worst_value, 0.0);
    assert_eq!(response.worst_residual, -1.0);
    assert!(response.physically_feasible());
    assert!(response.gradient.is_none());
    assert!(response.adjoints.is_empty());
    let cap = AffineLoadSpace3::new(2, 2.0 * evaluated.worst_compliance).unwrap();
    assert!(
        AffineLoadSpaceDesign3::new_with_response_limits(
            &mut study,
            &dead,
            &basis,
            &rho,
            cap,
            mass_goal(),
            AffineLoadSpaceDesignOptions3::default(),
            family,
            &mut control,
        )
        .is_err()
    );
    assert_eq!(study.operator().scales(), scales);
}
