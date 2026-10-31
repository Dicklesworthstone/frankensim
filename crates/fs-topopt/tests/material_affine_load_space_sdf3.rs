//! G1/G3: an affine uncertainty ball includes native material self-weight.
use std::ops::ControlFlow;

use fs_ascent::quadratic_ball::QuadraticBallBranch;
use fs_blake3::ContentHash;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::sdf3::robust::RobustOptions3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::design::{
    AffineLoadSpaceDesign3, AffineLoadSpaceDesignGoal3, AffineLoadSpaceDesignOptions3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::{
    AffineLoadSpace3, AffineLoadSpaceEvaluation3,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::{SimpParams, SolveBudget, SolveControl};

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
    let field =
        |cell: HexCell| {
            let id = if cell.lo()[0] < 0.5 { 2 } else { 8 };
            let tensor =
                std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
            Ok(ElasticAcousticState::try_new_mandel(
                tensor,
                f64::from(id),
                ContentHash([id; 32]),
                0.25,
            )
            .unwrap())
        };
    CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Top,
        ElasticMaterialField3::ByCell(&field),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut q,
    )
    .unwrap()
}
fn fixture() -> (CutDensityStudy3, Vec<f64>, Vec<f64>, Vec<f64>) {
    let op = operator();
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
