//! Native material self-weight through the complete CutFEM design derivatives.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{ElasticMaterialField3, MaterialBodyLoad3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::response::limits::binary::BinaryEndpointOptions3;
use fs_topopt::sdf3::response::limits::{
    BoundedResponse3, ResponseLimit3, ResponseLimitCase3, ResponseLimitSense3,
};
use fs_topopt::sdf3::robust::RobustOptions3;
use fs_topopt::sdf3::robust::envelope::{LoadEnvelope3, load_space::LoadSpace3};
use fs_topopt::sdf3::stress::{StressMeasure3, StressOptions3};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::{SimpParams, SolveBudget, SolveControl};

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
fn operator() -> AdaptiveSolveSpace3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
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
        &tree,
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
fn fixture(eta: f64, protected: bool) -> (Study, MaterialBodyLoad3) {
    let op = operator();
    let body = op
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let audit = op
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    let mut study = CutDensityStudy3::new(
        op,
        0.15,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.03,
            eta,
        },
    )
    .with_material_body_load(body)
    .unwrap();
    if protected {
        let mut regions = vec![PhysicalRegion3::Design; study.cells()];
        regions[0] = PhysicalRegion3::Solid;
        let last = regions.len() - 1;
        regions[last] = PhysicalRegion3::Void;
        study = study.with_physical_regions(regions).unwrap();
    }
    (study, audit)
}
fn raw(study: &Study) -> Vec<f64> {
    (0..study.cells())
        .map(|i| 0.33 + 0.055 * i as f64)
        .collect()
}
fn perturb(rho: &[f64], d: &[f64], h: f64) -> Vec<f64> {
    rho.iter().zip(d).map(|(r, d)| r + h * d).collect()
}
fn directions(n: usize) -> [Vec<f64>; 2] {
    [
        vec![1.0; n],
        (0..n).map(|i| ((i % 3) as f64 - 1.0) / 2.0).collect(),
    ]
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * a.abs().max(b.abs()).max(1e-10),
        "{a} != {b}"
    );
}
fn derivative(fd: f64, actual: f64, value: f64) {
    let scale = fd
        .abs()
        .max(actual.abs())
        .max(value.abs() * 1e-5)
        .max(1e-12);
    assert!(
        (fd - actual).abs() <= 4e-4 * scale,
        "FD {fd}, adjoint {actual}, scale {scale}"
    );
}

#[test]
fn compliance_has_both_body_force_derivatives_and_the_protected_projection() {
    let (mut s, audit) = fixture(0.5, true);
    assert!(s.has_material_body_load());
    let rho = raw(&s);
    let zero = vec![0.0; s.operator().n()];
    let loads = [LoadCase {
        force: &zero,
        weight: 1.0,
    }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let value = s.evaluate(&rho, &loads, &mut control).unwrap();
    assert!(value.objective.compliance > 1e-10);
    assert_eq!(value.projected_rho[0], 1.0);
    assert_eq!(value.projected_rho[rho.len() - 1], 0.0);
    let rhs = audit.load(&value.projected_rho, poll).unwrap();
    let u = &value.objective.displacements[0];
    close(value.objective.compliance, dot(&rhs, u));
    assert!(
        s.operator()
            .elasticity()
            .field_residual(u, &rhs, poll)
            .unwrap()
            < 1e-9
    );
    for i in 0..rho.len() {
        let mut d = vec![0.0; rho.len()];
        d[i] = 1.0;
        let h = 2e-5;
        let plus = s
            .evaluate(&perturb(&rho, &d, h), &loads, &mut control)
            .unwrap();
        let minus = s
            .evaluate(&perturb(&rho, &d, -h), &loads, &mut control)
            .unwrap();
        derivative(
            (plus.objective.compliance - minus.objective.compliance) / (2.0 * h),
            value.objective.gradient[i],
            value.objective.compliance,
        );
    }
    // Exact authored void has neither physical mass nor self-weight, despite
    // the retained positive regularization stiffness in every cell.
    let (empty, _) = fixture(0.5, false);
    let n = empty.cells();
    let mut empty = empty
        .with_physical_regions(vec![PhysicalRegion3::Void; n])
        .unwrap();
    let empty_value = empty.evaluate(&vec![0.7; n], &loads, &mut control).unwrap();
    assert_eq!(empty_value.objective.compliance, 0.0);
    assert!(empty_value.objective.gradient.iter().all(|g| *g == 0.0));
    assert!(
        empty_value.objective.displacements[0]
            .iter()
            .all(|u| *u == 0.0)
    );
}

#[test]
fn displacement_reaction_and_motion_stress_include_the_same_body_adjoint() {
    let (mut s, _) = fixture(0.5, true);
    let rho = raw(&s);
    let zero = vec![0.0; s.operator().n()];
    let q = s
        .operator()
        .elasticity()
        .body_load(&|p| [0.0, 0.0, 1.0 + p[1]], poll)
        .unwrap();
    let before = s.operator().elasticity().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    for reaction in [false, true] {
        let rows = [ResponseLimit3 {
            response: if reaction {
                BoundedResponse3::Reaction(&right)
            } else {
                BoundedResponse3::Displacement(&q)
            },
            sense: ResponseLimitSense3::AtMost,
            bound: 0.0,
            scale: 0.01,
        }];
        let cases = [ResponseLimitCase3 {
            force: &zero,
            prescribed: Some(&motion),
            limits: &rows,
        }];
        let value = s
            .evaluate_response_limits(&rho, &cases, Default::default(), &mut control)
            .unwrap();
        assert!(value.measurements[0][0].abs() > 1e-8);
        for d in directions(rho.len()) {
            let h = 2e-5;
            let a = s
                .evaluate_response_limits(
                    &perturb(&rho, &d, h),
                    &cases,
                    Default::default(),
                    &mut control,
                )
                .unwrap();
            let b = s
                .evaluate_response_limits(
                    &perturb(&rho, &d, -h),
                    &cases,
                    Default::default(),
                    &mut control,
                )
                .unwrap();
            derivative(
                (a.constraint - b.constraint) / (2.0 * h),
                dot(&d, &value.constraint_gradient),
                value.constraint,
            );
        }
    }
    let loads = [LoadCase {
        force: &zero,
        weight: 1.0,
    }];
    let options = StressOptions3 {
        aggregation_power: 4.0,
        ..Default::default()
    };
    let value = s
        .evaluate_stress_with_prescribed_motion(
            &rho,
            &loads,
            &[Some(&motion)],
            options,
            StressMeasure3::NormalizedAverage,
            &mut control,
        )
        .unwrap();
    assert!(value.aggregate > 1e-6);
    for d in directions(rho.len()) {
        let h = 2e-5;
        let a = s
            .evaluate_stress_with_prescribed_motion(
                &perturb(&rho, &d, h),
                &loads,
                &[Some(&motion)],
                options,
                StressMeasure3::NormalizedAverage,
                &mut control,
            )
            .unwrap();
        let b = s
            .evaluate_stress_with_prescribed_motion(
                &perturb(&rho, &d, -h),
                &loads,
                &[Some(&motion)],
                options,
                StressMeasure3::NormalizedAverage,
                &mut control,
            )
            .unwrap();
        derivative(
            (a.aggregate - b.aggregate) / (2.0 * h),
            dot(&d, &value.gradient),
            value.aggregate,
        );
    }
    assert_eq!(s.operator().elasticity().scales(), before);
}

#[test]
fn binary_self_weight_uses_zero_one_material_not_ersatz_stiffness() {
    let (mut s, _) = fixture(0.5, true);
    let rho = raw(&s);
    let zero = vec![0.0; s.operator().n()];
    let rows = [ResponseLimit3 {
        response: BoundedResponse3::Reaction(&right),
        sense: ResponseLimitSense3::AtMost,
        bound: 1.0,
        scale: 1.0,
    }];
    let cases = [ResponseLimitCase3 {
        force: &zero,
        prescribed: Some(&motion),
        limits: &rows,
    }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let before = s.operator().elasticity().scales().to_vec();
    let endpoint = s
        .evaluate_binary_response_limits(
            &rho,
            &cases,
            BinaryEndpointOptions3 {
                threshold: 0.5,
                void_stiffness: 0.1,
                ..Default::default()
            },
            &mut control,
        )
        .unwrap();
    assert!(endpoint.solid[0] && !endpoint.solid[rho.len() - 1]);
    let mut native = operator();
    let audit = native
        .elasticity()
        .prepare_material_body_load(&acceleration, 100_000, poll)
        .unwrap();
    native.set_scales(&endpoint.stiffness_scales()).unwrap();
    let occupancy: Vec<_> = endpoint
        .solid
        .iter()
        .map(|s| if *s { 1.0 } else { 0.0 })
        .collect();
    let lifting = native
        .elasticity()
        .prescribed_displacement_load(&motion, poll)
        .unwrap();
    let rhs: Vec<_> = audit
        .load(&occupancy, poll)
        .unwrap()
        .iter()
        .zip(&lifting)
        .map(|(b, g)| b + g)
        .collect();
    let u = &endpoint.displacements[0];
    assert!(native.elasticity().field_residual(u, &rhs, poll).unwrap() < 1e-9);
    close(
        endpoint.measurements[0][0],
        native
            .elasticity()
            .embedded_reaction(u, Some(&motion), &right, poll)
            .unwrap()
            .value,
    );
    let wrong: Vec<_> = audit
        .load(&endpoint.stiffness_scales(), poll)
        .unwrap()
        .iter()
        .zip(&lifting)
        .map(|(b, g)| b + g)
        .collect();
    assert!(native.elasticity().field_residual(u, &wrong, poll).unwrap() > 1e-5);
    assert_eq!(s.operator().elasticity().scales(), before);
}

#[test]
fn robust_erosion_and_finite_load_envelope_use_the_eroded_body_occupancy() {
    let (mut s, _) = fixture(0.5, true);
    let rho = raw(&s);
    let zero = vec![0.0; s.operator().n()];
    let loads = [LoadCase {
        force: &zero,
        weight: 1.0,
    }];
    let options = RobustOptions3::default();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let before = s.operator().elasticity().scales().to_vec();
    let robust = s
        .evaluate_robust(&rho, &loads, options, &mut control)
        .unwrap();
    let (mut eroded, audit) = fixture(0.5 + options.threshold_offset, true);
    let direct = eroded.evaluate(&rho, &loads, &mut control).unwrap();
    assert_eq!(robust.eroded.projected_rho, direct.projected_rho);
    close(
        robust.eroded.objective.compliance,
        direct.objective.compliance,
    );
    let rhs = audit.load(&direct.projected_rho, poll).unwrap();
    assert!(
        eroded
            .operator()
            .elasticity()
            .field_residual(&robust.eroded.objective.displacements[0], &rhs, poll)
            .unwrap()
            < 1e-9
    );
    let d = directions(rho.len())[1].clone();
    let h = 2e-5;
    let a = s
        .evaluate_robust(&perturb(&rho, &d, h), &loads, options, &mut control)
        .unwrap();
    let b = s
        .evaluate_robust(&perturb(&rho, &d, -h), &loads, options, &mut control)
        .unwrap();
    derivative(
        (a.eroded.objective.compliance - b.eroded.objective.compliance) / (2.0 * h),
        dot(&d, &robust.eroded.objective.gradient),
        robust.eroded.objective.compliance,
    );
    let fixed = s
        .operator()
        .elasticity()
        .body_load(&|_| [0.001, 0.0, 0.0], poll)
        .unwrap();
    let family = [&zero[..], &fixed[..]];
    let envelope = LoadEnvelope3::new(vec![1e-3, 2e-3], 0.03, 2).unwrap();
    let value = s
        .evaluate_robust_envelope(&rho, &family, options, &envelope, &mut control)
        .unwrap();
    let a = s
        .evaluate_robust_envelope(
            &perturb(&rho, &d, h),
            &family,
            options,
            &envelope,
            &mut control,
        )
        .unwrap();
    let b = s
        .evaluate_robust_envelope(
            &perturb(&rho, &d, -h),
            &family,
            options,
            &envelope,
            &mut control,
        )
        .unwrap();
    derivative(
        (a.envelope.value - b.envelope.value) / (2.0 * h),
        dot(&d, &value.gradient),
        value.envelope.value,
    );
    // An affine self-weight offset is not a linear spectral load basis.
    let spectral = LoadSpace3::new(2, 1e-3, 0.03).unwrap();
    let work = control.work();
    assert!(
        s.evaluate_robust_load_space(&rho, &family, options, &spectral, &mut control)
            .is_err()
    );
    assert_eq!(control.work(), work);
    assert_eq!(s.operator().elasticity().scales(), before);
}
