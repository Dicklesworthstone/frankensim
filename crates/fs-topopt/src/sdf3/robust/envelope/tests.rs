use super::*;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SimpParams, SolveBudget, SolveProgress};
use crate::sdf3::PhysicalRegion3;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64;3], _: [f64;3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
pub(super) fn fixture(beta: f64) -> (CutDensityStudy3, Vec<f64>, Vec<f64>) {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let op = CutElasticity3::build(HexCell::try_new([0.0;3], [1.0;3]).unwrap(), [2;3], &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap();
    let y = op.body_load(&|p| [0.0, -1.0 - p[1], 0.0], || ControlFlow::Continue(())).unwrap();
    let z = op.body_load(&|p| [0.0, 0.0, -1.0 - p[0]], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(op, 0.15, SimpParams { beta, penal: 2.0, e_min: 1e-3, ..Default::default() }), y, z)
}
pub(super) fn native_body_operator() -> CutElasticity3 {
    use fs_blake3::ContentHash;
    use fs_cutfem::elastic3::ElasticMaterialField3;
    use fs_material::elastic::acoustic::ElasticAcousticState;
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let field = |cell: HexCell| {
        let density = if cell.lo()[0] < 0.5 { 2_u8 } else { 9_u8 };
        let tensor = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
        Ok(ElasticAcousticState::try_new_mandel(tensor, f64::from(density), ContentHash([density; 32]), 1.0).unwrap())
    };
    CutElasticity3::build_with_material_field(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), [2; 3], &Slab,
        ElasticMaterialField3::ByCell(&field), &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut q).unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol * a.abs().max(b.abs()).max(1e-10), "{a} != {b}");
}

#[test]
fn shifted_envelope_has_a_controlled_bias_and_correct_compliance_derivatives() {
    let p = LoadEnvelope3::new(vec![2.0, 4.0, 8.0], 0.1, 3).unwrap();
    let c = [3.0, 6.0, 12.0];
    let e = p.aggregate(&c).unwrap();
    assert_eq!(e.worst_case, 0);
    close(e.value, 1.5 + p.bias_bound(), 1e-14);
    for i in 0..3 {
        close(e.compliance_weights[i], 1.0 / (3.0*p.references()[i]), 1e-14);
        let mut a = c; let mut b = c; a[i] += 1e-5; b[i] -= 1e-5;
        close((p.aggregate(&a).unwrap().value - p.aggregate(&b).unwrap().value)/2e-5,
            e.compliance_weights[i], 1e-8);
    }
    let extreme = LoadEnvelope3::new(vec![1.0;2], 1e-300, 2).unwrap().aggregate(&[1e300, 0.0]).unwrap();
    assert_eq!(extreme.value, 1e300); assert_eq!(extreme.compliance_weights, vec![1.0, 0.0]);
}

#[test]
fn single_case_and_duplicate_loads_keep_actual_independent_fields() {
    let (mut study, _, z) = fixture(2.0);
    let rho = vec![0.45; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    let robust = study.evaluate_robust(&rho, &[LoadCase { force: &z, weight: 1.0 }], RobustOptions3::default(), &mut c).unwrap();
    let policy = LoadEnvelope3::new(vec![1000.0], 0.1, 1).unwrap();
    let one = study.evaluate_robust_envelope(&rho, &[&z], RobustOptions3::default(), &policy, &mut c).unwrap();
    close(one.envelope.value, robust.eroded.objective.compliance/1000.0, 1e-13);
    for (g, reference) in one.gradient.iter().zip(&robust.eroded.objective.gradient) { close(*g, reference/1000.0, 1e-9); }
    let opposite: Vec<_> = z.iter().map(|v| -v).collect();
    let two = study.evaluate_robust_envelope(&rho, &[&z, &opposite], RobustOptions3::default(),
        &LoadEnvelope3::new(vec![1000.0;2], 0.1, 2).unwrap(), &mut c).unwrap();
    close(two.envelope.value, one.envelope.value + 0.1*fs_math::det::ln(2.0), 1e-13);
    assert!(two.envelope.maximum > 0.0);
    assert_eq!(two.physical.eroded.objective.displacements.len(), 2);
    for (a,b) in two.gradient.iter().zip(&one.gradient) { close(*a,*b,1e-8); }
}

#[test]
fn exact_design_gradient_and_dilated_volume_match_each_coordinate() {
    for beta in [1.0, 4.0] {
        let (mut study, y, z) = fixture(beta);
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.37 + 0.025*i as f64).collect();
        let p = LoadEnvelope3::new(vec![1200.0, 1700.0], 0.08, 2).unwrap();
        let mut poll = |_| ControlFlow::Continue(());
        let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
        let e = study.evaluate_robust_envelope(&rho, &[&y,&z], RobustOptions3::default(), &p, &mut c).unwrap();
        for i in 0..rho.len() {
            let h = 1e-4; let mut a = rho.clone(); let mut b = rho.clone(); a[i] += h; b[i] -= h;
            let plus = study.evaluate_robust_envelope(&a, &[&y,&z], RobustOptions3::default(), &p, &mut c).unwrap();
            let minus = study.evaluate_robust_envelope(&b, &[&y,&z], RobustOptions3::default(), &p, &mut c).unwrap();
            close((plus.envelope.value-minus.envelope.value)/(2.0*h), e.gradient[i], 3e-4);
            close((plus.physical.dilated_volume_fraction-minus.physical.dilated_volume_fraction)/(2.0*h),
                e.physical.dilated_volume_gradient[i], 3e-5);
        }
        let check = controlled_envelope_gradient_check(&mut study, &rho, &[&y,&z], RobustOptions3::default(),
            &p, GradientCheckOptions::default(), &mut c).unwrap();
        assert!(check.numerical.passed(), "{check:?}");
    }
}

#[test]
fn protected_regions_and_cancelled_pullback_leave_the_nominal_solution_unchanged() {
    let (study, y, z) = fixture(2.0);
    let mut mask = vec![PhysicalRegion3::Design; study.cells()]; mask[0] = PhysicalRegion3::Solid; mask[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(mask).unwrap();
    let rho = vec![0.4; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&rho, &[LoadCase { force: &z, weight: 1.0 }], &mut c).unwrap();
    let scales = study.operator().scales().to_vec();
    let p = LoadEnvelope3::new(vec![1000.0;2], 0.1, 2).unwrap();
    let e = study.evaluate_robust_envelope(&rho, &[&y,&z], RobustOptions3::default(), &p, &mut c).unwrap();
    for field in [&e.physical.eroded.projected_rho, &e.physical.nominal_rho, &e.physical.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
    assert_eq!(study.operator().scales(), scales);
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-envelope-adjoint" { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.evaluate_robust_envelope(&rho, &[&y,&z], RobustOptions3::default(), &p, &mut c), Err(EvaluationStop::Cancelled)));
    assert!(c.work().linear_iterations > 0);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn invalid_or_overflowing_cases_never_silently_disappear() {
    assert!(LoadEnvelope3::new(vec![], 0.1, 1).is_err());
    assert!(LoadEnvelope3::new(vec![1.0;2], 0.1, 1).is_err());
    for r in [0.0, -1.0, f64::NAN, f64::INFINITY] { assert!(LoadEnvelope3::new(vec![r], 0.1, 1).is_err()); }
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] { assert!(LoadEnvelope3::new(vec![1.0], t, 1).is_err()); }
    let p = LoadEnvelope3::new(vec![f64::MIN_POSITIVE], 0.1, 1).unwrap();
    assert!(p.aggregate(&[f64::MAX]).is_err());
    assert!(p.aggregate(&[f64::NAN]).is_err());
    assert!(p.aggregate(&[-1.0]).is_err());
    assert!(p.aggregate(&[]).is_err());
    let (mut study, _, z) = fixture(1.0);
    let raw = vec![0.4;study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(study.evaluate_robust_envelope(&raw, &[&z], RobustOptions3::default(), &p, &mut c).is_err());
    assert_eq!(c.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
}
