use super::*;
use crate::pipeline::LoadCase;
use crate::sdf3::PhysicalRegion3;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{SolveBudget, SolveProgress};
use std::ops::ControlFlow;

fn caps() -> ComplianceCaps3 { ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap() }
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 1e-8 * a.abs().max(b.abs()).max(1.0), "{a} != {b}"); }

#[test]
fn binary_material_body_force_matches_native_density_and_exact_occupancy() {
    use crate::sdf3::robust::envelope::tests::native_body_operator;
    let op = native_body_operator();
    let acceleration = [0.1, -0.2, 0.3];
    let prepared = op.prepare_material_body_load(&|_| acceleration, 100_000,
        || ControlFlow::Continue(())).unwrap();
    let dead = op.body_load(&|_| [0.0, 0.0, -0.05], || ControlFlow::Continue(())).unwrap();
    let opposite: Vec<_> = dead.iter().map(|f| -f).collect();
    let regions = op.cell_keys().iter().map(|key|
        if key[0] == 0 { PhysicalRegion3::Solid } else { PhysicalRegion3::Void }).collect();
    let params = SimpParams { beta: 2.0, penal: 3.0, e_min: 0.02, ..Default::default() };
    let mut study = CutDensityStudy3::new(op, 0.15, params)
        .with_physical_regions(regions).unwrap().with_material_body_load(prepared).unwrap();
    let incoming = study.operator().scales().to_vec();
    let raw = vec![0.6; study.cells()];
    let options = BinaryComplianceOptions3 { void_stiffness: 0.07, ..Default::default() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let binary = study.evaluate_binary_compliance_caps(&raw, &[&dead, &opposite],
        &caps(), options, &mut control).unwrap();
    assert!(binary.eroded.solid.iter().any(|s| *s) && binary.eroded.solid.iter().any(|s| !s));
    assert_eq!(study.operator().scales(), incoming);
    // Independent native integration uses full rho0 only where the binary
    // material is present. The heavy void half contributes no body force.
    let mut reference = native_body_operator();
    reference.set_scales(&binary.eroded_stiffness_scales()).unwrap();
    let body = reference.body_load(&|p| {
        let density = if p[0] < 0.5 { 2.0 } else { 0.0 };
        acceleration.map(|a| density * a)
    }, || ControlFlow::Continue(())).unwrap();
    for (case, external) in [&dead, &opposite].into_iter().enumerate() {
        let rhs: Vec<_> = external.iter().zip(&body).map(|(f, b)| f + b).collect();
        let solution = reference.solve_controlled(&rhs, 1e-11, 10_000, 16,
            |_| ControlFlow::Continue(())).unwrap();
        close(binary.case_compliances[case], solution.compliance());
        for (&a, &b) in binary.displacements[case].iter().zip(solution.coefficients()) { close(a, b); }
    }
    assert!((binary.case_compliances[0] - binary.case_compliances[1]).abs() > 1e-6,
        "a common body force must not be negated with each authored dead load");
}

#[test]
fn affine_body_force_refuses_continuous_spectral_endpoints_and_scans_before_work() {
    use crate::sdf3::robust::envelope::tests::native_body_operator;
    use crate::sdf3::robust::envelope::load_space::LoadSpace3;
    use crate::sdf3::robust::envelope::load_space::design::limits::{
        LoadSpaceCap3, binary::BinaryLoadSpaceScanOptions3,
    };
    let op = native_body_operator();
    let load = op.prepare_material_body_load(&|_| [0.0, 0.0, -1.0], 100_000,
        || ControlFlow::Continue(())).unwrap();
    let force = op.body_load(&|_| [0.0, -1.0, 0.0], || ControlFlow::Continue(())).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.15, SimpParams::default())
        .with_material_body_load(load).unwrap();
    let raw = vec![0.6; study.cells()];
    let basis = [&force[..], &force[..]];
    let scales = study.operator().scales().to_vec();
    let polls = std::cell::Cell::new(0);
    let mut callback = |_| { polls.set(polls.get() + 1); ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let work = control.work();
    let policy = LoadSpace3::new(2, 1e6, 0.02).unwrap();
    assert!(matches!(study.evaluate_robust_load_space(&raw, &basis,
        RobustOptions3::default(), &policy, &mut control), Err(EvaluationStop::Breakdown { .. })));
    let cap = LoadSpaceCap3::new(2, 1e6, 0.02).unwrap();
    assert!(matches!(study.evaluate_binary_load_space_cap(&raw, &basis, &cap,
        BinaryComplianceOptions3::default(), &mut control), Err(BinaryComplianceError3::Invalid(_))));
    assert!(matches!(study.scan_binary_load_space_cap(&raw, &basis, &cap, &[0.4, 0.6],
        BinaryLoadSpaceScanOptions3::default(), &mut control), Err(BinaryComplianceError3::Invalid(_))));
    assert!(matches!(study.scan_all_binary_load_space_cap(&raw, &basis, &cap,
        BinaryLoadSpaceScanOptions3::default(), &mut control), Err(BinaryComplianceError3::Invalid(_))));
    assert_eq!(control.work(), work);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(polls.get(), 0);
}

#[test]
fn binary_eroded_solves_match_an_independently_prescribed_material_image() {
    let (mut study, y, z) = fixture(2.0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.045 * i as f64).collect();
    let previous = study.operator().scales().to_vec();
    let options = BinaryComplianceOptions3 { threshold: 0.45, void_stiffness: 0.02, ..Default::default() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let binary = study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(), options, &mut control).unwrap();
    // One source filter, two physical equilibria, no fictitious adjoint solves.
    assert_eq!(binary.work.linear_solves, 3);
    assert_eq!(study.operator().scales(), previous);
    assert!(binary.eroded.solid.iter().any(|s| *s) && binary.eroded.solid.iter().any(|s| !s));
    for ((&e, &n), &d) in binary.eroded.solid.iter().zip(&binary.nominal.solid).zip(&binary.dilated.solid) {
        assert!(!e || n); assert!(!n || d);
    }
    let (mut reference, _, _) = fixture(2.0);
    reference.params.e_min = options.void_stiffness;
    let regions = binary.eroded.solid.iter().map(|&s|
        if s { PhysicalRegion3::Solid } else { PhysicalRegion3::Void }).collect();
    reference = reference.with_physical_regions(regions).unwrap();
    let independent = reference.evaluate(&raw, &[LoadCase { force: &y, weight: 0.5 },
        LoadCase { force: &z, weight: 0.5 }], &mut control).unwrap();
    for (a, b) in binary.case_compliances.iter().zip(&independent.objective.case_compliances) { close(*a, *b); }
    for (u, v) in binary.displacements.iter().zip(&independent.objective.displacements) {
        for (a, b) in u.iter().zip(v) { close(*a, *b); }
    }
    assert_eq!(binary.eroded_stiffness_scales(), reference.operator().scales());
    let tight = ComplianceCaps3::new(binary.case_compliances.iter().map(|c| c*0.999999).collect(), 0.02, 2).unwrap();
    let violating = study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &tight, options, &mut control).unwrap();
    assert!(violating.envelope.value - 1.0 < 0.1);
    assert!(!violating.feasible(0.1), "smooth tolerance must not waive a physical cap");
    assert!(binary.eroded.volume_fraction <= binary.nominal.volume_fraction);
    assert!(binary.nominal.volume_fraction <= binary.dilated.volume_fraction);
}

#[test]
fn a_gray_cap_pass_does_not_pass_a_rounded_endpoint() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.6; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let gray = study.evaluate_robust(&raw, &[LoadCase { force: &y, weight: 0.5 },
        LoadCase { force: &z, weight: 0.5 }], RobustOptions3::default(), &mut control).unwrap();
    let caps = ComplianceCaps3::new(gray.eroded.objective.case_compliances.iter().map(|c| 1.1*c).collect(), 0.02, 2).unwrap();
    assert!(caps.satisfied_by(&gray.eroded.objective.case_compliances));
    let binary = study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps,
        BinaryComplianceOptions3 { threshold: 1.0, ..Default::default() }, &mut control).unwrap();
    assert!(binary.eroded.solid.iter().all(|s| !s));
    assert!(binary.case_compliances.iter().zip(caps.caps()).all(|(c, cap)| c > cap));
    assert!(!binary.feasible(0.0)); assert!(!binary.feasible(0.5));
}

#[test]
fn ties_are_not_split_and_physical_regions_survive_all_realizations() {
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let tied = Design3 { projected: vec![0.5; 3], slope: vec![0.0; 3], scales: vec![1.0; 3], volume: 0.5 };
    let included = material(&tied, &[0.2, 0.3, 0.5], 0.5, &mut control).unwrap();
    assert_eq!(included.solid, vec![true; 3]); assert_eq!(included.volume_fraction, 1.0);
    let excluded = material(&tied, &[0.2, 0.3, 0.5], 0.5_f64.next_up(), &mut control).unwrap();
    assert_eq!(excluded.solid, vec![false; 3]); assert_eq!(excluded.volume_fraction, 0.0);
    let (study, y, z) = fixture(2.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let raw = vec![0.4; study.cells()];
    for threshold in [1e-6, 0.5, 1.0] {
        let b = study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(),
            BinaryComplianceOptions3 { threshold, ..Default::default() }, &mut control).unwrap();
        for field in [&b.eroded, &b.nominal, &b.dilated] { assert!(field.solid[0]); assert!(!field.solid[1]); }
    }
}

#[test]
fn zero_erosion_and_opposite_loads_retain_independent_results() {
    let (mut study, _, z) = fixture(2.0);
    let raw = vec![0.7; study.cells()]; let opposite: Vec<_> = z.iter().map(|f| -f).collect();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let b = study.evaluate_binary_compliance_caps(&raw, &[&z, &opposite], &caps(),
        BinaryComplianceOptions3 { robustness: RobustOptions3 { threshold_offset: 0.0 }, ..Default::default() }, &mut control).unwrap();
    assert_eq!(b.eroded, b.nominal); assert_eq!(b.nominal, b.dilated);
    assert!(b.case_compliances[0] > 0.0); close(b.case_compliances[0], b.case_compliances[1]);
    for (&u, &v) in b.displacements[0].iter().zip(&b.displacements[1]) { close(u, -v); }
}

#[test]
fn cancellation_or_exhaustion_preserves_a_preexisting_nominal_solution() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.6; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&raw, &[LoadCase { force: &y, weight: 1.0 }], &mut control).unwrap();
    let scales = study.operator().scales().to_vec();
    for stage in ["sdf3-binary-cap-elasticity", "sdf3-binary-cap-publish"] {
        let mut poll = |p: SolveProgress| if p.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        assert!(matches!(study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(),
            BinaryComplianceOptions3::default(), &mut control),
            Err(BinaryComplianceError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(study.operator().scales(), scales);
    }
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(),
        BinaryComplianceOptions3::default(), &mut control).is_err());
    assert_eq!(control.work().linear_iterations, 1); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn invalid_policy_and_load_shapes_are_refused_before_solves() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.6; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    for threshold in [0.0, -1.0, 1.1, f64::NAN] {
        let options = BinaryComplianceOptions3 { threshold, ..Default::default() };
        assert!(matches!(study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(), options, &mut control), Err(BinaryComplianceError3::Invalid(_))));
    }
    for void_stiffness in [0.0, 1.0, f64::INFINITY, f64::NAN] {
        let options = BinaryComplianceOptions3 { void_stiffness, ..Default::default() };
        assert!(study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(), options, &mut control).is_err());
    }
    assert!(study.evaluate_binary_compliance_caps(&raw, &[&y], &caps(), BinaryComplianceOptions3::default(), &mut control).is_err());
    assert!(study.evaluate_binary_compliance_caps(&raw, &[&y[..y.len()-1], &z], &caps(), BinaryComplianceOptions3::default(), &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
}

#[test]
fn adaptive_preconditioned_endpoint_matches_its_prescribed_binary_material() {
    use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
    use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
    use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
    use fs_cutfem::octree3::Octree3;
    use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
    use fs_ivl::Interval;
    use fs_material::IsotropicElastic;
    struct Slab;
    impl CutSdf3 for Slab {
        fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
        fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
            Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
        }
        fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
            let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
        }
    }
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree.refined(&[*tree.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap();
    let build = || {
        let mut poll = |_| ControlFlow::Continue(());
        let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
        AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), &tree, &Slab,
            &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
            ElasticityOptions3::default(), &mut q).unwrap()
    };
    let op = build();
    let force = op.body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(())).unwrap();
    let params = SimpParams { beta: 2.0, penal: 2.0, e_min: 0.02, ..Default::default() };
    let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15, params);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.3 + 0.4*(i as f64)/(study.cells() as f64)).collect();
    let scales = study.operator().scales().to_vec();
    let caps = ComplianceCaps3::new(vec![1e6], 0.02, 1).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let binary = study.evaluate_binary_compliance_caps(&raw, &[&force], &caps,
        BinaryComplianceOptions3 { threshold: 0.45, void_stiffness: 0.02, ..Default::default() }, &mut control).unwrap();
    let regions = binary.eroded.solid.iter().map(|&s| if s { PhysicalRegion3::Solid } else { PhysicalRegion3::Void }).collect();
    let mut reference = CutDensityStudy3::new(build(), 0.15, params).with_physical_regions(regions).unwrap();
    let fresh = reference.evaluate(&raw, &[LoadCase { force: &force, weight: 1.0 }], &mut control).unwrap();
    close(binary.case_compliances[0], fresh.objective.compliance);
    for (u, v) in binary.displacements[0].iter().zip(&fresh.objective.displacements[0]) { close(*u, *v); }
    assert_eq!(study.operator().scales(), scales);
}
