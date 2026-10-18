//! Sampled input reaches the actual native elasticity and optimizer owners.
use super::*;
use super::super::super as native;
use fs_topopt::pipeline::LoadCase;
use std::ops::ControlFlow;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-sampled.fsim"));
const LEGACY: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));

#[test]
fn sampled_input_binds_every_value_to_canonical_identity_and_replay() {
    let spec = native::spec::parse(FIXTURE).unwrap();
    assert!(spec.physical && spec.constructive.is_none());
    let sampled = &spec.sampled.as_ref().unwrap().field;
    assert_eq!(sampled.bounds(), spec.bounds);
    assert_eq!(sampled.value([0.5, 0.5, 0.5]), -0.25);
    assert_eq!(sampled.value([0.0, 0.0, 0.5]), -0.125);
    assert_eq!(native::spec::parse(&spec.canonical).unwrap().id, spec.id);
    let changed = native::spec::parse(&FIXTURE.replacen("-0.75", "-0.8", 1)).unwrap();
    assert_ne!(changed.id, spec.id, "an interior sample is part of geometry identity");
    let old = native::spec::parse(LEGACY).unwrap();
    assert!(old.sampled.is_none() && old.constructive.is_none() && !old.physical);
    assert_eq!(native::spec::parse(&old.canonical).unwrap().id, old.id);
}

#[test]
fn malformed_or_uncovered_grids_refuse_before_physics() {
    for bad in [
        FIXTURE.replace(":shape (3 3 3)", ":shape (3 3 2)"),
        FIXTURE.replace(":shape (3 3 3)", ":shape (1 3 9)"),
        FIXTURE.replace(":shape (3 3 3)", ":shape (3.5 3 3)"),
        FIXTURE.replace(":shape (3 3 3)", ":shape (1024 1024 1024)"),
        FIXTURE.replace(":spacing-m (0.5 0.5 0.5)", ":spacing-m (0.5 0.5 0.0)"),
        FIXTURE.replace(":spacing-m (0.5 0.5 0.5)", ":spacing-m (0.5 0.5 0.25)"),
        FIXTURE.replace(":origin-m (0.0 0.0 0.0)", ":origin-m (0.01 0.0 0.0)"),
        FIXTURE.replace("(1.0 1.0 1.0)", "(1.01 1.0 1.0)"),
        FIXTURE.replace(":values-m (", ":ignored 0 :values-m ("),
        FIXTURE.replacen("-0.75", "\"missing\"", 1),
        FIXTURE.replacen("-0.75", "1e101", 1),
        FIXTURE.replace(":memory-bytes 536870912", ":memory-bytes 134217728"),
    ] {
        assert!(native::spec::parse(&bad).is_err(), "invalid grid was accepted: {bad}");
    }
}

fn response(spec: &native::Spec, level: u8) -> (f64, f64, Vec<f64>) {
    let bounds = native::HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap();
    let tree = native::Octree3::uniform(level, 4, 1000).unwrap();
    let material = native::IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = native::QuadratureControl3::new(
        native::QuadratureOptions3::default(), &mut poll).unwrap();
    // Empty placeholder makes a missed geometry override fail immediately.
    // The same constructor serves stress, compliance, replay and refined grids.
    let placeholder = native::Domain { height: -1.0, curvature: 0.0 };
    let op = native::loading::build_operator(spec, bounds, &tree, &placeholder,
        &material, &mut quadrature).unwrap();
    let volume = op.volumes().iter().sum();
    let force = op.body_load(&|_| spec.loads[0].0, || ControlFlow::Continue(())).unwrap();
    let mut study = native::CutDensityStudy3::new(op, spec.radius, spec.schedule[0]);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = native::SolveControl::new(native::SolveBudget::default(), &mut poll);
    let result = study.evaluate(&vec![0.5; study.cells()],
        &[LoadCase { force: &force, weight: 1.0 }], &mut control).unwrap();
    (volume, result.objective.compliance, result.objective.displacements[0].clone())
}

#[test]
fn sampled_tent_changes_real_volume_and_has_repeatable_equilibrium() {
    let spec = native::spec::parse(FIXTURE).unwrap();
    let a = response(&spec, 1);
    // Height = 5/8 + (1/8) tent(x) tent(y); integral = 21/32.
    assert!((a.0 - 0.65625).abs() < 1e-8, "volume={}", a.0);
    assert!(a.1.is_finite() && a.1 > 0.0);
    assert!(a.2.iter().any(|u| *u != 0.0));
    let again = response(&native::spec::parse(&spec.canonical).unwrap(), 1);
    assert_eq!(a.1.to_bits(), again.1.to_bits());
    assert_eq!(a.2, again.2);
    let flat = FIXTURE.replace("-0.75", "-0.625")
        .replace("-0.25", "-0.125").replace(" 0.25 ", " 0.375 ");
    let flat = response(&native::spec::parse(&flat).unwrap(), 1);
    assert!((flat.0 - 0.625).abs() < 1e-8);
    assert!(a.0 > flat.0, "the interior sample must not be replaced by outer corners");
}

#[test]
fn sampled_equilibrium_preserves_material_scaling_and_refined_geometry() {
    let mut spec = native::spec::parse(FIXTURE).unwrap();
    let a = response(&spec, 1);
    spec.youngs *= 8.0;
    let stiff = response(&spec, 1);
    assert!((stiff.1 * 8.0 - a.1).abs() <= a.1 * 1e-6);
    spec.youngs /= 8.0;
    let refined = response(&spec, 2);
    assert!((refined.0 - a.0).abs() < 1e-8);
    assert!(refined.1.is_finite() && refined.1 > 0.0);
    assert!(refined.2.len() > a.2.len());
}

#[test]
fn sampled_study_runs_a_complete_optimizer_stage() {
    let mut spec = native::spec::parse(FIXTURE).unwrap();
    spec.schedule.truncate(1);
    spec.updates = 1;
    spec.max_level = spec.level;
    let result = native::compute(&spec, &native::CancelGate::new()).unwrap();
    assert_eq!(result.completed_stages, 1);
    assert_eq!(result.requested_stages, 1);
    assert!(result.study.cells() > 0);
}

#[test]
fn sampled_studies_obey_pre_cancellation_and_cumulative_geometry_limits() {
    let spec = native::spec::parse(FIXTURE).unwrap();
    let gate = native::CancelGate::new();
    gate.request();
    let error = match native::compute(&spec, &gate) {
        Err(error) => error,
        Ok(_) => panic!("cancelled sampled study published a result"),
    };
    assert_eq!(error.code, "cli-study-sdf3-cancelled");
    assert_eq!(error.exit, native::exit::CANCELLED);
    let mut limited = spec;
    limited.boxes = 1;
    let error = match native::compute(&limited, &native::CancelGate::new()) {
        Err(error) => error,
        Ok(_) => panic!("sampled construction ignored its geometry allowance"),
    };
    assert_eq!(error.exit, native::exit::BUDGET);
}
