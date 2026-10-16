use super::*;
use super::super::super as native;
use fs_topopt::pipeline::LoadCase;
use std::ops::ControlFlow;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));
const PLATE: &str = "(half-space :normal (0.0 0.0 1.0) :offset-m 0.63)";
const BORED: &str = "(difference :blend-m 0.0 :left (half-space :normal (0.0 0.0 1.0) :offset-m 0.63) :right (cylinder :axis y :center-m (0.55 0.0 0.28) :radius-m 0.12))";

fn shape(text: &str) -> Result<CsgDomain3> {
    super::parse(&fs_ir::sexpr::parse(text).map_err(|e| invalid(e.to_string()))?)
}

fn source(shape: &str) -> String {
    let mut source = FIXTURE.to_string();
    let begin = source.find("  (domain").unwrap();
    let end = source[begin..].find("  (physics").unwrap() + begin;
    source.replace_range(begin..end, &format!(
        "  (domain :type constructive-implicit :bounds ((0.0 0.0 0.0) (1.0 1.0 1.0)) :shape {shape})\n"));
    source
}

#[test]
fn grammar_admits_all_declared_primitives_and_preserves_shape_identity() {
    for text in [
        PLATE,
        BORED,
        "(sphere :center-m (0.1 0.2 0.3) :radius-m 0.2)",
        "(ellipsoid :center-m (0.1 0.2 0.3) :semi-axes-m (0.2 0.3 0.4))",
        "(box :center-m (0.1 0.2 0.3) :half-extents-m (0.2 0.3 0.4))",
        "(cylinder :axis x :center-m (0.1 0.2 0.3) :radius-m 0.2)",
    ] {
        assert!(shape(text).is_ok(), "declared primitive must be reachable: {text}");
    }
    let text = source(BORED);
    let spec = native::spec::parse(&text).unwrap();
    assert!(spec.physical);
    assert_eq!(spec.constructive.as_ref().unwrap().node_count(), 3);
    assert_eq!(native::spec::parse(&spec.canonical).unwrap().id, spec.id);
    let changed = native::spec::parse(&text.replace(":radius-m 0.12", ":radius-m 0.13")).unwrap();
    assert_ne!(changed.id, spec.id, "geometry must bind provenance and replay");
    assert!(native::spec::parse(FIXTURE).unwrap().constructive.is_none());
    assert!(native::spec::parse(&text.replace("(1.0 1.0 1.0)", "(1.0 1.0 -1.0)")).is_err());
}

#[test]
fn malformed_shapes_and_excessive_recipe_work_refuse_before_physics() {
    for bad in [
        "(sphere :center-m (0.0 0.0 0.0) :radius-m -0.2)",
        "(sphere :center-m (0.0 0.0 0.0) :radius-m 0.2 :ignored 1)",
        "(half-space :normal (0.0 0.0 0.0) :offset-m 0.2)",
        "(cylinder :axis diagonal :center-m (0.0 0.0 0.0) :radius-m 0.2)",
        "(mesh :file \"unread.stl\")",
        "(union :blend-m 0.0 :left 0.0 :right 1.0)",
    ] {
        assert!(shape(bad).is_err(), "must not silently ignore invalid geometry: {bad}");
    }
    let mut nested = PLATE.to_string();
    for _ in 1..MAX_DEPTH { nested = format!("(union :blend-m 0.0 :left {nested} :right {PLATE})"); }
    assert!(shape(&nested).is_ok(), "last admitted depth");
    nested = format!("(union :blend-m 0.0 :left {nested} :right {PLATE})");
    assert!(shape(&nested).is_err());
    let mut wide = PLATE.to_string();
    for _ in 0..6 { wide = format!("(union :blend-m 0.0 :left {wide} :right {wide})"); }
    assert_eq!(shape(&wide).unwrap().node_count(), 127);
    assert!(shape(&format!("(union :blend-m 0.0 :left {wide} :right {PLATE})")).is_err());
}

fn response(spec: &native::Spec) -> (f64, f64, Vec<f64>) {
    let bounds = native::HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap();
    let tree = native::Octree3::uniform(1, 4, 1000).unwrap();
    let material = native::IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = native::QuadratureControl3::new(
        native::QuadratureOptions3::default(), &mut poll).unwrap();
    // Deliberately empty legacy placeholder: substituting this for the authored
    // shape would fail. Both native optimizer modes use this exact build seam.
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
fn authored_bore_changes_actual_cut_volume_and_has_replayable_elastic_equilibrium() {
    let plate = native::spec::parse(&source(PLATE)).unwrap();
    let mut bored = native::spec::parse(&source(BORED)).unwrap();
    let full = response(&plate);
    let a = response(&bored);
    assert!((full.0 - 0.63).abs() < 1e-8);
    let removed = full.0 - a.0;
    let expected = std::f64::consts::PI * 0.12 * 0.12;
    assert!((removed - expected).abs() < 0.05 * expected, "{removed} vs {expected}");
    assert!(a.1.is_finite() && a.1 > 0.0);
    assert!(a.2.iter().any(|u| u.abs() > 0.0));
    let again = response(&bored);
    assert_eq!(a.1.to_bits(), again.1.to_bits());
    assert_eq!(a.2, again.2);
    bored.youngs *= 8.0;
    let stiff = response(&bored);
    assert!((stiff.1 * 8.0 - a.1).abs() <= a.1 * 1e-6);
}

#[test]
fn constructive_study_precancellation_never_starts_quadrature_or_equilibrium() {
    let spec = native::spec::parse(&source(BORED)).unwrap();
    let gate = native::CancelGate::new();
    gate.request();
    let error = match native::compute(&spec, &gate) {
        Err(error) => error,
        Ok(_) => panic!("a cancelled constructive study cannot publish a result"),
    };
    assert_eq!(error.code, "cli-study-sdf3-cancelled");
    assert_eq!(error.exit, native::exit::CANCELLED);
}
