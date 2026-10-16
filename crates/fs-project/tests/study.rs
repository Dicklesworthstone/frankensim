//! Integration tests for `.fsim` study schema parsing, validation, and refusal falsifiers.

use fs_project::{
    parse_study_json, parse_study_sexpr, print_study_json, print_study_sexpr,
    canonical_study_hash,
};

const VALID_STUDY_SEXPR: &str = r#"(fsim-study
  :version 1
  (metadata
    :name "bracket-2d"
    :created "2026-09-01"
    :context-of-use "2D topology optimization marquee benchmark"
    :intended-decision "select optimal material layout under compliance objective"
    :decision-gate scoping-estimate
    :consequence advisory)
  (versions
    :frankensim "0.0.1"
    :schema 1)
  (seeds
    :rng 1337)
  (budgets
    :wall-time 60 s
    :memory 1073741824 B
    :max-iterations 8)
  (capabilities
    "optimization.marquee-topopt"
    "geometry.sdf"
    "physics.cutfem")
  (units
    :system "SI")
  (domain
    :type sdf-plate-with-holes
    :bounds ((0.0 0.0) (1.0 1.0))
    :initial-holes (
      (hole :center (0.3 0.5) :radius 0.12)
      (hole :center (0.7 0.5) :radius 0.18)))
  (physics
    :type elasticity-2d
    :mesh-level 4)
  (scenario
    :fixed-boundary left
    :load-region right)
  (objective
    :type compliance
    :sense minimize
    :unit "J")
  (constraints
    :volume-fraction 0.5)
  (optimizer
    :type projected-gradient
    :step-size 1.0
    :r-min 0.08
    :r-max 0.20
    :steps 8))
"#;

#[test]
fn study_parses_and_validates_cleanly() {
    let spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse study sexpr");
    let violations = spec.validate();
    assert!(violations.is_empty(), "expected 0 violations, got: {violations:?}");

    let printed = print_study_sexpr(&spec);
    let roundtrip = parse_study_sexpr(&printed).expect("parse roundtrip");
    assert_eq!(spec.domain, roundtrip.domain);
    assert_eq!(spec.optimizer, roundtrip.optimizer);
    assert_eq!(spec.constraints, roundtrip.constraints);
}

#[test]
fn study_json_parses_and_validates() {
    let spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse study sexpr");
    let json_str = print_study_json(&spec);
    let roundtrip = parse_study_json(&json_str).expect("parse study json");
    let violations = roundtrip.validate();
    assert!(violations.is_empty(), "expected 0 violations for json, got: {violations:?}");
    assert_eq!(spec.domain, roundtrip.domain);
}

#[test]
fn study_refuses_undeclared_units() {
    let mut spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse");
    spec.units = None;
    let violations = spec.validate();
    assert!(violations.iter().any(|v| v.code == "project-undeclared-units"));
}

#[test]
fn study_refuses_load_on_interior() {
    let mut spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse");
    spec.scenario.as_mut().unwrap().load_region = "interior".to_string();
    let violations = spec.validate();
    assert!(violations.iter().any(|v| v.code == "study-load-non-boundary"));
}

#[test]
fn study_refuses_volume_fraction_out_of_bounds() {
    let mut spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse");
    spec.constraints.as_mut().unwrap().volume_fraction = 1.5;
    let violations = spec.validate();
    assert!(violations.iter().any(|v| v.code == "study-volume-fraction-out-of-bounds"));

    spec.constraints.as_mut().unwrap().volume_fraction = 0.0;
    let violations_zero = spec.validate();
    assert!(violations_zero.iter().any(|v| v.code == "study-volume-fraction-out-of-bounds"));
}

#[test]
fn study_refuses_objective_wrong_dimensions() {
    let mut spec = parse_study_sexpr(VALID_STUDY_SEXPR).expect("parse");
    spec.objective.as_mut().unwrap().unit = "W".to_string(); // Watts (Power) instead of Joules (Energy)
    let violations = spec.validate();
    assert!(violations.iter().any(|v| v.code == "study-objective-dimension-mismatch"));
}

#[test]
fn study_canonical_hash_is_deterministic() {
    let hash1 = canonical_study_hash(VALID_STUDY_SEXPR.as_bytes());
    let hash2 = canonical_study_hash(VALID_STUDY_SEXPR.as_bytes());
    assert_eq!(hash1, hash2);
}
