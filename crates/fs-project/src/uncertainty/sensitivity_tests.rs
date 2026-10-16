//! Global sensitivity is a fixed independent pick-freeze experiment, not QMC.
use super::*;

const SOURCE: &str = r#"(fsim-uncertainty-study
 :version 1 :project "cooling-reference.fsim" :samples 8 :seed 59
 :wall-time 120s :method sobol-sensitivity :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
 :materials ("aa6061.fsmcdpk") :interfaces ()
 :parameters (
  (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
  (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

#[test]
fn sensitivity_preserves_physical_binding_and_has_its_own_canonical_method() {
    let project = crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec;
    let source = UncertaintyStudy::parse(SOURCE).unwrap();
    assert!(source.sobol_sensitivity());
    assert_eq!(source.samples(), 8, "two rows cost eight real two-input solves");
    assert!(source.qmc().is_none() && source.mean_control().is_none() && source.compliance().is_none());
    assert_eq!(UncertaintyStudy::parse(source.canonical()).unwrap(), source);
    let mc = UncertaintyStudy::parse(&SOURCE.replace("sobol-sensitivity", "monte-carlo")).unwrap();
    assert!(!mc.sobol_sensitivity());
    assert_ne!(source.canonical(), mc.canonical());
    let sensitivity = source.bind(&project).unwrap();
    let ordinary = mc.bind(&project).unwrap();
    for values in [[4.2, 295.0], [5.8, 299.0], [4.2, 299.0]] {
        assert_eq!(sensitivity.sample_project(&values).unwrap(), ordinary.sample_project(&values).unwrap());
    }
    assert_eq!(sensitivity.base(), &project);
}

#[test]
fn sensitivity_refuses_incomplete_designs_constants_and_incompatible_probability_rules() {
    for total in [2, 4, 7, 9, 257] {
        assert!(UncertaintyStudy::parse(&SOURCE.replace(":samples 8", &format!(":samples {total}"))).is_err());
    }
    for changed in [
        SOURCE.replace(":high 6W", ":high 4W"),
        SOURCE.replace("independent", "(gaussian-copula :latent-correlation ((1.0 0.0) (0.0 1.0)))"),
        SOURCE.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)"),
        SOURCE.replace(":version 1", ":version 2 :compliance (bernoulli-mixture :required-probability 0.9 :alpha 0.05 :min-samples 2)"),
        SOURCE.replace(":method sobol-sensitivity", ":method sobol-sensitivity :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 4)"),
    ] { assert!(UncertaintyStudy::parse(&changed).is_err(), "{changed}"); }
    assert!(UncertaintyStudy::parse(&SOURCE.replace(":samples 8", ":samples 256")).is_ok());
}
