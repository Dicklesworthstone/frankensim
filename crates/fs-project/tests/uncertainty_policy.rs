//! Public admission of fixed-count and sequential native probability studies.

use fs_project::{
    ProjectSpec,
    uncertainty::{CompliancePolicy, UncertaintyStudy, VERSION},
};

const STUDY: &str = r#"(fsim-uncertainty-study
    :version 1 :project "cooling-reference.fsim" :samples 8 :seed 29
    :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :parameters (
        (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
        (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

fn base() -> ProjectSpec {
    fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-reference.fsim"
    )))
    .unwrap()
    .decoded
    .spec
}

const COMPLIANCE: &str =
    ":compliance (bernoulli-mixture :required-probability 0.5 :alpha 0.05 :min-samples 16)";

fn sequential_study() -> String {
    STUDY
        .replace(":version 1", &format!(":version 2 {COMPLIANCE}"))
        .replace(":samples 8", ":samples 64")
}

#[test]
fn version_two_binds_declared_bernoulli_policy_and_roundtrips() {
    let source = sequential_study();
    let study = UncertaintyStudy::parse(&source).unwrap();
    let policy = CompliancePolicy {
        required_probability: 0.5,
        alpha: 0.05,
        min_samples: 16,
    };
    assert_eq!(VERSION, 2);
    assert_eq!(study.compliance(), Some(&policy));
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    let canonical = study.canonical().to_string();
    let bound = study.bind(&base()).unwrap();
    assert_eq!(bound.study().compliance(), Some(&policy));
    assert_eq!(bound.threshold_k(), 348.15);
    assert!(
        bound
            .sample_project(&[5.0, 297.0])
            .unwrap()
            .validate()
            .is_empty()
    );
    for (from, to) in [
        ("0.5 :alpha", "0.6 :alpha"),
        (":alpha 0.05", ":alpha 0.1"),
        (":min-samples 16", ":min-samples 17"),
    ] {
        let changed = UncertaintyStudy::parse(&source.replace(from, to)).unwrap();
        assert_ne!(
            changed.canonical(),
            canonical,
            "stopping policy must bind retained study identity"
        );
    }
}

#[test]
fn study_version_requires_exactly_the_declared_policy() {
    assert!(UncertaintyStudy::parse(&STUDY.replace(":version 1", ":version 2")).is_err());
    assert!(
        UncertaintyStudy::parse(&STUDY.replace(":version 1", &format!(":version 1 {COMPLIANCE}")))
            .is_err()
    );
    let source = sequential_study();
    for (from, to) in [
        ("bernoulli-mixture", "gaussian-mixture"),
        (":compliance (", ":undeclared-policy ("),
        (":min-samples 16", ":min-samples 16 :half-width 0.1"),
        (":alpha 0.05", ":alpha 0.05 :alpha 0.1"),
        (":alpha 0.05 ", ""),
    ] {
        assert!(source.contains(from));
        assert!(
            UncertaintyStudy::parse(&source.replace(from, to)).is_err(),
            "admitted {to}"
        );
    }
    assert!(
        UncertaintyStudy::parse(&source.replace(COMPLIANCE, &format!("{COMPLIANCE} {COMPLIANCE}")))
            .is_err()
    );
}

#[test]
fn compliance_policy_rejects_wrong_types_endpoints_and_invalid_sample_counts() {
    let source = sequential_study();
    for value in [
        "false", "\"0.5\"", "0.5K", "0", "1", "-0.1", "1.1", "NaN", "1e999",
    ] {
        assert!(
            UncertaintyStudy::parse(&source.replace(
                ":required-probability 0.5",
                &format!(":required-probability {value}")
            ))
            .is_err(),
            "admitted probability {value}"
        );
        assert!(
            UncertaintyStudy::parse(&source.replace(":alpha 0.05", &format!(":alpha {value}")))
                .is_err(),
            "admitted alpha {value}"
        );
    }
    assert!(UncertaintyStudy::parse(&source.replace(":alpha 0.05", ":alpha 1e-320")).is_err());
    for value in ["0", "1", "65", "-1", "16.0", "true", "\"16\"", "16s"] {
        assert!(
            UncertaintyStudy::parse(
                &source.replace(":min-samples 16", &format!(":min-samples {value}"))
            )
            .is_err(),
            "admitted min-samples {value}"
        );
    }
    for value in ["2", "64"] {
        assert!(
            UncertaintyStudy::parse(
                &source.replace(":min-samples 16", &format!(":min-samples {value}"))
            )
            .is_ok()
        );
    }
}
