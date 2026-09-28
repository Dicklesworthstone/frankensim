//! Native Bernoulli decisions use completed physical QoIs, including when
//! deterministic input laws permit the same child solve to be reused.

use super::*;

// Use the numerical owner's actual source without adding a production or
// development dependency just for an independently replayed test oracle.
#[path = "../../../fs-eproc/src/bernoulli.rs"]
mod bernoulli_reference;
use bernoulli_reference::BernoulliMixtureCs;

fn policy_fixture(ambient: f64, cap: usize, minimum: usize, target: f64) -> Fixture {
    let fixture = Fixture::new();
    let source = STUDY
        .replace(":version 1", ":version 2")
        .replace(":samples 4", &format!(":samples {cap}"))
        .replace(":low 4W :high 6W", ":low 5W :high 5W")
        .replace(":low 294K :high 300K", &format!(":low {ambient}K :high {ambient}K"))
        .replace(":parameters (", &format!(
            ":compliance (bernoulli-mixture :required-probability {target} :alpha 0.05 :min-samples {minimum}) :parameters ("));
    std::fs::write(fixture.sources.join("study.fsim"), source).unwrap();
    fixture
}

fn policy_report(fixture: &Fixture, result: &JsonValue, cap: usize) -> (Vec<u8>, JsonValue) {
    assert_eq!(result.str_field("command"), Some("study"));
    let receipt = result.get("receipt").unwrap();
    assert_eq!(
        receipt.str_field("driver"),
        Some("native-cooling-uncertainty-v1")
    );
    assert_eq!(receipt.f64_field("samples_planned"), Some(cap as f64));
    let bytes = artifact(&fixture.ledger, receipt.str_field("report_json").unwrap());
    let report = JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(
        report.get("statistics"),
        Some(&JsonValue::Null),
        "a policy-stopped sample is not a completed fixed-count distribution"
    );
    (bytes, report)
}

/// Reconstruct the confidence process from published physical observations,
/// independently of the CLI policy, checkpoint and reporting implementations.
pub(super) fn verify_confidence(report: &JsonValue, minimum: usize, target: f64) -> Option<usize> {
    let compliance = report.get("compliance").unwrap();
    assert_eq!(
        compliance.str_field("method"),
        Some("bernoulli-beta-half-mixture")
    );
    assert_eq!(compliance.f64_field("required_probability"), Some(target));
    assert_eq!(compliance.f64_field("alpha"), Some(0.05));
    assert_eq!(
        compliance.f64_field("min_decision_samples"),
        Some(minimum as f64)
    );
    assert_eq!(
        compliance.f64_field("samples_evaluated"),
        Some(rows(report).len() as f64)
    );
    let threshold = report.f64_field("temperature_limit_k").unwrap();
    assert_eq!(
        threshold, 297.0,
        "the native margin stays in the physical event"
    );
    let mut oracle = BernoulliMixtureCs::new(0.05).unwrap();
    let mut first_decision = None;
    for (i, row) in rows(report).iter().enumerate() {
        oracle
            .observe(row.f64_field("value_k").unwrap() <= threshold)
            .unwrap();
        let interval = oracle.interval().unwrap().unwrap();
        if i + 1 >= minimum
            && (interval.lo >= target || interval.hi < target)
            && first_decision.is_none()
        {
            first_decision = Some(i + 1);
        }
    }
    let expected = oracle.interval().unwrap().unwrap();
    let interval = compliance
        .get("probability_confidence_sequence")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(interval.len(), 2);
    assert_eq!(
        interval[0].as_f64().unwrap().to_bits(),
        expected.lo.to_bits()
    );
    assert_eq!(
        interval[1].as_f64().unwrap().to_bits(),
        expected.hi.to_bits()
    );
    assert_eq!(
        compliance.f64_field("empirical_probability_of_compliance"),
        Some(expected.mean)
    );
    let decision = if rows(report).len() < minimum {
        "indeterminate"
    } else if expected.lo >= target {
        "meets-probability-target"
    } else if expected.hi < target {
        "below-probability-target"
    } else {
        "indeterminate"
    };
    assert_eq!(compliance.str_field("decision"), Some(decision));
    first_decision
}

#[test]
fn g1_native_all_pass_and_all_fail_decide_from_raw_qois_before_the_cap() {
    for (ambient, passes, decision) in [
        (294.0, true, "meets-probability-target"),
        (300.0, false, "below-probability-target"),
    ] {
        let fixture = policy_fixture(ambient, 16, 8, 0.5);
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        assert_eq!(result.str_field("status"), Some("decision-reached"));
        let (_, report) = policy_report(&fixture, &result, 16);
        assert_eq!(report.str_field("termination"), Some("probability-target"));
        assert_eq!(
            rows(&report).len(),
            8,
            "the declared minimum gates early decisions"
        );
        assert!(rows(&report).len() < 16);
        assert_eq!(
            verify_confidence(&report, 8, 0.5),
            Some(rows(&report).len())
        );
        let compliance = report.get("compliance").unwrap();
        assert_eq!(compliance.str_field("decision"), Some(decision));
        let first = &rows(&report)[0];
        let qoi = json_artifact(&fixture.ledger, first.str_field("qoi_receipt").unwrap());
        let physical = qoi.get("qoi").unwrap().as_array().unwrap()[0]
            .f64_field("value")
            .unwrap();
        assert_eq!(physical <= 297.0, passes);
        for row in rows(&report) {
            assert_eq!(
                row.f64_field("value_k").unwrap().to_bits(),
                physical.to_bits()
            );
            assert_eq!(row.str_field("run"), first.str_field("run"));
            assert_eq!(row.str_field("qoi_receipt"), first.str_field("qoi_receipt"));
        }
        let interval = compliance
            .get("probability_confidence_sequence")
            .unwrap()
            .as_array()
            .unwrap();
        assert!(
            interval[1].as_f64().unwrap() > interval[0].as_f64().unwrap(),
            "identical samples do not eliminate probability uncertainty"
        );
        let exported = command(
            &[
                "--json",
                "package",
                result.str_field("run").unwrap(),
                fixture.ledger.to_str().unwrap(),
            ],
            fs_cli::exit::SUCCESS,
        );
        let package = std::fs::read(exported.str_field("package").unwrap()).unwrap();
        assert_eq!(
            package,
            artifact(
                &fixture.ledger,
                result.get("receipt").unwrap().str_field("package").unwrap()
            )
        );
        let package = std::str::from_utf8(&package).unwrap();
        assert!(package.contains("cooling.uncertainty.compliance-probability"));
        assert!(package.contains(decision));
        let repeated = fixture.resume(
            result.str_field("run").unwrap(),
            None,
            fs_cli::exit::SUCCESS,
        );
        assert_eq!(repeated, result, "a resolved policy is terminal");
    }
}

#[test]
fn g4_g5_native_decision_resume_keeps_the_original_policy_and_exact_report() {
    let full = policy_fixture(294.0, 16, 8, 0.5);
    let completed = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, _) = policy_report(&full, &completed, 16);
    let split = policy_fixture(294.0, 16, 8, 0.5);
    let zero = split.study(Some("0"), fs_cli::exit::BUDGET);
    let (_, report) = policy_report(&split, &zero, 16);
    assert!(rows(&report).is_empty());
    let compliance = report.get("compliance").unwrap();
    assert_eq!(
        compliance.get("probability_confidence_sequence"),
        Some(&JsonValue::Null)
    );
    assert_eq!(
        compliance.get("empirical_probability_of_compliance"),
        Some(&JsonValue::Null)
    );
    assert_eq!(compliance.str_field("decision"), Some("indeterminate"));
    let first = split.resume(
        zero.str_field("run").unwrap(),
        Some("3"),
        fs_cli::exit::BUDGET,
    );
    let (_, report) = policy_report(&split, &first, 16);
    assert_eq!(rows(&report).len(), 3);
    assert_eq!(verify_confidence(&report, 8, 0.5), None);
    // Retain the original file but change the policy at its former path. Resume
    // is bound to the retained declaration, never the mutable source pathname.
    let source_path = split.sources.join("study.fsim");
    let changed = std::fs::read_to_string(&source_path)
        .unwrap()
        .replace(":required-probability 0.5", ":required-probability 0.99");
    std::fs::rename(&source_path, split.sources.join("original-policy.fsim")).unwrap();
    std::fs::write(&source_path, changed).unwrap();
    let second = split.resume(
        first.str_field("run").unwrap(),
        Some("4"),
        fs_cli::exit::BUDGET,
    );
    let (_, report) = policy_report(&split, &second, 16);
    assert_eq!(rows(&report).len(), 7);
    assert_eq!(verify_confidence(&report, 8, 0.5), None);
    let compliance = report.get("compliance").unwrap();
    assert!(
        compliance
            .get("probability_confidence_sequence")
            .unwrap()
            .as_array()
            .unwrap()[0]
            .as_f64()
            .unwrap()
            >= 0.5
    );
    assert_eq!(
        compliance.str_field("decision"),
        Some("indeterminate"),
        "the confidence boundary alone cannot bypass the declared minimum"
    );
    let finished = split.resume(
        second.str_field("run").unwrap(),
        None,
        fs_cli::exit::SUCCESS,
    );
    assert_eq!(finished.str_field("status"), Some("decision-reached"));
    let (actual, report) = policy_report(&split, &finished, 16);
    assert_eq!(actual, expected);
    assert_eq!(verify_confidence(&report, 8, 0.5), Some(8));
}

#[test]
fn g4_native_probability_target_unresolved_at_the_lifetime_cap_stays_partial() {
    let fixture = policy_fixture(294.0, 8, 8, 0.99);
    let result = fixture.study(None, fs_cli::exit::BUDGET);
    assert_eq!(result.str_field("status"), Some("budget-exhausted"));
    let (_, report) = policy_report(&fixture, &result, 8);
    assert_eq!(
        report.str_field("termination"),
        Some("lifetime-sample-budget")
    );
    assert_eq!(rows(&report).len(), 8);
    assert_eq!(verify_confidence(&report, 8, 0.99), None);
    let compliance = report.get("compliance").unwrap();
    assert_eq!(
        compliance.f64_field("empirical_probability_of_compliance"),
        Some(1.0)
    );
    assert_eq!(compliance.str_field("decision"), Some("indeterminate"));
    let repeated = fixture.resume(result.str_field("run").unwrap(), None, fs_cli::exit::BUDGET);
    assert_eq!(
        repeated, result,
        "resume cannot silently raise the lifetime sample cap"
    );
}
