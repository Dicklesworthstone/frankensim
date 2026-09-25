use super::*;
use super::super::super::super::parse as study_spec;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-projected-stress-2d.fsim"));
const PROTECTED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-protected-regions-2d.fsim"));
const FAMILY: &str = "    :load-family (independent\n      :aggregate weighted-sum\n      :max-solves 512\n      :max-recovery-solves 64\n      :additional (\n        (case :band (0.375 0.625) :traction-pa (0.0 1.0) :weight 1.0)\n      ))";
fn source() -> String {
    BASE.replace("    :stress-tolerance-pa 0.0)",
        &format!("    :stress-tolerance-pa 0.0\n{FAMILY}\n  )"))
}
fn spec() -> ElasticitySpec { study_spec(&source()).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn constraints(value: &JsonValue) -> &JsonValue { value.path(&["continuation", "constraints"]).unwrap() }
fn history(value: &JsonValue) -> &JsonValue { constraints(value).get("load_family").unwrap() }
fn restored(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, ConstraintEvidence) {
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (geometry, report) = decode(spec, &receipt, &design, &rows).unwrap();
    let evidence = ConstraintEvidence::read(constraints(&receipt), &report, spec.projected.as_ref().unwrap()).unwrap();
    (geometry, evidence)
}

#[test]
fn independent_load_declarations_are_complete_bounded_and_do_not_change_legacy_sources() {
    let legacy = study_spec(BASE).unwrap();
    assert!(stress_controls(&legacy).unwrap().family.is_none());
    assert!(!legacy.canonical.contains("load-family"));
    let original = spec();
    assert_eq!(study_spec(&original.canonical).unwrap().id, original.id);
    for (from, to) in [
        (":aggregate weighted-sum", ":aggregate summed-forces"),
        (":max-solves 512", ":max-solves 1"),
        (":max-recovery-solves 64", ":max-recovery-solves 100001"),
        (":traction-pa (0.0 1.0)", ":traction-pa (0.0 0.0)"),
        (":weight 1.0", ":weight -1.0"),
        (":weight 1.0", ":weight 1.0 :weight 2.0"),
        (":max-solves 512", ""),
    ] { assert!(study_spec(&source().replace(from, to)).is_err(), "{to}"); }
}

#[test]
fn opposite_operating_conditions_do_not_cancel_and_zero_weight_still_limits_stress() {
    for aggregate in ["weighted-sum", "worst-weighted-case"] {
        let spec = study_spec(&source().replace("weighted-sum", aggregate)).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let out = drive(&spec, &ledger, Some(0), &gate(), None).unwrap();
        let (_, evidence) = restored(&ledger, &out, &spec);
        let cases = &evidence.family.as_ref().unwrap().baseline;
        assert_eq!(cases.len(), 2);
        assert!(cases[0].compliance > 0.0);
        assert!((cases[0].compliance / cases[1].compliance - 1.0).abs() < 1e-8);
        let expected = if aggregate == "weighted-sum" { cases[0].compliance + cases[1].compliance }
            else { cases[0].compliance.max(cases[1].compliance) };
        assert_eq!(expected.to_bits(), evidence.baseline.compliance.to_bits());
        assert!(evidence.baseline.sampled_max_von_mises > 0.0);
        let limit = 1.5 * cases[0].sampled_max_von_mises;
        let bad = source().replace(":traction-pa (0.0 1.0)", ":traction-pa (0.0 2.0)")
            .replace(":weight 1.0", ":weight 0.0")
            .replace(":sampled-stress-limit-pa 1000000000000.0", &format!(":sampled-stress-limit-pa {limit:.17e}"));
        let error = drive(&study_spec(&bad).unwrap(), &ledger, None, &gate(), None).unwrap_err();
        assert!(error.message.contains("sampled stress limit exceeded"));
        assert!(!ledger.in_transaction());
    }
}

#[test]
fn multi_load_updates_preserve_regions_and_resume_without_refunding_study_work() {
    let text = PROTECTED.replace("    :design-regions (", &format!("{FAMILY}\n    :design-regions ("));
    let spec = study_spec(&text).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = super::super::drive(&spec, &ledger, Some(1), &gate(), None).unwrap();
    assert_eq!(first.status, "budget-exhausted", "{}", first.receipt);
    assert_eq!(integer(&json(&first), "iterations_completed").unwrap(), 1, "non-vacuous accepted update");
    let (phi, evidence) = restored(&ledger, &first, &spec);
    let policy = stress_controls(&spec).unwrap();
    let ControlFlow::Continue(prepared) = regions::prepare(&spec, &policy.regions,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("region preparation interrupted") };
    for &(node, expected) in &prepared.fixed_nodes { assert_eq!(phi.nodes()[node].to_bits(), expected.to_bits()); }
    assert!(evidence.current().compliance < evidence.baseline.compliance);
    let old = load(&ledger, &first.pointer).unwrap();
    let resumed = drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &whole_db, None, &gate(), None).unwrap();
    assert_eq!(whole.status, resumed.status);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(), linked(&whole_db, &json(&whole), key, kind).unwrap());
    }
    let a = json(&resumed); let b = json(&whole);
    for key in ["baseline_cases", "accepted_cases", "checkpoint_hex", "solves_started"] {
        assert_eq!(history(&a).get(key), history(&b).get(key), "{key}");
    }
    assert_eq!(integer(history(&a), "recovery_solves_used").unwrap(), 4);
    assert_eq!(integer(history(&b), "recovery_solves_used").unwrap(), 0);
    assert_eq!(load(&ledger, &first.pointer).unwrap().bytes, first.receipt);
}

#[test]
fn cancellation_in_late_case_keeps_the_field_but_charges_started_solves() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancelled = gate();
    let out = driver::drive_observed(&spec, &ledger, None, &cancelled, None, |stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseIterations { case: 1, iterations, .. } if iterations > 0) {
            cancelled.request();
        }
    }).unwrap();
    assert_eq!(out.status, "cancelled");
    let (phi, evidence) = restored(&ledger, &out, &spec);
    assert!(evidence.accepted.is_empty());
    assert_eq!(snapshot(&phi), evidence.baseline.snapshot);
    let spent = evidence.family.as_ref().unwrap().solves;
    assert_eq!(spent, 4, "two baseline solves and two started candidate solves");
    let old = load(&ledger, &out.pointer).unwrap();
    let retry = drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
    let value = json(&retry);
    assert_eq!(integer(&value, "iterations_completed").unwrap(), 1);
    assert!(integer(history(&value), "solves_started").unwrap() >= spent + 2);
    assert_eq!(integer(history(&value), "recovery_solves_used").unwrap(), 4);
}

#[test]
fn exhausted_solve_and_recovery_allowances_cannot_be_reset_by_resume() {
    let spec = study_spec(&source().replace(":max-solves 512", ":max-solves 2")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, None, &gate(), None).unwrap();
    assert_eq!(out.status, "budget-exhausted");
    assert_eq!(integer(&json(&out), "iterations_completed").unwrap(), 0);
    let old = load(&ledger, &out.pointer).unwrap();
    let repeated = drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap();
    assert_eq!(repeated.receipt, out.receipt);
    let spec = study_spec(&source().replace(":max-recovery-solves 64", ":max-recovery-solves 0")).unwrap();
    let seed = drive(&spec, &ledger, Some(0), &gate(), None).unwrap();
    let old = load(&ledger, &seed.pointer).unwrap();
    let error = drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap_err();
    assert_eq!(error.exit, exit::BUDGET);
    assert!(error.message.contains(&seed.pointer));
}
