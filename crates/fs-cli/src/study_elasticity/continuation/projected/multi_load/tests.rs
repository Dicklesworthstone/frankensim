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

fn assessed_source(text: &str) -> String {
    format!("{}\n  (assessment :type elasticity-dwr :max-solves-per-attempt 4)\n)\n",
        text.trim_end().strip_suffix(')').unwrap())
}

#[test]
fn weighted_dwr_requires_exact_family_solve_grant_and_rejects_nonsmooth_objectives() {
    let text = assessed_source(&source());
    let spec = study_spec(&text).unwrap();
    assert!(spec.final_dwr);
    assert_eq!(study_spec(&spec.canonical).unwrap().id, spec.id);
    for grant in [2, 3, 6, 33] {
        assert!(study_spec(&text.replace(":max-solves-per-attempt 4",
            &format!(":max-solves-per-attempt {grant}"))).is_err());
    }
    assert!(study_spec(&text.replace("weighted-sum", "worst-weighted-case")).is_err());
    assert!(!study_spec(&source()).unwrap().final_dwr);
}

#[test]
fn final_weighted_dwr_resumes_only_assessment_of_complete_or_zero_update_endpoints() {
    use fs_topols::WeightedComplianceDwrStage;
    for (stalled, boundary) in [(false, WeightedComplianceDwrStage::BeforeCase { case: 1 }),
        (true, WeightedComplianceDwrStage::BeforePublish)] {
        let text = source().replace(":steps 2", ":steps 1").replace(":max-iterations 2", ":max-iterations 1")
            .replace(":max-recovery-solves 64", ":max-recovery-solves 0");
        let text = if stalled {
            text.replace(":min-relative-improvement 0.00000001", ":min-relative-improvement 0.999999")
                .replace(":max-candidates 16", ":max-candidates 2")
        } else { text };
        let spec = study_spec(&assessed_source(&text)).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let cancel = gate();
        let mut reached = false;
        let interrupted = driver::drive_assessment_observed(&spec, &ledger, None, &cancel, None, |_| {}, |stage| {
            if stage == boundary { reached = true; cancel.request(); }
        }).unwrap();
        assert!(reached, "real case solves must reach the requested phase");
        assert_eq!(interrupted.status, "cancelled");
        let terminal = if stalled { "no-feasible-descent" } else { "completed" };
        let interrupted_json = json(&interrupted);
        assert_eq!(integer(&interrupted_json, "iterations_completed").unwrap(), if stalled { 0 } else { 1 });
        assert_eq!(history(&interrupted_json).str_field("optimizer_terminal"), Some(terminal));
        let charged = interrupted_json.f64_field("consumed_wall_s").unwrap();
        let pending = document(&linked(&ledger, &interrupted_json, "report_json", "study-report-json").unwrap()).unwrap();
        assert_eq!(pending.path(&["goal_error_assessment", "status"]).and_then(JsonValue::as_str), Some("pending"));
        assert_eq!(pending.path(&["goal_error_assessment", "max_solves_per_attempt"]).and_then(JsonValue::as_f64), Some(4.0));
        let old = load(&ledger, &interrupted.pointer).unwrap();
        let resumed = driver::drive_assessment_observed(&spec, &ledger, None, &gate(), Some(&old),
            |stage| panic!("assessment retry must not restart optimization: {stage:?}"), |_| {}).unwrap();
        assert_eq!(resumed.status, terminal);
        let receipt = json(&resumed);
        assert!(receipt.f64_field("consumed_wall_s").unwrap() >= charged);
        assert_eq!(history(&receipt), history(&interrupted_json), "no study/recovery work refunded or spent");
        assert_eq!(integer(history(&receipt), "recovery_solves_used").unwrap(), 0);
        for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
            assert_eq!(linked(&ledger, &receipt, key, kind).unwrap(), linked(&ledger, &interrupted_json, key, kind).unwrap());
        }
        let (_, retained) = restored(&ledger, &resumed, &spec);
        let goal = receipt.path(&["continuation", "goal_error_assessment"]).unwrap();
        assert_eq!(goal.str_field("status"), Some("estimated"));
        assert_eq!(goal.str_field("aggregate"), Some("weighted-sum"));
        assert_eq!(goal.f64_field("coarse_compliance_j").unwrap().to_bits(), retained.current().compliance.to_bits());
        let estimates = goal.get("cases").and_then(JsonValue::as_array).unwrap();
        let family = retained.family.as_ref().unwrap();
        let expected = family.accepted.last().unwrap_or(&family.baseline);
        assert_eq!(estimates.len(), expected.len());
        for (actual, expected) in estimates.iter().zip(expected) {
            assert_eq!(actual.f64_field("coarse_compliance_j").unwrap().to_bits(), expected.compliance.to_bits());
            assert!(actual.f64_field("coarse_relative_residual").unwrap() < 1e-12);
            assert!(actual.f64_field("enriched_relative_residual").unwrap() < 1e-12);
        }
        assert_eq!(goal.path(&["solver", "solves"]).and_then(JsonValue::as_f64), Some(4.0));
        let old = load(&ledger, &resumed.pointer).unwrap();
        let again = driver::drive_assessment_observed(&spec, &ledger, None, &gate(), Some(&old),
            |_| panic!("completed optimizer reused"), |_| panic!("completed DWR reused")).unwrap();
        assert_eq!(again.receipt, resumed.receipt);
        if !stalled {
            // Recover the real receipt immediately after the final update,
            // before the terminal marker was persisted (a possible crash gap).
            let pending = load(&ledger, &format!("study-{}",
                interrupted_json.str_field("predecessor").unwrap())).unwrap();
            let endpoint = load(&ledger, &format!("study-{}",
                pending.value.str_field("predecessor").unwrap())).unwrap();
            assert_eq!(endpoint.value.str_field("status"), Some("running"));
            assert!(history(&endpoint.value).get("optimizer_terminal").is_none());
            let recovered = driver::drive_assessment_observed(&spec, &ledger, None, &gate(), Some(&endpoint),
                |_| panic!("a fully updated feasible endpoint needs no optimizer recovery"), |_| {}).unwrap();
            assert_eq!(recovered.status, "completed");
            assert_eq!(history(&json(&recovered)), history(&receipt));
        }
    }
}

#[test]
fn unfinished_or_infeasible_restoration_never_runs_final_compliance_assessment() {
    let text = source().replace("      :max-recovery-solves 64\n",
        "      :max-recovery-solves 64\n      :stress-restoration-reduction 0.999999\n")
        .replace(":sampled-stress-limit-pa 1000000000000.0", ":sampled-stress-limit-pa 0.000000000001")
        .replace(":max-candidates 16", ":max-candidates 1");
    for exhausted in [true, false] {
        let text = if exhausted { text.replace(":max-solves 512", ":max-solves 2") } else { text.clone() };
        let spec = study_spec(&assessed_source(&text)).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let out = driver::drive_assessment_observed(&spec, &ledger, None, &gate(), None,
            |_| {}, |_| panic!("infeasible or unfinished endpoint must not be assessed")).unwrap();
        assert_eq!(out.status, if exhausted { "budget-exhausted" } else { "no-feasible-descent" });
        let (_, retained) = restored(&ledger, &out, &spec);
        assert!(!restoration::feasible(retained.current(), stress_controls(&spec).unwrap()));
        assert!(retained.family.as_ref().unwrap().optimizer_terminal.is_none());
        if !exhausted {
            assert_eq!(json(&out).path(&["continuation", "goal_error_assessment", "status"])
                .and_then(JsonValue::as_str), Some("refused"));
        }
    }
}
