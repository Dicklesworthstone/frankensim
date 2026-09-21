use super::*;
use fs_topols::evaluated::DesignEvaluationStage;

const ORIGINAL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/marquee/bracket-2d.fsim"));
const POLICY: &str = "    :constraint-mode projected-stress\n    :area-tolerance-m2 0.0001\n    :max-projection-shift 2.0\n    :max-area-evaluations 64\n    :max-candidates 16\n    :contraction 0.5\n    :min-relative-improvement 0.00000001\n    :cg-poll-iters 1\n    :sampled-stress-limit-pa 1000000000000.0\n    :stress-tolerance-pa 0.0)";

fn source() -> String {
    ORIGINAL.replace(":mesh-level 4", ":mesh-level 3")
        .replace(":youngs-modulus-pa 70000000000.0", ":youngs-modulus-pa 2.0")
        .replace(":load-traction-pa 1000000.0", ":load-traction-pa 1.0")
        .replace(":volume-fraction 0.45", ":volume-fraction 0.75")
        .replace(":max-iterations 8", ":max-iterations 2")
        .replace(":move-cells 0.35", ":move-cells 0.1")
        .replace(":nucleation-period 4", ":nucleation-period 0")
        .replace("    :steps 8)", &format!("    :steps 2\n{POLICY}"))
}
fn spec() -> ElasticitySpec { parse(&source()).expect("explicit projected native study") }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(output: &Outcome) -> JsonValue { document(output.receipt.as_bytes()).unwrap() }
fn constrained(value: &JsonValue) -> &JsonValue {
    value.path(&["continuation", "constraints"]).expect("retained constraints")
}
fn read(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, OptimizeReport, ConstraintEvidence) {
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &receipt, &design, &rows).unwrap();
    let state = ConstraintEvidence::read(constrained(&receipt), &report, spec.projected.unwrap()).unwrap();
    (phi, report, state)
}

#[test]
fn projected_study_controls_are_explicit_and_legacy_canonical_bytes_are_unchanged() {
    let legacy = parse(ORIGINAL).unwrap();
    assert!(legacy.projected.is_none());
    assert!(!legacy.canonical.contains("constraint-mode"));
    let declared = spec();
    let reparsed = parse(&declared.canonical).unwrap();
    assert_eq!(declared.canonical, reparsed.canonical);
    assert_eq!(declared.id, reparsed.id);
    for field in POLICY.lines().skip(1) {
        // Remove one required key/value but preserve the closing parenthesis.
        let key_value = field.trim_end_matches(')');
        assert!(parse(&source().replace(key_value, "")).is_err(), "missing {key_value}");
    }
    for (from, to) in [
        (":max-candidates 16", ":max-candidates 0"),
        (":cg-poll-iters 1", ":cg-poll-iters 0"),
        (":contraction 0.5", ":contraction 1.0"),
        (":area-tolerance-m2 0.0001", ":area-tolerance-m2 0.5"),
        (":stress-tolerance-pa 0.0", ":stress-tolerance-pa -1.0"),
        (":constraint-mode projected-stress", ":constraint-mode projected-stress :max-candidates 16"),
    ] { assert!(parse(&source().replace(from, to)).is_err(), "{to}"); }
    assert!(parse(&source().replace(":constraint-mode projected-stress", "")).is_err());
}

#[test]
fn projected_study_accepts_same_area_stress_feasible_design_and_resumes_exactly() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = super::super::drive(&spec, &ledger, Some(1), &gate(), None).unwrap();
    assert_eq!(first.status, "budget-exhausted", "{}", first.receipt);
    let (phi, rows, evidence) = read(&ledger, &first, &spec);
    assert_eq!(rows.rows.len(), 1, "the native fixture must produce a non-vacuous accepted update");
    assert!(rows.compliance[0] < evidence.baseline.compliance);
    assert!((rows.volume[0] - 0.75).abs() <= 1e-4);
    let oracle = fs_topols::evaluate_sampled_stress(&phi, fixture(&spec), settings(&spec, spec.steps)).unwrap();
    assert!(same(&oracle, evidence.current()));
    let prior_bytes = linked(&ledger, &json(&first), "design", "study-design").unwrap();
    let old = load(&ledger, &first.pointer).unwrap();
    let resumed = super::super::drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap();
    let (_, _, resumed_evidence) = read(&ledger, &resumed, &spec);
    let full_ledger = Ledger::open(":memory:").unwrap();
    let full = super::super::drive(&spec, &full_ledger, None, &gate(), None).unwrap();
    assert_eq!(full.status, resumed.status);
    let (_, _, full_evidence) = read(&full_ledger, &full, &spec);
    assert_eq!(resumed_evidence.json(), full_evidence.json());
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(),
            linked(&full_ledger, &json(&full), key, kind).unwrap());
    }
    assert_eq!(linked(&ledger, &json(&first), "design", "study-design").unwrap(), prior_bytes);
    assert_eq!(json(&resumed).path(&["continuation", "legacy_prefix_updates_replayed"])
        .and_then(JsonValue::as_f64), Some(0.0));
}

#[test]
fn projected_study_cancelled_final_solve_retains_feasible_seed_and_retry_is_exact() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancelled = gate();
    let out = drive_observed(&spec, &ledger, None, &cancelled, None, |stage| {
        if matches!(stage, Stage::Update(ProjectedStage::Evaluation(_, DesignEvaluationStage::Solve(n))) if n > 0) {
            cancelled.request();
        }
    }).unwrap();
    assert_eq!(out.status, "cancelled");
    let (phi, rows, state) = read(&ledger, &out, &spec);
    assert!(rows.rows.is_empty());
    assert_eq!(snapshot(&phi), state.baseline.snapshot);
    assert_ne!(snapshot(&phi), snapshot(&initial_phi(&spec)), "baseline must be projected before claiming feasibility");
    let old = load(&ledger, &out.pointer).unwrap();
    let retry = drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
    assert_eq!(integer(&json(&retry), "iterations_completed").unwrap(), 1);
    let summary = linked(&ledger, &json(&out), "report_json", "study-report-json").unwrap();
    let summary = document(&summary).unwrap();
    assert_eq!(summary.f64_field("final_compliance_j").unwrap().to_bits(), state.baseline.compliance.to_bits());
    assert!(summary.get("constraints").is_some());
}

#[test]
fn projected_study_infeasible_baseline_and_setup_cancellation_are_not_successes() {
    let bad = parse(&source().replace(":sampled-stress-limit-pa 1000000000000.0", ":sampled-stress-limit-pa 0.000000000001")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let error = drive(&bad, &ledger, None, &gate(), None).unwrap_err();
    assert!(error.message.contains("baseline refused"));
    assert!(!ledger.in_transaction());
    let cancelled = gate();
    let error = drive_observed(&spec(), &ledger, None, &cancelled, None, |stage| {
        if matches!(stage, Stage::Setup(ProjectedStressSetupStage::Stress(DesignEvaluationStage::StressCell(_)))) {
            cancelled.request();
        }
    }).unwrap_err();
    assert_eq!(error.exit, exit::CANCELLED);
    assert!(error.message.contains("no new feasible state published"));
}

#[test]
fn projected_study_no_descent_retains_refusals_and_never_reports_completion() {
    let spec = parse(&source().replace(":min-relative-improvement 0.00000001", ":min-relative-improvement 0.999999")
        .replace(":max-candidates 16", ":max-candidates 2")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, None, &gate(), None).unwrap();
    assert_eq!(out.status, "no-feasible-descent");
    let (phi, rows, evidence) = read(&ledger, &out, &spec);
    assert!(rows.rows.is_empty());
    assert_eq!(snapshot(&phi), evidence.baseline.snapshot);
    assert_eq!(evidence.refusals.len(), 2);
    let rendered = render(out, OutputMode::Json);
    assert_eq!(rendered.exit_code, exit::REFUSED);
    assert!(rendered.stdout.contains("no-feasible-descent"));
}
