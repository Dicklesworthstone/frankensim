use super::*;
use crate::study::elasticity::parse as parse_study;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-adaptive-dwr-2d.fsim"));

fn spec() -> ElasticitySpec { parse_study(EXAMPLE).unwrap() }
fn receipt(out: &Outcome) -> JsonValue { JsonValue::parse(&out.receipt).unwrap() }
fn summary(ledger: &Ledger, out: &Outcome) -> JsonValue {
    let bytes = linked(ledger, &receipt(out), "report_json", "study-report-json").unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}
fn assessed(value: &JsonValue) -> &JsonValue { value.get("goal_error_assessment").unwrap() }

#[test]
fn native_adaptive_assessment_requires_explicit_complete_bounded_controls() {
    let spec = spec();
    assert!(spec.final_dwr);
    assert_eq!(max_solves(&spec), 4);
    assert_eq!(parse_study(&spec.canonical).unwrap().id, spec.id);
    assert_eq!(spec.adaptive_dwr.unwrap().max_rounds, 2);
    for (from, to) in [(":max-solves-per-attempt 4", ":max-solves-per-attempt 3"),
        (":max-solves-per-attempt 4", ":max-solves-per-attempt 5"),
        (":max-cells 256", ":max-cells 63"), (":max-level 4", ":max-level 2"),
        (":max-rounds 2", ":max-rounds 0"), (":mark-fraction 0.3", ":mark-fraction 0.0"),
        (":relative-tolerance 0.0", ":relative-tolerance 1.0"),
        (":mark-fraction 0.3", ":mark-fraction 0.3 :mark-fraction 0.4"),
        (":max-cells 256", ""), (":memory 536870912 B", ":memory 134217728 B")] {
        assert!(parse_study(&EXAMPLE.replace(from, to)).is_err(), "must refuse {to}");
    }
    let changed = EXAMPLE.replace(":mark-fraction 0.3", ":mark-fraction 0.4");
    assert_ne!(parse_study(&changed).unwrap().id, spec.id);
}

#[test]
fn native_adaptive_assessment_refines_real_endpoint_and_keeps_optimizer_metrics() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = continuation::drive(&spec, &ledger, None, &CancelGate::new_clock_free(), None).unwrap();
    assert_eq!(out.status, "completed", "{}", out.receipt);
    let value = summary(&ledger, &out);
    let a = assessed(&value);
    assert_eq!(a.str_field("method"), Some("adaptive-elasticity-compliance-dwr"));
    let rows = a.get("rounds").and_then(JsonValue::as_array).unwrap();
    assert_eq!(rows.len(), 2, "exercise an actual refinement and new solve");
    assert_eq!(integer(&rows[0], "cells").unwrap(), 64);
    assert!(integer(&rows[1], "cells").unwrap() > 64);
    assert!(integer(&rows[1], "cells").unwrap() < 256, "refinement is local, not uniform");
    assert_eq!(integer(a, "solves").unwrap(), 4);
    assert_eq!(a.get("tolerance_satisfied"), Some(&JsonValue::Bool(false)));
    assert_eq!(a.str_field("stop"), Some("round-budget"));
    assert_eq!(a.str_field("snapshot"), value.str_field("snapshot"));
    assert_eq!(rows[0].f64_field("coarse_compliance_j").unwrap().to_bits(),
        value.f64_field("final_compliance_j").unwrap().to_bits());
    assert_eq!(a.f64_field("original_material_area_m2").unwrap().to_bits(),
        value.f64_field("final_material_area_m2").unwrap().to_bits());
    let cases = a.get("cases").and_then(JsonValue::as_array).unwrap();
    assert_eq!(cases.len(), 1);
    assert!(!cases[0].get("indicators").and_then(JsonValue::as_array).unwrap().is_empty());
    let html = linked(&ledger, &receipt(&out), "report_html", "study-report-html").unwrap();
    assert!(std::str::from_utf8(&html).unwrap().contains("Adaptive compliance assessment"));
    let old = load(&ledger, &out.pointer).unwrap();
    let count = ledger.table_count("ops").unwrap();
    let repeated = continuation::drive(&spec, &ledger, None, &CancelGate::new_clock_free(), Some(&old)).unwrap();
    assert_eq!(repeated.receipt, out.receipt);
    assert_eq!(ledger.table_count("ops").unwrap(), count, "sealed completion must not solve again");
}

fn beam(spec: &ElasticitySpec) -> (GridSdf, Vec<RobustLoadCase>, fs_topols::EvaluatedFinalState) {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let case = RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap();
    let measured = fs_topols::evaluate_compliance_design(&phi, fixture(spec), settings(spec, spec.steps)).unwrap();
    (phi, vec![case], measured)
}
fn design(m: fs_topols::EvaluatedFinalState) -> Design {
    Design { snapshot: m.snapshot, compliance: m.compliance, volume: m.volume }
}

#[test]
fn native_adaptive_assessment_cancellation_and_retry_preserve_the_exact_design() {
    let spec = spec();
    let (phi, cases, m) = beam(&spec);
    let before: Vec<_> = phi.nodes().iter().map(|x| x.to_bits()).collect();
    for stop in [AdaptiveDwrStage::BeforeCase { round: 0, case: 0 },
        AdaptiveDwrStage::RoundComplete(0), AdaptiveDwrStage::BeforeRefine(0),
        AdaptiveDwrStage::BeforeCase { round: 1, case: 0 }, AdaptiveDwrStage::BeforePublish] {
        let (status, result) = run(&spec, &phi, design(m), &cases, &[m.compliance], "completed", |stage| {
            if stage == stop { ControlFlow::Break("cancelled") } else { ControlFlow::Continue(()) }
        });
        assert_eq!(status, "cancelled", "must reach {stop:?}");
        assert!(result.is_none());
        assert_eq!(phi.nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(), before);
    }
    let fresh = || run(&spec, &phi, design(m), &cases, &[m.compliance], "completed", |_| ControlFlow::Continue(()));
    let (status, first) = fresh();
    assert_eq!(status, "completed");
    let (_, second) = fresh();
    assert_eq!(super::super::json(first.as_ref(), None), super::super::json(second.as_ref(), None));
}

#[test]
fn native_adaptive_assessment_rejects_mismatched_retained_response() {
    let spec = spec();
    let (phi, cases, m) = beam(&spec);
    let (status, result) = run(&spec, &phi, design(m), &cases, &[m.compliance * 1.01],
        "completed", |_| ControlFlow::Continue(()));
    assert_eq!(status, "numerical-failure");
    assert!(matches!(result, Some(FinalAssessment::Refused { max_solves: 4, .. })));
    assert_eq!(snapshot(&phi), m.snapshot);
}
