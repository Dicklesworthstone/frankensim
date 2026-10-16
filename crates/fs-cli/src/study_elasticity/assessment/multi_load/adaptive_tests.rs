use super::*;
use crate::study::elasticity::parse as parse_study;
use fs_topols::RobustAggregate;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-weighted-adaptive-dwr-2d.fsim"));

fn spec() -> ElasticitySpec { parse_study(EXAMPLE).unwrap() }
fn report(ledger: &Ledger, out: &Outcome) -> JsonValue {
    let receipt = JsonValue::parse(&out.receipt).unwrap();
    let bytes = linked(ledger, &receipt, "report_json", "study-report-json").unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

#[test]
fn native_adaptive_weighted_requires_every_declared_case_in_the_round_allowance() {
    let spec = spec();
    assert_eq!(max_solves(&spec), 8);
    assert_eq!(parse_study(&spec.canonical).unwrap().id, spec.id);
    for grant in [4, 7, 9] {
        assert!(parse_study(&EXAMPLE.replace(":max-solves-per-attempt 8",
            &format!(":max-solves-per-attempt {grant}"))).is_err());
    }
    assert!(parse_study(&EXAMPLE.replace(":aggregate weighted-sum",
        ":aggregate worst-weighted-case")).is_err());
    assert!(parse_study(&EXAMPLE.replace(":memory 536870912 B", ":memory 134217728 B")).is_err());
}

#[test]
fn native_adaptive_weighted_refines_all_loads_and_retains_zero_weight_response() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = continuation::drive(&spec, &ledger, None, &CancelGate::new_clock_free(), None).unwrap();
    assert_eq!(out.status, "completed", "{}", out.receipt);
    let value = report(&ledger, &out);
    let a = value.get("goal_error_assessment").unwrap();
    assert_eq!(a.str_field("method"), Some("adaptive-elasticity-compliance-dwr"));
    assert_eq!(integer(a, "solves").unwrap(), 8);
    let rows = a.get("rounds").and_then(JsonValue::as_array).unwrap();
    assert_eq!(rows.len(), 2);
    assert!(integer(&rows[1], "cells").unwrap() > integer(&rows[0], "cells").unwrap());
    let cases = a.get("cases").and_then(JsonValue::as_array).unwrap();
    assert_eq!(cases.len(), 2);
    assert_eq!(cases[1].f64_field("weight"), Some(0.0));
    for key in ["initial_compliance_j", "coarse_compliance_j", "enriched_compliance_j"] {
        let first = cases[0].f64_field(key).unwrap();
        let second = cases[1].f64_field(key).unwrap();
        assert!(first > 0.0 && (second / first - 9.0).abs() < 1e-7, "independent response {key}");
    }
    assert!(!cases[1].get("indicators").and_then(JsonValue::as_array).unwrap().is_empty());
    assert_eq!(cases[0].f64_field("initial_compliance_j").unwrap().to_bits(),
        value.f64_field("final_compliance_j").unwrap().to_bits());
    assert_eq!(a.str_field("snapshot"), value.str_field("snapshot"));
    let old = load(&ledger, &out.pointer).unwrap();
    let ops = ledger.table_count("ops").unwrap();
    let again = continuation::drive(&spec, &ledger, None, &CancelGate::new_clock_free(), Some(&old)).unwrap();
    assert_eq!(again.receipt, out.receipt);
    assert_eq!(ledger.table_count("ops").unwrap(), ops);
}

fn beam(spec: &ElasticitySpec) -> (GridSdf, fs_topols::RobustSampledStressEvaluation,
    Vec<SampledStressEvaluation>) {
    let Some(continuation::ProjectedControls::Stress(policy)) = &spec.projected else { panic!("stress policy") };
    let cases = policy.dwr_cases(spec).unwrap();
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let m = fs_topols::evaluate_robust_sampled_stress(&phi, &cases, settings(spec, spec.steps),
        RobustAggregate::WeightedSum).unwrap();
    let expected = (0..cases.len()).map(|i| SampledStressEvaluation {
        compliance: m.case_compliances[i], volume: m.volume, snapshot: m.snapshot,
        sampled_max_von_mises: m.case_sampled_max_von_mises[i],
        max_location: m.case_max_locations[i], sample_count: m.case_sample_counts[i],
    }).collect();
    (phi, m, expected)
}
fn design(m: &fs_topols::RobustSampledStressEvaluation) -> Design {
    Design { snapshot: m.snapshot, compliance: m.objective, volume: m.volume }
}

#[test]
fn native_adaptive_weighted_equal_aggregate_cannot_hide_changed_zero_weight_case() {
    let spec = spec();
    let (phi, m, mut expected) = beam(&spec);
    // This leaves the weighted objective, geometry identity and area unchanged.
    expected[1].compliance *= 1.01;
    let (status, result) = run(&spec, &phi, design(&m), &expected, "completed", |_| ControlFlow::Continue(()));
    assert_eq!(status, "numerical-failure");
    assert!(matches!(result, Some(FinalAssessment::Refused { max_solves: 8, .. })));
    assert_eq!(snapshot(&phi), m.snapshot);
}

#[test]
fn native_adaptive_weighted_late_case_cancellation_discards_partial_analysis_and_retries() {
    let spec = spec();
    let (phi, m, expected) = beam(&spec);
    let before: Vec<_> = phi.nodes().iter().map(|x| x.to_bits()).collect();
    let mut late_case_visits = 0;
    let (status, result) = run(&spec, &phi, design(&m), &expected, "completed", |stage| {
        if matches!(stage, WeightedComplianceDwrStage::BeforeCase { case: 1 }) {
            late_case_visits += 1;
            if late_case_visits == 2 { return ControlFlow::Break("cancelled"); }
        }
        ControlFlow::Continue(())
    });
    assert_eq!(late_case_visits, 2);
    assert_eq!(status, "cancelled");
    assert!(result.is_none());
    assert_eq!(phi.nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(), before);
    let fresh = || run(&spec, &phi, design(&m), &expected, "completed", |_| ControlFlow::Continue(()));
    let (status, first) = fresh();
    assert_eq!(status, "completed");
    let (_, second) = fresh();
    assert_eq!(super::super::json(first.as_ref(), None), super::super::json(second.as_ref(), None));
}
