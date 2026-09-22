use super::*;
use crate::study::elasticity::parse as parse_study;

const SOURCE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-multi-load-2d.fsim"));
const OPT_IN: &str = "      :stress-restoration-reduction 0.0\n";

fn source() -> String {
    SOURCE.replace("      :max-recovery-solves 64\n", &format!("      :max-recovery-solves 64\n{OPT_IN}"))
}
fn spec() -> ElasticitySpec { parse_study(&source()).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn constraints(out: &Outcome) -> JsonValue {
    json(out).path(&["continuation", "constraints"]).unwrap().clone()
}
fn read(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec)
    -> (GridSdf, OptimizeReport, ConstraintEvidence)
{
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let iterations = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &receipt, &design, &iterations).unwrap();
    let state = ConstraintEvidence::read(receipt.path(&["continuation", "constraints"]).unwrap(),
        &report, spec.projected.as_ref().unwrap()).unwrap();
    (phi, report, state)
}
fn owner(spec: &ElasticitySpec) -> MultiLoadProjectedOptimizer {
    let policy = spec.projected.as_ref().unwrap();
    let family = policy.family.as_ref().unwrap();
    let ControlFlow::Continue(prepared) = regions::prepare(spec, &policy.regions,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("preparation stopped") };
    MultiLoadProjectedOptimizer::new(prepared.geometry, &family.cases(spec).unwrap(), settings(spec, spec.steps),
        family.aggregate, prepared.fixed_nodes, policy.area, family.controls(policy)).unwrap()
        .with_stress_restoration(policy.stress, family.restoration_reduction.unwrap()).unwrap()
}
fn measured(stress: f64, compliance: f64) -> SampledStressEvaluation {
    SampledStressEvaluation { compliance, volume: 0.75, sampled_max_von_mises: stress,
        max_location: [0.5, 0.5], sample_count: 5, snapshot: 42 }
}

#[test]
fn explicit_restoration_policy_is_canonical_and_can_use_only_the_primary_load() {
    let ordinary = parse_study(SOURCE).unwrap();
    assert!(ordinary.projected.as_ref().unwrap().family.as_ref().unwrap().restoration_reduction.is_none());
    let declared = spec();
    assert_eq!(parse_study(&declared.canonical).unwrap().id, declared.id);
    for invalid in ["-0.1", "1.0", "2.0"] {
        assert!(parse_study(&source().replace("stress-restoration-reduction 0.0",
            &format!("stress-restoration-reduction {invalid}"))).is_err());
    }
    assert!(parse_study(&source().replace(OPT_IN, &format!("{OPT_IN}{OPT_IN}"))).is_err());
    let primary_only = source().replace(
        "        (case :band (0.375 0.625) :traction-pa (0.5 0.0) :weight 1.0)\n", "");
    let primary = parse_study(&primary_only).unwrap();
    assert_eq!(primary.projected.as_ref().unwrap().family.as_ref().unwrap().cases(&primary).unwrap().len(), 1);
    assert!(parse_study(&primary_only.replace(OPT_IN, "")).is_err());
}

#[test]
fn restoration_and_compliance_have_distinct_measured_gates_and_progress_baselines() {
    let mut policy = spec().projected.unwrap();
    policy.stress = SampledStressLimit::new(10.0, 0.0).unwrap();
    policy.family.as_mut().unwrap().restoration_reduction = Some(0.5);
    let initial = measured(12.0, 100.0);
    let partial = measured(10.9, 120.0); // worse compliance is permitted ONLY in restoration
    assert!(transition(&initial, &partial, &policy).unwrap());
    assert!(transition(&initial, &measured(11.0, 1.0), &policy).is_err()); // exact excess threshold
    assert!(transition(&initial, &measured(13.0, 1.0), &policy).is_err()); // compliance alone cannot repair stress
    assert_eq!(relative_reduction(&initial, &[partial.clone()], &policy), "null");
    let crossing = measured(10.0, 130.0);
    assert!(transition(&partial, &crossing, &policy).unwrap());
    assert!(transition(&crossing, &measured(9.0, 140.0), &policy).is_err());
    assert!(transition(&crossing, &measured(10.1, 1.0), &policy).is_err());
    let improved = measured(9.9, 117.0);
    assert!(!transition(&crossing, &improved, &policy).unwrap());
    let accepted = [partial, crossing.clone(), improved];
    let text = metadata(&initial, &accepted, &policy).unwrap();
    let status = document(text.as_bytes()).unwrap();
    assert_eq!(integer(&status, "restoration_updates").unwrap(), 2);
    assert_eq!(integer(&status, "compliance_updates").unwrap(), 1);
    assert_eq!(integer(&status, "first_feasible_update").unwrap(), 2);
    assert_eq!(status.get("feasible_baseline"), Some(&document(stress_json(&crossing).as_bytes()).unwrap()));
    assert_eq!(relative_reduction(&initial, &accepted, &policy).parse::<f64>().unwrap(), 0.1);
    let root = document(format!("{{\"stress_restoration\":{text}}}").as_bytes()).unwrap();
    assert!(check_retained(&root, &initial, &accepted, &policy).is_ok());
    let altered = document(format!("{{\"stress_restoration\":{text}}}")
        .replace("\"restoration_updates\":2", "\"restoration_updates\":0").as_bytes()).unwrap();
    assert!(check_retained(&altered, &initial, &accepted, &policy).is_err());
}

#[test]
fn unrepairable_budget_limited_inputs_are_retained_as_infeasible_not_completed() {
    let text = source().replace(":sampled-stress-limit-pa 1000000000000.0", ":sampled-stress-limit-pa 0.000000000001")
        .replace(":max-solves 512", ":max-solves 2")
        .replace(":traction-pa (0.5 0.0) :weight 1.0", ":traction-pa (0.0 -10.0) :weight 0.0");
    let spec = parse_study(&text).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let result = driver::drive(&spec, &ledger, None, &gate(), None).unwrap();
    assert_eq!(result.status, "budget-exhausted");
    let (_, rows, evidence) = read(&ledger, &result, &spec);
    assert!(rows.rows.is_empty());
    let cases = &evidence.family.as_ref().unwrap().baseline;
    assert!(cases[1].sampled_max_von_mises > cases[0].sampled_max_von_mises);
    let declared = constraints(&result);
    assert_eq!(declared.str_field("baseline_scope"), Some("area_feasible_study_start"));
    assert_eq!(declared.path(&["stress_restoration", "phase"]).and_then(JsonValue::as_str), Some("restoring-stress"));
    assert_eq!(declared.get("relative_reduction"), Some(&JsonValue::Null));
    assert!(evidence.html().contains("Stress remains infeasible"));
    let loaded = load(&ledger, &result.pointer).unwrap();
    let again = driver::drive(&spec, &ledger, None, &gate(), Some(&loaded)).unwrap();
    assert_eq!(again.receipt, result.receipt); // exhausted work is not restarted by resume
    let strict = parse_study(&text.replace(OPT_IN, "")).unwrap();
    assert!(driver::drive(&strict, &ledger, None, &gate(), None).is_err());
}

#[test]
fn native_restoration_replays_infeasible_state_and_matches_the_real_owner_search() {
    let pilot = owner(&spec());
    let bound = pilot.current_stress().unwrap().worst_sampled_von_mises * (1.0 - 1e-6);
    let text = source().replace(":sampled-stress-limit-pa 1000000000000.0",
        &format!(":sampled-stress-limit-pa {bound:.17e}"));
    let spec = parse_study(&text).unwrap();
    let mut reference = owner(&spec);
    assert!(reference.is_restoring_stress());
    let mut sampled = false;
    let ControlFlow::Continue(progress) = reference.advance_one_polling(1, |stage| {
        sampled |= matches!(stage, MultiLoadProjectedStage::StressCell { case: 1, .. });
        ControlFlow::<()>::Continue(())
    }).unwrap() else { panic!("reference interrupted") };
    assert!(sampled, "exercise complete candidate mechanics and stress, not only input parsing");
    let ledger = Ledger::open(":memory:").unwrap();
    let cancel = gate();
    let seed = driver::drive_observed(&spec, &ledger, None, &cancel, None, |stage| {
        if matches!(stage, MultiLoadProjectedStage::Direction) { cancel.request(); }
    }).unwrap();
    assert_eq!(seed.status, "cancelled");
    let loaded = load(&ledger, &seed.pointer).unwrap();
    let resumed = driver::drive(&spec, &ledger, Some(1), &gate(), Some(&loaded)).unwrap();
    let (phi, rows, state) = read(&ledger, &resumed, &spec);
    let history = state.family.as_ref().unwrap();
    assert_eq!(history.checkpoint, reference.checkpoint_bytes());
    assert_eq!(history.recovery_solves, 4);
    assert_eq!(phi.nodes(), reference.geometry().nodes());
    assert_eq!(rows.rows.len(), reference.next_iteration());
    match progress {
        MultiLoadProjectedProgress::Accepted(step) => {
            assert!(step.restoration);
            assert_eq!(updates(&state.baseline, &state.accepted, &state.policy), 1);
            assert!(state.current().sampled_max_von_mises < state.baseline.sampled_max_von_mises);
            assert_eq!(resumed.status, "budget-exhausted");
        }
        MultiLoadProjectedProgress::NoDescent(_) => assert_eq!(resumed.status, "no-feasible-descent"),
        _ => panic!("complete explicit search budget should be available"),
    }
    // Restoration may legitimately stall. The native result must be the actual
    // owner result, not fabricated feasibility, a different limit or another PDE.
    let replay = fs_topols::evaluate_robust_sampled_stress(&phi, reference.load_cases(),
        reference.settings(), reference.aggregate()).unwrap();
    assert_eq!(&replay, reference.current_stress().unwrap());
    let fixed = reference.fixed_nodes();
    for &(index, value) in fixed { assert_eq!(phi.nodes()[index].to_bits(), value.to_bits()); }
    assert_eq!(load(&ledger, &seed.pointer).unwrap().bytes, seed.receipt);
}

#[test]
fn cancellation_during_a_late_restoration_case_preserves_phase_and_charges_started_work() {
    let spec = parse_study(&source().replace(":sampled-stress-limit-pa 1000000000000.0",
        ":sampled-stress-limit-pa 0.000000000001")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancel = gate();
    let result = driver::drive_observed(&spec, &ledger, None, &cancel, None, |stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseIterations { case: 1, iterations, .. } if iterations > 0) {
            cancel.request();
        }
    }).unwrap();
    assert_eq!(result.status, "cancelled");
    let (_, rows, state) = read(&ledger, &result, &spec);
    assert!(rows.rows.is_empty());
    assert!(same(&state.baseline, state.current()));
    let history = state.family.as_ref().unwrap();
    assert_eq!(history.solves, 4); // two baseline and two started candidate solves
    let mut allowance = 4;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&history.checkpoint, &mut allowance).unwrap();
    assert!(restored.is_restoring_stress());
    assert_eq!(restored.restoration_updates(), 0);
    assert_eq!(allowance, 0);
    assert_eq!(restored.checkpoint_bytes(), history.checkpoint);
}
