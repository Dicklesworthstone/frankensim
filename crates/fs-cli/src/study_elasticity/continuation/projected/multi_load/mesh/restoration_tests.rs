//! The native reader and restarted owner use the same all-grid restoration gate.
use super::*;

fn restoring_source(strict: bool) -> String {
    source(strict).replace(":max-recovery-solves 64",
        ":max-recovery-solves 64 :stress-restoration-reduction 0.0")
        .replace(":sampled-stress-limit-pa 1000000000000.0", ":sampled-stress-limit-pa 0.00000001")
}

#[test]
fn native_mesh_restoration_keeps_infeasible_budget_terminals_and_numerical_refusals() {
    for strict in [false, true] {
        let text = restoring_source(strict).replace(":max-solves 512", ":max-solves 5");
        let spec = parse_study(&text).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let out = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), None,
            |_| panic!("baseline only"), |_| panic!("no final assessment"), |_| {}).unwrap();
        assert_eq!(out.status, if strict { "mesh-unresolved" } else { "budget-exhausted" });
        let (_, state) = history(&ledger, &out, &spec);
        let h = state.family.as_ref().unwrap();
        assert_eq!(h.solves, 4);
        assert!(state.accepted.is_empty());
        assert!(state.current().sampled_max_von_mises > state.policy.stress.admitted_max());
        let receipt = json(&out);
        assert_eq!(receipt.path(&["continuation", "constraints", "stress_restoration", "phase"])
            .and_then(JsonValue::as_str), Some("restoring-stress"));
        assert_eq!(receipt.path(&["continuation", "constraints", "relative_reduction"]), Some(&JsonValue::Null));
        let old = load(&ledger, &out.pointer).unwrap();
        let again = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), Some(&old),
            |_| panic!("terminal optimizer is inert"), |_| panic!("no final assessment"),
            |_| panic!("terminal mesh evidence is reused")).unwrap();
        assert_eq!(again.receipt, out.receipt);
    }
}

#[test]
fn native_mesh_restoration_cancel_resume_matches_the_real_checkpoint_owner() {
    let spec = parse_study(&restoring_source(false)).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let seed = driver::drive_mesh_observed(&spec, &ledger, None, &stop, None,
        |_| {}, |_| {}, |stage| {
            if matches!(stage, MultiLoadMeshStage::Optimizer(MultiLoadProjectedStage::Direction)) {
                stop.request();
            }
        }).unwrap();
    assert_eq!(seed.status, "cancelled");
    let (_, state) = history(&ledger, &seed, &spec);
    let h = state.family.as_ref().unwrap();
    assert_eq!(h.solves, 4);
    assert!(state.accepted.is_empty());
    let mut recovery = 2 * h.cases.len();
    let mut reference = MultiLoadProjectedOptimizer::restore_checkpoint(&h.checkpoint, &mut recovery).unwrap();
    let mut sampled = false;
    let ControlFlow::Continue(progress) = reference.advance_one_resolution_polling(
        state.policy.resolution.unwrap(), state.policy.search.poll_iters, |stage| {
            sampled |= matches!(stage, MultiLoadMeshStage::Optimizer(MultiLoadProjectedStage::StressCell { case: 1, .. }));
            ControlFlow::<()>::Continue(())
        }).unwrap() else { panic!("reference interrupted") };
    assert!(sampled, "exercise the real bounded restoration search, not only baseline admission");
    let old = load(&ledger, &seed.pointer).unwrap();
    let resumed = driver::drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
    let (phi, retained) = history(&ledger, &resumed, &spec);
    let actual = retained.family.as_ref().unwrap();
    assert_eq!(phi.nodes(), reference.geometry().nodes());
    assert_eq!(actual.checkpoint, reference.checkpoint_bytes());
    assert_eq!(actual.recovery_solves, 4);
    let (_, checked) = LastCheck::capture(progress).unwrap();
    assert_eq!(actual.mesh.as_ref().unwrap().json(), checked.as_ref().unwrap().json());
    assert_eq!(retained.accepted.len(), reference.next_iteration());
    if reference.next_iteration() > 0 {
        assert!(retained.current().sampled_max_von_mises < retained.baseline.sampled_max_von_mises);
    } else {
        assert_eq!(resumed.status, "no-feasible-descent");
    }
    assert_eq!(load(&ledger, &seed.pointer).unwrap().bytes, seed.receipt);
}

fn measurements(stresses: [[f64; 2]; 2], compliance: f64) -> MultiLoadResolutionReport {
    MultiLoadResolutionReport { rungs: stresses.iter().enumerate().map(|(i, grid)|
        MultiLoadResolutionRung { level: 3 + i as u32, cases: grid.iter().map(|&stress|
            SampledStressEvaluation { compliance, volume: 0.75,
                sampled_max_von_mises: stress, max_location: [0.5, 0.5], sample_count: 5,
                snapshot: 42 }).collect() }).collect() }
}
fn retained_check(check: &LastCheck, policy: &Controls, cases: &[RobustLoadCase]) -> Result<Option<LastCheck>> {
    let value = document(format!("{{{}}}", field(policy.resolution, Some(check)).trim_start_matches(',')).as_bytes()).unwrap();
    let accepted = vec![check.candidate.as_ref().unwrap().rungs[0].cases.clone()];
    super::super::read(&value, policy, cases, &check.baseline.rungs[0].cases, &accepted)
}

#[test]
fn native_retained_restoration_rejects_fine_only_regressions_and_premature_exit() {
    // Controlled public measurements isolate the retained-reader gate; these are
    // not substituted for the independent solves in the native execution tests.
    let spec = parse_study(&restoring_source(false)).unwrap();
    let mut policy = stress_controls(&spec).unwrap().clone();
    policy.stress = SampledStressLimit::new(10.0, 0.0).unwrap();
    policy.family.as_mut().unwrap().restoration_reduction = Some(0.1);
    let cases = policy.family.as_ref().unwrap().cases(&spec).unwrap();
    let mut check = LastCheck { outcome: "accepted",
        baseline: measurements([[2.0, 20.0], [4.0, 30.0]], 1.0),
        candidate: Some(measurements([[3.0, 17.0], [3.0, 24.0]], 2.0)),
        checked: 1, reason: None };
    assert!(retained_check(&check, &policy, &cases).is_ok()); // higher compliance is valid restoration
    check.candidate.as_mut().unwrap().rungs[1].cases[1].sampled_max_von_mises = 29.0;
    assert!(retained_check(&check, &policy, &cases).is_err()); // zero-weight fine case still counts
    check.candidate = Some(measurements([[3.0, 9.0], [3.0, 20.0]], 2.0));
    assert!(retained_check(&check, &policy, &cases).is_err()); // coarse-only exit
    check.candidate.as_mut().unwrap().rungs[1].cases[1].sampled_max_von_mises = 10.0;
    assert!(retained_check(&check, &policy, &cases).is_ok());
    policy.family.as_mut().unwrap().restoration_reduction = None;
    assert!(retained_check(&check, &policy, &cases).is_err()); // no undeclared restoration
}

#[test]
fn native_mesh_restoration_opt_in_preserves_already_feasible_compliance_descent() {
    let ordinary = parse_study(&source(false)).unwrap();
    let restoring = parse_study(&source(false).replace(":max-recovery-solves 64",
        ":max-recovery-solves 64 :stress-restoration-reduction 0.0")).unwrap();
    let db = Ledger::open(":memory:").unwrap();
    let a = driver::drive(&ordinary, &db, Some(1), &gate(), None).unwrap();
    let b = driver::drive(&restoring, &db, Some(1), &gate(), None).unwrap();
    let (phi_a, state_a) = history(&db, &a, &ordinary);
    let (phi_b, state_b) = history(&db, &b, &restoring);
    assert_eq!(state_a.accepted.len(), 1);
    assert_eq!(phi_a.nodes(), phi_b.nodes());
    assert_eq!(state_a.accepted, state_b.accepted);
    assert_eq!(restoration::updates(&state_b.baseline, &state_b.accepted, &state_b.policy), 0);
}
