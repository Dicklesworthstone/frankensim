use super::*;
use fs_uq::{CorrelationModel, ParameterUncertainty};

fn plan(n: usize) -> UqPlan {
    UqPlan::new("binary-event", PropagationMethod::MonteCarlo, n)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("x", 0.0, 1.0, "1"))
        .with_compliance_threshold(0.5)
}
fn policy(probability: f64) -> Policy {
    let mut policy = Policy::new(probability, 0.05).unwrap();
    policy.method = Method::BernoulliMixture;
    policy
}

#[test]
fn both_rare_event_directions_use_the_selected_owner_and_retain_actual_work() {
    for (value, p, decision) in [(0.0, 0.99, Decision::Satisfied), (1.0, 0.01, Decision::Violated)] {
        let mut calls = 0;
        let result = run(&plan(2048), Some(policy(p)), || false, |_| {
            calls += 1; Ok::<_, &str>(value)
        }).unwrap();
        let assessment = result.assessment.unwrap();
        assert_eq!(assessment.decision, decision);
        assert_eq!(calls, 1024); assert_eq!(result.result.samples_evaluated, calls);
        assert_eq!(result.result.status, UqStatus::BudgetTruncated);
        assert_eq!(assessment.policy.method, Method::BernoulliMixture);
        let mut direct = UqExecution::new(&plan(2048)).unwrap();
        direct.advance(calls, || false, |_| Ok::<_, &str>(value));
        let cs = direct.assess_bernoulli_compliance(0.05, 0.0).unwrap().unwrap();
        assert_eq!((assessment.lower, assessment.upper), (cs.lo, cs.hi));
        assert!(assessment.lower < assessment.upper);
    }
    let old = run(&plan(2048), Some(Policy::new(0.99, 0.05).unwrap()), || false,
        |_| Ok::<_, &str>(0.0)).unwrap();
    assert_eq!(old.result.samples_evaluated, 2048);
    assert_eq!(old.assessment.unwrap().decision, Decision::Inconclusive);
}

#[test]
fn failures_and_midrun_cancellation_cannot_reuse_a_previous_bernoulli_bound() {
    for terminal in [Err("solver refusal"), Ok(f64::NAN)] {
        let mut calls = 0;
        let result = run(&plan(2048), Some(policy(0.99)), || false, |_| {
            calls += 1; if calls == 33 { terminal } else { Ok(0.0) }
        }).unwrap();
        assert_eq!(calls, 33); assert_eq!(result.result.status, UqStatus::Refused);
        assert!(result.assessment.is_none());
    }
    let calls = std::cell::Cell::new(0);
    let result = run(&plan(2048), Some(policy(0.99)), || calls.get() == 33, |_| {
        calls.set(calls.get()+1); Ok::<_, &str>(0.0)
    }).unwrap();
    assert_eq!(result.result.samples_evaluated, 33);
    assert_eq!(result.result.status, UqStatus::Cancelled); assert!(result.assessment.is_none());
}

#[test]
fn bernoulli_budget_is_not_promoted_to_an_unsupported_probability_pass() {
    let run = run(&plan(19), Some(policy(0.999)), || false, |_| Ok::<_, &str>(0.0)).unwrap();
    assert_eq!(run.result.status, UqStatus::Complete);
    assert_eq!(run.result.samples_evaluated, 19);
    let a = run.assessment.unwrap();
    assert_eq!(a.decision, Decision::Inconclusive);
    assert!(a.lower < 0.999 && a.upper == 1.0);
}
