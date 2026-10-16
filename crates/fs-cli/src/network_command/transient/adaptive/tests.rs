use super::*;
use super::super::super::tests::{close, with_cx};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/adaptive-contact-pulse.json"));
fn field(doc: &J) -> Vec<f64> {
    doc.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}
fn config() -> Config {
    Config { absolute_tolerance_k: 0.1, relative_tolerance: 0.01, minimum_trial_step_s: 1e-5, max_trials: 4000 }
}

#[test]
fn full_field_estimator_scales_temperature_changes_not_kelvin_offsets() {
    with_cx(|cx| {
        let a = error_ratio(cx, &[300.0,300.0], &[301.0,310.0], &[301.0,310.5], config()).unwrap();
        let b = error_ratio(cx, &[400.0,400.0], &[401.0,410.0], &[401.0,410.5], config()).unwrap();
        close(a, 0.5/0.205, 1e-12);
        assert_eq!(a, b);
        assert!(a > 1.0, "the monitored first vertex is exact, but the full field must reject");
        assert!(error_ratio(cx, &[], &[], &[], config()).is_err());
        assert!(error_ratio(cx, &[300.0], &[300.0,301.0], &[300.0], config()).is_err());
        assert!(error_ratio(cx, &[300.0], &[f64::NAN], &[300.0], config()).is_err());
    });
}

#[test]
fn adaptive_contact_pulse_rejects_trials_without_adding_time_or_heat() {
    let r = Request::parse(FIXTURE).unwrap();
    let doc = J::parse(&execute(&r, &CancelGate::new_clock_free()).unwrap()).unwrap();
    let run = doc.get("transient").unwrap();
    let control = run.get("adaptive").unwrap();
    let trials = control.f64_field("trials").unwrap();
    let rejected = control.f64_field("rejected_trials").unwrap();
    assert!(rejected > 0.0);
    assert!(control.f64_field("largest_accepted_error_ratio").unwrap() <= 1.0);
    close(run.f64_field("steps").unwrap(), 2.0*(trials-rejected), 0.0);
    close(run.f64_field("input_energy_j").unwrap(), 600.0, 1e-7);
    close(run.f64_field("time_s").unwrap(), 150.0, 0.0);
    close(run.f64_field("sampled_peak_time_s").unwrap(), 30.0, 0.0);
    close(run.f64_field("sampled_peak_objective_k").unwrap(), 304.11908292812404, 5e-4);
    close(doc.path(&["objective","value_k"]).unwrap().as_f64().unwrap(), 301.74601435038426, 5e-4);
    assert!(run.f64_field("energy_residual_j").unwrap().abs() < r.limits.heat*150.0);
    let history = run.get("history").unwrap().as_array().unwrap();
    let mut previous = 0.0;
    let mut accepted_work = 0.0;
    let mut estimates = 0;
    for sample in &history[1..] {
        let time = sample.f64_field("time_s").unwrap();
        assert!(time > previous);
        assert!(!(previous < 30.0 && time > 30.0));
        close(sample.f64_field("source_w").unwrap(), if time <= 30.0 {20.0} else {0.0}, 1e-8);
        accepted_work += sample.f64_field("coupling_iterations").unwrap();
        if let Some(ratio) = sample.f64_field("estimated_local_error_ratio") { assert!(ratio <= 1.0); estimates += 1; }
        previous = time;
    }
    close(estimates as f64, trials-rejected, 0.0);
    assert!(run.f64_field("total_solid_solves").unwrap() > accepted_work, "discarded work must still count");
}

#[test]
fn tighter_local_tolerance_improves_this_resolved_reference_without_claiming_a_bound() {
    let mut loose = Request::parse(FIXTURE).unwrap();
    loose.transient.as_mut().unwrap().adaptive.as_mut().unwrap().relative_tolerance = 0.0;
    let loose_doc = J::parse(&execute(&loose, &CancelGate::new_clock_free()).unwrap()).unwrap();
    let mut tight = Request::parse(FIXTURE).unwrap();
    let control = tight.transient.as_mut().unwrap().adaptive.as_mut().unwrap();
    control.absolute_tolerance_k = 0.001; control.relative_tolerance = 0.0;
    let tight_doc = J::parse(&execute(&tight, &CancelGate::new_clock_free()).unwrap()).unwrap();
    let base = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/transient-contact-pulse.json"));
    let reference = Request::parse(&base.replace("\"max_step_s\": 2", "\"max_step_s\": 0.25")).unwrap();
    let reference = field(&J::parse(&execute(&reference, &CancelGate::new_clock_free()).unwrap()).unwrap());
    let error = |doc: &J| field(doc).iter().zip(&reference).map(|(a,b)| (a-b).abs()).fold(0.0_f64,f64::max);
    assert!(error(&tight_doc) < error(&loose_doc));
    assert!(tight_doc.path(&["transient","steps"]).unwrap().as_f64() > loose_doc.path(&["transient","steps"]).unwrap().as_f64());
}

#[test]
fn equilibrium_needs_no_refinement_and_preserves_zero_heat() {
    let mut r = Request::parse(FIXTURE).unwrap();
    for interval in &mut r.transient.as_mut().unwrap().intervals { interval.workload = Workload::Scale(0.0); }
    let doc = J::parse(&execute(&r, &CancelGate::new_clock_free()).unwrap()).unwrap();
    for t in field(&doc) { close(t, 300.0, 1e-8); }
    let run = doc.get("transient").unwrap();
    assert_eq!(run.path(&["adaptive","rejected_trials"]).unwrap().as_f64(), Some(0.0));
    close(run.f64_field("steps").unwrap(), 10.0, 0.0);
    close(run.f64_field("input_energy_j").unwrap(), 0.0, 0.0);
    close(run.f64_field("stored_energy_change_j").unwrap(), 0.0, 1e-6);
}

#[test]
fn adaptive_budgets_minimum_step_and_cancellation_refuse_without_a_trajectory() {
    for (key, replacement) in [("\"max_trials\": 4000", "\"max_trials\": 0"),
        ("\"relative_tolerance\": 0.001", "\"relative_tolerance\": -1"),
        ("\"minimum_trial_step_s\": 0.00001", "\"minimum_trial_step_s\": 31")] {
        assert!(FIXTURE.contains(key));
        assert!(Request::parse(&FIXTURE.replace(key,replacement)).is_err());
    }
    let r = Request::parse(&FIXTURE.replace("\"max_trials\": 4000", "\"max_trials\": 1")).unwrap();
    assert_eq!(execute(&r,&CancelGate::new_clock_free()).unwrap_err().code,"cooling-network-transient-budget");
    let r = Request::parse(&FIXTURE.replace("\"minimum_trial_step_s\": 0.00001", "\"minimum_trial_step_s\": 30")).unwrap();
    assert_eq!(execute(&r,&CancelGate::new_clock_free()).unwrap_err().code,"cooling-network-adaptive-resolution");
    let mut r = Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().max_steps = 2;
    assert_eq!(execute(&r,&CancelGate::new_clock_free()).unwrap_err().code,"cooling-network-transient-budget");
    let gate = CancelGate::new_clock_free(); gate.request();
    assert!(execute(&Request::parse(FIXTURE).unwrap(),&gate).is_err());
}
