use super::*;
use super::super::tests::{close,with_cx};

const FIXTURE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/transient-contact-pulse.json"));
fn request()->Request {Request::parse(FIXTURE).unwrap()}
fn final_field(doc:&J)->Vec<f64> {
    doc.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().map(|t|t.as_f64().unwrap()).collect()
}

#[test]
fn pulse_contact_and_fan_schedule_match_independent_transient_fem() {
    let r=request();
    let out=execute(&r,&CancelGate::new_clock_free()).unwrap();
    let doc=J::parse(&out).unwrap();
    let run=doc.get("transient").unwrap();
    close(run.f64_field("time_s").unwrap(),150.0,1e-12);
    assert_eq!(run.f64_field("steps"),Some(75.0));
    close(run.f64_field("sampled_peak_objective_k").unwrap(),304.11663335714235,3e-5);
    close(run.f64_field("sampled_peak_time_s").unwrap(),30.0,1e-12);
    close(run.f64_field("first_sampled_violation_s").unwrap(),24.0,1e-12);
    close(doc.get("objective").unwrap().f64_field("value_k").unwrap(),301.73054647705186,3e-5);
    close(run.f64_field("input_energy_j").unwrap(),600.0,1e-7);
    close(run.f64_field("stored_energy_change_j").unwrap(),533.2318718045553,3e-5);
    close(run.f64_field("air_energy_gain_j").unwrap(),66.76812788943381,3e-5);
    assert!(run.f64_field("energy_residual_j").unwrap().abs()<r.limits.heat*150.0);
    assert_eq!(doc.path(&["fan","speed_ratio"]).unwrap().as_f64(),Some(1.5));
    close(doc.path(&["fan","flow_m3_s"]).unwrap().as_f64().unwrap(),0.006,1e-9);
    let history=run.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.len(),76);
    assert!(history[1..16].iter().all(|s|s.f64_field("source_w").unwrap()>19.99999));
    assert!(history[16..].iter().all(|s|s.f64_field("source_w").unwrap()==0.0));
    assert_eq!(doc.get("adjoint_residual"),Some(&J::Null));
    assert_eq!(doc.get("contacts").unwrap().as_array().unwrap().len(),1);
}

#[test]
fn equilibrium_remains_uniform_and_keeps_zero_storage_and_exhaust() {
    let mut r=request();
    let schedule=r.transient.as_mut().unwrap();
    for interval in &mut schedule.intervals {interval.workload=Workload::Scale(0.0);}
    let doc=J::parse(&execute(&r,&CancelGate::new_clock_free()).unwrap()).unwrap();
    for t in final_field(&doc) {close(t,300.0,1e-8);}
    let run=doc.get("transient").unwrap();
    close(run.f64_field("stored_energy_change_j").unwrap(),0.0,1e-6);
    close(run.f64_field("air_energy_gain_j").unwrap(),0.0,1e-6);
    assert_eq!(run.get("first_sampled_violation_s"),Some(&J::Null));
}

#[test]
fn splitting_at_an_accepted_endpoint_preserves_the_numerical_tail() {
    let text=FIXTURE.replace("\"duration_s\": 30","\"duration_s\": 4").replace("\"duration_s\": 120","\"duration_s\": 4");
    let mut all=Request::parse(&text).unwrap();
    // A single constant interval of two steps versus two one-step invocations.
    let schedule=all.transient.as_mut().unwrap();
    schedule.intervals.truncate(1);schedule.total_steps=2;
    let complete=J::parse(&execute(&all,&CancelGate::new_clock_free()).unwrap()).unwrap();
    let mut first=Request::parse(&text).unwrap();
    let s=first.transient.as_mut().unwrap();s.intervals.truncate(1);s.intervals[0].duration=2.0;s.intervals[0].steps=1;s.total_steps=1;
    let middle=J::parse(&execute(&first,&CancelGate::new_clock_free()).unwrap()).unwrap();
    first.transient.as_mut().unwrap().initial=final_field(&middle);
    let resumed=J::parse(&execute(&first,&CancelGate::new_clock_free()).unwrap()).unwrap();
    for (a,b) in final_field(&complete).iter().zip(final_field(&resumed)) {close(*a,b,1e-8);}
}

#[test]
fn scheduled_fan_updates_flow_derived_coefficients_not_only_the_label() {
    let fixture=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/fan-correlated-hotspot.json"));
    let mut r=Request::parse(fixture).unwrap();r.gradient=false;
    let nominal=J::parse(&execute(&r,&CancelGate::new_clock_free()).unwrap()).unwrap();
    let spec=J::parse(r#"{"initial_temperature_k":300,"volumetric_heat_capacity_j_m3_k":2000000,"max_step_s":1,"max_steps":4,
        "intervals":[{"duration_s":2,"power_scale":1,"fan_speed_ratio":1},{"duration_s":2,"power_scale":1,"fan_speed_ratio":1.5}]}"#).unwrap();
    r.transient=Some(Schedule::parse(&spec,r.mesh.vertex_count(),r.mesh.element_count(),r.fan.as_ref()).unwrap());
    let result=J::parse(&execute(&r,&CancelGate::new_clock_free()).unwrap()).unwrap();
    close(result.path(&["fan","flow_m3_s"]).unwrap().as_f64().unwrap(),0.006,1e-9);
    let before=nominal.get("convection").unwrap().as_array().unwrap();
    let after=result.get("convection").unwrap().as_array().unwrap();
    assert!(!before.is_empty());
    assert!(before.iter().zip(after).any(|(a,b)|b.f64_field("htc_w_m2_k").unwrap()>a.f64_field("htc_w_m2_k").unwrap()));
}

#[test]
fn incompatible_steady_features_time_budgets_and_cancellation_refuse() {
    for (a,b) in [
        ("\"gradient\": false","\"gradient\": true"),
        ("\"max_step_s\": 2","\"max_step_s\": 0"),
        ("\"max_steps\": 1000","\"max_steps\": 1"),
        ("\"power_scale\": 1","\"power_scale\": -1"),
        ("\"fan_speed_ratio\": 1","\"fan_speed_ratio\": 0"),
        ("\"initial_temperature_k\": 300","\"initial_temperature_k\": 0"),
    ] {
        assert!(FIXTURE.contains(a),"fixture mutation did not apply: {a}");
        assert!(Request::parse(&FIXTURE.replace(a,b)).is_err(),"accepted {a}->{b}");
    }
    let r=request();let gate=CancelGate::new_clock_free();gate.request();
    assert!(execute(&r,&gate).is_err());
    let mut r=request();r.limits.coupling=1;
    with_cx(|cx|assert!(solve(&r,cx,r.transient.as_ref().unwrap()).is_err()));
}
