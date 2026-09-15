use super::*;
use super::super::super::tests::{close, with_cx};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/repeated-contact-pulse.json"));
fn document(r: &Request) -> J {
    J::parse(&execute(r,&CancelGate::new_clock_free()).expect("actual repeated FEM/air solve")).unwrap()
}
fn field(doc: &J) -> Vec<f64> {
    doc.get("solid_temperatures_k").unwrap().as_array().unwrap().iter()
        .map(|t|t.as_f64().unwrap()).collect()
}

#[test]
fn repeated_pulses_accumulate_heat_and_account_for_every_cycle() {
    let r=Request::parse(FIXTURE).unwrap();
    let doc=document(&r);
    let repeated=doc.get("repeated_cycles").unwrap();
    assert_eq!(repeated.f64_field("cycles_completed"),Some(10.0));
    close(repeated.f64_field("elapsed_time_s").unwrap(),1500.0,0.0);
    close(repeated.f64_field("total_accepted_steps").unwrap(),750.0,0.0);
    close(repeated.f64_field("input_energy_j").unwrap(),6000.0,1e-7);
    close(repeated.f64_field("sampled_peak_objective_k").unwrap(),310.748442023,1e-4);
    close(repeated.f64_field("sampled_peak_time_s").unwrap(),1380.0,0.0);
    close(field(&doc).into_iter().fold(f64::NEG_INFINITY,f64::max),304.948465080,1e-4);
    let cycles=repeated.get("cycles").unwrap().as_array().unwrap();
    assert!(cycles[0].f64_field("sampled_peak_objective_k").unwrap()<308.0);
    assert!(cycles[2].f64_field("sampled_peak_objective_k").unwrap()>308.0);
    let first=repeated.f64_field("first_sampled_violation_s").unwrap();
    assert!(first>300.0 && first<=330.0);
    assert!(repeated.f64_field("energy_residual_j").unwrap().abs()<r.limits.heat*1500.0);
    // The old transient record intentionally describes only the last cycle.
    close(doc.path(&["transient","time_s"]).unwrap().as_f64().unwrap(),150.0,0.0);
    close(repeated.f64_field("last_cycle_start_time_s").unwrap(),1350.0,0.0);
}

#[test]
fn one_repeat_preserves_the_original_numerical_result() {
    let mut r=Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().repeat.as_mut().unwrap().cycles=1;
    let repeated=document(&r);
    r.transient.as_mut().unwrap().repeat=None;
    let plain=document(&r);
    assert_eq!(plain.get("transient"),repeated.get("transient"));
    assert_eq!(plain.get("solid_temperatures_k"),repeated.get("solid_temperatures_k"));
    assert_eq!(plain.get("fan"),repeated.get("fan"));
}

#[test]
fn cycle_boundaries_pass_the_entire_field_to_the_next_solve() {
    let mut r=Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().repeat.as_mut().unwrap().cycles=3;
    with_cx(|cx| {
        let schedule=r.transient.as_ref().unwrap();
        let repeated=simulate(&r,cx,schedule,1.0,schedule.repeat.unwrap()).unwrap();
        let mut old=schedule.initial.clone();
        let mut manual_peak=f64::NEG_INFINITY;
        for _ in 0..3 {
            let cycle=simulate_cycle(&r,cx,schedule,1.0,&old,schedule.max_steps).unwrap();
            manual_peak=manual_peak.max(cycle.trajectory.peak_k);
            old=cycle.final_temperature;
        }
        assert_eq!(repeated.peak_k,manual_peak);
        let published=field(&J::parse(&repeated.output).unwrap());
        assert_eq!(old,published);
    });
}

#[test]
fn power_sizing_keeps_repetitions_and_does_not_reset_each_cycle() {
    let mut r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_mut().unwrap();
    schedule.repeat.as_mut().unwrap().cycles=3;
    schedule.limit=Some(307.0);
    schedule.power_design=Some(sizing::Config::parse_power(&J::parse(r#"{
        "min_power_multiplier":0,"max_power_multiplier":1,
        "power_multiplier_tolerance":0.001,"temperature_tolerance_k":0.001,
        "max_evaluations":32}"#).unwrap()).unwrap());
    let doc=document(&r);
    let design=doc.get("transient_power_design").unwrap();
    let repeated=doc.get("repeated_cycles").unwrap();
    assert_eq!(repeated.f64_field("cycles_completed"),Some(3.0));
    let multiplier=design.f64_field("selected_power_multiplier").unwrap();
    close(multiplier,7.0/(308.327173521872-300.0),0.001);
    assert!(design.f64_field("sampled_peak_objective_k").unwrap()<=307.0);
    assert!(design.path(&["failed_upper","sampled_peak_objective_k"]).unwrap().as_f64().unwrap()>307.0);
    close(repeated.f64_field("input_energy_j").unwrap(),1800.0*multiplier,1e-6);
    assert!(design.f64_field("total_accepted_steps").unwrap()>=2.0*225.0);
}

#[test]
fn invalid_counts_total_budget_and_cancellation_do_not_publish_partial_cycles() {
    for text in [r#"{"cycles":0,"max_total_steps":1000}"#,
        r#"{"cycles":2,"max_total_steps":149}"#,
        r#"{"cycles":2,"max_total_steps":299}"#,
        r#"{"cycles":2,"max_total_steps":1000,"unexpected":true}"#] {
        let adaptive=text.contains("299");
        assert!(Config::parse(&J::parse(text).unwrap(),75,adaptive).is_err());
    }
    let mut r=Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().repeat.as_mut().unwrap().max_total_steps=76;
    assert_eq!(execute(&r,&CancelGate::new_clock_free()).unwrap_err().code,"cooling-network-transient-budget");
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(execute(&r,&gate).is_err());
    with_cx(|cx| {
        assert!(field_residual(cx,&[300.0],&[f64::NAN]).is_err());
        assert!(field_residual(cx,&[],&[]).is_err());
    });
}
