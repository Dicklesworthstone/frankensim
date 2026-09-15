use super::*;
use super::super::super::tests::{close, with_cx};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/periodic-contact-pulse.json"));
fn document(r: &Request) -> J {
    J::parse(&execute(r,&CancelGate::new_clock_free()).expect("contact-coupled cycle convergence")).unwrap()
}

#[test]
fn periodic_closure_uses_every_node_not_just_peak_or_mean() {
    with_cx(|cx| {
        // Identical maxima and means do not imply identical thermal state.
        let residual=field_residual(cx,&[300.0,310.0],&[310.0,300.0]).unwrap();
        assert_eq!(residual,10.0);
        let rule=Periodic {tolerance_k:0.01,consecutive:2};
        let mut streak=0;
        assert!(!rule.observe(0.001,&mut streak));
        assert!(!rule.observe(residual,&mut streak));
        assert_eq!(streak,0);
        assert!(!rule.observe(0.001,&mut streak));
        assert!(rule.observe(0.002,&mut streak));
    });
}

#[test]
fn periodic_contact_cycle_retains_residual_and_nonzero_storage() {
    let r=Request::parse(FIXTURE).unwrap();
    let doc=document(&r);
    let run=doc.get("repeated_cycles").unwrap();
    assert_eq!(run.str_field("status"),Some("periodic-field-tolerance-met"));
    assert_eq!(run.f64_field("cycles_completed"),Some(52.0));
    let gate=run.get("periodic").unwrap();
    assert!(gate.f64_field("full_field_residual_k").unwrap()<=1e-4);
    assert_eq!(gate.f64_field("consecutive_cycles_met"),Some(2.0));
    close(run.f64_field("sampled_peak_objective_k").unwrap(),311.672519334,1e-3);
    close(run.f64_field("input_energy_j").unwrap(),31200.0,1e-6);
    close(run.f64_field("total_accepted_steps").unwrap(),3900.0,0.0);
    let cycles=run.get("cycles").unwrap().as_array().unwrap();
    assert!(cycles[49].f64_field("start_to_end_field_residual_k").unwrap()>1e-4);
    assert!(cycles[50].f64_field("start_to_end_field_residual_k").unwrap()<=1e-4);
    let last=cycles.last().unwrap();
    assert!(last.f64_field("stored_energy_change_j").unwrap()>0.0,
        "a small observed cycle residual must not be replaced by zero stored energy");
    let last_energy=last.f64_field("stored_energy_change_j").unwrap()
        -last.f64_field("input_energy_j").unwrap()+last.f64_field("air_energy_gain_j").unwrap();
    assert!(last_energy.abs()<r.limits.heat*150.0);
}

#[test]
fn converging_from_hot_initial_state_does_not_hide_warmup_peak() {
    let mut r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_mut().unwrap();
    schedule.initial.fill(350.0);
    for interval in &mut schedule.intervals {interval.workload=Workload::Scale(0.0);}
    let doc=document(&r);
    let run=doc.get("repeated_cycles").unwrap();
    // Consistent capacitance does not guarantee a discrete maximum principle:
    // preserve any early overshoot as well as the initial hot temperature.
    let first=&run.get("cycles").unwrap().as_array().unwrap()[0];
    assert!(run.f64_field("sampled_peak_objective_k").unwrap()>=350.0);
    assert!(run.f64_field("sampled_peak_objective_k").unwrap()>=first.f64_field("sampled_peak_objective_k").unwrap());
    assert!(run.f64_field("sampled_peak_time_s").unwrap()<=150.0);
    assert_eq!(run.f64_field("first_sampled_violation_s"),Some(0.0));
    assert!(doc.path(&["transient","sampled_peak_objective_k"]).unwrap().as_f64().unwrap()<301.0);
    close(run.f64_field("input_energy_j").unwrap(),0.0,0.0);
}

#[test]
fn periodic_power_sizing_waits_for_cycle_closure_in_each_candidate() {
    let mut r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_mut().unwrap();
    schedule.limit=Some(308.0);
    schedule.power_design=Some(sizing::Config::parse_power(&J::parse(r#"{
        "min_power_multiplier":0,"max_power_multiplier":1,
        "power_multiplier_tolerance":0.001,"temperature_tolerance_k":0.001,
        "max_evaluations":32}"#).unwrap()).unwrap());
    let doc=document(&r);
    let run=doc.get("repeated_cycles").unwrap();
    assert_eq!(run.str_field("status"),Some("periodic-field-tolerance-met"));
    assert!(run.f64_field("sampled_peak_objective_k").unwrap()<=308.0);
    let design=doc.get("transient_power_design").unwrap();
    let multiplier=design.f64_field("selected_power_multiplier").unwrap();
    // Reference is the direct affine fixed point of the independent FEM cycle.
    close(multiplier,8.0/(311.672953063822-300.0),0.001);
    let cycles=run.f64_field("cycles_completed").unwrap();
    assert!(cycles>2.0);
    close(run.f64_field("input_energy_j").unwrap(),600.0*multiplier*cycles,1e-5);
    assert!(design.path(&["failed_upper","sampled_peak_objective_k"]).unwrap().as_f64().unwrap()>308.0);
}

#[test]
fn equilibrium_adaptive_cycles_converge_without_suppressing_work_counts() {
    let mut r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_mut().unwrap();
    for interval in &mut schedule.intervals {interval.workload=Workload::Scale(0.0);}
    schedule.adaptive=Some(adaptive::Config {absolute_tolerance_k:0.01,relative_tolerance:0.0,
        minimum_trial_step_s:1e-5,max_trials:4000});
    let doc=document(&r);
    let run=doc.get("repeated_cycles").unwrap();
    assert_eq!(run.f64_field("cycles_completed"),Some(2.0));
    close(run.f64_field("sampled_peak_objective_k").unwrap(),300.0,1e-8);
    assert!(run.f64_field("total_solid_solves").unwrap()>run.f64_field("total_accepted_steps").unwrap());
    assert!(run.f64_field("energy_residual_j").unwrap().abs()<1e-5);
}

#[test]
fn invalid_periodic_contracts_and_unresolved_cycles_refuse() {
    for body in [r#"{"cycles":2,"until_periodic":{},"max_total_steps":1000}"#,
        r#"{"until_periodic":{"max_cycles":10,"temperature_tolerance_k":0,"consecutive_cycles":2},"max_total_steps":1000}"#,
        r#"{"until_periodic":{"max_cycles":10,"temperature_tolerance_k":0.01,"consecutive_cycles":1},"max_total_steps":1000}"#,
        r#"{"until_periodic":{"max_cycles":1,"temperature_tolerance_k":0.01,"consecutive_cycles":2},"max_total_steps":1000}"#,
        r#"{"until_periodic":{"max_cycles":10,"temperature_tolerance_k":0.01,"consecutive_cycles":2},"max_total_steps":149}"#] {
        assert!(Config::parse(&J::parse(body).unwrap(),75,false).is_err());
    }
    let mut r=Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().repeat.as_mut().unwrap().cycles=2;
    let error=execute(&r,&CancelGate::new_clock_free()).unwrap_err();
    assert_eq!(error.code,"cooling-network-transient-budget");
    assert!(error.message.contains("periodic cycle budget exhausted"));
    assert!(error.message.contains("full-field residual"));
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(execute(&r,&gate).is_err());
}
