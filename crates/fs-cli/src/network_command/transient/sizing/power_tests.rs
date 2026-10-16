use super::*;
use super::super::super::tests::{close,with_cx};

const FIXTURE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/size-transient-power.json"));
fn config()->Config {
    Config {control:Control::WorkloadPower,minimum:0.0,maximum:2.0,multiplier_tolerance:1e-4,
        temperature_tolerance_k:1e-5,max_evaluations:64}
}
fn numerical(peak:f64)->Trajectory {
    Trajectory {output:"{}\n".into(),peak_k:peak,peak_time_s:30.0,solid_solves:1,steps:1,design_gradient:None}
}

#[test]
fn power_search_keeps_the_passing_lower_side_and_never_returns_the_failed_upper() {
    with_cx(|cx| {
        let c=config();
        let chosen=search(cx,&c,306.0,|p|Ok(numerical(300.0+5.0*p))).unwrap();
        assert!(chosen.passing.peak_k<=306.0);
        assert!(306.0-chosen.passing.peak_k<=c.temperature_tolerance_k);
        assert!(chosen.failed.unwrap().multiplier>chosen.multiplier);
        assert!(chosen.failed.unwrap().peak_k>306.0);
        close(chosen.multiplier,1.2,3e-6);
        let boundary=search(cx,&c,320.0,|p|Ok(numerical(300.0+5.0*p))).unwrap();
        close(boundary.multiplier,2.0,0.0);
        assert!(boundary.failed.is_none());assert_eq!(boundary.history.len(),1);
        assert_eq!(search(cx,&c,299.0,|p|Ok(numerical(300.0+5.0*p))).err().unwrap().code,
            "cooling-network-transient-power-bracket");
    });
}

#[test]
fn scaling_preserves_footprints_zero_loads_initial_state_and_fan_schedule() {
    let r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_ref().unwrap();
    with_cx(|cx| {
        let scaled=power_schedule(cx,schedule,0.5).unwrap();
        assert_eq!(scaled.initial,schedule.initial);
        assert_eq!(scaled.capacities,schedule.capacities);
        assert!(scaled.power_design.is_none() && scaled.fan_speed_design.is_none());
        for (base,next) in schedule.intervals.iter().zip(&scaled.intervals) {
            assert_eq!(base.speed,next.speed);assert_eq!(base.duration,next.duration);assert_eq!(base.steps,next.steps);
            let b=base.workload.prepare(&r,cx).unwrap();
            let n=next.workload.prepare(&r,cx).unwrap();
            for v in 0..r.mesh.vertex_count() {close(n.source.at(v),0.5*b.source.at(v),1e-8);}
        }
        let zero=power_schedule(cx,schedule,0.0).unwrap();
        for interval in &zero.intervals {
            close(interval.workload.prepare(&r,cx).unwrap().expected_power_w.unwrap(),0.0,0.0);
        }
    });
    let Workload::Scale(s)=scaled_workload(&Workload::Scale(0.25),2.0).unwrap() else {panic!("wrong mode")};
    close(s,0.5,0.0);
    let powers=Workload::Components(BTreeMap::from([("cpu".into(),0.0),("gpu".into(),20.0)]));
    let Workload::Components(scaled)=scaled_workload(&powers,0.5).unwrap() else {panic!("wrong mode")};
    assert_eq!(scaled["cpu"],0.0);assert_eq!(scaled["gpu"],10.0);
    for invalid in [-1.0,f64::NAN,f64::INFINITY] {assert!(scaled_workload(&powers,invalid).is_err());}
}

#[test]
fn adaptive_power_design_returns_actual_scaled_watts_and_complete_passing_trajectory() {
    let r=Request::parse(FIXTURE).unwrap();
    let result=J::parse(&execute(&r,&CancelGate::new_clock_free()).unwrap()).unwrap();
    let d=result.get("transient_power_design").unwrap();
    assert!(result.get("transient_fan_speed_design").is_none());
    let multiplier=d.f64_field("selected_power_multiplier").unwrap();
    close(multiplier,0.7656970024108887,2e-4);
    assert!(d.f64_field("sampled_peak_objective_k").unwrap()<=315.0);
    assert!(315.0-d.f64_field("sampled_peak_objective_k").unwrap()<=1e-5);
    assert!(d.path(&["failed_upper","sampled_peak_objective_k"]).unwrap().as_f64().unwrap()>315.0);
    let t=result.get("transient").unwrap();
    close(t.f64_field("input_energy_j").unwrap(),600.0*multiplier,1e-6);
    assert!(t.f64_field("energy_residual_j").unwrap().abs()<r.limits.heat*150.0);
    assert_eq!(t.get("first_sampled_violation_s"),Some(&J::Null));
    for sample in &t.get("history").unwrap().as_array().unwrap()[1..] {
        let first=sample.f64_field("interval")==Some(0.0);
        close(sample.f64_field("source_w").unwrap(),if first {20.0*multiplier}else{0.0},1e-7);
        close(sample.path(&["component_powers_w","chip"]).unwrap().as_f64().unwrap(),if first {20.0*multiplier}else{0.0},1e-12);
        close(sample.f64_field("fan_speed_ratio").unwrap(),if first {1.0}else{1.5},0.0);
        assert!(sample.f64_field("objective_temperature_k").unwrap()<=315.0);
    }
    close(result.path(&["fan","speed_ratio"]).unwrap().as_f64().unwrap(),1.5,0.0);
    close(result.path(&["fan","flow_m3_s"]).unwrap().as_f64().unwrap(),0.006,1e-9);
    with_cx(|cx| {
        let applied=power_schedule(cx,r.transient.as_ref().unwrap(),multiplier).unwrap();
        let replay=J::parse(&simulate(&r,cx,&applied,1.0).unwrap().output).unwrap();
        assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    });
}

#[test]
fn power_schema_rejects_ambiguous_controls_and_admits_zero_lower_bound() {
    let r=Request::parse(FIXTURE).unwrap();
    assert_eq!(r.transient.as_ref().unwrap().power_design.as_ref().unwrap().minimum,0.0);
    for (from,to) in [("\"min_power_multiplier\": 0","\"min_power_multiplier\": -1"),
        ("\"max_power_multiplier\": 2","\"max_power_multiplier\": 0"),
        ("\"power_design\": {","\"fan_speed_design\": {},\"power_design\": {"),
        ("\"temperature_limit_k\": 315,",""),
        ("\"max_evaluations\": 64","\"max_evaluations\": 0")] {
        assert!(FIXTURE.contains(from));
        assert!(Request::parse(&FIXTURE.replace(from,to)).is_err());
    }
    // Workload sizing is also valid with a prescribed-pressure drive: it does
    // not require or invent a fan when the flow is already prescribed.
    assert!(config().validate(r.transient.as_ref().unwrap(),None).is_ok());
}
