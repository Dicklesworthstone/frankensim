use super::*;
use super::super::super::tests::{close,with_cx};

const FIXTURE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/size-transient-fan.json"));
fn config()->Config {
    Config {control:Control::FanSpeed,minimum:0.5,maximum:1.3,multiplier_tolerance:1e-4,temperature_tolerance_k:1e-5,max_evaluations:64}
}
fn numerical(peak:f64)->Trajectory {
    Trajectory {output:"{\"final_temperature_k\":300}\n".into(),peak_k:peak,peak_time_s:30.0,solid_solves:7,steps:3,design_gradient:None}
}
fn doc(r:&Request)->J {J::parse(&execute(r,&CancelGate::new_clock_free()).unwrap()).unwrap()}
fn values(doc:&J)->Vec<f64> {
    doc.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().map(|v|v.as_f64().unwrap()).collect()
}

#[test]
fn search_returns_the_evaluated_passing_peak_not_the_cool_final_field() {
    with_cx(|cx| {
        let c=config();
        // Two competing transient hotspots exchange control at multiplier 1.
        let response=|m:f64|(303.0-2.0*m).max(302.0-m);
        let selected=search(cx,&c,301.1,|m|Ok(numerical(response(m)))).unwrap();
        assert!(selected.passing.peak_k<=301.1);
        assert!(301.1-selected.passing.peak_k<=c.temperature_tolerance_k);
        assert!(selected.width<=c.multiplier_tolerance);
        assert!(selected.failed.unwrap().peak_k>301.1);
        close(selected.passing.peak_k,response(selected.multiplier),0.0);
        assert!(selected.history.iter().any(|t|t.peak_k>301.1));
    });
}

#[test]
fn minimum_bracket_budget_producer_and_cancellation_are_not_confused() {
    with_cx(|cx| {
        let mut c=config();c.max_evaluations=1;
        let selected=search(cx,&c,310.0,|_|Ok(numerical(300.0))).unwrap();
        assert!(selected.failed.is_none());assert_eq!(selected.history.len(),1);
        assert_eq!(search(cx,&c,310.0,|_|Ok(numerical(320.0))).err().unwrap().code,"cooling-network-transient-budget");
        c.max_evaluations=64;
        assert_eq!(search(cx,&c,310.0,|_|Ok(numerical(320.0))).err().unwrap().code,"cooling-network-transient-fan-bracket");
        assert_eq!(search(cx,&c,310.0,|_|Err(Failure {code:"upstream-test-refusal",message:"not feasibility".into()}))
            .err().unwrap().code,"upstream-test-refusal");
        let c=Config {minimum:1.0,maximum:f64::from_bits(1.0_f64.to_bits()+1),
            multiplier_tolerance:f64::EPSILON/4.0,..config()};
        assert_eq!(search(cx,&c,310.0,|m|Ok(numerical(if m==1.0 {320.0}else{310.0})))
            .err().unwrap().code,"cooling-network-transient-fan-resolution");
    });
    let gate=CancelGate::new_clock_free();gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:41,kernel_id:717,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let mut called=false;
        assert!(search(&cx,&config(),310.0,|_|{called=true;Ok(numerical(300.0))}).is_err());
        assert!(!called);
    });
}

#[test]
fn fixed_grid_design_replays_the_same_cold_initial_field_and_reports_actual_speeds() {
    let mut r=Request::parse(FIXTURE).unwrap();
    r.transient.as_mut().unwrap().adaptive=None;
    let designed=doc(&r);
    let d=designed.get("transient_fan_speed_design").unwrap();
    let multiplier=d.f64_field("selected_speed_multiplier").unwrap();
    close(multiplier,0.68603515625,0.001);
    assert!(d.f64_field("sampled_peak_objective_k").unwrap()<=320.85);
    assert!(d.path(&["failed_lower","sampled_peak_objective_k"]).unwrap().as_f64().unwrap()>320.85);
    let schedule=r.transient.as_ref().unwrap();
    assert_eq!(schedule.intervals[0].speed,Some(1.0));
    assert_eq!(schedule.intervals[1].speed,Some(1.5));
    let fresh=with_cx(|cx|simulate(&r,cx,schedule,multiplier).unwrap());
    assert_eq!(values(&designed),values(&J::parse(&fresh.output).unwrap()));
    let history=designed.path(&["transient","history"]).unwrap().as_array().unwrap();
    for sample in &history[1..] {
        assert!(sample.f64_field("objective_temperature_k").unwrap()<=320.85);
        let base=if sample.f64_field("interval")==Some(0.0){1.0}else{1.5};
        close(sample.f64_field("fan_speed_ratio").unwrap(),base*multiplier,0.0);
    }
    close(designed.path(&["fan","speed_ratio"]).unwrap().as_f64().unwrap(),1.5*multiplier,0.0);
    close(designed.path(&["fan","flow_m3_s"]).unwrap().as_f64().unwrap(),0.004*1.5*multiplier,1e-9);
    close(designed.path(&["transient","input_energy_j"]).unwrap().as_f64().unwrap(),600.0,1e-7);
    assert!(designed.get("contacts").unwrap().as_array().unwrap().len()==1);
}

#[test]
fn adaptive_design_checks_all_accepted_samples_including_initial_and_midpoints() {
    let r=Request::parse(FIXTURE).unwrap();
    let result=doc(&r);
    let d=result.get("transient_fan_speed_design").unwrap();
    close(d.f64_field("selected_speed_multiplier").unwrap(),1.1083984375,0.002);
    let run=result.get("transient").unwrap();
    let recorded=run.get("history").unwrap().as_array().unwrap();
    let peak=recorded.iter().map(|s|s.f64_field("objective_temperature_k").unwrap()).fold(f64::NEG_INFINITY,f64::max);
    close(peak,run.f64_field("sampled_peak_objective_k").unwrap(),0.0);
    assert!(peak<=320.85 && 320.85-peak<=1e-5);
    assert_eq!(run.get("first_sampled_violation_s"),Some(&J::Null));
    assert!(run.path(&["adaptive","rejected_trials"]).unwrap().as_f64().unwrap()>0.0);
    assert!(run.path(&["adaptive","largest_accepted_error_ratio"]).unwrap().as_f64().unwrap()<=1.0);
    assert!(result.path(&["objective","value_k"]).unwrap().as_f64().unwrap()<304.0);
    let trials=d.get("history").unwrap().as_array().unwrap();
    close(d.f64_field("total_accepted_steps").unwrap(),trials.iter().map(|t|t.f64_field("accepted_steps").unwrap()).sum(),0.0);
    close(d.f64_field("total_solid_solves").unwrap(),trials.iter().map(|t|t.f64_field("solid_solves").unwrap()).sum(),0.0);
    assert!(trials.iter().all(|t|t.f64_field("sampled_peak_time_s")==Some(30.0)));
}

#[test]
fn all_schedule_speed_bounds_and_target_are_admitted_before_physics() {
    for (from,to) in [
        ("\"max_speed_multiplier\": 1.3","\"max_speed_multiplier\": 1.4"),
        ("\"min_speed_multiplier\": 0.5","\"min_speed_multiplier\": 0.25"),
        ("\"max_evaluations\": 64","\"max_evaluations\": 0"),
        ("\"speed_multiplier_tolerance\": 0.0001","\"speed_multiplier_tolerance\": 0"),
        ("\"temperature_limit_k\": 320.85,",""),
    ] {
        assert!(FIXTURE.contains(from),"missing fixture mutation {from}");
        assert!(Request::parse(&FIXTURE.replace(from,to)).is_err(),"accepted {to}");
    }
    let r=Request::parse(FIXTURE).unwrap();
    let schedule=r.transient.as_ref().unwrap();
    assert!(config().validate(schedule,None).is_err());
    let mut initial_hot=Request::parse(FIXTURE).unwrap();
    let schedule=initial_hot.transient.as_mut().unwrap();
    schedule.adaptive=None;schedule.initial.fill(330.0);
    let run=with_cx(|cx|simulate(&initial_hot,cx,initial_hot.transient.as_ref().unwrap(),1.0).unwrap());
    assert!(run.peak_k>=330.0,"initial violation must not disappear after cooling");
}
