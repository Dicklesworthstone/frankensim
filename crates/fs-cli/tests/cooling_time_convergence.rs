//! Full command regressions: the ladder invokes the actual cooling producer,
//! and its final numerical grid is replayed without the study controller.
#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/radiative-contact-pulse.json"));
const NONMATCHING: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonmatching-contact-hotspot.json"));
fn num(value: f64) -> J { J::Number { value, raw: value.to_string() } }
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object") };
    &mut fields.iter_mut().find(|(name,_)| name == key).unwrap().1
}
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object") };
    if let Some((_,slot)) = fields.iter_mut().find(|(name,_)| name == key) { *slot = value; }
    else { fields.push((key.into(),value)); }
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object") };
    fields.retain(|(name,_)| name != key);
}
fn rows(root: &mut J) -> &mut Vec<J> { let J::Array(a) = root else { panic!("array") }; a }
fn encode(root: &J) -> String {
    match root {
        J::Null => "null".into(), J::Bool(b) => b.to_string(), J::Number { raw,.. } => raw.clone(),
        J::Str(s) => format!("\"{}\"",s.replace('\\',"\\\\").replace('"',"\\\"").replace('\n',"\\n")),
        J::Array(a) => format!("[{}]",a.iter().map(encode).collect::<Vec<_>>().join(",")),
        J::Object(a) => format!("{{{}}}",a.iter().map(|(k,v)|format!("{}:{}",encode(&J::Str(k.clone())),encode(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(input: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(encode(input).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(input: &J) -> J {
    let result = output(input);
    assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stderr));
    J::parse(std::str::from_utf8(&result.stdout).unwrap()).unwrap()
}
fn n(root: &J, key: &str) -> f64 { root.f64_field(key).unwrap() }
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol,"{a} versus {b}"); }
fn input(repeated: bool) -> J {
    let mut root = J::parse(BASE).unwrap();
    let schedule = member(&mut root,"transient");
    put(schedule,"max_step_s",num(1.0));
    put(schedule,"intervals",J::parse(r#"[{"duration_s":2.3,"power_scale":1,"fan_speed_ratio":1},{"duration_s":3.1,"power_scale":0,"fan_speed_ratio":1.5}]"#).unwrap());
    put(schedule,"time_convergence",J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":100,"max_total_steps":1000,"max_trace_bytes":1048576}"#).unwrap());
    if repeated { put(schedule,"repeat",J::parse(r#"{"cycles":2,"max_total_steps":1000}"#).unwrap()); }
    root
}
fn replay_request(input: &J, result: &J) -> J {
    let counts = result.path(&["time_convergence","final_steps_per_interval"]).unwrap().as_array().unwrap();
    let mut replay = input.clone();
    let schedule = member(&mut replay,"transient");
    remove(schedule,"time_convergence");
    let intervals = rows(member(schedule,"intervals"));
    assert_eq!(counts.len(),intervals.len());
    for (interval,count) in intervals.iter_mut().zip(counts) { put(interval,"steps",count.clone()); }
    replay
}
fn check_physical_replay(input: &J, result: &J) {
    let replay = run(&replay_request(input,result));
    let mut physical = result.clone();
    remove(&mut physical,"time_convergence");
    assert_eq!(physical,replay,"study observation must not change accepted trajectory fields, energy or work");
}

#[test]
fn radiating_single_and_repeated_grids_replay_exactly_and_count_only_final_physical_energy() {
    for repeated in [false,true] {
        let input = input(repeated);
        let result = run(&input);
        let study = result.get("time_convergence").unwrap();
        assert_eq!(study.str_field("status"),Some("successive-time-grid-tolerance-met"));
        assert_eq!(n(study,"trajectories_solved"),3.0);
        let history = study.get("history").unwrap().as_array().unwrap();
        for (row,expected) in history.iter().zip(["[3,4]","[6,8]","[12,16]"]) {
            assert_eq!(row.get("steps_per_interval"),Some(&J::parse(expected).unwrap()));
        }
        assert_eq!(n(study,"total_accepted_steps"),if repeated {98.0} else {49.0});
        let sum: f64 = history.iter().map(|row|n(row,"solid_solves")).sum();
        assert_eq!(n(study,"total_solid_solves"),sum);
        let trajectory = result.get(if repeated {"repeated_cycles"} else {"transient"}).unwrap();
        near(n(trajectory,"input_energy_j"),if repeated {92.0} else {46.0},1e-7);
        assert!(n(trajectory,"radiative_energy_loss_j") > 0.0);
        near(n(trajectory,"stored_energy_change_j")+n(trajectory,"air_energy_gain_j")
            +n(trajectory,"radiative_energy_loss_j"),n(trajectory,"input_energy_j"),1e-5);
        assert!(sum > n(trajectory,"total_solid_solves"));
        check_physical_replay(&input,&result);
        assert_eq!(result,run(&input));
    }
}

#[test]
fn initial_peak_agreement_cannot_hide_an_inaccurate_cooldown_field() {
    let mut input = input(true);
    put(&mut input,"objective",J::parse(r#"{"mean_wall_region":"first-face","gradient":false}"#).unwrap());
    let schedule = member(&mut input,"transient");
    put(schedule,"initial_temperature_k",num(350.0));
    for row in rows(member(schedule,"intervals")) { put(row,"power_scale",num(0.0)); }
    let result = run(&input);
    let history = result.path(&["time_convergence","history"]).unwrap().as_array().unwrap();
    for row in &history[1..] {
        near(n(row,"sampled_peak_change_k"),0.0,1e-9);
        assert!(n(row,"common_endpoint_field_change_k") > 1e-5);
    }
    // Tightening the FIELD criterion must refuse even though the selected
    // objective's initial maximum agrees throughout the grid ladder.
    put(member(member(&mut input,"transient"),"time_convergence"),"temperature_tolerance_k",num(1e-8));
    let refused = output(&input);
    assert_eq!(refused.status.code(),Some(6));
    assert!(refused.stdout.is_empty());
}

#[test]
fn independent_contact_meshes_remain_in_the_full_time_refinement_problem() {
    let mut input = J::parse(NONMATCHING).unwrap();
    put(member(&mut input,"objective"),"gradient",J::Bool(false));
    let template = self::input(false);
    let mut schedule = template.get("transient").unwrap().clone();
    remove(&mut schedule,"element_heat_capacities_j_m3_k");
    remove(&mut schedule,"nonlinear");
    put(&mut schedule,"volumetric_heat_capacity_j_m3_k",num(2_000_000.0));
    put(&mut input,"transient",schedule);
    put(&mut input,"radiation",template.get("radiation").unwrap().clone());
    let result = run(&input);
    assert_eq!(result.get("solid_temperatures_k").unwrap().as_array().unwrap().len(),32);
    let contact = &result.get("contacts").unwrap().as_array().unwrap()[0];
    assert_eq!(contact.str_field("discretization"),Some("planar-common-refinement-P1"));
    assert_eq!(contact.get("face_pairs"),Some(&J::Null));
    check_physical_replay(&input,&result);
}

#[test]
fn cumulative_and_per_trajectory_budgets_do_not_reset_with_each_grid() {
    let input = input(true);
    for (path,key,value) in [
        (vec!["transient","time_convergence"],"max_total_steps",41.0),
        (vec!["transient","time_convergence"],"max_trace_bytes",1.0),
        (vec!["transient","repeat"],"max_total_steps",14.0),
        (vec!["transient"],"max_steps",7.0),
        (vec!["budgets"],"wall_seconds",1e-12),
    ] {
        let mut trial = input.clone();
        let mut target = &mut trial;
        for part in path { target = member(target,part); }
        put(target,key,num(value));
        let result = output(&trial);
        assert_eq!(result.status.code(),Some(6),"{}",String::from_utf8_lossy(&result.stderr));
        assert!(result.stdout.is_empty());
    }
}

#[test]
fn unsupported_policies_and_nonrefining_explicit_steps_are_rejected() {
    let input = input(false);
    for (key,value) in [
        ("adaptive",r#"{"absolute_tolerance_k":1,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":100}"#),
        ("adjoint",r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#),
        ("repeat",r#"{"until_periodic":{"max_cycles":3,"temperature_tolerance_k":1,"consecutive_cycles":2},"max_total_steps":1000}"#),
        ("power_design",r#"{"min_power_multiplier":0,"max_power_multiplier":1,"power_multiplier_tolerance":0.01,"temperature_tolerance_k":0.01,"max_evaluations":10}"#),
    ] {
        let mut trial = input.clone();put(member(&mut trial,"transient"),key,J::parse(value).unwrap());
        let result = output(&trial);assert!(!result.status.success());assert!(result.stdout.is_empty());
    }
    for count in [0.0,2.0,3.5] {
        let mut trial = input.clone();remove(member(&mut trial,"transient"),"time_convergence");
        put(&mut rows(member(member(&mut trial,"transient"),"intervals"))[0],"steps",num(count));
        let result = output(&trial);assert!(!result.status.success());assert!(result.stdout.is_empty());
    }
}

#[test]
fn ordinary_runs_keep_the_same_fields_when_natural_counts_are_made_explicit() {
    let mut implicit = input(false);remove(member(&mut implicit,"transient"),"time_convergence");
    let mut explicit = implicit.clone();
    for (row,count) in rows(member(member(&mut explicit,"transient"),"intervals")).iter_mut().zip([3.0,4.0]) {
        put(row,"steps",num(count));
    }
    assert_eq!(run(&implicit),run(&explicit));
}
