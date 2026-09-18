//! Actual binary regressions for independent component allocation. The model
//! is never replaced by the scalar test seam used by the search's unit tests.
#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-component-contact-pulse.json"));
const SPEC: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/allocate-component-pulse.json"));
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path=std::env::temp_dir().join(format!("fs-component-allocation-{}-{}",
            std::process::id(),NEXT.fetch_add(1,Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap(); Self(path)
    }
    fn allocate(&self, base: &J, spec: &J) -> Output {
        std::fs::write(self.0.join("base.json"),text(base)).unwrap();
        std::fs::write(self.0.join("allocation.json"),text(spec)).unwrap();
        let out=Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json","cooling-component-design","base.json","allocation.json"]).output().unwrap();
        assert_eq!(std::fs::read_to_string(self.0.join("base.json")).unwrap(),text(base));
        assert_eq!(std::fs::read_to_string(self.0.join("allocation.json")).unwrap(),text(spec));
        out
    }
    fn cooling(&self, request: &J) -> J {
        std::fs::write(self.0.join("replay.json"),text(request)).unwrap();
        success(&Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json","cooling-network","replay.json"]).output().unwrap())
    }
}
impl Drop for Scratch { fn drop(&mut self) { let _=std::fs::remove_dir_all(&self.0); } }
fn number(value: f64) -> J { J::Number {value,raw:value.to_string()} }
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(rows)=root else {panic!("object")};
    &mut rows.iter_mut().find(|(name,_)|name==key).unwrap().1
}
fn put(root: &mut J,key: &str,value: J) {
    let J::Object(rows)=root else {panic!("object")};
    if let Some((_,v))=rows.iter_mut().find(|(name,_)|name==key) { *v=value; }
    else {rows.push((key.into(),value));}
}
fn remove(root: &mut J,key: &str) {
    let J::Object(rows)=root else {panic!("object")};rows.retain(|(name,_)|name!=key);
}
fn rows(root: &mut J) -> &mut Vec<J> {let J::Array(rows)=root else {panic!("array")};rows}
fn text(root: &J) -> String {
    match root {
        J::Null=>"null".into(),J::Bool(v)=>v.to_string(),J::Number {raw,..}=>raw.clone(),
        J::Str(s)=>format!("\"{}\"",s.replace('\\',"\\\\").replace('"',"\\\"").replace('\n',"\\n")),
        J::Array(a)=>format!("[{}]",a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(a)=>format!("{{{}}}",a.iter().map(|(k,v)|format!("{}:{}",text(&J::Str(k.clone())),text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn parsed(out: &Output) -> J {J::parse(std::str::from_utf8(&out.stdout).unwrap()).unwrap()}
fn success(out: &Output) -> J {
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));parsed(out)
}
fn n(root: &J,key: &str) -> f64 {root.f64_field(key).unwrap()}
fn near(a: f64,b: f64,t: f64) {assert!((a-b).abs()<=t,"{a} versus {b}");}
fn selected(result: &J,name: &str) -> f64 {
    let row=result.get("selected").unwrap().as_array().unwrap().iter()
        .find(|r|r.str_field("component")==Some(name)).unwrap();n(row,"selected_power_w")
}
fn verify_replay(dir: &Scratch,result: &J) {
    let cooling=result.get("cooling_result").unwrap();
    assert_eq!(&dir.cooling(result.get("resolved_request").unwrap()),cooling);
    let phase=cooling.get("repeated_cycles").unwrap_or_else(||cooling.get("transient").unwrap());
    near(n(result,"selected_sampled_peak_k"),n(phase,"sampled_peak_objective_k"),0.0);
    assert!(n(result,"selected_sampled_peak_k")<=n(result,"temperature_limit_k"));
    let values: Vec<J>=result.get("selected").unwrap().as_array().unwrap().iter()
        .map(|r|r.get("selected_power_w").unwrap().clone()).collect();
    assert!(result.get("history").unwrap().as_array().unwrap().iter()
        .any(|r|r.get("power_w")==Some(&J::Array(values.clone())) && r.get("passing")==Some(&J::Bool(true))));
}

#[test]
fn allocation_activates_a_dormant_footprint_and_replays_the_exact_repeated_result() {
    let dir=Scratch::new();let base=J::parse(BASE).unwrap();let spec=J::parse(SPEC).unwrap();
    let original=dir.allocate(&base,&spec);let result=success(&original);
    assert_eq!(result.str_field("status"),Some("priority-allocation-complete"));
    assert_eq!(n(&result,"completed_priorities"),3.0);
    assert_eq!(selected(&result,"standby"),3.0);assert_eq!(selected(&result,"memory"),10.0);
    near(selected(&result,"chip"),14.734144111050925,0.0011);
    assert!(n(&result,"newton_trials")>0.0);
    let resolved=result.get("resolved_request").unwrap();
    for key in ["solid","hydraulics","air","radiation","budgets","objective"] {
        assert_eq!(resolved.get(key),base.get(key));
    }
    assert_eq!(resolved.path(&["transient","intervals"]).unwrap().as_array().unwrap()[1],
        base.path(&["transient","intervals"]).unwrap().as_array().unwrap()[1]);
    let chosen=result.get("cooling_result").unwrap();
    assert_eq!(n(chosen.get("repeated_cycles").unwrap(),"cycles_completed"),2.0);
    near(n(chosen.get("repeated_cycles").unwrap(),"sampled_peak_time_s"),20.0,0.0);
    assert!(n(chosen.get("objective").unwrap(),"value_k")<n(&result,"selected_sampled_peak_k"));
    verify_replay(&dir,&result);
    let replay=dir.allocate(&base,&spec);success(&replay);assert_eq!(original.stdout,replay.stdout);
}

#[test]
fn derivative_free_and_reversed_priorities_change_only_the_declared_policy() {
    let dir=Scratch::new();let mut base=J::parse(BASE).unwrap();let spec=J::parse(SPEC).unwrap();
    remove(member(&mut base,"transient"),"adjoint");
    let result=success(&dir.allocate(&base,&spec));verify_replay(&dir,&result);
    assert_eq!(result.str_field("search_method"),Some("priority-component-bisection"));
    assert_eq!(n(&result,"newton_trials"),0.0);
    assert_eq!(result.path(&["cooling_result","repeated_cycles","adjoint"]),Some(&J::Null));
    near(selected(&result,"chip"),14.734144111050925,0.0011);
    let mut reverse=spec;rows(member(&mut reverse,"priority")).reverse();
    let changed=success(&dir.allocate(&base,&reverse));verify_replay(&dir,&changed);
    assert!(selected(&changed,"chip")>selected(&result,"chip")+1.0);
    assert!(selected(&changed,"standby")<0.01);
    assert!(selected(&changed,"memory")<0.03);
}

#[test]
fn partial_budgets_publish_only_the_last_evaluated_passing_allocation() {
    let dir=Scratch::new();let mut base=J::parse(BASE).unwrap();remove(member(&mut base,"transient"),"adjoint");
    for (key,value,completed,evaluations) in [
        ("max_evaluations",1.0,0.0,1.0),("max_evaluations",2.0,1.0,2.0),
        ("max_total_steps",14.0,0.0,1.0),
    ] {
        let mut spec=J::parse(SPEC).unwrap();put(&mut spec,key,number(value));
        let out=dir.allocate(&base,&spec);assert_eq!(out.status.code(),Some(6));
        let result=parsed(&out);assert_eq!(result.str_field("status"),Some("budget-exhausted"));
        assert_eq!(n(&result,"completed_priorities"),completed);
        assert_eq!(n(&result,"evaluations_completed"),evaluations);
        assert_eq!(n(&result,"total_completed_trajectory_steps"),14.0*evaluations);
        assert_eq!(selected(&result,"standby"),if completed==1.0 {3.0} else {0.0});
        assert_eq!(selected(&result,"chip"),0.0);verify_replay(&dir,&result);
    }
}

#[test]
fn model_failures_bad_controls_and_missing_peak_adjoints_never_become_hot_trials() {
    let dir=Scratch::new();let base=J::parse(BASE).unwrap();let spec=J::parse(SPEC).unwrap();
    let mut wrong=base.clone();
    put(member(member(&mut wrong,"transient"),"adjoint"),"qoi",J::Str("final".into()));
    let out=dir.allocate(&wrong,&spec);assert!(!out.status.success());assert!(out.stdout.is_empty());
    let mut bad_model=base.clone();put(member(&mut bad_model,"radiation"),"max_iterations",number(1.0));
    let out=dir.allocate(&bad_model,&spec);assert!(!out.status.success());assert!(out.stdout.is_empty());
    let mut duplicate=spec.clone();let priorities=rows(member(&mut duplicate,"priority"));priorities.push(priorities[0].clone());
    let out=dir.allocate(&base,&duplicate);assert!(!out.status.success());assert!(out.stdout.is_empty());
    let mut impossible=spec.clone();put(&mut impossible,"temperature_limit_k",number(299.0));
    let out=dir.allocate(&base,&impossible);assert!(!out.status.success());assert!(out.stdout.is_empty());
    let mut timed=spec;put(&mut timed,"wall_seconds",number(1e-9));
    let out=dir.allocate(&base,&timed);assert_eq!(out.status.code(),Some(6));assert!(out.stdout.is_empty());
}

#[test]
fn enclosure_reflections_remain_in_allocation_and_selected_trajectory_adjoints() {
    let dir=Scratch::new();
    let base=J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/adjoint-enclosure-pulse.json"))).unwrap();
    let mut spec=J::parse(SPEC).unwrap();put(&mut spec,"temperature_limit_k",number(316.0));
    put(&mut spec,"priority",J::parse(r#"[{"component":"heater","interval":0,"min_power_w":1,"max_power_w":8}]"#).unwrap());
    let result=success(&dir.allocate(&base,&spec));verify_replay(&dir,&result);
    assert!(selected(&result,"heater")>5.0 && selected(&result,"heater")<8.0);
    assert_eq!(result.path(&["cooling_result","radiation","model"]).unwrap().as_str(),Some("closed-gray-diffuse-enclosure"));
    near(n(result.path(&["cooling_result","radiation"]).unwrap(),"radiative_out_w"),0.0,1e-7);
    assert_eq!(result.path(&["resolved_request","radiation"]),base.get("radiation"));
    assert!(result.path(&["cooling_result","repeated_cycles","adjoint","radiation"]).unwrap().as_object().is_some());
}
