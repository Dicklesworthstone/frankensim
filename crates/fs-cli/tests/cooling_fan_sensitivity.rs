//! Actual-command tests; no replacement cooling evaluator is used.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const CORRELATED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-correlated-hotspot.json"));
const DECLARED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-hotspot.json"));
const PULSE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/nonlinear-contact-pulse.json"));

fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(rows) = root else { panic!("object required") };
    rows.iter_mut().find_map(|(name,value)| (name==key).then_some(value)).unwrap()
}
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(rows) = root else { panic!("object required") };
    if let Some((_,slot)) = rows.iter_mut().find(|(name,_)| name==key) { *slot=value; }
    else { rows.push((key.to_owned(),value)); }
}
fn remove(root: &mut J, key: &str) {
    let J::Object(rows) = root else { panic!("object required") };
    rows.retain(|(name,_)| name!=key);
}
fn n(value:f64)->J { J::Number { value, raw:value.to_string() } }
fn quote(text:&str)->String {
    let mut out=String::from("\"");
    for ch in text.chars() {
        match ch {
            '"'=>out.push_str("\\\""), '\\'=>out.push_str("\\\\"),
            '\n'=>out.push_str("\\n"), '\r'=>out.push_str("\\r"), '\t'=>out.push_str("\\t"),
            c if c<'\u{20}'=>out.push_str(&format!("\\u{:04x}",c as u32)), c=>out.push(c),
        }
    }
    out.push('"');out
}
fn encode(value:&J)->String {
    match value {
        J::Null=>"null".into(), J::Bool(v)=>v.to_string(), J::Number {raw,..}=>raw.clone(),
        J::Str(s)=>quote(s),
        J::Array(rows)=>format!("[{}]",rows.iter().map(encode).collect::<Vec<_>>().join(",")),
        J::Object(rows)=>format!("{{{}}}",rows.iter().map(|(k,v)|format!("{}:{}",quote(k),encode(v))).collect::<Vec<_>>().join(",")),
    }
}
fn request(text:&str,speed:f64,gradient:bool)->J {
    let mut root=J::parse(text).unwrap();
    put(member(member(&mut root,"hydraulics"),"fan"),"speed_ratio",n(speed));
    put(member(&mut root,"objective"),"gradient",J::Bool(gradient));
    root
}
fn invoke(root:&J)->Output {
    static NEXT:AtomicUsize=AtomicUsize::new(0);
    let dir:PathBuf=std::env::temp_dir().join(format!("frankensim-fan-gradient-{}-{}",
        std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    fs::create_dir_all(&dir).unwrap();
    let path=dir.join("request.json");fs::write(&path,encode(root)).unwrap();
    Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network"]).arg(path).output().unwrap()
}
fn run(root:&J)->J {
    let out=invoke(root);
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    J::parse(std::str::from_utf8(&out.stdout).unwrap()).unwrap()
}
fn value(root:&J,path:&[&str])->f64 {root.path(path).and_then(J::as_f64).unwrap()}
fn gradient(root:&J)->f64 {value(root,&["fan_speed_sensitivity","dobjective_dlog_speed_ratio_k"])}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a:.15e} != {b:.15e}");}

#[test]
fn total_fan_gradient_matches_perturbed_operating_points_with_and_without_correlations() {
    let eps=1e-4_f64;
    for (text,correlated) in [(DECLARED,false),(CORRELATED,true)] {
        for speed in [0.7,1.0,1.5] {
            let nominal=run(&request(text,speed,true));
            let plus=run(&request(text,speed*eps.exp(),false));
            let minus=run(&request(text,speed*(-eps).exp(),false));
            let fd=(value(&plus,&["objective","value_k"])-value(&minus,&["objective","value_k"]))/(2.0*eps);
            close(gradient(&nominal),fd,8e-5);
            let capacity=value(&nominal,&["fan_speed_sensitivity","capacity_contribution_k"]);
            let convection=value(&nominal,&["fan_speed_sensitivity","convection_contribution_k"]);
            close(capacity+convection,gradient(&nominal),1e-12);
            assert!(capacity.abs()>1e-3);
            if correlated {assert!(convection.abs()>1e-3);} else {close(convection,0.0,1e-12);}
            close(value(&nominal,&["fan","flow_m3_s"]),0.004*speed,1e-9);
            close(value(&nominal,&["source_w"]),1.0,1e-8);
            close(value(&nominal,&["robin_out_w"]),1.0,1e-7);
            assert_eq!(plus.get("fan_speed_sensitivity"),Some(&J::Null));
        }
    }
    // Independent dense P1 FEM and analytic air elimination, not the same
    // partitioned solve called a second time.
    let nominal=run(&request(CORRELATED,1.0,true));
    close(value(&nominal,&["objective","value_k"]),302.4575930842692,2e-5);
    close(gradient(&nominal),-0.6003053301825904,2e-5);
}

#[test]
fn nonlinear_contact_speed_gradient_retains_material_and_interface_feedback() {
    let mut base=request(PULSE,1.0,true);remove(&mut base,"transient");
    let nominal=run(&base);
    let eps=1e-4_f64;
    put(member(&mut base,"objective"),"gradient",J::Bool(false));
    put(member(member(&mut base,"hydraulics"),"fan"),"speed_ratio",n(eps.exp()));
    let plus=run(&base);
    put(member(member(&mut base,"hydraulics"),"fan"),"speed_ratio",n((-eps).exp()));
    let minus=run(&base);
    close(gradient(&nominal),(value(&plus,&["objective","value_k"])
        -value(&minus,&["objective","value_k"]))/(2.0*eps),5e-4);
    assert!(nominal.path(&["contact_sensitivities","rows"]).unwrap().as_array().unwrap().len()==1);
    close(value(&nominal,&["fan_speed_sensitivity","convection_contribution_k"]),0.0,1e-12);
    close(value(&nominal,&["source_w"]),20.0,1e-8);
    close(value(&nominal,&["robin_out_w"]),20.0,1e-7);
}

#[test]
fn reversing_an_edge_does_not_reverse_the_fan_derivative() {
    let root=request(CORRELATED,1.0,true);
    let nominal=run(&root);
    let mut reversed=root;
    let J::Array(branches)=member(member(&mut reversed,"hydraulics"),"branches") else {panic!()};
    put(&mut branches[0],"from",n(1.0));put(&mut branches[0],"to",n(0.0));
    let other=run(&reversed);
    close(gradient(&other),gradient(&nominal),1e-7);
    close(value(&other,&["objective","value_k"]),value(&nominal,&["objective","value_k"]),1e-7);
    assert!(other.get("branches").unwrap().as_array().unwrap()[0].f64_field("flow_m3_s").unwrap()<0.0);
}

#[test]
fn no_fan_or_adjoint_failure_cannot_publish_a_total_speed_gradient() {
    let mut base=request(DECLARED,1.0,true);
    let hyd=member(&mut base,"hydraulics");remove(hyd,"fan");
    put(hyd,"boundaries",J::parse(r#"[{"node":0,"pressure_pa":20,"temperature_k":300},{"node":2,"pressure_pa":0}]"#).unwrap());
    assert_eq!(run(&base).get("fan_speed_sensitivity"),Some(&J::Null));
    let mut bad=request(CORRELATED,1.0,true);
    put(member(&mut bad,"budgets"),"derivative_iterations",n(1.0));
    let out=invoke(&bad);assert!(!out.status.success());assert!(out.stdout.is_empty());
}

#[test]
fn gradient_guided_sizing_returns_an_evaluated_passing_bracket() {
    let mut base=request(CORRELATED,1.0,true);
    put(&mut base,"fan_speed_design",J::parse(r#"{
        "min_speed_ratio":0.5,"max_speed_ratio":2,"temperature_limit_k":302.4,
        "speed_ratio_tolerance":0.0001,"temperature_tolerance_k":0.00001,"max_evaluations":64
    }"#).unwrap());
    let result=run(&base);
    assert!(value(&result,&["fan_speed_design","newton_trials"])>0.0);
    let selected=value(&result,&["fan_speed_design","selected_speed_ratio"]);
    let temperature=value(&result,&["objective","value_k"]);
    assert!(temperature<=302.4 && 302.4-temperature<=1e-5);
    assert!(value(&result,&["fan_speed_design","speed_bracket_width"])<=1e-4);
    assert!(value(&result,&["fan_speed_design","failed_lower","temperature_k"])>302.4);
    let check=run(&request(CORRELATED,selected,false));
    close(value(&check,&["objective","value_k"]),temperature,1e-8);
    put(member(&mut base,"objective"),"gradient",J::Bool(false));
    let bisected=run(&base);
    close(value(&bisected,&["fan_speed_design","newton_trials"]),0.0,1e-12);
    close(value(&bisected,&["fan_speed_design","selected_speed_ratio"]),selected,1e-4);
    put(member(&mut base,"fan_speed_design"),"max_evaluations",n(2.0));
    let exhausted=invoke(&base);
    assert_eq!(exhausted.status.code(),Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(exhausted.stdout.is_empty());
}
