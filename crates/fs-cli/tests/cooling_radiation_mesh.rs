//! Actual-command regressions for radiation-preserving mesh studies. Tolerances
//! in the wiring tests are deliberately loose; none asserts a continuum bound.
#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/radiative-contact-hotspot.json"));
fn number(value: f64) -> J { J::Number { value, raw: value.to_string() } }
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object required") };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object required") };
    if let Some((_, slot)) = fields.iter_mut().find(|(name, _)| name == key) { *slot = value; }
    else { fields.push((key.to_string(), value)); }
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object required") };
    fields.retain(|(name, _)| name != key);
}
fn encode(value: &J) -> String {
    match value {
        J::Null => "null".into(), J::Bool(b) => b.to_string(), J::Number { raw, .. } => raw.clone(),
        J::Str(s) => format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"")
            .replace('\n', "\\n").replace('\r', "\\r").replace('\t', "\\t")),
        J::Array(a) => format!("[{}]", a.iter().map(encode).collect::<Vec<_>>().join(",")),
        J::Object(o) => format!("{{{}}}", o.iter().map(|(k,v)|
            format!("{}:{}", encode(&J::Str(k.clone())), encode(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(input: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(encode(input).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(input: &J) -> J {
    let result = output(input);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    J::parse(std::str::from_utf8(&result.stdout).unwrap()).unwrap()
}
fn n(value: &J, key: &str) -> f64 { value.f64_field(key).unwrap() }
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a} versus {b}"); }
fn input(adaptive: bool) -> J {
    let mut root = J::parse(BASE).unwrap();
    let mut policy = J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":100,"max_vertices":20000,"max_tetrahedra":100000}"#).unwrap();
    if adaptive {
        put(&mut policy, "strategy", J::Str("goal-recovery".into()));
        put(&mut policy, "marking_fraction", number(0.5));
        put(&mut policy, "max_refinements", number(3.0));
    }
    put(&mut root, "mesh_convergence", policy);
    root
}
fn verify(input: &J, result: &J, adaptive: bool) {
    let study = result.get("mesh_convergence").unwrap();
    let history = study.get("history").unwrap().as_array().unwrap();
    assert!(history.len() >= 3);
    assert_eq!(history.last().unwrap().str_field("arrived_by"), Some("uniform"));
    near(history.iter().map(|r| n(r,"solid_solves")).sum(), n(study,"total_solid_solves"), 0.0);
    for row in history { near(n(row,"source_w"), 20.0, 1e-7); }
    let radiation = result.get("radiation").unwrap();
    near(n(radiation,"radiative_out_w") + n(radiation,"convective_out_w"), 20.0, 1e-7);
    assert!(n(radiation,"radiative_out_w") > 1.0, "radiation must not be omitted");
    assert!(n(radiation,"energy_residual_w").abs() <= 1e-7);
    near(n(history.last().unwrap(),"solid_solves"), n(radiation,"solid_solves"), 0.0);
    assert!(n(radiation,"solid_solves") > n(result,"coupling_iterations"));
    let resolved = study.get("resolved_request").unwrap();
    assert_eq!(resolved.get("radiation"), input.get("radiation"), "keep the physical patch partition");
    assert_eq!(resolved.path(&["solid","materials"]), input.path(&["solid","materials"]));
    assert!(resolved.path(&["solid","component_power"]).is_none());
    assert!(resolved.path(&["solid","nodal_source_w_m3"]).is_some());
    assert!(resolved.get("mesh_convergence").is_none());
    let replay = run(resolved);
    assert_eq!(result.get("solid_temperatures_k"), replay.get("solid_temperatures_k"));
    assert_eq!(result.get("objective"), replay.get("objective"));
    assert_eq!(result.get("contacts"), replay.get("contacts"));
    assert_eq!(radiation.get("surfaces"), replay.path(&["radiation","surfaces"]));
    assert_eq!(result.get("adjoint_residual"), Some(&J::Null));
    assert_eq!(radiation.get("adjoint"), Some(&J::Null));
    if adaptive {
        assert_eq!(study.get("global_confirmation"), Some(&J::Bool(true)));
        assert!(n(study,"total_adjoint_sweeps") > 0.0);
        assert!(history.iter().any(|row| row.str_field("arrived_by") == Some("marked-edge-stars")));
        assert_eq!(n(radiation,"reconstruction_solid_solves"), 1.0);
        for row in history {
            assert!(n(row.get("marking").unwrap(),"adjoint_sweeps") > 0.0);
        }
    } else {
        assert_eq!(n(study,"total_adjoint_sweeps"), 0.0);
        assert_eq!(n(radiation,"reconstruction_solid_solves"), 0.0);
        assert_eq!(history.len(), 3);
    }
}

#[test]
fn uniform_refinement_replays_the_nonlinear_radiating_contact_problem() {
    let input = input(false); let result = run(&input); verify(&input,&result,false);
}

#[test]
fn local_refinement_uses_total_radiative_feedback_and_requires_global_confirmation() {
    let input = input(true); let result = run(&input); verify(&input,&result,true);
}

#[test]
fn insufficient_global_confirmation_or_radiation_budgets_never_publish_success() {
    let base = input(true);
    let mut local_only = base.clone();
    put(member(&mut local_only,"mesh_convergence"),"max_refinements",number(2.0));
    let result = output(&local_only); assert_eq!(result.status.code(),Some(6)); assert!(result.stdout.is_empty());
    for (section,key,value) in [("radiation","max_iterations",1.0),
        ("budgets","derivative_iterations",1.0),("mesh_convergence","max_tetrahedra",12.0),
        ("budgets","wall_seconds",1e-12)] {
        let mut bad = base.clone(); put(member(&mut bad,section),key,number(value));
        let result = output(&bad);
        assert!(!result.status.success(), "unexpected success for {section}.{key}");
        assert!(result.stdout.is_empty(), "partial output for {section}.{key}");
    }
}

#[test]
fn hot_reservoirs_remain_signed_energy_sources_after_refinement() {
    let mut input = input(false);
    let J::Array(patches) = member(member(&mut input,"radiation"),"surfaces") else { panic!() };
    for patch in patches { put(patch,"ambient_temperature_k",number(350.0)); }
    let result = run(&input);
    let radiation = result.get("radiation").unwrap();
    assert!(n(radiation,"radiative_out_w") < 0.0);
    near(n(radiation,"radiative_out_w") + n(radiation,"convective_out_w"),20.0,1e-7);
    assert!(n(radiation,"convective_out_w") > 20.0);
}

#[test]
fn changing_patch_order_does_not_change_the_adaptive_mesh_or_field() {
    let base = input(true); let first = output(&base);
    assert!(first.status.success(), "{}",String::from_utf8_lossy(&first.stderr));
    let mut reordered = base.clone();
    let J::Array(patches) = member(member(&mut reordered,"radiation"),"surfaces") else { panic!() };
    patches.reverse();
    let second = run(&reordered);
    let first = J::parse(std::str::from_utf8(&first.stdout).unwrap()).unwrap();
    assert_eq!(first.get("solid_temperatures_k"),second.get("solid_temperatures_k"));
    assert_eq!(first.path(&["mesh_convergence","history"]),second.path(&["mesh_convergence","history"]));
}

#[test]
fn uniform_studies_do_not_require_derivatives_and_vanishing_emissivity_recovers_convection() {
    let mut weak = input(false);
    put(member(&mut weak,"budgets"),"derivative_iterations",number(1.0));
    let J::Array(patches) = member(member(&mut weak,"radiation"),"surfaces") else { panic!() };
    for patch in patches { put(patch,"emissivity",number(1e-10)); }
    let radiating = run(&weak);
    let mut ordinary = weak.clone(); remove(&mut ordinary,"radiation");
    let plain = run(&ordinary);
    assert_eq!(n(radiating.get("mesh_convergence").unwrap(),"total_adjoint_sweeps"),0.0);
    for (a,b) in radiating.get("solid_temperatures_k").unwrap().as_array().unwrap().iter()
        .zip(plain.get("solid_temperatures_k").unwrap().as_array().unwrap()) {
        near(a.as_f64().unwrap(),b.as_f64().unwrap(),2e-6);
    }
}
