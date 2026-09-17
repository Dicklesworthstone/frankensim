#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const SLAB: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/mixed-slab.json"));
const CONTACT: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/radiative-contact-hotspot.json"));
const SIGMA: f64 = 5.670_374_419e-8;

fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object required") };
    if let Some((_, slot)) = fields.iter_mut().find(|(name, _)| name == key) { *slot = value; }
    else { fields.push((key.to_string(), value)); }
}
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object required") };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object required") };
    fields.retain(|(name, _)| name != key);
}
fn number(value: f64) -> J { J::Number { value, raw: value.to_string() } }
fn text(value: &J) -> String {
    match value {
        J::Null => "null".into(), J::Bool(b) => b.to_string(),
        J::Number { raw, .. } => raw.clone(),
        J::Str(s) => format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")),
        J::Array(a) => format!("[{}]", a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(o) => format!("{{{}}}", o.iter().map(|(k,v)|format!("{}:{}",text(&J::Str(k.clone())),text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(root: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(root).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(root: &J) -> J {
    let result = output(root);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    J::parse(&String::from_utf8(result.stdout).unwrap()).unwrap()
}
fn n(root: &J, key: &str) -> f64 { root.f64_field(key).unwrap() }
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a} versus {b}"); }
fn policy(ambient: f64) -> J {
    J::parse(&format!(r#"{{"max_iterations":128,"temperature_tolerance_k":1e-9,"relaxation":0.5,"surfaces":[{{"surface":"first-face","emissivity":0.85,"ambient_temperature_k":{ambient},"source":"analytic slab fixture"}},{{"surface":"last-face","emissivity":0.6,"ambient_temperature_k":{},"source":"analytic slab fixture"}}]}}"#, ambient+10.0)).unwrap()
}
fn slab(ambient: f64) -> J {
    let mut root = J::parse(SLAB).unwrap();
    put(member(&mut root, "objective"), "gradient", J::Bool(false));
    put(&mut root, "radiation", policy(ambient));
    root
}

/// Eliminate the two air streams and solve just the two slab-face temperatures.
/// This reference never assembles FEM matrices or invokes the product solver.
fn reference(ambient: f64) -> [f64; 2] {
    let capacity = [0.003*1.2*1007.0, 0.004*1.2*1007.0];
    let eps = [1.0-(-0.5_f64/capacity[0]).exp(), 1.0-(-0.8_f64/capacity[1]).exp()];
    let g = [capacity[0]*eps[0], capacity[1]*eps[1]];
    let mut w = [320.0_f64,320.0_f64];
    for _ in 0..64 {
        let mixing = 0.75*(330.0+eps[0]*(w[0]-330.0))+0.25*290.0;
        let rad = [0.85*SIGMA*0.01*(w[0].powi(4)-ambient.powi(4)),
            0.6*SIGMA*0.01*(w[1].powi(4)-(ambient+10.0).powi(4))];
        // k*A/L = 10*.01/.05 = 2 W/K.
        let r = [2.0*(w[0]-w[1])+g[0]*(w[0]-330.0)+rad[0],
            2.0*(w[1]-w[0])+g[1]*(w[1]-mixing)+rad[1]];
        if r[0].abs().max(r[1].abs()) < 1e-11 { return w; }
        let a = 2.0+g[0]+4.0*0.85*SIGMA*0.01*w[0].powi(3);
        let b = -2.0;
        let c = -2.0-g[1]*0.75*eps[0];
        let d = 2.0+g[1]+4.0*0.6*SIGMA*0.01*w[1].powi(3);
        let det = a*d-b*c;
        w[0] -= (d*r[0]-b*r[1])/det;
        w[1] -= (a*r[1]-c*r[0])/det;
    }
    panic!("independent slab equations did not converge")
}

#[test]
fn hot_and_cold_surroundings_match_the_independent_slab_solution() {
    for ambient in [270.0,300.0,400.0] {
        let input = slab(ambient); let result = run(&input); let expected = reference(ambient);
        let positions = input.path(&["solid","vertices_m"]).unwrap().as_array().unwrap();
        for (value,p) in result.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().zip(positions) {
            let x = p.as_array().unwrap()[0].as_f64().unwrap();
            near(value.as_f64().unwrap(), expected[0]+(expected[1]-expected[0])*x/0.05, 2e-6);
        }
        let report = result.get("radiation").unwrap();
        near(n(report,"radiative_out_w")+n(report,"convective_out_w"), 0.0, 1e-7);
        assert!(n(report,"energy_residual_w").abs() <= 1e-7);
        if ambient == 270.0 { assert!(n(report,"radiative_out_w")>0.0); }
        if ambient == 400.0 { assert!(n(report,"radiative_out_w")<0.0); }
        for row in report.get("surfaces").unwrap().as_array().unwrap() {
            let flux = n(row,"emissivity")*SIGMA*0.01*
                (n(row,"mean_temperature_k").powi(4)-n(row,"ambient_temperature_k").powi(4));
            near(n(row,"nonlinear_heat_w"),flux,1e-9);
            near(n(row,"applied_heat_w"),flux,1e-7);
        }
    }
}

#[test]
fn nonlinear_contact_source_closes_against_air_plus_radiation_not_air_alone() {
    let input = J::parse(CONTACT).unwrap(); let result = run(&input);
    let report = result.get("radiation").unwrap();
    let q = n(report,"radiative_out_w"); let air = n(report,"convective_out_w");
    assert!(q>1.0, "fixture must expose counting radiation as air heat");
    near(air+q,20.0,1e-7);near(n(&result,"robin_out_w"),20.0,1e-7);
    let walls=result.get("walls").unwrap().as_array().unwrap();
    near(walls.iter().map(|w|n(w,"outward_heat_w")).sum(),air,1e-8);
    let exhaust=result.get("branches").unwrap().as_array().unwrap().iter()
        .find(|b|b.str_field("name")==Some("mixed-stream")).unwrap();
    near(n(exhaust,"flow_m3_s")*1.2*1007.0*(n(exhaust,"outlet_k")-300.0),air,3e-7);
    assert_eq!(result.get("contact_sensitivities"),Some(&J::Null));
    assert!(!result.get("contacts").unwrap().as_array().unwrap().is_empty());
    assert!(n(report,"solid_solves")>n(&result,"coupling_iterations"));
    let mut ordinary=input.clone();remove(&mut ordinary,"radiation");
    let plain=run(&ordinary);
    assert!(n(result.get("objective").unwrap(),"value_k")<n(plain.get("objective").unwrap(),"value_k"));
}

#[test]
fn vanishing_emissivity_recovers_the_unmodified_cooling_problem() {
    let mut radiating=J::parse(CONTACT).unwrap();
    let J::Array(rows)=member(member(&mut radiating,"radiation"),"surfaces") else {panic!()};
    for row in rows {put(row,"emissivity",number(1e-10));}
    let mut plain=radiating.clone();remove(&mut plain,"radiation");
    let a=run(&radiating);let b=run(&plain);
    for (x,y) in a.get("solid_temperatures_k").unwrap().as_array().unwrap().iter()
        .zip(b.get("solid_temperatures_k").unwrap().as_array().unwrap()) {
        near(x.as_f64().unwrap(),y.as_f64().unwrap(),2e-6);
    }
    assert!(b.get("radiation").is_none());
}

#[test]
fn replay_and_patch_declaration_order_do_not_change_the_accepted_result() {
    let mut input=J::parse(CONTACT).unwrap();let a=output(&input);
    assert!(a.status.success(),"{}",String::from_utf8_lossy(&a.stderr));
    assert_eq!(a.stdout,output(&input).stdout);
    let J::Array(rows)=member(member(&mut input,"radiation"),"surfaces") else {panic!()};rows.reverse();
    let b=output(&input);assert!(b.status.success());assert_eq!(a.stdout,b.stdout);
}

#[test]
fn exhausted_or_unsupported_radiation_never_publishes_partial_cooling() {
    let input=slab(270.0);
    let mut exhausted=input.clone();put(member(&mut exhausted,"radiation"),"max_iterations",number(1.0));
    let result=output(&exhausted);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
    for (field,value) in [("emissivity",number(0.0)),("emissivity",number(1.01)),
        ("ambient_temperature_k",number(0.0)),("surface",J::Str("missing".into()))] {
        let mut bad=input.clone();
        let J::Array(rows)=member(member(&mut bad,"radiation"),"surfaces") else {panic!()};
        put(&mut rows[0],field,value);
        let result=output(&bad);assert!(!result.status.success());assert!(result.stdout.is_empty());
    }
    let mut derivative=input.clone();put(member(&mut derivative,"objective"),"gradient",J::Bool(true));
    put(member(&mut derivative,"budgets"),"derivative_iterations",number(1.0));
    let result=output(&derivative);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
    let mut cancelled=input.clone();put(member(&mut cancelled,"budgets"),"wall_seconds",number(1e-12));
    let result=output(&cancelled);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
}

#[path = "cooling_radiation/adjoint.rs"]
mod adjoint;
