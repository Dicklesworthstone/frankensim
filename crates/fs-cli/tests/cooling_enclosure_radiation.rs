//! These tests invoke the real command. The reference eliminates the two
//! radiosities and air references analytically; it does not repeat the solver.
#![cfg(unix)]
#[allow(dead_code)]
#[path="../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command,Output,Stdio};

const BASE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enclosure-radiation-gap.json"));
fn num(value:f64)->J {J::Number {value,raw:value.to_string()}}
fn member<'a>(root:&'a mut J,key:&str)->&'a mut J {
    let J::Object(fields)=root else {panic!("object required")};
    &mut fields.iter_mut().find(|(k,_)|k==key).unwrap().1
}
fn put(root:&mut J,key:&str,value:J) {
    let J::Object(fields)=root else {panic!("object required")};
    if let Some((_,v))=fields.iter_mut().find(|(k,_)|k==key){*v=value;}
    else {fields.push((key.into(),value));}
}
fn remove(root:&mut J,key:&str) {
    let J::Object(fields)=root else {panic!("object required")};fields.retain(|(k,_)|k!=key);
}
fn rows(root:&mut J)->&mut Vec<J> {let J::Array(r)=root else {panic!("array required")};r}
fn enclosure(root:&mut J)->&mut J {member(member(root,"radiation"),"enclosure")}
fn text(root:&J)->String {
    match root {
        J::Null=>"null".into(),J::Bool(v)=>v.to_string(),J::Number{raw,..}=>raw.clone(),
        J::Str(s)=>format!("\"{}\"",s.replace('\\',"\\\\").replace('"',"\\\"").replace('\n',"\\n")),
        J::Array(a)=>format!("[{}]",a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(a)=>format!("{{{}}}",a.iter().map(|(k,v)|format!("{}:{}",text(&J::Str(k.clone())),text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(input:&J)->Output {
    let mut child=Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(input).as_bytes()).unwrap();child.wait_with_output().unwrap()
}
fn success(output:&Output)->J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run(input:&J)->J {success(&output(input))}
fn n(root:&J,key:&str)->f64 {root.f64_field(key).unwrap()}
fn near(a:f64,b:f64,t:f64) {assert!((a-b).abs()<t,"{a} versus {b}");}
fn patch<'a>(result:&'a J,name:&str)->&'a J {
    result.path(&["radiation","surfaces"]).unwrap().as_array().unwrap().iter()
        .find(|p|p.str_field("surface")==Some(name)).unwrap()
}
fn finish(input:&mut J,emitter:f64,receiver:f64) {
    for (row,value) in rows(member(enclosure(input),"surfaces")).iter_mut().zip([emitter,receiver]) {
        put(row,"emissivity",num(value));
    }
}

fn reference(e0:f64,e1:f64)->(f64,f64,f64) {
    // Fan/system intersection is Q=0.004 m3/s. The heated branch carries
    // 0.003 m3/s; the cold bypass adds 0.001 before the downstream receiver.
    let c0=0.003*1.2*1007.0_f64;let c1=0.004*1.2*1007.0_f64;
    let g0=c0*(1.0-(-0.1/c0).exp());let g1=c1*(1.0-(-0.25/c1).exp());
    let temperatures=|q:f64| (300.0+(5.0-q)/g0,300.0+(5.0-q)/c1+q/g1);
    let mut lo=0.0;let mut hi=5.0;
    for _ in 0..80 {
        let q=0.5*(lo+hi);let (a,b)=temperatures(q);
        let emitted=0.01*5.670374419e-8*(a-b)*(a+b)*(a*a+b*b)/(1.0/e0+1.0/e1-1.0);
        if q<emitted {lo=q;} else {hi=q;}
    }
    let q=0.5*(lo+hi);let (a,b)=temperatures(q);(a,b,q)
}

fn check_reference(input:&J,e0:f64,e1:f64)->J {
    let result=run(input);let (a,b,q)=reference(e0,e1);
    assert_eq!(result.path(&["radiation","model"]).unwrap().as_str(),Some("closed-gray-diffuse-enclosure"));
    near(n(patch(&result,"emitter"),"mean_temperature_k"),a,2e-6);
    near(n(patch(&result,"receiver"),"mean_temperature_k"),b,2e-6);
    near(n(patch(&result,"emitter"),"outward_heat_w"),q,1e-7);
    near(n(patch(&result,"receiver"),"outward_heat_w"),-q,1e-7);
    near(n(result.get("radiation").unwrap(),"radiative_out_w"),0.0,1e-7);
    let air:f64=result.get("walls").unwrap().as_array().unwrap().iter().map(|r|n(r,"outward_heat_w")).sum();
    near(air,5.0,1e-7);near(n(&result,"robin_out_w"),5.0,1e-7);
    assert_eq!(result.get("dobjective_dinlet_k"),Some(&J::Null));
    result
}

#[test]
fn reflected_heat_crosses_the_gap_but_is_not_an_external_energy_sink() {
    let input=J::parse(BASE).unwrap();let result=check_reference(&input,0.8,0.6);
    near(n(patch(&result,"emitter"),"outward_heat_w"),1.263461440576686,1e-7);
    let mut black=input.clone();finish(&mut black,1.0,1.0);
    let black=check_reference(&black,1.0,1.0);
    assert!(n(patch(&black,"emitter"),"outward_heat_w")>n(patch(&result,"emitter"),"outward_heat_w"));
    let mut absent=input.clone();remove(&mut absent,"radiation");let absent=run(&absent);
    assert!(n(absent.get("objective").unwrap(),"value_k")>n(result.get("objective").unwrap(),"value_k")+5.0);
    // Near mirrors approach zero exchange, not blackbody exchange.
    let mut mirrors=input;finish(&mut mirrors,1e-5,1e-5);let mirrors=check_reference(&mirrors,1e-5,1e-5);
    near(n(mirrors.get("objective").unwrap(),"value_k"),n(absent.get("objective").unwrap(),"value_k"),1e-3);
}

#[test]
fn nonlinear_materials_keep_the_same_surface_balance_without_being_frozen() {
    let input=J::parse(BASE).unwrap();let baseline=run(&input);let mut nonlinear=input;
    for row in rows(member(member(&mut nonlinear,"solid"),"materials")) {
        remove(row,"conductivity_w_m_k");
        put(row,"conductivity_curve",J::parse(r#"{"temperature_k":[250,450],"conductivity_w_m_k":[0.5,4.5]}"#).unwrap());
    }
    let result=check_reference(&nonlinear,0.8,0.6);
    assert!((n(result.get("objective").unwrap(),"value_k")-n(baseline.get("objective").unwrap(),"value_k")).abs()>0.1);
}

#[test]
fn permuting_the_surface_list_and_both_matrix_axes_replays_exactly() {
    let input=J::parse(BASE).unwrap();let original=output(&input);success(&original);
    let mut reordered=input;
    rows(member(enclosure(&mut reordered),"surfaces")).reverse();
    let matrix=rows(member(enclosure(&mut reordered),"view_factors"));
    matrix.reverse();for row in matrix {rows(row).reverse();}
    let replay=output(&reordered);success(&replay);assert_eq!(original.stdout,replay.stdout);
}

#[test]
fn uniform_mesh_study_preserves_the_closed_model_and_final_field_replay() {
    let mut input=J::parse(BASE).unwrap();
    put(&mut input,"mesh_convergence",J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":2,"max_vertices":10000,"max_tetrahedra":10000}"#).unwrap());
    let result=run(&input);
    let resolved=result.path(&["mesh_convergence","resolved_request"]).unwrap();
    assert_eq!(resolved.get("radiation"),input.get("radiation"));
    let replay=run(resolved);
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(result.get("radiation"),replay.get("radiation"));
    near(n(patch(&result,"emitter"),"outward_heat_w"),reference(0.8,0.6).2,1e-7);
}

#[test]
fn incomplete_matrices_conflicting_models_and_unimplemented_gradients_refuse() {
    let input=J::parse(BASE).unwrap();let mut cases=Vec::new();
    let mut unclosed=input.clone();put(enclosure(&mut unclosed),"view_factors",J::parse("[[0,0.9],[1,0]]").unwrap());cases.push(unclosed);
    let mut nonreciprocal=input.clone();put(enclosure(&mut nonreciprocal),"view_factors",J::parse("[[0,1],[0.5,0.5]]").unwrap());cases.push(nonreciprocal);
    let mut both=input.clone();put(member(&mut both,"radiation"),"surfaces",J::Array(Vec::new()));cases.push(both);
    let mut gradient=input.clone();put(member(&mut gradient,"objective"),"gradient",J::Bool(true));cases.push(gradient);
    let mut marking=input.clone();put(&mut marking,"mesh_convergence",J::parse(r#"{"strategy":"goal-recovery","marking_fraction":0.5,"max_refinements":3,"consecutive_passes":2,"temperature_tolerance_k":2,"max_vertices":10000,"max_tetrahedra":10000}"#).unwrap());cases.push(marking);
    for input in cases {let output=output(&input);assert!(!output.status.success());assert!(output.stdout.is_empty());}
    let mut exhausted=input;
    put(member(&mut exhausted,"radiation"),"max_iterations",num(1.0));
    let output=output(&exhausted);assert_eq!(output.status.code(),Some(6));assert!(output.stdout.is_empty());
}

#[path="cooling_enclosure_radiation/uq.rs"]
mod uq;
