#![cfg(unix)]
#[allow(dead_code)]
#[path="../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command,Output,Stdio};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static FIXTURE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/mesh-convergence-contact-hotspot.json"))));

fn output(text:&str)->Output {
    let mut child=Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(text:&str)->J {
    let result=output(text);
    assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stderr));
    J::parse(&String::from_utf8(result.stdout).unwrap()).unwrap()
}
fn n(value:&J,key:&str)->f64{value.f64_field(key).unwrap()}
fn fields(value:&mut J)->&mut Vec<(String,J)>{let J::Object(fields)=value else{panic!("object")};fields}
fn member<'a>(value:&'a mut J,key:&str)->&'a mut J {
    fields(value).iter_mut().find_map(|(k,v)|(k==key).then_some(v)).unwrap()
}
fn quote(s:&str)->String {
    use std::fmt::Write as _;
    let mut out=String::from("\"");
    for c in s.chars(){match c{
        '"'=>out.push_str("\\\""),'\\'=>out.push_str("\\\\"),'\n'=>out.push_str("\\n"),
        '\r'=>out.push_str("\\r"),'\t'=>out.push_str("\\t"),
        c if c<'\u{20}'=>{write!(&mut out,"\\u{:04x}",c as u32).unwrap();},c=>out.push(c),
    }}out.push('"');out
}
fn encode(value:&J)->String {
    match value {
        J::Null=>"null".into(),J::Bool(v)=>v.to_string(),J::Number{raw,..}=>raw.clone(),J::Str(s)=>quote(s),
        J::Array(values)=>format!("[{}]",values.iter().map(encode).collect::<Vec<_>>().join(",")),
        J::Object(values)=>format!("{{{}}}",values.iter().map(|(k,v)|format!("{}:{}",quote(k),encode(v)))
            .collect::<Vec<_>>().join(",")),
    }
}

#[test]
fn contact_hotspot_study_solves_actual_meshes_and_replays_the_published_field() {
    let result=run(FIXTURE);
    let study=result.get("mesh_convergence").unwrap();
    assert_eq!(study.str_field("status"),Some("successive-mesh-tolerance-met"));
    let rows=study.get("history").unwrap().as_array().unwrap();
    assert!(rows.len()>=3);
    assert_eq!(n(study,"meshes_solved"),rows.len() as f64);
    assert!(n(study,"achieved_change_k")<=0.02);
    assert_eq!(n(rows.last().unwrap(),"consecutive_passes"),2.);
    for (level,row) in rows.iter().enumerate() {
        assert_eq!(n(row,"tetrahedra"),(12*8_usize.pow(level as u32)) as f64);
        assert!((n(row,"source_w")-1.).abs()<1e-7);
    }
    let resolved=study.get("resolved_request").unwrap();
    assert!(resolved.get("mesh_convergence").is_none());
    let replay=run(&encode(resolved));
    for key in ["objective","solid_temperatures_k","contacts","walls","branches","source_w"] {
        assert_eq!(result.get(key),replay.get(key),"replay changed {key}");
    }
    let solid=resolved.get("solid").unwrap();
    assert!(solid.get("component_power").is_none());
    let density=solid.get("nodal_source_w_m3").unwrap().as_array().unwrap();
    assert!((density[4].as_f64().unwrap()-48000.).abs()<1e-8);
    assert_eq!(density.len(),result.get("solid_temperatures_k").unwrap().as_array().unwrap().len());
    let original=J::parse(FIXTURE).unwrap();
    assert_eq!(solid.get("materials"),original.path(&["solid","materials"]));
    let contacts=result.get("contacts").unwrap().as_array().unwrap();
    assert!((n(&contacts[0],"area_m2")-0.01).abs()<1e-12);
    assert_eq!(n(&contacts[0],"face_pairs"),(2*4_usize.pow((rows.len()-1) as u32)) as f64);
}

#[test]
fn reusing_the_original_hot_vertex_is_detectably_a_different_physical_source() {
    let result=run(FIXTURE);
    let mut wrong=result.path(&["mesh_convergence","resolved_request"]).unwrap().clone();
    let original=J::parse(FIXTURE).unwrap();
    let component=original.path(&["solid","component_power"]).unwrap().clone();
    let solid=member(&mut wrong,"solid");
    fields(solid).retain(|(key,_)|key!="nodal_source_w_m3");
    fields(solid).push(("component_power".into(),component));
    let wrong=run(&encode(&wrong));
    assert!((n(&wrong,"source_w")-n(&result,"source_w")).abs()<1e-7);
    assert!(n(wrong.get("objective").unwrap(),"value_k")
        -n(result.get("objective").unwrap(),"value_k")>0.5,
        "equal total watts must not conceal the changed spatial source");
}

#[test]
fn refinement_count_mesh_size_and_wall_exhaustion_publish_no_convergence() {
    for (from,to) in [
        ("\"max_tetrahedra\":100000","\"max_tetrahedra\":12"),
        ("\"max_vertices\":20000","\"max_vertices\":16"),
        ("\"temperature_tolerance_k\":0.02","\"temperature_tolerance_k\":1e-14"),
        ("\"wall_seconds\":120","\"wall_seconds\":1e-12"),
    ] {
        assert!(FIXTURE.contains(from));
        let result=output(&FIXTURE.replace(from,to));
        assert_eq!(result.status.code(),Some(6),"{}",String::from_utf8_lossy(&result.stderr));
        assert!(result.stdout.is_empty(),"a budget is not mesh convergence");
        assert!(!result.stderr.is_empty());
    }
}

#[test]
fn unsupported_policy_or_ambiguous_source_never_falls_back_to_a_plain_solve() {
    for (from,to) in [
        ("\"consecutive_passes\":2","\"consecutive_passes\":1"),
        ("\"gradient\":false","\"gradient\":true"),
        ("\"mesh_convergence\":{","\"transient\":{},\"mesh_convergence\":{"),
        ("\"component_power\":{","\"nodal_source_w_m3\":[1],\"component_power\":{"),
    ] {
        assert!(FIXTURE.contains(from));
        let result=output(&FIXTURE.replace(from,to));
        assert!(!result.status.success());assert!(result.stdout.is_empty());assert!(!result.stderr.is_empty());
    }
}
