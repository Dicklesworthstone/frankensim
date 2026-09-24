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
    "/../../examples/cooling-network/adaptive-contact-hotspot.json"))));
fn output(text:&str)->Output {
    let mut child=Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();child.wait_with_output().unwrap()
}
fn run(text:&str)->J {
    let out=output(text);assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    J::parse(&String::from_utf8(out.stdout).unwrap()).unwrap()
}
fn n(j:&J,key:&str)->f64 {j.f64_field(key).unwrap()}
fn encode(j:&J)->String {
    fn quote(s:&str)->String {
        let mut out=String::from("\"");
        for c in s.chars() {match c {
            '"'=>out.push_str("\\\""),'\\'=>out.push_str("\\\\"),
            '\n'=>out.push_str("\\n"),'\r'=>out.push_str("\\r"),'\t'=>out.push_str("\\t"),
            c=>out.push(c),
        }}out.push('"');out
    }
    match j {
        J::Null=>"null".into(),J::Bool(b)=>b.to_string(),J::Number{raw,..}=>raw.clone(),J::Str(s)=>quote(s),
        J::Array(a)=>format!("[{}]",a.iter().map(encode).collect::<Vec<_>>().join(",")),
        J::Object(o)=>format!("{{{}}}",o.iter().map(|(k,v)|format!("{}:{}",quote(k),encode(v))).collect::<Vec<_>>().join(",")),
    }
}
#[test]
fn nonlinear_contact_study_uses_local_cells_then_a_global_check_and_replays_the_field() {
    let result=run(FIXTURE);let study=result.get("mesh_convergence").unwrap();
    assert_eq!(study.str_field("method"),Some("goal-recovery-edge-bisection"));
    assert_eq!(study.get("global_confirmation"),Some(&J::Bool(true)));
    assert!(n(study,"total_adjoint_sweeps")>0.0);
    let history=study.get("history").unwrap().as_array().unwrap();
    assert!(history.len()>=4);
    assert_eq!(history[1].str_field("arrived_by"),Some("marked-edge-stars"));
    assert!(n(&history[1],"tetrahedra")>12.0 && n(&history[1],"tetrahedra")<96.0);
    assert_eq!(history.last().unwrap().str_field("arrived_by"),Some("uniform"));
    assert!(n(history.last().unwrap(),"successive_change_k")<=0.02);
    for row in history {
        assert!((n(row,"source_w")-1.0).abs()<1e-7);
        let marking=row.get("marking").unwrap();
        assert_eq!(marking.get("score_is_error_bound"),Some(&J::Bool(false)));
        assert!(n(marking,"captured_fraction")>=0.5 || n(marking,"normalized_score_sum")==0.0);
    }
    let resolved=study.get("resolved_request").unwrap();
    assert!(resolved.get("mesh_convergence").is_none());
    assert!(resolved.path(&["solid","component_power"]).is_none());
    assert!(resolved.path(&["solid","nodal_source_w_m3"]).is_some());
    let original=J::parse(FIXTURE).unwrap();
    assert_eq!(resolved.path(&["solid","materials"]),original.path(&["solid","materials"]));
    let contact=&result.get("contacts").unwrap().as_array().unwrap()[0];
    assert!((n(contact,"area_m2")-0.01).abs()<1e-12);
    assert_eq!(n(contact,"resistance_m2_k_w"),0.01);
    let replay=run(&encode(resolved));
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(result.get("objective"),replay.get("objective"));
    assert_eq!(result.get("contacts"),replay.get("contacts"));
    assert_eq!(result.get("adjoint_residual"),Some(&J::Null));
}
#[test]
fn local_agreement_cannot_skip_the_budgeted_global_confirmation() {
    let too_short=FIXTURE.replace("\"max_refinements\":12","\"max_refinements\":2")
        .replace("\"temperature_tolerance_k\":0.02","\"temperature_tolerance_k\":100");
    let failure=output(&too_short);
    assert_eq!(failure.status.code(),Some(6));assert!(failure.stdout.is_empty());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("global confirmation"));
    let uniform=too_short.replace("\"strategy\":\"goal-recovery\",\"marking_fraction\":0.5,","");
    let uniform=run(&uniform);let study=uniform.get("mesh_convergence").unwrap();
    assert_eq!(n(study,"meshes_solved"),3.0);
    assert_eq!(study.str_field("method"),Some("uniform-red-tet-refinement"));
    assert_eq!(n(study,"total_adjoint_sweeps"),0.0);
}
#[test]
fn contact_vertex_order_does_not_change_the_adaptive_physical_problem() {
    let original=run(FIXTURE);
    let reordered=FIXTURE.replace("\"side_b\":[12,13,15]","\"side_b\":[15,12,13]")
        .replace("\"side_b\":[12,14,15]","\"side_b\":[14,15,12]");
    assert_ne!(reordered,FIXTURE);
    let reordered=run(&reordered);
    assert_eq!(original.get("solid_temperatures_k"),reordered.get("solid_temperatures_k"));
    assert_eq!(original.get("contacts"),reordered.get("contacts"));
    assert_eq!(original.path(&["mesh_convergence","history"]),reordered.path(&["mesh_convergence","history"]));
}
#[test]
fn marker_adjoint_and_resource_failures_never_publish_an_adaptive_result() {
    for (from,to) in [
        ("\"marking_fraction\":0.5","\"marking_fraction\":0"),
        ("\"marking_fraction\":0.5","\"marking_fraction\":1.1"),
        ("\"strategy\":\"goal-recovery\"","\"strategy\":\"unknown\""),
        ("\"max_tetrahedra\":100000","\"max_tetrahedra\":12"),
        ("\"derivative_iterations\":400","\"derivative_iterations\":1"),
        ("\"wall_seconds\":120","\"wall_seconds\":1e-12"),
        ("\"gradient\":false","\"gradient\":true"),
    ] {
        assert!(FIXTURE.contains(from));let result=output(&FIXTURE.replace(from,to));
        assert!(!result.status.success(),"accepted {to}");assert!(result.stdout.is_empty(),"published partial result for {to}");
        assert!(!result.stderr.is_empty());
    }
}
