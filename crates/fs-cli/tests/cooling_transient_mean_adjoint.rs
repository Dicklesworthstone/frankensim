#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Stdio};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static FIXTURE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-contact-pulse.json"))));
fn run(source: &str) -> J {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(source.as_bytes()).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(&String::from_utf8(output.stdout).unwrap()).unwrap()
}
fn mean(j: &J) -> f64 { j.path(&["objective","value_k"]).unwrap().as_f64().unwrap() }

#[test]
fn boundary_declaration_order_does_not_select_the_wrong_mean_wall_adjoint() {
    let source = FIXTURE.replace("\"duration_s\":30","\"duration_s\":6")
        .replace("\"duration_s\":120","\"duration_s\":8")
        .replace("\"qoi\":\"sampled-peak\"","\"qoi\":\"final\"")
        .replace("\"max_solid_temperature\":true","\"mean_wall_region\":\"first-face\"");
    let first = "{\"name\":\"first-face\",\"faces\":[[0,3,9],[0,6,9]],\"htc_w_m2_k\":50}";
    let last = "{\"name\":\"last-face\",\"faces\":[[2,5,11],[2,8,11]],\"htc_w_m2_k\":80}";
    assert!(source.contains(first) && source.contains(last));
    let reordered = source.replace(first,"SURFACE_SWAP").replace(last,first).replace("SURFACE_SWAP",last);
    let a=run(&source);let b=run(&reordered);
    assert_eq!(a.path(&["solid_temperatures_k"]),b.path(&["solid_temperatures_k"]));
    assert_eq!(a.path(&["transient","adjoint"]),b.path(&["transient","adjoint"]));
    let interval=&b.path(&["transient","adjoint","intervals"]).unwrap().as_array().unwrap()[0];
    let actual=interval.f64_field("dtemperature_dpower_multiplier_k").unwrap();
    let eps=1e-3;
    let plus=reordered.replace("\"power_scale\":1,",&format!("\"power_scale\":{},",1.0+eps));
    let minus=reordered.replace("\"power_scale\":1,",&format!("\"power_scale\":{},",1.0-eps));
    let expected=(mean(&run(&plus))-mean(&run(&minus)))/(2.0*eps);
    assert!((actual-expected).abs()<3e-4*expected.abs().max(1.0),"{actual} != {expected}");
}
