#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-contact-pulse.json"));
const ADJOINT: &str = "\"adjoint\":{\"qoi\":\"sampled-peak\",\"max_checkpoint_bytes\":1048576},";

fn output(text: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(text: &str) -> J {
    let out = output(text);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    J::parse(&String::from_utf8(out.stdout).unwrap()).unwrap()
}
fn small() -> String { FIXTURE.replace("\"duration_s\":30", "\"duration_s\":6")
    .replace("\"duration_s\":120", "\"duration_s\":8") }
fn plain() -> String { assert!(small().contains(ADJOINT)); small().replace(ADJOINT, "") }
fn at<'a>(j: &'a J, path: &[&str]) -> &'a J { j.path(path).unwrap() }
fn value(j: &J, qoi: &str) -> f64 {
    if qoi == "final" { at(j,&["objective","value_k"]).as_f64().unwrap() }
    else { at(j,&["transient","sampled_peak_objective_k"]).as_f64().unwrap() }
}
fn close(a:f64,b:f64) { assert!((a-b).abs() < 3e-4*b.abs().max(1.0), "{a:e} != {b:e}"); }

fn perturb(text: &str, axis: usize, eps: f64) -> String {
    match axis {
        0 => text.replace("\"initial_temperature_k\":300", &format!("\"initial_temperature_k\":{}",300.0+eps)),
        1 => text.replace("\"power_scale\":1,", &format!("\"power_scale\":{},",1.0+eps)),
        2 => {
            let original = "\"element_heat_capacities_j_m3_k\":[2000000,2000000,2000000,2000000,2000000,2000000,1000000,1000000,1000000,1000000,1000000,1000000]";
            assert!(text.contains(original));
            let capacities = (0..12).map(|i| ((if i<6 {2e6} else {1e6})*(1.0+eps)).to_string()).collect::<Vec<_>>().join(",");
            text.replace(original,&format!("\"element_heat_capacities_j_m3_k\":[{capacities}]"))
        }
        3 => text.replace("\"fan_speed_ratio\":1}", &format!("\"fan_speed_ratio\":{}}}",eps.exp())),
        4 => text.replace("\"fan_speed_ratio\":1.5}", &format!("\"fan_speed_ratio\":{}}}",1.5*eps.exp())),
        _ => panic!("unknown control"),
    }
}

#[test]
fn final_and_sampled_peak_gradients_match_complete_perturbed_trajectories() {
    for qoi in ["final","sampled-peak"] {
        let source = small().replace("\"qoi\":\"sampled-peak\"",&format!("\"qoi\":\"{qoi}\""));
        let result = run(&source);
        let gradient = at(&result,&["transient","adjoint"]);
        let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
        let actual = [
            gradient.f64_field("dtemperature_duniform_initial_k").unwrap(),
            intervals[0].f64_field("dtemperature_dpower_multiplier_k").unwrap(),
            gradient.f64_field("dtemperature_dcapacity_multiplier_k").unwrap(),
            intervals[0].f64_field("dtemperature_dlog_fan_speed_ratio_k").unwrap(),
            intervals[1].f64_field("dtemperature_dlog_fan_speed_ratio_k").unwrap(),
        ];
        let eps = 1e-3;
        for (axis, &derivative) in actual.iter().enumerate() {
            let plus = perturb(&plain(),axis,eps);
            let minus = perturb(&plain(),axis,-eps);
            assert_ne!(plus,plain()); assert_ne!(minus,plain());
            close(derivative,(value(&run(&plus),qoi)-value(&run(&minus),qoi))/(2.0*eps));
        }
        close(gradient.f64_field("value_k").unwrap(),value(&result,qoi));
        assert!(actual[1]>0.0);
    }
}

#[test]
fn requesting_adjoint_preserves_forward_bits_and_counts_reverse_reconstruction() {
    let no_gradient = run(&plain());
    let with_gradient = run(&small());
    for path in [&["solid_temperatures_k"][..], &["transient","history"][..]] {
        assert_eq!(at(&no_gradient,path),at(&with_gradient,path));
    }
    let transient = at(&with_gradient,&["transient"]);
    let gradient = transient.get("adjoint").unwrap();
    let rebuilt = gradient.f64_field("reconstructed_solid_endpoints").unwrap();
    assert_eq!(rebuilt,gradient.f64_field("state_index").unwrap());
    assert_eq!(transient.f64_field("total_solid_solves").unwrap(),
        transient.f64_field("forward_solid_solves").unwrap()+rebuilt);
    assert!(rebuilt>0.0 && rebuilt<=transient.f64_field("steps").unwrap());
    assert_eq!(at(&with_gradient,&["fan_speed_sensitivity"]),&J::Null,
        "steady gradient fields must stay separate from the trajectory adjoint");
    let again = run(&small());
    assert_eq!(at(&with_gradient,&["transient","adjoint"]),at(&again,&["transient","adjoint"]));
}

#[test]
fn future_cooldown_controls_cannot_change_an_earlier_peak() {
    let result = run(&small());
    let gradient = at(&result,&["transient","adjoint"]);
    assert!(gradient.f64_field("time_s").unwrap()<=6.0);
    let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
    assert_eq!(intervals[1].f64_field("dtemperature_dlog_fan_speed_ratio_k"),Some(0.0));
    assert_eq!(intervals[1].f64_field("dtemperature_dpower_multiplier_k"),Some(0.0));
}

#[test]
fn an_initial_state_maximum_needs_no_reverse_physics() {
    let source = FIXTURE.replace("\"initial_temperature_k\":300","\"initial_temperature_k\":350")
        .replace("\"power_scale\":1,","\"power_scale\":0,")
        .replace("\"duration_s\":30","\"duration_s\":600")
        .replace("\"duration_s\":120","\"duration_s\":600")
        .replace("\"max_step_s\":2,","\"max_step_s\":600,");
    let result = run(&source);
    let gradient = at(&result,&["transient","adjoint"]);
    assert_eq!(gradient.f64_field("state_index"),Some(0.0));
    assert_eq!(gradient.f64_field("reconstructed_solid_endpoints"),Some(0.0));
    assert_eq!(gradient.f64_field("adjoint_sweeps"),Some(0.0));
    assert_eq!(gradient.f64_field("dtemperature_duniform_initial_k"),Some(1.0));
    assert_eq!(gradient.f64_field("dtemperature_dcapacity_multiplier_k"),Some(0.0));
}

#[test]
fn unsupported_or_exhausted_adjoint_requests_never_publish_partial_gradients() {
    for source in [
        small().replace("\"max_checkpoint_bytes\":1048576","\"max_checkpoint_bytes\":1"),
        small().replace("\"qoi\":\"sampled-peak\"","\"qoi\":\"final\"")
            .replace("\"derivative_iterations\":400","\"derivative_iterations\":1"),
        small().replace("\"wall_seconds\":120","\"wall_seconds\":1e-12"),
        small().replace("\"adjoint\":{","\"repeat\":{},\"adjoint\":{"),
        small().replace("\"adjoint\":{","\"adaptive\":{},\"adjoint\":{"),
        small().replace("\"adjoint\":{","\"power_design\":{},\"adjoint\":{"),
    ] {
        let result = output(&source);
        assert!(!result.status.success());
        assert!(result.stdout.is_empty());
        assert!(!result.stderr.is_empty());
    }
}
