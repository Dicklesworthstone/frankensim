#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-repeated-contact-pulse.json"));
const ADJOINT: &str = "\"adjoint\":{\"qoi\":\"sampled-peak\",\"max_checkpoint_bytes\":1048576},";
const REPEAT: &str = "\"repeat\":{\"cycles\":3,\"max_total_steps\":21},";

fn output(text: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
        .spawn().unwrap();
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(text: &str) -> J {
    let result = output(text);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    J::parse(&String::from_utf8(result.stdout).unwrap()).unwrap()
}
fn base(repeated: bool) -> String {
    assert!(FIXTURE.contains(ADJOINT) && FIXTURE.contains(REPEAT));
    if repeated { FIXTURE.to_owned() } else { FIXTURE.replace(REPEAT, "") }
}
fn at<'a>(value: &'a J, path: &[&str]) -> &'a J { value.path(path).unwrap() }
fn n(value: &J, name: &str) -> f64 { value.f64_field(name).unwrap() }
fn trajectory(value: &J, repeated: bool) -> &J {
    value.get(if repeated { "repeated_cycles" } else { "transient" }).unwrap()
}
fn peak(value: &J, repeated: bool) -> f64 { n(trajectory(value,repeated), "sampled_peak_objective_k") }
fn name(power: bool) -> &'static str {
    if power { "transient_power_design" } else { "transient_fan_speed_design" }
}
fn multiplier_key(power: bool) -> &'static str {
    if power { "power_multiplier" } else { "speed_multiplier" }
}
fn scaled(text: &str, power: bool, multiplier: f64) -> String {
    if power {
        if text.contains("\"power_scale\":1,") {
            text.replace("\"power_scale\":1,", &format!("\"power_scale\":{multiplier},"))
        } else {
            assert!(text.contains("\"component_powers_w\":{\"chip\":20}"));
            text.replace("\"component_powers_w\":{\"chip\":20}",
                &format!("\"component_powers_w\":{{\"chip\":{}}}",20.0*multiplier))
        }
    } else {
        assert!(text.contains("\"fan_speed_ratio\":1}") && text.contains("\"fan_speed_ratio\":1.5}"));
        text.replace("\"fan_speed_ratio\":1}", &format!("\"fan_speed_ratio\":{multiplier}}}"))
            .replace("\"fan_speed_ratio\":1.5}", &format!("\"fan_speed_ratio\":{}}}",1.5*multiplier))
    }
}
fn design(text: &str, power: bool, limit: f64, gradient: bool) -> String {
    let controls = if power {
        "\"power_design\":{\"min_power_multiplier\":0,\"max_power_multiplier\":2,\"power_multiplier_tolerance\":0.0001,\"temperature_tolerance_k\":0.00001,\"max_evaluations\":64},"
    } else {
        "\"fan_speed_design\":{\"min_speed_multiplier\":0.5,\"max_speed_multiplier\":1.3,\"speed_multiplier_tolerance\":0.0001,\"temperature_tolerance_k\":0.00001,\"max_evaluations\":64},"
    };
    assert!(text.contains("\"temperature_limit_k\":305") && text.contains(ADJOINT));
    text.replace("\"temperature_limit_k\":305", &format!("\"temperature_limit_k\":{limit}"))
        .replace(ADJOINT, &format!("{controls}{}",if gradient {ADJOINT}else{""}))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn power_and_fan_designs_replay_the_actual_candidate_and_its_peak_derivative() {
    for repeated in [false, true] {
        for power in [false, true] {
            let source = base(repeated);
            let plain = source.replace(ADJOINT, "");
            // Define a reachable target using a separate complete forward run,
            // not a steady surrogate or an assumed monotone correlation.
            let target_multiplier = if power {0.8}else{0.85};
            let limit = peak(&run(&scaled(&plain,power,target_multiplier)),repeated);
            let result = run(&design(&source,power,limit,true));
            let search = result.get(name(power)).unwrap();
            let selected = n(search, &format!("selected_{}",multiplier_key(power)));
            assert!((selected-1.0).abs()>0.05, "fixture must test the non-unit chain rule");
            assert!(n(search,"newton_trials")>0.0);
            assert!(peak(&result,repeated)<=limit);
            assert!(limit-peak(&result,repeated)<=1e-5);
            assert!(n(search,"multiplier_bracket_width")<=1e-4);
            let failed = search.get(if power {"failed_upper"}else{"failed_lower"}).unwrap();
            assert!(n(failed,"sampled_peak_objective_k")>limit);
            let fresh = run(&scaled(&plain,power,selected));
            assert_eq!(result.get("solid_temperatures_k"),fresh.get("solid_temperatures_k"));
            assert_eq!(at(&result,&["transient","history"]),at(&fresh,&["transient","history"]));
            if repeated {
                assert_eq!(at(&result,&["repeated_cycles","cycles"]),at(&fresh,&["repeated_cycles","cycles"]));
                assert_eq!(at(&result,&["transient","adjoint"]),&J::Null);
            }
            let grad = trajectory(&result,repeated).get("adjoint").unwrap();
            close(n(grad,"value_k"),peak(&result,repeated),0.0);
            let relative: f64 = grad.get("intervals").unwrap().as_array().unwrap().iter()
                .map(|row| n(row,if power {"dtemperature_dpower_multiplier_k"}else{"dtemperature_dlog_fan_speed_ratio_k"})).sum();
            let slope = relative/selected;
            let trials = search.get("history").unwrap().as_array().unwrap();
            let accepted = trials.iter().find(|t| n(t,multiplier_key(power))==selected).unwrap();
            close(n(accepted,"dpeak_dmultiplier_k"),slope,1e-12);
            let eps = 1e-3;
            let plus = peak(&run(&scaled(&plain,power,selected+eps)),repeated);
            let minus = peak(&run(&scaled(&plain,power,selected-eps)),repeated);
            close(slope,(plus-minus)/(2.0*eps),3e-4);
            close(n(search,"total_solid_solves"),trials.iter().map(|t|n(t,"solid_solves")).sum(),0.0);
            let retained = trajectory(&result,repeated);
            assert!(n(retained,"total_solid_solves")>n(retained,"forward_solid_solves"));
        }
    }
}

#[test]
fn omitting_adjoint_keeps_bisection_and_does_no_hidden_reverse_work() {
    for power in [false,true] {
        let source = base(true);
        let plain = source.replace(ADJOINT,"");
        let limit = peak(&run(&scaled(&plain,power,0.85)),true);
        let result = run(&design(&source,power,limit,false));
        let search = result.get(name(power)).unwrap();
        assert_eq!(search.get("search_method").unwrap().as_str(),Some("bisection"));
        assert_eq!(n(search,"newton_trials"),0.0);
        assert!(search.get("history").unwrap().as_array().unwrap().iter()
            .all(|t|t.get("dpeak_dmultiplier_k")==Some(&J::Null)));
        let run = trajectory(&result,true);
        assert_eq!(run.get("adjoint"),Some(&J::Null));
        assert_eq!(n(run,"total_solid_solves"),n(run,"forward_solid_solves"));
    }
}

#[test]
fn named_workloads_and_zero_endpoint_keep_their_actual_source_semantics() {
    let source = base(true).replace("\"power_scale\":1,","\"component_powers_w\":{\"chip\":20},")
        .replace("\"power_scale\":0,","\"component_powers_w\":{\"chip\":0},");
    let plain = source.replace(ADJOINT,"");
    let limit = peak(&run(&scaled(&plain,true,0.8)),true);
    let result = run(&design(&source,true,limit,true));
    let search = result.get(name(true)).unwrap();
    let selected = n(search,"selected_power_multiplier");
    let trials = search.get("history").unwrap().as_array().unwrap();
    let zero = trials.iter().find(|row|n(row,"power_multiplier")==0.0).unwrap();
    assert_eq!(zero.get("dpeak_dmultiplier_k"),Some(&J::Null));
    let fresh = run(&scaled(&plain,true,selected));
    assert_eq!(result.get("solid_temperatures_k"),fresh.get("solid_temperatures_k"));
    for state in &at(&result,&["transient","history"]).as_array().unwrap()[1..] {
        let watts = at(state,&["component_powers_w","chip"]).as_f64().unwrap();
        close(watts,if n(state,"interval")==0.0 {20.0*selected}else{0.0},1e-12);
    }
}

#[test]
fn wrong_observable_unsupported_controls_and_exhaustion_publish_no_design() {
    let request = design(&base(true),true,303.5,true);
    for (from,to) in [
        ("\"qoi\":\"sampled-peak\"","\"qoi\":\"final\""),
        ("\"max_checkpoint_bytes\":1048576","\"max_checkpoint_bytes\":1"),
        ("\"derivative_iterations\":400","\"derivative_iterations\":1"),
        ("\"wall_seconds\":120","\"wall_seconds\":1e-12"),
        ("\"max_evaluations\":64","\"max_evaluations\":2"),
        ("\"adjoint\":{","\"adaptive\":{},\"adjoint\":{"),
        ("\"cycles\":3","\"until_periodic\":{}"),
        ("\"cycles\":3","\"cycles\":3,\"fan_controller\":{}"),
    ] {
        assert!(request.contains(from));
        let result=output(&request.replace(from,to));
        assert!(!result.status.success(),"accepted {to}");
        assert!(result.stdout.is_empty(),"partial design escaped on {to}");
        assert!(!result.stderr.is_empty());
    }
}
