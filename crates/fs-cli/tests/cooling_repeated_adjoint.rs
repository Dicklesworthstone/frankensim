#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static FIXTURE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-contact-pulse.json"))));
const ADJOINT: &str = "\"adjoint\":{\"qoi\":\"sampled-peak\",\"max_checkpoint_bytes\":1048576},";
const REPEAT: &str = "\"repeat\":{\"cycles\":3,\"max_total_steps\":21},";
const PHASES: &str = "{\"duration_s\":6,\"power_scale\":1,\"fan_speed_ratio\":1},{\"duration_s\":8,\"power_scale\":0,\"fan_speed_ratio\":1.5}";

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
fn small() -> String {
    assert!(FIXTURE.contains(ADJOINT));
    FIXTURE.replace("\"duration_s\":30", "\"duration_s\":6")
        .replace("\"duration_s\":120", "\"duration_s\":8")
}
fn repeated() -> String { small().replace(ADJOINT, &format!("{REPEAT}{ADJOINT}")) }
fn plain() -> String { repeated().replace(ADJOINT, "") }
fn at<'a>(j: &'a J, path: &[&str]) -> &'a J { j.path(path).unwrap() }
fn number(j: &J, key: &str) -> f64 { j.f64_field(key).unwrap() }
fn gradient(j: &J) -> &J { at(j,&["repeated_cycles","adjoint"]) }
fn value(j: &J, qoi: &str) -> f64 {
    if qoi=="final" { number(at(j,&["objective"]),"value_k") }
    else { number(at(j,&["repeated_cycles"]),"sampled_peak_objective_k") }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol*b.abs().max(1.0), "{a:e} != {b:e}");
}
fn perturb(text: &str, axis: usize, eps: f64) -> String {
    let changed = match axis {
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
    };
    assert_ne!(changed,text);
    changed
}

#[test]
fn all_cycle_gradients_match_full_nonlinear_contact_trajectory_perturbations() {
    for qoi in ["final","sampled-peak"] {
        let source = repeated().replace("\"qoi\":\"sampled-peak\"",&format!("\"qoi\":\"{qoi}\""));
        let result = run(&source);
        let g = gradient(&result);
        let rows = g.get("intervals").unwrap().as_array().unwrap();
        let actual = [number(g,"dtemperature_duniform_initial_k"),
            number(&rows[0],"dtemperature_dpower_multiplier_k"),
            number(g,"dtemperature_dcapacity_multiplier_k"),
            number(&rows[0],"dtemperature_dlog_fan_speed_ratio_k"),
            number(&rows[1],"dtemperature_dlog_fan_speed_ratio_k")];
        for (axis,&derivative) in actual.iter().enumerate() {
            let eps = 2e-3;
            let plus = run(&perturb(&plain(),axis,eps));
            let minus = run(&perturb(&plain(),axis,-eps));
            close(derivative,(value(&plus,qoi)-value(&minus,qoi))/(2.0*eps),4e-5);
        }
        close(number(g,"value_k"),value(&result,qoi),1e-12);
        assert_eq!(number(g,"cycles"),3.0);
        assert_eq!(number(g,"reconstructed_solid_endpoints"),number(g,"state_index"));
        if qoi=="sampled-peak" {
            assert!(number(g,"time_s")>28.0, "fixture peak must lie after two warm-up cycles");
            // The LAST cooldown is in the future, but the same base control
            // also changes earlier cooldowns. Dropping their carry gives zero.
            assert!(actual[4].abs()>1e-6);
        }
    }
}

#[test]
fn repeated_forward_history_is_unchanged_and_reverse_work_is_not_hidden() {
    let expected = run(&plain());
    let actual = run(&repeated());
    for path in [&["solid_temperatures_k"][..], &["transient","history"][..],
        &["repeated_cycles","cycles"][..]] {
        assert_eq!(at(&actual,path),at(&expected,path));
    }
    let r = at(&actual,&["repeated_cycles"]);
    assert_eq!(number(r,"total_solid_solves"),number(r,"forward_solid_solves")
        +number(gradient(&actual),"reconstructed_solid_endpoints"));
    assert_eq!(at(&actual,&["transient","adjoint"]),&J::Null,
        "a complete-history adjoint must not be mislabeled as a last-cycle derivative");
    assert_eq!(at(&actual,&["fan_speed_sensitivity"]),&J::Null);
    assert_eq!(gradient(&actual),gradient(&run(&repeated())));
}

#[test]
fn shared_phase_gradients_equal_the_sum_of_explicitly_unrolled_occurrences() {
    assert!(repeated().contains(PHASES));
    let unrolled = repeated().replace(REPEAT,"").replace(PHASES,&[PHASES;3].join(","));
    for spatial in ["\"max_solid_temperature\":true", "\"mean_wall_region\":\"last-face\""] {
        let actual = run(&repeated().replace("\"max_solid_temperature\":true",spatial));
        let expected = run(&unrolled.replace("\"max_solid_temperature\":true",spatial));
        assert_eq!(at(&actual,&["solid_temperatures_k"]),at(&expected,&["solid_temperatures_k"]));
        let ga = gradient(&actual);
        let ge = at(&expected,&["transient","adjoint"]);
        for key in ["value_k","time_s","state_index","dtemperature_duniform_initial_k",
            "dtemperature_dcapacity_multiplier_k"] { close(number(ga,key),number(ge,key),1e-10); }
        let a = ga.get("intervals").unwrap().as_array().unwrap();
        let e = ge.get("intervals").unwrap().as_array().unwrap();
        for phase in 0..2 {
            for key in ["dtemperature_dpower_multiplier_k","dtemperature_dlog_fan_speed_ratio_k"] {
                let sum: f64 = (0..3).map(|cycle| number(&e[2*cycle+phase],key)).sum();
                close(number(&a[phase],key),sum,1e-10);
            }
        }
    }
}

#[test]
fn one_repeated_cycle_is_the_same_discrete_adjoint_as_the_original_single_cycle() {
    let source = repeated().replace("\"cycles\":3","\"cycles\":1");
    let repeated = run(&source);
    let single = run(&small());
    assert_eq!(gradient(&repeated),at(&single,&["transient","adjoint"]));
}

#[test]
fn initial_all_cycle_peak_keeps_zero_future_influence_and_no_reverse_solves() {
    let source = repeated().replace("\"initial_temperature_k\":300","\"initial_temperature_k\":350")
        .replace("\"power_scale\":1,","\"power_scale\":0,")
        .replace("\"duration_s\":6,","\"duration_s\":600,")
        .replace("\"duration_s\":8,","\"duration_s\":600,")
        .replace("\"max_step_s\":2,","\"max_step_s\":600,");
    let result = run(&source);
    let g = gradient(&result);
    for key in ["time_s","state_index","reconstructed_solid_endpoints","adjoint_sweeps",
        "dtemperature_dcapacity_multiplier_k"] { assert_eq!(number(g,key),0.0); }
    assert_eq!(number(g,"dtemperature_duniform_initial_k"),1.0);
    assert!(number(at(&result,&["transient"]),"sampled_peak_objective_k")<value(&result,"sampled-peak"));
}

#[test]
fn checkpoint_budget_covers_all_cycles_and_state_dependent_horizons_refuse() {
    let single = run(&small());
    let one_cycle_bytes = number(at(&single,&["transient","adjoint"]),"checkpoint_bytes") as usize;
    let periodic = "\"repeat\":{\"until_periodic\":{\"max_cycles\":3,\"temperature_tolerance_k\":1,\"consecutive_cycles\":2},\"max_total_steps\":21},";
    let controlled = "\"repeat\":{\"cycles\":3,\"max_total_steps\":21,\"fan_controller\":{}},";
    for source in [
        repeated().replace("\"max_checkpoint_bytes\":1048576",&format!("\"max_checkpoint_bytes\":{one_cycle_bytes}")),
        repeated().replace(REPEAT,periodic),
        repeated().replace(REPEAT,controlled),
        repeated().replace("\"max_total_steps\":21","\"max_total_steps\":20"),
        repeated().replace("\"wall_seconds\":120","\"wall_seconds\":1e-12"),
    ] {
        let result = output(&source);
        assert!(!result.status.success());
        assert!(result.stdout.is_empty(),"no partial trajectory or derivative may escape");
        assert!(!result.stderr.is_empty());
    }
}
