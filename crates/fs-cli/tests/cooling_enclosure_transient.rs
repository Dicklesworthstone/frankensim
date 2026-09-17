//! Actual-command enclosure/storage regressions. Reference fields use direct
//! P1 equations with analytically eliminated two-plate radiation and mixed air.
#![cfg(unix)]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enclosure-radiation-pulse.json"));
const STEADY: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enclosure-radiation-gap.json"));
fn num(value: f64) -> J { J::Number { value, raw: value.to_string() } }
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object required") };
    &mut fields.iter_mut().find(|(k, _)| k == key).unwrap().1
}
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object required") };
    if let Some((_, v)) = fields.iter_mut().find(|(k, _)| k == key) { *v = value; }
    else { fields.push((key.into(), value)); }
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object required") };
    fields.retain(|(k, _)| k != key);
}
fn rows(root: &mut J) -> &mut Vec<J> { let J::Array(rows) = root else { panic!("array") }; rows }
fn text(root: &J) -> String {
    match root {
        J::Null => "null".into(), J::Bool(v) => v.to_string(), J::Number { raw, .. } => raw.clone(),
        J::Str(s) => format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")),
        J::Array(a) => format!("[{}]", a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(a) => format!("{{{}}}", a.iter().map(|(k, v)| format!("{}:{}", text(&J::Str(k.clone())), text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(input: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(input).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn success(out: &Output) -> J {
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    J::parse(std::str::from_utf8(&out.stdout).unwrap()).unwrap()
}
fn run(input: &J) -> J { success(&output(input)) }
fn n(root: &J, key: &str) -> f64 { root.f64_field(key).unwrap() }
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() < tolerance, "{a} versus {b}"); }
fn phase(result: &J) -> &J { result.get("repeated_cycles").unwrap_or_else(|| result.get("transient").unwrap()) }
fn peak(result: &J) -> f64 { n(phase(result), "sampled_peak_objective_k") }
fn fields(result: &J) -> Vec<f64> {
    result.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect()
}
fn compare_fields(a: &J, b: &J, tolerance: f64) {
    for (a, b) in fields(a).iter().zip(fields(b)) { near(*a, b, tolerance); }
}
fn energy(result: &J) {
    let p = phase(result);
    near(n(p, "stored_energy_change_j") + n(p, "air_energy_gain_j")
        + n(p, "radiative_energy_loss_j"), n(p, "input_energy_j"), 2e-5);
    near(n(p, "radiative_energy_loss_j"), 0.0, 2e-5);
    near(n(result.get("radiation").unwrap(), "radiative_out_w"), 0.0, 1e-7);
}
fn single() -> J {
    let mut input = J::parse(BASE).unwrap();
    remove(member(&mut input, "transient"), "repeat"); input
}

#[test]
fn manufactured_new_endpoints_use_implicit_radiation_in_both_heat_directions() {
    for (a, b, old) in [
        (325.0, 303.0, [320.47717363684444,331.0229474736791,323.3863526263162,326.6591789894716,
            323.3863526263156,326.6591789894722,324.47729474736707,325.9318842421039,
            302.73286427804163,302.91095475934753,302.6438190373891,303.0445226203262,
            302.64381903738905,303.0445226203263,302.10954759347203,303.4007035829375]),
        (303.0, 325.0, [302.2610922638762,300.6142394136088,301.8067880293193,301.2956957654436,
            301.8067880293192,301.29569576544407,301.6364239413605,301.40927182408274,
            334.177598580991,328.0591995269971,337.2367981079884,323.47040023650146,
            337.2367981079883,323.47040023650146,355.59199526997054,311.2336021285134]),
    ] {
        let mut input = J::parse(STEADY).unwrap();
        let mut schedule = J::parse(r#"{"max_step_s":1,"max_steps":1,"element_heat_capacities_j_m3_k":[20000,20000,20000,20000,20000,20000,10000,10000,10000,10000,10000,10000],"intervals":[{"duration_s":1,"power_scale":1,"fan_speed_ratio":1}]}"#).unwrap();
        put(&mut schedule, "initial_temperatures_k", J::Array(old.into_iter().map(num).collect()));
        put(&mut input, "transient", schedule);
        let result = run(&input); energy(&result);
        for (index, t) in fields(&result).iter().enumerate() { near(*t, if index < 8 { a } else { b }, 2e-6); }
        let expected = 0.01*5.670374419e-8*(a-b)*(a+b)*(a*a+b*b)/(1.0/0.8+1.0/0.6-1.0);
        for patch in result.path(&["radiation", "surfaces"]).unwrap().as_array().unwrap() {
            near(n(patch, "outward_heat_w"), if patch.str_field("surface") == Some("emitter") { expected } else { -expected }, 1e-7);
        }
        assert_eq!(result.path(&["radiation", "temporal_scope"]).unwrap().as_str(), Some("final-accepted-endpoint"));
        assert_eq!(result.path(&["radiation", "iterations"]), Some(&J::Null));
        assert!(n(phase(&result), "total_solid_solves") > 1.0);
    }
}

#[test]
fn nonlinear_repeated_pulse_matches_independent_fields_and_keeps_internal_heat_internal() {
    let result = run(&J::parse(BASE).unwrap()); energy(&result);
    let p = phase(&result);
    near(peak(&result), 315.42790320839714, 3e-5);
    near(n(p, "sampled_peak_time_s"), 40.0, 1e-12);
    near(n(p, "elapsed_time_s"), 60.0, 1e-12);
    assert_eq!(n(p, "total_accepted_steps"), 30.0);
    near(n(p, "input_energy_j"), 100.0, 1e-7);
    near(n(p, "stored_energy_change_j"), 35.006203990549764, 3e-5);
    near(n(p, "air_energy_gain_j"), 64.99379600945154, 3e-5);
    near(fields(&result)[0], 308.54369447940695, 3e-5);
    assert!(fields(&result)[8] > 301.0);
    assert_eq!(result.get("adjoint_residual"), Some(&J::Null));
}

#[test]
fn repeated_history_and_matrix_axis_permutations_replay_without_changing_physics() {
    let input = J::parse(BASE).unwrap(); let original = output(&input); let repeated = success(&original);
    let mut unrolled = input.clone();
    let schedule = member(&mut unrolled, "transient"); remove(schedule, "repeat");
    let intervals = rows(member(schedule, "intervals")); intervals.extend(intervals.clone());
    let unrolled = run(&unrolled); energy(&unrolled);
    assert_eq!(repeated.get("solid_temperatures_k"), unrolled.get("solid_temperatures_k"));
    assert_eq!(repeated.get("radiation"), unrolled.get("radiation"));
    near(peak(&repeated), peak(&unrolled), 1e-10);
    for key in ["stored_energy_change_j", "input_energy_j", "air_energy_gain_j"] {
        near(n(phase(&repeated), key), n(phase(&unrolled), key), 1e-7);
    }
    let mut reordered = input;
    let enclosure = member(member(&mut reordered, "radiation"), "enclosure");
    rows(member(enclosure, "surfaces")).reverse();
    let matrix = rows(member(enclosure, "view_factors")); matrix.reverse();
    for row in matrix { rows(row).reverse(); }
    let replay = output(&reordered); success(&replay); assert_eq!(original.stdout, replay.stdout);
}

#[test]
fn adaptive_trial_work_does_not_enter_accepted_storage_or_heat_history() {
    let mut input = single();
    put(member(&mut input, "transient"), "adaptive", J::parse(r#"{"absolute_tolerance_k":0.005,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":1000}"#).unwrap());
    put(member(&mut input, "transient"), "max_steps", num(1000.0));
    let adaptive = run(&input); energy(&adaptive);
    let history = adaptive.path(&["transient", "history"]).unwrap().as_array().unwrap();
    let intervals = history.iter().skip(1).map(|row| J::Object(vec![
        ("duration_s".into(), row.get("dt_s").unwrap().clone()),
        ("power_scale".into(), row.get("power_scale").unwrap().clone()),
        ("fan_speed_ratio".into(), row.get("fan_speed_ratio").unwrap().clone()),
        ("steps".into(), num(1.0)),
    ])).collect();
    remove(member(&mut input, "transient"), "adaptive");
    put(member(&mut input, "transient"), "intervals", J::Array(intervals));
    let replay = run(&input); energy(&replay); compare_fields(&adaptive, &replay, 3e-6);
    for key in ["stored_energy_change_j", "input_energy_j", "air_energy_gain_j"] {
        near(n(phase(&adaptive), key), n(phase(&replay), key), 1e-5);
    }
    assert!(n(phase(&adaptive), "total_solid_solves") > n(phase(&replay), "total_solid_solves"));
}

#[test]
fn time_grid_study_retains_enclosure_exchange_and_replays_the_final_grid() {
    let mut input = J::parse(BASE).unwrap();
    put(member(&mut input, "transient"), "time_convergence", J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":0.2,"max_total_steps":1000,"max_trace_bytes":1048576}"#).unwrap());
    let result = run(&input); energy(&result);
    let study = result.get("time_convergence").unwrap();
    assert_eq!(n(study, "trajectories_solved"), 3.0);
    assert_eq!(n(study, "total_accepted_steps"), 210.0);
    let counts = study.get("final_steps_per_interval").unwrap().as_array().unwrap();
    let schedule = member(&mut input, "transient"); remove(schedule, "time_convergence");
    for (row, count) in rows(member(schedule, "intervals")).iter_mut().zip(counts) { put(row, "steps", count.clone()); }
    let replay = run(&input);
    assert_eq!(result.get("solid_temperatures_k"), replay.get("solid_temperatures_k"));
    assert_eq!(result.get("radiation"), replay.get("radiation"));
    assert_eq!(result.get("repeated_cycles"), replay.get("repeated_cycles"));
}

#[test]
fn workload_sizing_rechecks_complete_radiating_candidates() {
    let mut input = single(); let reference = run(&input);
    let schedule = member(&mut input, "transient");
    put(schedule, "temperature_limit_k", num(peak(&reference)));
    put(schedule, "power_design", J::parse(r#"{"min_power_multiplier":0.8,"max_power_multiplier":1.2,"power_multiplier_tolerance":0.5,"temperature_tolerance_k":10,"max_evaluations":3}"#).unwrap());
    let designed = run(&input); energy(&designed);
    let decision = designed.get("transient_power_design").unwrap();
    assert_eq!(decision.str_field("status"), Some("target-bracketed"));
    let scale = n(decision, "selected_power_multiplier"); near(scale, 0.8, 1e-15);
    let schedule = member(&mut input, "transient"); remove(schedule, "power_design");
    for row in rows(member(schedule, "intervals")) {
        let value = n(row, "power_scale"); put(row, "power_scale", num(value*scale));
    }
    let replay = run(&input); compare_fields(&designed, &replay, 1e-8);
    assert_eq!(designed.get("radiation"), replay.get("radiation"));
    assert!(peak(&designed) <= peak(&reference));
}

#[test]
fn enclosure_adjoint_and_exhausted_radiation_step_or_time_budgets_never_publish() {
    let input = single();
    let mut gradient = input.clone();
    put(member(&mut gradient, "transient"), "adjoint", J::parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#).unwrap());
    let rejected = output(&gradient);
    assert!(!rejected.status.success()); assert!(rejected.stdout.is_empty());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("adjoint"));
    for (section, key, value) in [("radiation", "max_iterations", 1.0),
        ("transient", "max_steps", 1.0), ("budgets", "wall_seconds", 1e-9)] {
        let mut failed = input.clone(); put(member(&mut failed, section), key, num(value));
        let rejected = output(&failed); assert_eq!(rejected.status.code(), Some(6));
        assert!(rejected.stdout.is_empty());
    }
}

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("fs-enclosure-transient-{}-{}",
            std::process::id(), NEXT.fetch_add(1, Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap(); Self(path)
    }
    fn uq(&self, plan: &str, extra: &[&str]) -> Output {
        std::fs::write(self.0.join("uq.json"), plan).unwrap();
        Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json", "cooling-network-uq", "base.json", "uq.json"]).args(extra).output().unwrap()
    }
}
impl Drop for Scratch { fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); } }

#[test]
fn uncertain_finish_uses_whole_trajectories_and_resumes_exact_sample_history() {
    let dir = Scratch::new(); let input = J::parse(BASE).unwrap();
    std::fs::write(dir.0.join("base.json"), text(&input)).unwrap();
    let plan = r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":2,"wall_seconds":600,"qoi":{"kind":"transient-sampled-peak"},"temperature_limit_k":316,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"radiation-emissivity","surface":"emitter"},"distribution":{"kind":"uniform","lo":0.2,"hi":0.8}}]}"#;
    let zero = plan.replace("\"lo\":0.2", "\"lo\":0.8");
    let fixed = success(&dir.uq(&zero, &[])); near(n(&fixed, "mean_k"), peak(&run(&input)), 1e-9);
    assert_eq!(n(&fixed, "std_dev_k"), 0.0);
    let full = dir.uq(plan, &["--checkpoint", "full.uqcp"]); let parsed = success(&full);
    assert!(n(&parsed, "std_dev_k") > 1e-5);
    let part = dir.uq(plan, &["--checkpoint", "part.uqcp", "--max-new-samples", "1"]);
    assert_eq!(part.status.code(), Some(6));
    let original = std::fs::read(dir.0.join("part.uqcp")).unwrap();
    let resumed = dir.uq(plan, &["--resume", "part.uqcp", "--checkpoint", "done.uqcp"]); success(&resumed);
    assert_eq!(full.stdout, resumed.stdout);
    assert_eq!(std::fs::read(dir.0.join("full.uqcp")).unwrap(), std::fs::read(dir.0.join("done.uqcp")).unwrap());
    assert_eq!(original, std::fs::read(dir.0.join("part.uqcp")).unwrap());
}
