//! Real-binary trajectory derivatives across a latent plateau and a later
//! sensible endpoint. No temperature-to-enthalpy inverse supplies the history.
#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;

#[path = "cooling_enthalpy_adjoint/fan.rs"]
mod fan;
#[path = "cooling_enthalpy_adjoint/repeat.rs"]
mod repeat;

use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enthalpy-phase-pulse.json"
));
const INITIAL: [f64; 4] = [2000.0, 2030.0, 1970.0, 2010.0];

fn number(value: f64) -> J {
    J::Number {
        value,
        raw: value.to_string(),
    }
}
fn parse(value: &str) -> J {
    J::parse(value).unwrap()
}
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else {
        panic!("expected object")
    };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}
fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else {
        panic!("expected object")
    };
    if let Some((_, slot)) = fields.iter_mut().find(|(name, _)| name == key) {
        *slot = value;
    } else {
        fields.push((key.into(), value));
    }
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else {
        panic!("expected object")
    };
    fields.retain(|(name, _)| name != key);
}
fn array_mut(root: &mut J) -> &mut Vec<J> {
    let J::Array(values) = root else {
        panic!("expected array")
    };
    values
}
fn text(value: &J) -> String {
    match value {
        J::Null => "null".into(),
        J::Bool(v) => v.to_string(),
        J::Number { raw, .. } => raw.clone(),
        J::Str(v) => format!(
            "\"{}\"",
            v.replace('\\', "\\\\")
                .replace('"', "\\\"")
                .replace('\n', "\\n")
        ),
        J::Array(v) => format!("[{}]", v.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(v) => format!(
            "{{{}}}",
            v.iter()
                .map(|(k, v)| format!("{}:{}", text(&J::Str(k.clone())), text(v)))
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
}
fn output(request: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(text(request).as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}
fn run(request: &J) -> J {
    let result = output(request);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    J::parse(std::str::from_utf8(&result.stdout).unwrap()).unwrap()
}
fn n(root: &J, key: &str) -> f64 {
    root.f64_field(key).unwrap()
}
fn at<'a>(root: &'a J, path: &[&str]) -> &'a J {
    root.path(path).unwrap()
}
fn vector(root: &J, key: &str) -> Vec<f64> {
    root.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-5 * expected.abs().max(1.0),
        "adjoint {actual:.12e}, independent perturbed trajectory {expected:.12e}"
    );
}
fn final_temperature(result: &J) -> f64 {
    n(result.get("objective").unwrap(), "value_k")
}
fn peak(result: &J) -> f64 {
    n(result.get("transient").unwrap(), "sampled_peak_objective_k")
}

fn request(qoi: Option<&str>) -> J {
    let mut result = parse(FIXTURE);
    // At pressures [2,2,1,0], resistance=1/Q^2 gives capacity rates 5:2->7.
    // Only the first and last branches exchange heat; the cold bypass mixes
    // before the last port. These are numerical fixtures, not measured data.
    put(
        &mut result,
        "hydraulics",
        parse(
            r#"{
      "node_count":4,
      "boundaries":[{"node":0,"pressure_pa":2,"temperature_k":330},
                    {"node":1,"pressure_pa":2,"temperature_k":290},
                    {"node":3,"pressure_pa":0}],
      "branches":[
        {"name":"first","from":0,"to":2,"resistance_pa_s2_m6":0.04,"source":"synthetic first path","regions":["first"]},
        {"name":"bypass","from":1,"to":2,"resistance_pa_s2_m6":0.25,"source":"synthetic bypass","regions":[]},
        {"name":"last","from":2,"to":3,"resistance_pa_s2_m6":0.02040816326530612,"source":"synthetic last path","regions":["last"]}
      ]}"#,
        ),
    );
    put(
        member(&mut result, "solid"),
        "surfaces",
        parse(
            r#"[
      {"name":"first","faces":[[0,2,3]],"htc_w_m2_k":20},
      {"name":"last","faces":[[1,2,3],[0,1,3],[0,1,2]],"htc_w_m2_k":10}
    ]"#,
        ),
    );
    let radiation = member(&mut result, "radiation");
    let rows = array_mut(member(radiation, "surfaces"));
    put(&mut rows[0], "surface", J::Str("first".into()));
    let schedule = member(&mut result, "transient");
    put(
        schedule,
        "intervals",
        parse(
            r#"[
      {"duration_s":0.4,"power_scale":1},
      {"duration_s":0.4,"power_scale":6}
    ]"#,
        ),
    );
    let chart = member(schedule, "enthalpy");
    remove(chart, "initial_specific_enthalpy_j_kg");
    put(
        chart,
        "initial_specific_enthalpies_j_kg",
        J::Array(INITIAL.map(number).to_vec()),
    );
    let knots = array_mut(member(chart, "knots"));
    put(&mut knots[3], "specific_enthalpy_j_kg", number(10000.0));
    put(&mut knots[3], "temperature_k", number(1050.0));
    if let Some(qoi) = qoi {
        put(
            schedule,
            "adjoint",
            parse(&format!(
                r#"{{"qoi":"{qoi}","max_checkpoint_bytes":1048576}}"#
            )),
        );
    }
    result
}

fn scale_interval(request: &mut J, interval: usize, multiplier: f64) {
    let rows = array_mut(member(member(request, "transient"), "intervals"));
    let value = n(&rows[interval], "power_scale");
    put(
        &mut rows[interval],
        "power_scale",
        number(value * multiplier),
    );
}
fn shift_initial(request: &mut J, vertex: Option<usize>, delta: f64) {
    let chart = member(member(request, "transient"), "enthalpy");
    let values = array_mut(member(chart, "initial_specific_enthalpies_j_kg"));
    for (index, value) in values.iter_mut().enumerate() {
        if vertex.is_none_or(|v| v == index) {
            *value = number(value.as_f64().unwrap() + delta);
        }
    }
}
fn scaled(multiplier: f64) -> J {
    let mut base = request(None);
    for interval in 0..2 {
        scale_interval(&mut base, interval, multiplier);
    }
    base
}

#[test]
fn final_gradient_carries_latent_history_through_mixed_air_and_radiation() {
    let base = request(None);
    let plain = run(&base);
    let differentiated = run(&request(Some("final")));
    for path in [
        &["solid_specific_enthalpies_j_kg"][..],
        &["solid_temperatures_k"][..],
        &["transient", "history"][..],
    ] {
        assert_eq!(at(&plain, path), at(&differentiated, path));
    }
    let trajectory = plain.get("transient").unwrap();
    let history = trajectory.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.len(), 5);
    for row in &history[1..=2] {
        assert_eq!(n(row, "objective_temperature_k"), 350.0);
        assert!((1000.0..3000.0).contains(&n(row, "minimum_specific_enthalpy_j_kg")));
        assert!((1000.0..3000.0).contains(&n(row, "maximum_specific_enthalpy_j_kg")));
    }
    assert_ne!(
        n(&history[1], "mean_liquid_mass_fraction"),
        n(&history[2], "mean_liquid_mass_fraction")
    );
    assert!(final_temperature(&plain) > 360.0);
    assert!(n(trajectory, "radiative_energy_loss_j") > 0.0);
    let gradient = at(&differentiated, &["transient", "adjoint"]);
    assert_eq!(
        gradient.str_field("method"),
        Some("discrete-backward-euler-coupled-enthalpy-adjoint")
    );
    close(n(gradient, "value_k"), final_temperature(&plain));
    assert_eq!(n(gradient, "state_index"), 4.0);
    close(n(gradient, "time_s"), 0.8);
    let initial = vector(gradient, "dtemperature_dinitial_specific_enthalpies_k_kg_j");
    assert_eq!(initial.len(), 4);
    let uniform = n(
        gradient,
        "dtemperature_duniform_initial_specific_enthalpy_k_kg_j",
    );
    close(uniform, initial.iter().sum());
    assert!(
        uniform > 0.01,
        "latent history must affect later sensible temperature"
    );
    for (i, &derivative) in initial.iter().enumerate() {
        let mut plus = base.clone();
        let mut minus = base.clone();
        shift_initial(&mut plus, Some(i), 0.1);
        shift_initial(&mut minus, Some(i), -0.1);
        close(
            derivative,
            (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.2,
        );
    }
    let mut plus = base.clone();
    let mut minus = base.clone();
    shift_initial(&mut plus, None, 0.1);
    shift_initial(&mut minus, None, -0.1);
    close(
        uniform,
        (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.2,
    );
    let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
    assert_eq!(intervals.len(), 2);
    assert!(
        n(&intervals[0], "dtemperature_dpower_multiplier_k") > 1.0,
        "heat deposited entirely during the plateau must influence the final sensible field"
    );
    for (i, row) in intervals.iter().enumerate() {
        let mut plus = base.clone();
        let mut minus = base.clone();
        scale_interval(&mut plus, i, 1.001);
        scale_interval(&mut minus, i, 0.999);
        close(
            n(row, "dtemperature_dpower_multiplier_k"),
            (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.002,
        );
    }
    let inlets = vector(gradient, "dtemperature_dinlet_temperatures");
    assert_eq!(inlets.len(), 4);
    for i in 0..2 {
        let mut plus = base.clone();
        let mut minus = base.clone();
        for (request, delta) in [(&mut plus, 0.02), (&mut minus, -0.02)] {
            let rows = array_mut(member(member(request, "hydraulics"), "boundaries"));
            let value = n(&rows[i], "temperature_k");
            put(&mut rows[i], "temperature_k", number(value + delta));
        }
        close(
            inlets[i],
            (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.04,
        );
    }
    assert_eq!(inlets[2..], [0.0, 0.0]);
}

#[test]
fn power_sizing_uses_the_actual_sampled_peak_adjoint_and_returns_a_passing_replay() {
    let target = run(&scaled(0.9));
    let limit = peak(&target);
    assert!(
        limit > 350.0,
        "target must be beyond latent storage, away from a corner"
    );
    let mut request = request(Some("sampled-peak"));
    let schedule = member(&mut request, "transient");
    put(schedule, "temperature_limit_k", number(limit));
    put(
        schedule,
        "power_design",
        parse(
            r#"{
      "min_power_multiplier":0,"max_power_multiplier":1.6,
      "power_multiplier_tolerance":0.0001,"temperature_tolerance_k":0.0001,
      "max_evaluations":48
    }"#,
        ),
    );
    let result = run(&request);
    let design = result.get("transient_power_design").unwrap();
    assert_eq!(
        design.str_field("search_method"),
        Some("safeguarded-adjoint-newton-bisection")
    );
    assert!(n(design, "newton_trials") > 0.0);
    let selected = n(design, "selected_power_multiplier");
    assert!((selected - 0.9).abs() < 0.01);
    assert!(peak(&result) <= limit);
    assert!(limit - peak(&result) <= 0.0001);
    assert!(n(design, "multiplier_bracket_width") <= 0.0001);
    assert!(
        n(
            design.get("failed_upper").unwrap(),
            "sampled_peak_objective_k"
        ) > limit
    );
    let replay = run(&scaled(selected));
    assert_eq!(
        result.get("solid_specific_enthalpies_j_kg"),
        replay.get("solid_specific_enthalpies_j_kg")
    );
    assert_eq!(
        result.get("solid_temperatures_k"),
        replay.get("solid_temperatures_k")
    );
    assert_eq!(
        at(&result, &["transient", "history"]),
        at(&replay, &["transient", "history"])
    );
    let adjoint = at(&result, &["transient", "adjoint"]);
    close(n(adjoint, "value_k"), peak(&result));
    let relative: f64 = adjoint
        .get("intervals")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|row| n(row, "dtemperature_dpower_multiplier_k"))
        .sum();
    let derivative = relative / selected;
    let trials = design.get("history").unwrap().as_array().unwrap();
    let accepted = trials
        .iter()
        .find(|row| n(row, "power_multiplier") == selected)
        .unwrap();
    close(n(accepted, "dpeak_dmultiplier_k"), derivative);
    let epsilon = 1e-4;
    close(
        derivative,
        (peak(&run(&scaled(selected + epsilon))) - peak(&run(&scaled(selected - epsilon))))
            / (2.0 * epsilon),
    );
    let repeat = run(&request);
    assert_eq!(
        at(&result, &["transient", "adjoint"]),
        at(&repeat, &["transient", "adjoint"])
    );
}

#[test]
fn mean_wall_replay_and_an_initial_peak_use_the_selected_physical_branch() {
    let mut mean = request(None);
    put(
        &mut mean,
        "objective",
        parse(r#"{"mean_wall_region":"last","gradient":false}"#),
    );
    let forward = run(&mean);
    let mut differentiated = mean.clone();
    put(
        member(&mut differentiated, "transient"),
        "adjoint",
        parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#),
    );
    let result = run(&differentiated);
    let gradient = at(&result, &["transient", "adjoint"]);
    assert_eq!(
        n(gradient, "value_k").to_bits(),
        final_temperature(&forward).to_bits()
    );
    let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
    let mut plus = mean.clone();
    let mut minus = mean;
    scale_interval(&mut plus, 0, 1.001);
    scale_interval(&mut minus, 0, 0.999);
    close(
        n(&intervals[0], "dtemperature_dpower_multiplier_k"),
        (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.002,
    );

    let mut cooling = scaled(0.0);
    shift_initial(&mut cooling, None, 3000.0);
    put(
        member(&mut cooling, "transient"),
        "adjoint",
        parse(r#"{"qoi":"sampled-peak","max_checkpoint_bytes":1048576}"#),
    );
    let result = run(&cooling);
    let gradient = at(&result, &["transient", "adjoint"]);
    assert_eq!(n(gradient, "state_index"), 0.0);
    assert_eq!(n(gradient, "active_vertex"), 1.0);
    assert_eq!(n(gradient, "reconstructed_solid_endpoints"), 0.0);
    // The selected initial vertex is strictly inside the liquid segment:
    // (1050-350)/(10000-3000) = 0.1 K per J/kg, independent of later cooling.
    let initial = vector(gradient, "dtemperature_dinitial_specific_enthalpies_k_kg_j");
    assert_eq!(initial, [0.0, 0.1, 0.0, 0.0]);
    assert_eq!(
        vector(gradient, "dtemperature_dinlet_temperatures"),
        [0.0; 4]
    );
    for row in gradient.get("intervals").unwrap().as_array().unwrap() {
        assert_eq!(n(row, "dtemperature_dpower_multiplier_k"), 0.0);
    }
}

#[test]
fn exhausted_checkpoints_and_unsupported_controls_publish_no_gradient() {
    for control in ["checkpoint", "component_power", "contact_resistance"] {
        let mut request = request(Some("final"));
        let adjoint = member(member(&mut request, "transient"), "adjoint");
        if control == "checkpoint" {
            put(adjoint, "max_checkpoint_bytes", number(1.0));
        } else {
            put(adjoint, control, J::Bool(true));
        }
        let result = output(&request);
        assert!(
            !result.status.success(),
            "unexpected acceptance of {control}"
        );
        assert!(result.stdout.is_empty());
        assert!(!result.stderr.is_empty());
    }
}
