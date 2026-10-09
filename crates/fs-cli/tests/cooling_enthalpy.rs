//! G1/G3 actual-binary checks against an analytic latent-plateau balance.
//! The chart, geometry and air properties are synthetic numerical references.
#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;

use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const FIXTURE_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enthalpy-phase-pulse.json"
);
const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/enthalpy-phase-pulse.json"
));
const RHO: f64 = 10.0;
const SIGMA: f64 = 5.670_374_419e-8;
const EPSILON: f64 = 0.8;

fn number(value: f64) -> J {
    J::Number {
        value,
        raw: value.to_string(),
    }
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

fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else {
        panic!("expected object")
    };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}

fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else {
        panic!("expected object")
    };
    fields.retain(|(name, _)| name != key);
}

fn text(value: &J) -> String {
    match value {
        J::Null => "null".into(),
        J::Bool(value) => value.to_string(),
        J::Number { raw, .. } => raw.clone(),
        J::Str(value) => format!(
            "\"{}\"",
            value
                .replace('\\', "\\\\")
                .replace('"', "\\\"")
                .replace('\n', "\\n")
        ),
        J::Array(values) => format!(
            "[{}]",
            values.iter().map(text).collect::<Vec<_>>().join(",")
        ),
        J::Object(fields) => format!(
            "{{{}}}",
            fields
                .iter()
                .map(|(key, value)| format!("{}:{}", text(&J::Str(key.clone())), text(value)))
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

fn document(output: &Output) -> J {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}

fn run(request: &J) -> J {
    document(&output(request))
}
fn n(root: &J, key: &str) -> f64 {
    root.f64_field(key).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.14e} versus {b:.14e}");
}
fn values(root: &J, key: &str) -> Vec<f64> {
    root.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

// The right tetrahedron has V=1/6 and opposite face areas sqrt(3)/2,1/2,1/2,1/2.
// All T stay at 350 K inside the chart's latent plateau, even as h changes.
// For one unit-capacity-rate air stream, Q_air=50*(1-exp(-h*A)). Surface loads
// are (A-A_opposite)/3 times the uniform convective+radiative flux density.
// This oracle uses no FEM assembly, phase lookup, or cooling library routine.
fn heat(ambient: Option<f64>) -> (f64, f64, f64) {
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    let air = 50.0 * (1.0 - (-2.0 * area).exp());
    let radiation = ambient.map_or(0.0, |a| {
        EPSILON * SIGMA * area * (350.0_f64.powi(4) - a.powi(4))
    });
    (area, air, radiation)
}

fn enthalpies(time: f64, ambient: Option<f64>) -> [f64; 4] {
    let (area, air, radiation) = heat(ambient);
    let flux_density = (air + radiation) / area;
    let source_change = 6000.0 * time.min(0.4) / RHO;
    let mass = RHO / 24.0;
    [3.0_f64.sqrt() / 2.0, 0.5, 0.5, 0.5].map(|opposite| {
        2000.0 + source_change - time * ((area - opposite) / 3.0) * flux_density / mass
    })
}

fn check_plateau(result: &J, ambient: Option<f64>) {
    let (_, air, radiation) = heat(ambient);
    let trajectory = result.get("transient").unwrap();
    assert_eq!(
        trajectory.str_field("scheme"),
        Some("backward-euler-total-enthalpy")
    );
    near(n(trajectory, "time_s"), 0.8, 1e-14);
    assert_eq!(n(trajectory, "steps"), 4.0);
    near(n(trajectory, "sampled_peak_objective_k"), 350.0, 1e-8);
    near(n(trajectory, "input_energy_j"), 400.0, 1e-8);
    near(n(trajectory, "air_energy_gain_j"), 0.8 * air, 1e-7);
    near(
        n(trajectory, "stored_energy_change_j"),
        400.0 - 0.8 * (air + radiation),
        1e-6,
    );
    assert!(n(trajectory, "energy_residual_j").abs() <= 0.8e-7);
    if ambient.is_some() {
        near(
            n(trajectory, "radiative_energy_loss_j"),
            0.8 * radiation,
            1e-6,
        );
    } else {
        assert!(result.get("radiation").is_none());
    }
    let expected = enthalpies(0.8, ambient);
    let actual = values(result, "solid_specific_enthalpies_j_kg");
    let liquid = values(result, "solid_liquid_mass_fractions");
    assert_eq!(actual.len(), 4);
    assert_eq!(liquid.len(), 4);
    for i in 0..4 {
        assert!((1000.0..3000.0).contains(&expected[i]));
        near(actual[i], expected[i], 2e-6);
        near(liquid[i], (expected[i] - 1000.0) / 2000.0, 1e-9);
    }
    for temperature in values(result, "solid_temperatures_k") {
        near(temperature, 350.0, 1e-9);
    }
    let material = trajectory.get("enthalpy").unwrap();
    near(n(material, "reference_density_kg_m3"), RHO, 1e-14);
    near(
        n(material, "initial_total_enthalpy_j"),
        RHO / 6.0 * 2000.0,
        1e-9,
    );
    near(
        n(material, "final_total_enthalpy_j"),
        actual.iter().sum::<f64>() * RHO / 24.0,
        1e-8,
    );
    let history = trajectory.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.len(), 5);
    let mut storage = 0.0;
    for (index, row) in history.iter().enumerate() {
        let time = 0.2 * index as f64;
        let h = enthalpies(time, ambient);
        near(n(row, "time_s"), time, 1e-14);
        near(
            n(row, "minimum_specific_enthalpy_j_kg"),
            h.iter().copied().fold(f64::INFINITY, f64::min),
            2e-6,
        );
        near(
            n(row, "maximum_specific_enthalpy_j_kg"),
            h.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            2e-6,
        );
        near(
            n(row, "mean_liquid_mass_fraction"),
            h.iter().map(|v| (v - 1000.0) / 2000.0).sum::<f64>() / 4.0,
            1e-9,
        );
        if index > 0 {
            let source = if index <= 2 { 1000.0 } else { 0.0 };
            near(n(row, "dt_s"), 0.2, 1e-14);
            near(n(row, "source_w"), source, 1e-8);
            near(n(row, "air_heat_gain_w"), air, 1e-7);
            if ambient.is_some() {
                near(n(row, "radiative_heat_w"), radiation, 1e-6);
            }
            near(
                n(row, "stored_energy_change_j"),
                0.2 * (source - air - radiation),
                2e-7,
            );
            assert!(
                n(row, "coupling_iterations") > 1.0,
                "exercise multiple trials with one immutable history"
            );
            storage += n(row, "stored_energy_change_j");
        }
    }
    near(storage, n(trajectory, "stored_energy_change_j"), 1e-9);
    let wall = &result.get("walls").unwrap().as_array().unwrap()[0];
    near(n(wall, "outward_heat_w"), air, 1e-7);
    let branch = &result.get("branches").unwrap().as_array().unwrap()[0];
    near(n(branch, "flow_m3_s"), 1.0, 1e-12);
    near(n(branch, "outlet_k"), 300.0 + air, 1e-7);
}

#[test]
fn latent_storage_and_hot_cold_radiation_match_independent_endpoint_balances() {
    let cold = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", FIXTURE_PATH])
        .output()
        .unwrap();
    check_plateau(&document(&cold), Some(300.0));
    let input = J::parse(FIXTURE).unwrap();
    let replay = output(&input);
    document(&replay);
    assert_eq!(cold.stdout, replay.stdout);
    let mut hot = input.clone();
    let J::Array(patches) = member(member(&mut hot, "radiation"), "surfaces") else {
        panic!()
    };
    put(&mut patches[0], "ambient_temperature_k", number(400.0));
    check_plateau(&run(&hot), Some(400.0));
    let mut no_radiation = input;
    remove(&mut no_radiation, "radiation");
    check_plateau(&run(&no_radiation), None);
}

#[test]
fn final_nodal_enthalpy_restarts_the_same_physical_trajectory() {
    let input = J::parse(FIXTURE).unwrap();
    let full = run(&input);
    let mut first = input.clone();
    let J::Array(intervals) = member(member(&mut first, "transient"), "intervals") else {
        panic!()
    };
    intervals.truncate(1);
    let first = run(&first);
    let mut continuation = input;
    let schedule = member(&mut continuation, "transient");
    let J::Array(intervals) = member(schedule, "intervals") else {
        panic!()
    };
    intervals.remove(0);
    let storage = member(schedule, "enthalpy");
    remove(storage, "initial_specific_enthalpy_j_kg");
    put(
        storage,
        "initial_specific_enthalpies_j_kg",
        first.get("solid_specific_enthalpies_j_kg").unwrap().clone(),
    );
    let continued = run(&continuation);
    for key in [
        "solid_specific_enthalpies_j_kg",
        "solid_liquid_mass_fractions",
        "solid_temperatures_k",
    ] {
        for (a, b) in values(&full, key).iter().zip(values(&continued, key)) {
            near(*a, b, 1e-7);
        }
    }
    for key in [
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
        "radiative_energy_loss_j",
    ] {
        near(
            n(full.get("transient").unwrap(), key),
            n(first.get("transient").unwrap(), key) + n(continued.get("transient").unwrap(), key),
            1e-7,
        );
    }
    near(
        n(
            continued.path(&["transient", "enthalpy"]).unwrap(),
            "initial_total_enthalpy_j",
        ),
        n(
            first.path(&["transient", "enthalpy"]).unwrap(),
            "final_total_enthalpy_j",
        ),
        1e-8,
    );
}

fn refuses(input: &J) -> String {
    let output = output(input);
    assert!(!output.status.success());
    assert!(
        output.stdout.is_empty(),
        "a refused mode or endpoint must publish no partial trajectory"
    );
    String::from_utf8(output.stderr).unwrap()
}

#[test]
fn unsupported_storage_modes_and_invalid_or_exhausted_inputs_refuse() {
    let input = J::parse(FIXTURE).unwrap();
    for (key, value) in [
        (
            "adaptive",
            r#"{"absolute_tolerance_k":0.01,"relative_tolerance":0,"minimum_trial_step_s":0.01,"max_trials":100}"#,
        ),
        ("repeat", r#"{"cycles":2,"max_total_steps":16}"#),
        (
            "adjoint",
            r#"{"qoi":"sampled-peak","max_checkpoint_bytes":1048576}"#,
        ),
        (
            "power_design",
            r#"{"min_power_multiplier":0,"max_power_multiplier":1,"power_multiplier_tolerance":0.01,"temperature_tolerance_k":0.1,"max_evaluations":16}"#,
        ),
        (
            "fan_speed_design",
            r#"{"min_speed_multiplier":0.5,"max_speed_multiplier":1,"speed_multiplier_tolerance":0.01,"temperature_tolerance_k":0.1,"max_evaluations":16}"#,
        ),
        (
            "time_convergence",
            r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":0.01,"max_total_steps":100,"max_trace_bytes":1048576}"#,
        ),
    ] {
        let mut unsupported = input.clone();
        put(
            member(&mut unsupported, "transient"),
            key,
            J::parse(value).unwrap(),
        );
        let diagnostic = refuses(&unsupported);
        assert!(diagnostic.contains("enthalpy"), "{key}: {diagnostic}");
    }
    let mut gradient = input.clone();
    put(
        member(&mut gradient, "objective"),
        "gradient",
        J::Bool(true),
    );
    assert!(refuses(&gradient).contains("gradient"));
    for key in ["design", "fan_speed_design", "mesh_convergence"] {
        let mut unsupported = input.clone();
        put(&mut unsupported, key, J::Object(Vec::new()));
        assert!(
            refuses(&unsupported).contains("enthalpy"),
            "{key} must explicitly refuse the storage mode"
        );
    }
    let mut enclosure = input.clone();
    let policy = member(&mut enclosure, "radiation");
    remove(policy, "surfaces");
    put(policy, "enclosure", J::Object(Vec::new()));
    assert!(refuses(&enclosure).contains("enthalpy"));

    let mut duplicate = input.clone();
    put(
        member(&mut duplicate, "transient"),
        "initial_temperature_k",
        number(350.0),
    );
    assert!(refuses(&duplicate).contains("enthalpy"));
    let mut invalid = input.clone();
    put(
        member(member(&mut invalid, "transient"), "enthalpy"),
        "reference_density_kg_m3",
        number(-1.0),
    );
    refuses(&invalid);
    let mut work = input.clone();
    put(
        member(&mut work, "budgets"),
        "linear_iterations",
        number(32.0),
    );
    refuses(&work);
    let mut radiation = input.clone();
    put(
        member(&mut radiation, "radiation"),
        "max_iterations",
        number(1.0),
    );
    refuses(&radiation);
    let mut cancelled = input;
    put(
        member(&mut cancelled, "budgets"),
        "wall_seconds",
        number(1e-12),
    );
    let cancelled = output(&cancelled);
    assert_eq!(
        cancelled.status.code(),
        Some(i32::from(fs_cli::exit::BUDGET))
    );
    assert!(cancelled.stdout.is_empty());
}
