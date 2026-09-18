#![cfg(unix)]

//! Actual binary target searches checked against independent two-face equations.
//! The analytic enclosure below is a synthetic supplied exchange matrix; it is
//! not evidence that view factors have been inferred from the tetrahedral mesh.

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/size-radiative-slab.json"));
const SIGMA: f64 = 5.670_374_419e-8;
const LIMIT: f64 = 323.0;

fn put(root: &mut J, key: &str, value: J) {
    let J::Object(fields) = root else { panic!("object required") };
    if let Some((_, slot)) = fields.iter_mut().find(|(name, _)| name == key) { *slot = value; }
    else { fields.push((key.to_string(), value)); }
}
fn member<'a>(root: &'a mut J, key: &str) -> &'a mut J {
    let J::Object(fields) = root else { panic!("object required") };
    &mut fields.iter_mut().find(|(name, _)| name == key).unwrap().1
}
fn remove(root: &mut J, key: &str) {
    let J::Object(fields) = root else { panic!("object required") };
    fields.retain(|(name, _)| name != key);
}
fn number(value: f64) -> J { J::Number { value, raw: value.to_string() } }
fn text(value: &J) -> String {
    match value {
        J::Null => "null".into(), J::Bool(value) => value.to_string(),
        J::Number { raw, .. } => raw.clone(),
        J::Str(value) => format!("\"{}\"", value.replace('\\', "\\\\")
            .replace('"', "\\\"").replace('\n', "\\n")),
        J::Array(values) => format!("[{}]", values.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(values) => format!("{{{}}}", values.iter().map(|(key, value)|
            format!("{}:{}", text(&J::Str(key.clone())), text(value))).collect::<Vec<_>>().join(",")),
    }
}
fn output(root: &J) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", "/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(root).as_bytes()).unwrap();
    child.wait_with_output().unwrap()
}
fn run(root: &J) -> J {
    let result = output(root);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    J::parse(&String::from_utf8(result.stdout).unwrap()).unwrap()
}
fn n(root: &J, key: &str) -> f64 { root.f64_field(key).unwrap() }
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} versus {expected}");
}

/// Eliminate the air references analytically and Newton-solve two temperatures.
/// This uses neither the product's FEM, secant iteration, nor adjoint routines.
fn reference(h: f64, enclosure: bool) -> [f64; 2] {
    let capacity = [0.003 * 1.2 * 1007.0, 0.004 * 1.2 * 1007.0];
    let effectiveness = [1.0 - (-0.5_f64 / capacity[0]).exp(),
        1.0 - (-h * 0.01 / capacity[1]).exp()];
    let g = [capacity[0] * effectiveness[0], capacity[1] * effectiveness[1]];
    let mut t = [320.0_f64, 320.0_f64];
    for _ in 0..64 {
        let mixed = 0.75 * (330.0 + effectiveness[0] * (t[0] - 330.0)) + 0.25 * 290.0;
        let (q, j) = if enclosure {
            let exchange = SIGMA * 0.01 / (1.0 / 0.85 + 1.0 / 0.6 - 1.0);
            let heat = exchange * (t[0].powi(4) - t[1].powi(4));
            let a = 4.0 * exchange * t[0].powi(3);
            let b = -4.0 * exchange * t[1].powi(3);
            ([heat, -heat], [a, b, -a, -b])
        } else {
            ([0.85 * SIGMA * 0.01 * (t[0].powi(4) - 300.0_f64.powi(4)),
              0.6 * SIGMA * 0.01 * (t[1].powi(4) - 310.0_f64.powi(4))],
             [4.0 * 0.85 * SIGMA * 0.01 * t[0].powi(3), 0.0,
              0.0, 4.0 * 0.6 * SIGMA * 0.01 * t[1].powi(3)])
        };
        let residual = [2.0 * (t[0] - t[1]) + g[0] * (t[0] - 330.0) + q[0],
            2.0 * (t[1] - t[0]) + g[1] * (t[1] - mixed) + q[1]];
        if residual[0].abs().max(residual[1].abs()) < 1e-11 { return t; }
        let a = 2.0 + g[0] + j[0];
        let b = -2.0 + j[1];
        let c = -2.0 - g[1] * 0.75 * effectiveness[0] + j[2];
        let d = 2.0 + g[1] + j[3];
        let determinant = a * d - b * c;
        t[0] -= (d * residual[0] - b * residual[1]) / determinant;
        t[1] -= (a * residual[1] - c * residual[0]) / determinant;
    }
    panic!("independent two-face equations did not converge")
}

fn reference_target(enclosure: bool) -> f64 {
    let (mut low, mut high) = (10.0_f64, 1000.0_f64);
    assert!(reference(low, enclosure)[0] > LIMIT);
    assert!(reference(high, enclosure)[0] <= LIMIT);
    for _ in 0..64 {
        let middle = (low * high).sqrt();
        if middle <= low || middle >= high { break; }
        if reference(middle, enclosure)[0] <= LIMIT { high = middle; }
        else { low = middle; }
    }
    high
}

fn check_design(input: &J, result: &J, enclosure: bool) {
    let design = result.get("design").unwrap();
    let h = n(design, "selected_htc_w_m2_k");
    let objective = n(result.get("objective").unwrap(), "value_k");
    assert_eq!(design.str_field("status"), Some("target-bracketed"));
    assert!(objective <= LIMIT && LIMIT - objective <= 1e-5);
    assert!(n(design, "log_bracket_width") <= 1e-4);
    assert!(n(design.get("failed_lower").unwrap(), "temperature_k") > LIMIT);
    near(h, reference_target(enclosure), 0.01);
    let expected = reference(h, enclosure);
    near(objective, expected[0], 2e-6);
    let positions = input.path(&["solid", "vertices_m"]).unwrap().as_array().unwrap();
    let field = result.get("solid_temperatures_k").unwrap().as_array().unwrap();
    assert_eq!(field.len(), positions.len());
    for (value, position) in field.iter().zip(positions) {
        let x = position.as_array().unwrap()[0].as_f64().unwrap();
        near(value.as_f64().unwrap(), expected[0] + (expected[1] - expected[0]) * x / 0.05, 2e-6);
    }
    let wall = result.get("walls").unwrap().as_array().unwrap().iter()
        .find(|wall| wall.str_field("region") == Some("last-face")).unwrap();
    near(n(wall, "htc_w_m2_k"), h, 0.0);
    let delta = 1e-4_f64;
    let derivative = (reference(h * delta.exp(), enclosure)[0]
        - reference(h * (-delta).exp(), enclosure)[0]) / (2.0 * delta);
    near(n(wall, "dobjective_dlog_htc"), derivative, 2e-5);
    let report = result.get("radiation").unwrap();
    assert_ne!(report.get("adjoint"), Some(&J::Null));
    assert!(report.get("adjoint").is_some());
    let selected_work = n(report, if enclosure { "total_solid_solves" } else { "solid_solves" });
    assert!(n(design, "total_solid_solves") > selected_work);
    assert!(selected_work > n(result, "coupling_iterations"));
    near(n(design, "evaluations"), design.get("history").unwrap().as_array().unwrap().len() as f64, 0.0);
}

#[test]
fn ambient_target_search_matches_independent_temperature_field_and_total_adjoint() {
    let input = J::parse(FIXTURE).unwrap();
    let result = run(&input);
    check_design(&input, &result, false);
    let report = result.get("radiation").unwrap();
    near(n(report, "radiative_out_w") + n(report, "convective_out_w"), n(&result, "source_w"), 1e-7);
    assert!(n(report, "energy_residual_w").abs() <= 1e-7);
    assert!(n(report, "radiative_out_w") > 1.0);

    // Re-run the selected h without a search. Its field and mechanism report
    // must be the same evaluation, not a report retained from a rejected trial.
    let mut replay = input.clone();
    let selected = n(result.get("design").unwrap(), "selected_htc_w_m2_k");
    remove(&mut replay, "design");
    let J::Array(surfaces) = member(member(&mut replay, "solid"), "surfaces") else { panic!() };
    let target = surfaces.iter_mut().find(|surface| surface.str_field("name") == Some("last-face")).unwrap();
    put(target, "htc_w_m2_k", number(selected));
    let replayed = run(&replay);
    assert_eq!(result.get("solid_temperatures_k"), replayed.get("solid_temperatures_k"));
    assert_eq!(result.get("radiation"), replayed.get("radiation"));
}

#[test]
fn enclosure_target_search_keeps_reflection_and_internal_heat_accounting() {
    let mut input = J::parse(FIXTURE).unwrap();
    let radiation = member(&mut input, "radiation");
    remove(radiation, "surfaces");
    put(radiation, "enclosure", J::parse(r#"{
        "surfaces": [
            {"surface":"first-face","emissivity":0.85,"source":"synthetic equal-area exchange fixture"},
            {"surface":"last-face","emissivity":0.6,"source":"synthetic equal-area exchange fixture"}
        ],
        "view_factors": [[0,1],[1,0]],
        "row_sum_tolerance": 1e-12,
        "reciprocity_relative_tolerance": 1e-12,
        "evidence": {"kind":"analytic","geometry":"declared two-patch unit-view-factor reference; not inferred from the mesh"}
    }"#).unwrap());
    let result = run(&input);
    check_design(&input, &result, true);
    let report = result.get("radiation").unwrap();
    assert_eq!(report.str_field("model"), Some("closed-gray-diffuse-enclosure"));
    near(n(report, "radiative_out_w"), 0.0, 1e-7);
    let rows = report.get("surfaces").unwrap().as_array().unwrap();
    assert!(n(&rows[0], "outward_heat_w").abs() > 0.01);
    near(rows.iter().map(|row| n(row, "outward_heat_w")).sum(), 0.0, 1e-7);
    let air: f64 = result.get("walls").unwrap().as_array().unwrap().iter()
        .map(|wall| n(wall, "outward_heat_w")).sum();
    near(air, n(&result, "source_w"), 1e-7);
}

#[test]
fn radiative_minimum_feasibility_needs_only_one_complete_evaluation() {
    let mut input = J::parse(FIXTURE).unwrap();
    let design = member(&mut input, "design");
    put(design, "mean_temperature_limit_k", number(400.0));
    put(design, "max_evaluations", number(1.0));
    let result = run(&input);
    let design = result.get("design").unwrap();
    assert_eq!(design.str_field("status"), Some("minimum-feasible"));
    near(n(design, "selected_htc_w_m2_k"), 10.0, 0.0);
    near(n(design, "evaluations"), 1.0, 0.0);
    assert_eq!(design.get("failed_lower"), Some(&J::Null));
    near(n(design, "total_solid_solves"), n(result.get("radiation").unwrap(), "solid_solves"), 0.0);
}

#[test]
fn radiative_design_budget_and_missing_bracket_never_publish_partial_results() {
    for (object, key, value, exit, code) in [
        ("design", "max_evaluations", number(1.0), 6, "cooling-network-design-budget"),
        ("radiation", "max_iterations", number(1.0), 6, "cooling-network-radiation-budget"),
        ("design", "mean_temperature_limit_k", number(290.0), 4, "cooling-network-design-bracket"),
        ("objective", "gradient", J::Bool(false), 4, "cooling-network-input"),
        ("budgets", "wall_seconds", number(1e-12), 6, "cooling-network-time-budget"),
    ] {
        let mut input = J::parse(FIXTURE).unwrap();
        put(member(&mut input, object), key, value);
        let result = output(&input);
        assert_eq!(result.status.code(), Some(exit), "{}", String::from_utf8_lossy(&result.stderr));
        assert!(result.stdout.is_empty());
        let diagnostic = J::parse(&String::from_utf8(result.stderr).unwrap()).unwrap();
        assert_eq!(diagnostic.str_field("code"), Some(code));
    }
}

#[test]
fn radiative_search_replays_exactly_under_patch_permutation() {
    let mut input = J::parse(FIXTURE).unwrap();
    let first = output(&input);
    assert!(first.status.success(), "{}", String::from_utf8_lossy(&first.stderr));
    let J::Array(patches) = member(member(&mut input, "radiation"), "surfaces") else { panic!() };
    patches.reverse();
    let second = output(&input);
    assert!(second.status.success(), "{}", String::from_utf8_lossy(&second.stderr));
    assert_eq!(first.stdout, second.stdout);
}
