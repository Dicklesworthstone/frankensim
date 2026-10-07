//! `frankensim cooling-cht`: the voxel conjugate heat-transfer workflow
//! through the real binary — forced and natural convection on small scenes,
//! energy closure, and structured refusals.

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;

use json::JsonValue as J;
use std::path::PathBuf;
use std::process::Command;

fn scratch(name: &str, scene: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fs-cooling-cht-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, scene).unwrap();
    path
}

fn run(path: &PathBuf) -> (i32, J, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-cht"])
        .arg(path)
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    let document = if stdout.trim().is_empty() {
        J::parse(stderr.trim()).unwrap()
    } else {
        J::parse(stdout.trim()).unwrap()
    };
    (output.status.code().unwrap(), document, stderr)
}

fn f(document: &J, path: &[&str]) -> f64 {
    document.path(path).and_then(J::as_f64).unwrap()
}

const DUCT: &str = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.012, 0.004, 0.004], "voxel_m": 0.001,
 "materials": [{"name": "aluminium", "conductivity_w_m_k": 167.0}],
 "solids": [{"material": "aluminium", "min_m": [0.004, 0.0, 0.0], "max_m": [0.008, 0.002, 0.004]}],
 "sources": [{"name": "chip", "power_w": 0.1, "min_m": [0.004, 0.0, 0.0], "max_m": [0.008, 0.001, 0.004]}],
 "faces": {
  "x-": {"type": "inlet", "velocity_m_s": [0.1, 0.0, 0.0], "temperature_k": 300.0},
  "x+": {"type": "opening", "ambient_k": 300.0}
 },
 "solver": {"tolerance": 1e-8}
}"#;

#[test]
fn forced_convection_scene_closes_energy_through_the_binary() {
    let (code, result, stderr) = run(&scratch("duct.json", DUCT));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.str_field("status"), Some("completed"));
    assert_eq!(result.str_field("solver"), Some("fv-simplec"));
    // Adiabatic duct walls: every watt leaves through the open faces, almost
    // all by advection (a little diffuses back out of the Dirichlet inlet).
    let power = f(&result, &["energy", "source_w"]);
    assert!((power - 0.1).abs() < 1e-12);
    let leaving = f(&result, &["energy", "boundary_outflow_w"]);
    assert!(
        (leaving - power).abs() < 1e-9 * power,
        "{leaving} vs {power}"
    );
    let advected = f(&result, &["energy", "advective_outflow_w"]);
    assert!(
        advected > 0.999 * power && advected <= leaving,
        "{advected}"
    );
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    let (inflow, outflow) = (
        f(&result, &["flow", "inflow_m3_s"]),
        f(&result, &["flow", "outflow_m3_s"]),
    );
    assert!((inflow - 0.1 * 16e-6).abs() < 1e-15, "{inflow}");
    assert!((outflow - inflow).abs() < 1e-12 * inflow);
    let hottest = f(&result, &["max_solid_temperature_k"]);
    assert!(hottest > 300.0, "{hottest}");
    let sources = result.get("sources").and_then(J::as_array).unwrap();
    assert_eq!(sources[0].path(&["cells"]).and_then(J::as_f64), Some(16.0));
}

#[test]
fn natural_convection_scene_draws_air_through_its_openings() {
    let scene = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.004, 0.004, 0.016], "voxel_m": 0.001,
 "faces": {
  "x-": {"type": "wall", "temperature_k": 320.0},
  "z-": {"type": "opening", "ambient_k": 300.0},
  "z+": {"type": "opening", "ambient_k": 300.0}
 },
 "gravity_m_s2": [0.0, 0.0, -9.81],
 "solver": {"tolerance": 1e-7}
}"#;
    let (code, result, stderr) = run(&scratch("chimney.json", scene));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.str_field("solver"), Some("fv-simplec-boussinesq"));
    let outflow = f(&result, &["flow", "outflow_m3_s"]);
    assert!(outflow > 0.0, "a heated wall must draw air upward");
    assert!((f(&result, &["flow", "inflow_m3_s"]) - outflow).abs() < 1e-12 * outflow);
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    // The heated wall's heat leaves by advection through the top opening.
    assert!(f(&result, &["energy", "advective_outflow_w"]) > 0.0);
}

#[test]
fn malformed_scenes_refuse_with_structured_codes() {
    let cases = [
        (DUCT.replace("cooling-cht.v1", "cooling-cht.v9"), "schema"),
        (
            DUCT.replace("[0.012, 0.004", "[0.0125, 0.004"),
            "whole number",
        ),
        (
            DUCT.replace("\"material\": \"aluminium\"", "\"material\": \"copper\""),
            "unknown material",
        ),
        (
            DUCT.replace(
                "\"min_m\": [0.004, 0.0, 0.0], \"max_m\": [0.008, 0.001",
                "\"min_m\": [0.0, 0.0, 0.0], \"max_m\": [0.001, 0.001",
            ),
            "covers no solid",
        ),
        (
            DUCT.replace("\"type\": \"opening\"", "\"type\": \"vent\""),
            "type must be",
        ),
    ];
    for (index, (scene, needle)) in cases.iter().enumerate() {
        assert_ne!(
            scene.as_str(),
            DUCT,
            "case {index} did not mutate the scene"
        );
        let (code, diagnostic, _) = run(&scratch(&format!("bad-{index}.json"), scene));
        assert_eq!(code, 4, "case {index}");
        assert_eq!(diagnostic.str_field("code"), Some("cooling-cht-input"));
        let message = diagnostic.str_field("message").unwrap();
        assert!(message.contains(needle), "case {index}: {message}");
    }
}

/// ASCII STL of the unit cube [0, 1]^3, outward-wound.
fn unit_cube_stl() -> String {
    let v = |i: usize| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64];
    // Two triangles per face, counter-clockwise seen from outside.
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1], // z = 0
        [4, 5, 7, 6], // z = 1
        [0, 1, 5, 4], // y = 0
        [2, 6, 7, 3], // y = 1
        [0, 4, 6, 2], // x = 0
        [1, 3, 7, 5], // x = 1
    ];
    let mut out = String::from("solid cube\n");
    for [a, b, c, d] in faces {
        for tri in [[a, b, c], [a, c, d]] {
            out.push_str(" facet normal 0 0 0\n  outer loop\n");
            for i in tri {
                let p = v(i);
                out.push_str(&format!("   vertex {} {} {}\n", p[0], p[1], p[2]));
            }
            out.push_str("  endloop\n endfacet\n");
        }
    }
    out.push_str("endsolid cube\n");
    out
}

#[test]
fn stl_solids_voxelize_by_winding_number_with_scale_and_offset() {
    // A 4 mm cube from a unit-cube STL (scale 0.004, offset 2 mm) in an
    // 8 mm box at 1 mm voxels occupies exactly 4^3 voxels; a heated cube in
    // still air with one cold wall conducts its power out through that wall.
    scratch("cube.stl", &unit_cube_stl());
    let scene = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.008, 0.008, 0.008], "voxel_m": 0.001,
 "materials": [{"name": "copper", "conductivity_w_m_k": 400.0}],
 "solids": [{"material": "copper", "stl": "cube.stl", "scale": 0.004, "offset_m": [0.002, 0.002, 0.002]}],
 "sources": [{"name": "die", "power_w": 0.01, "min_m": [0.0, 0.0, 0.0], "max_m": [0.008, 0.008, 0.008]}],
 "faces": {"z-": {"type": "wall", "temperature_k": 300.0}}
}"#;
    let (code, result, stderr) = run(&scratch("cube.json", scene));
    assert_eq!(code, 0, "{stderr}");
    let materials = result.get("materials").and_then(J::as_array).unwrap();
    assert_eq!(
        materials[0].path(&["cells"]).and_then(J::as_f64),
        Some(64.0)
    );
    let sources = result.get("sources").and_then(J::as_array).unwrap();
    assert_eq!(sources[0].path(&["cells"]).and_then(J::as_f64), Some(64.0));
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    // No flow anywhere: the power leaves through the cold wall.
    assert!((f(&result, &["energy", "boundary_outflow_w"]) - 0.01).abs() < 1e-11);
    assert!(f(&result, &["max_solid_temperature_k"]) > 300.0);
}

#[test]
fn fan_face_finds_its_operating_point_on_the_curve() {
    // The duct scene driven by a fan curve instead of a fixed velocity: the
    // reported operating point lies on the curve and the delivered flow is
    // the inflow.
    let scene = DUCT.replace(
        r#""x-": {"type": "inlet", "velocity_m_s": [0.1, 0.0, 0.0], "temperature_k": 300.0}"#,
        r#""x-": {"type": "fan", "curve": [[0.0, 0.05], [4e-6, 0.0]], "temperature_k": 300.0}"#,
    );
    assert_ne!(scene, DUCT);
    let (code, result, stderr) = run(&scratch("fan.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    let q = f(&result, &["flow", "fan_flow_m3_s"]);
    let dp = f(&result, &["flow", "fan_pressure_pa"]);
    assert!(q > 0.0 && q < 4e-6, "{q}");
    // Linear curve: dp = 0.05 (1 - q / 4e-6).
    assert!((dp - 0.05 * (1.0 - q / 4e-6)).abs() < 1e-12, "{dp}");
    assert!(f(&result, &["flow", "fan_residual"]) < 1e-6);
    assert!((f(&result, &["flow", "inflow_m3_s"]) - q).abs() < 1e-12 * q);
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
}

#[test]
fn orthotropic_board_and_interface_resistance_reach_the_solver() {
    // A die (k 150) on a laminate (k 30 in-plane, 0.3 through) with a
    // 2e-4 m^2K/W interface, cooled from below: the joint and the weak
    // through-board conduction both raise the die temperature, and the
    // contact can only make it hotter.
    let base = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.008, 0.008, 0.004], "voxel_m": 0.001,
 "materials": [
  {"name": "pcb", "conductivity_w_m_k": [30.0, 30.0, 0.3]},
  {"name": "die", "conductivity_w_m_k": 150.0}
 ],
 "solids": [
  {"material": "pcb", "min_m": [0.0, 0.0, 0.0], "max_m": [0.008, 0.008, 0.002]},
  {"material": "die", "min_m": [0.003, 0.003, 0.002], "max_m": [0.005, 0.005, 0.003]}
 ],
 CONTACTS
 "sources": [{"name": "die", "power_w": 0.05, "min_m": [0.003, 0.003, 0.002], "max_m": [0.005, 0.005, 0.003]}],
 "faces": {"z-": {"type": "wall", "temperature_k": 300.0}}
}"#;
    let solve = |contacts: &str, name: &str| {
        let (code, result, stderr) = run(&scratch(name, &base.replace("CONTACTS", contacts)));
        assert_eq!(code, 0, "{stderr}");
        assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
        f(&result, &["max_solid_temperature_k"])
    };
    let bonded = solve("", "board-bonded.json");
    let joined = solve(
        r#""contacts": [{"between": ["die", "pcb"], "resistance_m2_k_w": 2e-4}],"#,
        "board-joint.json",
    );
    // All 0.05 W crossing the 4 mm^2 joint would add exactly
    // q'' R'' = 2.5 K; the still air around the die offers a parallel path
    // into the board surface, so the rise is bounded by that and measured
    // at 1.91 K.
    assert!(bonded > 302.0, "{bonded}");
    let rise = joined - bonded;
    assert!(rise > 1.5 && rise < 2.5, "{bonded} -> {joined}");
}

#[test]
fn transient_march_warms_toward_the_steady_junction_and_closes_energy() {
    // The duct scene with heat capacities and a transient block: a long
    // march from the inlet temperature settles on the steady solution's
    // hottest solid voxel; with the schedule at zero nothing heats.
    let with_capacity = DUCT.replace(
        r#""conductivity_w_m_k": 167.0}"#,
        r#""conductivity_w_m_k": 167.0, "volumetric_heat_capacity_j_m3_k": 2.4e6}"#,
    );
    assert_ne!(with_capacity, DUCT);
    let scene = with_capacity.replace(
        r#""solver": {"tolerance": 1e-8}"#,
        r#""solver": {"tolerance": 1e-8}, "transient": {"time_step_s": 2.0, "steps": 400}"#,
    );
    let (code, result, stderr) = run(&scratch("warm.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    let steady = f(&result, &["max_solid_temperature_k"]);
    let last = f(&result, &["transient", "final_max_solid_temperature_k"]);
    assert!(
        (last - steady).abs() < 1e-3 * (steady - 300.0),
        "{last} vs steady {steady}"
    );
    // Energy closure per step relative to the 0.2 J the die releases per
    // step (measured 4e-9 J: 2e-8 relative, the solver tolerance level).
    let closure = f(&result, &["transient", "worst_step_closure_j"]);
    assert!(closure < 1e-6 * 0.2, "{closure}");
    let records = result
        .path(&["transient", "records"])
        .and_then(J::as_array)
        .unwrap();
    let first = records[0]
        .path(&["max_solid_temperature_k"])
        .and_then(J::as_f64)
        .unwrap();
    assert!(first > 300.0 && first < last, "{first} {last}");
    // Power held off by the schedule: the body stays at the inlet temperature.
    let off = with_capacity.replace(
        r#""solver": {"tolerance": 1e-8}"#,
        r#""solver": {"tolerance": 1e-8}, "transient": {"time_step_s": 1.0, "steps": 5, "power_schedule": [[0.0, 0.0], [10.0, 0.0]]}"#,
    );
    let (code, result, stderr) = run(&scratch("off.json", &off));
    assert_eq!(code, 0, "{stderr}");
    let cold = f(&result, &["transient", "final_max_solid_temperature_k"]);
    assert!((cold - 300.0).abs() < 1e-9, "{cold}");
    // No heat capacity: refused, not guessed.
    let missing = DUCT.replace(
        r#""solver": {"tolerance": 1e-8}"#,
        r#""solver": {"tolerance": 1e-8}, "transient": {"time_step_s": 1.0, "steps": 5}"#,
    );
    let (code, diagnostic, _) = run(&scratch("nocap.json", &missing));
    assert_eq!(code, 4);
    assert!(
        diagnostic
            .str_field("message")
            .unwrap()
            .contains("heat capacity")
    );
}
