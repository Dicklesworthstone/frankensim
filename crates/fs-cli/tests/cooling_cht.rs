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
fn radiating_transient_settles_on_the_steady_radiating_solution() {
    // The duct scene with an emissive heatsink: the frozen-flow march
    // carries radiation through the inlet and opening, and a long march
    // settles on the steady radiating solve's temperature and radiated
    // power.
    let scene = DUCT
        .replace(
            r#""conductivity_w_m_k": 167.0}"#,
            r#""conductivity_w_m_k": 167.0, "volumetric_heat_capacity_j_m3_k": 2.4e6, "emissivity": 0.9}"#,
        )
        .replace(
            r#""solver": {"tolerance": 1e-8}"#,
            r#""solver": {"tolerance": 1e-8}, "transient": {"time_step_s": 2.0, "steps": 400}"#,
        );
    let (code, result, stderr) = run(&scratch("rad-warm.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    let steady = f(&result, &["max_solid_temperature_k"]);
    let steady_radiated = f(&result, &["radiation", "radiated_w"]);
    let last = f(&result, &["transient", "final_max_solid_temperature_k"]);
    assert!(
        (last - steady).abs() < 1e-3 * (steady - 300.0),
        "{last} vs steady {steady}"
    );
    let closure = f(&result, &["transient", "worst_step_closure_j"]);
    assert!(closure < 1e-6 * 0.2, "{closure}");
    let records = result
        .path(&["transient", "records"])
        .and_then(J::as_array)
        .unwrap();
    let radiated = records
        .last()
        .and_then(|r| r.path(&["radiated_w"]))
        .and_then(J::as_f64)
        .unwrap();
    assert!(steady_radiated > 0.0, "{steady_radiated}");
    assert!(
        (radiated - steady_radiated).abs() < 1e-3 * steady_radiated,
        "{radiated} vs steady {steady_radiated}"
    );
}

const CORNER_PLATE: &str = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.008, 0.008, 0.004], "voxel_m": 0.002,
 "fluid": {"density_kg_m3": 1.0, "specific_heat_j_kg_k": 1000.0,
           "conductivity_w_m_k": 1e-12, "kinematic_viscosity_m2_s": 1e-5},
 "materials": [{"name": "plate", "conductivity_w_m_k": 10.0}],
 "solids": [{"material": "plate", "min_m": [0.0, 0.0, 0.0], "max_m": [0.008, 0.008, 0.002]}],
 "sources": [{"name": "corner", "power_w": 0.05, "min_m": [0.006, 0.006, 0.0], "max_m": [0.008, 0.008, 0.002]}],
 "faces": {"x-": {"type": "wall", "temperature_k": 300.0}},
 "grid_convergence": {"splits": [1, 2, 3]}
}"#;

#[test]
fn grid_convergence_reports_the_observed_order_and_band() {
    // A plate heated in the corner farthest from its cold edge (the air
    // above it insulating): pure conduction, so the three grids (4, 8 and
    // 12 cells across) converge monotonically at second order and the
    // fine answer carries a GCI band that contains the extrapolation
    // (measured: 303.4112, 303.4496, 303.4576 K; order 1.80; extrapolated
    // 303.4650 K; GCI 0.0093 K).
    let (code, result, stderr) = run(&scratch("convergence.json", CORNER_PLATE));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        result.str_field("schema"),
        Some("frankensim.cooling-cht.convergence.v1")
    );
    let levels = result.get("levels").and_then(J::as_array).unwrap();
    let cells: Vec<f64> = levels
        .iter()
        .map(|l| l.path(&["cells"]).and_then(J::as_f64).unwrap())
        .collect();
    assert_eq!(cells, [32.0, 256.0, 864.0]);
    let quantities = result.get("quantities").and_then(J::as_array).unwrap();
    let hottest = quantities
        .iter()
        .find(|q| q.str_field("name") == Some("max_solid_temperature_k"))
        .unwrap();
    let values: Vec<f64> = hottest
        .path(&["values"])
        .and_then(J::as_array)
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let order = f(hottest, &["observed_order"]);
    let extrapolated = f(hottest, &["extrapolated"]);
    let gci = f(hottest, &["gci"]);
    eprintln!("values {values:?} order {order} extrapolated {extrapolated} gci {gci}");
    assert_eq!(hottest.str_field("convergence"), Some("monotone"));
    assert!((order - 2.0).abs() < 0.3, "observed order {order}");
    assert!(
        (extrapolated - values[2]).abs() < gci,
        "{extrapolated} {gci}"
    );
    // The fine level's own result rides along, and is the third value.
    assert_eq!(
        f(&result, &["fine_result", "max_solid_temperature_k"]),
        values[2]
    );
    assert!(
        quantities
            .iter()
            .any(|q| q.str_field("name") == Some("source:corner"))
    );
    // Malformed splits and a study alongside refuse.
    for (name, scene) in [
        (
            "splits-order.json",
            CORNER_PLATE.replace("[1, 2, 3]", "[2, 1, 3]"),
        ),
        (
            "splits-count.json",
            CORNER_PLATE.replace("[1, 2, 3]", "[1, 2]"),
        ),
        (
            "with-study.json",
            CORNER_PLATE.replace(
                r#""grid_convergence""#,
                r#""study": {"parameters": [{"name": "p", "path": ["sources", 0, "power_w"], "values": [0.05]}], "objective": {"minimize": "max_solid_temperature_k"}}, "grid_convergence""#,
            ),
        ),
    ] {
        let (code, diagnostic, _) = run(&scratch(name, &scene));
        assert_eq!(code, 4, "{name}");
        assert!(diagnostic.str_field("message").is_some());
    }
}

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
fn graded_grid_refines_the_block_and_closes_energy() {
    // The duct on a graded grid (0.5 mm around the block and near the
    // floor, 1 mm elsewhere): it solves, reports the grid as graded,
    // carries the declared inflow, closes energy, and lands near the
    // uniform 1 mm answer.
    let graded = DUCT.replace(
        r#""size_m": [0.012, 0.004, 0.004], "voxel_m": 0.001,"#,
        r#""grid": {
  "x": [{"to_m": 0.003, "voxel_m": 0.001}, {"to_m": 0.009, "voxel_m": 0.0005}, {"to_m": 0.012, "voxel_m": 0.001}],
  "y": [{"to_m": 0.003, "voxel_m": 0.0005}, {"to_m": 0.004, "voxel_m": 0.001}],
  "z": [{"to_m": 0.004, "voxel_m": 0.001}]},"#,
    );
    assert_ne!(graded, DUCT);
    let (code, uniform, stderr) = run(&scratch("duct-uniform-ref.json", DUCT));
    assert_eq!(code, 0, "{stderr}");
    let (code, result, stderr) = run(&scratch("duct-graded.json", &graded));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.get("graded"), Some(&J::Bool(true)));
    assert_eq!(uniform.get("graded"), Some(&J::Bool(false)));
    assert!((f(&result, &["flow", "inflow_m3_s"]) - 0.1 * 16e-6).abs() < 1e-15);
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    assert!((f(&result, &["energy", "source_w"]) - 0.1).abs() < 1e-12);
    let (a, b) = (
        f(&uniform, &["max_solid_temperature_k"]),
        f(&result, &["max_solid_temperature_k"]),
    );
    assert!(
        (a - b).abs() < 0.25 * (a - 300.0),
        "uniform {a} vs graded {b}"
    );
    // The refine-box form: 1 mm coarse, 0.5 mm in the block's slabs
    // (x 3..9 mm, y 0..3 mm, all of z): 18 x 7 x 8 cells.
    let refined = DUCT.replace(
        r#""voxel_m": 0.001,"#,
        r#""grid": {"voxel_m": 0.001, "refine": [
   {"min_m": [0.003, 0.0, 0.0], "max_m": [0.009, 0.003, 0.004], "voxel_m": 0.0005}]},"#,
    );
    assert_ne!(refined, DUCT);
    let (code, result, stderr) = run(&scratch("duct-refined.json", &refined));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.path(&["cells"]).and_then(J::as_f64), Some(1008.0));
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    // Both spacing declarations at once refuse.
    let both = graded.replace(r#""grid": {"#, r#""voxel_m": 0.001, "grid": {"#);
    let (code, diagnostic, _) = run(&scratch("duct-both.json", &both));
    assert_eq!(code, 4, "{diagnostic:?}");
}

#[test]
fn probes_and_vtk_export_report_the_fields() {
    // A probe in the heated block and one in the inlet air, and the fields
    // written as a VTK rectilinear grid next to the scene.
    let scene = DUCT.replace(
        r#""solver": {"tolerance": 1e-8}"#,
        r#""solver": {"tolerance": 1e-8},
 "probes": [{"name": "chip", "at_m": [0.006, 0.0005, 0.002]},
            {"name": "inlet-air", "at_m": [0.0005, 0.003, 0.002]}],
 "output": {"vtk": "duct-fields.vtr"}"#,
    );
    assert_ne!(scene, DUCT);
    let path = scratch("duct-probes.json", &scene);
    let (code, result, stderr) = run(&path);
    assert_eq!(code, 0, "{stderr}");
    let probes = result.get("probes").and_then(J::as_array).unwrap();
    assert_eq!(probes.len(), 2);
    assert_eq!(probes[0].get("solid"), Some(&J::Bool(true)));
    assert_eq!(probes[1].get("solid"), Some(&J::Bool(false)));
    let chip = probes[0]
        .path(&["temperature_k"])
        .and_then(J::as_f64)
        .unwrap();
    let air = probes[1]
        .path(&["temperature_k"])
        .and_then(J::as_f64)
        .unwrap();
    assert!(chip > air && air >= 300.0 - 1e-9, "{chip} {air}");
    let inflow_speed = probes[1]
        .path(&["velocity_m_s"])
        .and_then(J::as_array)
        .unwrap()[0]
        .as_f64()
        .unwrap();
    assert!(inflow_speed > 0.0);
    let vtk = std::fs::read_to_string(path.with_file_name("duct-fields.vtr")).unwrap();
    assert!(vtk.contains("<RectilinearGrid WholeExtent=\"0 12 0 4 0 4\">"));
    for name in ["temperature_k", "velocity_m_s", "pressure_pa", "material"] {
        assert!(vtk.contains(&format!("Name=\"{name}\"")), "{name}");
    }
    // 192 cells of temperature after its header line.
    let block: Vec<&str> = vtk
        .split("Name=\"temperature_k\" format=\"ascii\">\n")
        .nth(1)
        .unwrap()
        .split("</DataArray>")
        .next()
        .unwrap()
        .lines()
        .collect();
    assert_eq!(block.len(), 192);
}

#[test]
fn falling_conductivity_table_runs_the_chip_hotter() {
    // An aluminium whose conductivity falls steeply with temperature runs
    // hotter than the constant-k part, and the energy still closes.
    let table = DUCT.replace(
        r#""conductivity_w_m_k": 167.0}"#,
        r#""conductivity_w_m_k": 167.0, "conductivity_table": [[300.0, 167.0], [400.0, 0.5]]}"#,
    );
    assert_ne!(table, DUCT);
    let (code, constant, stderr) = run(&scratch("duct-k-constant.json", DUCT));
    assert_eq!(code, 0, "{stderr}");
    let (code, varying, stderr) = run(&scratch("duct-k-table.json", &table));
    assert_eq!(code, 0, "{stderr}");
    let (a, b) = (
        f(&constant, &["max_solid_temperature_k"]),
        f(&varying, &["max_solid_temperature_k"]),
    );
    assert!(b > a, "constant {a} vs k(T) {b}");
    assert!(f(&varying, &["energy", "balance_relative_residual"]) < 1e-9);
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
        (
            DUCT.replace(
                "\"solver\": {\"tolerance\": 1e-8}",
                "\"solver\": {\"tolerance\": 1e-8, \"turbulence\": \"k-epsilon\"}",
            ),
            "solver.turbulence",
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
fn lvel_turbulence_adds_eddy_viscosity_and_cools_the_chip() {
    // The duct at 4 m/s (local wall Reynolds numbers ~ 100; LVEL applies its
    // eddy viscosity wherever those are large, so this exercises the model's
    // plumbing, not a transition claim): the eddy viscosity is reported, the
    // energy still closes, and the added mixing cools the chip.
    let fast = DUCT.replace("[0.1, 0.0, 0.0]", "[4.0, 0.0, 0.0]");
    let lvel = fast.replace(
        "\"solver\": {\"tolerance\": 1e-8}",
        "\"solver\": {\"tolerance\": 1e-8, \"turbulence\": \"lvel\"}",
    );
    assert_ne!(fast, DUCT);
    assert_ne!(lvel, fast);
    let (code, laminar, stderr) = run(&scratch("duct-laminar.json", &fast));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        laminar.path(&["flow", "turbulence"]).and_then(J::as_str),
        Some("laminar")
    );
    let (code, turbulent, stderr) = run(&scratch("duct-lvel.json", &lvel));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        turbulent.path(&["flow", "turbulence"]).and_then(J::as_str),
        Some("lvel")
    );
    let ratio = f(&turbulent, &["flow", "max_eddy_viscosity_ratio"]);
    assert!(ratio > 0.0, "{ratio}");
    assert!(f(&turbulent, &["energy", "balance_relative_residual"]) < 1e-9);
    let (cold, hot) = (
        f(&turbulent, &["max_solid_temperature_k"]),
        f(&laminar, &["max_solid_temperature_k"]),
    );
    assert!(cold < hot && cold > 300.0, "lvel {cold} vs laminar {hot}");
}

const FAN_CHANNEL: &str = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.024, 0.004, 0.004], "voxel_m": 0.001,
 "faces": {
  "x-": {"type": "opening", "ambient_k": 300.0},
  "x+": {"type": "opening", "ambient_k": 300.0},
  "y-": {"type": "symmetry"}, "y+": {"type": "symmetry"},
  "z-": {"type": "symmetry"}, "z+": {"type": "symmetry"}
 },
 "internal_fans": [{"name": "axial", "axis": "x", "at_m": 0.008, "direction": "+",
   "min_m": [0.0, 0.0, 0.0], "max_m": [0.0, 0.004, 0.004],
   "curve": [[0.0, 2.0], [1.6e-5, 0.0]]}],
 "resistances": [{"type": "grille", "axis": "x", "at_m": 0.016,
   "min_m": [0.0, 0.0, 0.0], "max_m": [0.0, 0.004, 0.004], "free_area_ratio": 0.5}],
 "solver": {"tolerance": 1e-9}
}"#;

#[test]
fn internal_fan_and_grille_set_the_channel_operating_point() {
    // Plug flow between two openings: the fan's rise equals the grille's
    // 1/2 rho K U^2 (K from Idelchik's 50 % perforated plate), so the
    // operating point solves a quadratic.
    let (code, result, stderr) = run(&scratch("fan-channel.json", FAN_CHANNEL));
    assert_eq!(code, 0, "{stderr}");
    let fans = result
        .path(&["flow", "internal_fans"])
        .and_then(J::as_array)
        .unwrap();
    assert_eq!(fans[0].str_field("name"), Some("axial"));
    let q = fans[0].path(&["flow_m3_s"]).and_then(J::as_f64).unwrap();
    let rise = fans[0]
        .path(&["pressure_rise_pa"])
        .and_then(J::as_f64)
        .unwrap();
    let (rho, area, p0, q_max) = (1.1614, 16e-6, 2.0, 1.6e-5);
    let k = (1.0 + 0.707 * 0.5f64.sqrt() - 0.5).powi(2) / 0.25;
    let (qa, qb) = (0.5 * rho * k / (area * area), p0 / q_max);
    let exact = (-qb + (qb * qb + 4.0 * qa * p0).sqrt()) / (2.0 * qa);
    assert!((q - exact).abs() < 1e-6 * exact, "{q} vs {exact}");
    assert!((rise - p0 * (1.0 - q / q_max)).abs() < 1e-9, "{rise}");
    assert!((f(&result, &["flow", "inflow_m3_s"]) - q).abs() < 1e-9 * q);
    // A fan plane between voxel faces refuses.
    let off = FAN_CHANNEL.replace("\"at_m\": 0.008", "\"at_m\": 0.0085");
    assert_ne!(off, FAN_CHANNEL);
    let (code, diagnostic, _) = run(&scratch("fan-off-grid.json", &off));
    assert_eq!(code, 4);
    assert!(
        diagnostic
            .str_field("message")
            .unwrap()
            .contains("voxel face"),
        "{diagnostic:?}"
    );
}

#[test]
fn two_resistor_component_on_a_board_reports_its_junction() {
    // A 4 x 1 x 4 mm package (board side y-) on an unheated aluminium strip
    // in the duct: heat leaves through both resistors, they add up to its
    // power, and the energy closes through the binary.
    let scene = DUCT
        .replace(
            r#""max_m": [0.008, 0.002, 0.004]}],"#,
            r#""max_m": [0.008, 0.001, 0.004]}],
 "components": [{"name": "u1", "min_m": [0.004, 0.001, 0.0], "max_m": [0.008, 0.002, 0.004],
   "board_side": "y-", "power_w": 0.05, "junction_to_case_k_w": 40.0, "junction_to_board_k_w": 15.0}],"#,
        )
        .replace(
            r#" "sources": [{"name": "chip", "power_w": 0.1, "min_m": [0.004, 0.0, 0.0], "max_m": [0.008, 0.001, 0.004]}],
"#,
            "",
        );
    assert_ne!(scene, DUCT);
    assert!(!scene.contains("\"chip\""));
    let (code, result, stderr) = run(&scratch("component.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    let parts = result.get("components").and_then(J::as_array).unwrap();
    assert_eq!(parts[0].str_field("name"), Some("u1"));
    let tj = parts[0]
        .path(&["junction_temperature_k"])
        .and_then(J::as_f64)
        .unwrap();
    let case = parts[0].path(&["case_w"]).and_then(J::as_f64).unwrap();
    let board = parts[0].path(&["board_w"]).and_then(J::as_f64).unwrap();
    assert!((case + board - 0.05).abs() < 1e-9, "{case} + {board}");
    assert!(case > 0.0 && board > 0.0, "{case} {board}");
    assert!((f(&result, &["energy", "source_w"]) - 0.05).abs() < 1e-12);
    assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
    // The collapsed package cells report the junction temperature.
    assert!(f(&result, &["max_solid_temperature_k"]) >= tj - 1e-9);
    assert!(tj > 300.0, "{tj}");
}

#[test]
fn sealed_box_radiates_from_the_block_to_its_walls() {
    // A heated block inside a sealed box whose 2 mm walls are an emissive
    // solid held at 300 K from outside: no opening, so only surface-to-
    // surface exchange can carry radiation, and it must cool the block.
    let scene = |emissivity: f64| {
        format!(
            r#"{{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.02, 0.02, 0.02], "voxel_m": 0.002,
 "materials": [
  {{"name": "wall", "conductivity_w_m_k": 200.0, "emissivity": 0.9}},
  {{"name": "block", "conductivity_w_m_k": 50.0, "emissivity": {emissivity}}}
 ],
 "solids": [
  {{"material": "wall", "min_m": [0.0, 0.0, 0.0], "max_m": [0.002, 0.02, 0.02]}},
  {{"material": "wall", "min_m": [0.018, 0.0, 0.0], "max_m": [0.02, 0.02, 0.02]}},
  {{"material": "wall", "min_m": [0.0, 0.0, 0.0], "max_m": [0.02, 0.002, 0.02]}},
  {{"material": "wall", "min_m": [0.0, 0.018, 0.0], "max_m": [0.02, 0.02, 0.02]}},
  {{"material": "wall", "min_m": [0.0, 0.0, 0.0], "max_m": [0.02, 0.02, 0.002]}},
  {{"material": "wall", "min_m": [0.0, 0.0, 0.018], "max_m": [0.02, 0.02, 0.02]}},
  {{"material": "block", "min_m": [0.008, 0.008, 0.008], "max_m": [0.012, 0.012, 0.012]}}
 ],
 "sources": [{{"name": "block", "power_w": 0.2, "min_m": [0.008, 0.008, 0.008], "max_m": [0.012, 0.012, 0.012]}}],
 "faces": {{
  "x-": {{"type": "wall", "temperature_k": 300.0}}, "x+": {{"type": "wall", "temperature_k": 300.0}},
  "y-": {{"type": "wall", "temperature_k": 300.0}}, "y+": {{"type": "wall", "temperature_k": 300.0}},
  "z-": {{"type": "wall", "temperature_k": 300.0}}, "z+": {{"type": "wall", "temperature_k": 300.0}}
 }},
 "radiation": {{"rays_per_face": 512}}
}}"#
        )
    };
    let (code, dark, stderr) = run(&scratch("sealed-dark.json", &scene(0.0)));
    assert_eq!(code, 0, "{stderr}");
    let (code, bright, stderr) = run(&scratch("sealed-bright.json", &scene(0.9)));
    assert_eq!(code, 0, "{stderr}");
    let (hot, cool) = (
        f(&dark, &["max_solid_temperature_k"]),
        f(&bright, &["max_solid_temperature_k"]),
    );
    assert!(cool < hot - 1.0, "radiating {cool} vs dark {hot}");
    // Sealed: nothing reaches surroundings; the walls take it all.
    assert!(f(&bright, &["radiation", "radiated_w"]).abs() < 1e-9);
    assert_eq!(
        bright.path(&["radiation", "surface_exchange"]),
        Some(&J::Bool(true))
    );
    assert!(f(&bright, &["energy", "balance_relative_residual"]) < 1e-9);
    assert!((f(&bright, &["energy", "boundary_outflow_w"]) - 0.2).abs() < 1e-6);
}

#[test]
fn unsteady_march_ramps_the_duct_flow_and_closes_every_step() {
    // The duct started from rest with its inlet ramped over 0.2 s: flow and
    // energy march together; each step's energy closes, the chip warms, and
    // the final inflow is the declared one.
    let scene = DUCT
        .replace(
            r#""conductivity_w_m_k": 167.0}"#,
            r#""conductivity_w_m_k": 167.0, "volumetric_heat_capacity_j_m3_k": 2.4e6}"#,
        )
        .replace(
            r#""solver": {"tolerance": 1e-8}"#,
            r#""solver": {"tolerance": 1e-8},
 "transient": {"time_step_s": 0.05, "steps": 12, "flow": "unsteady",
   "inner_tolerance": 1e-8, "inlet_schedule": [[0.0, 0.0], [0.2, 1.0]]}"#,
        );
    assert_ne!(scene, DUCT);
    let (code, result, stderr) = run(&scratch("unsteady-duct.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.str_field("solver"), Some("fv-simplec-unsteady"));
    assert_eq!(
        result.path(&["transient", "scheme"]).and_then(J::as_str),
        Some("bdf2")
    );
    assert!(f(&result, &["transient", "worst_step_closure_j"]) < 1e-9);
    let inflow = f(&result, &["flow", "inflow_m3_s"]);
    assert!((inflow - 0.1 * 16e-6).abs() < 1e-12, "{inflow}");
    assert!(f(&result, &["flow", "final_kinetic_energy_j"]) > 0.0);
    let peak = f(&result, &["transient", "peak_solid_temperature_k"]);
    assert!(peak > 300.0, "{peak}");
    let records = result
        .path(&["transient", "records"])
        .and_then(J::as_array)
        .unwrap();
    assert_eq!(records.len(), 12);
    // The kinetic energy grows with the ramp and then holds.
    let ke = |i: usize| {
        records[i]
            .path(&["kinetic_energy_j"])
            .and_then(J::as_f64)
            .unwrap()
    };
    assert!(ke(0) < ke(3) && (ke(11) - ke(10)).abs() < 1e-6 * ke(11));
}

#[test]
fn steady_energy_on_the_mean_flow_matches_a_settled_flow() {
    // A flow that settles: the time average over 1..2 s of the started duct
    // is its steady flow, so the steady energy on that mean flow must give
    // the steady scene's chip temperature.
    let (code, steady, stderr) = run(&scratch("duct-steady-ref.json", DUCT));
    assert_eq!(code, 0, "{stderr}");
    let scene = DUCT.replace(
        r#""solver": {"tolerance": 1e-8}"#,
        r#""solver": {"tolerance": 1e-8},
 "transient": {"time_step_s": 0.1, "steps": 20, "flow": "unsteady",
   "energy": "steady-on-mean-flow", "inner_tolerance": 1e-9}"#,
    );
    assert_ne!(scene, DUCT);
    let (code, mean, stderr) = run(&scratch("duct-mean-flow.json", &scene));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        mean.str_field("solver"),
        Some("fv-simplec-unsteady-mean-flow")
    );
    assert_eq!(
        mean.path(&["flow", "averaged_steps"]).and_then(J::as_f64),
        Some(10.0)
    );
    assert!(f(&mean, &["energy", "balance_relative_residual"]) < 1e-9);
    let (a, b) = (
        f(&steady, &["max_solid_temperature_k"]),
        f(&mean, &["max_solid_temperature_k"]),
    );
    assert!(
        (a - b).abs() < 1e-4 * (a - 300.0),
        "steady {a} vs mean-flow {b}"
    );
}

#[test]
fn unsteady_buoyant_column_starts_its_own_draft() {
    // The heated-wall chimney from rest: buoyancy starts the draft without
    // any inlet, and the march reports it.
    let scene = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.004, 0.004, 0.016], "voxel_m": 0.001,
 "faces": {
  "x-": {"type": "wall", "temperature_k": 320.0},
  "z-": {"type": "opening", "ambient_k": 300.0},
  "z+": {"type": "opening", "ambient_k": 300.0}
 },
 "gravity_m_s2": [0.0, 0.0, -9.81],
 "transient": {"time_step_s": 0.02, "steps": 10, "flow": "unsteady", "inner_tolerance": 1e-8}
}"#;
    let (code, result, stderr) = run(&scratch("unsteady-chimney.json", scene));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        result.str_field("solver"),
        Some("fv-simplec-unsteady-boussinesq")
    );
    assert!(f(&result, &["flow", "outflow_m3_s"]) > 0.0);
    assert!(f(&result, &["flow", "final_kinetic_energy_j"]) > 0.0);
    assert!(f(&result, &["transient", "worst_step_closure_j"]) < 1e-9);
}

const FIN_STUDY: &str = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.016, 0.008, 0.006], "voxel_m": 0.001,
 "materials": [{"name": "aluminium", "conductivity_w_m_k": 167.0}],
 "solids": [{"material": "aluminium", "heatsink": {
   "base_min_m": [0.003, 0.0, 0.0], "base_size_m": [0.01, 0.008, 0.001],
   "fin_count": 2, "fin_thickness_m": 0.001, "fin_height_m": 0.003, "fins_along": "x"}}],
 "sources": [{"name": "chip", "power_w": 0.2, "min_m": [0.006, 0.0, 0.0], "max_m": [0.01, 0.008, 0.001]}],
 "faces": {
  "x-": {"type": "inlet", "velocity_m_s": [0.5, 0.0, 0.0], "temperature_k": 300.0},
  "x+": {"type": "opening", "ambient_k": 300.0}
 },
 "solver": {"tolerance": 1e-7},
 "study": {
  "parameters": [{"name": "fins", "path": ["solids", 0, "heatsink", "fin_count"], "values": [2, 3, 5]}],
  "objective": {"minimize": "source:chip"}
 }
}"#;

#[test]
fn fin_count_study_ranks_variants_and_records_refusals() {
    // Three fin counts on a 1 mm grid: 2 and 3 fins solve (more area runs
    // cooler at a fixed inflow); 5 fins merge on this grid and refuse.
    let (code, result, stderr) = run(&scratch("fin-study.json", FIN_STUDY));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(
        result.str_field("schema"),
        Some("frankensim.cooling-cht.study.v1")
    );
    let evaluations = result.get("evaluations").and_then(J::as_array).unwrap();
    assert_eq!(evaluations.len(), 3);
    let objective = |i: usize| evaluations[i].path(&["objective"]).and_then(J::as_f64);
    let (two, three) = (objective(0).unwrap(), objective(1).unwrap());
    assert!(three < two && two > 300.0, "{two} {three}");
    assert_eq!(evaluations[2].str_field("status"), Some("refused"));
    assert!(
        evaluations[2]
            .str_field("message")
            .unwrap()
            .contains("fins merge"),
        "{:?}",
        evaluations[2]
    );
    assert_eq!(
        result.path(&["best", "index"]).and_then(J::as_f64),
        Some(1.0)
    );
    // One thread or several: the same report, bit for bit.
    let serial = FIN_STUDY.replace(
        r#""objective": {"minimize": "source:chip"}"#,
        r#""objective": {"minimize": "source:chip"}, "parallelism": 1"#,
    );
    let (code, one, stderr) = run(&scratch("fin-study-serial.json", &serial));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(one.path(&["threads"]).and_then(J::as_f64), Some(1.0));
    assert_eq!(one.get("evaluations"), result.get("evaluations"));
    // An unreachable constraint leaves no feasible variant.
    let constrained = FIN_STUDY.replace(
        r#""objective": {"minimize": "source:chip"}"#,
        r#""objective": {"minimize": "source:chip"},
  "constraints": [{"quantity": "inflow_m3_s", "min": 1.0}]"#,
    );
    assert_ne!(constrained, FIN_STUDY);
    let (code, result, stderr) = run(&scratch("fin-study-constrained.json", &constrained));
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(result.get("best"), Some(&J::Null));
    // A path that does not exist refuses the study.
    let wrong = FIN_STUDY.replace(r#""heatsink", "fin_count""#, r#""heatsink", "fins""#);
    let (code, diagnostic, _) = run(&scratch("fin-study-wrong.json", &wrong));
    assert_eq!(code, 4, "{diagnostic:?}");
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

#[test]
fn emissive_block_in_a_vented_column_runs_cooler_and_closes_with_radiation() {
    // A heated block in a column open at both ends, cooled by natural
    // convection: declaring an emissivity adds radiation to the
    // surroundings through the openings, which can only cool it, and the
    // radiated power is the energy balance's sink line.
    let scene = r#"{
 "schema": "frankensim.cooling-cht.v1",
 "size_m": [0.006, 0.006, 0.016], "voxel_m": 0.001,
 "materials": [{"name": "aluminium", "conductivity_w_m_k": 167.0, "emissivity": EPS}],
 "solids": [{"material": "aluminium", "min_m": [0.002, 0.002, 0.006], "max_m": [0.004, 0.004, 0.010]}],
 "sources": [{"name": "block", "power_w": 0.05, "min_m": [0.002, 0.002, 0.006], "max_m": [0.004, 0.004, 0.010]}],
 "faces": {
  "z-": {"type": "opening", "ambient_k": 300.0},
  "z+": {"type": "opening", "ambient_k": 300.0}
 },
 "gravity_m_s2": [0.0, 0.0, -9.81],
 "radiation": {"rays_per_face": 512},
 "solver": {"tolerance": 1e-7}
}"#;
    let solve = |eps: &str, name: &str| {
        let (code, result, stderr) = run(&scratch(name, &scene.replace("EPS", eps)));
        assert_eq!(code, 0, "{stderr}");
        assert!(f(&result, &["energy", "balance_relative_residual"]) < 1e-9);
        result
    };
    let dark = solve("0.0", "column-dark.json");
    let grey = solve("0.9", "column-grey.json");
    assert!(dark.get("radiation").is_none());
    let (t_dark, t_grey) = (
        f(&dark, &["max_solid_temperature_k"]),
        f(&grey, &["max_solid_temperature_k"]),
    );
    assert!(
        t_grey < t_dark - 0.1,
        "radiation must cool: {t_dark} -> {t_grey}"
    );
    let radiated = f(&grey, &["radiation", "radiated_w"]);
    let sink = f(&grey, &["energy", "sink_outflow_w"]);
    assert!(radiated > 0.0 && radiated < 0.05, "{radiated}");
    assert!(
        (radiated - sink).abs() < 1e-6 * 0.05,
        "{radiated} vs {sink}"
    );
    // Everything generated leaves by advection, conduction or radiation.
    let leaving = f(&grey, &["energy", "boundary_outflow_w"]) + sink;
    assert!((leaving - 0.05).abs() < 1e-9, "{leaving}");
    // The buoyant march from rest carries the same radiation: every step
    // closes with the radiated power inside its outflow.
    let marched = scene
        .replace(
            r#""emissivity": EPS}"#,
            r#""emissivity": 0.9, "volumetric_heat_capacity_j_m3_k": 2.4e6}"#,
        )
        .replace(
            r#""solver": {"tolerance": 1e-7}"#,
            r#""solver": {"tolerance": 1e-7}, "transient": {"time_step_s": 0.5, "steps": 4, "flow": "unsteady", "inner_tolerance": 1e-8}"#,
        );
    let (code, result, stderr) = run(&scratch("rad-unsteady.json", &marched));
    assert_eq!(code, 0, "{stderr}");
    assert!(f(&result, &["transient", "worst_step_closure_j"]) < 1e-9);
    let records = result
        .path(&["transient", "records"])
        .and_then(J::as_array)
        .unwrap();
    let radiated: Vec<f64> = records
        .iter()
        .map(|r| r.path(&["radiated_w"]).and_then(J::as_f64).unwrap())
        .collect();
    // The block warms from 300 K, so it radiates more every step.
    assert!(
        radiated[0] > 0.0 && radiated.windows(2).all(|w| w[1] > w[0]),
        "{radiated:?}"
    );
}
