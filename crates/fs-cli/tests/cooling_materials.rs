//! Actual product-boundary regressions. No surrogate thermal evaluator is used.
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;

use json::JsonValue as J;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}

fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-materials-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}

fn command(kind: &str, base: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.arg("--json").arg(kind).arg(base);
    command
}

fn document(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) { assert!((a - b).abs() <= tolerance, "{a} != {b}"); }

fn fixed_request(path: &Path, axial_k: f64) {
    fs::write(path, format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":300,"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"material-principal-conductivity","material":"substrate","axis":0}},"distribution":{{"kind":"uniform","lo":{axial_k},"hi":{axial_k}}}}}]}}"#)).unwrap();
}

#[test]
fn oriented_principal_values_and_the_same_full_tensor_produce_the_same_field() {
    let base = example("orthotropic-hotspot.json");
    let dir = scratch("tensor");
    let text = json::compact(&fs::read_to_string(&base).unwrap());
    let from = r#""orthotropic":{"principal_axes":[[0.6,0.8,0],[-0.8,0.6,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}"#;
    assert!(text.contains(from));
    let tensor = dir.join("tensor.json");
    fs::write(&tensor, text.replace(from, r#""conductivity_tensor_w_m_k":[[8.48,8.64,0],[8.64,13.52,0],[0,0,1]]"#)).unwrap();
    let a = document(&command("cooling-network", &base).output().unwrap());
    let b = document(&command("cooling-network", &tensor).output().unwrap());
    let a_field = a.get("solid_temperatures_k").unwrap().as_array().unwrap();
    let b_field = b.get("solid_temperatures_k").unwrap().as_array().unwrap();
    assert_eq!(a_field.len(), b_field.len());
    for (a, b) in a_field.iter().zip(b_field) { close(a.as_f64().unwrap(), b.as_f64().unwrap(), 1e-6); }
    close(a.f64_field("source_w").unwrap(), 1.0, 1e-9);
    close(a.f64_field("robin_out_w").unwrap(), 1.0, 1e-7);
    let materials = a.path(&["solid_inputs", "constitutive", "materials"]).unwrap().as_array().unwrap();
    let substrate = materials.iter().find(|row| row.str_field("name") == Some("substrate")).unwrap();
    assert_eq!(substrate.str_field("coordinate_frame"), Some("mesh-cartesian"));
    assert!(substrate.get("resolved_conductivity_tensor_w_m_k").is_some());
}

#[test]
fn zero_material_uncertainty_matches_the_nominal_hotspot_and_changed_k_changes_it() {
    let base = example("orthotropic-hotspot.json");
    let dir = scratch("degenerate");
    let nominal = document(&command("cooling-network", &base).output().unwrap());
    let qoi = nominal.path(&["objective", "value_k"]).and_then(J::as_f64).unwrap();
    let request = dir.join("fixed-k.json");
    fixed_request(&request, 20.0);
    let full = document(&command("cooling-network-uq", &base).arg(&request).output().unwrap());
    close(full.f64_field("mean_k").unwrap(), qoi, 1e-8);
    close(full.f64_field("std_dev_k").unwrap(), 0.0, 1e-10);
    fixed_request(&request, 2.0);
    let changed = document(&command("cooling-network-uq", &base).arg(&request).output().unwrap());
    assert!((changed.f64_field("mean_k").unwrap() - qoi).abs() > 1e-5,
        "changing an active material coefficient must change the real thermal answer");
}

#[test]
fn directional_material_uncertainty_resumes_to_the_identical_real_solve_result() {
    let base = example("orthotropic-hotspot.json");
    let request = example("uq-orthotropic-hotspot.json");
    let dir = scratch("resume");
    let full_path = dir.join("full.bin");
    let full = command("cooling-network-uq", &base).arg(&request)
        .arg("--checkpoint").arg(&full_path).output().unwrap();
    let report = document(&full);
    assert!(report.f64_field("std_dev_k").unwrap() > 0.0);
    let prefix = dir.join("prefix.bin");
    let partial = command("cooling-network-uq", &base).arg(&request)
        .arg("--checkpoint").arg(&prefix).args(["--max-new-samples", "3"]).output().unwrap();
    assert_eq!(partial.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&partial.stderr));
    let final_path = dir.join("final.bin");
    let resumed = command("cooling-network-uq", &base).arg(&request)
        .arg("--resume").arg(&prefix).arg("--checkpoint").arg(&final_path).output().unwrap();
    document(&resumed);
    assert_eq!(full.stdout, resumed.stdout);
    assert_eq!(fs::read(full_path).unwrap(), fs::read(final_path).unwrap());
}

#[test]
fn scalarizing_a_directional_material_refuses_before_creating_a_checkpoint() {
    let dir = scratch("refusal");
    let base = example("orthotropic-hotspot.json");
    let request = dir.join("invalid-target.json");
    fixed_request(&request, 20.0);
    let text = json::compact(&fs::read_to_string(&request).unwrap());
    fs::write(&request, text.replace(
        r#""kind":"material-principal-conductivity","material":"substrate","axis":0"#,
        r#""kind":"material-conductivity","material":"substrate""#,
    )).unwrap();
    let checkpoint = dir.join("must-not-exist.bin");
    let result = command("cooling-network-uq", &base).arg(&request)
        .arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    assert!(!checkpoint.exists());
    assert!(String::from_utf8_lossy(&result.stderr).contains("do not scalarize a tensor"));
}

#[test]
fn nonlinear_material_reaches_the_real_solver_and_refuses_extrapolation() {
    let dir = scratch("nonlinear");
    let text = json::compact(&fs::read_to_string(example("orthotropic-hotspot.json")).unwrap());
    let scalar = r#""name":"spreader","conductivity_w_m_k":20"#;
    assert!(text.contains(scalar));
    let curve = dir.join("curve.json");
    fs::write(&curve, text.replace(scalar,
        r#""name":"spreader","conductivity_curve":{"temperature_k":[250,400],"conductivity_w_m_k":[17,2]}"#,
    )).unwrap();
    let nonlinear = document(&command("cooling-network", &curve).output().unwrap());
    close(nonlinear.f64_field("source_w").unwrap(), 1.0, 1e-9);
    close(nonlinear.f64_field("robin_out_w").unwrap(), 1.0, 1e-7);
    let materials = nonlinear.path(&["solid_inputs", "constitutive", "materials"]).unwrap().as_array().unwrap();
    let spreader = materials.iter().find(|row| row.str_field("name") == Some("spreader")).unwrap();
    assert_eq!(spreader.str_field("temperature_extrapolation"), Some("refused"));
    let frozen = dir.join("frozen.json");
    fs::write(&frozen, text.replace(scalar, r#""name":"spreader","conductivity_w_m_k":17"#)).unwrap();
    let constant = document(&command("cooling-network", &frozen).output().unwrap());
    let qoi = |doc: &J| doc.path(&["objective", "value_k"]).and_then(J::as_f64).unwrap();
    assert!((qoi(&nonlinear) - qoi(&constant)).abs() > 1e-5,
        "the inactive first-knot scalar must not replace the temperature law");
    let outside = dir.join("outside-span.json");
    fs::write(&outside, text.replace(scalar,
        r#""name":"spreader","conductivity_curve":{"temperature_k":[250,260],"conductivity_w_m_k":[17,16]}"#,
    )).unwrap();
    let refused = command("cooling-network", &outside).output().unwrap();
    assert!(!refused.status.success());
    assert!(refused.stdout.is_empty());
}
