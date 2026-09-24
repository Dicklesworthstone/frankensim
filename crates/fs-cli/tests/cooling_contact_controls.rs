//! Actual-binary tests for contact controls. Closed-form slab and independent
//! dense nonlinear FEM values are reference calculations, not prior Rust runs.
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static SLAB: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/size-contact-slab.json"))));
// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static PULSE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonlinear-contact-pulse.json"))));

fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-contact-controls-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name); fs::write(&path, text).unwrap(); path
}
fn command(kind: &str, path: &Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    cmd.args(["--json", kind]).arg(path); cmd
}
fn document(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn direct(dir: &Path, name: &str, text: &str) -> J {
    let path = write(dir, name, text);
    document(&command("cooling-network", &path).output().unwrap())
}
fn number(doc: &J, path: &[&str]) -> f64 { doc.path(path).and_then(J::as_f64).unwrap() }
fn close(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a} != {b}"); }
fn without_last_section(text: &str, section: &str) -> String {
    let delimiter = format!(",\"{section}\":");
    let (prefix, _) = text.split_once(&delimiter).expect("fixture's final section");
    format!("{prefix}}}")
}
fn slab(resistance: f64, gradient: bool) -> String {
    without_last_section(SLAB, "design")
        .replace("\"gradient\":true", &format!("\"gradient\":{gradient}"))
        .replace("\"resistance_m2_k_w\":0.01", &format!("\"resistance_m2_k_w\":{resistance}"))
}
fn nonlinear_steady(resistance: f64, gradient: bool) -> String {
    without_last_section(PULSE, "transient")
        .replace("\"gradient\":false", &format!("\"gradient\":{gradient}"))
        .replace("\"resistance_m2_k_w\":0.01", &format!("\"resistance_m2_k_w\":{resistance}"))
}
fn sensitivity_rows(doc: &J) -> &[J] {
    doc.path(&["contact_sensitivities", "rows"]).unwrap().as_array().unwrap()
}
fn gradient(doc: &J) -> f64 {
    sensitivity_rows(doc)[0].f64_field("dobjective_dlog_resistance_k").unwrap()
}
fn uq_plan(lo: f64, hi: f64, samples: usize, transient: bool) -> String {
    let qoi = if transient { ",\"qoi\":{\"kind\":\"transient-sampled-peak\"}" } else { "" };
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{samples},"wall_seconds":300,"temperature_limit_k":305,"correlation":{{"kind":"independent"}}{qoi},"parameters":[{{"target":{{"kind":"contact-resistance","contact":"bondline"}},"distribution":{{"kind":"uniform","lo":{lo},"hi":{hi}}}}}]}}"#)
}

#[test]
fn contact_gradient_matches_closed_form_slab_and_perturbed_coupled_solves() {
    let dir = scratch("slab");
    let nominal = direct(&dir, "nominal.json", &slab(0.01, true));
    let c1 = 1.2 * 0.003 * 1007.0_f64;
    let ct = 1.2 * 0.004 * 1007.0_f64;
    let f1 = c1 * -(-0.5 / c1).exp_m1();
    let f2 = ct * -(-0.8 / ct).exp_m1();
    let resistance = 0.5 + 1.0 + 1.0 / f1 + 1.0 / f2 - 1.0 / ct;
    let expected = 10.0 / (f1 * resistance * resistance);
    close(gradient(&nominal), expected, 2e-6);
    close(number(&nominal, &["objective", "value_k"]), 330.0 - 10.0 / resistance / f1, 1e-5);
    close(sensitivity_rows(&nominal)[0].f64_field("dobjective_dresistance_w_m2").unwrap(), expected / 0.01, 2e-4);
    let delta = 1e-4_f64;
    let plus = direct(&dir, "plus.json", &slab(0.01 * delta.exp(), false));
    let minus = direct(&dir, "minus.json", &slab(0.01 * (-delta).exp(), false));
    close(gradient(&nominal), (number(&plus, &["objective", "value_k"])
        - number(&minus, &["objective", "value_k"])) / (2.0 * delta), 5e-5);
    assert_eq!(plus.get("contact_sensitivities"), Some(&J::Null));
    let reversed = slab(0.01, true).replace("\"side_a\":", "\"temporary_side\":")
        .replace("\"side_b\":", "\"side_a\":").replace("\"temporary_side\":", "\"side_b\":");
    let reversed = direct(&dir, "reversed.json", &reversed);
    close(gradient(&nominal), gradient(&reversed), 1e-7);
}

#[test]
fn nonlinear_contact_gradient_includes_the_temperature_dependent_material_jacobian() {
    let dir = scratch("nonlinear");
    let nominal = direct(&dir, "nominal.json", &nonlinear_steady(0.01, true));
    // Independent 16-node P1 FEM with air-network elimination and K'(T).
    close(number(&nominal, &["objective", "value_k"]), 332.4354972844827, 3e-5);
    close(gradient(&nominal), 3.0602605735196953, 3e-5);
    let delta = 1e-4_f64;
    let plus = direct(&dir, "plus.json", &nonlinear_steady(0.01 * delta.exp(), false));
    let minus = direct(&dir, "minus.json", &nonlinear_steady(0.01 * (-delta).exp(), false));
    close(gradient(&nominal), (number(&plus, &["objective", "value_k"])
        - number(&minus, &["objective", "value_k"])) / (2.0 * delta), 5e-4);
    close(number(&nominal, &["source_w"]), 20.0, 1e-8);
    close(number(&nominal, &["robin_out_w"]), 20.0, 1e-7);
}

#[test]
fn splitting_the_interface_preserves_total_derivative_and_orders_controls_by_name() {
    let dir = scratch("partition");
    let text = slab(0.01, true);
    let start = text.find("\"contacts\":[").unwrap();
    let end = text.find("},\"objective\"").unwrap();
    let replacement = r#""contacts":[
        {"name":"z-first","source":"fixture","side_a_material":"left","side_b_material":"right","resistance_m2_k_w":0.01,
         "face_pairs":[{"side_a":[1,4,10],"side_b":[12,13,15]}]},
        {"name":"a-second","source":"fixture","side_a_material":"left","side_b_material":"right","resistance_m2_k_w":0.01,
         "face_pairs":[{"side_a":[1,7,10],"side_b":[12,14,15]}]}]"#;
    let split = format!("{}{replacement}{}", &text[..start], &text[end..]);
    let original = direct(&dir, "whole.json", &text);
    let split = direct(&dir, "split.json", &split);
    let rows = sensitivity_rows(&split);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].str_field("contact"), Some("a-second"));
    assert_eq!(rows[1].str_field("contact"), Some("z-first"));
    close(rows.iter().map(|r| r.f64_field("dobjective_dlog_resistance_k").unwrap()).sum(), gradient(&original), 1e-7);
    close(number(&split, &["objective", "value_k"]), number(&original, &["objective", "value_k"]), 1e-7);
}

#[test]
fn fixed_contact_uncertainty_matches_real_steady_and_transient_observables() {
    let dir = scratch("degenerate");
    for (index, transient) in [false, true].into_iter().enumerate() {
        let text = if transient { PULSE.to_string() } else { slab(0.01, false) };
        let base = write(&dir, &format!("base-{index}.json"), &text);
        let nominal = document(&command("cooling-network", &base).output().unwrap());
        let path: &[&str] = if transient { &["transient", "sampled_peak_objective_k"] } else { &["objective", "value_k"] };
        let expected = number(&nominal, path);
        let uq = write(&dir, &format!("uq-{index}.json"), &uq_plan(0.01, 0.01, 2, transient));
        let report = document(&command("cooling-network-uq", &base).arg(&uq).output().unwrap());
        close(number(&report, &["mean_k"]), expected, 1e-8);
        close(number(&report, &["std_dev_k"]), 0.0, 1e-10);
        assert_eq!(number(&report, &["samples_evaluated"]), 2.0);
        assert_eq!(nominal.get("contact_sensitivities"), Some(&J::Null));
        if transient {
            close(expected, 306.1651033635, 3e-5);
            assert_eq!(number(&report, &["empirical_probability_of_compliance"]), 0.0);
        }
    }
}

#[test]
fn uncertain_contact_trajectories_resume_with_identical_results_and_checkpoint_bytes() {
    let dir = scratch("replay");
    let text = PULSE.replace("\"duration_s\":30", "\"duration_s\":4")
        .replace("\"duration_s\":120", "\"duration_s\":4");
    let base = write(&dir, "base.json", &text);
    let uq = write(&dir, "uq.json", &uq_plan(0.005, 0.02, 4, true));
    let full_path = dir.join("full.bin");
    let full = command("cooling-network-uq", &base).arg(&uq)
        .arg("--checkpoint").arg(&full_path).output().unwrap();
    assert!(number(&document(&full), &["std_dev_k"]) > 0.0);
    let prefix = dir.join("prefix.bin");
    let first = command("cooling-network-uq", &base).arg(&uq)
        .arg("--checkpoint").arg(&prefix).args(["--max-new-samples", "2"]).output().unwrap();
    assert_eq!(first.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    let final_path = dir.join("final.bin");
    let resumed = command("cooling-network-uq", &base).arg(&uq)
        .arg("--resume").arg(&prefix).arg("--checkpoint").arg(&final_path).output().unwrap();
    document(&resumed);
    assert_eq!(full.stdout, resumed.stdout);
    assert_eq!(fs::read(full_path).unwrap(), fs::read(final_path).unwrap());
}

#[test]
fn invalid_contacts_and_nonphysical_contact_samples_cannot_publish_statistics() {
    let dir = scratch("refused");
    let base = write(&dir, "base.json", &slab(0.01, false));
    let plan = uq_plan(0.005, 0.02, 2, false);
    let missing = write(&dir, "missing.json", &plan.replace("\"contact\":\"bondline\"", "\"contact\":\"missing\""));
    let destination = dir.join("must-not-exist.bin");
    let output = command("cooling-network-uq", &base).arg(missing)
        .arg("--checkpoint").arg(&destination).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(!destination.exists());
    let invalid = write(&dir, "invalid.json", &plan.replace(
        "\"kind\":\"uniform\",\"lo\":0.005,\"hi\":0.02",
        "\"kind\":\"gaussian\",\"mean\":-1,\"std_dev\":0"));
    let failed_path = dir.join("failed.bin");
    let output = command("cooling-network-uq", &base).arg(invalid)
        .arg("--checkpoint").arg(&failed_path).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(fs::read(failed_path).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
}
