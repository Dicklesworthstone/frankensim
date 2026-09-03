//! Integration tests for C1 STEP import path parity with STL (bead `frankensim-rc-root-q61wp.51`).
//!
//! Verifies:
//! 1. Full 7-stage solve and QoI parity between STL and STEP heatsink bodies.
//! 2. Import receipt structure and STEP entity admission naming.
//! 3. Refusal of non-manifold STEP with missing face (hole detection).
//! 4. Refusal of unit mismatch at conduction with `cli-solve-conduction-length-unit`.
//! 5. Provenance audit and face count matching.

use fs_cli::{exit, run};
use fs_io::{NamedFaceGroup, STEP_FACETED_DECODER_VERSION, parse_step};
use fs_ledger::Ledger;
use std::path::{Path, PathBuf};

fn args(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fs-cli-step-heatsink-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

fn extract_temperature_max(json_str: &str) -> Option<f64> {
    // Parse the JSON output for "temperature-max" or find "value":<float>
    for line in json_str.lines() {
        if line.contains("\"name\":\"temperature-max\"") || line.contains("\"temperature-max\"") {
            if let Some(pos) = line.find("\"value\":") {
                let rest = &line[pos + 8..];
                let num_str: String = rest
                    .chars()
                    .take_while(|c| c.is_ascii_digit() || *c == '.' || *c == '-' || *c == 'e' || *c == 'E')
                    .collect();
                if let Ok(v) = num_str.parse::<f64>() {
                    return Some(v);
                }
            }
        }
    }
    // Fallback: search anywhere for `"value":` near temperature
    if let Some(pos) = json_str.find("\"value\":") {
        let rest = &json_str[pos + 8..];
        let num_str: String = rest
            .chars()
            .take_while(|c| c.is_ascii_digit() || *c == '.' || *c == '-' || *c == 'e' || *c == 'E' || *c == '+')
            .collect();
        if let Ok(v) = num_str.parse::<f64>() {
            return Some(v);
        }
    }
    None
}

#[test]
fn step_heatsink_001_qoi_parity_with_stl_across_seven_stages() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim_stl = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let fsim_step = root.join("examples/heatsink-fan/heatsink-fan-step.fsim");
    let step = root.join("examples/heatsink-fan/heatsink.step");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");

    // 1. Validate STEP project
    let validated = run(args(&[
        "--json",
        "validate",
        fsim_step.to_string_lossy().as_ref(),
    ]));
    assert_eq!(validated.exit_code, exit::SUCCESS, "validate stderr: {}", validated.stderr);
    assert!(validated.stdout.contains("\"status\":\"ok\""));
    assert!(validated.stdout.contains("\"finding_count\":0"));

    // 2. Run STL pipeline (baseline)
    let dir_stl = scratch("stl-baseline");
    let ledger_stl = dir_stl.join("heatsink-stl.db");
    let import_stl = run(args(&[
        "--json",
        "import",
        fsim_stl.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger_stl.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(import_stl.exit_code, exit::SUCCESS, "stl import stderr: {}", import_stl.stderr);

    let run_stl = run(args(&[
        "--json",
        "run",
        fsim_stl.to_string_lossy().as_ref(),
        ledger_stl.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(run_stl.exit_code, exit::SUCCESS, "stl run stdout: {} / stderr: {}", run_stl.stdout, run_stl.stderr);
    assert!(run_stl.stdout.contains("\"stages_completed\":7"));
    assert!(run_stl.stdout.contains("\"status\":\"completed\""));
    let qoi_stl = extract_temperature_max(&run_stl.stdout).expect("temperature-max in STL run output");

    // 3. Run STEP pipeline
    let dir_step = scratch("step-run");
    let ledger_step = dir_step.join("heatsink-step.db");
    let import_step = run(args(&[
        "--json",
        "import",
        fsim_step.to_string_lossy().as_ref(),
        step.to_string_lossy().as_ref(),
        ledger_step.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--step-root",
        "5001",
        "--target-h",
        "0.005",
    ]));
    assert_eq!(import_step.exit_code, exit::SUCCESS, "step import stderr: {}", import_step.stderr);

    // Verify import receipt in STEP ledger
    let ledger = Ledger::open(ledger_step.to_str().unwrap()).expect("open step ledger");
    assert!(ledger.lint().expect("lint").is_clean());

    let run_step = run(args(&[
        "--json",
        "run",
        fsim_step.to_string_lossy().as_ref(),
        ledger_step.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(run_step.exit_code, exit::SUCCESS, "step run stdout: {} / stderr: {}", run_step.stdout, run_step.stderr);
    assert!(run_step.stdout.contains("\"stages_completed\":7"));
    assert!(run_step.stdout.contains("\"status\":\"completed\""));
    let qoi_step = extract_temperature_max(&run_step.stdout).expect("temperature-max in STEP run output");

    // 4. Verify QoI parity: within discretization term
    let rel_diff = ((qoi_step - qoi_stl) / qoi_stl).abs();
    assert!(
        rel_diff < 1e-4,
        "QoI parity failure: STL qoi={qoi_stl}, STEP qoi={qoi_step}, rel_diff={rel_diff:.6e} exceeds 1e-4"
    );
}

#[test]
fn step_heatsink_002_import_receipt_records_step_entities_and_mesh_counts() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim_step = root.join("examples/heatsink-fan/heatsink-fan-step.fsim");
    let step = root.join("examples/heatsink-fan/heatsink.step");

    let dir = scratch("step-receipt-test");
    let ledger_path = dir.join("receipt.db");
    let imported = run(args(&[
        "--json",
        "import",
        fsim_step.to_string_lossy().as_ref(),
        step.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--step-root",
        "5001",
        "--target-h",
        "0.005",
    ]));
    assert_eq!(imported.exit_code, exit::SUCCESS, "import failed: {}", imported.stderr);

    let ledger = Ledger::open(ledger_path.to_str().unwrap()).expect("open ledger");
    let all_artifacts = ledger.all_artifacts().expect("read artifacts");
    let mut found_receipt = false;
    for (_hash, bytes) in all_artifacts {
        if let Ok(text) = String::from_utf8(bytes) {
            if text.contains("frankensim.cli.faceted-step-import-receipt.v1") {
                found_receipt = true;
                assert!(text.contains("\"kind\":\"step-triangular-faceted-brep-receipt\""));
                assert!(text.contains("\"root_id\":5001"));
                assert!(text.contains("\"repaired_mesh\":{\"vertices\":56,\"triangles\":108}"));
                assert!(text.contains("CONFIG_CONTROL_DESIGN"));
            }
        }
    }
    assert!(found_receipt, "faceted-step-import-receipt not found in ledger artifacts");
}

#[test]
fn step_heatsink_003_missing_face_refuses_at_import_with_defect_diagnostics() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let step_bytes = std::fs::read_to_string(root.join("examples/heatsink-fan/heatsink.step"))
        .expect("read heatsink.step");

    // Remove one face (#3000) from the CLOSED_SHELL list (#5000)
    // This creates an open boundary hole in what was a closed manifold.
    let bad_step = step_bytes.replace("#3000,", "");
    let dir = scratch("missing-face");
    let bad_step_path = dir.join("heatsink-open.step");
    std::fs::write(&bad_step_path, &bad_step).expect("write bad step");

    let parsed = parse_step(bad_step.as_bytes()).expect("syntax is still valid Part-21");
    let fp = parsed.receipt().source_fingerprint();

    // Create a project matching the new source fingerprint
    let fsim_content = std::fs::read_to_string(root.join("examples/heatsink-fan/heatsink-fan-step.fsim"))
        .expect("read fsim");
    let bad_fsim_content = fsim_content.replace(
        "2ad6dfc2e3e7ff92",
        &format!("{:016x}", fp),
    );
    let bad_fsim_path = dir.join("heatsink-open.fsim");
    std::fs::write(&bad_fsim_path, bad_fsim_content).expect("write bad fsim");

    let ledger_path = dir.join("open.db");
    let imported = run(args(&[
        "--json",
        "import",
        bad_fsim_path.to_string_lossy().as_ref(),
        bad_step_path.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--step-root",
        "5001",
        "--target-h",
        "0.005",
    ]));

    // Must refuse at import because of residual boundary defects (open shell / hole)
    assert_eq!(
        imported.exit_code,
        exit::REFUSED,
        "importing non-manifold / open STEP must refuse; stdout={}, stderr={}",
        imported.stdout,
        imported.stderr
    );
    assert!(
        imported.stderr.contains("cli-import-step-tessellation") || imported.stderr.contains("boundary"),
        "stderr should name the tessellation defect / boundary edges: {}",
        imported.stderr
    );
}

#[test]
fn step_heatsink_004_millimetre_unit_refuses_at_conduction() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim_step = root.join("examples/heatsink-fan/heatsink-fan-step.fsim");
    let step = root.join("examples/heatsink-fan/heatsink.step");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");

    // Modify the assignment in the project to declare length-unit "mm"
    let fsim_content = std::fs::read_to_string(&fsim_step).expect("read fsim");
    let mm_fsim_content = fsim_content.replace(
        ":length-unit \"m\"",
        ":length-unit \"mm\"",
    );
    let dir = scratch("unit-mismatch");
    let mm_fsim_path = dir.join("heatsink-mm.fsim");
    std::fs::write(&mm_fsim_path, mm_fsim_content).expect("write mm fsim");

    let ledger_path = dir.join("mm.db");
    // Import with --unit mm
    let imported = run(args(&[
        "--json",
        "import",
        mm_fsim_path.to_string_lossy().as_ref(),
        step.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--unit",
        "mm",
        "--step-root",
        "5001",
        "--target-h",
        "5.0",
    ]));
    assert_eq!(imported.exit_code, exit::SUCCESS, "import with mm succeeds: {}", imported.stderr);

    // Now run solve / conduction
    let run_res = run(args(&[
        "--json",
        "run",
        mm_fsim_path.to_string_lossy().as_ref(),
        ledger_path.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));

    // Must refuse at conduction stage with code cli-solve-conduction-length-unit
    assert_eq!(
        run_res.exit_code,
        exit::REFUSED,
        "conduction with non-metre unit must refuse: stdout={}, stderr={}",
        run_res.stdout,
        run_res.stderr
    );
    assert!(
        run_res.stderr.contains("cli-solve-conduction-length-unit"),
        "stderr must contain cli-solve-conduction-length-unit; got: {}",
        run_res.stderr
    );
    assert!(
        run_res.stderr.contains("metres"),
        "stderr must explain metres requirement; got: {}",
        run_res.stderr
    );
}

#[test]
fn step_heatsink_005_provenance_record_matches_cad_model() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let prov_path = root.join("examples/heatsink-fan/heatsink.step.provenance.json");
    let step_path = root.join("examples/heatsink-fan/heatsink.step");

    let prov_bytes = std::fs::read(&prov_path).expect("read provenance");
    let prov_str = String::from_utf8(prov_bytes).expect("provenance utf8");
    assert!(prov_str.contains("\"face_count\": 108"));
    assert!(prov_str.contains("\"vertex_count\": 56"));
    assert!(prov_str.contains("\"root_id\": 5001"));
    assert!(prov_str.contains("\"source_fingerprint\": \"2ad6dfc2e3e7ff92\""));
    assert!(prov_str.contains("\"parser_version\": \"step-triangular-faceted-brep-v2\""));

    let step_bytes = std::fs::read(&step_path).expect("read step");
    let parsed = parse_step(&step_bytes).expect("parse step");
    assert_eq!(
        format!("{:016x}", parsed.receipt().source_fingerprint()),
        "2ad6dfc2e3e7ff92",
        "step file source fingerprint must match recorded provenance"
    );
}
