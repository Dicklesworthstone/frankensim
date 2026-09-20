//! Actual command: exact fit does not establish identifiable parameters.
use std::{path::PathBuf, process::{Command, Output}};
fn examples() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/equilibrium-uncertainty") }
fn controls() -> Vec<String> {
    "--rank-relative-tolerance 0.00001 --max-observations 16 --max-adjoints 16 --point-x left 0 --point-x right 0"
        .split_whitespace().map(str::to_owned).collect()
}
fn run(kind: &str, controls: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .arg(examples().join(format!("sensitivity-{kind}.model")))
        .arg(examples().join(format!("sensitivity-{kind}.fit")))
        .args(controls).output().unwrap()
}
fn text(output: &Output) -> String {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    String::from_utf8(output.stdout.clone()).unwrap()
}
fn scalar(text: &str, key: &str) -> f64 {
    text.split(&format!("\"{key}\":")).nth(1).unwrap().split([',', '}']).next().unwrap().parse().unwrap()
}
fn temporary(name: &str, content: &str) -> PathBuf {
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-observations-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&dir).unwrap(); let path = dir.join(name); std::fs::write(&path, content).unwrap(); path
}

#[test]
fn exact_fits_distinguish_independent_stiffnesses_from_a_shared_observable_sum() {
    for (kind, rank, cases) in [("independent", 2.0, 1.0), ("parallel", 1.0, 2.0)] {
        let json = text(&run(kind, &controls()));
        assert!(json.contains("frankensim-equilibrium-sensitivity-v1"));
        assert!(scalar(&json, "objective") < 1e-20);
        assert_eq!(scalar(&json, "numerical_rank"), rank);
        assert_eq!(scalar(&json, "case_solves"), cases);
        assert_eq!(scalar(&json, "evaluations"), 1.0);
        assert_eq!(scalar(&json, "observation_adjoints"), 2.0);
        assert!(scalar(&json, "jacobian_scale") > 0.1);
        assert!(json.contains("\"parameter_directions\":") && json.contains("\"d_displacement_d_x\":"));
        if rank == 1.0 { assert!(json.contains("\"condition_number\":null")); }
        else { assert!((scalar(&json, "condition_number") - 2.0).abs() < 1e-9); }
    }
}

#[test]
fn named_point_order_replay_and_input_files_are_preserved() {
    let model_path = examples().join("sensitivity-independent.model");
    let design_path = examples().join("sensitivity-independent.fit");
    let model = std::fs::read(&model_path).unwrap(); let design = std::fs::read(&design_path).unwrap();
    let first = run("independent", &controls()); text(&first);
    let mut reverse = controls(); reverse[7] = "right".into(); reverse[10] = "left".into();
    assert_eq!(first.stdout, run("independent", &reverse).stdout);
    assert_eq!(first.stdout, run("independent", &controls()).stdout);
    assert_eq!(model, std::fs::read(model_path).unwrap()); assert_eq!(design, std::fs::read(design_path).unwrap());
}

#[test]
fn missing_points_or_adjoint_capacity_and_unsupported_precision_refuse() {
    for (flag, value) in [("--max-adjoints", "1"), ("--max-observations", "1"), ("--rank-relative-tolerance", "0.00000001")] {
        let mut c = controls(); let i = c.iter().position(|s| s == flag).unwrap(); c[i+1] = value.into();
        let out = run("independent", &c); assert!(!out.status.success()); assert!(out.stdout.is_empty());
        if flag != "--rank-relative-tolerance" { assert!(String::from_utf8_lossy(&out.stderr).contains("attempted 0 cases and 0 observation adjoints")); }
    }
    let mut c = controls(); c.truncate(c.len()-3);
    let out = run("independent", &c); assert!(!out.status.success()); assert!(out.stdout.is_empty());
    let mut c = controls(); c.extend(["--point-x", "left", "1"].map(str::to_owned));
    let out = run("independent", &c); assert!(!out.status.success()); assert!(out.stdout.is_empty());
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .args(["missing.model", "missing.fit", "--rank-relative-tolerance", "0.00000001", "--max-observations", "16", "--max-adjoints", "16"])
        .output().unwrap();
    assert!(!out.status.success()); assert!(String::from_utf8_lossy(&out.stderr).contains("local information requires"));
}

#[test]
fn zero_weight_rows_still_have_raw_derivatives_but_no_fitting_information() {
    let design = std::fs::read_to_string(examples().join("sensitivity-independent.fit")).unwrap();
    let path = temporary("zero-weights.fit", &design.replace("0.01 1", "0.01 0"));
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .arg(examples().join("sensitivity-independent.model")).arg(path).args(controls()).output().unwrap();
    let json = text(&out);
    assert_eq!(scalar(&json, "numerical_rank"), 0.0); assert_eq!(scalar(&json, "jacobian_scale"), 0.0);
    assert_eq!(scalar(&json, "observation_adjoints"), 2.0);
    let row = json.split("\"d_displacement_d_x\":[").nth(1).unwrap();
    assert!(row.split(',').next().unwrap().parse::<f64>().unwrap().abs() > 0.01);
}

#[test]
fn late_physical_refusal_and_contact_switches_publish_no_partial_analysis() {
    let original = std::fs::read_to_string(examples().join("sensitivity-parallel.fit")).unwrap();
    let path = temporary("failed-case.fit", &original.replace("load 0 0 2\n", "load 0 0 1000000000\n"));
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .arg(examples().join("sensitivity-parallel.model")).arg(path).args(controls()).output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("attempted 2 cases and 1 observation adjoints"));
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .arg(examples().join("contact-onset.model")).arg(examples().join("contact-onset.fit"))
        .args(["--rank-relative-tolerance", "0.00001", "--max-observations", "16", "--max-adjoints", "16", "--point-x", "force-N", "0"])
        .output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
}
