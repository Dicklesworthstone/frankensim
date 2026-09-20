//! Actual-binary tests: arbitrary files reach the existing bounded physics study.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const MODEL: &str = include_str!("../../fs-couple/examples/equilibrium-design.model");
const DESIGN: &str = include_str!("../../fs-couple/examples/equilibrium-design.fit");
static SERIAL: AtomicUsize = AtomicUsize::new(0);
fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("equilibrium-fit-{}-{stamp}-{}", std::process::id(), SERIAL.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap(); path
}
fn inputs(dir: &Path, model: &str, design: &str) -> (PathBuf, PathBuf) {
    let m = dir.join("model.performance"); let d = dir.join("design.fit");
    std::fs::write(&m, model).unwrap(); std::fs::write(&d, design).unwrap(); (m,d)
}
fn run(model: &Path, design: &Path, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_equilibrium_fit")).arg(model).arg(design).args(extra).output().unwrap()
}
fn text(result: Output) -> String {
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    String::from_utf8(result.stdout).unwrap()
}
fn number_after(text: &str, key: &str) -> f64 {
    text.split_once(key).unwrap().1.split([',','}',']']).next().unwrap().parse().unwrap()
}
fn parameter(text: &str, name: &str) -> f64 {
    number_after(text.split_once(&format!("\"name\":\"{name}\"")).unwrap().1, "\"value\":")
}

#[test]
fn complete_files_fit_shared_contact_parameters_and_relocate_without_changing_results() {
    let dir = directory(); let (m,d) = inputs(&dir, MODEL, DESIGN);
    let result = text(run(&m,&d,&[]));
    assert_eq!(result.lines().count(),1);
    assert!(result.contains("\"schema\":\"frankensim-equilibrium-fit-v1\""));
    assert!(number_after(&result,"\"objective\":") < 1e-8);
    for (name,expected) in [("support-N-per-m",600.0),("contact-N-per-m2",1.8e8),("gap-m",0.0002)] {
        assert!((parameter(&result,name)-expected).abs() < 1e-3*expected);
    }
    assert_eq!(result.matches("\"predicted_m\":").count(),6);
    assert_eq!(result.matches("\"active_contacts\":1").count(),3);
    assert!(number_after(&result,"\"evaluations_including_audit\":") <= 256.0);
    let (moved_m,moved_d) = inputs(&directory(), MODEL, DESIGN);
    assert_eq!(text(run(&moved_m,&moved_d,&[])),result);
    assert_eq!(std::fs::read_to_string(m).unwrap(),MODEL);
    assert_eq!(std::fs::read_to_string(d).unwrap(),DESIGN);
}

#[test]
fn unrelated_model_reaches_a_real_bound_optimum_with_a_nonzero_objective_gradient() {
    // One elastic coordinate, different frequency and shape, no normal contacts.
    // u = F/64 m; target requires F=2 N, but the design caps it at 1 N.
    let model = "frankensim-modal-performance-v2\nsample_rate_hz 24000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 1 100 1000 1000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\nmode 16 0.1 0 0 0 0\nport 0 2\ncoupling_limits 0 10000 0.9 1000 1000 1000 1e-10 1e-10 1e-8\nconnections 0\nevents 0\n";
    let design = "frankensim-equilibrium-design-v1\npreload_limits 0 1 10000\nsensitivity_limits 0 10000 10000 0\ndesign_limits 1 1 1 2\ncases 1\ncase different-rig 1 1\nload 0 0 0.5\ntarget 0 0 0.03125 0.015625 1\nvariables 1\nvariable force-N 0.5 0.5 0 1 1\nbind actuator-force 0 0\n";
    let (m,d) = inputs(&directory(),model,design);
    let result = text(run(&m,&d,&["--tolerance","1e-9"]));
    assert!(result.contains("\"converged\":true"));
    assert!(result.contains("\"stop\":\"Converged\""));
    assert!((parameter(&result,"force-N")-1.0).abs() < 1e-10);
    assert!((number_after(&result,"\"objective\":")-0.5).abs() < 1e-10);
    assert!(number_after(&result,"\"upper_multiplier_decision\":") > 0.49);
    assert!(number_after(&result,"\"stationarity\":") < 1e-9);
}

#[test]
fn budget_stop_reserves_a_complete_final_re_solve_and_is_not_labelled_converged() {
    let (m,d) = inputs(&directory(),MODEL,DESIGN);
    let result = text(run(&m,&d,&["--evaluations","2"]));
    assert!(result.contains("\"stop\":\"EvaluationLimit\""));
    assert!(result.contains("\"converged\":false"));
    assert_eq!(number_after(&result,"\"evaluations_including_audit\":"),2.0);
    assert_eq!(number_after(&result,"\"case_solves\":"),6.0);
    assert_eq!(number_after(&result,"\"iterations\":"),0.0);
    assert_eq!(parameter(&result,"support-N-per-m"),400.0);
    assert_eq!(result.matches("\"predicted_m\":").count(),6);
}

#[test]
fn malformed_inputs_fail_without_a_partial_result_and_names_are_json_escaped() {
    let (m,d) = inputs(&directory(),MODEL,&DESIGN.replace("target 1 0","target 99 0"));
    let bad = run(&m,&d,&[]);
    assert!(!bad.status.success()); assert!(bad.stdout.is_empty()); assert!(!bad.stderr.is_empty());
    let (m,d) = inputs(&directory(),MODEL,DESIGN);
    for flags in [vec!["--evaluations","1"],vec!["--iterations","1","--iterations","2"],vec!["--max-kkt-dimension","1"]] {
        let bad = run(&m,&d,&flags); assert!(!bad.status.success()); assert!(bad.stdout.is_empty());
    }
    let escaped = DESIGN.replace("load-0.8N","load-\"\\N");
    let (m,d) = inputs(&directory(),MODEL,&escaped);
    let result = text(run(&m,&d,&["--evaluations","2"]));
    assert!(result.contains("\"name\":\"load-\\\"\\\\N\""));
}

#[path = "equilibrium_fit/response_constraints.rs"]
mod response_constraints;
