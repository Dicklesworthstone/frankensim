//! Actual-command tests. Reference values use an independent dense P1 FEM
//! calculation with air mixing eliminated analytically, not a second call to
//! the production staggered solver. The materials in the fixture are synthetic.
#[allow(dead_code)]
#[path="../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::fs;
use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::time::{SystemTime,UNIX_EPOCH};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static FIXTURE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/nonlinear-contact-pulse.json"))));
const POLICY:&str=r#""nonlinear":{"max_iterations":32,"residual_rtol":1e-10,"residual_atol_j":1e-10,"armijo_c":1e-4,"shrink":0.5,"max_backtracks":24},"#;
fn scratch(name:&str)->PathBuf {
    let nonce=SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("fs-nonlinear-transient-{name}-{}-{nonce}",std::process::id()));
    fs::create_dir(&path).unwrap();path
}
fn run(path:&Path)->Output {
    Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json","cooling-network"])
        .arg(path).output().unwrap()
}
fn document(output:&Output)->J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn close(a:f64,b:f64,tol:f64){assert!((a-b).abs()<=tol,"{a} != {b}");}
fn write(dir:&Path,name:&str,text:&str)->PathBuf {
    let path=dir.join(name);fs::write(&path,text).unwrap();path
}

#[test]
fn nonlinear_contact_pulse_matches_independent_endpoint_fem_and_energy() {
    let dir=scratch("oracle");
    let path=write(&dir,"pulse.json",FIXTURE);
    let output=run(&path);let doc=document(&output);
    let trajectory=doc.get("transient").unwrap();
    close(trajectory.f64_field("sampled_peak_objective_k").unwrap(),306.1651033635,3e-5);
    close(trajectory.f64_field("sampled_peak_time_s").unwrap(),30.0,1e-10);
    close(doc.path(&["objective","value_k"]).and_then(J::as_f64).unwrap(),301.4152325108374,3e-5);
    close(trajectory.f64_field("stored_energy_change_j").unwrap(),536.0649695201166,3e-5);
    close(trajectory.f64_field("air_energy_gain_j").unwrap(),63.93503047989178,3e-5);
    close(trajectory.f64_field("input_energy_j").unwrap(),600.0,1e-7);
    assert!(trajectory.f64_field("energy_residual_j").unwrap().abs()<=150.0e-7);
    let nonlinear=trajectory.get("nonlinear").unwrap();
    assert_eq!(nonlinear.str_field("method"),Some("endpoint-newton-fgmres"));
    assert!(nonlinear.f64_field("newton_updates").unwrap()>75.0);
    assert!(nonlinear.f64_field("worst_accepted_residual_ratio").unwrap()<=1.0);
    assert_eq!(nonlinear.f64_field("solid_solves"),trajectory.f64_field("total_solid_solves"));
    assert_eq!(doc.get("contacts").unwrap().as_array().unwrap().len(),1);
    assert_eq!(run(&path).stdout,output.stdout,"same deterministic schedule must replay");
}

#[test]
fn a_constant_curve_matches_the_original_linear_material_trajectory() {
    let dir=scratch("linear-control");
    let constant=FIXTURE.replace("\"conductivity_w_m_k\":[2,110]","\"conductivity_w_m_k\":[20,20]")
        .replace("\"conductivity_w_m_k\":[1,7]","\"conductivity_w_m_k\":[2,2]");
    assert_ne!(constant,FIXTURE);
    let doc=document(&run(&write(&dir,"constant.json",&constant)));
    let trajectory=doc.get("transient").unwrap();
    close(trajectory.f64_field("sampled_peak_objective_k").unwrap(),306.3431585002163,3e-5);
    close(doc.path(&["objective","value_k"]).and_then(J::as_f64).unwrap(),301.46723362812094,3e-5);
    close(trajectory.f64_field("stored_energy_change_j").unwrap(),536.2871880825311,3e-5);
    // Even a flat sampled table retains its finite validity span. The old
    // linear API intentionally refuses tables; use real scalar declarations.
    let legacy=constant.replace(POLICY,"")
        .replace(r#""conductivity_curve":{"temperature_k":[280,400],"conductivity_w_m_k":[20,20]}"#,r#""conductivity_w_m_k":20"#)
        .replace(r#""conductivity_curve":{"temperature_k":[280,400],"conductivity_w_m_k":[2,2]}"#,r#""conductivity_w_m_k":2"#);
    assert!(!legacy.contains("conductivity_curve"));
    let linear=document(&run(&write(&dir,"legacy.json",&legacy)));
    let a=doc.get("solid_temperatures_k").unwrap().as_array().unwrap();
    let b=linear.get("solid_temperatures_k").unwrap().as_array().unwrap();
    for (a,b) in a.iter().zip(b){close(a.as_f64().unwrap(),b.as_f64().unwrap(),1e-7);}
}

#[test]
fn missing_invalid_and_exhausted_policies_never_publish_a_trajectory() {
    let dir=scratch("refusals");
    assert!(FIXTURE.contains(POLICY));
    let cases=[
        FIXTURE.replace(POLICY,""),
        FIXTURE.replace("\"max_iterations\":32","\"max_iterations\":0"),
        FIXTURE.replace("\"max_iterations\":32","\"max_iterations\":1"),
        FIXTURE.replace("\"shrink\":0.5","\"shrink\":1.0"),
        FIXTURE.replace("\"residual_atol_j\":1e-10","\"residual_atol_j\":-1.0"),
        FIXTURE.replace("\"temperature_k\":[280,400]","\"temperature_k\":[100,200]"),
    ];
    for (i,text) in cases.iter().enumerate(){
        assert_ne!(text,FIXTURE);
        let output=run(&write(&dir,&format!("invalid-{i}.json"),text));
        assert!(!output.status.success());
        assert!(output.stdout.is_empty(),"a rejected endpoint must not publish a partial history");
        if i==0 {assert!(String::from_utf8_lossy(&output.stderr).contains("transient.nonlinear"));}
        if i==2 {assert_eq!(output.status.code(),Some(i32::from(fs_cli::exit::BUDGET)));}
    }
}

#[test]
fn adaptive_trials_use_endpoint_materials_and_count_discarded_solve_work() {
    let dir=scratch("adaptive");
    let text=FIXTURE.replace("\"duration_s\":30","\"duration_s\":4")
        .replace("\"duration_s\":120","\"duration_s\":4")
        .replace("\"max_step_s\":2",r#""adaptive":{"absolute_tolerance_k":0.01,"relative_tolerance":0.001,"minimum_trial_step_s":0.00001,"max_trials":1000},"max_step_s":2"#);
    let doc=document(&run(&write(&dir,"adaptive.json",&text)));
    let trajectory=doc.get("transient").unwrap();
    assert_ne!(trajectory.get("adaptive"),Some(&J::Null));
    close(trajectory.f64_field("input_energy_j").unwrap(),80.0,1e-7);
    assert!(trajectory.f64_field("energy_residual_j").unwrap().abs()<=8.0e-7);
    let nonlinear=trajectory.get("nonlinear").unwrap();
    assert_eq!(nonlinear.f64_field("solid_solves"),trajectory.f64_field("total_solid_solves"));
    assert!(nonlinear.f64_field("solid_solves").unwrap()>trajectory.f64_field("steps").unwrap());
    assert!(nonlinear.f64_field("worst_accepted_residual_ratio").unwrap()<=1.0);
}
