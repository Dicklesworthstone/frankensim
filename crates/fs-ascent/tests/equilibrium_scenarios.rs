#![cfg(feature = "equilibrium-design")]
//! Actual-command regressions for finite-scenario physical inverse design.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const MODEL: &str = include_str!("../../fs-couple/examples/equilibrium-design.model");
const DESIGN: &str = include_str!("../../fs-couple/examples/equilibrium-response-limits.fit");
const SCENARIOS: &str = include_str!("../../fs-couple/examples/equilibrium-response-tolerances.scenarios");
static SERIAL: AtomicUsize = AtomicUsize::new(0);
fn inputs(model: &str, design: &str, scenarios: &str) -> (PathBuf,PathBuf,PathBuf) {
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("equilibrium-scenarios-{}-{stamp}-{}",std::process::id(),SERIAL.fetch_add(1,Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();
    let m=dir.join("model");let d=dir.join("design");let s=dir.join("scenarios");
    std::fs::write(&m,model).unwrap();std::fs::write(&d,design).unwrap();std::fs::write(&s,scenarios).unwrap();(m,d,s)
}
fn run(m:&Path,d:&Path,s:Option<&Path>,extra:&[&str])->Output {
    let mut command=Command::new(env!("CARGO_BIN_EXE_equilibrium_fit"));command.arg(m).arg(d);
    if let Some(s)=s {command.arg("--scenarios").arg(s);}
    command.args(extra).output().unwrap()
}
fn success(output:Output)->String {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    String::from_utf8(output.stdout).unwrap()
}
fn number(text:&str,key:&str)->f64 {
    text.split_once(key).unwrap().1.split([',','}',']']).next().unwrap().parse().unwrap()
}
fn load_parameter(text:&str)->f64 {number(text.split_once("\"name\":\"load-N\"").unwrap().1,"\"value\":")}

#[test]
fn tolerance_scenarios_change_the_design_and_keep_every_realization_under_the_force_cap() {
    let (m,d,s)=inputs(MODEL,DESIGN,SCENARIOS);
    let nominal=success(run(&m,&d,None,&[]));
    let robust=success(run(&m,&d,Some(&s),&[]));
    assert!(robust.contains("\"scope\":\"local-static-finite-scenario-minimax\""));
    assert!(robust.contains("\"converged\":true"));
    assert_eq!(robust.lines().count(),1);
    assert!((load_parameter(&nominal)-0.9077770876399966).abs()<1e-6);
    assert!((load_parameter(&robust)-0.8077770876399966).abs()<1e-6);
    assert!((number(&robust,"\"worst_objective\":")-0.9632655963309).abs()<1e-6);
    assert!(number(&robust,"\"epigraph_violation\":")<=1e-8);
    assert_eq!(robust.matches("\"quantity\":\"contact-force\"").count(),3);
    assert_eq!(robust.matches("\"predicted_m\":").count(),3);
    let high=robust.split_once("\"name\":\"high-load\"").unwrap().1;
    let cap=high.split_once("\"name\":\"normal-cap\"").unwrap().1;
    assert!((number(cap,"\"value\":")-0.8).abs()<1e-7);
    assert!(number(cap,"\"multiplier_normalized\":")>1.0);
    assert!(number(&robust,"\"physical_evaluations_including_audit\":")<=256.0);
    assert_eq!(success(run(&m,&d,None,&[])),nominal,"nominal invocation remains unchanged");
}

#[test]
fn unrelated_linear_model_reports_both_active_worst_case_multipliers() {
    let model="frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 10 1000 1000 1000000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\nmode 2 0.1 0 0 0 0\nport 0 1\ncoupling_limits 0 16384 0.9 1000 1000000 10000 1e-10 1e-11 1e-9\nconnections 0\nevents 0\n";
    let design="frankensim-equilibrium-design-v1\npreload_limits 0 1 16384\nsensitivity_limits 0 16384 16384 0\ndesign_limits 1 1 1 2\ncases 1\ncase response 1 1\nload 0 0 4\ntarget 0 0 1 1 1\nvariables 1\nvariable load-N 4 4 1 8 1\nbind actuator-force 0 0\n";
    let scenarios="frankensim-equilibrium-scenarios-v1\nscenarios 2\nscenario low\noffset load-N -0.8\nscenario high\noffset load-N 2.4\n";
    let (m,d,s)=inputs(model,design,scenarios);
    let result=success(run(&m,&d,Some(&s),&["--tolerance","1e-9"]));
    assert!((load_parameter(&result)-3.2).abs()<1e-7);
    assert!((number(&result,"\"worst_objective\":")-0.08).abs()<1e-9);
    assert_eq!(result.matches("\"epigraph_multiplier\":").count(),2);
    for segment in result.split("\"epigraph_multiplier\":").skip(1) {
        let dual=segment.split(',').next().unwrap().parse::<f64>().unwrap();assert!((dual-0.5).abs()<1e-7);
    }
}

#[test]
fn final_complete_family_fits_inside_the_physical_budget_and_invalid_inputs_publish_nothing() {
    let (m,d,s)=inputs(MODEL,DESIGN,SCENARIOS);
    let stopped=success(run(&m,&d,Some(&s),&["--evaluations","6"]));
    assert!(stopped.contains("\"stop\":\"EvaluationLimit\""));assert!(stopped.contains("\"converged\":false"));
    assert_eq!(number(&stopped,"\"iterations\":"),0.0);
    assert_eq!(number(&stopped,"\"physical_evaluations_including_audit\":"),6.0);
    assert_eq!(number(&stopped,"\"case_solves\":"),6.0);
    assert_eq!(number(&stopped,"\"ensemble_evaluations_including_audit\":"),2.0);
    assert!(number(&stopped,"\"violation\":")>0.0);
    for flags in [vec!["--evaluations","5"],vec!["--max-kkt-dimension","5"],vec!["--scenarios","duplicate"]] {
        let bad=run(&m,&d,Some(&s),&flags);assert!(!bad.status.success());assert!(bad.stdout.is_empty());
    }
    for malformed in [SCENARIOS.replace("scenarios 3","scenarios 33"),SCENARIOS.replace("load-N","unknown"),
        SCENARIOS.replace("-0.1","NaN"),SCENARIOS.replace("high-load","low-load"),
        SCENARIOS.replace("offset load-N 0.1\n",""),format!("{SCENARIOS}ignored\n")] {
        let (m,d,s)=inputs(MODEL,DESIGN,&malformed);
        let bad=run(&m,&d,Some(&s),&[]);assert!(!bad.status.success());assert!(bad.stdout.is_empty());
    }
}

#[test]
fn complete_realizations_relocate_and_user_labels_are_escaped_without_input_mutation() {
    let scenarios=SCENARIOS.replace("high-load","high-\"\\load");
    let (m,d,s)=inputs(MODEL,DESIGN,&scenarios);
    let first=success(run(&m,&d,Some(&s),&["--evaluations","6"]));
    assert!(first.contains("\"name\":\"high-\\\"\\\\load\""));
    let (other_m,other_d,other_s)=inputs(MODEL,DESIGN,&scenarios);
    assert_eq!(first,success(run(&other_m,&other_d,Some(&other_s),&["--evaluations","6"])));
    assert_eq!(std::fs::read_to_string(m).unwrap(),MODEL);
    assert_eq!(std::fs::read_to_string(d).unwrap(),DESIGN);
    assert_eq!(std::fs::read_to_string(s).unwrap(),scenarios);
}
