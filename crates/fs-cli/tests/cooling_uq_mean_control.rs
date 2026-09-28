//! Public-command regressions: actual coupled child solves, not fabricated QoIs.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::ffi::OsString;
use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::sync::atomic::{AtomicU64,Ordering};
use std::time::{SystemTime,UNIX_EPOCH};

const BASE: &str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/fan-correlated-hotspot.json"));
static NEXT: AtomicU64=AtomicU64::new(0);
fn directory()->PathBuf {
    let path=std::env::temp_dir().join(format!("fs-uq-mean-{}-{}-{}",std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos(),NEXT.fetch_add(1,Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap(); path
}
fn replace_once(text:&str,from:&str,to:&str)->String {
    assert_eq!(text.matches(from).count(),1,"fixture mutation must match exactly once: {from}");
    text.replacen(from,to,1)
}
fn command(args:Vec<OsString>)->Output {
    Command::new(env!("CARGO_BIN_EXE_frankensim")).arg("--json").args(args).output().unwrap()
}
fn document(output:Output)->J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn uq(dir:&Path,base:&str,plan:&str,control:bool)->J {
    let base_path=dir.join(if control {"controlled-base.json"}else{"raw-base.json"});
    let plan_path=dir.join(if control {"controlled-plan.json"}else{"raw-plan.json"});
    std::fs::write(&base_path,base).unwrap();std::fs::write(&plan_path,plan).unwrap();
    let mut args=vec!["cooling-network-uq".into(),base_path.into_os_string(),plan_path.into_os_string()];
    if control {args.extend(["--mean-control".into(),"adjoint".into()]);}
    document(command(args))
}
fn solve(dir:&Path,name:&str,base:&str)->f64 {
    let path=dir.join(name);std::fs::write(&path,base).unwrap();
    document(command(vec!["cooling-network".into(),path.into_os_string()]))
        .path(&["objective","value_k"]).unwrap().as_f64().unwrap()
}
fn plan(parameters:&str)->String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":16,"wall_seconds":120,"temperature_limit_k":312,"correlation":{{"kind":"independent"}},"parameters":[{parameters}]}}"#)
}
fn number(value:&J,key:&str)->f64 {value.f64_field(key).unwrap()}

#[test]
fn inlet_control_uses_distribution_mean_and_preserves_every_raw_result_field() {
    let dir=directory();
    let p=plan(r#"{"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"uniform","lo":308,"hi":312}}"#);
    let raw=uq(&dir,BASE,&p,false);
    let controlled=uq(&dir,BASE,&p,true);
    let cv=controlled.get("mean_control_variate").unwrap();
    let nominal=solve(&dir,"mean-inlet.json",&replace_once(BASE,"\"temperature_k\": 300","\"temperature_k\": 310"));
    assert!((number(cv,"mean_k")-nominal).abs()<1e-6);
    assert!((number(cv,"nominal_objective_k")-nominal).abs()<1e-6);
    assert!(number(cv,"adjusted_std_dev_k")<1e-6);
    assert!(number(cv,"adjusted_to_raw_variance_ratio")<1e-10);
    assert!(number(&raw,"std_dev_k")>0.1);
    assert_eq!(cv.get("parameter_means").unwrap().as_array().unwrap()[0].as_f64(),Some(310.0));
    assert_eq!(number(cv,"sample_model_evaluations"),16.0);
    assert_eq!(number(cv,"total_model_evaluations"),17.0);
    let mut untouched=controlled.clone();
    let J::Object(rows)=&mut untouched else {panic!("result object required")};
    rows.retain(|(name,_)|name!="mean_control_variate");
    assert_eq!(untouched,raw,"probability, bounds, quantiles and raw standard error must remain unchanged");
    let retry_dir=directory();
    assert_eq!(controlled,uq(&retry_dir,BASE,&p,true),"same seed and model must replay the complete result");
}

#[test]
fn fan_control_uses_the_actual_total_adjoint_in_declared_parameter_order() {
    let dir=directory();
    let p=plan(r#"{"target":{"kind":"fan-speed-ratio"},"distribution":{"kind":"uniform","lo":1.09,"hi":1.11}},{"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"uniform","lo":302,"hi":304}}"#);
    let controlled=uq(&dir,BASE,&p,true);
    let cv=controlled.get("mean_control_variate").unwrap();
    let g=cv.get("gradient_k_per_parameter_unit").unwrap().as_array().unwrap();
    let center=replace_once(BASE,"\"temperature_k\": 300","\"temperature_k\": 303");
    let plus=replace_once(&center,"\"speed_ratio\": 1,","\"speed_ratio\": 1.1001,");
    let minus=replace_once(&center,"\"speed_ratio\": 1,","\"speed_ratio\": 1.0999,");
    let fd=(solve(&dir,"fan-plus.json",&plus)-solve(&dir,"fan-minus.json",&minus))/0.0002;
    assert!((g[0].as_f64().unwrap()-fd).abs()<1e-4*fd.abs().max(1.0),"actual fan/flow/convection chain must match independent solves");
    assert!((g[1].as_f64().unwrap()-1.0).abs()<1e-7);
    assert!(number(cv,"adjusted_to_raw_variance_ratio")<0.01);
}

#[test]
fn unsupported_modes_refuse_before_reading_or_running_a_model() {
    for extra in [vec!["--qmc-replicates","2"],vec!["--sensitivity","sobol"],
        vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "8"]] {
        let dir=directory();
        let mut args:Vec<OsString>=vec!["cooling-network-uq".into(),dir.join("missing-base.json").into_os_string(),
            dir.join("missing-plan.json").into_os_string(),"--mean-control".into(),"adjoint".into()];
        args.extend(extra.into_iter().map(|value| {
            if value == "must-not-be-created.bin" { dir.join(value).into_os_string() }
            else { OsString::from(value) }
        }));
        let output=command(args);
        assert_eq!(output.status.code(),Some(2),"option admission must precede file/physics work");
        assert!(output.stdout.is_empty());
        assert!(!dir.join("must-not-be-created.bin").exists());
    }
}

fn recoverable_command(dir: &Path, extra: &[&str]) -> Output {
    let mut args: Vec<OsString> = vec!["cooling-network-uq".into(),
        dir.join("base.json").into_os_string(), dir.join("plan.json").into_os_string(),
        "--mean-control".into(), "adjoint".into()];
    for &arg in extra {
        args.push(if arg.ends_with(".bin") { dir.join(arg).into_os_string() } else { arg.into() });
    }
    command(args)
}
fn recovery_inputs(dir: &Path) -> String {
    let p = plan(r#"{"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"uniform","lo":308,"hi":312}}"#);
    std::fs::write(dir.join("base.json"), BASE).unwrap();
    std::fs::write(dir.join("plan.json"), &p).unwrap();
    p
}
fn progress(output: Output) -> J {
    assert_eq!(output.status.code(), Some(6), "{}", String::from_utf8_lossy(&output.stderr));
    let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert!(result.get("mean_control_variate").is_none(), "a partial prefix must not masquerade as a completed mean estimate");
    assert_eq!(result.str_field("mean_control_state"), Some("frozen-adjoint-retained"));
    result
}

#[test]
fn actual_adjoint_and_samples_survive_zero_chunk_multichunk_and_terminal_recovery() {
    let dir = directory();
    recovery_inputs(&dir);
    let zero = progress(recoverable_command(&dir, &["--checkpoint", "zero.bin", "--max-new-samples", "0"]));
    assert_eq!(number(&zero, "samples_evaluated"), 0.0);
    assert_eq!(number(&zero, "total_model_evaluations_this_run"), 1.0);
    assert_eq!(number(&zero, "nominal_forward_adjoint_evaluations_this_run"), 1.0);
    let zero_bytes = std::fs::read(dir.join("zero.bin")).unwrap();
    let partial = progress(recoverable_command(&dir, &["--resume", "zero.bin", "--checkpoint", "part.bin", "--max-new-samples", "5"]));
    assert_eq!(number(&partial, "samples_evaluated"), 5.0);
    assert_eq!(number(&partial, "total_model_evaluations"), 6.0);
    assert_eq!(number(&partial, "total_model_evaluations_this_run"), 5.0);
    assert_eq!(number(&partial, "nominal_forward_adjoint_evaluations_this_run"), 0.0);
    assert_eq!(partial.get("control_restored"), Some(&J::Bool(true)));
    assert_eq!(std::fs::read(dir.join("zero.bin")).unwrap(), zero_bytes);
    let part_bytes = std::fs::read(dir.join("part.bin")).unwrap();
    let resumed = document(recoverable_command(&dir, &["--resume", "part.bin", "--checkpoint", "complete.bin"]));
    assert_eq!(std::fs::read(dir.join("part.bin")).unwrap(), part_bytes);
    let whole = document(recoverable_command(&dir, &["--checkpoint", "whole.bin"]));
    assert_eq!(resumed, whole, "same raw samples, nominal record and complete controlled estimate");
    assert_eq!(std::fs::read(dir.join("complete.bin")).unwrap(), std::fs::read(dir.join("whole.bin")).unwrap());
    let terminal = document(recoverable_command(&dir, &["--resume", "complete.bin"]));
    assert_eq!(terminal, whole);
    assert_eq!(number(terminal.get("mean_control_variate").unwrap(), "total_model_evaluations"), 17.0);
}

#[test]
fn corrupted_or_rebound_adjoint_prefixes_refuse_before_creating_an_output() {
    let dir = directory();
    let p = recovery_inputs(&dir);
    progress(recoverable_command(&dir, &["--checkpoint", "zero.bin", "--max-new-samples", "0"]));
    let original = std::fs::read(dir.join("zero.bin")).unwrap();
    for (ordinal, offset) in [8, 16, 40, original.len() - 1].into_iter().enumerate() {
        let mut changed = original.clone(); changed[offset] ^= 1;
        let input = format!("bad-{ordinal}.bin");
        let output = format!("not-created-{ordinal}.bin");
        std::fs::write(dir.join(&input), changed).unwrap();
        let rejected = recoverable_command(&dir, &["--resume", &input, "--checkpoint", &output]);
        assert_eq!(rejected.status.code(), Some(4), "{}", String::from_utf8_lossy(&rejected.stderr));
        assert!(rejected.stdout.is_empty());
        assert!(!dir.join(output).exists());
    }
    // The same core plan does not legitimize a different fixed model input.
    let changed_base = replace_once(BASE, "\"temperature_k\": 300", "\"temperature_k\": 301");
    std::fs::write(dir.join("base.json"), changed_base).unwrap();
    assert_eq!(recoverable_command(&dir, &["--resume", "zero.bin", "--checkpoint", "different.bin"]).status.code(), Some(4));
    assert!(!dir.join("different.bin").exists());
    std::fs::write(dir.join("base.json"), BASE).unwrap();
    let changed_plan = replace_once(&p, "\"samples\":16", "\"samples\":17");
    std::fs::write(dir.join("plan.json"), changed_plan).unwrap();
    assert_eq!(recoverable_command(&dir, &["--resume", "zero.bin", "--checkpoint", "budget-changed.bin"]).status.code(), Some(4));
    assert!(!dir.join("budget-changed.bin").exists());
    // A per-invocation wall allowance can change without altering the sampled law.
    std::fs::write(dir.join("plan.json"), replace_once(&p, "\"wall_seconds\":120", "\"wall_seconds\":240")).unwrap();
    let copied = progress(recoverable_command(&dir, &["--resume", "zero.bin", "--checkpoint", "copied.bin", "--max-new-samples", "0"]));
    assert_eq!(number(&copied, "total_model_evaluations_this_run"), 0.0);
    assert_eq!(std::fs::read(dir.join("copied.bin")).unwrap(), original);
}

#[test]
fn controlled_recovery_never_overwrites_input_or_silently_changes_execution_mode() {
    let dir = directory(); recovery_inputs(&dir);
    progress(recoverable_command(&dir, &["--checkpoint", "zero.bin", "--max-new-samples", "0"]));
    let original = std::fs::read(dir.join("zero.bin")).unwrap();
    assert_eq!(recoverable_command(&dir, &["--resume", "zero.bin", "--checkpoint", "zero.bin"]).status.code(), Some(4));
    assert_eq!(std::fs::read(dir.join("zero.bin")).unwrap(), original);
    let raw = command(vec!["cooling-network-uq".into(), dir.join("base.json").into_os_string(),
        dir.join("plan.json").into_os_string(), "--resume".into(), dir.join("zero.bin").into_os_string(),
        "--checkpoint".into(), dir.join("raw-output.bin").into_os_string()]);
    assert_eq!(raw.status.code(), Some(4));
    assert!(!dir.join("raw-output.bin").exists());
    let raw_zero = command(vec!["cooling-network-uq".into(), dir.join("base.json").into_os_string(),
        dir.join("plan.json").into_os_string(), "--checkpoint".into(), dir.join("raw.bin").into_os_string(),
        "--max-new-samples".into(), "0".into()]);
    assert_eq!(raw_zero.status.code(), Some(6));
    assert_eq!(recoverable_command(&dir, &["--resume", "raw.bin", "--checkpoint", "late-control.bin"]).status.code(), Some(4));
    assert!(!dir.join("late-control.bin").exists());
    std::fs::write(dir.join("occupied.bin"), b"existing user data").unwrap();
    assert_eq!(recoverable_command(&dir, &["--checkpoint", "occupied.bin"]).status.code(), Some(4));
    assert_eq!(std::fs::read(dir.join("occupied.bin")).unwrap(), b"existing user data".to_vec());
}
