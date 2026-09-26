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
    for extra in [vec!["--checkpoint","must-not-be-created.bin"],vec!["--resume","missing.bin"],
        vec!["--qmc-replicates","2"],vec!["--sensitivity","sobol"]] {
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
