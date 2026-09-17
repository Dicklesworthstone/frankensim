#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/radiative-contact-hotspot.json"));
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("fs-radiation-uq-{}-{}", std::process::id(), NEXT.fetch_add(1,Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap();
        std::fs::write(path.join("base.json"),BASE).unwrap();
        Self(path)
    }
    fn file(&self, name: &str) -> PathBuf { self.0.join(name) }
}
impl Drop for Scratch {
    fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
}
fn plan(vary: bool) -> String {
    let (elo,ehi,tlo,thi)=if vary {(0.75,0.95,285.0,305.0)}else{(0.85,0.85,290.0,290.0)};
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":600,"temperature_limit_k":335,"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"radiation-emissivity","surface":"first-face"}},"distribution":{{"kind":"uniform","lo":{elo},"hi":{ehi}}}},{{"target":{{"kind":"radiation-ambient-temperature","surface":"last-face"}},"distribution":{{"kind":"uniform","lo":{tlo},"hi":{thi}}}}]}}"#)
}
fn evaluate(dir: &Scratch, request: &str, extras: &[&str]) -> Output {
    let spec=dir.file("uq.json");std::fs::write(&spec,request).unwrap();
    Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json","cooling-network-uq"])
        .arg(dir.file("base.json")).arg(spec).args(extras).output().unwrap()
}
fn success(output: &Output) -> J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn s(path: &Path) -> &str { path.to_str().unwrap() }

#[test]
fn radiation_parameter_draws_reach_the_real_cooling_equations() {
    let dir=Scratch::new();
    let direct=Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json","cooling-network"])
        .arg(dir.file("base.json")).output().unwrap();
    let direct=success(&direct);
    let result=success(&evaluate(&dir,&plan(false),&[]));
    let expected=direct.path(&["objective","value_k"]).unwrap().as_f64().unwrap();
    assert!((result.f64_field("mean_k").unwrap()-expected).abs()<1e-10);
    assert_eq!(result.f64_field("std_dev_k"),Some(0.0));
    assert_eq!(result.f64_field("samples_evaluated"),Some(4.0));
    let variable=success(&evaluate(&dir,&plan(true),&[]));
    assert!(variable.f64_field("std_dev_k").unwrap()>1e-5,
        "varying radiation must not leave the nominal child input unchanged");
    let parameters=variable.get("parameters").unwrap().as_array().unwrap();
    assert_eq!(parameters[0].str_field("unit"),Some("1"));
    assert_eq!(parameters[1].str_field("unit"),Some("K"));
}

#[test]
fn complete_and_resumed_radiative_samples_have_identical_observation_bytes() {
    let dir=Scratch::new();let spec=plan(true);
    let full_path=dir.file("full.uqcp");let part_path=dir.file("part.uqcp");let done_path=dir.file("done.uqcp");
    let full=evaluate(&dir,&spec,&["--checkpoint",s(&full_path)]);success(&full);
    let part=evaluate(&dir,&spec,&["--checkpoint",s(&part_path),"--max-new-samples","2"]);
    assert_eq!(part.status.code(),Some(6));
    let before=std::fs::read(&part_path).unwrap();
    let done=evaluate(&dir,&spec,&["--resume",s(&part_path),"--checkpoint",s(&done_path)]);success(&done);
    assert_eq!(full.stdout,done.stdout);
    assert_eq!(std::fs::read(full_path).unwrap(),std::fs::read(done_path).unwrap());
    assert_eq!(before,std::fs::read(part_path).unwrap());
}

#[test]
fn out_of_range_emissivity_is_a_terminal_model_failure_not_a_filtered_draw() {
    let dir=Scratch::new();let checkpoint=dir.file("failed.uqcp");
    let invalid=r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":600,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"radiation-emissivity","surface":"first-face"},"distribution":{"kind":"gaussian","mean":1.1,"std_dev":0}}]}"#;
    let output=evaluate(&dir,invalid,&["--checkpoint",s(&checkpoint)]);
    assert!(!output.status.success());assert!(output.stdout.is_empty());
    assert!(!output.stderr.is_empty());
    assert!(std::fs::read(checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
}
