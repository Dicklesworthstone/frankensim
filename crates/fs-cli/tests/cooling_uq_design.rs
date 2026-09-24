#![cfg(unix)]

#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::ffi::OsString;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static BASE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/fan-hotspot.json"))));
// Compacted so fixture edits are independent of the example's formatting
// (a 2026-09-22 reformat silently turned every compact-spelled edit into a no-op).
static TRANSIENT: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-repeated-contact-pulse.json"))));
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn dir() -> PathBuf {
    let path = std::env::temp_dir().join(format!("frankensim-uq-design-{}-{}",
        std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    fs::create_dir(&path).unwrap();
    path
}
fn decoded(output: &Output) -> J {
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn cooling(text: &str) -> J {
    let mut child = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    decoded(&output)
}
fn n(value: &J, key: &str) -> f64 { value.f64_field(key).unwrap() }
fn objective(value: &J) -> f64 { value.path(&["objective","value_k"]).unwrap().as_f64().unwrap() }
fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }
fn policy(power: bool) -> Vec<OsString> {
    args(&[if power {"--power-candidates"}else{"--fan-speed-candidates"},"0.5,1,1.5",
        "--compliance-probability","0.5","--confidence-alpha","0.05","--min-decision-samples","32"])
}
fn uq(limit: f64, transient: bool) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":64,"wall_seconds":180,"temperature_limit_k":{limit},"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"{}"{} }},"distribution":{{"kind":"uniform","lo":{},"hi":{} }} }}] {} }}"#,
        if transient {"initial-temperature"}else{"component-power"},
        if transient {""}else{",\"component\":\"chip\""},
        if transient {300.0}else{1.0},if transient {300.0}else{1.0},
        if transient {",\"qoi\":{\"kind\":\"transient-sampled-peak\"}"}else{""})
}
fn run(base: &str, plan: &str, options: &[OsString]) -> Output {
    let root = dir();
    let base_path = root.join("base.json"); let plan_path = root.join("uq.json");
    fs::write(&base_path,base).unwrap(); fs::write(&plan_path,plan).unwrap();
    Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network-uq"]).arg(base_path).arg(plan_path)
        .args(options).output().unwrap()
}
fn with_paths(mut options: Vec<OsString>, saved: Option<&PathBuf>, output: &PathBuf,
    chunk: Option<usize>) -> Vec<OsString> {
    if let Some(saved) = saved { options.push("--resume".into());options.push(saved.as_os_str().to_owned()); }
    options.push("--checkpoint".into());options.push(output.as_os_str().to_owned());
    if let Some(chunk) = chunk { options.push("--max-new-samples".into());options.push(chunk.to_string().into()); }
    options
}
fn scaled_power(text: &str, value: f64) -> String {
    assert!(text.contains("\"total_w\":1") && text.contains("\"watts\":1"));
    text.replace("\"total_w\":1",&format!("\"total_w\":{value}"))
        .replace("\"watts\":1",&format!("\"watts\":{value}"))
}
fn power_limit() -> f64 { objective(&cooling(&scaled_power(BASE,1.25))) }

#[test]
fn actual_workload_and_fan_candidates_use_family_confidence_and_sampled_inputs() {
    for power in [true,false] {
        let limit = if power { power_limit() } else {
            assert!(BASE.contains("\"speed_ratio\":1,"));
            let low = objective(&cooling(&BASE.replace("\"speed_ratio\":1,","\"speed_ratio\":0.5,")));
            let high = objective(&cooling(BASE));
            assert!(low > high); 0.5*(low+high)
        };
        let output = run(BASE,&uq(limit,false),&policy(power));
        assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
        let result = decoded(&output);
        assert_eq!(result.str_field("status"),Some("selected"));
        assert_eq!(n(&result,"selected_multiplier"),1.0);
        assert!(n(&result,"alpha_per_candidate")*3.0 <= n(&result,"family_alpha"));
        let rows = result.get("candidates").unwrap().as_array().unwrap();
        let rejected = if power {2}else{0}; let unused = if power {0}else{2};
        assert_eq!(rows[rejected].str_field("decision"),Some("below-probability-target"));
        assert_eq!(rows[1].str_field("decision"),Some("meets-probability-target"));
        assert_eq!(n(&rows[unused],"samples_evaluated"),0.0);
        assert_eq!(n(&result,"samples_evaluated"),64.0);
    }
}

#[test]
fn cross_candidate_chunks_and_terminal_resume_reproduce_result_and_checkpoint_bytes() {
    let plan = uq(power_limit(),false); let root = dir();
    let full = root.join("full.bin"); let first = root.join("first.bin");
    let second = root.join("second.bin"); let last = root.join("last.bin");
    let complete = run(BASE,&plan,&with_paths(policy(true),None,&full,None));
    assert!(complete.status.success(),"{}",String::from_utf8_lossy(&complete.stderr));
    let a = run(BASE,&plan,&with_paths(policy(true),None,&first,Some(11)));
    assert_eq!(a.status.code(),Some(6)); assert_eq!(n(&decoded(&a),"samples_evaluated"),11.0);
    let b = run(BASE,&plan,&with_paths(policy(true),Some(&first),&second,Some(26)));
    assert_eq!(b.status.code(),Some(6));
    let rows = decoded(&b); let rows = rows.get("candidates").unwrap().as_array().unwrap();
    assert_eq!(n(&rows[2],"samples_evaluated"),32.0);
    assert_eq!(n(&rows[1],"samples_evaluated"),5.0);
    let c = run(BASE,&plan,&with_paths(policy(true),Some(&second),&last,None));
    assert!(c.status.success(),"{}",String::from_utf8_lossy(&c.stderr));
    assert_eq!(c.stdout,complete.stdout); assert_eq!(fs::read(&last).unwrap(),fs::read(&full).unwrap());
    let terminal = root.join("terminal.bin");
    let tiny_time = plan.replace("\"wall_seconds\":180","\"wall_seconds\":1e-12");
    let d = run(BASE,&tiny_time,&with_paths(policy(true),Some(&last),&terminal,Some(0)));
    assert!(d.status.success()); assert_eq!(d.stdout,complete.stdout);
    assert_eq!(fs::read(&terminal).unwrap(),fs::read(&full).unwrap());

    let mut changed = policy(true); changed[1] = "0.5,1,1.6".into();
    let refused = root.join("wrong-family.bin");
    let e = run(BASE,&plan,&with_paths(changed,Some(&first),&refused,None));
    assert!(!e.status.success()); assert!(e.stdout.is_empty()); assert!(!refused.exists());
    let before = fs::read(&last).unwrap();
    let e = run(BASE,&plan,&with_paths(policy(true),Some(&last),&last,None));
    assert!(!e.status.success()); assert_eq!(before,fs::read(&last).unwrap());
}

#[test]
fn repeated_trajectory_selection_uses_the_all_cycle_peak_not_the_final_temperature() {
    let adjoint = "\"adjoint\":{\"qoi\":\"sampled-peak\",\"max_checkpoint_bytes\":1048576},";
    assert!(TRANSIENT.contains(adjoint));
    let base = TRANSIENT.replace(adjoint,"");
    let higher = base.replace("\"power_scale\":1,","\"power_scale\":1.25,");
    let reference = cooling(&higher);
    let limit = reference.path(&["repeated_cycles","sampled_peak_objective_k"]).unwrap().as_f64().unwrap();
    assert!(limit > objective(&reference));
    let output = run(&base,&uq(limit,true),&policy(true));
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    let result = decoded(&output);
    assert_eq!(n(&result,"selected_multiplier"),1.0);
    assert_eq!(result.path(&["qoi","cycles"]).unwrap().as_f64(),Some(3.0));
    assert_eq!(result.path(&["qoi","observation"]).unwrap().as_str(),Some("one-completed-trajectory"));
}

#[test]
fn unresolved_sample_budgets_and_model_failures_never_publish_a_selected_design() {
    let plan = uq(power_limit(),false);
    let small = plan.replace("\"samples\":64","\"samples\":2");
    let mut options = policy(true); *options.last_mut().unwrap() = "2".into();
    let out = run(BASE,&small,&options);
    assert_eq!(out.status.code(),Some(6));
    let result = decoded(&out);
    assert_eq!(result.str_field("status"),Some("inconclusive"));
    assert_eq!(result.get("selected_multiplier"),Some(&J::Null));
    assert_eq!(n(&result,"samples_evaluated"),6.0);

    let all_bad = run(BASE,&uq(299.0,false),&policy(true));
    assert!(all_bad.status.success());
    assert_eq!(decoded(&all_bad).str_field("status"),Some("no-qualified-candidate"));
    assert_eq!(decoded(&all_bad).get("selected_multiplier"),Some(&J::Null));

    let failed = uq(400.0,false).replace("\"kind\":\"uniform\",\"lo\":1,\"hi\":1",
        "\"kind\":\"gaussian\",\"mean\":-1,\"std_dev\":0");
    assert_ne!(failed,uq(400.0,false));
    let checkpoint = dir().join("failed.bin");
    let out = run(BASE,&failed,&with_paths(policy(true),None,&checkpoint,None));
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(fs::read(&checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED"));
}
