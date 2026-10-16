//! Public subprocess tests for the normalized thermal study CLI.
#![deny(unsafe_code)]

use std::process::{Command, Output};

fn arguments(command: &str, steps: usize) -> Vec<String> {
    format!("{command} --units normalized --model-version 1 --hole 0.35,0.5,0.09 --hole 0.65,0.5,0.09 --level 4 --steps {steps} --step-size 0 --area 0.95 --r-min 0.03 --r-max 0.12 --max-base-cells 256 --max-evaluations {steps}")
        .split_whitespace().map(str::to_string).collect()
}

fn call(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee")).args(args).output().expect("run fs-marquee")
}

fn stdout(output: &Output) -> &str {
    std::str::from_utf8(&output.stdout).expect("UTF-8 stdout")
}

fn replace(args: &mut [String], key: &str, value: &str) {
    let index = args.iter().position(|arg| arg == key).expect("fixture option");
    args[index + 1] = value.to_string();
}

fn trace(output: &Output) -> String {
    let prefix = "\"trace_hash\":\"";
    let text = stdout(output);
    let start = text.find(prefix).expect("complete trace hash") + prefix.len();
    text[start..].split('"').next().expect("trace field").to_string()
}

#[test]
fn check_projects_without_a_solve_or_a_complete_claim() {
    let output = call(&arguments("check", 0));
    assert!(output.status.success(), "{:?}", output);
    let text = stdout(&output);
    assert!(text.contains("\"event\":\"admitted\""));
    assert!(!text.contains("\"event\":\"iteration\""));
    assert!(!text.contains("\"event\":\"complete\""));
    assert!(text.contains("\"rng\":\"none\""));
    assert!(text.contains("\"units\":\"normalized\""));
}

#[test]
fn zero_step_run_has_no_fabricated_objective() {
    let output = call(&arguments("run", 0));
    assert!(output.status.success(), "{:?}", output);
    assert!(stdout(&output).contains("\"accepted_compliance\":null"));
    assert_eq!(stdout(&output).lines().count(), 2);
    assert_eq!(trace(&output).len(), 64);
}

#[test]
fn invalid_options_and_geometry_refuse_before_output() {
    for (key, value) in [
        ("--hole", "0.35,0.5,NaN"), ("--hole", "0.01,0.5,0.1"),
        ("--hole", "0.35,0.5,0.1,1"), ("--level", "64"),
        ("--steps", "257"), ("--r-min", "0"), ("--r-max", "inf"),
        ("--area", "1"), ("--units", "SI"), ("--model-version", "2"),
    ] {
        let mut args = arguments("check", 0);
        replace(&mut args, key, value);
        let output = call(&args);
        assert_eq!(output.status.code(), Some(2), "{key}={value}: {output:?}");
        assert!(output.stdout.is_empty());
    }
    let mut args = arguments("check", 0);
    args.extend(["--level".to_string(), "4".to_string()]);
    assert_eq!(call(&args).status.code(), Some(2));
}

#[test]
fn resource_budgets_refuse_before_solver_work() {
    let mut args = arguments("run", 1);
    replace(&mut args, "--max-base-cells", "255");
    let output = call(&args);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cell_budget_exceeded"));
    replace(&mut args, "--max-base-cells", "256");
    replace(&mut args, "--max-evaluations", "0");
    let output = call(&args);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("evaluation_budget_exceeded"));
}

#[test]
fn armijo_budget_accounts_for_the_current_solve_and_all_trials() {
    let mut args = arguments("check", 1);
    replace(&mut args, "--step-size", "0.01");
    replace(&mut args, "--max-evaluations", "9");
    assert_eq!(call(&args).status.code(), Some(2));
    replace(&mut args, "--max-evaluations", "10");
    assert!(call(&args).status.success());
}

#[test]
fn actual_pde_run_replays_and_detects_a_wrong_trace() {
    let args = arguments("run", 1);
    let first = call(&args);
    assert!(first.status.success(), "{first:?}");
    assert!(stdout(&first).contains("\"event\":\"iteration\""));
    assert!(!stdout(&first).contains("\"accepted_compliance\":null"));
    assert!(stdout(&first).contains("\"evidence\":\"estimated\""));
    let mut replay_args = args.clone();
    replay_args.extend(["--expect-trace".to_string(), trace(&first)]);
    let second = call(&replay_args);
    assert!(second.status.success(), "{second:?}");
    assert_eq!(first.stdout, second.stdout);
    replay_args.pop();
    replay_args.push("0".repeat(64));
    let mismatch = call(&replay_args);
    assert_eq!(mismatch.status.code(), Some(5));
    assert!(!stdout(&mismatch).contains("\"event\":\"complete\""));
}

#[test]
fn affine_load_is_part_of_replay_identity() {
    let base = call(&arguments("run", 0));
    assert!(base.status.success());
    let mut args = arguments("run", 0);
    args.extend(["--source".to_string(), "1,0.5,0".to_string()]);
    let changed = call(&args);
    assert!(changed.status.success());
    assert_ne!(trace(&base), trace(&changed));
    replace(&mut args, "--source", "1,-2,0");
    let invalid = call(&args);
    assert_eq!(invalid.status.code(), Some(2));
    assert!(invalid.stdout.is_empty());
}
