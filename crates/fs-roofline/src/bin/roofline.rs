//! Roofline harness CLI (plan §14.4 nightly lane).
//!
//! Usage:
//!   roofline [--n <elements>] [--warmup <k>] [--reps <k>] [--ledger <db>]
//!
//! Probes the machine axes, runs the default kernel registry, prints one
//! JSON line per kernel (plus the axes line and the §14.1 coverage table),
//! and — when `--ledger` is given — records the run as ledger provenance
//! and reports staleness for every registered kernel.

use fs_roofline::kernels::default_registry;
use fs_roofline::{MachineAxes, SECTION_14_1_TARGETS, run_registry, staleness};

fn fail(detail: &str) -> std::process::ExitCode {
    eprintln!("{{\"error\":\"Roofline\",\"detail\":\"{detail}\"}}");
    std::process::ExitCode::FAILURE
}

fn parse_flag(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1).cloned())
}

fn main() -> std::process::ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let n = match parse_flag(&args, "--n").map(|v| v.parse::<usize>()) {
        None => 1 << 22, // 32 MiB per f64 buffer: streams past every L2/L3
        Some(Ok(v)) if v > 0 => v,
        Some(_) => return fail("--n must be a positive integer"),
    };
    let warmup = match parse_flag(&args, "--warmup").map(|v| v.parse::<usize>()) {
        None => 2,
        Some(Ok(v)) => v,
        Some(Err(_)) => return fail("--warmup must be an integer"),
    };
    let reps = match parse_flag(&args, "--reps").map(|v| v.parse::<usize>()) {
        None => 9,
        Some(Ok(v)) if v > 0 => v,
        Some(_) => return fail("--reps must be a positive integer"),
    };

    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());

    let mut registry = default_registry(n);
    let results = run_registry(&mut registry, warmup, reps, &axes);
    for r in &results {
        println!("{}", r.to_jsonl());
    }
    for row in SECTION_14_1_TARGETS {
        println!(
            "{{\"target\":\"{}\",\"statement\":\"{}\",\"landed\":{}}}",
            row.kernel, row.statement, row.landed
        );
    }

    if let Some(db) = parse_flag(&args, "--ledger") {
        let ledger = match fs_ledger::Ledger::open(&db) {
            Ok(l) => l,
            Err(e) => return fail(&e.to_string().replace('"', "'")),
        };
        match fs_roofline::record_run(&ledger, &axes, &results) {
            Ok(op) => println!("{{\"ledgered\":true,\"op\":{op},\"db\":\"{db}\"}}"),
            Err(e) => return fail(&e.to_string().replace('"', "'")),
        }
        for r in &results {
            match staleness(&ledger, &r.kernel, axes.fingerprint) {
                Ok(s) => println!("{{\"kernel\":\"{}\",\"staleness\":\"{s:?}\"}}", r.kernel),
                Err(e) => return fail(&e.to_string().replace('"', "'")),
            }
        }
    }
    std::process::ExitCode::SUCCESS
}
