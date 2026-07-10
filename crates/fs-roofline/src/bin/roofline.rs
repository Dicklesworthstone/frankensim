//! Roofline harness CLI (plan §14.4 nightly lane).
//!
//! Usage:
//!   roofline [--n <elements>] [--warmup <k>] [--reps <k>] [--ledger <db>]
//!
//! Probes the machine axes, runs the default kernel registry, prints one
//! JSON line per kernel (plus the axes line and the §14.1 coverage table),
//! and — when `--ledger` is given — records the run as ledger provenance
//! and reports staleness for every registered kernel.

use fs_roofline::kernels::production_registry_with_ledger;
use fs_roofline::{MachineAxes, SECTION_14_1_TARGETS, run_is_citable, run_registry, staleness};

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(c));
            }
            c => escaped.push(c),
        }
    }
    escaped
}

fn fail(detail: &str) -> std::process::ExitCode {
    eprintln!(
        "{{\"error\":\"Roofline\",\"detail\":\"{}\"}}",
        json_escape(detail)
    );
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
        Some(Ok(v)) if v > 0 => v,
        Some(_) => return fail("--warmup must be a positive integer"),
    };
    let reps = match parse_flag(&args, "--reps").map(|v| v.parse::<usize>()) {
        None => 9,
        Some(Ok(v)) if v > 0 => v,
        Some(_) => return fail("--reps must be a positive integer"),
    };

    let ledger_path = parse_flag(&args, "--ledger");
    let tune_ledger = match ledger_path.as_deref() {
        Some(path) => match fs_ledger::Ledger::open(path) {
            Ok(ledger) => Some(ledger),
            Err(error) => return fail(&error.to_string()),
        },
        None => None,
    };

    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());

    let mut registry = production_registry_with_ledger(n, &axes, tune_ledger);
    let results = run_registry(&mut registry, warmup, reps, &axes);
    let post_axes = MachineAxes::probe();
    println!("{}", post_axes.to_jsonl());
    let citable = run_is_citable(&axes, &post_axes, &results);
    for r in &results {
        println!("{}", r.to_jsonl());
    }
    for row in SECTION_14_1_TARGETS {
        println!(
            "{{\"target\":\"{}\",\"statement\":\"{}\",\"landed\":{}}}",
            json_escape(row.kernel),
            json_escape(row.statement),
            row.landed
        );
    }

    // The registry owns fsqlite's deliberately !Send tune connection. Drop it
    // before reopening the same database for the atomic evidence transaction;
    // run_registry is synchronous, so no kernel work survives this point.
    drop(registry);
    if let Some(db) = ledger_path {
        let ledger = match fs_ledger::Ledger::open(&db) {
            Ok(l) => l,
            Err(e) => return fail(&e.to_string()),
        };
        match fs_roofline::record_run(&ledger, &axes, &post_axes, &results) {
            Ok(op) => {
                println!(
                    "{{\"ledgered\":true,\"citable\":{citable},\"op\":{op},\"db\":\"{}\"}}",
                    json_escape(&db)
                );
            }
            Err(e) => return fail(&e.to_string()),
        }
        for r in &results {
            match staleness(&ledger, &r.kernel, &r.version, axes.fingerprint) {
                Ok(s) => println!(
                    "{{\"kernel\":\"{}\",\"staleness\":\"{s:?}\"}}",
                    json_escape(&r.kernel)
                ),
                Err(e) => return fail(&e.to_string()),
            }
        }
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::json_escape;

    #[test]
    fn manual_json_fields_escape_hostile_paths_and_diagnostics() {
        assert_eq!(
            json_escape("ledger\\\"row\n\t\u{0001}.db"),
            "ledger\\\\\\\"row\\n\\t\\u0001.db"
        );
    }
}
