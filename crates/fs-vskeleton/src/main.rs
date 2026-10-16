//! CLI for the PV vertical skeleton:
//! `fs-vskeleton run <study.fsir> <ledger.db>` — execute + ledger + report.
//! `fs-vskeleton replay <ledger.db>` — re-execute and compare hashes.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("run") if args.len() == 4 => {
            let study = match std::fs::read_to_string(&args[2]) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("cannot read study {}: {e}", args[2]);
                    return ExitCode::FAILURE;
                }
            };
            match fs_vskeleton::run_study(&study, &args[3]) {
                Ok(outcome) => {
                    println!("{}", outcome.report);
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("run failed: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Some("replay") if args.len() == 3 => match fs_vskeleton::replay(&args[2]) {
            Ok(()) => {
                println!("replay OK: ledger reproduces bit-for-bit");
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("replay FAILED: {e}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("usage: fs-vskeleton run <study.fsir> <ledger.db> | replay <ledger.db>");
            ExitCode::FAILURE
        }
    }
}
