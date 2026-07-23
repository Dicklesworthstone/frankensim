//! Deterministic CLI adapter for the complete wedge decision audit.

use fs_govern::{WedgeAuditLog, WedgeDecisionAuditRequest, build_wedge_decision_audit};
use std::process::ExitCode;

const USAGE: &str = "usage: wedge-audit --measured-days <positive-days> --cycle-time-evidence <locator> [--seed-fault missing-cycle-time-evidence]";

#[derive(Debug, Clone, PartialEq)]
struct RunArgs {
    measured_days: f64,
    cycle_time_evidence: String,
    seed_missing_evidence: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CliRefusal {
    code: &'static str,
    detail: String,
    fix: &'static str,
}

enum Action {
    Help,
    Run(RunArgs),
}

fn next_value(
    args: &mut impl Iterator<Item = String>,
    flag: &'static str,
) -> Result<String, CliRefusal> {
    args.next().ok_or_else(|| CliRefusal {
        code: "missing-flag-value",
        detail: format!("{flag} requires a value"),
        fix: USAGE,
    })
}

fn parse_args() -> Result<Action, CliRefusal> {
    let mut measured_days = None;
    let mut cycle_time_evidence = None;
    let mut seed_missing_evidence = false;
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--help" | "-h" => return Ok(Action::Help),
            "--measured-days" => {
                if measured_days.is_some() {
                    return Err(CliRefusal {
                        code: "duplicate-flag",
                        detail: "--measured-days was supplied more than once".to_string(),
                        fix: USAGE,
                    });
                }
                let value = next_value(&mut args, "--measured-days")?;
                measured_days = Some(value.parse::<f64>().map_err(|error| CliRefusal {
                    code: "invalid-measured-days",
                    detail: format!("--measured-days value `{value}` is not a number: {error}"),
                    fix: "supply one positive finite working-day value",
                })?);
            }
            "--cycle-time-evidence" => {
                if cycle_time_evidence.is_some() {
                    return Err(CliRefusal {
                        code: "duplicate-flag",
                        detail: "--cycle-time-evidence was supplied more than once".to_string(),
                        fix: USAGE,
                    });
                }
                cycle_time_evidence = Some(next_value(&mut args, "--cycle-time-evidence")?);
            }
            "--seed-fault" => {
                if seed_missing_evidence {
                    return Err(CliRefusal {
                        code: "duplicate-flag",
                        detail: "--seed-fault was supplied more than once".to_string(),
                        fix: USAGE,
                    });
                }
                let fault = next_value(&mut args, "--seed-fault")?;
                if fault != "missing-cycle-time-evidence" {
                    return Err(CliRefusal {
                        code: "unknown-seeded-fault",
                        detail: format!("unsupported seeded fault `{fault}`"),
                        fix: "use --seed-fault missing-cycle-time-evidence",
                    });
                }
                seed_missing_evidence = true;
            }
            unknown => {
                return Err(CliRefusal {
                    code: "unknown-argument",
                    detail: format!("unknown argument `{unknown}`"),
                    fix: USAGE,
                });
            }
        }
    }
    Ok(Action::Run(RunArgs {
        measured_days: measured_days.ok_or_else(|| CliRefusal {
            code: "missing-required-flag",
            detail: "--measured-days is required".to_string(),
            fix: USAGE,
        })?,
        cycle_time_evidence: cycle_time_evidence.ok_or_else(|| CliRefusal {
            code: "missing-required-flag",
            detail: "--cycle-time-evidence is required".to_string(),
            fix: USAGE,
        })?,
        seed_missing_evidence,
    }))
}

fn warn(code: &'static str, detail: impl Into<String>, fix: &'static str) {
    eprintln!(
        "{}",
        WedgeAuditLog::warning("wedge-audit-refusal", code, detail, fix).to_json()
    );
}

fn main() -> ExitCode {
    let action = match parse_args() {
        Ok(action) => action,
        Err(refusal) => {
            warn(refusal.code, refusal.detail, refusal.fix);
            return ExitCode::from(2);
        }
    };
    let Action::Run(mut args) = action else {
        println!("{USAGE}");
        return ExitCode::SUCCESS;
    };
    if args.seed_missing_evidence {
        args.cycle_time_evidence.clear();
    }
    let request = match WedgeDecisionAuditRequest::new(args.measured_days, args.cycle_time_evidence)
    {
        Ok(request) => request,
        Err(error) => {
            warn(error.code(), error.to_string(), error.fix());
            return ExitCode::from(3);
        }
    };
    let audit = match build_wedge_decision_audit(&request) {
        Ok(audit) => audit,
        Err(error) => {
            warn(error.code(), error.to_string(), error.fix());
            return ExitCode::from(3);
        }
    };
    for log in audit.logs() {
        eprintln!("{}", log.to_json());
    }
    println!("{}", audit.artifact());
    ExitCode::SUCCESS
}
