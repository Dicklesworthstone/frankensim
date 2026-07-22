//! Command-line entry point for the isolated high-precision oracle lane.

use std::process::ExitCode;

use frankensim_high_precision_oracle::{AuditConfig, run_audit};

fn parse_value<T: std::str::FromStr>(flag: &str, value: Option<String>) -> Result<T, String> {
    let raw = value.ok_or_else(|| format!("{flag} requires a value"))?;
    raw.parse()
        .map_err(|_| format!("{flag} has an invalid value: {raw}"))
}

fn parse_config() -> Result<Option<AuditConfig>, String> {
    let mut config = AuditConfig::default();
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--samples" => config.samples = parse_value("--samples", args.next())?,
            "--precision-bits" => {
                config.precision_bits = parse_value("--precision-bits", args.next())?;
            }
            "-h" | "--help" => return Ok(None),
            _ => return Err(format!("unknown argument: {argument}")),
        }
    }
    config.validate().map(Some)
}

fn main() -> ExitCode {
    let config = match parse_config() {
        Ok(Some(config)) => config,
        Ok(None) => {
            println!("usage: frankensim-high-precision-oracle [--samples N] [--precision-bits N]");
            return ExitCode::SUCCESS;
        }
        Err(error) => {
            eprintln!("oracle admission refused: {error}");
            return ExitCode::FAILURE;
        }
    };
    let report = match run_audit(config) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("oracle execution refused: {error}");
            return ExitCode::FAILURE;
        }
    };
    for line in report.render_json_lines() {
        println!("{line}");
    }
    if report.passed() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
