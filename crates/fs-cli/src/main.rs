//! `frankensim` command-line entry point.

mod network_command;
mod uq_command;

use std::ffi::OsString;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<OsString> = std::env::args_os().skip(1).collect();
    let json = args.first().is_some_and(|arg| arg == "--json");
    let command = usize::from(json);
    let mut output = if args.get(command).is_some_and(|arg| arg == "cooling-network") {
        network_command::run(&args[command + 1..], json)
    } else if args.get(command).is_some_and(|arg| arg == "cooling-network-uq") {
        uq_command::run(&args[command + 1..], json)
    } else {
        fs_cli::run_os(args.clone())
    };
    if !json && output.exit_code == 0
        && (args.is_empty() || (args.len() == 1 && matches!(args[0].to_str(), Some("--help" | "-h" | "help"))))
    {
        output.stdout.push_str("\nExperimental coupled workflows:\n  cooling-network <request.json>\n  cooling-network-uq <base-request.json> <uq-request.json>\nUse each command's --help for its file-driven interface.\n");
    }
    print!("{}", output.stdout);
    eprint!("{}", output.stderr);
    ExitCode::from(output.exit_code)
}
