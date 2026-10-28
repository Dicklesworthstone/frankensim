//! `frankensim` command-line entry point.

use std::process::ExitCode;

fn main() -> ExitCode {
    // On Linux, retain the running image even if an installation replaces its
    // pathname while a checkpointed study is running.
    #[cfg(target_os = "linux")]
    let worker = Some(std::path::PathBuf::from("/proc/self/exe"));
    #[cfg(not(target_os = "linux"))]
    let worker = std::env::current_exe().ok();
    let output = fs_cli::run_os_with_context(
        std::env::args_os().skip(1),
        fs_cli::CommandContext { cooling_worker: worker.as_deref() },
    );
    print!("{}", output.stdout);
    eprint!("{}", output.stderr);
    ExitCode::from(output.exit_code)
}
