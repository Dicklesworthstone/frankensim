use super::*;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

const MAX_CHILD_OUTPUT_BYTES: usize = 32 * 1024 * 1024;

#[derive(Debug)]
pub(super) enum EvaluationError { Budget, Child(String) }
impl fmt::Display for EvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Budget => f.write_str("overall UQ wall-time budget expired during a cooling sample"),
            Self::Child(message) => write!(f, "{message}"),
        }
    }
}

/// On Linux this names the running image even if an installation replaces its
/// pathname during a long study. Hashing and launch use the same image.
pub(super) fn executable_path() -> std::io::Result<PathBuf> {
    #[cfg(target_os = "linux")]
    { Ok(PathBuf::from("/proc/self/exe")) }
    #[cfg(not(target_os = "linux"))]
    { std::env::current_exe() }
}

struct RunningChild(Child);
impl Drop for RunningChild {
    fn drop(&mut self) {
        // Best-effort cleanup also covers early pipe/write/read errors. A
        // reaped child simply returns an error from kill; never mask the
        // original numerical or I/O diagnosis with a cleanup error.
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

pub(super) fn evaluate_sample(request: &str, deadline: Instant) -> std::result::Result<f64, EvaluationError> {
    let executable = executable_path().map_err(|error| EvaluationError::Child(format!("cannot locate frankensim executable: {error}")))?;
    let mut command = Command::new(executable);
    command.arg("--json").arg("cooling-network").arg("/dev/stdin");
    evaluate_process(command, request, deadline)
}

fn evaluate_process(mut command: Command, request: &str, deadline: Instant) -> std::result::Result<f64, EvaluationError> {
    if Instant::now() >= deadline { return Err(EvaluationError::Budget); }
    let mut child = RunningChild(command
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
        .spawn().map_err(|error| EvaluationError::Child(format!("cannot launch cooling sample: {error}")))?);
    let mut stdin = child.0.stdin.take().ok_or_else(|| EvaluationError::Child("cooling child stdin unavailable".into()))?;
    let stdout = child.0.stdout.take().ok_or_else(|| EvaluationError::Child("cooling child stdout unavailable".into()))?;
    let stderr = child.0.stderr.take().ok_or_else(|| EvaluationError::Child("cooling child stderr unavailable".into()))?;
    let stdout_reader = std::thread::spawn(move || drain(stdout, MAX_CHILD_OUTPUT_BYTES));
    let stderr_reader = std::thread::spawn(move || drain(stderr, MAX_CHILD_OUTPUT_BYTES));
    // A request can be much larger than a pipe. Never write it on the watchdog
    // thread: a blocked write would make the declared timeout unenforceable.
    let request = request.as_bytes().to_vec();
    let input_writer = std::thread::spawn(move || stdin.write_all(&request));
    let status = loop {
        match child.0.try_wait() {
            Ok(Some(status)) => break Ok(status),
            Ok(None) => {
                if Instant::now() >= deadline { break Err(EvaluationError::Budget); }
                std::thread::sleep(Duration::from_millis(5));
            }
            Err(error) => break Err(EvaluationError::Child(format!("cannot poll cooling child: {error}"))),
        }
    };
    // Termination closes all pipes in this same-binary producer; join every
    // I/O thread only AFTER killing/reaping on timeout or polling failure.
    drop(child);
    let input = input_writer.join();
    let stdout = stdout_reader.join();
    let stderr = stderr_reader.join();
    let status = status?;
    let stdout = stdout.map_err(|_| EvaluationError::Child("cooling stdout reader panicked".into()))?.map_err(EvaluationError::Child)?;
    let stderr = stderr.map_err(|_| EvaluationError::Child("cooling stderr reader panicked".into()))?.map_err(EvaluationError::Child)?;
    if !status.success() {
        return Err(EvaluationError::Child(format!("cooling sample exited with {status}: {}", String::from_utf8_lossy(&stderr).trim())));
    }
    input.map_err(|_| EvaluationError::Child("cooling stdin writer panicked".into()))?
        .map_err(|error| EvaluationError::Child(format!("cannot send cooling sample: {error}")))?;
    let text = std::str::from_utf8(&stdout).map_err(|_| EvaluationError::Child("cooling sample output is not UTF-8".into()))?;
    let document = J::parse(text).map_err(|error| EvaluationError::Child(format!("cooling sample emitted invalid JSON: {error}")))?;
    document.path(&["objective", "value_k"]).and_then(J::as_f64).filter(|value| value.is_finite())
        .ok_or_else(|| EvaluationError::Child("cooling sample has no finite objective.value_k".into()))
}

fn drain(mut reader: impl Read, cap: usize) -> std::result::Result<Vec<u8>, String> {
    let mut kept = Vec::new();
    let mut total = 0usize;
    let mut buffer = [0u8; 8192];
    loop {
        let count = reader.read(&mut buffer).map_err(|error| format!("cannot read cooling child output: {error}"))?;
        if count == 0 { break; }
        total = total.checked_add(count).ok_or_else(|| "cooling child output length overflowed".to_string())?;
        if kept.len() < cap + 1 {
            let room = cap + 1 - kept.len();
            kept.extend_from_slice(&buffer[..count.min(room)]);
        }
    }
    if total > cap { Err(format!("cooling child output exceeds {cap} bytes")) } else { Ok(kept) }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;

    #[test]
    fn a_child_that_never_reads_a_large_request_cannot_block_the_watchdog() {
        let mut command = Command::new("/bin/sh");
        // exec avoids a grandchild inheriting the pipes after the child dies.
        command.arg("-c").arg("exec sleep 2");
        let request = "x".repeat(2 * 1024 * 1024);
        let result = evaluate_process(command, &request, Instant::now() + Duration::from_millis(30));
        assert!(matches!(result, Err(EvaluationError::Budget)));
    }

    #[test]
    fn a_real_child_refusal_stays_a_model_failure() {
        let mut command = Command::new("/bin/sh");
        command.arg("-c").arg("printf 'solver-refusal-sentinel' >&2; exit 7");
        let result = evaluate_process(command, "", Instant::now() + Duration::from_secs(2));
        match result {
            Err(EvaluationError::Child(message)) => assert!(message.contains("solver-refusal-sentinel")),
            other => panic!("expected model refusal, got {other:?}"),
        }
    }
}
