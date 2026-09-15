use super::*;
use std::io::{Read, Write};
use std::process::{Command, Stdio};
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

pub(super) fn evaluate_sample(request: &str, deadline: Instant) -> std::result::Result<f64, EvaluationError> {
    if Instant::now() >= deadline { return Err(EvaluationError::Budget); }
    let executable = std::env::current_exe().map_err(|error| EvaluationError::Child(format!("cannot locate frankensim executable: {error}")))?;
    let mut child = Command::new(executable)
        .arg("--json").arg("cooling-network").arg("/dev/stdin")
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
        .spawn().map_err(|error| EvaluationError::Child(format!("cannot launch cooling sample: {error}")))?;
    {
        let mut stdin = child.stdin.take().ok_or_else(|| EvaluationError::Child("cooling child stdin unavailable".into()))?;
        stdin.write_all(request.as_bytes()).map_err(|error| EvaluationError::Child(format!("cannot send cooling sample: {error}")))?;
    }
    let stdout = child.stdout.take().ok_or_else(|| EvaluationError::Child("cooling child stdout unavailable".into()))?;
    let stderr = child.stderr.take().ok_or_else(|| EvaluationError::Child("cooling child stderr unavailable".into()))?;
    let stdout_reader = std::thread::spawn(move || drain(stdout, MAX_CHILD_OUTPUT_BYTES));
    let stderr_reader = std::thread::spawn(move || drain(stderr, MAX_CHILD_OUTPUT_BYTES));
    let status = loop {
        match child.try_wait() {
            Ok(Some(status)) => break status,
            Ok(None) => {
                if Instant::now() >= deadline {
                    let _ = child.kill(); let _ = child.wait();
                    let _ = stdout_reader.join(); let _ = stderr_reader.join();
                    return Err(EvaluationError::Budget);
                }
                std::thread::sleep(Duration::from_millis(5));
            }
            Err(error) => {
                let _ = child.kill(); let _ = child.wait();
                let _ = stdout_reader.join(); let _ = stderr_reader.join();
                return Err(EvaluationError::Child(format!("cannot poll cooling child: {error}")));
            }
        }
    };
    let stdout = stdout_reader.join().map_err(|_| EvaluationError::Child("cooling stdout reader panicked".into()))?.map_err(EvaluationError::Child)?;
    let stderr = stderr_reader.join().map_err(|_| EvaluationError::Child("cooling stderr reader panicked".into()))?.map_err(EvaluationError::Child)?;
    if !status.success() {
        return Err(EvaluationError::Child(format!("cooling sample exited with {status}: {}", String::from_utf8_lossy(&stderr).trim())));
    }
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
