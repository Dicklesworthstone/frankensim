//! Bounded empirical Monte Carlo through the actual `cooling-network` product.
//!
//! Every sample is serialized as a complete cooling request and executed by a
//! child invocation of this same `frankensim` binary. This deliberately pays
//! process-launch overhead to keep one parser/physics/product boundary instead
//! of maintaining a second cooling implementation inside UQ.

#[path = "json_read.rs"]
mod json;
mod child;
mod execute;
mod model;

use fs_cli::{CommandOutput, exit};
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan, UqStatus,
};
use json::JsonValue as J;
use std::ffi::{OsStr, OsString};
use std::fmt;
use std::fs::File;
use std::io::Read;

const SCHEMA: &str = "frankensim.cooling-network-uq.v1";
const RESULT_SCHEMA: &str = "frankensim.cooling-network-uq.result.v1";
const MAX_BASE_BYTES: u64 = 16 * 1024 * 1024;
const MAX_UQ_BYTES: u64 = 4 * 1024 * 1024;
const MAX_PRODUCT_SAMPLES: usize = 10_000;
const HELP: &str = "Usage: frankensim [--json] cooling-network-uq <base-request.json> <uq-request.json>\n\nRun fixed-count empirical Monte Carlo by invoking the actual steady cooling-network\nproducer once per sample. See examples/cooling-network/COOLING_UQ.md.\n";

#[derive(Debug)]
struct Failure {
    code: &'static str,
    message: String,
}
impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}
impl std::error::Error for Failure {}
type Result<T> = std::result::Result<T, Failure>;

fn bad(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-input", message: message.into() }
}
fn budget(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-budget", message: message.into() }
}
fn model_failure(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-model", message: message.into() }
}

pub(super) fn run(args: &[OsString], json_mode: bool) -> CommandOutput {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return CommandOutput {
            exit_code: exit::SUCCESS,
            stdout: if json_mode {
                format!("{{\"schema\":{},\"help\":{}}}\n", quote(RESULT_SCHEMA), quote(HELP))
            } else {
                HELP.into()
            },
            stderr: String::new(),
        };
    }
    if args.len() != 2 {
        return diagnostic(exit::USAGE, bad(HELP), json_mode);
    }
    let base = match read(&args[0], MAX_BASE_BYTES) {
        Ok(value) => value,
        Err(error) => return diagnostic(exit::INPUT, error, json_mode),
    };
    let uq = match read(&args[1], MAX_UQ_BYTES) {
        Ok(value) => value,
        Err(error) => return diagnostic(exit::INPUT, error, json_mode),
    };
    match execute::execute(&base, &uq) {
        Ok(stdout) => CommandOutput { exit_code: exit::SUCCESS, stdout, stderr: String::new() },
        Err(error) => {
            let class = if error.code == "cooling-network-uq-budget" {
                exit::BUDGET
            } else {
                exit::REFUSED
            };
            diagnostic(class, error, json_mode)
        }
    }
}

fn diagnostic(code: u8, failure: Failure, json_mode: bool) -> CommandOutput {
    let stderr = if json_mode {
        format!(
            "{{\"schema\":\"frankensim.cooling-network-uq.diagnostic.v1\",\"code\":{},\"message\":{}}}\n",
            quote(failure.code),
            quote(&failure.message)
        )
    } else {
        format!("{failure}\n")
    };
    CommandOutput { exit_code: code, stdout: String::new(), stderr }
}

fn read(path: &OsStr, cap: u64) -> Result<String> {
    let mut text = String::new();
    File::open(std::path::Path::new(path))
        .map_err(|error| bad(error.to_string()))?
        .take(cap + 1)
        .read_to_string(&mut text)
        .map_err(|error| bad(error.to_string()))?;
    if text.len() as u64 > cap {
        return Err(bad("input exceeds admitted byte cap"));
    }
    Ok(text)
}

fn object<'a>(value: &'a J, allowed: &[&str], path: &str) -> Result<&'a J> {
    let members = value
        .as_object()
        .ok_or_else(|| bad(format!("{path} must be an object")))?;
    for (name, _) in members {
        if !allowed.contains(&name.as_str()) {
            return Err(bad(format!("unknown field {path}.{name}")));
        }
    }
    Ok(value)
}
fn field<'a>(value: &'a J, key: &str) -> Result<&'a J> {
    value.get(key).ok_or_else(|| bad(format!("missing required field {key}")))
}
fn string(value: &J, path: &str) -> Result<String> {
    let text = value.as_str().ok_or_else(|| bad(format!("{path} must be a string")))?;
    if text.is_empty() || text.trim() != text || text.len() > 256 || text.chars().any(char::is_control) {
        return Err(bad(format!("{path} must be a nonempty trimmed control-free string <=256 bytes")));
    }
    Ok(text.into())
}
fn number(value: &J, path: &str) -> Result<f64> {
    value
        .as_f64()
        .filter(|number| number.is_finite())
        .ok_or_else(|| bad(format!("{path} must be a finite number")))
}
fn positive(value: &J, path: &str) -> Result<f64> {
    let value = number(value, path)?;
    if value > 0.0 { Ok(value) } else { Err(bad(format!("{path} must be positive"))) }
}
fn integer_raw(value: &J, path: &str) -> Result<usize> {
    value
        .number_raw()
        .and_then(|raw| raw.parse::<usize>().ok())
        .ok_or_else(|| bad(format!("{path} must use a nonnegative integer JSON spelling")))
}
fn integer(value: &J, path: &str, maximum: usize) -> Result<usize> {
    let value = integer_raw(value, path)?;
    if value <= maximum { Ok(value) } else { Err(bad(format!("{path} exceeds {maximum}"))) }
}
fn array<'a>(value: &'a J, path: &str, maximum: usize) -> Result<&'a [J]> {
    let values = value.as_array().ok_or_else(|| bad(format!("{path} must be an array")))?;
    if values.len() > maximum {
        Err(bad(format!("{path} has more than {maximum} entries")))
    } else {
        Ok(values)
    }
}
fn quote(text: &str) -> String {
    use std::fmt::Write as _;
    let mut output = String::from("\"");
    for ch in text.chars() {
        match ch {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            c if c < '\u{20}' => { let _ = write!(output, "\\u{:04x}", c as u32); }
            c => output.push(c),
        }
    }
    output.push('"');
    output
}
fn number_json(value: f64) -> Result<String> {
    if value.is_finite() { Ok(value.to_string()) } else { Err(model_failure("nonfinite result cannot be published")) }
}
fn optional_number(value: Option<f64>) -> Result<String> {
    value.map_or_else(|| Ok("null".into()), number_json)
}
