//! Bounded Monte Carlo through the actual `cooling-network` product, with
//! steady or completed-trajectory observables and optional probability-target
//! stopping using the existing UQ confidence sequence.
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
pub(super) mod component_design;

use fs_cli::{CommandOutput, exit};
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan, UqStatus,
};
use json::JsonValue as J;
use std::ffi::{OsStr, OsString};
use std::fmt;
use std::fs::File;
use std::io::Read as _;

const SCHEMA: &str = "frankensim.cooling-network-uq.v1";
const RESULT_SCHEMA: &str = "frankensim.cooling-network-uq.result.v1";
const MAX_BASE_BYTES: u64 = 16 * 1024 * 1024;
const MAX_UQ_BYTES: u64 = 4 * 1024 * 1024;
const MAX_PRODUCT_SAMPLES: usize = 10_000;
const HELP: &str = "Usage: frankensim [--json] cooling-network-uq <base-request.json> <uq-request.json> [options]\n\nRun Monte Carlo through actual cooling-network solves. Default steady mode uses\nthe declared objective. Transient bases require qoi.kind=transient-sampled-peak\nin the UQ document: one observation is the initial/accepted-endpoint peak of one\nCOMPLETE trajectory, including all fixed repeated cycles, not its final value.\nActive schedule loads, fan speeds, initial temperature and heat capacities may\nbe sampled. Nested design searches and variable periodic horizons are refused.\n\nRecovery options:\n--checkpoint NEW-PATH writes a new model-bound checkpoint and updates it atomically\nafter every completed sample. Existing destinations are refused, never overwritten.\n--resume SAVED-PATH restores a trusted checkpoint under the identical base request,\nplan, executable and stopping policy. Use a fresh --checkpoint for further progress.\n--max-new-samples N limits this invocation (zero allowed) and requires --checkpoint.\nThe request's wall_seconds is a fresh evaluation-time allowance per invocation;\nthe original lifetime sample count never resets. Interrupted samples retry the\nsame ordinal; an unfinished transient restarts from its declared initial state.\nFixed-count partial runs emit progress, not a distribution.\n\nRandomized quasi-Monte Carlo (fixed layout):\n--qmc-replicates R uses R independently keyed Owen-scrambled Sobol nets.\nR must be in 2..256; samples/R must be a power of two >=2, with no remainder.\nAt most ten uncertain parameters are supported. Checkpoint/resume and chunk\nlimits retain complete evaluations, including an unfinished net; repeat the\nidentical --qmc-replicates on resume. Partial runs publish progress only.\nOnly the complete fixed layout publishes between-replicate standard errors.\nSequential compliance and candidate-selection flags are unavailable for QMC.\nSee examples/cooling-network/QMC_UQ.md.\n\nDirect-model variance sensitivity:\n--sensitivity sobol estimates main and total input effects through actual solves.\nFor d independent varying inputs, samples must equal N*(d+2), N >= 2.\nCheckpoint/resume retains every A/B/hybrid evaluation, including partial rows.\nRepeat --sensitivity sobol on resume; --max-new-samples counts physical solves.\nOnly the complete fixed design publishes indices. Partial runs publish progress.\nQMC, sequential compliance and candidate-selection flags are incompatible.\nSee examples/cooling-network/SENSITIVITY.md.\n\nSequential compliance options (supply ALL three together):\n--compliance-probability P declares the required probability in (0,1).\n--confidence-alpha A declares the confidence-sequence error level in (0,1).\n--min-decision-samples N requires at least N observations before a decision.\nThe UQ request must declare temperature_limit_k. Each accepted cooling solve\nupdates the existing Bernoulli-indicator Gaussian-mixture confidence sequence.\nStop when its interval resolves P(QoI <= limit) >= P; samples is the lifetime cap.\nThe compliance.v1 result distinguishes meets-probability-target, below-probability-\ntarget and indeterminate. Either resolved decision exits SUCCESS; an unresolved\nsample/time/chunk budget exits BUDGET. This is inference about the declared\nnumerical model, NOT physical safety or continuous-time peak certification.\nWithout these flags, the fixed-count result.v1 behavior is retained.\n\nUncertainty-aware design (requires the complete compliance policy above):\n--fan-speed-candidates 0.75,1,1.25 selects the lowest qualified fan multiplier.\n--power-candidates 0.5,0.75,1 selects the highest qualified workload multiplier.\nChoose one control and 2..64 strictly increasing predeclared candidates. Fan\nmultipliers must be positive; workload multipliers can include zero. Each\ncandidate scales the actual sampled input AFTER drawing uncertainty. Transient\ncontrols scale the complete declared schedule, including every repeated cycle.\nSteady workload selection requires solid.component_power. The design.v1 result\nselects only when all more-preferred candidates are below target; unresolved\npreferred candidates give inconclusive, even if qualified_multiplier is present.\nNo interpolation, monotonicity or continuous optimum is inferred. Alpha is a\nFAMILY error budget divided across all predeclared candidates; common draws do\nnot require independence between candidates. samples is the cap PER candidate;\ncandidate count times samples must be <=10000. The one wall deadline and\n--max-new-samples apply across the whole invocation. Recovery retains all\ncandidate prefixes; changing the family, control or policy on resume refuses.\nA resolved no-qualified-candidate exits SUCCESS as a calculation, not approval.\nSee examples/cooling-network/UNCERTAINTY_AWARE_DESIGN.md.\nSee examples/cooling-network/COOLING_UQ.md and TRANSIENT_UQ.md.\n";

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
    if args.len() < 2 {
        return diagnostic(exit::USAGE, bad(HELP), json_mode);
    }
    let options = match execute::Options::parse(&args[2..]) {
        Ok(options) => options,
        Err(error) => return diagnostic(exit::USAGE, error, json_mode),
    };
    let base = match read(&args[0], MAX_BASE_BYTES) {
        Ok(value) => value,
        Err(error) => return diagnostic(exit::INPUT, error, json_mode),
    };
    let uq = match read(&args[1], MAX_UQ_BYTES) {
        Ok(value) => value,
        Err(error) => return diagnostic(exit::INPUT, error, json_mode),
    };
    let result = if args.len() == 2 {
        execute::execute(&base, &uq).map(|stdout| execute::ExecutionOutput {
            stdout, exit_code: exit::SUCCESS,
        })
    } else {
        execute::execute_with_options(&base, &uq, &options)
    };
    match result {
        Ok(output) => CommandOutput { exit_code: output.exit_code, stdout: output.stdout, stderr: String::new() },
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
