//! Command-line access to the existing normalized thermal marquee study.
#![deny(unsafe_code)]

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::process::ExitCode;

use fs_marquee::study::{
    IterRecord, MAX_ARMIJO_BACKTRACKS, PlateWithHoles, StudyConfig, StudyReport,
    StudyRunner, ThermalSource,
};

const MODEL: &str = "thermal-poisson-unit-plate-v1";
const MAX_HOLES: usize = 64;
const MAX_STEPS: usize = 256;
const MAX_BASE_CELLS: u64 = 65_536;
const MAX_ARGUMENT_BYTES: usize = 32_768;
const MAX_ARGUMENTS: usize = 256;
const HELP: &str = "fs-marquee: normalized thermal-compliance radius optimization\n\
Usage:\n\
  fs-marquee check [study options]\n\
  fs-marquee run [study options] [--expect-trace HEX]\n\
Required study options:\n\
  --units normalized --model-version 1\n\
  --hole X,Y,R                         repeat for each cooling hole (maximum 64)\n\
  --level N --steps N --step-size X     base-grid level 1..=8, steps 0..=256\n\
  --area X --r-min X --r-max X          material-area target and radius bounds\n\
  --max-base-cells N --max-evaluations N\n\
Optional study option:\n\
  --source CONSTANT,X_SLOPE,Y_SLOPE     default 1,0,0\n\
Output is deterministic JSONL. check admits/projects without solving a PDE.\n\
Budgets bound base-grid cells and worst-case solve-and-grade calls, not peak\n\
RAM or wall time. Each grade includes a state solve and DWR estimation.\n\
This is a deterministic, RNG-free, normalized Poisson model, not SI-calibrated\n\
elasticity. Error estimates are Estimated, not physical validation.\n\
Exit codes: 0 success, 2 input/admission, 3 solver, 4 I/O, 5 replay mismatch.\n";

type CliResult<T> = Result<T, CliError>;

#[derive(Debug)]
struct CliError {
    code: &'static str,
    message: String,
    exit: u8,
}

impl CliError {
    fn new(code: &'static str, message: impl Into<String>, exit: u8) -> Self {
        Self { code, message: message.into(), exit }
    }

    fn input(message: impl Into<String>) -> Self {
        Self::new("invalid_input", message, 2)
    }
}

impl From<io::Error> for CliError {
    fn from(error: io::Error) -> Self {
        Self::new("io_error", error.to_string(), 4)
    }
}

#[derive(Debug, Clone)]
struct Request {
    design: PlateWithHoles,
    config: StudyConfig,
    source: ThermalSource,
    max_base_cells: u64,
    max_evaluations: u64,
}

#[derive(Debug, Default)]
struct RunOptions {
    expected_trace: Option<String>,
}

fn finite_number(text: &str, name: &str) -> CliResult<f64> {
    let value = text.parse::<f64>().map_err(|_| CliError::input(format!("{name} must be a number")))?;
    if !value.is_finite() {
        return Err(CliError::input(format!("{name} must be finite")));
    }
    Ok(value)
}

fn unsigned_number(text: &str, name: &str) -> CliResult<u64> {
    text.parse::<u64>().map_err(|_| CliError::input(format!("{name} must be an unsigned integer")))
}

fn triple(text: &str, name: &str) -> CliResult<[f64; 3]> {
    let mut parts = text.split(',');
    let mut values = [0.0; 3];
    for value in &mut values {
        let part = parts.next().ok_or_else(|| CliError::input(format!("{name} requires three comma-separated numbers")))?;
        *value = finite_number(part, name)?;
    }
    if parts.next().is_some() {
        return Err(CliError::input(format!("{name} requires exactly three numbers")));
    }
    Ok(values)
}

fn valid_trace(text: &str) -> bool {
    text.len() == 64 && text.bytes().all(|byte| byte.is_ascii_hexdigit())
}

#[allow(clippy::too_many_lines)]
fn parse_request(args: &[String], allow_run_options: bool) -> CliResult<(Request, RunOptions)> {
    if args.len() > MAX_ARGUMENTS || args.iter().map(String::len).sum::<usize>() > MAX_ARGUMENT_BYTES {
        return Err(CliError::input("argument budget exceeded"));
    }
    let mut fields: BTreeMap<&str, &str> = BTreeMap::new();
    let mut design = PlateWithHoles { centers: Vec::new(), radii: Vec::new() };
    let mut options = RunOptions::default();
    let mut pairs = args.chunks_exact(2);
    for pair in &mut pairs {
        let (key, value) = (pair[0].as_str(), pair[1].as_str());
        if key == "--hole" {
            if design.radii.len() == MAX_HOLES {
                return Err(CliError::input("at most 64 holes are supported"));
            }
            let [x, y, radius] = triple(value, "--hole")?;
            design.centers.push([x, y]);
            design.radii.push(radius);
            continue;
        }
        match key {
            "--units" | "--model-version" | "--level" | "--steps" | "--step-size"
            | "--area" | "--r-min" | "--r-max" | "--max-base-cells"
            | "--max-evaluations" | "--source" => {},
            "--expect-trace" if allow_run_options => {},
            _ => return Err(CliError::input(format!("unknown or misplaced option: {key}"))),
        }
        if fields.insert(key, value).is_some() {
            return Err(CliError::input(format!("duplicate option: {key}")));
        }
    }
    if !pairs.remainder().is_empty() {
        return Err(CliError::input("every option requires a value"));
    }
    let required = |name: &str| -> CliResult<&str> {
        fields.get(name).copied().ok_or_else(|| CliError::input(format!("missing {name}")))
    };
    if required("--units")? != "normalized" || required("--model-version")? != "1" {
        return Err(CliError::input("only --units normalized --model-version 1 is supported"));
    }
    let level = unsigned_number(required("--level")?, "--level")?;
    let steps = unsigned_number(required("--steps")?, "--steps")?;
    if !(1..=8).contains(&level) || steps > MAX_STEPS as u64 {
        return Err(CliError::input("level must be in 1..=8 and steps in 0..=256"));
    }
    let [constant, x_slope, y_slope] = triple(fields.get("--source").copied().unwrap_or("1,0,0"), "--source")?;
    if let Some(expected) = fields.get("--expect-trace") {
        if !valid_trace(expected) {
            return Err(CliError::input("--expect-trace requires exactly 64 hexadecimal digits"));
        }
        options.expected_trace = Some(expected.to_ascii_lowercase());
    }
    let request = Request {
        design,
        config: StudyConfig {
            level: u32::try_from(level).map_err(|_| CliError::input("level overflow"))?,
            steps: usize::try_from(steps).map_err(|_| CliError::input("steps overflow"))?,
            step_size: finite_number(required("--step-size")?, "--step-size")?,
            area_target: finite_number(required("--area")?, "--area")?,
            r_min: finite_number(required("--r-min")?, "--r-min")?,
            r_max: finite_number(required("--r-max")?, "--r-max")?,
        },
        source: ThermalSource { constant, x_slope, y_slope },
        max_base_cells: unsigned_number(required("--max-base-cells")?, "--max-base-cells")?,
        max_evaluations: unsigned_number(required("--max-evaluations")?, "--max-evaluations")?,
    };
    Ok((request, options))
}

impl Request {
    fn required_evaluations(&self) -> u64 {
        let per_step = if self.config.step_size == 0.0 { 1 } else { MAX_ARMIJO_BACKTRACKS as u64 + 2 };
        self.config.steps as u64 * per_step
    }

    fn admit(&self) -> CliResult<StudyRunner> {
        let base_cells = 1_u64 << (2 * self.config.level);
        if self.max_base_cells > MAX_BASE_CELLS || base_cells > self.max_base_cells {
            return Err(CliError::new("cell_budget_exceeded", "base-grid cells exceed the declared budget or the 65536-cell CLI ceiling", 2));
        }
        if self.required_evaluations() > self.max_evaluations {
            return Err(CliError::new("evaluation_budget_exceeded", format!("the configured worst case requires {} solve-and-grade calls", self.required_evaluations()), 2));
        }
        StudyRunner::new_with_source(self.design.clone(), self.config.clone(), self.source)
            .map_err(|error| CliError::new("admission_failed", format!("{error:?}"), 2))
    }
}

fn json_string(text: &str) -> String {
    use std::fmt::Write as _;
    let mut result = String::with_capacity(text.len() + 2);
    result.push('"');
    for character in text.chars() {
        match character {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\r' => result.push_str("\\r"),
            '\t' => result.push_str("\\t"),
            control if control <= '\u{1f}' => { let _ = write!(result, "\\u{:04x}", u32::from(control)); },
            other => result.push(other),
        }
    }
    result.push('"');
    result
}

fn numbers(values: &[f64]) -> String {
    format!("[{}]", values.iter().map(|value| format!("{value:.17e}")).collect::<Vec<_>>().join(","))
}

fn geometry_json(design: &PlateWithHoles) -> String {
    let centers = design.centers.iter().map(|center| numbers(center)).collect::<Vec<_>>().join(",");
    format!("{{\"centers\":[{centers}],\"radii\":{},\"material_area\":{:.17e}}}", numbers(&design.radii), design.area())
}

fn admission_row(request: &Request, runner: &StudyRunner) -> String {
    format!(
        "{{\"schema\":1,\"event\":\"admitted\",\"model\":{},\"version\":{},\"units\":\"normalized\",\"rng\":\"none\",\"capabilities\":[\"thermal-poisson-radius-study\"],\"base_cells\":{},\"max_base_cells\":{},\"worst_case_evaluations\":{},\"max_evaluations\":{},\"steps\":{},\"source\":{},\"initial_projected_design\":{}}}\n",
        json_string(MODEL), json_string(fs_marquee::VERSION), 1_u64 << (2 * request.config.level),
        request.max_base_cells, request.required_evaluations(), request.max_evaluations, request.config.steps,
        numbers(&[request.source.constant, request.source.x_slope, request.source.y_slope]), geometry_json(runner.design()),
    )
}

fn iteration_row(record: &IterRecord) -> String {
    format!(
        "{{\"schema\":1,\"event\":\"iteration\",\"record\":{},\"radii\":{},\"accepted_radii\":{}}}\n",
        record.jsonl_row().trim_end(), numbers(&record.radii), numbers(&record.accepted_radii),
    )
}

fn result_row(report: &StudyReport) -> String {
    let objective = report.iterations.last().map_or_else(|| "null".to_string(), |record| format!("{:.17e}", record.accepted_compliance));
    format!(
        "{{\"schema\":1,\"event\":\"complete\",\"model\":{},\"trace_hash\":{},\"iterations\":{},\"accepted_compliance\":{objective},\"evidence\":\"estimated\",\"final_design\":{}}}\n",
        json_string(MODEL), json_string(&report.trace_hash), report.iterations.len(), geometry_json(&report.design),
    )
}

fn emit(output: &mut impl Write, text: &str) -> CliResult<()> {
    output.write_all(text.as_bytes())?;
    output.flush()?;
    Ok(())
}

fn run_study_cli(request: &Request, options: &RunOptions, output: &mut impl Write) -> CliResult<()> {
    let mut runner = request.admit()?;
    emit(output, &admission_row(request, &runner))?;
    while runner.advance().map_err(|error| CliError::new("solver_failed", format!("{error:?}"), 3))? {
        if let Some(record) = runner.iterations().last() {
            emit(output, &iteration_row(record))?;
        }
    }
    let report = runner.report();
    if let Some(expected) = &options.expected_trace {
        if *expected != report.trace_hash {
            return Err(CliError::new("trace_mismatch", format!("expected {expected}, obtained {}", report.trace_hash), 5));
        }
    }
    emit(output, &result_row(&report))
}

fn execute(args: &[String], output: &mut impl Write) -> CliResult<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(CliError::input("missing command; use --help"));
    };
    match command {
        "--help" | "help" if args.len() == 1 => emit(output, HELP),
        "--version" if args.len() == 1 => emit(output, &format!("fs-marquee {}\n", fs_marquee::VERSION)),
        "check" => {
            let (request, _) = parse_request(&args[1..], false)?;
            let runner = request.admit()?;
            emit(output, &admission_row(&request, &runner))
        },
        "run" => {
            let (request, options) = parse_request(&args[1..], true)?;
            run_study_cli(&request, &options, output)
        },
        _ => Err(CliError::input("unknown command; use --help")),
    }
}

fn arguments() -> CliResult<Vec<String>> {
    let mut args = Vec::new();
    let mut bytes = 0_usize;
    for arg in std::env::args_os().skip(1) {
        let text = arg.into_string().map_err(|_| CliError::input("arguments must be UTF-8"))?;
        bytes = bytes.checked_add(text.len()).ok_or_else(|| CliError::input("argument size overflow"))?;
        if bytes > MAX_ARGUMENT_BYTES || args.len() == MAX_ARGUMENTS {
            return Err(CliError::input("argument budget exceeded"));
        }
        args.push(text);
    }
    Ok(args)
}

fn main() -> ExitCode {
    let result = arguments().and_then(|args| execute(&args, &mut io::stdout().lock()));
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let _ = writeln!(io::stderr().lock(), "{{\"schema\":1,\"event\":\"error\",\"code\":{},\"message\":{}}}", json_string(error.code), json_string(&error.message));
            ExitCode::from(error.exit)
        },
    }
}
