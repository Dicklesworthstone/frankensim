//! Deterministic command-line contract for FrankenSim's Cooling 0.1 product
//! workflow.
//!
//! The v0 surface exposes `validate`, one-source `import`, `solve`, `report`,
//! and `package` while keeping authority honest. Project validation delegates
//! to the strict [`fs_project`] readers. Product stages whose producing Beads
//! remain open fail before side effects with `cli-stage-unavailable`; a
//! CLI-shaped mock is not substituted for an integrated workflow.

mod import;

use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};

pub use import::{
    GeometryImportLimits, GeometryImportRefusal, GeometryImportRun, RawGeometryLibrary,
    RecordedImportRefusal, RetainedGeometryImport, import_project_geometry,
};

/// Maximum project source accepted by the CLI.
pub const MAX_PROJECT_BYTES: u64 = 16 * 1024 * 1024;

/// Stable process exit classes.
pub mod exit {
    /// Command completed successfully.
    pub const SUCCESS: u8 = 0;
    /// Arguments did not match the documented grammar.
    pub const USAGE: u8 = 2;
    /// Input could not be read, decoded, or admitted by the CLI resource cap.
    pub const INPUT: u8 = 3;
    /// The project was read but refused by its schema or semantic validator.
    pub const REFUSED: u8 = 4;
    /// The command is reserved but its authoritative producer is not shipped.
    pub const UNAVAILABLE: u8 = 5;
}

const RESULT_SCHEMA: &str = "frankensim.cli.result.v1";
const DIAGNOSTIC_SCHEMA: &str = "frankensim.cli.diagnostic.v1";
const VALIDATION_AUTHORITY: &str = "structural-project-admission";
const VALIDATION_NO_CLAIM: &str =
    "does not prove artifact existence, capability availability, solvability, or physical validity";
const USAGE: &str = "frankensim [--json] validate <project.fsim|project.json> | import <project> <source> <ledger.db> --unit <unit> (--max-hole-edges <n> | --step-root <id> --target-h <spacing>) | solve <project> | solve --resume <run-id> | report <run-id> | package <run-id>";

/// Captured command output. Final result records are on stdout; diagnostics
/// are on stderr.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommandOutput {
    /// Stable process exit code.
    pub exit_code: u8,
    /// Final result records only.
    pub stdout: String,
    /// Diagnostics, one JSON object per line in JSON mode.
    pub stderr: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputMode {
    Text,
    Json,
}

#[derive(Debug, Clone, PartialEq)]
enum Command {
    Help,
    Validate(PathBuf),
    Import(ImportCommand),
    SolveProject(PathBuf),
    Resume(String),
    Report(String),
    Package(String),
}

#[derive(Debug, Clone, PartialEq)]
struct ImportCommand {
    project: PathBuf,
    source: PathBuf,
    ledger: PathBuf,
    unit: String,
    policy: ImportPolicy,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum ImportPolicy {
    Mesh { max_hole_edges: usize },
    FacetedStep { root_id: u64, target_h: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProjectSyntax {
    Sexpr,
    Json,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Diagnostic {
    command: &'static str,
    code: String,
    message: String,
    fix: String,
    subject: Option<String>,
}

impl Diagnostic {
    fn new(
        command: &'static str,
        code: impl Into<String>,
        message: impl Into<String>,
        fix: impl Into<String>,
    ) -> Self {
        Self {
            command,
            code: code.into(),
            message: message.into(),
            fix: fix.into(),
            subject: None,
        }
    }

    fn with_subject(mut self, subject: impl Into<String>) -> Self {
        self.subject = Some(subject.into());
        self
    }
}

/// Run the CLI from UTF-8 arguments excluding the executable name.
#[must_use]
pub fn run(args: impl IntoIterator<Item = String>) -> CommandOutput {
    let (mode, command) = match parse_args(args) {
        Ok(parsed) => parsed,
        Err((mode, diagnostic)) => return refusal(mode, exit::USAGE, &diagnostic, None),
    };

    match command {
        Command::Help => help(mode),
        Command::Validate(path) => validate_path(&path, mode),
        Command::Import(command) => import_path(&command, mode),
        Command::SolveProject(path) => unavailable(
            mode,
            "solve",
            &path.to_string_lossy(),
            "frankensim-extreal-program-f85xj.6.5",
        ),
        Command::Resume(run_id) => unavailable(
            mode,
            "solve",
            &run_id,
            "frankensim-extreal-program-f85xj.6.5",
        ),
        Command::Report(run_id) => unavailable(
            mode,
            "report",
            &run_id,
            "frankensim-extreal-program-f85xj.6.9",
        ),
        Command::Package(run_id) => unavailable(
            mode,
            "package",
            &run_id,
            "frankensim-extreal-program-f85xj.6.10",
        ),
    }
}

/// Run the CLI from platform arguments excluding the executable name.
/// Non-UTF-8 arguments are a stable usage refusal rather than a panic.
#[must_use]
pub fn run_os(args: impl IntoIterator<Item = OsString>) -> CommandOutput {
    let mut utf8 = Vec::new();
    for argument in args {
        match argument.into_string() {
            Ok(argument) => utf8.push(argument),
            Err(_) => {
                return refusal(
                    OutputMode::Text,
                    exit::USAGE,
                    &Diagnostic::new(
                        "arguments",
                        "cli-argument-encoding",
                        "an argument is not valid UTF-8",
                        "pass UTF-8 command names, paths, and run identifiers",
                    ),
                    None,
                );
            }
        }
    }
    run(utf8)
}

/// Validate already-loaded project bytes. This is the pure validation seam
/// used by conformance tests and embedders that own their own bounded I/O.
#[must_use]
pub fn validate_source(
    project: &str,
    source: &str,
    json_input: bool,
    json_output: bool,
) -> CommandOutput {
    let syntax = if json_input {
        ProjectSyntax::Json
    } else {
        ProjectSyntax::Sexpr
    };
    let mode = if json_output {
        OutputMode::Json
    } else {
        OutputMode::Text
    };
    validate_loaded(project, source, syntax, mode)
}

fn parse_args(
    args: impl IntoIterator<Item = String>,
) -> Result<(OutputMode, Command), (OutputMode, Diagnostic)> {
    let mut mode = OutputMode::Text;
    let mut saw_json = false;
    let mut positional = Vec::new();
    for argument in args {
        if argument == "--json" {
            if saw_json {
                return Err((
                    OutputMode::Json,
                    Diagnostic::new(
                        "arguments",
                        "cli-duplicate-flag",
                        "`--json` was provided more than once",
                        "provide `--json` at most once",
                    ),
                ));
            }
            saw_json = true;
            mode = OutputMode::Json;
        } else {
            positional.push(argument);
        }
    }

    let command = match positional.as_slice() {
        [flag] if flag == "--help" || flag == "help" => Command::Help,
        [verb, project] if verb == "validate" && is_operand(project) => {
            Command::Validate(PathBuf::from(project))
        }
        [verb, rest @ ..] if verb == "import" => {
            Command::Import(parse_import_args(rest).map_err(|diagnostic| (mode, diagnostic))?)
        }
        [verb, project] if verb == "solve" && is_operand(project) => {
            Command::SolveProject(PathBuf::from(project))
        }
        [verb, resume, run_id] if verb == "solve" && resume == "--resume" && is_operand(run_id) => {
            Command::Resume(run_id.clone())
        }
        [verb, run_id] if verb == "report" && is_operand(run_id) => Command::Report(run_id.clone()),
        [verb, run_id] if verb == "package" && is_operand(run_id) => {
            Command::Package(run_id.clone())
        }
        _ => {
            return Err((
                mode,
                Diagnostic::new(
                    "arguments",
                    "cli-usage",
                    "arguments do not match the v0 command grammar",
                    USAGE,
                ),
            ));
        }
    };
    Ok((mode, command))
}

fn parse_import_args(args: &[String]) -> Result<ImportCommand, Diagnostic> {
    if args.len() < 5 || !is_operand(&args[0]) || !is_operand(&args[1]) || !is_operand(&args[2]) {
        return Err(import_usage_diagnostic());
    }
    let mut unit = None;
    let mut max_hole_edges = None;
    let mut step_root = None;
    let mut target_h = None;
    let mut index = 3usize;
    while index < args.len() {
        let flag = &args[index];
        let Some(value) = args.get(index + 1) else {
            return Err(import_usage_diagnostic());
        };
        match flag.as_str() {
            "--unit" if unit.is_none() && is_operand(value) => unit = Some(value.clone()),
            "--max-hole-edges" if max_hole_edges.is_none() => {
                max_hole_edges = value.parse::<usize>().ok();
                if max_hole_edges.is_none() {
                    return Err(Diagnostic::new(
                        "import",
                        "cli-import-argument",
                        format!(
                            "`--max-hole-edges` requires a non-negative integer; got `{value}`"
                        ),
                        USAGE,
                    ));
                }
            }
            "--step-root" if step_root.is_none() => {
                step_root = value.parse::<u64>().ok().filter(|root| *root > 0);
                if step_root.is_none() {
                    return Err(Diagnostic::new(
                        "import",
                        "cli-import-argument",
                        format!("`--step-root` requires a positive integer; got `{value}`"),
                        USAGE,
                    ));
                }
            }
            "--target-h" if target_h.is_none() => {
                target_h = value
                    .parse::<f64>()
                    .ok()
                    .filter(|spacing| spacing.is_finite() && *spacing > 0.0);
                if target_h.is_none() {
                    return Err(Diagnostic::new(
                        "import",
                        "cli-import-argument",
                        format!("`--target-h` requires a finite positive number; got `{value}`"),
                        USAGE,
                    ));
                }
            }
            _ => return Err(import_usage_diagnostic()),
        }
        index += 2;
    }
    let unit = unit.ok_or_else(import_usage_diagnostic)?;
    let policy = match (max_hole_edges, step_root, target_h) {
        (Some(max_hole_edges), None, None) => ImportPolicy::Mesh { max_hole_edges },
        (None, Some(root_id), Some(target_h)) => ImportPolicy::FacetedStep { root_id, target_h },
        _ => return Err(import_usage_diagnostic()),
    };
    Ok(ImportCommand {
        project: PathBuf::from(&args[0]),
        source: PathBuf::from(&args[1]),
        ledger: PathBuf::from(&args[2]),
        unit,
        policy,
    })
}

fn import_usage_diagnostic() -> Diagnostic {
    Diagnostic::new(
        "import",
        "cli-import-usage",
        "import requires one project, one raw source, one ledger, one unit, and exactly one format policy",
        USAGE,
    )
}

fn is_operand(value: &str) -> bool {
    !value.is_empty() && !value.starts_with('-')
}

fn help(mode: OutputMode) -> CommandOutput {
    let stdout = match mode {
        OutputMode::Text => format!("{USAGE}\n"),
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"help\",\"status\":\"ok\",\"usage\":");
            push_json_string(&mut out, USAGE);
            out.push_str("}\n");
            out
        }
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

fn import_path(command: &ImportCommand, mode: OutputMode) -> CommandOutput {
    let project_label = command.project.to_string_lossy();
    let decoded = match read_project_for_import(&command.project, mode) {
        Ok(decoded) => decoded,
        Err(output) => return output,
    };
    let findings = decoded.findings();
    if !findings.is_empty() {
        let mut stderr = String::new();
        for finding in &findings {
            stderr.push_str(&format_diagnostic(
                mode,
                &Diagnostic::new(
                    "import",
                    finding.code,
                    finding.what.clone(),
                    finding.fix.clone(),
                )
                .with_subject(project_label.as_ref()),
            ));
        }
        return CommandOutput {
            exit_code: exit::REFUSED,
            stdout: format_result(
                mode,
                "import",
                "refused",
                &project_label,
                None,
                findings.len(),
            ),
            stderr,
        };
    }
    let Some(geometry) = decoded.spec.geometry.as_ref() else {
        return import_refusal(
            mode,
            &project_label,
            "project-geometry-missing",
            "project has no geometry section",
            "declare exactly one imported geometry receipt row for this command",
        );
    };
    if geometry.len() != 1 {
        return import_refusal(
            mode,
            &project_label,
            "cli-import-source-count",
            format!(
                "the v0 import command requires exactly one geometry row; the project declares {}",
                geometry.len()
            ),
            "import one reference enclosure, or use the library orchestration surface for a multi-source project",
        );
    }
    let artifact = &geometry[0];
    let declared_memory = decoded
        .spec
        .budgets
        .as_ref()
        .and_then(|budgets| usize::try_from(budgets.memory_bytes).ok())
        .unwrap_or(usize::MAX);
    let source_cap = GeometryImportLimits::DEFAULT
        .max_source_bytes
        .min(declared_memory);
    let source_bytes =
        match read_raw_import_source(&command.source, source_cap, mode, &project_label) {
            Ok(bytes) => bytes,
            Err(output) => return output,
        };
    let ledger_path = match command.ledger.to_str() {
        Some(path) if !path.is_empty() => path,
        _ => {
            return import_input_refusal(
                mode,
                &project_label,
                "cli-import-ledger-path",
                "ledger path is not valid non-empty UTF-8",
                "provide a UTF-8 SQLite ledger path",
            );
        }
    };

    let mut raw = RawGeometryLibrary::new();
    match command.policy {
        ImportPolicy::Mesh { max_hole_edges } => {
            raw.insert_mesh(
                artifact,
                command.source.to_string_lossy(),
                source_bytes,
                command.unit.clone(),
                max_hole_edges,
                Vec::new(),
            );
        }
        ImportPolicy::FacetedStep { root_id, target_h } => {
            raw.insert_faceted_step(
                artifact,
                command.source.to_string_lossy(),
                source_bytes,
                root_id,
                command.unit.clone(),
                target_h,
                Vec::new(),
            );
        }
    }
    let ledger = match fs_ledger::Ledger::open(ledger_path) {
        Ok(ledger) => ledger,
        Err(error) => {
            return import_input_refusal(
                mode,
                &project_label,
                "cli-import-ledger-open",
                format!("cannot open import ledger `{ledger_path}`: {error}"),
                "provide a writable ledger path whose parent directory exists, then retry",
            );
        }
    };

    let mut limits = GeometryImportLimits::DEFAULT;
    limits.max_sources = 1;
    if let Some(budgets) = decoded.spec.budgets.as_ref() {
        let memory_bytes = usize::try_from(budgets.memory_bytes).unwrap_or(usize::MAX);
        limits.max_source_bytes = limits.max_source_bytes.min(memory_bytes);
        limits.max_total_source_bytes = limits.max_total_source_bytes.min(memory_bytes);
    }
    let seed = decoded.spec.seeds.as_ref().map_or(0, |seeds| seeds.master);
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed,
                kernel_id: 0x66_73_63_6c_69_69_6d_70,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        import_project_geometry(&decoded.spec, &raw, &ledger, limits, &cx)
    });
    match result {
        Ok(run) => format_import_success(mode, &project_label, ledger_path, &run),
        Err(error) => {
            let mut message = error.what;
            if let Some(recorded) = error.recorded {
                let _ = write!(
                    message,
                    "; refusal retained as ledger operation {} and diagnostic artifact {}",
                    recorded.op_id,
                    recorded.diagnostic_artifact.to_hex()
                );
            }
            import_refusal(mode, &project_label, error.code, message, error.fix)
        }
    }
}

fn read_project_for_import(
    path: &Path,
    mode: OutputMode,
) -> Result<fs_project::DecodedProject, CommandOutput> {
    let label = path.to_string_lossy();
    let syntax = match path.extension().and_then(|extension| extension.to_str()) {
        Some("fsim") => ProjectSyntax::Sexpr,
        Some("json") => ProjectSyntax::Json,
        _ => {
            return Err(import_input_refusal(
                mode,
                &label,
                "cli-input-format",
                format!("project `{label}` has no admitted .fsim or .json extension"),
                "name canonical s-expression projects *.fsim and canonical JSON projects *.json",
            ));
        }
    };
    let metadata = std::fs::metadata(path).map_err(|error| {
        import_input_refusal(
            mode,
            &label,
            "cli-input-read",
            format!("cannot inspect project: {error}"),
            "provide a readable regular UTF-8 project file",
        )
    })?;
    if !metadata.is_file() || metadata.len() > MAX_PROJECT_BYTES {
        return Err(import_input_refusal(
            mode,
            &label,
            "cli-input-too-large",
            format!("project must be a regular file no larger than {MAX_PROJECT_BYTES} bytes"),
            "provide a stable bounded canonical project file",
        ));
    }
    let file = std::fs::File::open(path).map_err(|error| {
        import_input_refusal(
            mode,
            &label,
            "cli-input-read",
            format!("cannot open project: {error}"),
            "provide a readable regular UTF-8 project file",
        )
    })?;
    let mut source = String::new();
    file.take(MAX_PROJECT_BYTES.saturating_add(1))
        .read_to_string(&mut source)
        .map_err(|error| {
            import_input_refusal(
                mode,
                &label,
                "cli-input-read",
                format!("cannot read UTF-8 project: {error}"),
                "provide a readable regular UTF-8 project file",
            )
        })?;
    if u64::try_from(source.len()).map_or(true, |length| length > MAX_PROJECT_BYTES) {
        return Err(import_input_refusal(
            mode,
            &label,
            "cli-input-too-large",
            format!("project exceeded the {MAX_PROJECT_BYTES}-byte cap while being read"),
            "retry against a stable file no larger than the documented cap",
        ));
    }
    let decoded = match syntax {
        ProjectSyntax::Sexpr => fs_project::parse_sexpr(&source),
        ProjectSyntax::Json => fs_project::parse_json(&source),
    }
    .map_err(|error| import_refusal(mode, &label, error.code, error.detail, error.hint))?;
    Ok(decoded)
}

fn read_raw_import_source(
    path: &Path,
    source_cap: usize,
    mode: OutputMode,
    project: &str,
) -> Result<Vec<u8>, CommandOutput> {
    let label = path.to_string_lossy();
    let cap = u64::try_from(source_cap).unwrap_or(u64::MAX);
    let metadata = std::fs::metadata(path).map_err(|error| {
        import_input_refusal(
            mode,
            project,
            "cli-import-source-read",
            format!("cannot inspect source `{label}`: {error}"),
            "provide the readable regular raw geometry file named by this import attempt",
        )
    })?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err(import_input_refusal(
            mode,
            project,
            "cli-import-source-size",
            format!("raw source must be a regular file no larger than {cap} bytes"),
            "reduce the source or use the library surface with an explicitly budgeted envelope",
        ));
    }
    let file = std::fs::File::open(path).map_err(|error| {
        import_input_refusal(
            mode,
            project,
            "cli-import-source-read",
            format!("cannot open source `{label}`: {error}"),
            "provide the readable regular raw geometry file named by this import attempt",
        )
    })?;
    let mut bytes = Vec::new();
    file.take(cap.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|error| {
            import_input_refusal(
                mode,
                project,
                "cli-import-source-read",
                format!("cannot read source `{label}`: {error}"),
                "retry against a stable regular raw geometry file",
            )
        })?;
    if u64::try_from(bytes.len()).map_or(true, |length| length > cap) {
        return Err(import_input_refusal(
            mode,
            project,
            "cli-import-source-size",
            format!("raw source exceeded the {cap}-byte cap while being read"),
            "retry against a stable source no larger than the documented cap",
        ));
    }
    Ok(bytes)
}

fn format_import_success(
    mode: OutputMode,
    project: &str,
    ledger: &str,
    run: &GeometryImportRun,
) -> CommandOutput {
    let project_hash = run.project_hash.to_hex();
    let summary = run.summary_artifact.to_hex();
    let stdout = match mode {
        OutputMode::Text => format!(
            "status=ok\ncommand=import\nproject={}\nproject_hash={project_hash}\nledger={}\nop_id={}\nsummary_artifact={summary}\nartifact_count={}\nassignment_table={}\nauthority=retained-import-and-assignment-evidence\nno_claim={}\n",
            escape_text(project),
            escape_text(ledger),
            run.op_id,
            run.artifacts.len(),
            escape_text(&run.assignment_table),
            escape_text(GeometryImportRun::no_claim()),
        ),
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"import\",\"status\":\"ok\",\"project\":");
            push_json_string(&mut out, project);
            out.push_str(",\"project_hash\":");
            push_json_string(&mut out, &project_hash);
            out.push_str(",\"ledger\":");
            push_json_string(&mut out, ledger);
            let _ = write!(out, ",\"op_id\":{},\"summary_artifact\":", run.op_id);
            push_json_string(&mut out, &summary);
            let _ = write!(out, ",\"artifact_count\":{}", run.artifacts.len());
            out.push_str(",\"assignment_table\":");
            push_json_string(&mut out, &run.assignment_table);
            out.push_str(
                ",\"authority\":\"retained-import-and-assignment-evidence\",\"no_claim\":",
            );
            push_json_string(&mut out, GeometryImportRun::no_claim());
            out.push_str("}\n");
            out
        }
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

fn import_refusal(
    mode: OutputMode,
    project: &str,
    code: impl Into<String>,
    message: impl Into<String>,
    fix: impl Into<String>,
) -> CommandOutput {
    refusal(
        mode,
        exit::REFUSED,
        &Diagnostic::new("import", code, message, fix).with_subject(project),
        Some(("import", "refused", project, None, 1)),
    )
}

fn import_input_refusal(
    mode: OutputMode,
    project: &str,
    code: impl Into<String>,
    message: impl Into<String>,
    fix: impl Into<String>,
) -> CommandOutput {
    refusal(
        mode,
        exit::INPUT,
        &Diagnostic::new("import", code, message, fix).with_subject(project),
        Some(("import", "refused", project, None, 0)),
    )
}

fn validate_path(path: &Path, mode: OutputMode) -> CommandOutput {
    let label = path.to_string_lossy();
    let syntax = match path.extension().and_then(|extension| extension.to_str()) {
        Some("fsim") => ProjectSyntax::Sexpr,
        Some("json") => ProjectSyntax::Json,
        _ => {
            return refusal(
                mode,
                exit::INPUT,
                &Diagnostic::new(
                    "validate",
                    "cli-input-format",
                    format!("project `{label}` has no admitted .fsim or .json extension"),
                    "name canonical s-expression projects *.fsim and canonical JSON projects *.json",
                )
                .with_subject(label.as_ref()),
                Some(("validate", "refused", label.as_ref(), None, 0)),
            );
        }
    };

    let metadata = match std::fs::metadata(path) {
        Ok(metadata) => metadata,
        Err(error) => {
            return input_error(mode, &label, format!("cannot inspect project: {error}"));
        }
    };
    if !metadata.is_file() {
        return input_error(
            mode,
            &label,
            "project path is not a regular file".to_string(),
        );
    }
    if metadata.len() > MAX_PROJECT_BYTES {
        return refusal(
            mode,
            exit::INPUT,
            &Diagnostic::new(
                "validate",
                "cli-input-too-large",
                format!(
                    "project is {} bytes; the CLI cap is {MAX_PROJECT_BYTES}",
                    metadata.len()
                ),
                "reduce the project description or split external payloads into content-addressed artifacts",
            )
            .with_subject(label.as_ref()),
            Some(("validate", "refused", label.as_ref(), None, 0)),
        );
    }
    let file = match std::fs::File::open(path) {
        Ok(file) => file,
        Err(error) => return input_error(mode, &label, format!("cannot open project: {error}")),
    };
    let mut source = String::new();
    match file
        .take(MAX_PROJECT_BYTES.saturating_add(1))
        .read_to_string(&mut source)
    {
        Ok(_) => {}
        Err(error) => {
            return input_error(mode, &label, format!("cannot read UTF-8 project: {error}"));
        }
    }
    if u64::try_from(source.len()).map_or(true, |length| length > MAX_PROJECT_BYTES) {
        return refusal(
            mode,
            exit::INPUT,
            &Diagnostic::new(
                "validate",
                "cli-input-too-large",
                format!("project exceeded the {MAX_PROJECT_BYTES}-byte cap while being read"),
                "retry against a stable file no larger than the documented cap",
            )
            .with_subject(label.as_ref()),
            Some(("validate", "refused", label.as_ref(), None, 0)),
        );
    }
    validate_loaded(&label, &source, syntax, mode)
}

fn input_error(mode: OutputMode, project: &str, message: String) -> CommandOutput {
    refusal(
        mode,
        exit::INPUT,
        &Diagnostic::new(
            "validate",
            "cli-input-read",
            message,
            "provide a readable regular UTF-8 project file",
        )
        .with_subject(project),
        Some(("validate", "refused", project, None, 0)),
    )
}

fn validate_loaded(
    project: &str,
    source: &str,
    syntax: ProjectSyntax,
    mode: OutputMode,
) -> CommandOutput {
    let decoded = match syntax {
        ProjectSyntax::Sexpr => fs_project::parse_sexpr(source),
        ProjectSyntax::Json => fs_project::parse_json(source),
    };
    let decoded = match decoded {
        Ok(decoded) => decoded,
        Err(error) => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new("validate", error.code, error.detail, error.hint)
                    .with_subject(project),
                Some(("validate", "refused", project, None, 0)),
            );
        }
    };

    let findings = decoded.findings();
    if !findings.is_empty() {
        let mut stderr = String::new();
        for finding in &findings {
            let diagnostic = Diagnostic::new(
                "validate",
                finding.code,
                finding.what.clone(),
                finding.fix.clone(),
            )
            .with_subject(project);
            stderr.push_str(&format_diagnostic(mode, &diagnostic));
        }
        return CommandOutput {
            exit_code: exit::REFUSED,
            stdout: format_result(mode, "validate", "refused", project, None, findings.len()),
            stderr,
        };
    }

    let hash = decoded.hash().to_hex();
    let stdout = match mode {
        OutputMode::Text => {
            let project = escape_text(project);
            format!(
                "status=ok\ncommand=validate\nproject={project}\nproject_hash={hash}\nfsim_version={}\nfinding_count=0\nauthority={VALIDATION_AUTHORITY}\nno_claim={VALIDATION_NO_CLAIM}\n",
                fs_project::FSIM_VERSION
            )
        }
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"validate\",\"status\":\"ok\",\"project\":");
            push_json_string(&mut out, project);
            out.push_str(",\"project_hash\":");
            push_json_string(&mut out, &hash);
            let _ = write!(
                out,
                ",\"fsim_version\":{},\"finding_count\":0,\"default_receipt_count\":0,\"canonicalization_receipt\":false,\"authority\":",
                fs_project::FSIM_VERSION
            );
            push_json_string(&mut out, VALIDATION_AUTHORITY);
            out.push_str(",\"no_claim\":");
            push_json_string(&mut out, VALIDATION_NO_CLAIM);
            out.push_str("}\n");
            out
        }
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

fn unavailable(
    mode: OutputMode,
    command: &'static str,
    subject: &str,
    dependency: &'static str,
) -> CommandOutput {
    let message = format!(
        "`{command}` is reserved but cannot execute until `{dependency}` supplies its authoritative product stage"
    );
    let fix = format!(
        "complete and verify `{dependency}`; do not substitute a skeleton run or placeholder artifact"
    );
    let mut output = refusal(
        mode,
        exit::UNAVAILABLE,
        &Diagnostic::new(command, "cli-stage-unavailable", message, fix).with_subject(subject),
        Some((command, "unavailable", subject, None, 0)),
    );
    if mode == OutputMode::Json {
        let marker = "}\n";
        if let Some(at) = output.stdout.rfind(marker) {
            output
                .stdout
                .insert_str(at, &format!(",\"dependency\":\"{dependency}\""));
        }
    } else {
        let _ = writeln!(output.stdout, "dependency={dependency}");
    }
    output
}

fn refusal(
    mode: OutputMode,
    exit_code: u8,
    diagnostic: &Diagnostic,
    result: Option<(&str, &str, &str, Option<&str>, usize)>,
) -> CommandOutput {
    let stdout = result.map_or_else(String::new, |(command, status, subject, hash, findings)| {
        format_result(mode, command, status, subject, hash, findings)
    });
    CommandOutput {
        exit_code,
        stdout,
        stderr: format_diagnostic(mode, diagnostic),
    }
}

fn format_result(
    mode: OutputMode,
    command: &str,
    status: &str,
    subject: &str,
    hash: Option<&str>,
    finding_count: usize,
) -> String {
    match mode {
        OutputMode::Text => {
            let subject = escape_text(subject);
            let mut out = format!(
                "status={status}\ncommand={command}\nsubject={subject}\nfinding_count={finding_count}\n"
            );
            if let Some(hash) = hash {
                let _ = writeln!(out, "project_hash={hash}");
            }
            out
        }
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":");
            push_json_string(&mut out, command);
            out.push_str(",\"status\":");
            push_json_string(&mut out, status);
            out.push_str(",\"subject\":");
            push_json_string(&mut out, subject);
            if let Some(hash) = hash {
                out.push_str(",\"project_hash\":");
                push_json_string(&mut out, hash);
            }
            let _ = writeln!(out, ",\"finding_count\":{finding_count}}}");
            out
        }
    }
}

fn format_diagnostic(mode: OutputMode, diagnostic: &Diagnostic) -> String {
    match mode {
        OutputMode::Text => format!(
            "ERROR {}: {}\nFIX: {}\n",
            escape_text(&diagnostic.code),
            escape_text(&diagnostic.message),
            escape_text(&diagnostic.fix)
        ),
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, DIAGNOSTIC_SCHEMA);
            out.push_str(",\"command\":");
            push_json_string(&mut out, diagnostic.command);
            out.push_str(",\"severity\":\"error\",\"code\":");
            push_json_string(&mut out, &diagnostic.code);
            out.push_str(",\"message\":");
            push_json_string(&mut out, &diagnostic.message);
            out.push_str(",\"fix\":");
            push_json_string(&mut out, &diagnostic.fix);
            if let Some(subject) = &diagnostic.subject {
                out.push_str(",\"subject\":");
                push_json_string(&mut out, subject);
            }
            out.push_str("}\n");
            out
        }
    }
}

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character <= '\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out.push('"');
}

fn escape_text(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(out, "\\u{{{:x}}}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out
}
