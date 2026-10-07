//! FrankenScript executor v0 (bead frankensim-rc-root-q61wp.38): run an
//! admitted `(study …)` program through the same stage drivers the `.fsim`
//! verbs use, so a program and the equivalent command sequence produce the
//! same run, receipts, report and package.
//!
//! # Pipeline
//!
//! 1. Parse (`*.fs` s-expression or `*.fs.json` canonical JSON) and admit
//!    through `fs_ir::admission::admit`: unit, budget, capability and
//!    explicit-pillar errors refuse here, before anything executes.
//! 2. Re-recognize the LOWERED program (`Study::from_node` on the lowering
//!    receipt's canonical text) and bind every clause statically against
//!    [`EXECUTABLE_VERBS`]. An unbound verb, malformed argument, unreadable
//!    project, project-hash mismatch, or explicit that disagrees with the
//!    project refuses with structured diagnostics before any stage runs.
//! 3. Execute the bound steps in program order. Each step calls the CLI's own
//!    stage driver (`import`, `solve`, or the validate → solve → report →
//!    package `run` workflow); the first non-success step stops the program
//!    and its output is returned with the steps that completed.
//!
//! # Program shape
//!
//! ```text
//! (study "heatsink-fan-journey-a"
//!   (seed 0x7) (versions (constellation :lock "2026-07"))
//!   (budget (wall 60s) (mem 64MiB))
//!   (capability :cores 1 :mem 64MiB :wall 60s :ops (cooling.*))
//!   (let project (cooling.project "heatsink-fan.fsim" :hash "…"))
//!   (cooling.import project :sources ("heatsink.stl") :unit "m" :max-hole-edges 0)
//!   (cooling.run project :materials ("aa6061.fsmcdpk")))
//!
//! (study "marquee-bracket-2d-journey-b"
//!   (seed 0x539) (versions (constellation :lock "2026-07"))
//!   (budget (wall 300s) (mem 512MiB))
//!   (capability :cores 1 :mem 512MiB :wall 300s :ops (study.*))
//!   (study.run "bracket-2d.fsim" :hash "…"))
//! ```
//!
//! Paths resolve against the program file's directory. The study seed must
//! equal the project's `seeds.root`, and a declared `(wall …)` / `(mem …)`
//! budget must cover the project's solve-time and memory budgets: the program
//! cannot silently run under different explicits than it states.
//!
//! # No-claim boundaries
//!
//! v0 binds the cooling project pipeline (`cooling.project`, `cooling.import`,
//! `cooling.solve`, `cooling.run`), native uncertainty/sensitivity studies
//! (`cooling.study project :source "study.fsim" [:hash "…"] [:budget N]`), and
//! the canonical study driver on any `.fsim` study file (`study.run
//! "study.fsim" [:hash "…"] [:budget N]`, e.g. the 2-D marquee). Native study
//! assets resolve against the study source; its sampling seed is distinct from
//! the physical project seed, while `study.run` requires the program seed to
//! equal the study file's `seeds.root` and leaves budgets to the study driver.
//! Budget stops retain the study run ID for ordinary `study --resume`;
//! whole-program resume is not implemented. The catalog's physics operators
//! (`flux.*`, `ascent.*`, …) are admitted by fs-ir but have no stage binding
//! here and refuse as not executable. Project and study files are re-read by
//! each stage driver exactly as the CLI verbs read them; `:hash` pins are
//! checked once during binding, not atomically with subsequent file reads.
//! Per-step wall admission is not aggregate program wall-time metering.

use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use fs_ir::admission::{AdmissionContext, RegimePolicy, Severity, admit};
use fs_ir::ast::{CountUnit, Node, NodeKind};
use fs_ir::study::Study;

#[path = "frankenscript/native_study.rs"]
mod native_study;

use crate::cards::CardPackKind;
use crate::{
    CommandOutput, DIAGNOSTIC_SCHEMA, Diagnostic, ImportCommand, ImportPolicy, MAX_PROJECT_BYTES,
    OutputMode, RESULT_SCHEMA, escape_text, exit, format_diagnostic, import_path, push_json_string,
    read_project_for_solve, run_workflow_path, solve_path,
};

/// Verbs with a stage binding, in catalog order.
pub const EXECUTABLE_VERBS: [&str; 6] = [
    "cooling.project",
    "cooling.import",
    "cooling.solve",
    "cooling.run",
    "cooling.study",
    "study.run",
];

const LENGTH: fs_qty::Dims = fs_qty::Dims([1, 0, 0, 0, 0, 0]);
const TIME: fs_qty::Dims = fs_qty::Dims([0, 0, 1, 0, 0, 0]);

/// Domain separator of the program identity hash.
const PROGRAM_DOMAIN: &str = "frankensim.frankenscript.program.v0";

/// True when `path` names a FrankenScript program (`*.fs` or `*.fs.json`).
pub(crate) fn is_program(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".fs") || name.ends_with(".fs.json")
}

/// One bound, not yet executed step.
#[derive(Debug)]
enum Step {
    Study(native_study::StudyStep),
    Import(ImportCommand),
    Solve {
        project: PathBuf,
        cards: Vec<(CardPackKind, PathBuf)>,
    },
    Run {
        project: PathBuf,
        cards: Vec<(CardPackKind, PathBuf)>,
    },
    StudyFile {
        path: PathBuf,
        budget: Option<String>,
    },
}

impl Step {
    const fn verb(&self) -> &'static str {
        match self {
            Self::Study(_) => "cooling.study",
            Self::Import(_) => "cooling.import",
            Self::Solve { .. } => "cooling.solve",
            Self::Run { .. } => "cooling.run",
            Self::StudyFile { .. } => "study.run",
        }
    }
}

/// A `cooling.project` binding.
#[derive(Debug, Clone)]
struct ProjectBinding {
    path: PathBuf,
    hash: String,
}

struct Refusals {
    mode: OutputMode,
    label: String,
    stderr: String,
    count: usize,
}

impl Refusals {
    fn push(&mut self, code: &str, message: impl Into<String>, fix: impl Into<String>) {
        let diagnostic =
            Diagnostic::new("run", code.to_string(), message, fix).with_subject(self.label.clone());
        self.stderr
            .push_str(&format_diagnostic(self.mode, &diagnostic));
        self.count += 1;
    }

    fn finish(self, exit_code: u8) -> CommandOutput {
        CommandOutput {
            exit_code,
            stdout: crate::format_result(
                self.mode,
                "run",
                "refused",
                &self.label,
                None,
                self.count,
            ),
            stderr: self.stderr,
        }
    }
}

/// Keyword arguments of a verb form after its positional operands.
fn keywords<'a>(
    items: &'a [Node],
    refusals: &mut Refusals,
    verb: &str,
) -> Option<(Vec<&'a Node>, Vec<(&'a str, &'a Node)>)> {
    let mut positional = Vec::new();
    let mut named = Vec::new();
    let mut index = 1;
    while index < items.len() {
        if let NodeKind::Keyword(key) = &items[index].kind {
            let Some(value) = items.get(index + 1) else {
                refusals.push(
                    "frankenscript-argument",
                    format!("`{verb}` keyword :{key} has no value"),
                    "supply a value after every keyword",
                );
                return None;
            };
            if named.iter().any(|(k, _)| *k == key.as_str()) {
                refusals.push(
                    "frankenscript-argument",
                    format!("`{verb}` repeats keyword :{key}"),
                    "state each keyword once",
                );
                return None;
            }
            named.push((key.as_str(), value));
            index += 2;
        } else {
            if !named.is_empty() {
                refusals.push(
                    "frankenscript-argument",
                    format!("`{verb}` has a positional operand after its keywords"),
                    "put positional operands before keywords",
                );
                return None;
            }
            positional.push(&items[index]);
            index += 1;
        }
    }
    Some((positional, named))
}

fn string_of(node: &Node) -> Option<&str> {
    match &node.kind {
        NodeKind::Str(s) => Some(s),
        _ => None,
    }
}

fn string_list(node: &Node) -> Option<Vec<&str>> {
    node.items()?.iter().map(string_of).collect()
}

/// The positive integer value of `node`.
fn count_of(node: &Node) -> Option<u64> {
    match node.kind {
        NodeKind::Int(v) => u64::try_from(v).ok(),
        NodeKind::Seed(v) => Some(v),
        _ => None,
    }
}

struct Binder<'p> {
    base: PathBuf,
    ledger: &'p Path,
    seed: Option<u64>,
    wall_seconds: Option<f64>,
    projects: Vec<(String, ProjectBinding)>,
    studies: Vec<ProjectBinding>,
}

impl Binder<'_> {
    fn resolve(&self, relative: &str) -> PathBuf {
        let path = Path::new(relative);
        if path.is_absolute() {
            path.to_path_buf()
        } else {
            self.base.join(path)
        }
    }

    fn project_operand(
        &self,
        verb: &str,
        positional: &[&Node],
        refusals: &mut Refusals,
    ) -> Option<ProjectBinding> {
        let name = match positional {
            [only] => match &only.kind {
                NodeKind::Symbol(name) => name.as_str(),
                _ => "",
            },
            _ => "",
        };
        let found = self
            .projects
            .iter()
            .find(|(bound, _)| bound == name)
            .map(|(_, binding)| binding.clone());
        if found.is_none() {
            refusals.push(
                "frankenscript-unbound-project",
                format!("`{verb}` must name exactly one (let <name> (cooling.project …)) binding"),
                "bind the project with (let project (cooling.project \"p.fsim\")) and pass `project`",
            );
        }
        found
    }

    fn cards(
        &self,
        verb: &str,
        named: &[(&str, &Node)],
        refusals: &mut Refusals,
    ) -> Option<Vec<(CardPackKind, PathBuf)>> {
        let mut cards = Vec::new();
        for (key, value) in named {
            let kind = match *key {
                "materials" => CardPackKind::Material,
                "interfaces" => CardPackKind::Interface,
                other => {
                    refusals.push(
                        "frankenscript-argument",
                        format!("`{verb}` has no keyword :{other}"),
                        "use :materials (\"pack\" …) and :interfaces (\"pack\" …)",
                    );
                    return None;
                }
            };
            let Some(paths) = string_list(value) else {
                refusals.push(
                    "frankenscript-argument",
                    format!("`{verb}` :{key} must be a list of path strings"),
                    format!("write :{key} (\"a.fsmcdpk\" …)"),
                );
                return None;
            };
            cards.extend(paths.into_iter().map(|p| (kind, self.resolve(p))));
        }
        Some(cards)
    }

    /// `(study.run "study.fsim" [:budget <n>] [:hash "<blake3 of bytes>"])`:
    /// the canonical `frankensim study` driver on a study file whose root
    /// seed must equal the program's.
    fn bind_study(
        &mut self,
        positional: &[&Node],
        named: &[(&str, &Node)],
        refusals: &mut Refusals,
    ) -> Option<Step> {
        let path = match positional {
            [only] => string_of(only),
            _ => None,
        };
        let mut budget = None;
        let mut pin = None;
        for (key, value) in named {
            match (*key, &value.kind) {
                ("budget", NodeKind::Int(n)) if *n > 0 => budget = Some(n.to_string()),
                ("budget", NodeKind::Str(text)) => budget = Some(text.clone()),
                ("hash", NodeKind::Str(text)) => pin = Some(text.as_str()),
                _ => {
                    refusals.push(
                        "frankenscript-argument",
                        format!("`study.run` keyword :{key} is not admitted here"),
                        "use :budget <positive integer> and :hash \"<study hash>\"",
                    );
                    return None;
                }
            }
        }
        let Some(path) = path else {
            refusals.push(
                "frankenscript-argument",
                "`study.run` takes one study path string",
                "write (study.run \"study.fsim\")",
            );
            return None;
        };
        let path = self.resolve(path);
        if path.extension().and_then(|e| e.to_str()) != Some("fsim") {
            refusals.push(
                "frankenscript-argument",
                format!(
                    "`study.run` binds canonical .fsim studies only, not `{}`",
                    path.display()
                ),
                "run JSON studies with `frankensim study` directly",
            );
            return None;
        }
        let bytes = match std::fs::metadata(&path).and_then(|m| {
            if m.len() > MAX_PROJECT_BYTES {
                Err(std::io::Error::other("study exceeds the project size cap"))
            } else {
                std::fs::read(&path)
            }
        }) {
            Ok(bytes) => bytes,
            Err(error) => {
                refusals.push(
                    "frankenscript-input",
                    format!("cannot read study `{}`: {error}", path.display()),
                    "provide a readable .fsim study next to the program",
                );
                return None;
            }
        };
        let hash = fs_blake3::hash_bytes(&bytes).to_hex();
        if let Some(pin) = pin
            && pin != hash
        {
            refusals.push(
                "frankenscript-study-hash",
                format!(
                    "study `{}` hashes to {hash}, not the pinned {pin}",
                    path.display()
                ),
                "re-pin :hash after reviewing the study change, or restore the pinned study",
            );
        }
        let root_seed = std::str::from_utf8(&bytes)
            .ok()
            .and_then(|text| fs_ir::sexpr::parse(text).ok())
            .and_then(|node| {
                node.items()?
                    .iter()
                    .find(|clause| clause.head() == Some("seeds"))
                    .and_then(|seeds| {
                        let items = seeds.items()?;
                        let at = items
                            .iter()
                            .position(|i| matches!(&i.kind, NodeKind::Keyword(k) if k == "root"))?;
                        items.get(at + 1).and_then(count_of)
                    })
            });
        if root_seed.is_none() || root_seed != self.seed {
            refusals.push(
                "frankenscript-explicit-seed",
                format!(
                    "the study seed {:?} differs from the study file's seeds.root {root_seed:?}",
                    self.seed
                ),
                "state the study file's root seed in (seed …) so the program records the seed it runs",
            );
        }
        self.studies.push(ProjectBinding {
            path: path.clone(),
            hash,
        });
        Some(Step::StudyFile { path, budget })
    }

    #[allow(clippy::too_many_lines)] // one keyword grammar per verb
    fn bind(&mut self, clause: &Node, refusals: &mut Refusals) -> Option<Step> {
        let verb = clause.head().unwrap_or("");
        let items = clause.items().unwrap_or(&[]);
        if !EXECUTABLE_VERBS.contains(&verb) || verb == "cooling.project" {
            refusals.push(
                "frankenscript-not-executable",
                format!("not executable: verb `{verb}` has no stage binding"),
                format!(
                    "use one of the executable verbs ({}); `cooling.project` binds through (let …)",
                    EXECUTABLE_VERBS.join(", ")
                ),
            );
            return None;
        }
        let (positional, named) = keywords(items, refusals, verb)?;
        if verb == "study.run" {
            return self.bind_study(&positional, &named, refusals);
        }
        let project = self.project_operand(verb, &positional, refusals)?;
        match verb {
            "cooling.study" => native_study::bind(self, &project, &named, refusals).map(Step::Study),
            "cooling.import" => {
                let mut sources = None;
                let mut unit = None;
                let mut max_hole_edges = None;
                let mut step_root = None;
                let mut target_h = None;
                for (key, value) in &named {
                    match *key {
                        "sources" => sources = string_list(value),
                        "unit" => unit = string_of(value),
                        "max-hole-edges" => {
                            max_hole_edges = count_of(value).and_then(|v| usize::try_from(v).ok());
                        }
                        "step-root" => step_root = count_of(value),
                        "target-h" => {
                            target_h = match value.kind {
                                NodeKind::Qty { value, dims, .. } if dims == LENGTH => Some(value),
                                NodeKind::Float(v) => Some(v),
                                _ => None,
                            };
                        }
                        other => {
                            refusals.push(
                                "frankenscript-argument",
                                format!("`cooling.import` has no keyword :{other}"),
                                "use :sources, :unit, and :max-hole-edges or :step-root with :target-h",
                            );
                            return None;
                        }
                    }
                }
                let (Some(sources), Some(unit)) = (sources, unit) else {
                    refusals.push(
                        "frankenscript-argument",
                        "`cooling.import` needs :sources (\"path\" …) and :unit \"<unit>\"",
                        "e.g. (cooling.import project :sources (\"part.stl\") :unit \"m\" :max-hole-edges 0)",
                    );
                    return None;
                };
                let policy = match (max_hole_edges, step_root, target_h) {
                    (Some(max_hole_edges), None, None) => ImportPolicy::Mesh { max_hole_edges },
                    (None, Some(root_id), Some(target_h)) => {
                        ImportPolicy::FacetedStep { root_id, target_h }
                    }
                    _ => {
                        refusals.push(
                            "frankenscript-argument",
                            "`cooling.import` needs exactly one policy: :max-hole-edges <n>, or \
                             :step-root <id> with :target-h <length>",
                            "declare the mesh hole policy or the faceted-STEP root and spacing",
                        );
                        return None;
                    }
                };
                Some(Step::Import(ImportCommand {
                    project: project.path,
                    sources: sources.iter().map(|s| self.resolve(s)).collect(),
                    ledger: self.ledger.to_path_buf(),
                    unit: unit.to_string(),
                    policy,
                }))
            }
            "cooling.solve" => Some(Step::Solve {
                project: project.path,
                cards: self.cards(verb, &named, refusals)?,
            }),
            _ => Some(Step::Run {
                project: project.path,
                cards: self.cards(verb, &named, refusals)?,
            }),
        }
    }
}

/// Seconds and bytes a `(budget …)` clause declares, when stated.
fn declared_budget(budget: &Node) -> (Option<f64>, Option<u64>) {
    let (mut wall, mut mem) = (None, None);
    for clause in budget.items().unwrap_or(&[]).iter().skip(1) {
        let value = clause.items().and_then(|items| items.get(1));
        match (clause.head(), value.map(|v| &v.kind)) {
            (Some("wall"), Some(NodeKind::Qty { value, dims, .. })) if *dims == TIME => {
                wall = Some(*value);
            }
            (Some("mem"), Some(NodeKind::Count { value, unit })) if *unit != CountUnit::Cores => {
                mem = value.integral_bytes(*unit);
            }
            _ => {}
        }
    }
    (wall, mem)
}

fn read_program(path: &Path, refusals: &mut Refusals) -> Option<String> {
    let mut source = String::new();
    let read = std::fs::metadata(path)
        .map_err(|e| e.to_string())
        .and_then(|metadata| {
            if metadata.is_file() && metadata.len() <= MAX_PROJECT_BYTES {
                Ok(())
            } else {
                Err(format!(
                    "not a regular file within {MAX_PROJECT_BYTES} bytes"
                ))
            }
        })
        .and_then(|()| std::fs::File::open(path).map_err(|e| e.to_string()))
        .and_then(|file| {
            file.take(MAX_PROJECT_BYTES.saturating_add(1))
                .read_to_string(&mut source)
                .map_err(|e| e.to_string())
        });
    match read {
        Ok(_) if source.len() as u64 <= MAX_PROJECT_BYTES => Some(source),
        Ok(_) => {
            refusals.push(
                "frankenscript-input",
                "program grew past the size cap while being read",
                "retry against a stable program file",
            );
            None
        }
        Err(error) => {
            refusals.push(
                "frankenscript-input",
                format!("cannot read program: {error}"),
                "provide a readable UTF-8 *.fs or *.fs.json program",
            );
            None
        }
    }
}

/// Parse, admit, bind, and execute one program file against `ledger`.
#[allow(clippy::too_many_lines)] // the four pipeline phases in order
pub(crate) fn run_program_path(program: &Path, ledger: &Path, mode: OutputMode) -> CommandOutput {
    let label = program.to_string_lossy().into_owned();
    let mut refusals = Refusals {
        mode,
        label: label.clone(),
        stderr: String::new(),
        count: 0,
    };
    let Some(source) = read_program(program, &mut refusals) else {
        return refusals.finish(exit::INPUT);
    };
    let json = label.ends_with(".fs.json");
    let parsed = if json {
        fs_ir::json::parse(&source)
    } else {
        fs_ir::sexpr::parse(&source)
    };
    let node = match parsed {
        Ok(node) => node,
        Err(error) => {
            refusals.push("frankenscript-parse", error.detail, error.hint);
            return refusals.finish(exit::INPUT);
        }
    };

    // 1. Admission: the only authority for units, budgets and capabilities.
    let context = AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: std::collections::BTreeMap::new(),
        capability: None,
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };
    let report = admit(&node, &context);
    let mut warnings = String::new();
    for finding in &report.findings {
        let fix = finding
            .fixes
            .first()
            .map_or_else(String::new, |fix| fix.action.clone());
        let code = format!("admission-{}", finding.check);
        if finding.severity == Severity::Reject {
            refusals.push(&code, finding.what.clone(), fix);
        } else {
            match mode {
                OutputMode::Text => {
                    let _ = writeln!(
                        warnings,
                        "WARN {}: {}",
                        escape_text(&code),
                        escape_text(&finding.what)
                    );
                }
                OutputMode::Json => {
                    warnings.push_str("{\"schema\":");
                    push_json_string(&mut warnings, DIAGNOSTIC_SCHEMA);
                    warnings.push_str(",\"command\":\"run\",\"severity\":\"warn\",\"code\":");
                    push_json_string(&mut warnings, &code);
                    warnings.push_str(",\"message\":");
                    push_json_string(&mut warnings, &finding.what);
                    warnings.push_str(",\"fix\":");
                    push_json_string(&mut warnings, &fix);
                    warnings.push_str("}\n");
                }
            }
        }
    }
    if !report.admitted {
        if refusals.count == 0 {
            refusals.push(
                "admission-rejected",
                report.diagnosis(),
                "apply the ranked fixes above",
            );
        }
        return refusals.finish(exit::REFUSED);
    }
    let Some(lowered_text) = report.lowering.lowered_canonical() else {
        refusals.push(
            "admission-identity",
            "admission returned no lowered identity",
            "report this as an fs-ir defect",
        );
        return refusals.finish(exit::REFUSED);
    };
    let program_hash = fs_blake3::hash_domain(PROGRAM_DOMAIN, lowered_text.as_bytes()).to_hex();
    let lowered = match fs_ir::VersionedProgram::parse_sexpr(lowered_text) {
        Ok(program) => program,
        Err(error) => {
            refusals.push("admission-identity", error.detail, error.hint);
            return refusals.finish(exit::REFUSED);
        }
    };
    let study = match Study::from_node(lowered.program()) {
        Ok(study) => study,
        Err(error) => {
            refusals.push("frankenscript-study", error.detail, error.hint);
            return refusals.finish(exit::REFUSED);
        }
    };

    // 2. Static binding: nothing executes unless every clause binds.
    let mut binder = Binder {
        base: program
            .parent()
            .map_or_else(PathBuf::new, Path::to_path_buf),
        ledger,
        seed: study.seed,
        wall_seconds: study.budget.and_then(|budget| declared_budget(budget).0),
        projects: Vec::new(),
        studies: Vec::new(),
    };
    for (name, value) in &study.lets {
        if value.head() != Some("cooling.project") {
            refusals.push(
                "frankenscript-not-executable",
                format!(
                    "not executable: (let {name} …) binds `{}`, which has no stage binding",
                    value.head().unwrap_or("a literal")
                ),
                "v0 executes (let <name> (cooling.project \"p.fsim\" [:hash \"…\"])) bindings only",
            );
            continue;
        }
        let items = value.items().unwrap_or(&[]);
        let Some((positional, named)) = keywords(items, &mut refusals, "cooling.project") else {
            continue;
        };
        let path = match positional.as_slice() {
            [only] => string_of(only),
            _ => None,
        };
        let pin = named
            .iter()
            .find(|(key, _)| *key == "hash")
            .and_then(|(_, v)| string_of(v));
        let (Some(path), true) = (path, named.iter().all(|(key, _)| *key == "hash")) else {
            refusals.push(
                "frankenscript-argument",
                format!(
                    "(let {name} (cooling.project …)) takes one path string and an optional :hash"
                ),
                "write (cooling.project \"p.fsim\" :hash \"<project hash>\")",
            );
            continue;
        };
        let path = binder.resolve(path);
        let decoded = match read_project_for_solve(&path, mode) {
            Ok(decoded) => decoded,
            Err(output) => {
                refusals.stderr.push_str(&output.stderr);
                refusals.count += 1;
                continue;
            }
        };
        let hash = decoded.hash().to_hex();
        if let Some(pin) = pin
            && pin != hash
        {
            refusals.push(
                "frankenscript-project-hash",
                format!(
                    "project `{}` hashes to {hash}, not the pinned {pin}",
                    path.display()
                ),
                "re-pin :hash after reviewing the project change, or restore the pinned project",
            );
        }
        let project_seed = decoded.spec.seeds.as_ref().map(|seeds| seeds.root);
        if project_seed != study.seed {
            refusals.push(
                "frankenscript-explicit-seed",
                format!(
                    "the study seed {:?} differs from the project's seeds.root {project_seed:?}",
                    study.seed
                ),
                "state the project's root seed in (seed …) so the program records the seed it runs",
            );
        }
        if let (Some(budget), Some(project_budget)) = (study.budget, decoded.spec.budgets.as_ref())
        {
            let (wall, mem) = declared_budget(budget);
            if let Some(wall) = wall
                && wall < project_budget.solve_time.value
            {
                refusals.push(
                    "frankenscript-explicit-budget",
                    format!(
                        "the study wall budget {wall} s is below the project's solve-time {} s",
                        project_budget.solve_time.value
                    ),
                    "raise (budget (wall …)) to cover the project's declared solve time",
                );
            }
            if let Some(mem) = mem
                && mem < project_budget.memory_bytes
            {
                refusals.push(
                    "frankenscript-explicit-budget",
                    format!(
                        "the study memory budget {mem} B is below the project's {} B",
                        project_budget.memory_bytes
                    ),
                    "raise (budget (mem …)) to cover the project's declared memory",
                );
            }
        }
        binder
            .projects
            .push(((*name).to_string(), ProjectBinding { path, hash }));
    }
    let mut steps = Vec::new();
    for clause in &study.body {
        if let Some(step) = binder.bind(clause, &mut refusals) {
            steps.push(step);
        }
    }
    if steps.is_empty() && refusals.count == 0 {
        refusals.push(
            "frankenscript-empty",
            "the study has no executable step",
            "add a (cooling.run project …) clause",
        );
    }
    if refusals.count > 0 {
        let mut output = refusals.finish(exit::REFUSED);
        output.stderr = format!("{warnings}{}", output.stderr);
        return output;
    }

    // 3. Execute in program order through the CLI's own stage drivers.
    let mut stderr = warnings;
    let mut records = Vec::with_capacity(steps.len());
    let mut status = exit::SUCCESS;
    for step in &steps {
        let output = match step {
            Step::Study(study) => study.run(ledger, mode),
            Step::Import(command) => import_path(command, mode),
            Step::Solve { project, cards } => solve_path(project, ledger, cards, mode),
            Step::Run { project, cards } => run_workflow_path(project, ledger, cards, mode),
            Step::StudyFile { path, budget } => {
                crate::study::study_path(path, ledger, budget.as_deref(), mode)
            }
        };
        stderr.push_str(&output.stderr);
        records.push((step.verb(), output.exit_code, output.stdout));
        if output.exit_code != exit::SUCCESS {
            status = output.exit_code;
            break;
        }
    }
    let completed = status == exit::SUCCESS;
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"run\",\"status\":");
            push_json_string(&mut out, if completed { "completed" } else { "stopped" });
            out.push_str(",\"subject\":");
            push_json_string(&mut out, &label);
            out.push_str(",\"study\":");
            push_json_string(&mut out, study.name);
            out.push_str(",\"program_hash\":");
            push_json_string(&mut out, &program_hash);
            out.push_str(",\"projects\":[");
            for (index, (name, binding)) in binder.projects.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                out.push_str("{\"name\":");
                push_json_string(&mut out, name);
                out.push_str(",\"path\":");
                push_json_string(&mut out, &binding.path.to_string_lossy());
                out.push_str(",\"hash\":");
                push_json_string(&mut out, &binding.hash);
                out.push('}');
            }
            out.push_str("],\"studies\":[");
            for (index, binding) in binder.studies.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                out.push_str("{\"path\":");
                push_json_string(&mut out, &binding.path.to_string_lossy());
                out.push_str(",\"hash\":");
                push_json_string(&mut out, &binding.hash);
                out.push('}');
            }
            out.push_str("],\"steps\":[");
            for (index, (verb, code, result)) in records.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                out.push_str("{\"verb\":");
                push_json_string(&mut out, verb);
                let _ = write!(out, ",\"exit\":{code},\"result\":");
                let trimmed = result.trim();
                out.push_str(if trimmed.is_empty() { "null" } else { trimmed });
                out.push('}');
            }
            let _ = writeln!(out, "],\"steps_planned\":{}}}", steps.len());
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status={}\ncommand=run\nsubject={}\nstudy={}\nprogram_hash={program_hash}\n",
                if completed { "completed" } else { "stopped" },
                escape_text(&label),
                escape_text(study.name),
            );
            for (name, binding) in &binder.projects {
                let _ = writeln!(
                    out,
                    "project={} path={} hash={}",
                    escape_text(name),
                    escape_text(&binding.path.to_string_lossy()),
                    binding.hash
                );
            }
            for binding in &binder.studies {
                let _ = writeln!(
                    out,
                    "study_file={} hash={}",
                    escape_text(&binding.path.to_string_lossy()),
                    binding.hash
                );
            }
            for (index, (verb, code, result)) in records.iter().enumerate() {
                let _ = writeln!(out, "step={} verb={verb} exit={code}", index + 1);
                for line in result.lines() {
                    let _ = writeln!(out, "  {line}");
                }
            }
            let _ = writeln!(out, "steps_planned={}", steps.len());
            out
        }
    };
    CommandOutput {
        exit_code: status,
        stdout,
        stderr,
    }
}
