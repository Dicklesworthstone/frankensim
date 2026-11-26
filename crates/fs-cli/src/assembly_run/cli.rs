//! CLI handlers for `frankensim run <assembly>`, `frankensim run --resume`,
//! and `frankensim explain <run-id> [<ledger.db>]` (`.11.7` Stage 3).

use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use fs_blake3::{ContentHash, hash_bytes};
use fs_ir::ast::{CountUnit, Node, NodeKind};
use fs_ledger::Ledger;

use super::{
    AssemblyAnalysisMode, AssemblyCancelPlan, explain_assembly_run, resume_assembly, run_assembly,
};
use crate::cards::CardPackKind;
use crate::json_read::JsonValue;
use crate::solve::{SolveOutcome, SolveRefusal, SolveRunId, SolveRunStatus};
use crate::{
    CommandOutput, Diagnostic, MAX_PROJECT_BYTES, OutputMode, RESULT_SCHEMA, admit_card_packs,
    escape_text, exit, format_diagnostic, open_solve_ledger, push_json_string, refusal,
};

fn skip_sexpr_ws_and_comments(mut s: &str) -> &str {
    loop {
        s = s.trim_start();
        if let Some(rest) = s.strip_prefix(';') {
            s = rest.find('\n').map_or("", |i| &rest[i + 1..]);
        } else {
            return s;
        }
    }
}

fn push_kv(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, ",\"{key}\":");
    push_json_string(out, val);
}

/// Detect whether `path` is a canonical `(assembly ...)` `.fsim` or isomorphic
/// `fs-ir` `.json` assembly description.
pub(crate) fn is_assembly_file(path: &Path) -> bool {
    let ext = path.extension().and_then(|e| e.to_str());
    if !matches!(ext, Some("fsim" | "json")) {
        return false;
    }
    let Ok(meta) = std::fs::metadata(path) else {
        return false;
    };
    if !meta.is_file() || meta.len() > MAX_PROJECT_BYTES {
        return false;
    }
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let mut buf = Vec::new();
    if file
        .take(MAX_PROJECT_BYTES.min(65_536))
        .read_to_end(&mut buf)
        .is_err()
    {
        return false;
    }
    let text = match std::str::from_utf8(&buf) {
        Ok(s) => s,
        Err(e) if e.valid_up_to() > 0 => match std::str::from_utf8(&buf[..e.valid_up_to()]) {
            Ok(s) => s,
            Err(_) => return false,
        },
        Err(_) => return false,
    };
    let trimmed = skip_sexpr_ws_and_comments(text);
    if let Some(rest) = trimmed
        .strip_prefix('(')
        .map(str::trim_start)
        .and_then(|s| s.strip_prefix("assembly"))
    {
        return rest.is_empty() || rest.starts_with(|c: char| c.is_ascii_whitespace() || c == '"');
    }
    ext == Some("json")
        && trimmed.starts_with(['{', '['])
        && std::fs::read_to_string(path)
            .ok()
            .and_then(|full| fs_ir::json::parse(&full).ok())
            .is_some_and(|ast| ast.head() == Some("assembly"))
}

/// Bounded read of an assembly `.fsim` or `.json` file.
pub(crate) fn read_assembly_bytes(
    path: &Path,
    command: &'static str,
    mode: OutputMode,
) -> Result<Vec<u8>, CommandOutput> {
    let label = path.to_string_lossy();
    let refuse = |code: &str, msg: String, fix: &str| {
        let d = Diagnostic::new(command, code, msg, fix).with_subject(label.to_string());
        refusal(
            mode,
            exit::INPUT,
            &d,
            Some((command, "refused", &label, None, 0)),
        )
    };
    let meta = std::fs::metadata(path).map_err(|e| {
        refuse(
            "cli-input-read",
            format!("cannot inspect assembly `{label}`: {e}"),
            "provide a readable regular UTF-8 assembly file",
        )
    })?;
    if !meta.is_file() || meta.len() > MAX_PROJECT_BYTES {
        return Err(refuse(
            "cli-input-too-large",
            format!("assembly `{label}` must be a regular file within {MAX_PROJECT_BYTES} bytes"),
            "provide a bounded canonical assembly file",
        ));
    }
    let file = std::fs::File::open(path).map_err(|e| {
        refuse(
            "cli-input-read",
            format!("cannot open assembly `{label}`: {e}"),
            "provide a readable regular UTF-8 assembly file",
        )
    })?;
    let mut bytes = Vec::new();
    file.take(MAX_PROJECT_BYTES.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|e| {
            refuse(
                "cli-input-read",
                format!("cannot read assembly `{label}`: {e}"),
                "provide a readable regular UTF-8 assembly file",
            )
        })?;
    if !u64::try_from(bytes.len()).is_ok_and(|n| n <= MAX_PROJECT_BYTES) {
        return Err(refuse(
            "cli-input-too-large",
            format!("assembly `{label}` exceeded {MAX_PROJECT_BYTES} bytes while being read"),
            "retry against a stable file inside the documented cap",
        ));
    }
    Ok(bytes)
}

/// Canonical assembly header metadata used by FrankenScript `assembly.run`.
pub(crate) struct AssemblyHeader {
    pub(crate) source_hash: ContentHash,
    pub(crate) seed: Option<u64>,
    pub(crate) wall_budget_s: Option<f64>,
    pub(crate) mem_budget_bytes: Option<u64>,
}

const TIME_DIMS: fs_qty::Dims = fs_qty::Dims([0, 0, 1, 0, 0, 0]);

fn mem_bytes_of(node: &Node) -> Option<u64> {
    let NodeKind::Count { value, unit } = node.kind else {
        return None;
    };
    if unit == CountUnit::Cores {
        return None;
    }
    value.integral_bytes(unit)
}

pub(crate) fn parse_assembly_header(source: &[u8]) -> Result<AssemblyHeader, String> {
    let text = std::str::from_utf8(source).map_err(|e| format!("non-UTF-8 source: {e}"))?;
    let ast = if text.trim_start().starts_with(['{', '[']) {
        fs_ir::json::parse(text)
    } else {
        fs_ir::sexpr::parse(text)
    }
    .map_err(|e| format!("assembly parse failed: {e}"))?;
    let canonical =
        fs_ir::sexpr::print_checked(&ast).map_err(|e| format!("canonicalization failed: {e}"))?;
    let items = ast
        .items()
        .ok_or_else(|| "assembly must be an s-expression list".to_string())?;
    if ast.head() != Some("assembly") || items.len() < 3 {
        return Err("top-level form must be (assembly \"name\" ...)".to_string());
    }
    let (mut seed, mut wall_budget_s, mut mem_budget_bytes) = (None, None, None);
    for clause in &items[2..] {
        let Some(sub) = clause.items() else { continue };
        match clause.head() {
            Some("seed") if sub.len() == 2 => {
                if let NodeKind::Seed(s) = sub[1].kind {
                    seed = Some(s);
                }
            }
            Some("budget") => {
                for entry in &sub[1..] {
                    let Some(pair) = entry.items() else { continue };
                    match (entry.head(), pair.get(1).map(|n| &n.kind)) {
                        (Some("wall"), Some(NodeKind::Qty { value, dims, .. }))
                            if *dims == TIME_DIMS && value.is_finite() =>
                        {
                            wall_budget_s = Some(*value);
                        }
                        (Some("mem"), Some(_)) => mem_budget_bytes = mem_bytes_of(&pair[1]),
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
    Ok(AssemblyHeader {
        source_hash: hash_bytes(canonical.as_bytes()),
        seed,
        wall_budget_s,
        mem_budget_bytes,
    })
}

fn render_assembly_refusal(
    command: &'static str,
    mode: OutputMode,
    subject: &str,
    err: &SolveRefusal,
) -> CommandOutput {
    let (status, exit_code) = if err.dependency.is_some() {
        ("unavailable", exit::UNAVAILABLE)
    } else if err.code == "cli-solve-cancelled" {
        ("cancelled", exit::CANCELLED)
    } else if err.code == "cli-solve-budget" {
        ("budget_exceeded", exit::BUDGET)
    } else {
        ("refused", exit::REFUSED)
    };
    let mut message = err.what.clone();
    if let Some(op) = err.recorded_op {
        let _ = write!(message, "; refusal retained as ledger operation {op}");
    }
    let diag = Diagnostic::new(command, err.code, message, err.fix.clone()).with_subject(subject);
    let stderr = format_diagnostic(mode, &diag);
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", command);
            push_kv(&mut out, "status", status);
            push_kv(&mut out, "subject", subject);
            if let Some(run) = &err.run {
                push_kv(&mut out, "run", run);
            }
            if let Some(stage) = err.stage {
                push_kv(&mut out, "stage", stage);
            }
            out.push_str("}\n");
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status={status}\ncommand={command}\nsubject={}\n",
                escape_text(subject)
            );
            if let Some(run) = &err.run {
                let _ = writeln!(out, "run={run}");
            }
            if let Some(stage) = err.stage {
                let _ = writeln!(out, "stage={stage}");
            }
            out
        }
    };
    CommandOutput {
        exit_code,
        stdout,
        stderr,
    }
}

fn read_report_fields(
    ledger: &Ledger,
    outcome: &SolveOutcome,
) -> (Option<String>, Option<String>, Option<String>) {
    let Some(rep) = outcome.stages.iter().find(|s| s.stage == "report") else {
        return (None, None, None);
    };
    let Some(val) = ContentHash::from_hex(&rep.receipt)
        .and_then(|h| ledger.get_artifact(&h).ok().flatten())
        .and_then(|b| String::from_utf8(b).ok())
        .and_then(|t| JsonValue::parse(&t).ok())
    else {
        return (None, None, None);
    };
    (
        val.str_field("plan_hash").map(str::to_string),
        val.str_field("report_content_hash").map(str::to_string),
        val.str_field("package_root").map(str::to_string),
    )
}

fn render_assembly_outcome(
    ledger: &Ledger,
    mode: OutputMode,
    subject: &str,
    analysis: Option<AssemblyAnalysisMode>,
    outcome: &SolveOutcome,
    progress: &[String],
) -> CommandOutput {
    let mut stderr = String::new();
    for line in progress {
        stderr.push_str(line);
        stderr.push('\n');
    }
    let (status, exit_code) = match outcome.status {
        SolveRunStatus::Completed => ("completed", exit::SUCCESS),
        SolveRunStatus::BudgetExceeded { .. } => ("budget_exceeded", exit::BUDGET),
        SolveRunStatus::Cancelled => ("cancelled", exit::CANCELLED),
    };
    let stages_completed = outcome.prior_stages as usize + outcome.stages.len();
    let (plan_hash, report_hash, pkg_root) = read_report_fields(ledger, outcome);
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "run");
            push_kv(&mut out, "kind", "assembly");
            push_kv(&mut out, "status", status);
            push_kv(&mut out, "subject", subject);
            if let Some(a) = analysis {
                push_kv(&mut out, "analysis", a.as_str());
            }
            push_kv(&mut out, "run", &outcome.run);
            let _ = write!(out, ",\"stages_completed\":{stages_completed}");
            if let Some(ph) = &plan_hash {
                push_kv(&mut out, "plan_hash", ph);
            }
            if let Some(rh) = &report_hash {
                push_kv(&mut out, "report_content_hash", rh);
            }
            if let Some(pr) = &pkg_root {
                push_kv(&mut out, "package_root", pr);
                push_kv(&mut out, "checker", "pass");
            }
            out.push_str("}\n");
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status={status}\ncommand=run\nkind=assembly\nsubject={}\nrun={}\nstages_completed={stages_completed}\n",
                escape_text(subject),
                outcome.run,
            );
            if let Some(a) = analysis {
                let _ = writeln!(out, "analysis={}", a.as_str());
            }
            if let Some(ph) = &plan_hash {
                let _ = writeln!(out, "plan_hash={ph}");
            }
            if let Some(rh) = &report_hash {
                let _ = writeln!(out, "report_content_hash={rh}");
            }
            if let Some(pr) = &pkg_root {
                let _ = writeln!(out, "package_root={pr}\nchecker=pass");
            }
            out
        }
    };
    CommandOutput {
        exit_code,
        stdout,
        stderr,
    }
}

/// Execute `frankensim run <assembly.fsim|assembly.json> <ledger.db> ...`
/// against already-read canonical bytes (used both by CLI and FrankenScript so
/// the exact hash-pinned bytes execute with no TOCTOU re-read).
pub(crate) fn run_assembly_bytes(
    path: &Path,
    bytes: &[u8],
    ledger_path: &Path,
    declared_cards: &[(CardPackKind, PathBuf)],
    analysis: AssemblyAnalysisMode,
    mode: OutputMode,
) -> CommandOutput {
    let label = path.to_string_lossy();
    let packs = match admit_card_packs(declared_cards, mode) {
        Ok(p) => p,
        Err(out) => return out,
    };
    let ledger = match open_solve_ledger(ledger_path, &label, mode) {
        Ok(l) => l,
        Err(out) => return out,
    };
    let started = std::time::Instant::now();
    let mut clock = move || started.elapsed().as_secs_f64();
    match run_assembly(
        &ledger,
        bytes,
        &packs,
        analysis,
        &mut clock,
        AssemblyCancelPlan::default(),
    ) {
        Ok((outcome, progress)) => {
            render_assembly_outcome(&ledger, mode, &label, Some(analysis), &outcome, &progress)
        }
        Err(err) => render_assembly_refusal("run", mode, &label, &err),
    }
}

/// Execute `frankensim run <assembly.fsim|assembly.json> <ledger.db> ...`.
pub(crate) fn run_assembly_path(
    path: &Path,
    ledger_path: &Path,
    declared_cards: &[(CardPackKind, PathBuf)],
    analysis: AssemblyAnalysisMode,
    mode: OutputMode,
) -> CommandOutput {
    let bytes = match read_assembly_bytes(path, "run", mode) {
        Ok(b) => b,
        Err(out) => return out,
    };
    run_assembly_bytes(path, &bytes, ledger_path, declared_cards, analysis, mode)
}

/// Execute `frankensim run --resume <run-id> <ledger.db>`.
pub(crate) fn resume_assembly_path(
    run_id: &str,
    ledger_path: &Path,
    mode: OutputMode,
) -> CommandOutput {
    let Some(parsed) = SolveRunId::parse_hex(run_id) else {
        let diag = Diagnostic::new(
            "run",
            "cli-run-resume-id",
            format!("invalid assembly run id `{run_id}`"),
            "pass a 64-character lowercase hex run id",
        )
        .with_subject(run_id);
        return refusal(
            mode,
            exit::USAGE,
            &diag,
            Some(("run", "refused", run_id, None, 1)),
        );
    };
    let ledger = match open_solve_ledger(ledger_path, run_id, mode) {
        Ok(l) => l,
        Err(out) => return out,
    };
    let started = std::time::Instant::now();
    let mut clock = move || started.elapsed().as_secs_f64();
    match resume_assembly(&ledger, parsed, &mut clock, AssemblyCancelPlan::default()) {
        Ok((outcome, progress)) => {
            render_assembly_outcome(&ledger, mode, run_id, None, &outcome, &progress)
        }
        Err(err) => render_assembly_refusal("run", mode, run_id, &err),
    }
}

fn render_explanation(run_id: &str, exp: &super::AssemblyExplanation, mode: OutputMode) -> String {
    match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            push_kv(&mut out, "command", "explain");
            push_kv(&mut out, "status", "ok");
            push_kv(&mut out, "subject", run_id);
            push_kv(&mut out, "run", &exp.run);
            let _ = write!(out, ",\"stages_completed\":{}", exp.stages_completed);
            match &exp.plan_hash {
                Some(ph) => push_kv(&mut out, "plan_hash", ph),
                None => out.push_str(",\"plan_hash\":null"),
            }
            let _ = write!(
                out,
                ",\"lowering_receipt\":{},\"execution_plan\":{},\"causalization_audit\":{},\"window_ledger\":{},\"probe_series\":{},\"qoi_budget\":{},\"stages\":[",
                exp.lowering_receipt_json
                    .as_deref()
                    .map_or("null", str::trim),
                exp.plan_json.as_deref().map_or("null", str::trim),
                exp.causalization_audit_json
                    .as_deref()
                    .map_or("null", str::trim),
                exp.window_ledger_json.as_deref().map_or("null", str::trim),
                exp.probe_series_json.as_deref().map_or("null", str::trim),
                exp.qoi_budget_json.as_deref().map_or("null", str::trim),
            );
            for (i, (stage, op_id, receipt)) in exp.stages.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                out.push_str("{\"stage\":");
                push_json_string(&mut out, stage);
                let _ = write!(out, ",\"op_id\":{op_id}");
                push_kv(&mut out, "receipt", receipt);
                out.push('}');
            }
            out.push_str("]}\n");
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status=ok\ncommand=explain\nsubject={}\nrun={}\nstages_completed={}\n",
                escape_text(run_id),
                exp.run,
                exp.stages_completed,
            );
            if let Some(ph) = &exp.plan_hash {
                let _ = writeln!(out, "plan_hash={ph}");
            }
            let _ = writeln!(
                out,
                "lowering_receipt_present={}\nexecution_plan_present={}\ncausalization_audit_present={}\nwindow_ledger_present={}\nprobe_series_present={}\nqoi_budget_present={}",
                exp.lowering_receipt_json.is_some(),
                exp.plan_json.is_some(),
                exp.causalization_audit_json.is_some(),
                exp.window_ledger_json.is_some(),
                exp.probe_series_json.is_some(),
                exp.qoi_budget_json.is_some(),
            );
            if let Some(raw) = &exp.qoi_budget_json
                && let Ok(qb) = super::budget::EvaluatedQoiBudget::from_json(raw)
            {
                for r in &qb.requirements {
                    let unm = r.unmeasured_terms.join(",");
                    let _ = writeln!(
                        out,
                        "requirement={} qoi={} outcome={} limit={:.6e} margin={:.6e} unit={} unmeasured=[{unm}]",
                        r.id, r.qoi, r.outcome, r.effective_limit, r.nominal_margin, r.unit
                    );
                }
                for t in &qb.terms {
                    let val = t
                        .value
                        .map_or_else(|| "NO-DATA".into(), |v| format!("{v:.6e}"));
                    let _ = writeln!(
                        out,
                        "budget_term qoi={} kind={} state={} value={val} reason={}",
                        t.qoi,
                        t.kind,
                        t.state,
                        escape_text(&t.reason)
                    );
                }
            }
            for (stage, op_id, receipt) in &exp.stages {
                let _ = writeln!(out, "stage={stage} op_id={op_id} receipt={receipt}");
            }
            out
        }
    }
}

/// Execute `frankensim explain <run-id> [<ledger.db>]`.
pub(crate) fn explain_assembly_path(
    run_id: &str,
    ledger_path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let err_out = |code: &str, exit_code: u8, msg: String, fix: &str, errs: usize| {
        let d = Diagnostic::new("explain", code, msg, fix).with_subject(run_id);
        refusal(
            mode,
            exit_code,
            &d,
            Some(("explain", "refused", run_id, None, errs)),
        )
    };
    let Some(path) = ledger_path else {
        return err_out(
            "cli-export-ledger-required",
            exit::INPUT,
            "no ledger path was given".to_string(),
            "pass the ledger the run wrote into: `frankensim explain <run-id> <ledger.db>`",
            0,
        );
    };
    if !path.is_file() {
        return err_out(
            "cli-explain-ledger-missing",
            exit::INPUT,
            format!("ledger `{}` does not exist", path.display()),
            "pass an existing ledger file written by `frankensim run`",
            0,
        );
    }
    let Some(parsed) = SolveRunId::parse_hex(run_id) else {
        return err_out(
            "cli-explain-run-id",
            exit::USAGE,
            format!("invalid assembly run id `{run_id}`"),
            "pass a 64-character lowercase hex run id",
            1,
        );
    };
    let ledger = match open_solve_ledger(path, run_id, mode) {
        Ok(l) => l,
        Err(out) => return out,
    };
    let exp = match explain_assembly_run(&ledger, parsed) {
        Ok(e) => e,
        Err(err) => return render_assembly_refusal("explain", mode, run_id, &err),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout: render_explanation(run_id, &exp, mode),
        stderr: String::new(),
    }
}
