//! `report` verb: export the engineering report a completed solve run sealed.
//!
//! The report is produced by the solve driver's seventh stage and retained in
//! the ledger as three artifacts (HTML, JSON twin, evidence package). This verb
//! never renders anything itself: it locates the run through the same loader
//! `solve --resume` uses, copies the exact retained bytes next to the ledger,
//! and prints their content hashes. An export therefore cannot disagree with
//! what the solver sealed, and an incomplete run cannot mint a report.
//!
//! 2026-08-25: the previous body of this file fabricated report values; do not
//! restore them. Every value shown to a user comes from a retained receipt.

use std::path::{Path, PathBuf};

use fs_ledger::Ledger;

use crate::json_read::JsonValue;
use crate::solve::{CompletedRunExport, SolveRefusal, load_completed_run};
use crate::{
    CommandOutput, Diagnostic, OutputMode, RESULT_SCHEMA, escape_text, exit, push_json_string,
    refusal,
};

/// A completed run located in a ledger, with its parsed report receipt and the
/// directory exports are written into.
pub(crate) struct LoadedExport {
    pub(crate) export: CompletedRunExport,
    pub(crate) receipt: JsonValue,
    pub(crate) dir: PathBuf,
}

/// What `report` exported.
pub(crate) struct ReportExport {
    pub(crate) html_path: PathBuf,
    pub(crate) json_path: PathBuf,
    pub(crate) content_hash: String,
    pub(crate) verdict: String,
    pub(crate) stages_completed: usize,
}

/// A typed refusal for an export verb.
pub(crate) fn refuse(
    mode: OutputMode,
    command: &'static str,
    exit_code: u8,
    code: &str,
    subject: &str,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> CommandOutput {
    let diagnostic = Diagnostic::new(command, code.to_string(), what.into(), fix.into())
        .with_subject(subject.to_string());
    refusal(
        mode,
        exit_code,
        &diagnostic,
        Some((command, "refused", subject, None, 0)),
    )
}

/// Render a solve-layer refusal (unknown run, incomplete run, ledger fault)
/// under the export verb's name, keeping the stable code and exit class.
pub(crate) fn refuse_solve(
    mode: OutputMode,
    command: &'static str,
    subject: &str,
    solve_refusal: &SolveRefusal,
) -> CommandOutput {
    let (status, exit_code) = if solve_refusal.dependency.is_some() {
        ("unavailable", exit::UNAVAILABLE)
    } else {
        ("refused", exit::REFUSED)
    };
    let mut what = solve_refusal.what.clone();
    if let Some(dependency) = solve_refusal.dependency {
        what.push_str(&format!(" (owning bead: {dependency})"));
    }
    let diagnostic = Diagnostic::new(
        command,
        solve_refusal.code.to_string(),
        what,
        solve_refusal.fix.clone(),
    )
    .with_subject(subject.to_string());
    refusal(
        mode,
        exit_code,
        &diagnostic,
        Some((command, status, subject, None, 0)),
    )
}

/// Open the ledger read path and locate the completed run.
///
/// Refuses (never creates a ledger) when the path is missing, and refuses
/// through the solve loader's own codes when the run is unknown or incomplete.
pub(crate) fn load_export(
    command: &'static str,
    run_id: &str,
    ledger: Option<&Path>,
    mode: OutputMode,
) -> Result<LoadedExport, CommandOutput> {
    let Some(ledger_path) = ledger else {
        return Err(refuse(
            mode,
            command,
            exit::INPUT,
            "cli-export-ledger-required",
            run_id,
            "no ledger path was given",
            format!("pass the ledger the solve wrote into: `frankensim {command} <run-id> <ledger.db>`"),
        ));
    };
    if !ledger_path.is_file() {
        return Err(refuse(
            mode,
            command,
            exit::INPUT,
            "cli-export-ledger-missing",
            run_id,
            format!("ledger `{}` does not exist", ledger_path.display()),
            "pass the ledger the solve wrote into; exports never create a ledger",
        ));
    }
    let Some(path_str) = ledger_path.to_str() else {
        return Err(refuse(
            mode,
            command,
            exit::INPUT,
            "cli-export-ledger-path",
            run_id,
            "the ledger path is not valid UTF-8",
            "pass a UTF-8 ledger path",
        ));
    };
    let opened = Ledger::open(path_str).map_err(|error| {
        refuse(
            mode,
            command,
            exit::INPUT,
            "cli-export-ledger-open",
            run_id,
            format!("cannot open ledger `{path_str}`: {error}"),
            "pass the ledger the solve wrote into",
        )
    })?;
    let export = load_completed_run(&opened, run_id)
        .map_err(|solve_refusal| refuse_solve(mode, command, run_id, &solve_refusal))?;
    let receipt = JsonValue::parse(&export.report_receipt).map_err(|error| {
        refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-receipt-shape",
            run_id,
            format!("the retained report receipt does not parse: {error}"),
            "regenerate the run; exports never repair a receipt",
        )
    })?;
    let dir = ledger_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
    Ok(LoadedExport {
        export,
        receipt,
        dir,
    })
}

/// Write retained bytes to `path`. An existing identical file is left alone
/// (exports are idempotent); an existing differing file is a conflict, never
/// overwritten.
pub(crate) fn write_retained(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if path.exists() {
        let existing = std::fs::read(path)
            .map_err(|error| format!("cannot read existing `{}`: {error}", path.display()))?;
        if existing == bytes {
            return Ok(());
        }
        return Err(format!(
            "`{}` exists and differs from the retained artifact",
            path.display()
        ));
    }
    std::fs::write(path, bytes)
        .map_err(|error| format!("cannot write `{}`: {error}", path.display()))
}

/// Export the retained report (HTML + JSON twin) for `run_id` next to the ledger.
pub(crate) fn export_report(
    command: &'static str,
    run_id: &str,
    ledger: Option<&Path>,
    mode: OutputMode,
) -> Result<ReportExport, CommandOutput> {
    let loaded = load_export(command, run_id, ledger, mode)?;
    let run = loaded.export.run.clone();
    let html_path = loaded.dir.join(format!("{run}.report.html"));
    let json_path = loaded.dir.join(format!("{run}.report.json"));
    for (path, bytes) in [
        (&html_path, &loaded.export.report_html),
        (&json_path, &loaded.export.report_json),
    ] {
        write_retained(path, bytes).map_err(|why| {
            refuse(
                mode,
                command,
                exit::REFUSED,
                "cli-export-output-conflict",
                &run,
                why,
                "move the existing file aside; exports never overwrite a differing artifact",
            )
        })?;
    }
    Ok(ReportExport {
        html_path,
        json_path,
        content_hash: loaded
            .receipt
            .str_field("report_content_hash")
            .unwrap_or_default()
            .to_string(),
        verdict: loaded
            .receipt
            .str_field("verdict")
            .unwrap_or_default()
            .to_string(),
        stages_completed: loaded.export.stages.len(),
    })
}

/// Execute the `report` verb.
#[must_use]
pub fn report_path(run_id: &str, ledger_path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    let export = match export_report("report", run_id, ledger_path, mode) {
        Ok(export) => export,
        Err(output) => return output,
    };
    let run = run_id;
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"report\",\"status\":\"ok\",\"subject\":");
            push_json_string(&mut out, run);
            out.push_str(",\"run\":");
            push_json_string(&mut out, run);
            out.push_str(",\"report_html\":");
            push_json_string(&mut out, &export.html_path.to_string_lossy());
            out.push_str(",\"report_json\":");
            push_json_string(&mut out, &export.json_path.to_string_lossy());
            out.push_str(",\"content_hash\":");
            push_json_string(&mut out, &export.content_hash);
            out.push_str(",\"verdict\":");
            push_json_string(&mut out, &export.verdict);
            out.push_str(&format!(
                ",\"stages_completed\":{},\"authority\":\"projection-of-retained-receipts\"}}\n",
                export.stages_completed
            ));
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=report\nsubject={run}\nrun={run}\nreport_html={}\nreport_json={}\ncontent_hash={}\nverdict={}\nstages_completed={}\nauthority=projection-of-retained-receipts\n",
            escape_text(&export.html_path.to_string_lossy()),
            escape_text(&export.json_path.to_string_lossy()),
            export.content_hash,
            escape_text(&export.verdict),
            export.stages_completed,
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}
