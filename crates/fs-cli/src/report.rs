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

use std::io::Write;
use std::path::{Path, PathBuf};

use fs_blake3::{hash_bytes, hash_domain};
use fs_ledger::Ledger;
use fs_package::EvidencePackage;

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
    /// How the retained run was proven before export; exports never replay
    /// physics, they re-emit sealed evidence and say so.
    pub(crate) verification: &'static str,
}

const REPORT_RECEIPT_SCHEMA: &str = "frankensim.cli.solve-report.v1";
const REPORT_CONTENT_HASH_DOMAIN: &str = "org.frankensim.report.engineering.v1";

fn receipt_string<'a>(receipt: &'a JsonValue, key: &str) -> Result<&'a str, String> {
    receipt
        .str_field(key)
        .ok_or_else(|| format!("report receipt is missing string field `{key}`"))
}

fn source_receipt<'a>(receipt: &'a JsonValue, key: &str) -> Result<&'a str, String> {
    receipt
        .get("sources")
        .and_then(|sources| sources.str_field(key))
        .ok_or_else(|| format!("report receipt is missing sources.{key}"))
}

fn stage_receipt<'a>(export: &'a CompletedRunExport, stage: &str) -> Result<&'a str, String> {
    let mut receipts = export
        .stages
        .iter()
        .filter(|(name, _, _)| *name == stage)
        .map(|(_, _, receipt)| receipt.as_str());
    let Some(receipt) = receipts.next() else {
        return Err(format!("loaded run has no `{stage}` stage receipt"));
    };
    if receipts.next().is_some() {
        return Err(format!("loaded run has multiple `{stage}` stage receipts"));
    }
    Ok(receipt)
}

/// Bind the receipt's self-described identities to the already re-attested
/// completed run before an export can write any retained byte.
fn validate_export_receipt(export: &CompletedRunExport, receipt: &JsonValue) -> Result<(), String> {
    let html_hash = hash_bytes(&export.report_html).to_hex();
    let json_hash = hash_bytes(&export.report_json).to_hex();
    let package_hash = hash_bytes(&export.package_json).to_hex();
    for (field, expected) in [
        ("schema", REPORT_RECEIPT_SCHEMA),
        ("run", export.run.as_str()),
        ("stage", "report"),
        ("project_hash", export.project_hash.as_str()),
        ("report_html", html_hash.as_str()),
        ("report_json", json_hash.as_str()),
        ("package", package_hash.as_str()),
    ] {
        if receipt_string(receipt, field)? != expected {
            return Err(format!(
                "report receipt `{field}` does not match the loaded run"
            ));
        }
    }
    let report_content_hash = hash_domain(REPORT_CONTENT_HASH_DOMAIN, &export.report_html).to_hex();
    if receipt_string(receipt, "report_content_hash")? != report_content_hash.as_str() {
        return Err(
            "report receipt `report_content_hash` does not match the loaded report".to_string(),
        );
    }
    let report_text = std::str::from_utf8(&export.report_json)
        .map_err(|_| "loaded report JSON is not UTF-8".to_string())?;
    let report = JsonValue::parse(report_text)
        .map_err(|error| format!("loaded report JSON does not parse: {error}"))?;
    let verdict = report
        .get("requirements")
        .and_then(JsonValue::as_array)
        .and_then(<[JsonValue]>::first)
        .and_then(|requirement| requirement.str_field("outcome"))
        .ok_or_else(|| "loaded report JSON is missing requirements[0].outcome".to_string())?;
    if receipt_string(receipt, "verdict")? != verdict {
        return Err("report receipt `verdict` does not match the loaded report JSON".to_string());
    }
    let checker = receipt
        .get("checker")
        .ok_or_else(|| "report receipt is missing object field `checker`".to_string())?;
    if !matches!(checker.get("passed"), Some(JsonValue::Bool(true))) {
        return Err("report receipt checker did not record a passing decision".to_string());
    }
    if checker.f64_field("protocol") != Some(f64::from(fs_checker::CHECKER_PROTOCOL_VERSION)) {
        return Err("report receipt checker protocol is not current".to_string());
    }
    let package_text = std::str::from_utf8(&export.package_json)
        .map_err(|_| "loaded evidence package is not UTF-8".to_string())?;
    let package = EvidencePackage::from_json(package_text)
        .map_err(|error| format!("loaded evidence package does not parse: {error}"))?;
    let package_root = package
        .try_merkle_root()
        .map_err(|error| format!("loaded evidence package does not seal: {error}"))?
        .to_hex();
    if receipt_string(receipt, "package_root")? != package_root {
        return Err("report receipt `package_root` does not match the loaded package".to_string());
    }
    for (field, stage) in [
        ("qoi_receipt", "qoi"),
        ("conduction_receipt", "conduction"),
        ("material_receipt", "material-resolve"),
    ] {
        if source_receipt(receipt, field)? != stage_receipt(export, stage)? {
            return Err(format!(
                "report receipt sources.{field} does not match the loaded `{stage}` receipt"
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        CompletedRunExport, EvidencePackage, JsonValue, REPORT_CONTENT_HASH_DOMAIN,
        REPORT_RECEIPT_SCHEMA, hash_bytes, hash_domain, validate_export_receipt,
    };

    fn receipt_fixture() -> (CompletedRunExport, String) {
        let package = EvidencePackage::new(fs_package::Provenance::new("test", "test"));
        let package_json = package.to_json().expect("fixture package serializes");
        let export = CompletedRunExport {
            verification: "sealed-evidence",
            run: "1".repeat(64),
            project_hash: "2".repeat(64),
            stages: vec![
                ("import-verify", 1, "3".repeat(64)),
                ("assign", 2, "4".repeat(64)),
                ("material-resolve", 3, "5".repeat(64)),
                ("flow-network", 4, "6".repeat(64)),
                ("conduction", 5, "7".repeat(64)),
                ("qoi", 6, "8".repeat(64)),
                ("report", 7, "9".repeat(64)),
            ],
            report_receipt: String::new(),
            report_html: b"<html>fixture</html>".to_vec(),
            report_json: b"{\"requirements\":[{\"outcome\":\"indeterminate\"}]}".to_vec(),
            package_json: package_json.into_bytes(),
        };
        let package_root = EvidencePackage::from_json(
            std::str::from_utf8(&export.package_json).expect("fixture package is UTF-8"),
        )
        .expect("fixture package parses")
        .try_merkle_root()
        .expect("fixture package seals")
        .to_hex();
        let receipt = format!(
            "{{\"schema\":\"{REPORT_RECEIPT_SCHEMA}\",\"run\":\"{}\",\"stage\":\"report\",\"project_hash\":\"{}\",\"report_html\":\"{}\",\"report_json\":\"{}\",\"report_content_hash\":\"{}\",\"package\":\"{}\",\"package_root\":\"{package_root}\",\"checker\":{{\"passed\":true,\"protocol\":{}}},\"verdict\":\"indeterminate\",\"sources\":{{\"qoi_receipt\":\"{}\",\"conduction_receipt\":\"{}\",\"material_receipt\":\"{}\"}}}}",
            export.run,
            export.project_hash,
            hash_bytes(&export.report_html).to_hex(),
            hash_bytes(&export.report_json).to_hex(),
            hash_domain(REPORT_CONTENT_HASH_DOMAIN, &export.report_html).to_hex(),
            hash_bytes(&export.package_json).to_hex(),
            fs_checker::CHECKER_PROTOCOL_VERSION,
            export.stages[5].2,
            export.stages[4].2,
            export.stages[2].2,
        );
        (export, receipt)
    }

    #[test]
    fn receipt_binding_refuses_every_export_identity_substitution() {
        let (export, receipt) = receipt_fixture();
        let parsed = JsonValue::parse(&receipt).expect("fixture receipt parses");
        assert!(validate_export_receipt(&export, &parsed).is_ok());

        let html_hash = hash_bytes(&export.report_html).to_hex();
        let json_hash = hash_bytes(&export.report_json).to_hex();
        let content_hash = hash_domain(REPORT_CONTENT_HASH_DOMAIN, &export.report_html).to_hex();
        let package_hash = hash_bytes(&export.package_json).to_hex();
        let substitutions = vec![
            (
                format!("\"schema\":\"{REPORT_RECEIPT_SCHEMA}\""),
                "\"schema\":\"wrong.schema\"".to_string(),
            ),
            (
                format!("\"run\":\"{}\"", export.run),
                "\"run\":\"a\"".to_string(),
            ),
            (
                "\"stage\":\"report\"".to_string(),
                "\"stage\":\"wrong-stage\"".to_string(),
            ),
            (
                format!("\"project_hash\":\"{}\"", export.project_hash),
                "\"project_hash\":\"b\"".to_string(),
            ),
            (
                format!("\"report_html\":\"{html_hash}\""),
                "\"report_html\":\"c\"".to_string(),
            ),
            (
                format!("\"report_json\":\"{json_hash}\""),
                "\"report_json\":\"d\"".to_string(),
            ),
            (
                format!("\"report_content_hash\":\"{content_hash}\""),
                "\"report_content_hash\":\"wrong\"".to_string(),
            ),
            (
                format!("\"package\":\"{package_hash}\""),
                "\"package\":\"e\"".to_string(),
            ),
            (
                format!("\"qoi_receipt\":\"{}\"", export.stages[5].2),
                "\"qoi_receipt\":\"f\"".to_string(),
            ),
            (
                format!("\"conduction_receipt\":\"{}\"", export.stages[4].2),
                "\"conduction_receipt\":\"g\"".to_string(),
            ),
            (
                format!("\"material_receipt\":\"{}\"", export.stages[2].2),
                "\"material_receipt\":\"h\"".to_string(),
            ),
            (
                "\"verdict\":\"indeterminate\"".to_string(),
                "\"verdict\":\"pass\"".to_string(),
            ),
            (
                "\"passed\":true".to_string(),
                "\"passed\":false".to_string(),
            ),
            (
                format!("\"protocol\":{}", fs_checker::CHECKER_PROTOCOL_VERSION),
                format!("\"protocol\":{}", fs_checker::CHECKER_PROTOCOL_VERSION + 1),
            ),
        ];
        for (expected, replacement) in substitutions {
            let hostile = receipt.replacen(&expected, &replacement, 1);
            let parsed = JsonValue::parse(&hostile).expect("hostile receipt still parses");
            assert!(
                validate_export_receipt(&export, &parsed).is_err(),
                "receipt substitution must refuse"
            );
        }
        for required in [
            format!("\"report_content_hash\":\"{content_hash}\","),
            "\"verdict\":\"indeterminate\",".to_string(),
            format!(
                "\"checker\":{{\"passed\":true,\"protocol\":{}}},",
                fs_checker::CHECKER_PROTOCOL_VERSION
            ),
        ] {
            let hostile = receipt.replacen(&required, "", 1);
            let parsed = JsonValue::parse(&hostile).expect("omitted-field receipt still parses");
            assert!(
                validate_export_receipt(&export, &parsed).is_err(),
                "required receipt field omission must refuse"
            );
        }

        let package_root = EvidencePackage::from_json(
            std::str::from_utf8(&export.package_json).expect("fixture package is UTF-8"),
        )
        .expect("fixture package parses")
        .try_merkle_root()
        .expect("fixture package seals")
        .to_hex();
        let hostile = receipt.replacen(&package_root, "i", 1);
        let parsed = JsonValue::parse(&hostile).expect("hostile package-root receipt parses");
        assert!(validate_export_receipt(&export, &parsed).is_err());
    }
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
    let (opened, dir) = open_export_ledger(command, run_id, ledger, mode)?;
    load_export_from(command, run_id, &opened, dir, mode)
}

/// Open the ledger read path for an export verb without locating a run.
///
/// Returns the opened ledger and the directory exports are written into.
/// Refuses (never creates a ledger) when the path is missing.
pub(crate) fn open_export_ledger(
    command: &'static str,
    run_id: &str,
    ledger: Option<&Path>,
    mode: OutputMode,
) -> Result<(Ledger, PathBuf), CommandOutput> {
    let Some(ledger_path) = ledger else {
        return Err(refuse(
            mode,
            command,
            exit::INPUT,
            "cli-export-ledger-required",
            run_id,
            "no ledger path was given",
            format!(
                "pass the ledger the solve wrote into: `frankensim {command} <run-id> <ledger.db>`"
            ),
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
    let trace = std::env::var_os("FS_CLI_TRACE_EXPORT").is_some();
    let t_open = std::time::Instant::now();
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
    if trace {
        eprintln!(
            "TRACE export: ledger open {:.3}s",
            t_open.elapsed().as_secs_f64()
        );
    }
    let dir = ledger_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
    Ok((opened, dir))
}

/// Locate the completed run in an already opened ledger and bind its report
/// receipt to the loaded artifacts.
pub(crate) fn load_export_from(
    command: &'static str,
    run_id: &str,
    opened: &Ledger,
    dir: PathBuf,
    mode: OutputMode,
) -> Result<LoadedExport, CommandOutput> {
    let trace = std::env::var_os("FS_CLI_TRACE_EXPORT").is_some();
    let t_load = std::time::Instant::now();
    let export = load_completed_run(opened, run_id)
        .map_err(|solve_refusal| refuse_solve(mode, command, run_id, &solve_refusal))?;
    if trace {
        eprintln!(
            "TRACE export: load_completed_run {:.3}s",
            t_load.elapsed().as_secs_f64()
        );
    }
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
    validate_export_receipt(&export, &receipt).map_err(|why| {
        refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-receipt-binding",
            run_id,
            why,
            "regenerate the run; exports never write bytes named by an unbound receipt",
        )
    })?;
    Ok(LoadedExport {
        export,
        receipt,
        dir,
    })
}

/// Check whether `path` is absent or already contains the retained bytes.
///
/// The boolean is true only for an existing identical file. A differing file
/// is always a conflict; read failures are not treated as absence.
fn retained_compatible(path: &Path, bytes: &[u8]) -> Result<bool, String> {
    match std::fs::read(path) {
        Ok(existing) if existing == bytes => Ok(true),
        Ok(_) => Err(format!(
            "`{}` exists and differs from the retained artifact",
            path.display()
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(error) => Err(format!(
            "cannot read existing `{}`: {error}",
            path.display()
        )),
    }
}

/// Write retained bytes to `path`. An existing identical file is left alone
/// (exports are idempotent); an existing differing file is a conflict, never
/// overwritten.
pub(crate) fn write_retained(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if retained_compatible(path, bytes)? {
        return Ok(());
    }
    match std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
    {
        Ok(mut file) => file
            .write_all(bytes)
            .map_err(|error| format!("cannot write `{}`: {error}", path.display())),
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            match retained_compatible(path, bytes)? {
                true => Ok(()),
                false => Err(format!(
                    "`{}` disappeared during retained-artifact conflict checking",
                    path.display()
                )),
            }
        }
        Err(error) => Err(format!("cannot write `{}`: {error}", path.display())),
    }
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
    let content_hash = receipt_string(&loaded.receipt, "report_content_hash")
        .map_err(|why| {
            refuse(
                mode,
                command,
                exit::REFUSED,
                "cli-export-receipt-binding",
                &run,
                why,
                "regenerate the run; exports never write bytes named by an unbound receipt",
            )
        })?
        .to_string();
    let verdict = receipt_string(&loaded.receipt, "verdict")
        .map_err(|why| {
            refuse(
                mode,
                command,
                exit::REFUSED,
                "cli-export-receipt-binding",
                &run,
                why,
                "regenerate the run; exports never write bytes named by an unbound receipt",
            )
        })?
        .to_string();
    let html_path = loaded.dir.join(format!("{run}.report.html"));
    let json_path = loaded.dir.join(format!("{run}.report.json"));
    // The two representations are one product. Refuse every conflict visible
    // at invocation time before publishing either path; write_retained uses
    // create_new so a concurrent creator is never overwritten. Two independent
    // filesystem paths are not a transactional pair: a concurrent creator can
    // still make the second path conflict after the first path is published.
    for (path, bytes) in [
        (&html_path, &loaded.export.report_html),
        (&json_path, &loaded.export.report_json),
    ] {
        retained_compatible(path, bytes).map_err(|why| {
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
        verification: loaded.export.verification,
        content_hash,
        verdict,
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
                ",\"stages_completed\":{},\"authority\":\"projection-of-retained-receipts\",\"verification\":\"{}\"}}\n",
                export.stages_completed, export.verification
            ));
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=report\nsubject={run}\nrun={run}\nreport_html={}\nreport_json={}\ncontent_hash={}\nverdict={}\nstages_completed={}\nauthority=projection-of-retained-receipts\nverification={}\n",
            escape_text(&export.html_path.to_string_lossy()),
            escape_text(&export.json_path.to_string_lossy()),
            export.content_hash,
            escape_text(&export.verdict),
            export.stages_completed,
            export.verification,
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}
