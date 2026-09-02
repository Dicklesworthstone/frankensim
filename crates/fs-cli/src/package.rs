//! `package` verb: export the evidence package a completed solve run sealed,
//! re-parse it, and re-run the solver-free checker on the exact bytes.
//!
//! The package is built by the solve driver's report stage from retained
//! receipts and retained as a ledger artifact. This verb copies those bytes
//! next to the ledger, proves they still parse as a format-9 package, and
//! prints the checker verdict and Merkle root. It mints no claim of its own.
//!
//! 2026-08-25: the previous body of this file fabricated package claims; do
//! not restore them.

use std::path::{Path, PathBuf};

use fs_package::EvidencePackage;

use crate::report::{load_export, refuse, write_retained};
use crate::{CommandOutput, OutputMode, RESULT_SCHEMA, escape_text, exit, push_json_string};

/// What `package` exported.
pub(crate) struct PackageExport {
    pub(crate) path: PathBuf,
    pub(crate) merkle_root: String,
    pub(crate) claim_count: usize,
}

/// Export the retained evidence package for `run_id` next to the ledger and
/// re-check it with the solver-free checker.
pub(crate) fn export_package(
    command: &'static str,
    run_id: &str,
    ledger: Option<&Path>,
    mode: OutputMode,
) -> Result<PackageExport, CommandOutput> {
    let loaded = load_export(command, run_id, ledger, mode)?;
    let run = loaded.export.run.clone();
    let text = std::str::from_utf8(&loaded.export.package_json).map_err(|_| {
        refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-package-parse",
            &run,
            "the retained package is not UTF-8",
            "regenerate the run; exports never repair a package",
        )
    })?;
    let package = EvidencePackage::from_json(text).map_err(|error| {
        refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-package-parse",
            &run,
            format!("the retained package does not parse: {error}"),
            "regenerate the run; exports never repair a package",
        )
    })?;
    let check = fs_checker::check(&package);
    if !check.passed() {
        return Err(refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-package-check",
            &run,
            "the solver-free checker refused the retained package",
            "regenerate the run and report the packaging defect; a refused package is never exported",
        ));
    }
    let merkle_root = check.merkle_root().to_hex();
    if let Some(recorded) = loaded.receipt.str_field("package_root")
        && recorded != merkle_root
    {
        return Err(refuse(
            mode,
            command,
            exit::REFUSED,
            "cli-export-package-check",
            &run,
            format!(
                "the report receipt records package root {recorded} but the retained bytes hash to {merkle_root}"
            ),
            "regenerate the run; the receipt and the retained package disagree",
        ));
    }
    let path = loaded.dir.join(format!("{run}.fspkg"));
    write_retained(&path, &loaded.export.package_json).map_err(|why| {
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
    Ok(PackageExport {
        path,
        merkle_root,
        claim_count: package.declared_claims_unverified().len(),
    })
}

/// Execute the `package` verb.
#[must_use]
pub fn package_path(run_id: &str, ledger_path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    let export = match export_package("package", run_id, ledger_path, mode) {
        Ok(export) => export,
        Err(output) => return output,
    };
    let stdout = match mode {
        OutputMode::Json => {
            let mut out = String::from("{\"schema\":");
            push_json_string(&mut out, RESULT_SCHEMA);
            out.push_str(",\"command\":\"package\",\"status\":\"ok\",\"subject\":");
            push_json_string(&mut out, run_id);
            out.push_str(",\"run\":");
            push_json_string(&mut out, run_id);
            out.push_str(",\"package\":");
            push_json_string(&mut out, &export.path.to_string_lossy());
            out.push_str(",\"merkle_root\":");
            push_json_string(&mut out, &export.merkle_root);
            out.push_str(&format!(
                ",\"claim_count\":{},\"checker\":\"pass\",\"authority\":\"structural-integrity-only\"}}\n",
                export.claim_count
            ));
            out
        }
        OutputMode::Text => format!(
            "status=ok\ncommand=package\nsubject={run_id}\nrun={run_id}\npackage={}\nmerkle_root={}\nclaim_count={}\nchecker=pass\nauthority=structural-integrity-only\n",
            escape_text(&export.path.to_string_lossy()),
            export.merkle_root,
            export.claim_count,
        ),
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}
