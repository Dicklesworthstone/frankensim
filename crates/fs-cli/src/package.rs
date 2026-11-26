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

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use fs_package::{
    CoverageStatus, EvidencePackage, PackageCoverageReport, Standard, package_coverage,
    verified_package_coverage,
};

use crate::report::{load_export, refuse, write_retained};
use crate::{
    CommandOutput, Diagnostic, OutputMode, RESULT_SCHEMA, escape_text, exit, is_operand,
    push_json_string,
};

/// What `package` exported.
pub(crate) struct PackageExport {
    pub(crate) path: PathBuf,
    pub(crate) merkle_root: String,
    pub(crate) claim_count: usize,
    pub(crate) package: EvidencePackage,
}

/// Parse a CLI `--coverage <standard>` slug matching `Standard::label()`.
pub(crate) fn parse_coverage_standard(slug: &str) -> Option<Standard> {
    match slug {
        "asme-vv-10" => Some(Standard::AsmeVvV10),
        "asme-vv-20" => Some(Standard::AsmeVvV20),
        "asme-vv-40" => Some(Standard::AsmeVvV40),
        "faa-easa-cba" => Some(Standard::FaaEasaCbA),
        _ => None,
    }
}

/// Parse `[<ledger.db>] [--coverage <standard>]` operands after `package <run-id>`.
pub(crate) fn parse_package_args(
    args: &[String],
) -> Result<(Option<PathBuf>, Option<Standard>), Diagnostic> {
    let (mut ledger, mut coverage) = (None, None);
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--coverage" {
            if coverage.is_some() {
                return Err(Diagnostic::new(
                    "package",
                    "cli-package-coverage-usage",
                    "`--coverage` may be stated at most once",
                    "pass a single standard: asme-vv-10, asme-vv-20, asme-vv-40, or faa-easa-cba",
                ));
            }
            let Some(slug) = args.get(i + 1).filter(|v| is_operand(v)) else {
                return Err(Diagnostic::new(
                    "package",
                    "cli-package-coverage-usage",
                    "`--coverage` requires one of asme-vv-10, asme-vv-20, asme-vv-40, or faa-easa-cba",
                    "pass `--coverage <standard>`",
                ));
            };
            let Some(parsed) = parse_coverage_standard(slug) else {
                return Err(Diagnostic::new(
                    "package",
                    "cli-package-coverage-usage",
                    format!("unknown coverage standard `{slug}`"),
                    "choose one of asme-vv-10, asme-vv-20, asme-vv-40, or faa-easa-cba",
                ));
            };
            coverage = Some(parsed);
            i += 2;
        } else if is_operand(&args[i]) && ledger.is_none() {
            ledger = Some(PathBuf::from(&args[i]));
            i += 1;
        } else {
            return Err(Diagnostic::new(
                "package",
                "cli-package-usage",
                format!("unexpected `package` argument `{}`", args[i]),
                "use `frankensim [--json] package <run-id> [<ledger.db>] [--coverage <standard>]`",
            ));
        }
    }
    Ok((ledger, coverage))
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
    let claim_count = package.declared_claims_unverified().len();
    Ok(PackageExport {
        path,
        merkle_root,
        claim_count,
        package,
    })
}

fn coverage_status_slug(status: &CoverageStatus) -> &'static str {
    match status {
        CoverageStatus::Covered => "covered",
        CoverageStatus::MappedButAbsent => "mapped-but-absent",
        CoverageStatus::NoClaim => "no-claim",
    }
}

fn coverage_counts(report: &PackageCoverageReport) -> (usize, usize, usize) {
    let (mut covered, mut absent, mut no_claim) = (0usize, 0usize, 0usize);
    for (_, status, _) in report.rows() {
        match status {
            CoverageStatus::Covered => covered += 1,
            CoverageStatus::MappedButAbsent => absent += 1,
            CoverageStatus::NoClaim => no_claim += 1,
        }
    }
    (covered, absent, no_claim)
}

fn push_coverage_json(out: &mut String, report: &PackageCoverageReport) {
    let (covered, absent, no_claim) = coverage_counts(report);
    out.push_str(",\"coverage\":{\"standard\":");
    push_json_string(out, report.standard().label());
    out.push_str(",\"standard_name\":");
    push_json_string(out, report.standard().full_name());
    let _ = write!(
        out,
        ",\"crosswalk_version\":{},\"package_format\":{},\"decision_hash\":",
        report.crosswalk_version(),
        report.package_format()
    );
    push_json_string(out, &report.decision_hash().to_hex());
    out.push_str(",\"verified_receipt\":");
    match report.receipt() {
        Some(r) => push_json_string(out, &r.receipt_hash().to_hex()),
        None => out.push_str("null"),
    }
    let _ = write!(
        out,
        ",\"summary\":{{\"covered\":{covered},\"mapped_but_absent\":{absent},\"no_claim\":{no_claim}}},\"rows\":["
    );
    for (i, (concept, status, why)) in report.rows().iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"concept\":");
        push_json_string(out, concept.label());
        out.push_str(",\"status\":");
        push_json_string(out, coverage_status_slug(status));
        out.push_str(",\"why\":");
        push_json_string(out, why);
        out.push('}');
    }
    out.push_str("]}");
}

fn push_coverage_text(out: &mut String, report: &PackageCoverageReport) {
    let (covered, absent, no_claim) = coverage_counts(report);
    let _ = writeln!(
        out,
        "coverage_standard={}\ncoverage_crosswalk_version={}\ncoverage_package_format={}\ncoverage_decision_hash={}\ncoverage_covered={covered}\ncoverage_mapped_but_absent={absent}\ncoverage_no_claim={no_claim}",
        report.standard().label(),
        report.crosswalk_version(),
        report.package_format(),
        report.decision_hash().to_hex()
    );
    for (concept, status, why) in report.rows() {
        let _ = writeln!(
            out,
            "coverage_row={} status={} why={}",
            concept.label(),
            coverage_status_slug(status),
            escape_text(why)
        );
    }
}

/// Execute the `package` verb with optional `--coverage <standard>` reporting.
#[must_use]
pub(crate) fn package_path(
    run_id: &str,
    ledger_path: Option<&Path>,
    coverage: Option<Standard>,
    mode: OutputMode,
) -> CommandOutput {
    if run_id.starts_with("study-") {
        if coverage.is_some() {
            return refuse(
                mode,
                "package",
                exit::REFUSED,
                "cli-package-study-coverage-unsupported",
                run_id,
                "standards coverage is wired for sealed solve EvidencePackage runs, not study-* exports",
                "omit `--coverage` when exporting a study-* package, or run `--coverage` on a sealed solve run",
            );
        }
        return crate::study::export("package", run_id, ledger_path, mode);
    }
    let export = match export_package("package", run_id, ledger_path, mode) {
        Ok(export) => export,
        Err(output) => return output,
    };
    let report = match coverage {
        Some(standard) => {
            let report = match export.package.clone().into_verified() {
                Ok(verified) => verified_package_coverage(&verified, standard),
                Err(_) => package_coverage(&export.package, standard),
            };
            if !report.validate_decision_hash() {
                return refuse(
                    mode,
                    "package",
                    exit::REFUSED,
                    "cli-package-coverage-hash",
                    run_id,
                    "coverage decision hash failed self-validation",
                    "report this coverage integrity defect; no unvalidated coverage is emitted",
                );
            }
            Some(report)
        }
        None => None,
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
            let _ = write!(
                out,
                ",\"claim_count\":{},\"checker\":\"pass\",\"authority\":\"structural-integrity-only\"",
                export.claim_count
            );
            if let Some(cov) = &report {
                push_coverage_json(&mut out, cov);
            }
            out.push_str("}\n");
            out
        }
        OutputMode::Text => {
            let mut out = format!(
                "status=ok\ncommand=package\nsubject={run_id}\nrun={run_id}\npackage={}\nmerkle_root={}\nclaim_count={}\nchecker=pass\nauthority=structural-integrity-only\n",
                escape_text(&export.path.to_string_lossy()),
                export.merkle_root,
                export.claim_count,
            );
            if let Some(cov) = &report {
                push_coverage_text(&mut out, cov);
            }
            out
        }
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}
