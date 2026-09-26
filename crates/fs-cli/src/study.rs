//! Canonical `frankensim study` routing.
//!
//! Thermal, free-boundary elasticity and native uncertainty studies share one
//! command while retaining distinct numerical producers and authority claims.

use std::path::Path;

use crate::{CommandOutput, OutputMode};

#[path = "study_elasticity.rs"]
mod elasticity;
#[path = "study_thermal.rs"]
mod thermal;
#[path = "study_uncertainty.rs"]
mod uncertainty;

/// Shared retained-study receipt envelope. Numerical producer identity remains
/// an explicit `driver` field inside every receipt.
pub const STUDY_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.study-run-receipt.v1";

pub(crate) fn study_path(
    path: &Path,
    ledger_path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    if uncertainty::looks_like(path) {
        uncertainty::study_path(path, ledger_path, override_text, mode)
    } else if elasticity::looks_like(path) {
        elasticity::study_path(path, ledger_path, override_text, mode)
    } else {
        thermal::study_path(path, ledger_path, override_text, mode)
    }
}

pub(crate) fn resume_path(
    pointer: &str,
    path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    if uncertainty::owns_run(pointer, path) {
        uncertainty::resume_path(pointer, path, override_text, mode)
    } else if elasticity::owns_run(pointer, path) {
        elasticity::resume_path(pointer, path, override_text, mode)
    } else {
        thermal::resume_path(pointer, path, override_text, mode)
    }
}

pub(crate) fn export(
    command: &'static str,
    pointer: &str,
    path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    if path.is_some_and(|ledger| uncertainty::owns_run(pointer, ledger)) {
        uncertainty::export(command, pointer, path, mode)
    } else if path.is_some_and(|ledger| elasticity::owns_run(pointer, ledger)) {
        elasticity::export(command, pointer, path, mode)
    } else {
        thermal::export(command, pointer, path, mode)
    }
}
