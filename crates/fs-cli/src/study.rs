//! Canonical `frankensim study` routing.
//!
//! Thermal and free-boundary elasticity studies share one receipt envelope and
//! user-facing command while retaining distinct numerical producers and claims.

use std::path::Path;

use crate::{CommandOutput, OutputMode};

#[path = "study_elasticity.rs"]
mod elasticity;
#[path = "study_thermal.rs"]
mod thermal;

/// Shared retained-study receipt envelope. Numerical producer identity remains
/// an explicit `driver` field inside every receipt.
pub const STUDY_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.study-run-receipt.v1";

pub(crate) fn study_path(
    path: &Path,
    ledger_path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    if elasticity::looks_like(path) {
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
    if elasticity::owns_run(pointer, path) {
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
    if path.is_some_and(|ledger| elasticity::owns_run(pointer, ledger)) {
        elasticity::export(command, pointer, path, mode)
    } else {
        thermal::export(command, pointer, path, mode)
    }
}
