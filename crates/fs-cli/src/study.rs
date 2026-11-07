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
pub(crate) use uncertainty::{PreparedStudy, StudyPins};

enum StudyFileModel {
    Thermal(thermal::PreparedStudy),
    Elasticity(elasticity::PreparedStudy),
}

/// The actual self-contained study model admitted by `study.run`.
/// The source path is consumed only while loading; execution borrows this state.
pub(crate) struct PreparedStudyFile {
    source_hash: fs_blake3::ContentHash,
    model: StudyFileModel,
}

impl std::fmt::Debug for PreparedStudyFile {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.debug_struct("PreparedStudyFile")
            .field("source_hash", &self.source_hash.to_hex())
            .field("seed", &self.seed())
            .finish_non_exhaustive()
    }
}

impl PreparedStudyFile {
    pub(crate) fn load(path: &Path, input_limit: u64) -> Result<Self, String> {
        let bytes = uncertainty::read_input_snapshot(path, input_limit.min(crate::MAX_PROJECT_BYTES))?;
        let source_hash = fs_blake3::hash_bytes(&bytes);
        let source = std::str::from_utf8(&bytes).map_err(|error| error.to_string())?;
        let model = if elasticity::looks_like_source(source) {
            StudyFileModel::Elasticity(elasticity::PreparedStudy::parse(source)?)
        } else {
            StudyFileModel::Thermal(thermal::PreparedStudy::parse(source)?)
        };
        let prepared = Self { source_hash, model };
        if prepared.input_bytes() > input_limit {
            return Err("canonical study exceeds the remaining input-storage allowance".into());
        }
        Ok(prepared)
    }

    pub(crate) fn source_hash(&self) -> fs_blake3::ContentHash { self.source_hash }

    /// Retained canonical input bytes, not a peak-RSS or decoded-state bound.
    pub(crate) fn input_bytes(&self) -> u64 {
        match &self.model {
            StudyFileModel::Thermal(study) => study.source().len() as u64,
            StudyFileModel::Elasticity(study) => study.source().len() as u64,
        }
    }

    pub(crate) fn seed(&self) -> u64 {
        match &self.model {
            StudyFileModel::Thermal(study) => study.seed(),
            StudyFileModel::Elasticity(study) => study.seed(),
        }
    }

    pub(crate) fn check_budget(&self, text: Option<&str>) -> Result<(), String> {
        match &self.model {
            StudyFileModel::Thermal(study) => study.check_budget(text),
            StudyFileModel::Elasticity(study) => study.check_budget(text),
        }
    }

    pub(crate) fn run(
        &self, ledger_path: &Path, override_text: Option<&str>, mode: OutputMode,
    ) -> CommandOutput {
        match &self.model {
            StudyFileModel::Thermal(study) => study.run(ledger_path, override_text, mode),
            StudyFileModel::Elasticity(study) => study.run(ledger_path, override_text, mode),
        }
    }
}

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
