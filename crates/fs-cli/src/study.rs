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

/// Logical begin/end timestamps of one study persistence op (carries no other metrics).
pub(crate) const STUDY_OP_START_T: i64 = 0;
pub(crate) const STUDY_OP_END_T: i64 = 1;
/// Reserved metric key recording cumulative study wall-seconds outside content-addressed
/// `ops.budget` and `study-run-receipt` bytes (bead `frankensim-4fs7g`).
pub(crate) const STUDY_CONSUMED_WALL_S_METRIC: &str = "study.consumed_wall_s";

/// Atomically persist cumulative study wall-seconds on `op` at `STUDY_OP_END_T`.
#[allow(clippy::result_large_err)]
pub(crate) fn record_study_wall_s(
    ledger: &fs_ledger::Ledger,
    op: i64,
    wall_s: f64,
) -> Result<(), fs_ledger::LedgerError> {
    if !wall_s.is_finite() || wall_s < 0.0 {
        return Err(fs_ledger::LedgerError::Invalid {
            field: STUDY_CONSUMED_WALL_S_METRIC.to_string(),
            problem: format!("study wall charge {wall_s} must be finite and nonnegative"),
        });
    }
    ledger.record_metric(op, STUDY_OP_END_T, STUDY_CONSUMED_WALL_S_METRIC, wall_s)
}

/// Resolve and validate cumulative wall-seconds for an attested study receipt from its
/// producing op's atomic `metrics` row (`study.consumed_wall_s` at `STUDY_OP_END_T`) or
/// legacy pre-`4fs7g` JSON field. Fail-closed on unsealed, missing, duplicate, mistimed,
/// conflicting, negative, or non-finite values.
#[allow(clippy::float_cmp)]
pub(crate) fn resolve_study_wall_s(
    ledger: &fs_ledger::Ledger,
    receipt_hash: &fs_blake3::ContentHash,
    receipt: &crate::json_read::JsonValue,
) -> Result<f64, String> {
    let op_id = ledger
        .artifact_output_seal(receipt_hash)
        .map_err(|e| e.to_string())?
        .ok_or_else(|| "unsealed study receipt lacks producing operation".to_string())?;
    let mut metric_wall: Option<f64> = None;
    for (t, name, value) in ledger.metrics_for_op(op_id).map_err(|e| e.to_string())? {
        if name == STUDY_CONSUMED_WALL_S_METRIC
            && (t != STUDY_OP_END_T || metric_wall.replace(value).is_some())
        {
            return Err(format!(
                "study op {op_id} carries duplicate or mistimed `{STUDY_CONSUMED_WALL_S_METRIC}` metric"
            ));
        }
    }
    let legacy_wall = match receipt.get("consumed_wall_s") {
        Some(v) => Some(v.as_f64().ok_or_else(|| {
            format!("study op {op_id} receipt has non-numeric `consumed_wall_s` field")
        })?),
        None => None,
    };
    let wall_s = match (metric_wall, legacy_wall) {
        (Some(m), Some(l)) if m != l => {
            return Err(format!(
                "study op {op_id} wall metric conflicts with legacy receipt field"
            ));
        }
        (Some(m), _) => m,
        (None, Some(l)) => l,
        (None, None) => {
            return Err(format!(
                "study op {op_id} lacks retained `{STUDY_CONSUMED_WALL_S_METRIC}` metric"
            ));
        }
    };
    if !wall_s.is_finite() || wall_s < 0.0 {
        return Err(format!("study op {op_id} carries invalid wall charge"));
    }
    Ok(wall_s)
}

/// Return retained `receipt_kind` output artifact hashes for completed ops of `study_id`
/// on `MAIN_BRANCH`, newest first, for same-ledger idempotent rerun re-attestation.
#[allow(clippy::result_large_err)]
pub(crate) fn retained_study_receipts(
    ledger: &fs_ledger::Ledger,
    study_id: fs_blake3::ContentHash,
    receipt_kind: &str,
) -> Result<Vec<fs_blake3::ContentHash>, fs_ledger::LedgerError> {
    let mut out = Vec::new();
    for id in ledger
        .visible_op_ids(fs_ledger::MAIN_BRANCH, None)?
        .into_iter()
        .rev()
    {
        let Some(row) = ledger.op(id)? else { continue };
        if row.session.as_deref() != Some(study_id.as_bytes().as_slice())
            || row.outcome.as_deref() != Some("ok")
        {
            continue;
        }
        for edge in ledger.op_artifact_edges_bounded(id, 64)?.edges {
            if edge.role == fs_ledger::EdgeRole::Out
                && ledger
                    .artifact_info(&edge.artifact)?
                    .is_some_and(|info| info.kind == receipt_kind)
                && !out.contains(&edge.artifact)
            {
                out.push(edge.artifact);
            }
        }
    }
    Ok(out)
}

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
