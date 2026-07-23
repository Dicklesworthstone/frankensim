//! Deterministic end-to-end wedge decision audit (bead
//! `frankensim-extreal-program-f85xj.1.5`).
//!
//! This module composes the measured inputs, explicit comparison,
//! sensitivity tables, cycle-time baseline, kill-criterion evaluation, and
//! fail-closed ratification record into one content-addressed artifact.  The
//! caller must supply a positive measured FrankenSim cycle time and a
//! non-empty evidence locator for that measurement.  The locator is retained
//! as provenance data; this pure module does not authenticate or dereference
//! it.

use crate::{RatificationError, json_escape, ratification_json, ratified_vertical};
use core::fmt;
use core::fmt::Write as _;
use fs_blake3::{ContentHash, hash_domain};
use fs_wedge::{
    CHT_BASELINE, KillCriterionError, KillCriterionEvaluation, KillVerdict, ScoringError, audit,
    default_recommendation, render_comparison_report, to_json as wedge_manifest_json,
};

/// Stable schema label for the complete wedge decision audit artifact.
pub const WEDGE_DECISION_AUDIT_SCHEMA: &str = "frankensim-wedge-decision-audit-v1";
/// Domain separating the artifact payload identity from every other content
/// identity in the constellation.
pub const WEDGE_DECISION_AUDIT_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.wedge-decision-audit.v1";
/// Maximum UTF-8 bytes admitted for the caller's cycle-time evidence locator.
pub const MAX_WEDGE_AUDIT_EVIDENCE_BYTES: usize = 4_096;
/// Maximum bytes published in one audit artifact.
pub const MAX_WEDGE_AUDIT_ARTIFACT_BYTES: usize = 1_048_576;

/// Structured severity for deterministic audit progress and refusal logs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WedgeAuditLogLevel {
    /// A successful derivation step.
    Info,
    /// A refusal that requires caller action.
    Warn,
}

impl WedgeAuditLogLevel {
    const fn label(self) -> &'static str {
        match self {
            Self::Info => "info",
            Self::Warn => "warn",
        }
    }
}

/// One deterministic structured log row from the wedge audit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WedgeAuditLog {
    level: WedgeAuditLogLevel,
    step: &'static str,
    code: Option<&'static str>,
    detail: String,
    fix: Option<&'static str>,
}

impl WedgeAuditLog {
    fn info(step: &'static str, detail: String) -> Self {
        Self {
            level: WedgeAuditLogLevel::Info,
            step,
            code: None,
            detail,
            fix: None,
        }
    }

    /// Build an actionable warning row for a refused CLI or audit request.
    #[must_use]
    pub fn warning(
        step: &'static str,
        code: &'static str,
        detail: impl Into<String>,
        fix: &'static str,
    ) -> Self {
        Self {
            level: WedgeAuditLogLevel::Warn,
            step,
            code: Some(code),
            detail: detail.into(),
            fix: Some(fix),
        }
    }

    /// Stable JSON-lines rendering.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"level\":\"{}\",\"step\":\"{}\",\"code\":{},\"detail\":\"{}\",\"fix\":{}}}",
            self.level.label(),
            json_escape(self.step),
            self.code.map_or_else(
                || "null".to_string(),
                |code| format!("\"{}\"", json_escape(code))
            ),
            json_escape(&self.detail),
            self.fix.map_or_else(
                || "null".to_string(),
                |fix| format!("\"{}\"", json_escape(fix))
            ),
        )
    }
}

/// Validated caller input for one wedge decision audit.
#[derive(Debug, Clone, PartialEq)]
pub struct WedgeDecisionAuditRequest {
    measured_days: f64,
    cycle_time_evidence: String,
}

impl WedgeDecisionAuditRequest {
    /// Validate a measured cycle time and its evidence locator.
    pub fn new(
        measured_days: f64,
        cycle_time_evidence: impl Into<String>,
    ) -> Result<Self, WedgeDecisionAuditError> {
        if !measured_days.is_finite() || measured_days <= 0.0 {
            return Err(WedgeDecisionAuditError::KillCriterion(
                KillCriterionError::NonMeasurableCycleTime { measured_days },
            ));
        }
        let cycle_time_evidence = cycle_time_evidence.into();
        if cycle_time_evidence.trim().is_empty() {
            return Err(WedgeDecisionAuditError::MissingCycleTimeEvidence);
        }
        if cycle_time_evidence.len() > MAX_WEDGE_AUDIT_EVIDENCE_BYTES {
            return Err(WedgeDecisionAuditError::CycleTimeEvidenceTooLong {
                observed: cycle_time_evidence.len(),
            });
        }
        Ok(Self {
            measured_days,
            cycle_time_evidence,
        })
    }

    /// Measured FrankenSim cycle time in working days.
    #[must_use]
    pub const fn measured_days(&self) -> f64 {
        self.measured_days
    }

    /// Caller-supplied evidence locator retained by the audit.
    #[must_use]
    pub fn cycle_time_evidence(&self) -> &str {
        &self.cycle_time_evidence
    }
}

/// Refusal from the complete wedge decision audit boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum WedgeDecisionAuditError {
    /// The measured cycle time carried no provenance locator.
    MissingCycleTimeEvidence,
    /// The cycle-time provenance locator exceeded its resource cap.
    CycleTimeEvidenceTooLong {
        /// Observed UTF-8 byte length.
        observed: usize,
    },
    /// The underlying fs-wedge self-audit found drift or incompleteness.
    WedgeSelfAudit {
        /// Deterministically ordered human-readable gaps.
        gaps: Vec<String>,
    },
    /// The measured comparison or sensitivity render refused.
    Scoring(ScoringError),
    /// The program ratification record refused current inputs.
    Ratification(RatificationError),
    /// The cycle-time kill criterion refused its input or baseline.
    KillCriterion(KillCriterionError),
    /// The complete artifact exceeded its hard byte cap.
    ArtifactTooLarge {
        /// Observed artifact byte length.
        observed: usize,
    },
}

impl WedgeDecisionAuditError {
    /// Stable machine-readable refusal code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::MissingCycleTimeEvidence => "missing-cycle-time-evidence",
            Self::CycleTimeEvidenceTooLong { .. } => "cycle-time-evidence-too-long",
            Self::WedgeSelfAudit { .. } => "wedge-self-audit-failed",
            Self::Scoring(_) => "measured-scoring-refused",
            Self::Ratification(_) => "ratification-refused",
            Self::KillCriterion(_) => "kill-criterion-refused",
            Self::ArtifactTooLarge { .. } => "wedge-audit-artifact-too-large",
        }
    }

    /// Stable actionable remedy suitable for warning logs.
    #[must_use]
    pub const fn fix(&self) -> &'static str {
        match self {
            Self::MissingCycleTimeEvidence => {
                "provide --cycle-time-evidence with the retained run or timing receipt locator"
            }
            Self::CycleTimeEvidenceTooLong { .. } => {
                "use a bounded content-addressed locator rather than embedding the evidence"
            }
            Self::WedgeSelfAudit { .. } => {
                "repair the named fs-wedge drift gaps before rendering a decision audit"
            }
            Self::Scoring(_) => {
                "repair the measured factor table, normalized weights, or sensitivity inputs"
            }
            Self::Ratification(_) => {
                "repair or supersede the ratification record against the current measured inputs"
            }
            Self::KillCriterion(_) => {
                "provide a positive finite measured cycle time and a complete non-placeholder baseline"
            }
            Self::ArtifactTooLarge { .. } => {
                "reduce diagnostic payload size or version the artifact resource envelope"
            }
        }
    }
}

impl fmt::Display for WedgeDecisionAuditError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingCycleTimeEvidence => {
                formatter.write_str("measured cycle time has no evidence locator")
            }
            Self::CycleTimeEvidenceTooLong { observed } => write!(
                formatter,
                "cycle-time evidence locator is {observed} bytes; limit is {MAX_WEDGE_AUDIT_EVIDENCE_BYTES}"
            ),
            Self::WedgeSelfAudit { gaps } => {
                write!(formatter, "fs-wedge self-audit failed: {}", gaps.join("; "))
            }
            Self::Scoring(source) => write!(formatter, "measured scoring refused: {source:?}"),
            Self::Ratification(source) => {
                write!(formatter, "vertical ratification refused: {source:?}")
            }
            Self::KillCriterion(source) => {
                write!(formatter, "cycle-time kill criterion refused: {source:?}")
            }
            Self::ArtifactTooLarge { observed } => write!(
                formatter,
                "wedge audit artifact is {observed} bytes; limit is {MAX_WEDGE_AUDIT_ARTIFACT_BYTES}"
            ),
        }
    }
}

impl std::error::Error for WedgeDecisionAuditError {}

/// Complete deterministic artifact plus the logs that reconstruct its
/// derivation.
#[derive(Debug, Clone, PartialEq)]
pub struct WedgeDecisionAudit {
    identity: ContentHash,
    artifact: String,
    logs: Vec<WedgeAuditLog>,
    kill_verdict: KillVerdict,
}

impl WedgeDecisionAudit {
    /// Domain-separated identity of the exact artifact payload.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Complete deterministic JSON artifact.
    #[must_use]
    pub fn artifact(&self) -> &str {
        &self.artifact
    }

    /// Ordered structured logs. The constituent-data rows plus the final
    /// artifact row permit reconstruction without rerunning the command.
    #[must_use]
    pub fn logs(&self) -> &[WedgeAuditLog] {
        &self.logs
    }

    /// Conservative kill-criterion verdict carried by the artifact.
    #[must_use]
    pub const fn kill_verdict(&self) -> KillVerdict {
        self.kill_verdict
    }
}

fn audit_checks_json(report: &fs_wedge::WedgeAudit) -> String {
    let mut out = String::from("[");
    for (index, check) in report.checks.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"name\":\"{}\",\"passed\":{}}}",
            json_escape(check.name),
            check.passed
        )
        .expect("write to String");
    }
    out.push(']');
    out
}

fn kill_evaluation_json(evaluation: &KillCriterionEvaluation) -> String {
    format!(
        "{{\"vertical\":\"{}\",\"assembled_on\":\"{}\",\"baseline_provenance\":\"{}\",\"baseline_days_low\":{},\"baseline_days_high\":{},\"target_reduction\":{},\"measured_days\":{},\"measured_days_bits\":\"{:016x}\",\"reduction_low\":{},\"reduction_high\":{},\"verdict\":\"{}\",\"derivation\":\"{}\"}}",
        json_escape(evaluation.vertical),
        json_escape(evaluation.assembled_on),
        evaluation.provenance.label(),
        evaluation.baseline_days_low,
        evaluation.baseline_days_high,
        evaluation.target_reduction,
        evaluation.measured_days,
        evaluation.measured_days.to_bits(),
        evaluation.reduction_low,
        evaluation.reduction_high,
        evaluation.verdict.label(),
        json_escape(&evaluation.derivation),
    )
}

/// Render the complete measured wedge decision state and every derivation log
/// row without reading or writing external state.
pub fn build_wedge_decision_audit(
    request: &WedgeDecisionAuditRequest,
) -> Result<WedgeDecisionAudit, WedgeDecisionAuditError> {
    let wedge_audit = audit();
    if !wedge_audit.ok() {
        return Err(WedgeDecisionAuditError::WedgeSelfAudit {
            gaps: wedge_audit.gaps,
        });
    }
    let checks_json = audit_checks_json(&wedge_audit);
    let recommendation = default_recommendation().map_err(WedgeDecisionAuditError::Scoring)?;
    let comparison_report = render_comparison_report().map_err(WedgeDecisionAuditError::Scoring)?;
    let ratification = ratified_vertical().map_err(WedgeDecisionAuditError::Ratification)?;
    let ratification_json = ratification_json().map_err(WedgeDecisionAuditError::Ratification)?;
    let evaluation = CHT_BASELINE
        .evaluate_kill_criterion(request.measured_days)
        .map_err(WedgeDecisionAuditError::KillCriterion)?;
    let kill_json = kill_evaluation_json(&evaluation);
    let wedge_manifest = wedge_manifest_json();

    let payload = format!(
        "{{\"authority\":\"commercial-decision-audit-only\",\"cycle_time_input\":{{\"measured_days\":{},\"measured_days_bits\":\"{:016x}\",\"evidence\":\"{}\"}},\"wedge_self_audit\":{},\"wedge_manifest\":{},\"comparison_report_tsv\":\"{}\",\"kill_evaluation\":{},\"ratification\":{},\"no_claim\":\"The caller-supplied cycle-time locator is retained but not authenticated or dereferenced; this artifact does not prove market demand, technical maturity, scientific validity, or that the kill criterion was met by a production run.\"}}",
        request.measured_days,
        request.measured_days.to_bits(),
        json_escape(&request.cycle_time_evidence),
        checks_json,
        wedge_manifest,
        json_escape(&comparison_report),
        kill_json,
        ratification_json,
    );
    let identity = hash_domain(WEDGE_DECISION_AUDIT_IDENTITY_DOMAIN, payload.as_bytes());
    let artifact = format!(
        "{{\"schema\":\"{WEDGE_DECISION_AUDIT_SCHEMA}\",\"identity\":\"{identity}\",\"payload\":{payload}}}"
    );
    if artifact.len() > MAX_WEDGE_AUDIT_ARTIFACT_BYTES {
        return Err(WedgeDecisionAuditError::ArtifactTooLarge {
            observed: artifact.len(),
        });
    }

    let input_json = format!(
        "{{\"measured_days\":{},\"measured_days_bits\":\"{:016x}\",\"cycle_time_evidence\":\"{}\"}}",
        request.measured_days,
        request.measured_days.to_bits(),
        json_escape(&request.cycle_time_evidence),
    );
    let ranking_json = format!(
        "{{\"recommended\":\"{}\",\"runner_up\":\"{}\",\"report_tsv\":\"{}\"}}",
        recommendation.recommended,
        recommendation.runner_up,
        json_escape(&comparison_report),
    );
    let logs = vec![
        WedgeAuditLog::info("request-validated", input_json),
        WedgeAuditLog::info("wedge-self-audit", checks_json),
        WedgeAuditLog::info("measured-inputs-and-baseline", wedge_manifest),
        WedgeAuditLog::info("scoring-and-sensitivity", ranking_json),
        WedgeAuditLog::info(
            "ratification-validated",
            format!(
                "{{\"record_id\":\"{}\",\"record\":{}}}",
                json_escape(ratification.id),
                ratification_json
            ),
        ),
        WedgeAuditLog::info("kill-criterion-evaluated", kill_json),
        WedgeAuditLog::info(
            "artifact-rendered",
            format!(
                "{{\"identity\":\"{identity}\",\"bytes\":{},\"artifact\":{artifact}}}",
                artifact.len()
            ),
        ),
    ];

    Ok(WedgeDecisionAudit {
        identity,
        artifact,
        logs,
        kill_verdict: evaluation.verdict,
    })
}
