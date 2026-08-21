//! E10.0 registry-conformance AUDIT (bead wf-root-guzez.11.1; the
//! freeze itself is E1.7). Confirms that no partition, metric, prior,
//! or band changed after protected results were observed: every
//! artifact the E1.7 registry froze by content hash is re-hashed from
//! the exact bytes on disk, and the result is a RECEIPT — conformant
//! rows, or violations enumerated WITH their identity history (the
//! frozen hash vs the observed hash IS the history at this tier).
//!
//! The audit distrusts the registry file too: callers who hold the
//! EvidenceRegistryId (the registry's own content hash, cited by
//! consumers) pass it in, and a doctored registry refuses before any
//! row is judged.

use crate::{Refusal, refuse};
use fs_blake3::{hash_bytes, hash_domain};
use std::fs;
use std::path::Path;

/// Frozen-file-count cap (a registry with more is a different schema).
pub const MAX_FROZEN_FILES: usize = 64;

/// One audited artifact row.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AuditRow {
    /// Artifact file name.
    pub file: String,
    /// Hash the freeze pinned.
    pub frozen_hex: String,
    /// Hash of today's bytes.
    pub current_hex: String,
    /// Conformant?
    pub conformant: bool,
}

/// The committed audit receipt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegistryAuditReceiptV1 {
    /// Schema id.
    pub schema: &'static str,
    /// The registry's own content hash (EvidenceRegistryId).
    pub evidence_registry_id: String,
    /// Registry status string at audit time.
    pub registry_status: String,
    /// Per-artifact rows (every frozen file, never a summary only).
    pub rows: Vec<AuditRow>,
    /// Violation count.
    pub violations: usize,
    /// Verdict.
    pub verdict: &'static str,
    /// Digest over the receipt payload.
    pub receipt_digest: String,
}

impl RegistryAuditReceiptV1 {
    /// Bounded JSONL rendering (the committed audit trail).
    #[must_use]
    pub fn to_json(&self) -> String {
        let rows: Vec<String> = self
            .rows
            .iter()
            .map(|r| {
                format!(
                    "  {{\"file\": \"{}\", \"frozen\": \"{}\", \"current\": \"{}\", \"conformant\": {}}}",
                    r.file, r.frozen_hex, r.current_hex, r.conformant
                )
            })
            .collect();
        format!(
            "{{\n \"schema\": \"{}\",\n \"evidence_registry_id\": \"{}\",\n \"registry_status\": \"{}\",\n \"rows\": [\n{}\n ],\n \"violations\": {},\n \"verdict\": \"{}\",\n \"receipt_digest\": \"{}\"\n}}\n",
            self.schema,
            self.evidence_registry_id,
            self.registry_status,
            rows.join(",\n"),
            self.violations,
            self.verdict,
            self.receipt_digest
        )
    }
}

/// Dependency-free extraction of the `files` block and status (the
/// same honest tier as the E1.7 guard: parse only what is asserted).
fn parse_registry(text: &str) -> Result<(Vec<(String, String)>, String), Refusal> {
    let files_block = text
        .split("\"files\": {")
        .nth(1)
        .and_then(|rest| rest.split('}').next())
        .ok_or_else(|| {
            refuse(
                "registry-shape-invalid",
                "no files block".into(),
                "the E1.7 registry carries frozen_artifacts.files",
            )
        })?;
    let mut files = Vec::new();
    for line in files_block.lines() {
        let line = line.trim().trim_end_matches(',');
        if let Some((k, v)) = line.split_once("\": \"") {
            files.push((
                k.trim_start_matches('"').to_string(),
                v.trim_end_matches('"').to_string(),
            ));
        }
    }
    let status = text
        .split("\"status\": \"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .unwrap_or("MISSING")
        .to_string();
    Ok((files, status))
}

/// Run the conformance audit against a registry directory.
///
/// # Errors
/// `registry-unreadable`; `registry-shape-invalid`; `registry-files-
/// invalid` (0 files, or beyond the cap — AT the cap admits);
/// `registry-identity-mismatch` (the registry file itself does not
/// hash to the id the caller cites — a doctored registry never gets
/// to judge rows).
pub fn audit_registry(
    data_dir: &Path,
    registry_file: &str,
    expected_registry_id: Option<&str>,
) -> Result<RegistryAuditReceiptV1, Refusal> {
    let reg_path = data_dir.join(registry_file);
    let reg_bytes = fs::read(&reg_path).map_err(|e| {
        refuse(
            "registry-unreadable",
            format!("{}: {e}", reg_path.display()),
            "audit runs against the tracked data directory",
        )
    })?;
    let evidence_registry_id = hash_bytes(&reg_bytes).to_hex().to_string();
    if let Some(expected) = expected_registry_id {
        if expected != evidence_registry_id {
            return Err(refuse(
                "registry-identity-mismatch",
                format!("registry hashes to {evidence_registry_id}, caller cites {expected}"),
                "a doctored registry never judges rows; recover the cited registry",
            ));
        }
    }
    let text = String::from_utf8_lossy(&reg_bytes);
    let (files, registry_status) = parse_registry(&text)?;
    if files.is_empty() || files.len() > MAX_FROZEN_FILES {
        return Err(refuse(
            "registry-files-invalid",
            format!(
                "{} frozen files outside [1, {MAX_FROZEN_FILES}]",
                files.len()
            ),
            "the v1 freeze covers eight artifacts",
        ));
    }
    let mut rows = Vec::with_capacity(files.len());
    let mut violations = 0usize;
    for (file, frozen_hex) in files {
        let current_hex = match fs::read(data_dir.join(&file)) {
            Ok(b) => hash_bytes(&b).to_hex().to_string(),
            Err(e) => format!("UNREADABLE:{e}"),
        };
        let conformant = current_hex == frozen_hex;
        if !conformant {
            violations += 1;
        }
        rows.push(AuditRow {
            file,
            frozen_hex,
            current_hex,
            conformant,
        });
    }
    let verdict = if violations == 0 {
        "CONFORMANT"
    } else {
        "VIOLATED"
    };
    let mut b = evidence_registry_id.as_bytes().to_vec();
    for r in &rows {
        b.extend_from_slice(r.file.as_bytes());
        b.extend_from_slice(r.frozen_hex.as_bytes());
        b.extend_from_slice(r.current_hex.as_bytes());
        b.push(u8::from(r.conformant));
    }
    let receipt_digest = hash_domain("org.frankensim.wf.registry-audit.v1", &b).to_hex();
    Ok(RegistryAuditReceiptV1 {
        schema: "org.frankensim.wf.registry-audit.v1",
        evidence_registry_id,
        registry_status,
        rows,
        violations,
        verdict,
        receipt_digest,
    })
}
