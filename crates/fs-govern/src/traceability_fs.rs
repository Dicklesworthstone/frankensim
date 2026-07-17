//! Bounded filesystem adapter for traceability source snapshots.
//!
//! The declaration registry remains pure data. This module is the explicit
//! world-boundary adapter that reads root-relative Beads, contract, and
//! registry files, hashes the exact bytes actually read, and feeds those
//! identities into the existing sealed source snapshot. Beads JSONL receives
//! a conservative lexical-envelope and unique top-level `id` audit only; this
//! module does not interpret tracker status, dependencies, or scientific proof.

use crate::traceability::{
    MAX_TRACEABILITY_FIELD_BYTES, MAX_TRACEABILITY_SOURCE_LOCATOR_BYTES, MAX_TRACEABILITY_SOURCES,
    TraceabilitySource, TraceabilitySourceKind, TraceabilitySourceSnapshot,
};
use fs_blake3::{ContentHash, hash_bytes};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::Read as _;
use std::path::{Component, Path, PathBuf};

/// Versioned meaning of the bounded filesystem source adapter.
pub const TRACEABILITY_FILESYSTEM_ADAPTER_VERSION: &str =
    "frankensim-traceability-filesystem-adapter-v1";

/// Hard maximum bytes read from one source artifact.
pub const MAX_TRACEABILITY_SOURCE_FILE_BYTES: usize = 64 * 1024 * 1024;
/// Hard maximum bytes read across one source snapshot.
pub const MAX_TRACEABILITY_TOTAL_SOURCE_BYTES: usize = 256 * 1024 * 1024;
/// Hard maximum nonblank Beads JSONL records.
pub const MAX_TRACEABILITY_BEADS_RECORDS: usize = 100_000;
/// Hard maximum bytes in one Beads JSONL record.
pub const MAX_TRACEABILITY_BEADS_LINE_BYTES: usize = 1024 * 1024;
/// Hard maximum JSON container nesting in one Beads record.
pub const MAX_TRACEABILITY_BEADS_JSON_NESTING: usize = 256;

/// Defensive limits for one filesystem source-snapshot load.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TraceabilityFilesystemLimits {
    /// Maximum source specs.
    pub max_sources: usize,
    /// Maximum bytes read from one regular file.
    pub max_source_bytes: usize,
    /// Maximum bytes read across all files.
    pub max_total_bytes: usize,
    /// Maximum nonblank Beads JSONL records.
    pub max_beads_records: usize,
    /// Maximum bytes in one Beads JSONL record.
    pub max_beads_line_bytes: usize,
    /// Maximum JSON container nesting in one Beads record.
    pub max_beads_json_nesting: usize,
}

impl Default for TraceabilityFilesystemLimits {
    fn default() -> Self {
        Self {
            max_sources: MAX_TRACEABILITY_SOURCES,
            max_source_bytes: MAX_TRACEABILITY_SOURCE_FILE_BYTES,
            max_total_bytes: MAX_TRACEABILITY_TOTAL_SOURCE_BYTES,
            max_beads_records: MAX_TRACEABILITY_BEADS_RECORDS,
            max_beads_line_bytes: MAX_TRACEABILITY_BEADS_LINE_BYTES,
            max_beads_json_nesting: MAX_TRACEABILITY_BEADS_JSON_NESTING,
        }
    }
}

/// One root-relative source file to load and bind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TraceabilityFileSpec<'a> {
    /// Artifact class.
    pub kind: TraceabilitySourceKind,
    /// Portable root-relative path; absolute, parent, and current-directory
    /// components refuse.
    pub relative_path: &'a Path,
}

/// Field named by one filesystem-adapter diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum TraceabilityFilesystemField {
    /// Adapter limit configuration.
    Limits,
    /// Snapshot root directory.
    Root,
    /// Portable root-relative source path.
    RelativePath,
    /// Regular-file metadata or root-containment check.
    Metadata,
    /// Bounded exact file bytes.
    Content,
    /// Beads JSONL envelope and canonical ID index.
    BeadsJsonl,
    /// Existing Beads/contract/registry source-snapshot admission.
    Snapshot,
}

impl TraceabilityFilesystemField {
    /// Stable diagnostic field name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Limits => "limits",
            Self::Root => "root",
            Self::RelativePath => "relative_path",
            Self::Metadata => "metadata",
            Self::Content => "content",
            Self::BeadsJsonl => "beads_jsonl",
            Self::Snapshot => "snapshot",
        }
    }
}

/// One deterministic filesystem source refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityFilesystemDiagnostic {
    /// Portable source locator or aggregate scope.
    pub source: String,
    /// Exact failed adapter field.
    pub field: TraceabilityFilesystemField,
    /// Actionable refusal reason.
    pub reason: String,
}

/// Complete fail-closed filesystem adapter audit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityFilesystemAudit {
    /// Number of offered source specs.
    pub total: usize,
    /// Every deterministic refusal.
    pub diagnostics: Vec<TraceabilityFilesystemDiagnostic>,
}

impl TraceabilityFilesystemAudit {
    /// Whether all files and the resulting source snapshot were admitted.
    #[must_use]
    pub fn ok(&self) -> bool {
        self.total > 0 && self.diagnostics.is_empty()
    }
}

/// Exact read receipt for one filesystem source artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityFilesystemSourceReceipt {
    kind: TraceabilitySourceKind,
    locator: String,
    content_identity: ContentHash,
    byte_count: usize,
    beads_record_count: Option<usize>,
}

impl TraceabilityFilesystemSourceReceipt {
    /// Artifact class.
    #[must_use]
    pub const fn kind(&self) -> TraceabilitySourceKind {
        self.kind
    }

    /// Canonical portable root-relative locator.
    #[must_use]
    pub fn locator(&self) -> &str {
        &self.locator
    }

    /// BLAKE3 identity of the exact bytes read.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.content_identity
    }

    /// Exact byte count read and hashed.
    #[must_use]
    pub const fn byte_count(&self) -> usize {
        self.byte_count
    }

    /// Number of admitted nonblank Beads records, only for Beads sources.
    #[must_use]
    pub const fn beads_record_count(&self) -> Option<usize> {
        self.beads_record_count
    }
}

/// Complete filesystem read receipt bound to the sealed source snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceabilityFilesystemReceipt {
    adapter_version: &'static str,
    limits: TraceabilityFilesystemLimits,
    total_bytes: usize,
    sources: Vec<TraceabilityFilesystemSourceReceipt>,
    source_snapshot_identity: ContentHash,
}

impl TraceabilityFilesystemReceipt {
    /// Filesystem adapter semantics version.
    #[must_use]
    pub const fn adapter_version(&self) -> &'static str {
        self.adapter_version
    }

    /// Exact defensive limits.
    #[must_use]
    pub const fn limits(&self) -> TraceabilityFilesystemLimits {
        self.limits
    }

    /// Total bytes read and hashed.
    #[must_use]
    pub const fn total_bytes(&self) -> usize {
        self.total_bytes
    }

    /// Canonically ordered per-source receipts.
    #[must_use]
    pub fn sources(&self) -> &[TraceabilityFilesystemSourceReceipt] {
        &self.sources
    }

    /// Identity of the existing sealed source snapshot.
    #[must_use]
    pub const fn source_snapshot_identity(&self) -> ContentHash {
        self.source_snapshot_identity
    }
}

/// Sealed source snapshot paired with the concrete filesystem read receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoadedTraceabilitySourceSnapshot {
    snapshot: TraceabilitySourceSnapshot,
    receipt: TraceabilityFilesystemReceipt,
}

impl LoadedTraceabilitySourceSnapshot {
    /// Existing pure source snapshot for declaration-ledger binding.
    #[must_use]
    pub const fn snapshot(&self) -> &TraceabilitySourceSnapshot {
        &self.snapshot
    }

    /// Concrete filesystem read receipt.
    #[must_use]
    pub const fn receipt(&self) -> &TraceabilityFilesystemReceipt {
        &self.receipt
    }

    /// Consume into the sealed snapshot and filesystem receipt.
    #[must_use]
    pub fn into_parts(self) -> (TraceabilitySourceSnapshot, TraceabilityFilesystemReceipt) {
        (self.snapshot, self.receipt)
    }
}

#[derive(Debug)]
struct PlannedSource {
    kind: TraceabilitySourceKind,
    relative_path: PathBuf,
    locator: String,
}

fn diagnostic(
    source: impl Into<String>,
    field: TraceabilityFilesystemField,
    reason: impl Into<String>,
) -> TraceabilityFilesystemDiagnostic {
    TraceabilityFilesystemDiagnostic {
        source: source.into(),
        field,
        reason: reason.into(),
    }
}

/// Read, hash, audit, and seal concrete Beads/contract/registry source files.
///
/// Paths are interpreted relative to one canonical root. Parent traversal,
/// absolute paths, final symlinks, non-regular files, root escapes after
/// canonicalization, oversized/racing reads, malformed Beads envelopes, and
/// incomplete source-class coverage all refuse without a partial snapshot.
///
/// # Errors
/// Returns a deterministically ordered [`TraceabilityFilesystemAudit`] naming
/// every source that could not be admitted.
#[allow(clippy::too_many_lines)] // One ordered world-boundary read and seal transaction.
pub fn load_traceability_source_snapshot(
    root: &Path,
    specs: &[TraceabilityFileSpec<'_>],
    limits: TraceabilityFilesystemLimits,
) -> Result<LoadedTraceabilitySourceSnapshot, TraceabilityFilesystemAudit> {
    let mut diagnostics = validate_limits(limits);
    if specs.is_empty() {
        diagnostics.push(diagnostic(
            "<snapshot>",
            TraceabilityFilesystemField::Snapshot,
            "filesystem source snapshot is empty",
        ));
    }
    if specs.len() > limits.max_sources {
        diagnostics.push(diagnostic(
            "<snapshot>",
            TraceabilityFilesystemField::Limits,
            format!(
                "offered {} source specs; configured maximum is {}",
                specs.len(),
                limits.max_sources
            ),
        ));
    }
    if !diagnostics.is_empty() {
        return Err(finish_audit(specs.len(), diagnostics));
    }

    let canonical_root = match std::fs::canonicalize(root) {
        Ok(root) => root,
        Err(error) => {
            return Err(finish_audit(
                specs.len(),
                vec![diagnostic(
                    "<root>",
                    TraceabilityFilesystemField::Root,
                    format!("cannot canonicalize source root: {error}"),
                )],
            ));
        }
    };
    match std::fs::metadata(&canonical_root) {
        Ok(metadata) if metadata.is_dir() => {}
        Ok(_) => {
            return Err(finish_audit(
                specs.len(),
                vec![diagnostic(
                    "<root>",
                    TraceabilityFilesystemField::Root,
                    "source root is not a directory",
                )],
            ));
        }
        Err(error) => {
            return Err(finish_audit(
                specs.len(),
                vec![diagnostic(
                    "<root>",
                    TraceabilityFilesystemField::Root,
                    format!("cannot inspect source root: {error}"),
                )],
            ));
        }
    }

    let mut planned = Vec::new();
    if let Err(error) = planned.try_reserve_exact(specs.len()) {
        return Err(finish_audit(
            specs.len(),
            vec![diagnostic(
                "<snapshot>",
                TraceabilityFilesystemField::Limits,
                format!("allocation refused for source plan: {error}"),
            )],
        ));
    }
    for spec in specs {
        match portable_locator(spec.relative_path) {
            Ok(locator) => planned.push(PlannedSource {
                kind: spec.kind,
                relative_path: spec.relative_path.to_path_buf(),
                locator,
            }),
            Err(reason) => diagnostics.push(diagnostic(
                "<invalid-relative-path>",
                TraceabilityFilesystemField::RelativePath,
                reason,
            )),
        }
    }
    planned.sort_unstable_by(|left, right| {
        left.locator
            .cmp(&right.locator)
            .then_with(|| left.kind.cmp(&right.kind))
    });
    for pair in planned.windows(2) {
        if pair[0].locator == pair[1].locator {
            diagnostics.push(diagnostic(
                &pair[0].locator,
                TraceabilityFilesystemField::RelativePath,
                "duplicate source locator; one file cannot be relabeled across source classes",
            ));
        }
    }
    if !diagnostics.is_empty() {
        return Err(finish_audit(specs.len(), diagnostics));
    }

    let mut receipts = Vec::new();
    if let Err(error) = receipts.try_reserve_exact(planned.len()) {
        return Err(finish_audit(
            specs.len(),
            vec![diagnostic(
                "<snapshot>",
                TraceabilityFilesystemField::Limits,
                format!("allocation refused for source receipts: {error}"),
            )],
        ));
    }
    let mut total_bytes = 0usize;
    for source in &planned {
        match load_one_source(&canonical_root, source, limits, total_bytes) {
            Ok(receipt) => {
                total_bytes = match total_bytes.checked_add(receipt.byte_count) {
                    Some(total) => total,
                    None => {
                        diagnostics.push(diagnostic(
                            &source.locator,
                            TraceabilityFilesystemField::Content,
                            "total source byte count overflowed",
                        ));
                        continue;
                    }
                };
                receipts.push(receipt);
            }
            Err(source_diagnostics) => diagnostics.extend(source_diagnostics),
        }
    }
    if !diagnostics.is_empty() {
        return Err(finish_audit(specs.len(), diagnostics));
    }

    let mut references = Vec::new();
    if let Err(error) = references.try_reserve_exact(receipts.len()) {
        return Err(finish_audit(
            specs.len(),
            vec![diagnostic(
                "<snapshot>",
                TraceabilityFilesystemField::Limits,
                format!("allocation refused for source references: {error}"),
            )],
        ));
    }
    references.extend(receipts.iter().map(|receipt| TraceabilitySource {
        kind: receipt.kind,
        locator: &receipt.locator,
        content_identity: receipt.content_identity,
    }));
    let snapshot = match TraceabilitySourceSnapshot::new(&references) {
        Ok(snapshot) => snapshot,
        Err(audit) => {
            let diagnostics = audit
                .diagnostics
                .into_iter()
                .map(|source| {
                    diagnostic(
                        source.source,
                        TraceabilityFilesystemField::Snapshot,
                        format!("{}: {}", source.field.name(), source.reason),
                    )
                })
                .collect();
            return Err(finish_audit(specs.len(), diagnostics));
        }
    };
    let receipt = TraceabilityFilesystemReceipt {
        adapter_version: TRACEABILITY_FILESYSTEM_ADAPTER_VERSION,
        limits,
        total_bytes,
        sources: receipts,
        source_snapshot_identity: snapshot.identity(),
    };
    Ok(LoadedTraceabilitySourceSnapshot { snapshot, receipt })
}

fn validate_limits(limits: TraceabilityFilesystemLimits) -> Vec<TraceabilityFilesystemDiagnostic> {
    let mut diagnostics = Vec::new();
    for (name, value, hard_max) in [
        ("max_sources", limits.max_sources, MAX_TRACEABILITY_SOURCES),
        (
            "max_source_bytes",
            limits.max_source_bytes,
            MAX_TRACEABILITY_SOURCE_FILE_BYTES,
        ),
        (
            "max_total_bytes",
            limits.max_total_bytes,
            MAX_TRACEABILITY_TOTAL_SOURCE_BYTES,
        ),
        (
            "max_beads_records",
            limits.max_beads_records,
            MAX_TRACEABILITY_BEADS_RECORDS,
        ),
        (
            "max_beads_line_bytes",
            limits.max_beads_line_bytes,
            MAX_TRACEABILITY_BEADS_LINE_BYTES,
        ),
        (
            "max_beads_json_nesting",
            limits.max_beads_json_nesting,
            MAX_TRACEABILITY_BEADS_JSON_NESTING,
        ),
    ] {
        if value == 0 || value > hard_max {
            diagnostics.push(diagnostic(
                "<limits>",
                TraceabilityFilesystemField::Limits,
                format!("{name} must be in 1..={hard_max}; got {value}"),
            ));
        }
    }
    if limits.max_source_bytes > limits.max_total_bytes {
        diagnostics.push(diagnostic(
            "<limits>",
            TraceabilityFilesystemField::Limits,
            "max_source_bytes exceeds max_total_bytes",
        ));
    }
    diagnostics
}

fn portable_locator(path: &Path) -> Result<String, String> {
    if path.as_os_str().is_empty() {
        return Err("source path is empty".to_string());
    }
    let Some(raw_path) = path.to_str() else {
        return Err("source path contains non-UTF-8 bytes".to_string());
    };
    let mut locator = String::new();
    for component in path.components() {
        let Component::Normal(component) = component else {
            return Err(format!(
                "source path {:?} is not a strict root-relative normal path",
                path
            ));
        };
        let Some(component) = component.to_str() else {
            return Err("source path contains non-UTF-8 bytes".to_string());
        };
        if component.is_empty()
            || component.contains('\\')
            || component.chars().any(char::is_control)
        {
            return Err(
                "source path contains an empty, backslash, or control-bearing component"
                    .to_string(),
            );
        }
        if !locator.is_empty() {
            locator.push('/');
        }
        locator.push_str(component);
        if locator.len() > MAX_TRACEABILITY_SOURCE_LOCATOR_BYTES {
            return Err(format!(
                "source locator exceeds {MAX_TRACEABILITY_SOURCE_LOCATOR_BYTES} bytes"
            ));
        }
    }
    if locator.is_empty() {
        return Err("source path has no normal component".to_string());
    }
    if raw_path != locator {
        return Err(
            "source path must use canonical '/' separators without '.' or redundant separators"
                .to_string(),
        );
    }
    Ok(locator)
}

fn load_one_source(
    canonical_root: &Path,
    source: &PlannedSource,
    limits: TraceabilityFilesystemLimits,
    current_total: usize,
) -> Result<TraceabilityFilesystemSourceReceipt, Vec<TraceabilityFilesystemDiagnostic>> {
    let requested = canonical_root.join(&source.relative_path);
    let metadata = match std::fs::symlink_metadata(&requested) {
        Ok(metadata) => metadata,
        Err(error) => {
            return Err(vec![diagnostic(
                &source.locator,
                TraceabilityFilesystemField::Metadata,
                format!("cannot inspect source file: {error}"),
            )]);
        }
    };
    if metadata.file_type().is_symlink() {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Metadata,
            "final source path is a symlink; exact regular files are required",
        )]);
    }
    if !metadata.is_file() {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Metadata,
            "source path is not a regular file",
        )]);
    }
    let canonical_file = match std::fs::canonicalize(&requested) {
        Ok(path) => path,
        Err(error) => {
            return Err(vec![diagnostic(
                &source.locator,
                TraceabilityFilesystemField::Metadata,
                format!("cannot canonicalize source file: {error}"),
            )]);
        }
    };
    if !canonical_file.starts_with(canonical_root) {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Metadata,
            "canonical source file escapes the configured root",
        )]);
    }
    let expected_len = match usize::try_from(metadata.len()) {
        Ok(length) => length,
        Err(_) => {
            return Err(vec![diagnostic(
                &source.locator,
                TraceabilityFilesystemField::Content,
                "source file length does not fit this platform",
            )]);
        }
    };
    if expected_len == 0 {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Content,
            "source file is empty",
        )]);
    }
    if expected_len > limits.max_source_bytes {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Content,
            format!(
                "source file is {expected_len} bytes; configured maximum is {}",
                limits.max_source_bytes
            ),
        )]);
    }
    if current_total
        .checked_add(expected_len)
        .is_none_or(|total| total > limits.max_total_bytes)
    {
        return Err(vec![diagnostic(
            &source.locator,
            TraceabilityFilesystemField::Content,
            format!(
                "source would exceed configured total-byte maximum {}",
                limits.max_total_bytes
            ),
        )]);
    }
    let bytes = match read_bounded(&canonical_file, expected_len, limits.max_source_bytes) {
        Ok(bytes) => bytes,
        Err(reason) => {
            return Err(vec![diagnostic(
                &source.locator,
                TraceabilityFilesystemField::Content,
                reason,
            )]);
        }
    };
    let beads_record_count = match source.kind {
        TraceabilitySourceKind::Beads => match audit_beads_jsonl(&bytes, limits) {
            Ok(records) => Some(records),
            Err(reason) => {
                return Err(vec![diagnostic(
                    &source.locator,
                    TraceabilityFilesystemField::BeadsJsonl,
                    reason,
                )]);
            }
        },
        TraceabilitySourceKind::Contract | TraceabilitySourceKind::Registry => {
            match core::str::from_utf8(&bytes) {
                Ok(text) if !text.contains('\0') => None,
                Ok(_) => {
                    return Err(vec![diagnostic(
                        &source.locator,
                        TraceabilityFilesystemField::Content,
                        "text source contains a NUL byte",
                    )]);
                }
                Err(error) => {
                    return Err(vec![diagnostic(
                        &source.locator,
                        TraceabilityFilesystemField::Content,
                        format!("text source is not UTF-8: {error}"),
                    )]);
                }
            }
        }
    };
    Ok(TraceabilityFilesystemSourceReceipt {
        kind: source.kind,
        locator: source.locator.clone(),
        content_identity: hash_bytes(&bytes),
        byte_count: bytes.len(),
        beads_record_count,
    })
}

fn read_bounded(path: &Path, expected_len: usize, limit: usize) -> Result<Vec<u8>, String> {
    let file = File::open(path).map_err(|error| format!("cannot open source file: {error}"))?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(expected_len)
        .map_err(|error| format!("allocation refused for {expected_len} source bytes: {error}"))?;
    let read_limit = limit
        .checked_add(1)
        .and_then(|value| u64::try_from(value).ok())
        .ok_or_else(|| "bounded read limit overflowed".to_string())?;
    file.take(read_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read source file: {error}"))?;
    if bytes.len() > limit {
        return Err(format!("source exceeded {limit} bytes while being read"));
    }
    if bytes.len() != expected_len {
        return Err(format!(
            "source changed while being read: metadata length {expected_len}, read length {}",
            bytes.len()
        ));
    }
    Ok(bytes)
}

fn audit_beads_jsonl(bytes: &[u8], limits: TraceabilityFilesystemLimits) -> Result<usize, String> {
    core::str::from_utf8(bytes).map_err(|error| format!("Beads JSONL is not UTF-8: {error}"))?;
    let mut ids = BTreeSet::new();
    let mut record_count = 0usize;
    for (line_offset, raw_line) in bytes.split(|byte| *byte == b'\n').enumerate() {
        let raw_line = raw_line.strip_suffix(b"\r").unwrap_or(raw_line);
        let line = trim_ascii_whitespace(raw_line);
        if line.is_empty() {
            continue;
        }
        if line.len() > limits.max_beads_line_bytes {
            return Err(format!(
                "Beads line {} is {} bytes; configured maximum is {}",
                line_offset + 1,
                line.len(),
                limits.max_beads_line_bytes
            ));
        }
        record_count = record_count
            .checked_add(1)
            .ok_or_else(|| "Beads record count overflowed".to_string())?;
        if record_count > limits.max_beads_records {
            return Err(format!(
                "Beads JSONL exceeds configured {}-record maximum",
                limits.max_beads_records
            ));
        }
        validate_json_object_envelope(line, limits.max_beads_json_nesting)
            .map_err(|reason| format!("Beads line {}: {reason}", line_offset + 1))?;
        let id = top_level_bead_id(line)
            .map_err(|reason| format!("Beads line {}: {reason}", line_offset + 1))?;
        if id.len() > MAX_TRACEABILITY_FIELD_BYTES {
            return Err(format!(
                "Beads line {} id is {} bytes; maximum is {MAX_TRACEABILITY_FIELD_BYTES}",
                line_offset + 1,
                id.len()
            ));
        }
        if !ids.insert(id.to_string()) {
            return Err(format!(
                "Beads line {} repeats canonical id {id:?}",
                line_offset + 1
            ));
        }
    }
    if record_count == 0 {
        return Err("Beads JSONL contains no nonblank records".to_string());
    }
    Ok(record_count)
}

fn trim_ascii_whitespace(mut bytes: &[u8]) -> &[u8] {
    while bytes.first().is_some_and(|byte| byte.is_ascii_whitespace()) {
        bytes = &bytes[1..];
    }
    while bytes.last().is_some_and(|byte| byte.is_ascii_whitespace()) {
        bytes = &bytes[..bytes.len() - 1];
    }
    bytes
}

fn validate_json_object_envelope(line: &[u8], max_nesting: usize) -> Result<(), String> {
    if line.first() != Some(&b'{') || line.last() != Some(&b'}') {
        return Err("record must be exactly one JSON object".to_string());
    }
    let mut stack = [0u8; MAX_TRACEABILITY_BEADS_JSON_NESTING];
    let mut depth = 0usize;
    let mut index = 0usize;
    let mut in_string = false;
    while index < line.len() {
        let byte = line[index];
        if in_string {
            match byte {
                b'"' => in_string = false,
                b'\\' => {
                    index += 1;
                    let Some(escaped) = line.get(index).copied() else {
                        return Err("record ends inside a JSON escape".to_string());
                    };
                    match escaped {
                        b'"' | b'\\' | b'/' | b'b' | b'f' | b'n' | b'r' | b't' => {}
                        b'u' => {
                            let end = index.checked_add(5).ok_or_else(|| {
                                "JSON Unicode escape index overflowed".to_string()
                            })?;
                            let Some(hex) = line.get(index + 1..end) else {
                                return Err("short JSON Unicode escape".to_string());
                            };
                            if !hex.iter().all(u8::is_ascii_hexdigit) {
                                return Err("invalid JSON Unicode escape".to_string());
                            }
                            index += 4;
                        }
                        _ => return Err("invalid JSON string escape".to_string()),
                    }
                }
                0x00..=0x1f => {
                    return Err("unescaped control byte in JSON string".to_string());
                }
                _ => {}
            }
            index += 1;
            continue;
        }
        match byte {
            b'"' => in_string = true,
            b'{' | b'[' => {
                if depth >= max_nesting {
                    return Err(format!(
                        "JSON nesting exceeds configured maximum {max_nesting}"
                    ));
                }
                stack[depth] = byte;
                depth += 1;
            }
            b'}' | b']' => {
                if depth == 0 {
                    return Err("JSON container closes without an opener".to_string());
                }
                let opener = stack[depth - 1];
                if !matches!((opener, byte), (b'{', b'}') | (b'[', b']')) {
                    return Err("JSON container delimiters are mismatched".to_string());
                }
                depth -= 1;
                if depth == 0 && index + 1 != line.len() {
                    return Err("record has bytes after its outer JSON object".to_string());
                }
            }
            0x00..=0x1f if !byte.is_ascii_whitespace() => {
                return Err("control byte outside JSON string".to_string());
            }
            _ => {}
        }
        index += 1;
    }
    if in_string {
        return Err("record ends inside a JSON string".to_string());
    }
    if depth != 0 {
        return Err("record ends inside a JSON container".to_string());
    }
    Ok(())
}

fn top_level_bead_id(line: &[u8]) -> Result<&str, String> {
    let mut depth = 0usize;
    let mut index = 0usize;
    let mut found = None;
    while index < line.len() {
        match line[index] {
            b'{' | b'[' => {
                depth += 1;
                index += 1;
            }
            b'}' | b']' => {
                depth = depth.saturating_sub(1);
                index += 1;
            }
            b'"' => {
                let (end, escaped) = json_string_end(line, index)?;
                let previous = previous_non_whitespace(line, index);
                if depth == 1 && matches!(previous, Some(b'{') | Some(b',')) {
                    let mut colon = end + 1;
                    while line
                        .get(colon)
                        .is_some_and(|byte| byte.is_ascii_whitespace())
                    {
                        colon += 1;
                    }
                    if line.get(colon) == Some(&b':') && !escaped && &line[index + 1..end] == b"id"
                    {
                        let mut value = colon + 1;
                        while line
                            .get(value)
                            .is_some_and(|byte| byte.is_ascii_whitespace())
                        {
                            value += 1;
                        }
                        if line.get(value) != Some(&b'"') {
                            return Err("top-level id value must be a JSON string".to_string());
                        }
                        let (value_end, value_escaped) = json_string_end(line, value)?;
                        if value_escaped {
                            return Err(
                                "top-level id must use canonical unescaped UTF-8".to_string()
                            );
                        }
                        let raw = &line[value + 1..value_end];
                        let id = core::str::from_utf8(raw)
                            .map_err(|error| format!("top-level id is not UTF-8: {error}"))?;
                        if id.trim().is_empty() || id.chars().any(char::is_control) {
                            return Err("top-level id is blank or control-bearing".to_string());
                        }
                        if found.replace(id).is_some() {
                            return Err("record repeats the top-level id key".to_string());
                        }
                    }
                }
                index = end + 1;
            }
            _ => index += 1,
        }
    }
    found.ok_or_else(|| "record has no canonical top-level string id".to_string())
}

fn json_string_end(line: &[u8], start: usize) -> Result<(usize, bool), String> {
    let mut index = start + 1;
    let mut escaped_any = false;
    while index < line.len() {
        match line[index] {
            b'"' => return Ok((index, escaped_any)),
            b'\\' => {
                escaped_any = true;
                index += 1;
                let Some(escaped) = line.get(index).copied() else {
                    return Err("record ends inside a JSON escape".to_string());
                };
                if escaped == b'u' {
                    let end = index
                        .checked_add(5)
                        .ok_or_else(|| "JSON Unicode escape index overflowed".to_string())?;
                    let Some(hex) = line.get(index + 1..end) else {
                        return Err("short JSON Unicode escape".to_string());
                    };
                    if !hex.iter().all(u8::is_ascii_hexdigit) {
                        return Err("invalid JSON Unicode escape".to_string());
                    }
                    index += 4;
                }
            }
            0x00..=0x1f => return Err("unescaped control byte in JSON string".to_string()),
            _ => {}
        }
        index += 1;
    }
    Err("record ends inside a JSON string".to_string())
}

fn previous_non_whitespace(line: &[u8], index: usize) -> Option<u8> {
    line[..index]
        .iter()
        .rev()
        .copied()
        .find(|byte| !byte.is_ascii_whitespace())
}

fn finish_audit(
    total: usize,
    mut diagnostics: Vec<TraceabilityFilesystemDiagnostic>,
) -> TraceabilityFilesystemAudit {
    diagnostics.sort_by(|left, right| {
        left.source
            .cmp(&right.source)
            .then_with(|| left.field.cmp(&right.field))
            .then_with(|| left.reason.cmp(&right.reason))
    });
    diagnostics.dedup();
    TraceabilityFilesystemAudit { total, diagnostics }
}
