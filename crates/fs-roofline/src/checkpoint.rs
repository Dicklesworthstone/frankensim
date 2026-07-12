//! Bounded, authenticated staleness-history checkpoints (bead vm3i).
//!
//! The exhaustive verifier ([`crate::staleness_at`]) re-validates every
//! matching tune row on every probe: op fetch, protocol parse, dependency
//! receipt, manifest membership, then a per-sibling `tune_get` walk — for H
//! retained runs of a K-kernel registry that is O(H·(K+4)) ledger queries,
//! and sustained CI turns it into an unbounded resource sink.
//!
//! A [`StalenessCheckpoint`] is one exhaustive pass SEALED: per covered row
//! it retains a content hash (over key + params + measured bytes), the
//! build identity, the op-bound dependency-receipt digests, the recorded
//! timestamp, and the verdict. Checkpoints live in the ledger's own tune
//! table under a reserved kernel name, machine-keyed exactly like the
//! production rows they cover, and chain through a domain-separated digest:
//! `digest_i = H(domain, prev_digest_bytes ‖ body_bytes)`. Authentication
//! is the chain over verifier-produced content — operator-observed,
//! tamper-EVIDENT, no cryptographic unforgeability claim (the workspace
//! no-crypto no-claim applies: anyone who can write tune rows can mint a
//! chain; what they cannot do is alter covered history under an EXISTING
//! chain without detection).
//!
//! The fast path ([`staleness_at_checkpointed`]) costs a constant number of
//! ledger queries plus the exhaustive cost of only the DELTA (rows newer
//! than the checkpoint):
//! 1. fetch checkpoint rows (1 bounded query), verify ordinal contiguity
//!    and the digest chain — any break FAILS CLOSED to the exhaustive path;
//! 2. fetch the kernel's rows (1 bounded query) and classify the cheap
//!    lattice prefix identically to the exhaustive path;
//! 3. every covered row must hash to its checkpointed entry (tamper and
//!    rollback → `CorruptEvidence`); every checkpointed row must still
//!    exist (removal → `CorruptEvidence`); corrupt verdicts are TOMBSTONES
//!    (once corrupt, forever corrupt — re-checkpointing preserves them, so
//!    compaction can refuse but never un-corrupt);
//! 4. covered valid rows replay the build/dependency scan from checkpointed
//!    metadata against the CURRENT build and dependency binding — no op
//!    fetches; delta rows run the full exhaustive validator.

use fs_ledger::{Ledger, LedgerError};

use crate::{
    BuildRowScan, DependencyReceiptBinding, RowSelection, Staleness, ValidatedRooflineRow,
    classify_scanned_rows, executable_build_identity, roofline_machine_key, select_matching_rows,
    validate_roofline_row,
};

/// Reserved tune-table kernel prefix for checkpoint rows. Distinct from any
/// production kernel name, so production row queries never see checkpoints
/// and vice versa.
const CHECKPOINT_KERNEL_PREFIX: &str = "roofline-staleness-checkpoint:";
const CHECKPOINT_SHAPE_PREFIX: &str = "roofline-ckpt-v1:";
const CHECKPOINT_SCHEMA: &str = "fs-roofline-staleness-checkpoint-v1";
const CHECKPOINT_CHAIN_DOMAIN: &str = "org.frankensim.fs-roofline.staleness-checkpoint-chain.v1";
const ROW_CONTENT_DOMAIN: &str = "org.frankensim.fs-roofline.staleness-row-content.v1";

/// Per-row verdict sealed into a checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckpointVerdict {
    /// The row passed the exhaustive verifier at checkpoint time.
    Valid,
    /// The row failed verification — a permanent tombstone.
    Corrupt,
}

impl CheckpointVerdict {
    fn as_str(self) -> &'static str {
        match self {
            CheckpointVerdict::Valid => "valid",
            CheckpointVerdict::Corrupt => "corrupt",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        match text {
            "valid" => Some(CheckpointVerdict::Valid),
            "corrupt" => Some(CheckpointVerdict::Corrupt),
            _ => None,
        }
    }
}

/// One covered row inside a checkpoint body.
#[derive(Debug, Clone)]
struct CheckpointEntry {
    shape_class: String,
    row_hash: fs_blake3::ContentHash,
    build: fs_blake3::ContentHash,
    dep_digest: fs_blake3::ContentHash,
    dep_artifact: fs_blake3::ContentHash,
    recorded_at_ns: i64,
    verdict: CheckpointVerdict,
}

/// A verified checkpoint chain head plus its covered entries.
#[derive(Debug)]
struct VerifiedCheckpoint {
    entries: Vec<CheckpointEntry>,
    /// Digest of the newest checkpoint (the chain tip).
    tip_digest: fs_blake3::ContentHash,
    /// Number of checkpoints in the chain (== next free ordinal).
    len: u64,
}

/// Receipt returned by [`checkpoint_staleness_history`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckpointReceipt {
    /// Ordinal of the checkpoint written (0-based, contiguous per chain).
    pub ordinal: u64,
    /// Chained digest sealing this checkpoint.
    pub digest: fs_blake3::ContentHash,
    /// Rows covered.
    pub rows: usize,
    /// Rows sealed as corrupt tombstones (including inherited ones).
    pub corrupt_rows: usize,
}

/// All-zero placeholder hash for tombstone entries (no verified metadata).
fn zero_hash() -> fs_blake3::ContentHash {
    fs_blake3::ContentHash::from_hex(&"0".repeat(64)).expect("zero hash parses")
}

fn checkpoint_kernel(kernel: &str) -> String {
    format!("{CHECKPOINT_KERNEL_PREFIX}{kernel}")
}

/// One sealed entry from an exhaustive-verifier result. `None` (failed or
/// inherited-tombstone) seals a permanent Corrupt tombstone with zeroed
/// metadata.
fn seal_entry(
    row: &fs_ledger::TuneRow,
    validated: Option<ValidatedRooflineRow>,
) -> CheckpointEntry {
    match validated {
        Some(v) => CheckpointEntry {
            shape_class: row.shape_class.clone(),
            row_hash: row_content_hash(row),
            build: v.build_identity,
            dep_digest: v.dependency_receipt_digest,
            dep_artifact: v.dependency_receipt_artifact,
            recorded_at_ns: v.recorded_at_ns,
            verdict: CheckpointVerdict::Valid,
        },
        None => CheckpointEntry {
            shape_class: row.shape_class.clone(),
            row_hash: row_content_hash(row),
            build: zero_hash(),
            dep_digest: zero_hash(),
            dep_artifact: zero_hash(),
            recorded_at_ns: 0,
            verdict: CheckpointVerdict::Corrupt,
        },
    }
}

fn checkpoint_shape_class(version: &str, ordinal: u64) -> String {
    format!("{CHECKPOINT_SHAPE_PREFIX}{version}:{ordinal:08}")
}

fn checkpoint_params(ordinal: u64, digest: fs_blake3::ContentHash) -> String {
    format!("{{\"schema\":\"{CHECKPOINT_SCHEMA}\",\"ordinal\":{ordinal},\"digest\":\"{digest}\"}}")
}

/// Content hash binding one tune row's full stored identity.
fn row_content_hash(row: &fs_ledger::TuneRow) -> fs_blake3::ContentHash {
    let mut material = Vec::new();
    for part in [
        row.kernel.as_bytes(),
        row.shape_class.as_bytes(),
        row.machine.as_slice(),
        row.params.as_bytes(),
        row.measured.as_bytes(),
    ] {
        material.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
        material.extend_from_slice(part);
    }
    fs_blake3::hash_domain(ROW_CONTENT_DOMAIN, &material)
}

fn entry_json(entry: &CheckpointEntry) -> String {
    format!(
        "{{\"shape_class\":\"{}\",\"row_hash\":\"{}\",\"build\":\"{}\",\"dep_digest\":\"{}\",\"dep_artifact\":\"{}\",\"recorded_at_ns\":{},\"verdict\":\"{}\"}}",
        entry.shape_class,
        entry.row_hash,
        entry.build,
        entry.dep_digest,
        entry.dep_artifact,
        entry.recorded_at_ns,
        entry.verdict.as_str(),
    )
}

fn body_json(
    kernel: &str,
    version: &str,
    ordinal: u64,
    prev: Option<fs_blake3::ContentHash>,
    entries: &[CheckpointEntry],
) -> String {
    let prev_text = prev.map_or_else(|| "null".to_string(), |p| format!("\"{p}\""));
    let rows = entries.iter().map(entry_json).collect::<Vec<_>>().join(",");
    format!(
        "{{\"schema\":\"{CHECKPOINT_SCHEMA}\",\"kernel\":\"{kernel}\",\"version\":\"{version}\",\"ordinal\":{ordinal},\"prev\":{prev_text},\"rows\":[{rows}]}}",
    )
}

fn chain_digest(
    prev: Option<fs_blake3::ContentHash>,
    body: &str,
) -> fs_blake3::ContentHash {
    let mut material = Vec::new();
    if let Some(prev) = prev {
        material.extend_from_slice(prev.as_bytes());
    }
    material.extend_from_slice(body.as_bytes());
    fs_blake3::hash_domain(CHECKPOINT_CHAIN_DOMAIN, &material)
}

/// Parsed checkpoint body fields.
struct ParsedBody {
    kernel: String,
    version: String,
    ordinal: u64,
    prev: Option<fs_blake3::ContentHash>,
    entries: Vec<CheckpointEntry>,
}

/// Strict parse of one checkpoint body; byte-exact round trip enforced.
fn parse_body(text: &str) -> Option<ParsedBody> {
    let rest = text.strip_prefix(&format!("{{\"schema\":\"{CHECKPOINT_SCHEMA}\",\"kernel\":\""))?;
    let (kernel, rest) = rest.split_once("\",\"version\":\"")?;
    let (version, rest) = rest.split_once("\",\"ordinal\":")?;
    let (ordinal_text, rest) = rest.split_once(",\"prev\":")?;
    let ordinal: u64 = ordinal_text.parse().ok()?;
    let (prev_text, rest) = rest.split_once(",\"rows\":[")?;
    let prev = if prev_text == "null" {
        None
    } else {
        Some(fs_blake3::ContentHash::from_hex(
            prev_text.strip_prefix('"')?.strip_suffix('"')?,
        )?)
    };
    let rows_text = rest.strip_suffix("]}")?;
    let mut entries = Vec::new();
    if !rows_text.is_empty() {
        for raw in rows_text.split("},") {
            let raw = if raw.ends_with('}') {
                raw.to_string()
            } else {
                format!("{raw}}}")
            };
            let inner = raw.strip_prefix("{\"shape_class\":\"")?;
            let (shape_class, inner) = inner.split_once("\",\"row_hash\":\"")?;
            let (row_hash, inner) = inner.split_once("\",\"build\":\"")?;
            let (build, inner) = inner.split_once("\",\"dep_digest\":\"")?;
            let (dep_digest, inner) = inner.split_once("\",\"dep_artifact\":\"")?;
            let (dep_artifact, inner) = inner.split_once("\",\"recorded_at_ns\":")?;
            let (recorded_text, inner) = inner.split_once(",\"verdict\":\"")?;
            let verdict_text = inner.strip_suffix("\"}")?;
            entries.push(CheckpointEntry {
                shape_class: shape_class.to_string(),
                row_hash: fs_blake3::ContentHash::from_hex(row_hash)?,
                build: fs_blake3::ContentHash::from_hex(build)?,
                dep_digest: fs_blake3::ContentHash::from_hex(dep_digest)?,
                dep_artifact: fs_blake3::ContentHash::from_hex(dep_artifact)?,
                recorded_at_ns: recorded_text.parse().ok()?,
                verdict: CheckpointVerdict::parse(verdict_text)?,
            });
        }
    }
    // Byte-exact round trip: non-canonical spellings are refused.
    (body_json(kernel, version, ordinal, prev, &entries) == text).then_some(ParsedBody {
        kernel: kernel.to_string(),
        version: version.to_string(),
        ordinal,
        prev,
        entries,
    })
}

/// Chain-load outcome: absent, present-but-unverifiable, or verified.
enum ChainState {
    /// No checkpoint rows exist for this (kernel, version, machine).
    Empty,
    /// Rows exist but the chain fails verification (parse failure, ordinal
    /// gap, digest mismatch, params/shape inconsistency).
    Broken,
    /// The full chain verified; carries the newest checkpoint's entries.
    Verified(VerifiedCheckpoint),
}

/// Load and verify the checkpoint chain for `(kernel, version, machine)`.
/// Callers FAIL CLOSED on both `Empty` and `Broken`.
fn load_chain(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    machine_key: [u8; 40],
) -> Result<ChainState, LedgerError> {
    let rows = ledger.tune_rows(&checkpoint_kernel(kernel))?;
    let shape_prefix = format!("{CHECKPOINT_SHAPE_PREFIX}{version}:");
    let mut chain: Vec<&fs_ledger::TuneRow> = rows
        .iter()
        .filter(|r| r.machine == machine_key && r.shape_class.starts_with(&shape_prefix))
        .collect();
    if chain.is_empty() {
        return Ok(ChainState::Empty);
    }
    chain.sort_by(|a, b| a.shape_class.cmp(&b.shape_class));
    let mut prev: Option<fs_blake3::ContentHash> = None;
    let mut latest_entries: Option<Vec<CheckpointEntry>> = None;
    for (index, row) in chain.iter().enumerate() {
        let Some(body) = parse_body(&row.measured) else {
            return Ok(ChainState::Broken);
        };
        let expected_digest = chain_digest(prev, &row.measured);
        let params_ok = row.params == checkpoint_params(body.ordinal, expected_digest);
        if body.kernel != kernel
            || body.version != version
            || body.ordinal != index as u64
            || body.prev != prev
            || !params_ok
            || row.shape_class != checkpoint_shape_class(version, body.ordinal)
        {
            return Ok(ChainState::Broken);
        }
        prev = Some(expected_digest);
        latest_entries = Some(body.entries);
    }
    match (latest_entries, prev) {
        (Some(entries), Some(tip_digest)) => Ok(ChainState::Verified(VerifiedCheckpoint {
            entries,
            tip_digest,
            len: chain.len() as u64,
        })),
        _ => Ok(ChainState::Broken),
    }
}

/// Run the exhaustive verifier once over the current matching rows and seal
/// the results as the next checkpoint in the chain. Prior tombstones are
/// PRESERVED: a row sealed corrupt in any earlier checkpoint remains corrupt
/// here even if its stored bytes now re-verify (never un-corrupt).
///
/// # Errors
/// Ledger errors propagate. A lattice-prefix verdict (no rows for this
/// machine/baseline) is an [`LedgerError::Invalid`] refusal — there is
/// nothing to seal.
pub fn checkpoint_staleness_history(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
) -> Result<CheckpointReceipt, LedgerError> {
    checkpoint_history_with_dependency(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        DependencyReceiptBinding::current().ok(),
    )
}

pub(crate) fn checkpoint_history_with_dependency(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<CheckpointReceipt, LedgerError> {
    let machine_key = roofline_machine_key(current_fingerprint, current_baseline);
    let matching = match select_matching_rows(
        ledger,
        kernel,
        version,
        current_fingerprint,
        Some(current_baseline),
    )? {
        RowSelection::Rows(rows) => rows,
        RowSelection::Verdict(v) => {
            return Err(LedgerError::Invalid {
                field: "staleness-checkpoint".to_string(),
                problem: format!(
                    "nothing to seal: row selection classified {v:?} before per-row verification"
                ),
            });
        }
    };

    // Inherit tombstones and the chain tip from the existing chain. A chain
    // that EXISTS but fails verification blocks sealing entirely: appending
    // to a broken chain would launder it.
    let (prior_tombstones, next_ordinal, prev_digest) =
        match load_chain(ledger, kernel, version, machine_key)? {
            ChainState::Verified(verified) => (
                verified
                    .entries
                    .iter()
                    .filter(|e| e.verdict == CheckpointVerdict::Corrupt)
                    .map(|e| e.shape_class.clone())
                    .collect::<Vec<String>>(),
                verified.len,
                Some(verified.tip_digest),
            ),
            ChainState::Empty => (Vec::new(), 0, None),
            ChainState::Broken => {
                return Err(LedgerError::Invalid {
                    field: "staleness-checkpoint".to_string(),
                    problem: "existing checkpoint chain fails verification; refusing to \
                              extend it (a broken chain is permanent evidence, never \
                              overwritten)"
                        .to_string(),
                });
            }
        };

    let mut entries = Vec::with_capacity(matching.len());
    let mut corrupt_rows = 0usize;
    for row in &matching {
        let inherited_corrupt = prior_tombstones.contains(&row.shape_class);
        let validated = if inherited_corrupt {
            None
        } else {
            validate_roofline_row(
                ledger,
                row,
                kernel,
                version,
                current_fingerprint,
                current_baseline,
                expected_dependency,
            )?
        };
        if validated.is_none() {
            corrupt_rows += 1;
        }
        entries.push(seal_entry(row, validated));
    }
    entries.sort_by(|a, b| a.shape_class.cmp(&b.shape_class));

    let body = body_json(kernel, version, next_ordinal, prev_digest, &entries);
    let digest = chain_digest(prev_digest, &body);
    let params = checkpoint_params(next_ordinal, digest);
    // Append-only: a colliding ordinal must never overwrite sealed history.
    ledger.tune_put_if_absent(
        &checkpoint_kernel(kernel),
        &checkpoint_shape_class(version, next_ordinal),
        &machine_key,
        &params,
        &body,
    )?;
    let stored = ledger
        .tune_get(
            &checkpoint_kernel(kernel),
            &checkpoint_shape_class(version, next_ordinal),
            &machine_key,
        )?
        .ok_or_else(|| LedgerError::Invalid {
            field: "staleness-checkpoint".to_string(),
            problem: "checkpoint row missing immediately after insert".to_string(),
        })?;
    if stored.measured != body || stored.params != params {
        return Err(LedgerError::Invalid {
            field: "staleness-checkpoint".to_string(),
            problem: format!(
                "checkpoint ordinal {next_ordinal} already sealed with different content; \
                 refusing to overwrite chain history"
            ),
        });
    }
    Ok(CheckpointReceipt {
        ordinal: next_ordinal,
        digest,
        rows: entries.len(),
        corrupt_rows,
    })
}

/// Checkpoint-accelerated staleness: identical verdict lattice to
/// [`crate::staleness_at`], with covered rows verified against the sealed
/// chain instead of re-walking their manifests. Falls back to the exhaustive
/// path whenever no verified chain exists (fail closed).
///
/// # Errors
/// Ledger and executable-identity errors propagate.
pub fn staleness_at_checkpointed(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
) -> Result<Staleness, LedgerError> {
    let current_build = executable_build_identity()?;
    staleness_at_checkpointed_with(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        observed_wall_ns,
        current_build,
        DependencyReceiptBinding::current().ok(),
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn staleness_at_checkpointed_with(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Staleness, LedgerError> {
    let matching = match select_matching_rows(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
    )? {
        RowSelection::Verdict(v) => return Ok(v),
        RowSelection::Rows(rows) => rows,
    };
    let baseline = current_baseline.expect("baseline present when rows match");
    let machine_key = roofline_machine_key(current_fingerprint, baseline);

    let verified = match load_chain(ledger, kernel, version, machine_key)? {
        ChainState::Verified(v) => v,
        // No chain, or a chain that fails verification: FAIL CLOSED to the
        // exhaustive per-row path over the already-selected rows.
        ChainState::Empty | ChainState::Broken => {
            return exhaustive_scan(
                ledger,
                &matching,
                kernel,
                version,
                current_fingerprint,
                baseline,
                observed_wall_ns,
                current_build,
                expected_dependency,
            );
        }
    };

    let mut build_scan = BuildRowScan::default();
    let mut covered: std::collections::BTreeMap<&str, &CheckpointEntry> = verified
        .entries
        .iter()
        .map(|e| (e.shape_class.as_str(), e))
        .collect();
    for row in &matching {
        if let Some(entry) = covered.remove(row.shape_class.as_str()) {
            if entry.verdict == CheckpointVerdict::Corrupt {
                // Tombstone: permanently corrupt.
                return Ok(Staleness::CorruptEvidence);
            }
            if row_content_hash(row) != entry.row_hash {
                // Tampered since sealing.
                return Ok(Staleness::CorruptEvidence);
            }
            let dependency_matches_current = expected_dependency.is_some_and(|expected| {
                entry.dep_artifact == expected.artifact_hash
                    && entry.dep_digest == expected.domain_digest
            });
            let replayed = ValidatedRooflineRow {
                build_identity: entry.build,
                recorded_at_ns: entry.recorded_at_ns,
                dependency_matches_current,
                dependency_receipt_digest: entry.dep_digest,
                dependency_receipt_artifact: entry.dep_artifact,
            };
            if !build_scan.observe(&replayed, current_build) {
                return Ok(Staleness::CorruptEvidence);
            }
        } else {
            // Delta row (newer than the checkpoint): full validation.
            let Some(validated) = validate_roofline_row(
                ledger,
                row,
                kernel,
                version,
                current_fingerprint,
                baseline,
                expected_dependency,
            )?
            else {
                return Ok(Staleness::CorruptEvidence);
            };
            if !build_scan.observe(&validated, current_build) {
                return Ok(Staleness::CorruptEvidence);
            }
        }
    }
    if !covered.is_empty() {
        // A sealed row vanished from the ledger: rollback/removal.
        return Ok(Staleness::CorruptEvidence);
    }
    Ok(classify_scanned_rows(build_scan, observed_wall_ns))
}

/// The exhaustive per-row scan over pre-selected rows (shared shape with
/// `staleness_at_with_build_and_dependency`'s tail, reused by the fallback).
#[allow(clippy::too_many_arguments)]
fn exhaustive_scan(
    ledger: &Ledger,
    matching: &[fs_ledger::TuneRow],
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    baseline: fs_blake3::ContentHash,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Staleness, LedgerError> {
    let mut build_scan = BuildRowScan::default();
    for row in matching {
        let Some(validated) = validate_roofline_row(
            ledger,
            row,
            kernel,
            version,
            current_fingerprint,
            baseline,
            expected_dependency,
        )?
        else {
            return Ok(Staleness::CorruptEvidence);
        };
        if !build_scan.observe(&validated, current_build) {
            return Ok(Staleness::CorruptEvidence);
        }
    }
    Ok(classify_scanned_rows(build_scan, observed_wall_ns))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::default_registry;
    use crate::production::{ProductionProbe, ProductionRunConfig};
    use crate::{
        AxisBaselinePolicy, BaselineAxes, BaselineCandidate, BaselineIdentity, MachineAxes,
        STALENESS_MAX_AGE_NS, promote_baseline,
    };

    const FINGERPRINT: u64 = 0xC4EC;

    fn temp_db(tag: &str) -> String {
        static NEXT: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-roofline-ckpt-{tag}-{}-{n}.db",
                std::process::id()
            ))
            .display()
            .to_string()
    }

    fn cleanup_db(path: &str) {
        for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
            let _ = std::fs::remove_file(format!("{path}{suffix}"));
        }
    }

    fn synthetic_axes(fingerprint: u64) -> MachineAxes {
        // Roofs far above any real machine (bead xjhz): cache-resident test
        // kernels must never outrun the fixture roof.
        MachineAxes {
            fingerprint,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100_000.0,
            bandwidth_all_core_gbs: 400_000.0,
            peak_single_gflops: 50_000.0,
            peak_all_core_gflops: 300_000.0,
        }
    }

    fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        let identity =
            BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.checkpoint-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic checkpoint fixture",
            20_000,
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    const CONFIG: ProductionRunConfig = ProductionRunConfig {
        n: 1 << 10,
        warmup: 0,
        reps: 1,
    };

    const TEST_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"checkpoint-battery\"}";

    fn test_binding() -> DependencyReceiptBinding {
        let digest = fs_blake3::hash_domain(
            fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
            TEST_DEPGRAPH_RECEIPT.as_bytes(),
        );
        DependencyReceiptBinding::from_parts(TEST_DEPGRAPH_RECEIPT, digest)
            .expect("test receipt digest agrees")
    }

    struct Fixture {
        ledger: Ledger,
        baseline: BaselineAxes,
        kernels: Vec<(String, String)>,
        recorded_at: i64,
    }

    /// Record one sealed receipt-backed production run into `ledger`.
    fn record_one_run(ledger: &Ledger) -> (Vec<(String, String)>, i64) {
        let axes = synthetic_axes(FINGERPRINT);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_test_receipt(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                TEST_DEPGRAPH_RECEIPT,
            )
            .expect("sealed checkpoint fixture");
        assert!(run.citation_eligible());
        let kernels = run
            .results()
            .iter()
            .map(|result| (result.kernel.clone(), result.version.clone()))
            .collect();
        let op = run.record(ledger).expect("record checkpoint fixture");
        let recorded_at = ledger
            .op(op)
            .unwrap()
            .expect("recorded op")
            .t_end
            .expect("finished op");
        (kernels, recorded_at)
    }

    fn fixture(db: &str) -> Fixture {
        let ledger = Ledger::open(db).expect("open ledger");
        let (kernels, recorded_at) = record_one_run(&ledger);
        let axes = synthetic_axes(FINGERPRINT);
        let (baseline, _) = trusted_baseline(&axes);
        Fixture {
            ledger,
            baseline,
            kernels,
            recorded_at,
        }
    }

    fn seal(fx: &Fixture) -> CheckpointReceipt {
        checkpoint_history_with_dependency(
            &fx.ledger,
            &fx.kernels[0].0,
            &fx.kernels[0].1,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect("seal checkpoint")
    }

    fn fast(fx: &Fixture, kernel: &str, version: &str, at: i64) -> Staleness {
        staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            at,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("checkpointed probe")
    }

    fn exhaustive(fx: &Fixture, kernel: &str, version: &str, at: i64) -> Staleness {
        crate::staleness_at_with_build_and_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            at,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("exhaustive probe")
    }

    /// The stored roofline row for one production kernel.
    fn production_row(ledger: &Ledger, kernel: &str) -> fs_ledger::TuneRow {
        let mut rows: Vec<_> = ledger
            .tune_rows(kernel)
            .expect("tune rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert!(!rows.is_empty(), "expected roofline rows for {kernel}");
        rows.pop().expect("row")
    }

    #[test]
    fn checkpointed_verdicts_match_exhaustive_across_the_lattice() {
        let db = temp_db("equivalence");
        let fx = fixture(&db);
        let receipt = seal(&fx);
        assert_eq!(receipt.ordinal, 0);
        assert_eq!(receipt.corrupt_rows, 0);
        assert!(receipt.rows >= 1);

        let (kernel, version) = &fx.kernels[0];
        // Fresh / Expired / ClockRollback classification parity.
        for at in [
            fx.recorded_at + 1,
            fx.recorded_at + STALENESS_MAX_AGE_NS,
            fx.recorded_at + STALENESS_MAX_AGE_NS + 1,
            fx.recorded_at - 1,
        ] {
            assert_eq!(
                fast(&fx, kernel, version, at),
                exhaustive(&fx, kernel, version, at),
                "fast/exhaustive divergence at offset {}",
                at - fx.recorded_at
            );
        }
        assert_eq!(fast(&fx, kernel, version, fx.recorded_at + 1), Staleness::Fresh);
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + STALENESS_MAX_AGE_NS + 1),
            Staleness::Expired
        );
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at - 1),
            Staleness::ClockRollback
        );

        // Lattice-prefix verdicts decided before the chain is even loaded.
        assert_eq!(
            fast(&fx, "never-measured-kernel", version, fx.recorded_at + 1),
            Staleness::NeverMeasured
        );
        let wrong_fp = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT + 1,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("wrong-fingerprint probe");
        assert_eq!(wrong_fp, Staleness::FingerprintDrift);
        let no_baseline = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            None,
            fx.recorded_at + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("no-baseline probe");
        assert_eq!(no_baseline, Staleness::BaselineUnavailable);

        // Foreign build (injected): both paths must agree on BuildDrift.
        let foreign_build = fs_blake3::hash_domain("fs-roofline.ckpt-test-build.v1", b"other");
        let fast_foreign = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            foreign_build,
            Some(test_binding()),
        )
        .expect("foreign-build fast probe");
        let exhaustive_foreign = crate::staleness_at_with_build_and_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            foreign_build,
            Some(test_binding()),
        )
        .expect("foreign-build exhaustive probe");
        assert_eq!(fast_foreign, exhaustive_foreign);
        assert_eq!(fast_foreign, Staleness::BuildDrift);
        cleanup_db(&db);
    }

    #[test]
    fn fast_path_costs_two_reads_and_undercuts_exhaustive() {
        let db = temp_db("budget");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];

        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, fx.recorded_at + 1), Staleness::Fresh);
        let fast_reads = fx.ledger.read_queries() - before;

        let before = fx.ledger.read_queries();
        assert_eq!(
            exhaustive(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        let exhaustive_reads = fx.ledger.read_queries() - before;

        // Fully covered fast path: one tune_rows for the production kernel,
        // one for the checkpoint chain. Nothing per-row.
        assert_eq!(
            fast_reads, 2,
            "covered fast path must not scale with history"
        );
        assert!(
            exhaustive_reads > fast_reads,
            "exhaustive ({exhaustive_reads}) must cost more than checkpointed ({fast_reads})"
        );
        cleanup_db(&db);
    }

    #[test]
    fn tampered_production_row_is_corrupt_under_the_checkpoint() {
        let db = temp_db("tamper");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];
        assert_eq!(fast(&fx, kernel, version, fx.recorded_at + 1), Staleness::Fresh);

        let row = production_row(&fx.ledger, kernel);
        let forged = row.measured.replace("\"dispersion\":", "\"dispersion\": ");
        assert_ne!(forged, row.measured);
        fx.ledger
            .tune_put(&row.kernel, &row.shape_class, &row.machine, &row.params, &forged)
            .expect("overwrite row");

        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence,
            "content-hash mismatch against the sealed entry must be corrupt"
        );
        cleanup_db(&db);
    }

    #[test]
    fn covered_row_missing_from_the_ledger_is_corrupt() {
        // Rollback simulation: rebuild a ledger holding the checkpoint chain
        // and the second run's rows but MISSING the first run's row — as if
        // history was truncated underneath the seal.
        let db = temp_db("rollback-src");
        let fx = fixture(&db);
        let (kernels2, recorded_at2) = record_one_run(&fx.ledger);
        assert_eq!(fx.kernels, kernels2);
        seal(&fx); // covers both runs' rows
        let (kernel, version) = &fx.kernels[0];
        assert_eq!(fast(&fx, kernel, version, recorded_at2 + 1), Staleness::Fresh);

        let replay_db = temp_db("rollback-dst");
        let replay = Ledger::open(&replay_db).expect("open replay ledger");
        // Copy the checkpoint chain verbatim.
        for row in fx
            .ledger
            .tune_rows(&checkpoint_kernel(kernel))
            .expect("chain rows")
        {
            replay
                .tune_put(&row.kernel, &row.shape_class, &row.machine, &row.params, &row.measured)
                .expect("copy chain row");
        }
        // Copy only the NEWER of the two production rows.
        let mut rows: Vec<_> = fx
            .ledger
            .tune_rows(kernel)
            .expect("production rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert_eq!(rows.len(), 2, "two runs leave two rows");
        rows.sort_by_key(|row| row.shape_class.clone());
        let kept = rows.pop().expect("newest row");
        replay
            .tune_put(&kept.kernel, &kept.shape_class, &kept.machine, &kept.params, &kept.measured)
            .expect("copy surviving row");

        let verdict = staleness_at_checkpointed_with(
            &replay,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            recorded_at2 + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("rollback probe");
        assert_eq!(
            verdict,
            Staleness::CorruptEvidence,
            "a sealed row vanishing from history must be corrupt, not silently fresh"
        );
        cleanup_db(&db);
        cleanup_db(&replay_db);
    }

    #[test]
    fn tampered_chain_fails_closed_to_exhaustive_and_blocks_sealing() {
        let db = temp_db("chain-tamper");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];

        // Flip bytes inside the sealed body. The chain no longer verifies.
        let chain_kernel = checkpoint_kernel(kernel);
        let row = fx
            .ledger
            .tune_rows(&chain_kernel)
            .expect("chain rows")
            .pop()
            .expect("chain row");
        let forged = row.measured.replace("\"verdict\":\"valid\"", "\"verdict\":\"corrupt\"");
        assert_ne!(forged, row.measured, "fixture must have a valid entry to flip");
        fx.ledger
            .tune_put(&row.kernel, &row.shape_class, &row.machine, &row.params, &forged)
            .expect("tamper chain");

        // Fail closed: the fast path falls back to the exhaustive verdict.
        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, fx.recorded_at + 1), Staleness::Fresh);
        let fallback_reads = fx.ledger.read_queries() - before;
        assert!(
            fallback_reads > 2,
            "broken chain must trigger the exhaustive fallback, not the covered path"
        );

        // And extending a broken chain is refused outright.
        let refusal = checkpoint_history_with_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("sealing over a broken chain must refuse");
        assert!(matches!(refusal, LedgerError::Invalid { .. }));
        cleanup_db(&db);
    }

    #[test]
    fn ordinal_gaps_and_duplicate_ordinals_break_the_chain() {
        let db = temp_db("gap");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];
        let chain_kernel = checkpoint_kernel(kernel);
        let row = fx
            .ledger
            .tune_rows(&chain_kernel)
            .expect("chain rows")
            .pop()
            .expect("chain row");

        // Insert a structurally plausible row at a gapped ordinal (2 with no 1).
        fx.ledger
            .tune_put(
                &row.kernel,
                &checkpoint_shape_class(version, 2),
                &row.machine,
                &row.params,
                &row.measured,
            )
            .expect("gapped row");
        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, fx.recorded_at + 1), Staleness::Fresh);
        assert!(
            fx.ledger.read_queries() - before > 2,
            "gapped chain must fail closed to exhaustive"
        );
        cleanup_db(&db);
    }

    #[test]
    fn tombstones_survive_row_restoration_and_resealing() {
        let db = temp_db("tombstone");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let original = production_row(&fx.ledger, kernel);

        // Tamper, then seal: the exhaustive pass records tombstones.
        let forged = original.measured.replace("\"dispersion\":", "\"dispersion\": ");
        assert_ne!(forged, original.measured);
        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &forged,
            )
            .expect("tamper row");
        let sealed = seal(&fx);
        assert!(sealed.corrupt_rows >= 1, "tampered history must seal tombstones");

        // Restore the original bytes: the row would re-verify exhaustively...
        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &original.measured,
            )
            .expect("restore row");
        assert_eq!(
            exhaustive(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh,
            "precondition: restored bytes verify exhaustively"
        );
        // ...but the tombstone is permanent on the fast path...
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence,
            "restoration must never un-corrupt a sealed tombstone"
        );
        // ...and re-sealing inherits it rather than forgiving it.
        let resealed = seal(&fx);
        assert_eq!(resealed.ordinal, sealed.ordinal + 1);
        assert!(
            resealed.corrupt_rows >= sealed.corrupt_rows,
            "re-checkpointing must preserve every prior tombstone"
        );
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    #[test]
    fn delta_rows_recorded_after_the_seal_still_classify() {
        let db = temp_db("delta");
        let fx = fixture(&db);
        seal(&fx);
        let (kernels2, recorded_at2) = record_one_run(&fx.ledger);
        assert_eq!(fx.kernels, kernels2);
        let (kernel, version) = &fx.kernels[0];

        // The second run's row is not covered by the checkpoint; the fast
        // path validates it exhaustively and classifies over BOTH rows.
        assert_eq!(
            fast(&fx, kernel, version, recorded_at2 + 1),
            exhaustive(&fx, kernel, version, recorded_at2 + 1)
        );
        assert_eq!(fast(&fx, kernel, version, recorded_at2 + 1), Staleness::Fresh);

        // Delta validation costs more than covered but the next seal
        // re-covers everything.
        let before = fx.ledger.read_queries();
        let _ = fast(&fx, kernel, version, recorded_at2 + 1);
        let delta_reads = fx.ledger.read_queries() - before;
        assert!(delta_reads > 2, "delta rows must be exhaustively validated");

        let receipt = seal(&fx);
        assert_eq!(receipt.ordinal, 1);
        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, recorded_at2 + 1), Staleness::Fresh);
        assert_eq!(
            fx.ledger.read_queries() - before,
            2,
            "after re-sealing, the fast path is fully covered again"
        );
        cleanup_db(&db);
    }

    #[test]
    #[ignore = "scale fixture: ~40 production runs; run explicitly"]
    fn fast_path_read_cost_stays_constant_as_history_grows() {
        let db = temp_db("scale");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let mut newest = fx.recorded_at;
        let mut exhaustive_read_counts = Vec::new();
        for _ in 0..39 {
            let (_, at) = record_one_run(&fx.ledger);
            newest = at;
        }
        seal(&fx);

        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, newest + 1), Staleness::Fresh);
        let fast_reads = fx.ledger.read_queries() - before;
        assert_eq!(fast_reads, 2, "40 covered rows must still cost two reads");

        let before = fx.ledger.read_queries();
        assert_eq!(exhaustive(&fx, kernel, version, newest + 1), Staleness::Fresh);
        exhaustive_read_counts.push(fx.ledger.read_queries() - before);
        assert!(
            exhaustive_read_counts[0] > 40,
            "exhaustive cost must scale with history ({} reads)",
            exhaustive_read_counts[0]
        );
        cleanup_db(&db);
    }
}
