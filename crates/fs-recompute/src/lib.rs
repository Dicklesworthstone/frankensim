//! fs-recompute — Proposal 2's STORE (bead lmp4.6). Layer: L6.
//!
//! A content-addressed Merkle DAG where every node records
//! `(op_id, input_hashes, params, code_version_hash, rng_seed,
//! achieved_error, required_tolerance)` and the gap
//! `required_tolerance − achieved_error` is the node's SLACK — the
//! resource incremental recompute spends. The Error Ledger becomes a
//! build graph with a SOUNDNESS CERTIFICATE for every skip:
//! [`Store::can_skip`] answers "is the cached artifact still good
//! enough?" with the slack attached, and a tolerance tightened past
//! the achieved error forces recomputation with the deficit named.
//!
//! DETERMINISM IS THE CERTIFIED CONTRACT here, not a nicety:
//! tolerance-level memoization requires bit-stable recomputation, so
//! [`Store::put`] TRIPS ([`StoreError::DeterminismViolation`]) when
//! the same node record arrives with different artifact bytes — the
//! write path itself polices the contract, and the conformance battery
//! certifies bit-identical artifacts across worker counts and
//! adversarial completion orders (risk R2, owned here).
//!
//! Pinning: nodes referenced by evidence packages or contracts are
//! NEVER evicted — the eviction pass can only touch unpinned nodes.

#[cfg(feature = "tolerance-invalidation")]
pub mod api;
#[cfg(feature = "tolerance-invalidation")]
pub mod invalidate;

use fs_ledger::{ContentHash, hash_bytes};
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A canonical parameter value (floats travel as bits).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ParamValue {
    /// A float, keyed by its bit pattern.
    F64(u64),
    /// An integer.
    Int(i64),
    /// A string.
    Str(String),
}

impl ParamValue {
    /// Convenience: from a float.
    #[must_use]
    pub fn f(v: f64) -> ParamValue {
        ParamValue::F64(v.to_bits())
    }
}

fn push_u64(buf: &mut Vec<u8>, value: u64) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn push_bytes(buf: &mut Vec<u8>, value: &[u8]) {
    push_u64(buf, value.len() as u64);
    buf.extend_from_slice(value);
}

fn push_string(buf: &mut Vec<u8>, value: &str) {
    push_bytes(buf, value.as_bytes());
}

pub(crate) fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{00}'..='\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(ch));
            }
            _ => out.push(ch),
        }
    }
    out.push('"');
    out
}

pub(crate) fn json_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.17e}")
    } else {
        "null".to_string()
    }
}

/// The seven-field node record (the Merkle DAG schema).
#[derive(Debug, Clone, PartialEq)]
pub struct NodeRecord {
    /// Operator identity.
    pub op_id: String,
    /// Content hashes of the inputs (edges of the DAG).
    pub input_hashes: Vec<ContentHash>,
    /// Canonical parameters (sorted by key at hash time).
    pub params: Vec<(String, ParamValue)>,
    /// The code version that computed it.
    pub code_version_hash: ContentHash,
    /// The seed (P2: seeds are data).
    pub rng_seed: u64,
    /// The error the computation ACHIEVED.
    pub achieved_error: f64,
    /// The tolerance the consumer REQUIRED.
    pub required_tolerance: f64,
}

impl NodeRecord {
    /// The node's SLACK: `required_tolerance − achieved_error`. May be
    /// NEGATIVE (an over-budget node) — representable on purpose, and
    /// a negative-slack node never satisfies a skip.
    #[must_use]
    pub fn slack(&self) -> f64 {
        self.required_tolerance - self.achieved_error
    }

    /// Stable content hash of the record (canonical serialization,
    /// floats as bits, params sorted by key).
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        // Versioned, length-prefixed binary encoding. Delimiter-based text
        // encoding is not injective when caller-controlled strings can contain
        // newlines or field-looking prefixes.
        let mut buf = Vec::new();
        buf.extend_from_slice(b"fs-recompute-node-v2\0");
        push_string(&mut buf, &self.op_id);
        push_u64(&mut buf, self.input_hashes.len() as u64);
        for h in &self.input_hashes {
            buf.extend_from_slice(h.as_bytes());
        }
        let mut params = self.params.clone();
        params.sort();
        push_u64(&mut buf, params.len() as u64);
        for (k, v) in &params {
            push_string(&mut buf, k);
            match v {
                ParamValue::F64(bits) => {
                    buf.push(0);
                    push_u64(&mut buf, *bits);
                }
                ParamValue::Int(i) => {
                    buf.push(1);
                    buf.extend_from_slice(&i.to_le_bytes());
                }
                ParamValue::Str(s) => {
                    buf.push(2);
                    push_string(&mut buf, s);
                }
            }
        }
        buf.extend_from_slice(self.code_version_hash.as_bytes());
        push_u64(&mut buf, self.rng_seed);
        push_u64(&mut buf, self.achieved_error.to_bits());
        push_u64(&mut buf, self.required_tolerance.to_bits());
        hash_bytes(&buf)
    }

    /// Canonical ledger row (node fields + slack).
    #[must_use]
    pub fn to_row(&self, artifact: &ContentHash) -> String {
        let input_hashes = self
            .input_hashes
            .iter()
            .map(|h| format!("\"{}\"", h.to_hex()))
            .collect::<Vec<_>>()
            .join(",");
        let mut params = self.params.clone();
        params.sort();
        let params = params
            .iter()
            .map(|(key, value)| match value {
                ParamValue::F64(bits) => format!(
                    "{{\"key\":{},\"kind\":\"f64\",\"bits\":\"{bits:016X}\"}}",
                    json_string(key)
                ),
                ParamValue::Int(value) => format!(
                    "{{\"key\":{},\"kind\":\"int\",\"value\":{value}}}",
                    json_string(key)
                ),
                ParamValue::Str(value) => format!(
                    "{{\"key\":{},\"kind\":\"string\",\"value\":{}}}",
                    json_string(key),
                    json_string(value)
                ),
            })
            .collect::<Vec<_>>()
            .join(",");
        let slack = self.slack();
        format!(
            "{{\"op\":{},\"node\":\"{}\",\"artifact\":\"{}\",\
             \"input_hashes\":[{input_hashes}],\"params\":[{params}],\
             \"code_version\":\"{}\",\"seed\":{},\"achieved\":{},\
             \"achieved_bits\":\"{:016X}\",\"required\":{},\
             \"required_bits\":\"{:016X}\",\"slack\":{},\
             \"slack_bits\":\"{:016X}\"}}",
            json_string(&self.op_id),
            self.content_hash().to_hex(),
            artifact.to_hex(),
            self.code_version_hash.to_hex(),
            self.rng_seed,
            json_f64(self.achieved_error),
            self.achieved_error.to_bits(),
            json_f64(self.required_tolerance),
            self.required_tolerance.to_bits(),
            json_f64(slack),
            slack.to_bits()
        )
    }
}

/// Why a node is pinned (never evicted).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum PinReason {
    /// Referenced by an evidence package (Proposal 12).
    EvidencePackage(String),
    /// Referenced by a contract (Proposal E).
    Contract(String),
}

/// A stored node. The RECORD is immutable (it IS the content
/// identity); absorbed perturbations accumulate in `burned`, mutable
/// runtime state that never touches the hash.
#[derive(Debug, Clone)]
pub struct StoredNode {
    /// The record (immutable identity).
    pub record: NodeRecord,
    /// Hash of the artifact bytes this record produced.
    pub artifact_hash: ContentHash,
    /// Pins (empty = evictable).
    pub pins: Vec<PinReason>,
    /// Insertion order (deterministic eviction).
    pub seq: u64,
    /// Slack burned by absorbed perturbations (runtime state).
    pub burned: f64,
}

impl StoredNode {
    /// Slack remaining after burns: `record.slack() − burned`.
    #[must_use]
    pub fn effective_slack(&self) -> f64 {
        self.record.slack() - self.burned
    }
}

/// Outcome of a put.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PutOutcome {
    /// New node stored.
    Inserted(ContentHash),
    /// Identical record + identical artifact already present (the
    /// memoization hit at write time).
    Deduped(ContentHash),
}

/// Skip-soundness decision.
#[derive(Debug, Clone, PartialEq)]
pub enum SkipDecision {
    /// The cached artifact satisfies the new tolerance: skipping is
    /// SOUND, with this much slack left.
    Hit {
        /// The cached node.
        node: ContentHash,
        /// `new_tolerance − achieved_error` (≥ 0).
        slack: f64,
    },
    /// The tolerance tightened past what the cached run achieved:
    /// recompute, and by this much.
    ToleranceTightened {
        /// `achieved_error − new_tolerance` (> 0).
        deficit: f64,
    },
    /// The requested tolerance was not a finite, non-negative error
    /// magnitude, so no skip certificate can be issued.
    InvalidTolerance {
        /// Exact bits supplied by the caller.
        tolerance_bits: u64,
    },
    /// No node with this identity exists.
    Miss,
}

/// Store errors (the determinism trip-wire lives here).
#[derive(Debug, Clone, PartialEq)]
pub enum StoreError {
    /// THE CONTRACT TRIP-WIRE: the same node record produced different
    /// artifact bytes — determinism is broken and memoization would be
    /// UNSOUND. This is a stop-the-line error, not a warning.
    DeterminismViolation {
        /// The node whose recomputation diverged.
        node: ContentHash,
        /// The artifact hash on record.
        expected: String,
        /// The artifact hash just produced.
        got: String,
    },
    /// Unknown node.
    UnknownNode {
        /// The hash asked for.
        node: ContentHash,
    },
    /// The cache's PINNED population alone exceeds the requested
    /// capacity — a structured refusal, never an OOM or a deadlock.
    CacheFullOfPins {
        /// How many nodes are pinned.
        pinned: usize,
        /// The capacity requested.
        capacity: usize,
    },
    /// Error bounds are magnitudes and must be finite and non-negative.
    InvalidErrorBudget {
        /// Supplied achieved-error bits.
        achieved_bits: u64,
        /// Supplied required-tolerance bits.
        required_bits: u64,
    },
    /// A slack burn was malformed or did not fit strictly inside the
    /// currently available slack.
    InvalidSlackBurn {
        /// The node whose slack would be changed.
        node: ContentHash,
        /// Requested burn bits.
        amount_bits: u64,
        /// Available slack bits.
        available_bits: u64,
    },
    /// A pure plan was committed after the store state it certified had
    /// changed.
    StalePlan {
        /// Store revision captured by the plan.
        planned_revision: u64,
        /// Store revision at commit time.
        current_revision: u64,
        /// State fingerprint captured by the plan.
        planned_state: ContentHash,
        /// State fingerprint at commit time.
        current_state: ContentHash,
    },
}

impl core::fmt::Display for StoreError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            StoreError::DeterminismViolation {
                node,
                expected,
                got,
            } => write!(
                f,
                "DETERMINISM CONTRACT VIOLATION at node {}: the same \
                 (op, inputs, params, code, seed) produced artifact {got} where \
                 {expected} is on record — tolerance-level memoization is UNSOUND \
                 until the op is fixed (unordered reduction? unstable sort? \
                 uninitialized padding?); this is stop-the-line, not a warning",
                node.to_hex()
            ),
            StoreError::UnknownNode { node } => {
                write!(f, "node {} is not in the store", node.to_hex())
            }
            StoreError::CacheFullOfPins { pinned, capacity } => write!(
                f,
                "{pinned} pinned nodes exceed the requested capacity {capacity};                  pins are re-verifiability PROMISES (evidence packages, contracts)                  and cannot be evicted — raise the capacity or retire the promises"
            ),
            StoreError::InvalidErrorBudget {
                achieved_bits,
                required_bits,
            } => write!(
                f,
                "invalid error budget: achieved={} required={}; both must be finite, non-negative magnitudes",
                f64::from_bits(*achieved_bits),
                f64::from_bits(*required_bits)
            ),
            StoreError::InvalidSlackBurn {
                node,
                amount_bits,
                available_bits,
            } => write!(
                f,
                "invalid slack burn at node {}: amount={} must be finite, non-negative, and strictly below available slack {}",
                node.to_hex(),
                f64::from_bits(*amount_bits),
                f64::from_bits(*available_bits)
            ),
            StoreError::StalePlan {
                planned_revision,
                current_revision,
                planned_state,
                current_state,
            } => write!(
                f,
                "stale recompute plan: certified store revision {planned_revision} state {}, current revision {current_revision} state {}; re-plan before committing",
                planned_state.to_hex(),
                current_state.to_hex()
            ),
        }
    }
}

impl std::error::Error for StoreError {}

/// The content-addressed store.
#[derive(Debug, Default)]
pub struct Store {
    nodes: BTreeMap<[u8; 32], StoredNode>,
    seq: u64,
    revision: u64,
    rows: Vec<String>,
}

fn key(h: &ContentHash) -> [u8; 32] {
    *h.as_bytes()
}

impl Store {
    /// Empty store.
    #[must_use]
    pub fn new() -> Self {
        Store::default()
    }

    /// Number of stored nodes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// True when empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// The canonical ledger rows written so far.
    #[must_use]
    pub fn rows(&self) -> &[String] {
        &self.rows
    }

    /// Monotonic mutation revision used to bind pure plans to the state
    /// against which their certificates were computed.
    #[must_use]
    pub fn revision(&self) -> u64 {
        self.revision
    }

    fn bump_revision(&mut self) {
        self.revision = self.revision.saturating_add(1);
    }

    /// Deterministic fingerprint of every field that can affect an
    /// invalidation certificate. This catches cross-store plans and remains
    /// authoritative if the diagnostic revision counter saturates.
    pub(crate) fn state_fingerprint(&self) -> ContentHash {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"fs-recompute-store-state-v1\0");
        push_u64(&mut buf, self.nodes.len() as u64);
        for (node_key, node) in &self.nodes {
            buf.extend_from_slice(node_key);
            buf.extend_from_slice(node.artifact_hash.as_bytes());
            push_u64(&mut buf, node.burned.to_bits());
        }
        hash_bytes(&buf)
    }

    /// Store a computed node. Re-putting the identical record with the
    /// identical artifact dedupes; the identical record with DIFFERENT
    /// artifact bytes trips the determinism contract.
    ///
    /// # Errors
    /// [`StoreError::DeterminismViolation`] — stop the line.
    pub fn put(
        &mut self,
        record: NodeRecord,
        artifact_bytes: &[u8],
    ) -> Result<PutOutcome, StoreError> {
        if !record.achieved_error.is_finite()
            || record.achieved_error < 0.0
            || !record.required_tolerance.is_finite()
            || record.required_tolerance < 0.0
        {
            return Err(StoreError::InvalidErrorBudget {
                achieved_bits: record.achieved_error.to_bits(),
                required_bits: record.required_tolerance.to_bits(),
            });
        }
        let node_hash = record.content_hash();
        let artifact_hash = hash_bytes(artifact_bytes);
        if let Some(existing) = self.nodes.get(&key(&node_hash)) {
            if existing.artifact_hash == artifact_hash {
                return Ok(PutOutcome::Deduped(node_hash));
            }
            return Err(StoreError::DeterminismViolation {
                node: node_hash,
                expected: existing.artifact_hash.to_hex(),
                got: artifact_hash.to_hex(),
            });
        }
        self.rows.push(record.to_row(&artifact_hash));
        self.nodes.insert(
            key(&node_hash),
            StoredNode {
                record,
                artifact_hash,
                pins: Vec::new(),
                seq: self.seq,
                burned: 0.0,
            },
        );
        self.seq += 1;
        self.bump_revision();
        Ok(PutOutcome::Inserted(node_hash))
    }

    /// The stored node for a record identity, if any.
    #[must_use]
    pub fn lookup(&self, record: &NodeRecord) -> Option<&StoredNode> {
        self.nodes.get(&key(&record.content_hash()))
    }

    /// The stored node by hash.
    #[must_use]
    pub fn get(&self, node: &ContentHash) -> Option<&StoredNode> {
        self.nodes.get(&key(node))
    }

    /// Skip soundness: is the cached artifact for this identity (op,
    /// inputs, params, code, seed — tolerance EXCLUDED from identity
    /// here) still good enough for `new_tolerance`? The certificate is
    /// the returned slack.
    #[must_use]
    pub fn can_skip(&self, record: &NodeRecord, new_tolerance: f64) -> SkipDecision {
        if !new_tolerance.is_finite() || new_tolerance < 0.0 {
            return SkipDecision::InvalidTolerance {
                tolerance_bits: new_tolerance.to_bits(),
            };
        }
        // Identity for skip purposes: the record with its tolerance
        // fields normalized out.
        let mut probe = record.clone();
        probe.achieved_error = 0.0;
        probe.required_tolerance = 0.0;
        let probe_hash = probe.content_hash();
        // Scan for a node with the same normalized identity.
        let mut best: Option<(&StoredNode, f64)> = None;
        for stored in self.nodes.values() {
            let mut norm = stored.record.clone();
            norm.achieved_error = 0.0;
            norm.required_tolerance = 0.0;
            if norm.content_hash() == probe_hash {
                let effective_error = stored.record.achieved_error + stored.burned;
                if best.is_none_or(|(_, current)| effective_error < current) {
                    best = Some((stored, effective_error));
                }
            }
        }
        let Some((stored, effective_error)) = best else {
            return SkipDecision::Miss;
        };
        let slack = new_tolerance - effective_error;
        if slack >= 0.0 {
            SkipDecision::Hit {
                node: stored.record.content_hash(),
                slack,
            }
        } else {
            SkipDecision::ToleranceTightened { deficit: -slack }
        }
    }

    /// Pin a node (evidence package / contract reference): pinned
    /// nodes are NEVER evicted.
    ///
    /// # Errors
    /// [`StoreError::UnknownNode`].
    pub fn pin(&mut self, node: &ContentHash, reason: PinReason) -> Result<(), StoreError> {
        let entry = self
            .nodes
            .get_mut(&key(node))
            .ok_or(StoreError::UnknownNode { node: *node })?;
        if !entry.pins.contains(&reason) {
            entry.pins.push(reason);
            entry.pins.sort();
            self.bump_revision();
        }
        Ok(())
    }

    /// Evict unpinned nodes (oldest first, deterministic) until at
    /// most `keep` UNPINNED nodes remain. Returns how many were
    /// evicted. Pinned nodes are untouchable by construction.
    pub fn evict_unpinned(&mut self, keep: usize) -> u32 {
        let mut unpinned: Vec<([u8; 32], u64)> = self
            .nodes
            .iter()
            .filter(|(_, n)| n.pins.is_empty())
            .map(|(k, n)| (*k, n.seq))
            .collect();
        unpinned.sort_by_key(|&(_, seq)| seq);
        let excess = unpinned.len().saturating_sub(keep);
        let mut evicted = 0;
        for &(k, _) in unpinned.iter().take(excess) {
            self.nodes.remove(&k);
            evicted += 1;
        }
        if evicted > 0 {
            self.bump_revision();
        }
        evicted
    }

    /// Iterate stored nodes (BTree key order; deterministic).
    pub fn iter(&self) -> impl Iterator<Item = ([u8; 32], &StoredNode)> {
        self.nodes.iter().map(|(k, v)| (*k, v))
    }

    /// Burn absorbed perturbation into a node's achieved error (the
    /// slack is a SPENDABLE resource: repeat perturbations see the
    /// reduced remainder).
    ///
    /// # Errors
    /// [`StoreError::UnknownNode`].
    pub fn burn_slack(&mut self, node: &ContentHash, amount: f64) -> Result<(), StoreError> {
        self.commit_burns(self.revision, self.state_fingerprint(), &[(*node, amount)])
    }

    pub(crate) fn commit_burns(
        &mut self,
        planned_revision: u64,
        planned_state: ContentHash,
        burns: &[(ContentHash, f64)],
    ) -> Result<(), StoreError> {
        let current_state = self.state_fingerprint();
        if self.revision != planned_revision || current_state != planned_state {
            return Err(StoreError::StalePlan {
                planned_revision,
                current_revision: self.revision,
                planned_state,
                current_state,
            });
        }
        let mut aggregated = BTreeMap::<[u8; 32], (ContentHash, f64)>::new();
        for (node, amount) in burns {
            let entry = self
                .nodes
                .get(&key(node))
                .ok_or(StoreError::UnknownNode { node: *node })?;
            let available = entry.effective_slack();
            if !amount.is_finite() || *amount < 0.0 {
                return Err(StoreError::InvalidSlackBurn {
                    node: *node,
                    amount_bits: amount.to_bits(),
                    available_bits: available.to_bits(),
                });
            }
            let total = &mut aggregated.entry(key(node)).or_insert((*node, 0.0)).1;
            *total += *amount;
            if !total.is_finite() {
                return Err(StoreError::InvalidSlackBurn {
                    node: *node,
                    amount_bits: total.to_bits(),
                    available_bits: available.to_bits(),
                });
            }
        }
        for (node_key, (node, amount)) in &aggregated {
            let available = self.nodes[node_key].effective_slack();
            if *amount >= available {
                return Err(StoreError::InvalidSlackBurn {
                    node: *node,
                    amount_bits: amount.to_bits(),
                    available_bits: available.to_bits(),
                });
            }
        }
        for (node_key, (_, amount)) in aggregated {
            self.nodes
                .get_mut(&node_key)
                .expect("burns prevalidated")
                .burned += amount;
        }
        self.bump_revision();
        Ok(())
    }

    /// Remove a node by raw key (the eviction path).
    #[allow(dead_code)] // wired by the eviction path as it lands; keeping the seam
    pub(crate) fn remove_by_key(&mut self, k: [u8; 32]) {
        if self.nodes.remove(&k).is_some() {
            self.bump_revision();
        }
    }

    /// Serialize the store to its canonical text form (round-trips;
    /// "hash stability under fork").
    #[must_use]
    pub fn snapshot(&self) -> String {
        let mut out = String::from("fsrecompute v2\n");
        for node in self.nodes.values() {
            let _ = writeln!(out, "{}", node.record.to_row(&node.artifact_hash));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(required_tolerance: f64) -> NodeRecord {
        NodeRecord {
            op_id: "aggregate-burn".to_string(),
            input_hashes: Vec::new(),
            params: Vec::new(),
            code_version_hash: hash_bytes(b"test-code"),
            rng_seed: 1,
            achieved_error: 0.0,
            required_tolerance,
        }
    }

    #[test]
    fn duplicate_burns_are_aggregated_before_mutation() {
        let mut store = Store::new();
        let PutOutcome::Inserted(node) = store.put(record(1.0), b"artifact").expect("put") else {
            unreachable!("fresh store");
        };
        let revision = store.revision();
        let state = store.state_fingerprint();
        let refused = store.commit_burns(revision, state, &[(node, 0.6), (node, 0.6)]);
        assert!(matches!(refused, Err(StoreError::InvalidSlackBurn { .. })));
        assert_eq!(
            store.get(&node).expect("node").effective_slack().to_bits(),
            1.0f64.to_bits(),
            "aggregate overflow must refuse before any partial burn"
        );

        store.revision = u64::MAX;
        store.bump_revision();
        assert_eq!(store.revision(), u64::MAX, "revision must never wrap");
    }
}
