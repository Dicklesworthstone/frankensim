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
        let mut buf = String::new();
        let _ = writeln!(buf, "op:{}", self.op_id);
        for h in &self.input_hashes {
            let _ = writeln!(buf, "in:{}", h.to_hex());
        }
        let mut params = self.params.clone();
        params.sort();
        for (k, v) in &params {
            match v {
                ParamValue::F64(bits) => {
                    let _ = writeln!(buf, "pf:{k}={bits:016X}");
                }
                ParamValue::Int(i) => {
                    let _ = writeln!(buf, "pi:{k}={i}");
                }
                ParamValue::Str(s) => {
                    let _ = writeln!(buf, "ps:{k}={s}");
                }
            }
        }
        let _ = writeln!(buf, "code:{}", self.code_version_hash.to_hex());
        let _ = writeln!(buf, "seed:{}", self.rng_seed);
        let _ = writeln!(buf, "ach:{:016X}", self.achieved_error.to_bits());
        let _ = writeln!(buf, "req:{:016X}", self.required_tolerance.to_bits());
        hash_bytes(buf.as_bytes())
    }

    /// Canonical ledger row (node fields + slack).
    #[must_use]
    pub fn to_row(&self, artifact: &ContentHash) -> String {
        format!(
            "{{\"op\":\"{}\",\"node\":\"{}\",\"artifact\":\"{}\",\"inputs\":{},\
             \"seed\":{},\"achieved\":{:.6e},\"required\":{:.6e},\"slack\":{:.6e}}}",
            self.op_id,
            self.content_hash().to_hex(),
            artifact.to_hex(),
            self.input_hashes.len(),
            self.rng_seed,
            self.achieved_error,
            self.required_tolerance,
            self.slack()
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
        }
    }
}

impl std::error::Error for StoreError {}

/// The content-addressed store.
#[derive(Debug, Default)]
pub struct Store {
    nodes: BTreeMap<[u8; 32], StoredNode>,
    seq: u64,
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
        // Identity for skip purposes: the record with its tolerance
        // fields normalized out.
        let mut probe = record.clone();
        probe.achieved_error = 0.0;
        probe.required_tolerance = 0.0;
        let probe_hash = probe.content_hash();
        // Scan for a node with the same normalized identity.
        for stored in self.nodes.values() {
            let mut norm = stored.record.clone();
            norm.achieved_error = 0.0;
            norm.required_tolerance = 0.0;
            if norm.content_hash() == probe_hash {
                let slack = new_tolerance - (stored.record.achieved_error + stored.burned);
                if slack >= 0.0 {
                    return SkipDecision::Hit {
                        node: stored.record.content_hash(),
                        slack,
                    };
                }
                return SkipDecision::ToleranceTightened { deficit: -slack };
            }
        }
        SkipDecision::Miss
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
        let entry = self
            .nodes
            .get_mut(&key(node))
            .ok_or(StoreError::UnknownNode { node: *node })?;
        entry.burned += amount;
        Ok(())
    }

    /// Remove a node by raw key (the eviction path).
    #[allow(dead_code)] // wired by the eviction path as it lands; keeping the seam
    pub(crate) fn remove_by_key(&mut self, k: [u8; 32]) {
        self.nodes.remove(&k);
    }

    /// Serialize the store to its canonical text form (round-trips;
    /// "hash stability under fork").
    #[must_use]
    pub fn snapshot(&self) -> String {
        let mut out = String::from("fsrecompute v1\n");
        for node in self.nodes.values() {
            let _ = writeln!(out, "{}", node.record.to_row(&node.artifact_hash));
        }
        out
    }
}
