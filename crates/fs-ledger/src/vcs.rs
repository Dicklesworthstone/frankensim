//! VERSION CONTROL FOR PHYSICS (addendum Proposal 10, the base verbs):
//! COMMITS, BRANCHES, CHECKOUT over Merkle roots — free-riding on the
//! forkable-worlds machinery (`travel`) that already implements
//! hash-shared branches. PLM systems version files; this versions the
//! design-plus-ledger STATE: a commit is the Merkle root of a branch's
//! visible frozen ops (Five Explicits + outcome + linked artifact
//! hashes; wall-clock times and rowids EXCLUDED so logically identical
//! histories produce identical roots). Branches are pointers; a checkout
//! between nearby commits costs the SYMMETRIC-DIFFERENCE frontier — the
//! `perturb()`-style delta the incremental-recompute store consumes —
//! not a full re-solve. Garbage collection respects every live branch by
//! construction (content addressing + the travel GC's reachability walk).

use crate::hash::{Blake3, ContentHash, hash_bytes};
use crate::{EventRow, Ledger, LedgerError};
use fsqlite::SqliteValue;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// One recorded commit.
#[derive(Debug, Clone, PartialEq)]
pub struct CommitInfo {
    /// The Merkle root (the commit's identity).
    pub root: ContentHash,
    /// The branch it snapshots.
    pub branch: i64,
    /// The frontier op id captured (None for an empty history).
    pub frontier_op: Option<i64>,
    /// The parent commit on the same branch, if any.
    pub parent: Option<ContentHash>,
}

/// The checkout delta between two commits: the ops a delta-solver must
/// reconcile (everything else is hash-shared and untouched).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckoutDelta {
    /// Op ids visible in the target but not the source.
    pub added: Vec<i64>,
    /// Op ids visible in the source but not the target.
    pub removed: Vec<i64>,
    /// Op ids shared by both views (the hash-shared bulk).
    pub shared: usize,
}

/// Merge-view bookkeeping for the diff/bisect/merge consumers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergeViews {
    /// The merge base: ops visible from BOTH branches.
    pub base: Vec<i64>,
    /// Ops only on branch A.
    pub only_a: Vec<i64>,
    /// Ops only on branch B.
    pub only_b: Vec<i64>,
}

/// Artifact-sharing audit: the storage story for N branches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StorageAudit {
    /// Physical artifact rows in the store.
    pub physical_artifacts: u64,
    /// Sum over branches of visible artifact references.
    pub logical_references: u64,
    /// Live branch count.
    pub branches: u64,
}

/// The in-session commit registry (commits are also persisted as
/// `vcs-commit` events; the registry is the fast lookup surface).
#[derive(Debug, Default)]
pub struct Vcs {
    commits: BTreeMap<[u8; 32], CommitInfo>,
    heads: BTreeMap<i64, ContentHash>,
}

impl Ledger {
    /// Sorted artifact hashes linked to an op (lineage edges, both
    /// roles) — the artifact content folded into commit leaves.
    ///
    /// # Errors
    /// Engine errors.
    pub fn op_artifact_hashes(&self, op: i64) -> Result<Vec<ContentHash>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT artifact FROM edges WHERE op = ?1 ORDER BY artifact",
                &[SqliteValue::Integer(op)],
            )
            .map_err(|e| LedgerError::Sql {
                context: "op_artifact_hashes".to_string(),
                detail: e.to_string(),
            })?;
        let mut out = Vec::with_capacity(rows.len());
        for row in &rows {
            match row.get(0) {
                Some(SqliteValue::Blob(b)) if b.len() == 32 => {
                    let mut h = [0u8; 32];
                    h.copy_from_slice(b);
                    out.push(ContentHash(h));
                }
                other => {
                    return Err(LedgerError::Sql {
                        context: "op_artifact_hashes".to_string(),
                        detail: format!("artifact: expected 32-byte BLOB, got {other:?}"),
                    });
                }
            }
        }
        Ok(out)
    }

    /// The commit LEAF hash of one op: canonical frozen content (Five
    /// Explicits, outcome, diag, exec mode via the op row) + sorted
    /// linked-artifact hashes. Wall times and rowids are EXCLUDED —
    /// logically identical histories hash identically.
    ///
    /// # Errors
    /// Engine errors; unknown op ids.
    pub fn commit_leaf(&self, op: i64) -> Result<ContentHash, LedgerError> {
        let row = self.op(op)?.ok_or_else(|| LedgerError::Invalid {
            field: "op".to_string(),
            problem: format!("op {op} does not exist"),
        })?;
        let mut canon = String::new();
        let _ = write!(
            canon,
            "op;ir={};versions={};budget={};capability={};outcome={};diag={}",
            row.ir,
            row.versions,
            row.budget,
            row.capability,
            row.outcome.as_deref().unwrap_or("<in-flight>"),
            row.diag.as_deref().unwrap_or("")
        );
        let mut hasher = Blake3::new();
        hasher.update(canon.as_bytes());
        hasher.update(&row.seed);
        for a in self.op_artifact_hashes(op)? {
            hasher.update(a.as_bytes());
        }
        Ok(hasher.finalize())
    }
}

/// Binary Merkle fold over leaves (odd tails promote; empty list hashes
/// a fixed domain tag so the empty commit is well-defined).
fn merkle_root(mut level: Vec<ContentHash>) -> ContentHash {
    if level.is_empty() {
        return hash_bytes(b"vcs-empty-commit");
    }
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        for pair in level.chunks(2) {
            if pair.len() == 2 {
                let mut h = Blake3::new();
                h.update(pair[0].as_bytes());
                h.update(pair[1].as_bytes());
                next.push(h.finalize());
            } else {
                next.push(pair[0]);
            }
        }
        level = next;
    }
    level[0]
}

impl Vcs {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Self {
        Vcs::default()
    }

    /// COMMIT: snapshot a branch's visible state as a Merkle root,
    /// chain it to the branch's previous commit, and persist a
    /// `vcs-commit` event. Reproducible: the same visible state yields
    /// the same root, on any run.
    ///
    /// # Errors
    /// Unknown branches and engine errors.
    pub fn commit(&mut self, ledger: &Ledger, branch: i64) -> Result<CommitInfo, LedgerError> {
        if ledger.branch(branch)?.is_none() {
            return Err(LedgerError::Invalid {
                field: "branch".to_string(),
                problem: format!("branch {branch} does not exist"),
            });
        }
        let ops = ledger.visible_op_ids(branch, None)?;
        let mut leaves = Vec::with_capacity(ops.len());
        for op in &ops {
            leaves.push(ledger.commit_leaf(*op)?);
        }
        let root = merkle_root(leaves);
        let info = CommitInfo {
            root,
            branch,
            frontier_op: ops.last().copied(),
            parent: self.heads.get(&branch).copied(),
        };
        let payload = format!(
            "{{\"kind\":\"vcs-commit\",\"root\":\"{}\",\"branch\":{},\"frontier\":{},\
             \"parent\":{}}}",
            root.to_hex(),
            branch,
            info.frontier_op
                .map_or("null".to_string(), |o| o.to_string()),
            info.parent
                .map_or("null".to_string(), |p| format!("\"{}\"", p.to_hex())),
        );
        ledger.append_event(&EventRow {
            session: None,
            t: info.frontier_op.unwrap_or(0),
            kind: "vcs-commit",
            payload: Some(&payload),
        })?;
        self.commits.insert(root.0, info.clone());
        self.heads.insert(branch, root);
        Ok(info)
    }

    /// A commit by root.
    #[must_use]
    pub fn lookup(&self, root: &ContentHash) -> Option<&CommitInfo> {
        self.commits.get(&root.0)
    }

    /// The current head of a branch.
    #[must_use]
    pub fn head(&self, branch: i64) -> Option<ContentHash> {
        self.heads.get(&branch).copied()
    }

    /// CHECKOUT: materialize a commit's view (ops + finished artifacts).
    ///
    /// # Errors
    /// A structured error for an unknown root (the boundary case);
    /// engine errors.
    pub fn checkout(
        &self,
        ledger: &Ledger,
        root: &ContentHash,
    ) -> Result<crate::travel::ViewSnapshot, LedgerError> {
        let info = self.lookup(root).ok_or_else(|| LedgerError::Invalid {
            field: "commit".to_string(),
            problem: format!("no commit with root {} is known", root.to_hex()),
        })?;
        ledger.at_time(info.branch, i64::MAX).map(|mut snap| {
            // Trim to the committed frontier (later ops are not part of
            // the commit).
            if let Some(frontier) = info.frontier_op {
                snap.ops.retain(|o| o.id <= frontier);
            } else {
                snap.ops.clear();
            }
            snap
        })
    }

    /// The CHECKOUT DELTA between two commits: the symmetric difference
    /// of visible op sets — the frontier a delta-solver reconciles. A
    /// nearby checkout is CHEAP because `shared` dominates.
    ///
    /// # Errors
    /// Unknown roots; engine errors.
    pub fn checkout_delta(
        &self,
        ledger: &Ledger,
        from: &ContentHash,
        to: &ContentHash,
    ) -> Result<CheckoutDelta, LedgerError> {
        let visible = |root: &ContentHash| -> Result<BTreeSet<i64>, LedgerError> {
            let info = self.lookup(root).ok_or_else(|| LedgerError::Invalid {
                field: "commit".to_string(),
                problem: format!("no commit with root {} is known", root.to_hex()),
            })?;
            Ok(ledger
                .visible_op_ids(info.branch, info.frontier_op)?
                .into_iter()
                .collect())
        };
        let a = visible(from)?;
        let b = visible(to)?;
        Ok(CheckoutDelta {
            added: b.difference(&a).copied().collect(),
            removed: a.difference(&b).copied().collect(),
            shared: a.intersection(&b).count(),
        })
    }

    /// The divergence point of two branches: the deepest common branch
    /// segment's cap (the shared history boundary the merge base builds
    /// on).
    ///
    /// # Errors
    /// Unknown branches; engine errors.
    pub fn merge_views(
        &self,
        ledger: &Ledger,
        branch_a: i64,
        branch_b: i64,
    ) -> Result<MergeViews, LedgerError> {
        let a: BTreeSet<i64> = ledger.visible_op_ids(branch_a, None)?.into_iter().collect();
        let b: BTreeSet<i64> = ledger.visible_op_ids(branch_b, None)?.into_iter().collect();
        Ok(MergeViews {
            base: a.intersection(&b).copied().collect(),
            only_a: a.difference(&b).copied().collect(),
            only_b: b.difference(&a).copied().collect(),
        })
    }

    /// STORAGE AUDIT: physical artifact rows vs logical per-branch
    /// references — the "N branches ≈ 1× + deltas" claim, measured.
    ///
    /// # Errors
    /// Engine errors.
    pub fn storage_audit(&self, ledger: &Ledger) -> Result<StorageAudit, LedgerError> {
        let physical = ledger.table_count("artifacts")?;
        let branches = ledger.branches()?;
        let mut logical = 0u64;
        for b in &branches {
            let snap = ledger.at_time(b.id, i64::MAX)?;
            logical += snap.artifacts.len() as u64;
        }
        Ok(StorageAudit {
            physical_artifacts: physical,
            logical_references: logical,
            branches: branches.len() as u64,
        })
    }
}
