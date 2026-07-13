//! VERSION CONTROL FOR PHYSICS (addendum Proposal 10, the base verbs):
//! COMMITS, BRANCHES, CHECKOUT over Merkle roots — free-riding on the
//! forkable-worlds machinery (`travel`) that already implements
//! hash-shared branches. PLM systems version files; this versions the
//! design-plus-ledger STATE: a commit is the Merkle root of a branch's
//! visible frozen ops (Five Explicits + outcome + execution mode +
//! role-qualified linked artifact hashes; wall-clock times and rowids
//! EXCLUDED so logically identical histories produce identical roots).
//! Branches are pointers; a checkout
//! between nearby commits costs the SYMMETRIC-DIFFERENCE frontier — the
//! `perturb()`-style delta the incremental-recompute store consumes —
//! not a full re-solve. Garbage collection respects every live branch by
//! construction (content addressing + the travel GC's reachability walk).

use crate::hash::{Blake3, ContentHash};
use crate::travel::{ExecMode, ViewSnapshot};
use crate::{EventRow, Ledger, LedgerError, VCS_IDENTITY_EVENT_KIND};
use fsqlite::SqliteValue;
use std::collections::{BTreeMap, BTreeSet};

const COMMIT_LEAF_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.commit-leaf.v2";
const MERKLE_PAIR_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.merkle-pair.v2";
const MERKLE_ODD_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.merkle-odd.v2";
const COMMIT_ROOT_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.commit-root.v2";
const LEDGER_IDENTITY_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.ledger-identity.v1";

fn hash_frame(hasher: &mut Blake3, bytes: &[u8]) {
    let len = u64::try_from(bytes.len()).expect("frame length fits in u64");
    hasher.update(&len.to_le_bytes());
    hasher.update(bytes);
}

fn domain_hasher(domain: &[u8]) -> Blake3 {
    let mut hasher = Blake3::new();
    hash_frame(&mut hasher, b"domain");
    hash_frame(&mut hasher, domain);
    hasher
}

fn hash_field(hasher: &mut Blake3, name: &[u8], value: &[u8]) {
    hash_frame(hasher, name);
    hash_frame(hasher, value);
}

fn hash_optional_field(hasher: &mut Blake3, name: &[u8], value: Option<&[u8]>) {
    hash_frame(hasher, name);
    match value {
        Some(value) => {
            hash_frame(hasher, b"present");
            hash_frame(hasher, value);
        }
        None => hash_frame(hasher, b"absent"),
    }
}

fn framed_hash(domain: &[u8], fields: &[(&[u8], &[u8])]) -> ContentHash {
    let mut hasher = domain_hasher(domain);
    for (name, value) in fields {
        hash_field(&mut hasher, name, value);
    }
    hasher.finalize()
}

/// The full ENVELOPE identity of a commit (bead gp3.17): WHICH ledger,
/// WHICH branch, WHICH semantic root. The Merkle root alone is the
/// SEMANTIC-STATE identity — two branches (or two ledgers) can reach
/// the same semantic state, and their commit envelopes must not
/// overwrite each other.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CommitId {
    /// The persisted ledger identity ([`Ledger::vcs_identity`]).
    pub ledger: ContentHash,
    /// The branch within that ledger.
    pub branch: i64,
    /// The semantic Merkle root.
    pub root: ContentHash,
}

/// One recorded commit.
#[derive(Debug, Clone, PartialEq)]
pub struct CommitInfo {
    /// The persisted identity of the ledger this commit came from.
    pub ledger: ContentHash,
    /// The Merkle root (the SEMANTIC identity; the envelope identity
    /// is [`CommitInfo::id`]).
    pub root: ContentHash,
    /// The branch it snapshots.
    pub branch: i64,
    /// The frontier op id captured (None for an empty history);
    /// LEDGER-LOCAL, never compared across ledgers.
    pub frontier_op: Option<i64>,
    /// The parent commit root on the same branch, if any.
    pub parent: Option<ContentHash>,
}

impl CommitInfo {
    /// The full envelope identity.
    #[must_use]
    pub fn id(&self) -> CommitId {
        CommitId {
            ledger: self.ledger,
            branch: self.branch,
            root: self.root,
        }
    }
}

/// One op in a checkout delta: identified SEMANTICALLY by its commit
/// leaf hash (portable across ledgers), with the op id THAT SIDE's
/// ledger uses locally (added ops carry target-local ids, removed ops
/// source-local ids).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeltaOp {
    /// The stable semantic leaf identity ([`Ledger::commit_leaf`]).
    pub leaf: ContentHash,
    /// The op id in the ledger of the side this entry came from.
    pub local_op: i64,
}

/// The checkout delta between two commits: the ops a delta-solver must
/// reconcile (everything else is hash-shared and untouched). Computed
/// from SEMANTIC leaf identities (bead gp3.17), so commits from
/// different ledger instances — or histories imported in a different
/// row order — compare correctly; local row ids are never compared
/// across ledgers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckoutDelta {
    /// Ops visible in the target but not the source (target-local ids).
    pub added: Vec<DeltaOp>,
    /// Ops visible in the source but not the target (source-local ids).
    pub removed: Vec<DeltaOp>,
    /// Leaf instances shared by both views (the hash-shared bulk).
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

type EnvelopeKey = (ContentHash, i64, [u8; 32]); // (ledger, branch, root)

#[derive(Debug, Clone)]
struct CommitRecord {
    info: CommitInfo,
    /// Semantic leaf hashes aligned with `snapshot.ops` (the portable
    /// identities cross-ledger deltas compare).
    leaves: Vec<ContentHash>,
    snapshot: ViewSnapshot,
}

/// The in-session commit registry (commits are also persisted as
/// `vcs-commit` events; the registry is the fast lookup surface).
/// Keyed by the FULL envelope identity (ledger, branch, root) — equal
/// semantic roots from different branches or ledgers coexist (bead
/// gp3.17).
#[derive(Debug, Default)]
pub struct Vcs {
    commits: BTreeMap<EnvelopeKey, CommitRecord>,
    heads: BTreeMap<(ContentHash, i64), CommitInfo>,
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

    fn op_artifact_edges(&self, op: i64) -> Result<Vec<(String, ContentHash)>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT role, artifact FROM edges WHERE op = ?1 ORDER BY role, artifact",
                &[SqliteValue::Integer(op)],
            )
            .map_err(|e| LedgerError::Sql {
                context: "op_artifact_edges".to_string(),
                detail: e.to_string(),
            })?;
        let mut edges = Vec::with_capacity(rows.len());
        for row in &rows {
            let role = match row.get(0) {
                Some(SqliteValue::Text(role)) if matches!(role.as_str(), "in" | "out") => {
                    role.as_str().to_string()
                }
                other => {
                    return Err(LedgerError::Corrupt {
                        hash_hex: String::new(),
                        detail: format!("op {op}: invalid edge role {other:?}"),
                    });
                }
            };
            let hash = match row.get(1) {
                Some(SqliteValue::Blob(bytes)) => {
                    ContentHash::from_slice(bytes).ok_or_else(|| LedgerError::Corrupt {
                        hash_hex: String::new(),
                        detail: format!("op {op}: edge contains a malformed artifact hash"),
                    })?
                }
                other => {
                    return Err(LedgerError::Sql {
                        context: "op_artifact_edges".to_string(),
                        detail: format!("artifact: expected 32-byte BLOB, got {other:?}"),
                    });
                }
            };
            edges.push((role, hash));
        }
        Ok(edges)
    }

    fn commit_exec_mode(&self, op: i64) -> Result<ExecMode, LedgerError> {
        let mode = self.bounded_op_exec_mode(op)?;
        ExecMode::parse(&mode).ok_or_else(|| LedgerError::OpCorrupt {
            op,
            detail: "execution mode passed storage guard but not enum parsing".to_string(),
        })
    }

    /// The commit LEAF hash of one op: canonical frozen content (Five
    /// Explicits, outcome, diagnostic, execution mode) + sorted,
    /// role-qualified linked-artifact hashes. The encoding is domain
    /// separated and length framed. Wall times, rowids, branch ids, and
    /// session envelopes are EXCLUDED so logically identical histories
    /// hash identically.
    ///
    /// # Errors
    /// Engine errors; unknown op ids.
    pub fn commit_leaf(&self, op: i64) -> Result<ContentHash, LedgerError> {
        let row = self.op(op)?.ok_or_else(|| LedgerError::Invalid {
            field: "op".to_string(),
            problem: format!("op {op} does not exist"),
        })?;
        let exec_mode = self.commit_exec_mode(op)?;
        let edges = self.op_artifact_edges(op)?;
        let mut hasher = domain_hasher(COMMIT_LEAF_DOMAIN);
        hash_field(&mut hasher, b"ir", row.ir.as_bytes());
        hash_field(&mut hasher, b"seed", &row.seed);
        hash_field(&mut hasher, b"versions", row.versions.as_bytes());
        hash_field(&mut hasher, b"budget", row.budget.as_bytes());
        hash_field(&mut hasher, b"capability", row.capability.as_bytes());
        hash_optional_field(
            &mut hasher,
            b"outcome",
            row.outcome.as_deref().map(str::as_bytes),
        );
        hash_optional_field(&mut hasher, b"diag", row.diag.as_deref().map(str::as_bytes));
        hash_field(&mut hasher, b"exec_mode", exec_mode.as_str().as_bytes());
        let edge_count = u64::try_from(edges.len())
            .expect("edge count fits in u64")
            .to_le_bytes();
        hash_field(&mut hasher, b"edge_count", &edge_count);
        for (role, artifact) in edges {
            hash_field(&mut hasher, b"edge_role", role.as_bytes());
            hash_field(&mut hasher, b"artifact_hash", artifact.as_bytes());
        }
        Ok(hasher.finalize())
    }

    /// The PERSISTED identity of this ledger database (bead gp3.17):
    /// read from the first `vcs-identity` event, minted and recorded on
    /// first use (a hash of the path and the wall clock at minting).
    /// Copies of the file share the identity — they ARE the same ledger
    /// lineage; independent databases get distinct identities, so
    /// commit envelopes and local op ids are never conflated across
    /// instances.
    ///
    /// # Errors
    /// Engine errors.
    pub fn vcs_identity(&self) -> Result<ContentHash, LedgerError> {
        let rows = self
            .conn
            .query(
                "SELECT payload FROM events WHERE kind = 'vcs-identity' \
                 ORDER BY id LIMIT 1",
            )
            .map_err(|e| LedgerError::Sql {
                context: "vcs_identity read".to_string(),
                detail: e.to_string(),
            })?;
        if let Some(row) = rows.first() {
            let payload = match row.get(0) {
                Some(SqliteValue::Text(t)) => t.as_str().to_string(),
                other => {
                    return Err(LedgerError::Corrupt {
                        hash_hex: String::new(),
                        detail: format!("vcs-identity payload: expected TEXT, got {other:?}"),
                    });
                }
            };
            let hex = payload
                .split("\"identity\":\"")
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .ok_or_else(|| LedgerError::Corrupt {
                    hash_hex: String::new(),
                    detail: format!("vcs-identity payload is malformed: {payload}"),
                })?;
            return ContentHash::from_hex(hex).ok_or_else(|| LedgerError::Corrupt {
                hash_hex: String::new(),
                detail: format!("vcs-identity carries a malformed hash: {hex}"),
            });
        }
        let minted = framed_hash(
            LEDGER_IDENTITY_DOMAIN,
            &[
                (b"path", self.path().as_bytes()),
                (b"minted_ns", &crate::now_wall_ns().to_le_bytes()),
            ],
        );
        let payload = format!(
            "{{\"kind\":\"vcs-identity\",\"identity\":\"{}\"}}",
            minted.to_hex()
        );
        self.append_vcs_identity_event(&EventRow {
            session: None,
            t: 0,
            kind: VCS_IDENTITY_EVENT_KIND,
            payload: Some(&payload),
        })?;
        // Re-read: if a concurrent handle minted first, the FIRST event by
        // rowid is the authority for every handle.
        let rows = self
            .conn
            .query(
                "SELECT payload FROM events WHERE kind = 'vcs-identity' \
                 ORDER BY id LIMIT 1",
            )
            .map_err(|e| LedgerError::Sql {
                context: "vcs_identity reread".to_string(),
                detail: e.to_string(),
            })?;
        let payload = match rows.first().and_then(|row| row.get(0)) {
            Some(SqliteValue::Text(t)) => t.as_str().to_string(),
            other => {
                return Err(LedgerError::Corrupt {
                    hash_hex: String::new(),
                    detail: format!("vcs-identity reread: expected TEXT, got {other:?}"),
                });
            }
        };
        let hex = payload
            .split("\"identity\":\"")
            .nth(1)
            .and_then(|rest| rest.split('"').next())
            .ok_or_else(|| LedgerError::Corrupt {
                hash_hex: String::new(),
                detail: format!("vcs-identity payload is malformed: {payload}"),
            })?;
        ContentHash::from_hex(hex).ok_or_else(|| LedgerError::Corrupt {
            hash_hex: String::new(),
            detail: format!("vcs-identity carries a malformed hash: {hex}"),
        })
    }
}

/// Binary Merkle fold over leaves. Pair nodes, odd nodes, and the final root
/// have distinct domains; the root also binds the leaf count so neither tree
/// shape nor a leaf can be confused with a commit root.
fn merkle_root(mut level: Vec<ContentHash>) -> ContentHash {
    let leaf_count = u64::try_from(level.len())
        .expect("leaf count fits in u64")
        .to_le_bytes();
    if level.is_empty() {
        return framed_hash(COMMIT_ROOT_DOMAIN, &[(b"leaf_count", &leaf_count)]);
    }
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        for pair in level.chunks(2) {
            if pair.len() == 2 {
                next.push(framed_hash(
                    MERKLE_PAIR_DOMAIN,
                    &[
                        (b"left", pair[0].as_bytes()),
                        (b"right", pair[1].as_bytes()),
                    ],
                ));
            } else {
                next.push(framed_hash(
                    MERKLE_ODD_DOMAIN,
                    &[(b"child", pair[0].as_bytes())],
                ));
            }
        }
        level = next;
    }
    framed_hash(
        COMMIT_ROOT_DOMAIN,
        &[(b"leaf_count", &leaf_count), (b"tree", level[0].as_bytes())],
    )
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
        let ledger_id = ledger.vcs_identity()?;
        if ledger.branch(branch)?.is_none() {
            return Err(LedgerError::Invalid {
                field: "branch".to_string(),
                problem: format!("branch {branch} does not exist"),
            });
        }
        let ops = ledger.visible_op_ids(branch, None)?;
        let mut leaves = Vec::with_capacity(ops.len());
        let mut frozen_ops = Vec::with_capacity(ops.len());
        let mut artifact_set = BTreeSet::new();
        let mut artifacts = Vec::new();
        for op in &ops {
            let row = ledger.op(*op)?.ok_or_else(|| LedgerError::Corrupt {
                hash_hex: String::new(),
                detail: format!("visible op {op} disappeared while creating a commit"),
            })?;
            if row.outcome.is_none() || row.t_end.is_none() {
                return Err(LedgerError::Invalid {
                    field: "commit".to_string(),
                    problem: format!(
                        "op {op} is still in flight; drain and finalize every op before commit"
                    ),
                });
            }
            leaves.push(ledger.commit_leaf(*op)?);
            for (role, artifact) in ledger.op_artifact_edges(*op)? {
                if role == "out" && artifact_set.insert(artifact) {
                    artifacts.push(artifact);
                }
            }
            frozen_ops.push(row);
        }
        let root = merkle_root(leaves.clone());
        let cutoff_ns = frozen_ops
            .iter()
            .filter_map(|op| op.t_end)
            .max()
            .unwrap_or(0);
        let snapshot = ViewSnapshot {
            branch,
            cutoff_ns,
            ops: frozen_ops,
            in_flight: 0,
            artifacts,
        };
        if let Some(head) = self.heads.get(&(ledger_id, branch))
            && head.root == root
        {
            // Commits identify state, not button presses. Recommitting an
            // unchanged branch is idempotent; recording `root` as its own
            // parent would create a cycle in the commit graph.
            return Ok(head.clone());
        }
        let info = CommitInfo {
            ledger: ledger_id,
            root,
            branch,
            frontier_op: ops.last().copied(),
            parent: self.heads.get(&(ledger_id, branch)).map(|head| head.root),
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
            t: cutoff_ns,
            kind: "vcs-commit",
            payload: Some(&payload),
        })?;
        self.commits.insert(
            (ledger_id, branch, root.0),
            CommitRecord {
                info: info.clone(),
                leaves,
                snapshot,
            },
        );
        self.heads.insert((ledger_id, branch), info.clone());
        Ok(info)
    }

    /// A commit by its FULL envelope identity (bead gp3.17).
    #[must_use]
    pub fn lookup(&self, id: &CommitId) -> Option<&CommitInfo> {
        self.commits
            .get(&(id.ledger, id.branch, id.root.0))
            .map(|record| &record.info)
    }

    /// Every commit envelope sharing a SEMANTIC root — possibly several
    /// branches or ledgers (deterministic order: ledger, then branch).
    #[must_use]
    pub fn lookup_semantic(&self, root: &ContentHash) -> Vec<&CommitInfo> {
        self.commits
            .values()
            .filter(|record| record.info.root == *root)
            .map(|record| &record.info)
            .collect()
    }

    /// The current head of a branch WITHIN a ledger (bead gp3.17: heads
    /// are envelope-scoped; equal roots on other branches or ledgers do
    /// not clobber this pointer).
    ///
    /// # Errors
    /// Engine errors while resolving the ledger identity.
    pub fn head(&self, ledger: &Ledger, branch: i64) -> Result<Option<ContentHash>, LedgerError> {
        let ledger_id = ledger.vcs_identity()?;
        Ok(self.heads.get(&(ledger_id, branch)).map(|head| head.root))
    }

    /// CHECKOUT: return the exact frozen commit view (ops + finished output
    /// artifacts) for a commit ON THIS LEDGER AND BRANCH (bead gp3.17:
    /// snapshots carry ledger-local op ids, so checkout is
    /// envelope-scoped — an equal semantic root from another branch or
    /// ledger is a different commit). Later ops and edges cannot enter
    /// an older snapshot.
    ///
    /// # Errors
    /// A structured error for an unknown commit (the boundary case);
    /// engine errors.
    pub fn checkout(
        &self,
        ledger: &Ledger,
        branch: i64,
        root: &ContentHash,
    ) -> Result<crate::travel::ViewSnapshot, LedgerError> {
        let ledger_id = ledger.vcs_identity()?;
        self.commits
            .get(&(ledger_id, branch, root.0))
            .map(|record| record.snapshot.clone())
            .ok_or_else(|| LedgerError::Invalid {
                field: "commit".to_string(),
                problem: format!(
                    "no commit with root {} is known on this ledger's branch {branch}",
                    root.to_hex()
                ),
            })
    }

    /// The CHECKOUT DELTA between two commits: the symmetric difference
    /// of visible SEMANTIC LEAF multisets — the frontier a delta-solver
    /// reconciles (bead gp3.17: leaf hashes are portable, so the two
    /// commits may come from DIFFERENT ledger instances or from
    /// histories imported in a different row order; local op ids are
    /// reported per side and never compared across ledgers). A nearby
    /// checkout is CHEAP because `shared` dominates.
    ///
    /// # Errors
    /// Unknown commits; engine errors.
    pub fn checkout_delta(
        &self,
        from: &CommitId,
        to: &CommitId,
    ) -> Result<CheckoutDelta, LedgerError> {
        let record = |id: &CommitId| -> Result<&CommitRecord, LedgerError> {
            self.commits
                .get(&(id.ledger, id.branch, id.root.0))
                .ok_or_else(|| LedgerError::Invalid {
                    field: "commit".to_string(),
                    problem: format!(
                        "no commit with root {} is known on ledger {} branch {}",
                        id.root.to_hex(),
                        id.ledger.to_hex(),
                        id.branch
                    ),
                })
        };
        let source = record(from)?;
        let target = record(to)?;
        // Leaf -> that side's local op ids (a multiset: semantically
        // identical ops are distinct instances).
        let index = |rec: &CommitRecord| -> BTreeMap<ContentHash, Vec<i64>> {
            let mut map: BTreeMap<ContentHash, Vec<i64>> = BTreeMap::new();
            for (leaf, op) in rec.leaves.iter().zip(rec.snapshot.ops.iter()) {
                map.entry(*leaf).or_default().push(op.id);
            }
            map
        };
        let a = index(source);
        let b = index(target);
        let mut added = Vec::new();
        let mut removed = Vec::new();
        let mut shared = 0usize;
        let keys: BTreeSet<&ContentHash> = a.keys().chain(b.keys()).collect();
        for leaf in keys {
            let from_ids = a.get(leaf).map_or(&[][..], Vec::as_slice);
            let to_ids = b.get(leaf).map_or(&[][..], Vec::as_slice);
            let common = from_ids.len().min(to_ids.len());
            shared += common;
            for &local_op in &to_ids[common..] {
                added.push(DeltaOp {
                    leaf: *leaf,
                    local_op,
                });
            }
            for &local_op in &from_ids[common..] {
                removed.push(DeltaOp {
                    leaf: *leaf,
                    local_op,
                });
            }
        }
        Ok(CheckoutDelta {
            added,
            removed,
            shared,
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
