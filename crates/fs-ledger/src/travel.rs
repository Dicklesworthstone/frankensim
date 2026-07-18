//! Time travel, forkable worlds, `explain()`, and the replay audit
//! (plan §11.2; Decalogue P9) — the features that become nearly free
//! because sessions are event-sourced over content-addressed state.
//!
//! - A **fork** is a new branch of the op log sharing every artifact by
//!   hash: N forks cost 1× artifacts + deltas (the storage audit test
//!   proves it). Branch visibility: a branch sees its own ops plus each
//!   ancestor's ops up to the fork point.
//! - **`at_time`** materializes a consistent view at any instant, including
//!   mid-sweep: ops begun by then appear, outcomes not yet written are
//!   masked, and outputs of unfinished ops are not yet visible.
//! - **`explain`** walks the lineage DAG backward from an artifact and
//!   renders the full causal tree (structured + human-readable).
//! - **`replay_verdict`** compares two ledgers op-by-op: `deterministic`
//!   ops must reproduce artifact hashes EXACTLY; `fast` ops may diverge and
//!   are reported separately (the G5-at-study-scale primitive; golden
//!   ledgers in CI use exactly this).

use std::collections::BTreeSet;

use fsqlite::SqliteValue;

use crate::{
    ContentHash, FiveExplicits, Ledger, LedgerError, OpRow, blob_param, opt_blob_param, row_i64,
    sql_err, text_param,
};

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// The root branch every ledger is born with (created by the v2 migration).
pub const MAIN_BRANCH: i64 = 1;

/// Ancestry depth guard: a parent chain longer than this is corruption.
const MAX_BRANCH_DEPTH: usize = 64;

/// Recorded execution mode of an op (plan §5.4): part of provenance, and
/// the contract the replay audit enforces per op.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecMode {
    /// Bit-stable: a replay must reproduce artifact hashes exactly.
    Deterministic,
    /// Throughput mode: replays may diverge; the divergence class is
    /// reported, never silently absorbed.
    Fast,
}

impl ExecMode {
    /// Stable lowercase name (the `ops.exec_mode` column value).
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            ExecMode::Deterministic => "deterministic",
            ExecMode::Fast => "fast",
        }
    }

    /// Parse a column value.
    #[must_use]
    pub fn parse(s: &str) -> Option<ExecMode> {
        match s {
            "deterministic" => Some(ExecMode::Deterministic),
            "fast" => Some(ExecMode::Fast),
            _ => None,
        }
    }
}

/// One branch of the op log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchInfo {
    /// Branch id (rowid; `MAIN_BRANCH` = 1).
    pub id: i64,
    /// Unique human name.
    pub name: String,
    /// Parent branch (None for main).
    pub parent: Option<i64>,
    /// Divergence point: the last parent op id visible from this branch
    /// (None for main).
    pub fork_op: Option<i64>,
    /// Wall-clock creation time (provenance envelope).
    pub created_at: i64,
}

/// A consistent historical view of one branch at a time cutoff.
#[derive(Debug, Clone)]
pub struct ViewSnapshot {
    /// The branch viewed.
    pub branch: i64,
    /// The cutoff instant (wall ns, caller-supplied logical time allowed).
    pub cutoff_ns: i64,
    /// Ops begun by the cutoff, in log order; outcomes not yet written at
    /// the cutoff are masked to in-flight.
    pub ops: Vec<OpRow>,
    /// Ops still in flight at the cutoff.
    pub in_flight: usize,
    /// Artifacts produced (out-edges) by ops FINISHED at the cutoff —
    /// an unfinished op's outputs do not exist yet at that instant.
    pub artifacts: Vec<ContentHash>,
}

/// Set difference between two branches' visible op logs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchDiff {
    /// Op ids visible only on the first branch.
    pub only_a: Vec<i64>,
    /// Op ids visible only on the second branch.
    pub only_b: Vec<i64>,
    /// Ops visible on both (the shared prefix through fork points).
    pub shared: usize,
}

/// One op-pair divergence found by the replay audit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayMismatch {
    /// Position in the visible op sequence (0-based).
    pub position: usize,
    /// Op id in this ledger.
    pub op_a: i64,
    /// Op id in the other ledger.
    pub op_b: i64,
    /// Output hashes (hex) present only in this ledger.
    pub only_a: Vec<String>,
    /// Output hashes (hex) present only in the other ledger.
    pub only_b: Vec<String>,
}

/// Replay audit result (see [`Ledger::replay_verdict`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReplayVerdict {
    /// Op pairs compared.
    pub compared: usize,
    /// Hash mismatches on `deterministic` ops — replay FAILURES.
    pub deterministic_mismatches: Vec<ReplayMismatch>,
    /// Hash divergences on `fast` ops — allowed, reported.
    pub fast_divergences: Vec<ReplayMismatch>,
    /// Set when the op sequences themselves disagree (different length or
    /// any frozen operation-semantic mismatch) — the ledgers do not record
    /// the same study.
    pub structure_mismatch: Option<String>,
}

impl ReplayVerdict {
    /// True when the replay reproduced every deterministic op exactly and
    /// the op sequences agree (fast divergences do not fail a replay).
    #[must_use]
    pub fn is_replay_clean(&self) -> bool {
        self.deterministic_mismatches.is_empty() && self.structure_mismatch.is_none()
    }
}

/// Garbage-collection report for unreferenced artifacts.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GcReport {
    /// Hex hashes of artifacts unreachable from every supported durable root.
    pub candidates: Vec<String>,
    /// Rows actually deleted (0 on dry runs).
    pub deleted: usize,
}

/// One op in an explain tree.
#[derive(Debug, Clone)]
pub struct ExplainOp {
    /// Op id.
    pub id: i64,
    /// Frozen IR (JSON).
    pub ir: String,
    /// Seed bytes as hex.
    pub seed_hex: String,
    /// Frozen versions (JSON — carries the constellation lock hash).
    pub versions: String,
    /// Frozen budget (JSON).
    pub budget: String,
    /// Frozen capability (JSON).
    pub capability: String,
    /// Recorded execution mode.
    pub exec_mode: String,
    /// Final operation outcome, or `None` while the producer is in flight.
    pub outcome: Option<String>,
    /// Structured diagnostic JSON, when the operation recorded one.
    pub diag: Option<String>,
    /// Branch the op ran on.
    pub branch: i64,
    /// The artifacts this op consumed, recursively explained.
    pub inputs: Vec<ExplainNode>,
}

/// One artifact in an explain tree: the full causal ancestry of a number.
#[derive(Debug, Clone)]
pub struct ExplainNode {
    /// Content hash (hex).
    pub hash_hex: String,
    /// Artifact kind.
    pub kind: String,
    /// True when traversal stopped here (depth limit or already expanded
    /// elsewhere in the tree — content addressing makes lineage a DAG).
    pub truncated: bool,
    /// Producing ops (normally one; content addressing dedupes identical
    /// results from independent ops into one artifact row).
    pub produced_by: Vec<ExplainOp>,
}

impl ExplainNode {
    /// Structured rendering (one JSON object; the `explain()` payload).
    #[must_use]
    pub fn to_json(&self) -> String {
        let ops: Vec<String> = self
            .produced_by
            .iter()
            .map(|op| {
                let inputs: Vec<String> = op.inputs.iter().map(ExplainNode::to_json).collect();
                let outcome = op
                    .outcome
                    .as_deref()
                    .map_or_else(|| "null".to_string(), json_string);
                let diag = op.diag.as_deref().unwrap_or("null");
                format!(
                    "{{\"op\":{},\"ir\":{},\"seed\":{},\"versions\":{},\"budget\":{},\
                     \"capability\":{},\"exec_mode\":{},\"outcome\":{},\"diag\":{},\
                     \"branch\":{},\"inputs\":[{}]}}",
                    op.id,
                    op.ir,
                    json_string(&op.seed_hex),
                    op.versions,
                    op.budget,
                    op.capability,
                    json_string(&op.exec_mode),
                    outcome,
                    diag,
                    op.branch,
                    inputs.join(",")
                )
            })
            .collect();
        format!(
            "{{\"artifact\":{},\"kind\":{},\"truncated\":{},\"produced_by\":[{}]}}",
            json_string(&self.hash_hex),
            json_string(&self.kind),
            self.truncated,
            ops.join(",")
        )
    }

    /// Human-readable indented rendering.
    #[must_use]
    pub fn render_text(&self) -> String {
        let mut out = String::new();
        self.render_into(&mut out, 0);
        out
    }

    fn render_into(&self, out: &mut String, depth: usize) {
        use std::fmt::Write as _;
        let pad = "  ".repeat(depth);
        let mark = if self.truncated { " …" } else { "" };
        let _ = writeln!(
            out,
            "{pad}{} [{}]{mark}",
            &self.hash_hex[..16.min(self.hash_hex.len())],
            self.kind
        );
        for op in &self.produced_by {
            let _ = writeln!(
                out,
                "{pad}  <- op {} ({}, branch {}, seed {})",
                op.id, op.exec_mode, op.branch, op.seed_hex
            );
            for input in &op.inputs {
                input.render_into(out, depth + 2);
            }
        }
    }
}

/// One segment of a branch's ancestry: ops on `branch` with id ≤ `cap`
/// (no cap for the branch's own segment).
struct ChainSegment {
    branch: i64,
    cap: Option<i64>,
}

impl Ledger {
    // -- branches -----------------------------------------------------------

    /// Fork a new branch off `parent` at its current op frontier. The fork
    /// shares every artifact by hash; only new ops (and their new
    /// artifacts) cost storage.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown parent;
    /// [`LedgerError::Invalid`] for an empty or duplicate name.
    pub fn fork(&self, name: &str, parent: i64) -> Result<i64, LedgerError> {
        if name.is_empty() {
            return Err(LedgerError::Invalid {
                field: "name".to_string(),
                problem: "empty; branches are addressed by name".to_string(),
            });
        }
        if self.branch(parent)?.is_none() {
            return Err(LedgerError::NotFound {
                what: format!("branch {parent}"),
            });
        }
        let fork_op = self.visible_frontier(parent)?;
        let insert = self
            .conn
            .prepare(
                "INSERT INTO branches(name, parent, fork_op, created_at) \
                 VALUES (?1, ?2, ?3, ?4)",
            )
            .map_err(|e| sql_err("fork prepare", &e))?
            .execute_with_params(&[
                text_param(name),
                SqliteValue::Integer(parent),
                fork_op.map_or(SqliteValue::Null, SqliteValue::Integer),
                SqliteValue::Integer(crate::now_wall_ns()),
            ]);
        match insert {
            Ok(_) => Ok(self.conn.last_insert_rowid()),
            Err(e) if crate::is_duplicate_key(&e) => Err(LedgerError::Invalid {
                field: "name".to_string(),
                problem: format!("branch {name:?} already exists; branch names are unique"),
            }),
            Err(e) => Err(sql_err("fork insert", &e)),
        }
    }

    /// The newest op id visible on `branch` (None when the branch's view is
    /// empty). This is the divergence point recorded by [`Ledger::fork`].
    ///
    /// # Errors
    /// Engine errors; unknown branches error via the ancestry walk.
    pub fn visible_frontier(&self, branch: i64) -> Result<Option<i64>, LedgerError> {
        Ok(self.visible_op_ids(branch, None)?.last().copied())
    }

    /// Fetch one branch, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn branch(&self, id: i64) -> Result<Option<BranchInfo>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT id, name, parent, fork_op, created_at FROM branches WHERE id = ?1",
                &[SqliteValue::Integer(id)],
            )
            .map_err(|e| sql_err("branch fetch", &e))?;
        rows.first().map(branch_from_row).transpose()
    }

    /// Fetch one branch by name, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn branch_by_name(&self, name: &str) -> Result<Option<BranchInfo>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT id, name, parent, fork_op, created_at FROM branches WHERE name = ?1",
                &[text_param(name)],
            )
            .map_err(|e| sql_err("branch fetch", &e))?;
        rows.first().map(branch_from_row).transpose()
    }

    /// All branches, in creation order.
    ///
    /// # Errors
    /// Engine errors.
    pub fn branches(&self) -> Result<Vec<BranchInfo>, LedgerError> {
        let rows = self
            .conn
            .query("SELECT id, name, parent, fork_op, created_at FROM branches ORDER BY id")
            .map_err(|e| sql_err("branch scan", &e))?;
        rows.iter().map(branch_from_row).collect()
    }

    /// Begin an op on a specific branch with an explicit execution mode
    /// (see [`Ledger::begin_op`] for the main-branch deterministic default).
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown branch;
    /// [`LedgerError::MissingExplicit`] naming an offending frozen explicit;
    /// [`LedgerError::Invalid`] for oversized session or IR fields.
    pub fn begin_op_on(
        &self,
        branch: i64,
        mode: ExecMode,
        session: Option<&[u8]>,
        ir: &str,
        explicits: &FiveExplicits<'_>,
        t_start_ns: i64,
    ) -> Result<i64, LedgerError> {
        if let Some(session) = session {
            self.require_op_field_bound(
                "session",
                session.len(),
                crate::MAX_OP_SESSION_BYTES,
                false,
            )?;
        }
        if explicits.seed.is_empty() {
            return Err(LedgerError::MissingExplicit {
                field: "seed".to_string(),
                problem: "empty; record the RNG seed bytes that reproduce this op".to_string(),
            });
        }
        self.require_op_field_bound("seed", explicits.seed.len(), crate::MAX_OP_SEED_BYTES, true)?;
        self.require_bounded_op_json("ir", ir, crate::MAX_OP_IR_BYTES, false)?;
        self.require_bounded_op_json(
            "versions",
            explicits.versions,
            crate::MAX_OP_VERSIONS_BYTES,
            true,
        )?;
        self.require_bounded_op_json("budget", explicits.budget, crate::MAX_OP_BUDGET_BYTES, true)?;
        self.require_bounded_op_json(
            "capability",
            explicits.capability,
            crate::MAX_OP_CAPABILITY_BYTES,
            true,
        )?;
        if self.branch(branch)?.is_none() {
            return Err(LedgerError::NotFound {
                what: format!("branch {branch}"),
            });
        }
        let owns_transaction = !self.in_transaction();
        if owns_transaction {
            self.begin()?;
        }
        let inserted = (|| {
            self.conn
                .prepare(
                    "INSERT INTO ops(session, ir, seed, versions, budget, capability, t_start, \
                     branch, exec_mode) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
                )
                .map_err(|e| sql_err("op insert prepare", &e))?
                .execute_with_params(&[
                    opt_blob_param(session),
                    text_param(ir),
                    blob_param(explicits.seed),
                    text_param(explicits.versions),
                    text_param(explicits.budget),
                    text_param(explicits.capability),
                    SqliteValue::Integer(t_start_ns),
                    SqliteValue::Integer(branch),
                    text_param(mode.as_str()),
                ])
                .map_err(|e| sql_err("op insert", &e))?;
            let op = self.conn.last_insert_rowid();
            self.insert_op_content_identity(op, session, ir, explicits)?;
            Ok(op)
        })();
        match inserted {
            Ok(op) if owns_transaction => match self.commit() {
                Ok(()) => Ok(op),
                Err(error) => {
                    let _ = self.rollback();
                    Err(error)
                }
            },
            Ok(op) => Ok(op),
            Err(error) => {
                if owns_transaction {
                    let _ = self.rollback();
                }
                Err(error)
            }
        }
    }

    // -- visibility ---------------------------------------------------------

    /// The branch's ancestry as (branch, cap) segments: itself uncapped,
    /// each ancestor capped at the fork point below it.
    fn branch_chain(&self, branch: i64) -> Result<Vec<ChainSegment>, LedgerError> {
        let mut segments = Vec::new();
        let mut current = branch;
        let mut cap: Option<i64> = None;
        for _ in 0..MAX_BRANCH_DEPTH {
            let Some(info) = self.branch(current)? else {
                return Err(LedgerError::NotFound {
                    what: format!("branch {current}"),
                });
            };
            segments.push(ChainSegment {
                branch: info.id,
                cap,
            });
            match info.parent {
                None => return Ok(segments),
                Some(parent) => {
                    // The parent's ops are visible only up to the fork
                    // point, further capped by any cap already in force. A
                    // NULL fork point means the parent had no visible ops
                    // at fork time: nothing of it is visible, ever.
                    cap = match (info.fork_op, cap) {
                        (Some(f), Some(c)) => Some(f.min(c)),
                        (Some(f), None) => Some(f),
                        (None, _) => Some(0),
                    };
                    current = parent;
                }
            }
        }
        Err(LedgerError::Corrupt {
            hash_hex: String::new(),
            detail: format!(
                "branch {branch}: ancestry deeper than {MAX_BRANCH_DEPTH} — parent cycle?"
            ),
        })
    }

    /// Op ids visible on `branch` (its own ops plus ancestors' up to each
    /// fork point), ascending, optionally capped at `upto_op`.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown branch.
    pub fn visible_op_ids(
        &self,
        branch: i64,
        upto_op: Option<i64>,
    ) -> Result<Vec<i64>, LedgerError> {
        let mut ids = Vec::new();
        for segment in self.branch_chain(branch)? {
            let cap = match (segment.cap, upto_op) {
                (Some(c), Some(u)) => c.min(u),
                (Some(c), None) => c,
                (None, Some(u)) => u,
                (None, None) => i64::MAX,
            };
            let rows = self
                .conn
                .query_with_params(
                    "SELECT id FROM ops WHERE branch = ?1 AND id <= ?2 ORDER BY id",
                    &[
                        SqliteValue::Integer(segment.branch),
                        SqliteValue::Integer(cap),
                    ],
                )
                .map_err(|e| sql_err("visible ops", &e))?;
            for row in &rows {
                ids.push(row_i64(row, 0, "visible ops")?);
            }
        }
        ids.sort_unstable();
        Ok(ids)
    }

    /// A consistent view of `branch` at instant `t_ns` (mid-sweep points
    /// included): ops begun by then, outcomes masked if not yet written,
    /// outputs visible only for ops finished by then.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown branch.
    pub fn at_time(&self, branch: i64, t_ns: i64) -> Result<ViewSnapshot, LedgerError> {
        let mut ops = Vec::new();
        let mut artifacts = Vec::new();
        let mut seen = BTreeSet::new();
        let mut in_flight = 0usize;
        for id in self.visible_op_ids(branch, None)? {
            let Some(mut op) = self.op(id)? else { continue };
            if op.t_start > t_ns {
                continue;
            }
            let finished_by_cutoff = op.t_end.is_some_and(|t| t <= t_ns);
            if finished_by_cutoff {
                for h in self.op_output_hashes(id)? {
                    if seen.insert(h) {
                        artifacts.push(h);
                    }
                }
            } else {
                // In flight at the cutoff: its recorded ending is the
                // future from this view's perspective.
                op.t_end = None;
                op.outcome = None;
                op.diag = None;
                in_flight += 1;
            }
            ops.push(op);
        }
        Ok(ViewSnapshot {
            branch,
            cutoff_ns: t_ns,
            ops,
            in_flight,
            artifacts,
        })
    }

    /// Output artifact hashes of one op (edges with role `out`).
    fn op_output_hashes(&self, op: i64) -> Result<Vec<ContentHash>, LedgerError> {
        self.op_edge_hashes(op, "out")
    }

    /// Input artifact hashes of one op (edges with role `in`).
    fn op_input_hashes(&self, op: i64) -> Result<Vec<ContentHash>, LedgerError> {
        self.op_edge_hashes(op, "in")
    }

    fn op_edge_hashes(&self, op: i64, role: &'static str) -> Result<Vec<ContentHash>, LedgerError> {
        let (sql, context) = match role {
            "in" => (
                "SELECT artifact FROM edges WHERE op = ?1 AND role = 'in' ORDER BY artifact",
                "op inputs",
            ),
            "out" => (
                "SELECT artifact FROM edges WHERE op = ?1 AND role = 'out' ORDER BY artifact",
                "op outputs",
            ),
            _ => {
                return Err(LedgerError::Invalid {
                    field: "edge.role".to_string(),
                    problem: format!("unsupported internal edge role {role:?}"),
                });
            }
        };
        let rows = self
            .conn
            .query_with_params(sql, &[SqliteValue::Integer(op)])
            .map_err(|e| sql_err(context, &e))?;
        let mut out = Vec::with_capacity(rows.len());
        for row in &rows {
            match row.get(0) {
                Some(SqliteValue::Blob(b)) => match ContentHash::from_slice(b) {
                    Some(h) => out.push(h),
                    None => {
                        return Err(LedgerError::Corrupt {
                            hash_hex: String::new(),
                            detail: format!("op {op}: malformed artifact hash in edges"),
                        });
                    }
                },
                other => {
                    return Err(LedgerError::Sql {
                        context: context.to_string(),
                        detail: format!("artifact: expected BLOB, got {other:?}"),
                    });
                }
            }
        }
        Ok(out)
    }

    /// Set difference between two branches' visible op logs.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown branch.
    pub fn branch_diff(&self, a: i64, b: i64) -> Result<BranchDiff, LedgerError> {
        let sa: BTreeSet<i64> = self.visible_op_ids(a, None)?.into_iter().collect();
        let sb: BTreeSet<i64> = self.visible_op_ids(b, None)?.into_iter().collect();
        Ok(BranchDiff {
            only_a: sa.difference(&sb).copied().collect(),
            only_b: sb.difference(&sa).copied().collect(),
            shared: sa.intersection(&sb).count(),
        })
    }

    // -- explain ------------------------------------------------------------

    /// Walk the lineage DAG backward from an artifact and render its full
    /// causal tree. `max_depth` bounds recursion for enormous campaigns;
    /// repeated sub-lineages are expanded once and marked truncated after
    /// (content addressing makes lineage a DAG, not a tree).
    ///
    /// # Errors
    /// Engine errors; an unknown artifact is `Ok(None)`.
    pub fn explain(
        &self,
        artifact: &ContentHash,
        max_depth: usize,
    ) -> Result<Option<ExplainNode>, LedgerError> {
        let mut expanded = BTreeSet::new();
        self.explain_inner(artifact, max_depth, &mut expanded)
    }

    fn explain_inner(
        &self,
        artifact: &ContentHash,
        depth_left: usize,
        expanded: &mut BTreeSet<ContentHash>,
    ) -> Result<Option<ExplainNode>, LedgerError> {
        let Some(info) = self.artifact_info(artifact)? else {
            return Ok(None);
        };
        let first_visit = expanded.insert(*artifact);
        if depth_left == 0 || !first_visit {
            return Ok(Some(ExplainNode {
                hash_hex: artifact.to_hex(),
                kind: info.kind,
                truncated: true,
                produced_by: Vec::new(),
            }));
        }
        let producer_rows = self
            .conn
            .query_with_params(
                "SELECT op FROM edges WHERE artifact = ?1 AND role = 'out' ORDER BY op",
                &[blob_param(artifact.as_bytes())],
            )
            .map_err(|e| sql_err("explain producers", &e))?;
        let mut produced_by = Vec::new();
        for row in &producer_rows {
            let op_id = row_i64(row, 0, "explain producers")?;
            let Some(op) = self.op(op_id)? else {
                return Err(LedgerError::Corrupt {
                    hash_hex: artifact.to_hex(),
                    detail: format!("edge references missing op {op_id} — lineage is broken"),
                });
            };
            let input_rows = self
                .conn
                .query_with_params(
                    "SELECT artifact FROM edges WHERE op = ?1 AND role = 'in' ORDER BY artifact",
                    &[SqliteValue::Integer(op_id)],
                )
                .map_err(|e| sql_err("explain inputs", &e))?;
            let mut inputs = Vec::new();
            for input_row in &input_rows {
                let b = match input_row.get(0) {
                    Some(SqliteValue::Blob(bytes)) => bytes,
                    other => {
                        return Err(LedgerError::Corrupt {
                            hash_hex: artifact.to_hex(),
                            detail: format!(
                                "op {op_id} has a non-blob input-artifact identity {other:?}"
                            ),
                        });
                    }
                };
                let Some(h) = ContentHash::from_slice(b) else {
                    return Err(LedgerError::Corrupt {
                        hash_hex: artifact.to_hex(),
                        detail: format!(
                            "op {op_id} has a malformed input-artifact identity: expected 32 \
                             bytes, got {}",
                            b.len()
                        ),
                    });
                };
                match self.explain_inner(&h, depth_left - 1, expanded)? {
                    Some(node) => inputs.push(node),
                    None => {
                        return Err(LedgerError::Corrupt {
                            hash_hex: h.to_hex(),
                            detail: format!(
                                "op {op_id} consumed an artifact that does not exist — \
                                 an orphan input in the lineage"
                            ),
                        });
                    }
                }
            }
            let seed_hex: String = op.seed.iter().fold(String::new(), |mut s, b| {
                use std::fmt::Write as _;
                let _ = write!(s, "{b:02x}");
                s
            });
            produced_by.push(ExplainOp {
                id: op.id,
                ir: op.ir,
                seed_hex,
                versions: op.versions,
                budget: op.budget,
                capability: op.capability,
                exec_mode: self.op_exec_mode(op_id)?,
                outcome: op.outcome,
                diag: op.diag,
                branch: self.op_branch(op_id)?,
                inputs,
            });
        }
        Ok(Some(ExplainNode {
            hash_hex: artifact.to_hex(),
            kind: info.kind,
            truncated: false,
            produced_by,
        }))
    }

    fn op_exec_mode(&self, op: i64) -> Result<String, LedgerError> {
        self.bounded_op_exec_mode(op)
    }

    fn op_branch(&self, op: i64) -> Result<i64, LedgerError> {
        let row = self
            .conn
            .query_with_params(
                "SELECT branch FROM ops WHERE id = ?1",
                &[SqliteValue::Integer(op)],
            )
            .map_err(|e| sql_err("op branch", &e))?;
        match row.first() {
            Some(r) => row_i64(r, 0, "op branch"),
            None => Err(LedgerError::NotFound {
                what: format!("op {op}"),
            }),
        }
    }

    // -- replay audit ---------------------------------------------------------

    /// Compare this ledger's branch against another ledger's branch as a
    /// replay: op sequences must agree semantically (IR, all frozen
    /// explicits, execution mode, input lineage, outcome, and diagnostic,
    /// in order); `deterministic` ops must reproduce output hashes exactly;
    /// `fast` divergences are reported without failing the audit. Row ids,
    /// branch ids, sessions, and wall-clock timestamps are provenance
    /// envelopes rather than study semantics and are deliberately excluded.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for unknown branches; engine errors.
    #[allow(clippy::too_many_lines)] // one fail-closed comparison keeps replay semantics visibly ordered
    pub fn replay_verdict(
        &self,
        branch_self: i64,
        other: &Ledger,
        branch_other: i64,
    ) -> Result<ReplayVerdict, LedgerError> {
        let ids_a = self.visible_op_ids(branch_self, None)?;
        let ids_b = other.visible_op_ids(branch_other, None)?;
        let mut verdict = ReplayVerdict::default();
        if ids_a.is_empty() && ids_b.is_empty() {
            verdict.structure_mismatch = Some(
                "both branches contain no operations — no executed study exists to replay"
                    .to_string(),
            );
            return Ok(verdict);
        }
        if ids_a.len() != ids_b.len() {
            verdict.structure_mismatch = Some(format!(
                "op count differs: {} here vs {} in the replay — not the same study",
                ids_a.len(),
                ids_b.len()
            ));
            return Ok(verdict);
        }
        for (position, (&ia, &ib)) in ids_a.iter().zip(ids_b.iter()).enumerate() {
            let (Some(op_a), Some(op_b)) = (self.op(ia)?, other.op(ib)?) else {
                verdict.structure_mismatch = Some(format!("missing op row at position {position}"));
                return Ok(verdict);
            };
            if op_a.outcome.is_none()
                || op_a.t_end.is_none()
                || op_b.outcome.is_none()
                || op_b.t_end.is_none()
            {
                verdict.structure_mismatch = Some(format!(
                    "op at position {position} is still in flight in at least one ledger — \
                     drain and finalize both studies before replay admission"
                ));
                return Ok(verdict);
            }
            let mode_a_raw = self.op_exec_mode(ia)?;
            let mode_b_raw = other.op_exec_mode(ib)?;
            let mode_a = ExecMode::parse(&mode_a_raw).ok_or_else(|| LedgerError::Corrupt {
                hash_hex: String::new(),
                detail: format!("op {ia}: invalid execution mode {mode_a_raw:?}"),
            })?;
            let mode_b = ExecMode::parse(&mode_b_raw).ok_or_else(|| LedgerError::Corrupt {
                hash_hex: String::new(),
                detail: format!("replay op {ib}: invalid execution mode {mode_b_raw:?}"),
            })?;
            let inputs_a = self.op_input_hashes(ia)?;
            let inputs_b = other.op_input_hashes(ib)?;
            let mut changed = Vec::new();
            if op_a.ir != op_b.ir {
                changed.push("ir");
            }
            if op_a.seed != op_b.seed {
                changed.push("seed");
            }
            if op_a.versions != op_b.versions {
                changed.push("versions");
            }
            if op_a.budget != op_b.budget {
                changed.push("budget");
            }
            if op_a.capability != op_b.capability {
                changed.push("capability");
            }
            if mode_a != mode_b {
                changed.push("exec_mode");
            }
            if inputs_a != inputs_b {
                changed.push("input_lineage");
            }
            if op_a.outcome != op_b.outcome {
                changed.push("outcome");
            }
            if op_a.diag != op_b.diag {
                changed.push("diag");
            }
            if !changed.is_empty() {
                verdict.structure_mismatch = Some(format!(
                    "op semantics diverge at position {position}: {} differ \
                     — the replay ran a different study",
                    changed.join(", ")
                ));
                return Ok(verdict);
            }
            let out_a: BTreeSet<String> = self
                .op_output_hashes(ia)?
                .iter()
                .map(ContentHash::to_hex)
                .collect();
            let out_b: BTreeSet<String> = other
                .op_output_hashes(ib)?
                .iter()
                .map(ContentHash::to_hex)
                .collect();
            verdict.compared += 1;
            if out_a != out_b {
                let mismatch = ReplayMismatch {
                    position,
                    op_a: ia,
                    op_b: ib,
                    only_a: out_a.difference(&out_b).cloned().collect(),
                    only_b: out_b.difference(&out_a).cloned().collect(),
                };
                if mode_a == ExecMode::Fast {
                    verdict.fast_divergences.push(mismatch);
                } else {
                    verdict.deterministic_mismatches.push(mismatch);
                }
            }
        }
        Ok(verdict)
    }

    // -- garbage collection ---------------------------------------------------

    /// Artifacts with neither a lineage edge on ANY branch, an immutable
    /// solver-checkpoint receipt, either side of an immutable quantity
    /// dimension crosswalk, an immutable semantic state-checkpoint receipt,
    /// nor an explicit artifact-semantic binding are unreachable from every
    /// supported root and safe to reclaim.
    /// Referenced artifacts are immortal (Decalogue P9). Dry runs only report.
    ///
    /// # Errors
    /// Engine errors; on failure during deletion the transaction rolls back.
    pub fn gc_unreferenced_artifacts(&self, dry_run: bool) -> Result<GcReport, LedgerError> {
        let rows = self
            .conn
            .query(
                "SELECT a.hash FROM artifacts a \
                 LEFT JOIN edges e ON a.hash = e.artifact \
                 LEFT JOIN session_checkpoint_receipts c \
                    ON a.hash = c.solver_state_artifact \
                 LEFT JOIN qty_dimension_crosswalks q_old \
                    ON a.hash = q_old.old_hash \
                 LEFT JOIN qty_dimension_crosswalks q_new \
                    ON a.hash = q_new.new_hash \
                 LEFT JOIN semantic_state_checkpoint_receipts s \
                    ON a.hash = s.runtime_state_artifact \
                 LEFT JOIN artifact_semantic_bindings identity_binding \
                    ON a.hash = identity_binding.artifact_hash \
                 WHERE e.artifact IS NULL \
                   AND c.solver_state_artifact IS NULL \
                   AND q_old.old_hash IS NULL \
                   AND q_new.new_hash IS NULL \
                   AND s.runtime_state_artifact IS NULL \
                   AND identity_binding.artifact_hash IS NULL \
                 ORDER BY a.hash",
            )
            .map_err(|e| sql_err("gc scan", &e))?;
        let mut candidates = Vec::with_capacity(rows.len());
        for row in &rows {
            if let Some(SqliteValue::Blob(b)) = row.get(0)
                && let Some(h) = ContentHash::from_slice(b)
            {
                candidates.push(h.to_hex());
            }
        }
        if dry_run || candidates.is_empty() {
            return Ok(GcReport {
                candidates,
                deleted: 0,
            });
        }
        let owns_txn = !self.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        let result = (|| -> Result<usize, LedgerError> {
            let mut deleted = 0usize;
            for hex in &candidates {
                let Some(h) = ContentHash::from_hex(hex) else {
                    continue;
                };
                self.conn
                    .prepare("DELETE FROM artifact_chunks WHERE hash = ?1")
                    .map_err(|e| sql_err("gc chunks", &e))?
                    .execute_with_params(&[blob_param(h.as_bytes())])
                    .map_err(|e| sql_err("gc chunks", &e))?;
                deleted += self
                    .conn
                    .prepare("DELETE FROM artifacts WHERE hash = ?1")
                    .map_err(|e| sql_err("gc artifact", &e))?
                    .execute_with_params(&[blob_param(h.as_bytes())])
                    .map_err(|e| sql_err("gc artifact", &e))?;
            }
            Ok(deleted)
        })();
        match (result, owns_txn) {
            (Ok(deleted), true) => {
                self.commit()?;
                Ok(GcReport {
                    candidates,
                    deleted,
                })
            }
            (Ok(deleted), false) => Ok(GcReport {
                candidates,
                deleted,
            }),
            (Err(e), owns) => {
                if owns {
                    let _ = self.rollback();
                }
                Err(e)
            }
        }
    }
}

fn branch_from_row(row: &fsqlite::Row) -> Result<BranchInfo, LedgerError> {
    let name = match row.get(1) {
        Some(SqliteValue::Text(t)) => t.as_str().to_string(),
        other => {
            return Err(LedgerError::Sql {
                context: "branch fetch".to_string(),
                detail: format!("name: expected TEXT, got {other:?}"),
            });
        }
    };
    let opt_int = |idx: usize| -> Option<i64> {
        match row.get(idx) {
            Some(SqliteValue::Integer(v)) => Some(*v),
            _ => None,
        }
    };
    Ok(BranchInfo {
        id: row_i64(row, 0, "branch.id")?,
        name,
        parent: opt_int(2),
        fork_op: opt_int(3),
        created_at: row_i64(row, 4, "branch.created_at")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EdgeRole, Ledger, OpOutcome};
    use fs_blake3::hash_bytes;

    const FX: FiveExplicits<'static> = FiveExplicits {
        seed: &[0xAB],
        versions: "{}",
        budget: "{}",
        capability: "{}",
    };

    fn mem() -> Ledger {
        Ledger::open(":memory:").expect("open")
    }

    /// One complete op producing one artifact with deterministic bytes.
    fn unit(l: &Ledger, branch: i64, mode: ExecMode, tag: u64, t: i64) -> ContentHash {
        let op = l
            .begin_op_on(branch, mode, None, &format!("{{\"tag\":{tag}}}"), &FX, t)
            .expect("begin");
        let bytes: Vec<u8> = (0..64u64).map(|i| ((tag * 131 + i) % 251) as u8).collect();
        let receipt = l.put_artifact("fixture", &bytes, None).expect("put");
        l.link(op, &receipt.hash, EdgeRole::Out).expect("link");
        l.finish_op(op, OpOutcome::Ok, None, t + 5).expect("finish");
        receipt.hash
    }

    #[derive(Clone, Copy)]
    struct ReplayFixture<'a> {
        session: Option<&'a [u8]>,
        ir: &'a str,
        seed: &'a [u8],
        versions: &'a str,
        budget: &'a str,
        capability: &'a str,
        mode: ExecMode,
        input: &'a [u8],
        outcome: OpOutcome,
        diag: Option<&'a str>,
        t: i64,
    }

    const REPLAY_FIXTURE: ReplayFixture<'static> = ReplayFixture {
        session: Some(b"session-a"),
        ir: "{\"op\":\"solve\"}",
        seed: b"seed-a",
        versions: "{\"solver\":1}",
        budget: "{\"wall_ns\":10}",
        capability: "{\"cores\":1}",
        mode: ExecMode::Deterministic,
        input: b"input-a",
        outcome: OpOutcome::Ok,
        diag: None,
        t: 10,
    };

    fn replay_fixture(spec: ReplayFixture<'_>) -> Ledger {
        let ledger = mem();
        let explicits = FiveExplicits {
            seed: spec.seed,
            versions: spec.versions,
            budget: spec.budget,
            capability: spec.capability,
        };
        let op = ledger
            .begin_op_on(
                MAIN_BRANCH,
                spec.mode,
                spec.session,
                spec.ir,
                &explicits,
                spec.t,
            )
            .expect("begin replay fixture");
        let input = ledger
            .put_artifact("fixture-input", spec.input, None)
            .expect("put replay input");
        ledger
            .link(op, &input.hash, EdgeRole::In)
            .expect("link replay input");
        let output = ledger
            .put_artifact("fixture-output", b"stable-output", None)
            .expect("put replay output");
        ledger
            .link(op, &output.hash, EdgeRole::Out)
            .expect("link replay output");
        ledger
            .finish_op(op, spec.outcome, spec.diag, spec.t + 1)
            .expect("finish replay fixture");
        ledger
    }

    #[test]
    fn main_branch_exists_after_migration() {
        let l = mem();
        let main = l.branch(MAIN_BRANCH).unwrap().expect("main branch");
        assert_eq!(main.name, "main");
        assert_eq!(main.parent, None);
        assert!(l.branch_by_name("main").unwrap().is_some());
    }

    #[test]
    fn fork_visibility_and_independence() {
        let l = mem();
        unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 1, 10);
        unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 2, 20);
        let a = l.fork("branch-a", MAIN_BRANCH).unwrap();
        let b = l.fork("branch-b", MAIN_BRANCH).unwrap();
        unit(&l, a, ExecMode::Deterministic, 3, 30);
        unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 4, 40);
        // A sees the shared prefix + its own op, NOT main's post-fork op.
        assert_eq!(l.visible_op_ids(a, None).unwrap(), vec![1, 2, 3]);
        // B sees only the shared prefix.
        assert_eq!(l.visible_op_ids(b, None).unwrap(), vec![1, 2]);
        // Main sees its own log, not A's op.
        assert_eq!(l.visible_op_ids(MAIN_BRANCH, None).unwrap(), vec![1, 2, 4]);
        let diff = l.branch_diff(a, MAIN_BRANCH).unwrap();
        assert_eq!(diff.only_a, vec![3]);
        assert_eq!(diff.only_b, vec![4]);
        assert_eq!(diff.shared, 2);
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn fork_rejects_bad_inputs() {
        let l = mem();
        assert_eq!(l.fork("x", 999).unwrap_err().code(), "LedgerNotFound");
        assert_eq!(l.fork("", MAIN_BRANCH).unwrap_err().code(), "LedgerInvalid");
        l.fork("dup", MAIN_BRANCH).unwrap();
        assert_eq!(
            l.fork("dup", MAIN_BRANCH).unwrap_err().code(),
            "LedgerInvalid"
        );
        assert_eq!(
            l.begin_op_on(999, ExecMode::Fast, None, "{}", &FX, 1)
                .unwrap_err()
                .code(),
            "LedgerNotFound"
        );
    }

    #[test]
    fn replay_binds_every_operation_semantic_field() {
        let original = replay_fixture(REPLAY_FIXTURE);
        let mutations = [
            (
                "ir",
                ReplayFixture {
                    ir: "{\"op\":\"different\"}",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "seed",
                ReplayFixture {
                    seed: b"seed-b",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "versions",
                ReplayFixture {
                    versions: "{\"solver\":2}",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "budget",
                ReplayFixture {
                    budget: "{\"wall_ns\":11}",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "capability",
                ReplayFixture {
                    capability: "{\"cores\":2}",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "exec_mode",
                ReplayFixture {
                    mode: ExecMode::Fast,
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "input_lineage",
                ReplayFixture {
                    input: b"input-b",
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "outcome",
                ReplayFixture {
                    outcome: OpOutcome::Error,
                    ..REPLAY_FIXTURE
                },
            ),
            (
                "diag",
                ReplayFixture {
                    diag: Some("{\"code\":\"changed\"}"),
                    ..REPLAY_FIXTURE
                },
            ),
        ];
        for (field, mutation) in mutations {
            let replay = replay_fixture(mutation);
            let verdict = original
                .replay_verdict(MAIN_BRANCH, &replay, MAIN_BRANCH)
                .expect("replay verdict");
            assert!(
                !verdict.is_replay_clean(),
                "semantic mutation {field} passed replay"
            );
            assert!(
                verdict
                    .structure_mismatch
                    .as_deref()
                    .is_some_and(|message| message.contains(field)),
                "semantic mutation {field} was not diagnosed: {verdict:?}"
            );
        }
    }

    #[test]
    fn replay_excludes_run_envelope_fields() {
        let original = replay_fixture(REPLAY_FIXTURE);
        let replay = replay_fixture(ReplayFixture {
            session: Some(b"session-b"),
            t: 99_000,
            ..REPLAY_FIXTURE
        });
        let verdict = original
            .replay_verdict(MAIN_BRANCH, &replay, MAIN_BRANCH)
            .expect("replay verdict");
        assert!(verdict.is_replay_clean(), "envelope drift: {verdict:?}");
        assert_eq!(verdict.compared, 1);
    }

    #[test]
    fn replay_never_admits_two_unfinished_studies() {
        let original = mem();
        let replay = mem();
        original
            .begin_op_on(
                MAIN_BRANCH,
                ExecMode::Deterministic,
                None,
                "{\"op\":\"unfinished\"}",
                &FX,
                1,
            )
            .expect("begin original");
        replay
            .begin_op_on(
                MAIN_BRANCH,
                ExecMode::Deterministic,
                None,
                "{\"op\":\"unfinished\"}",
                &FX,
                2,
            )
            .expect("begin replay");

        let verdict = original
            .replay_verdict(MAIN_BRANCH, &replay, MAIN_BRANCH)
            .expect("replay verdict");
        assert!(!verdict.is_replay_clean(), "unfinished replay was admitted");
        assert!(
            verdict
                .structure_mismatch
                .as_deref()
                .is_some_and(|message| message.contains("still in flight")),
            "unexpected refusal: {verdict:?}"
        );
        assert_eq!(verdict.compared, 0);
    }

    #[test]
    fn replay_never_admits_two_empty_branches() {
        let original = mem();
        let replay = mem();

        let verdict = original
            .replay_verdict(MAIN_BRANCH, &replay, MAIN_BRANCH)
            .expect("replay verdict");
        assert!(!verdict.is_replay_clean(), "empty replay was admitted");
        assert_eq!(verdict.compared, 0);
        assert!(
            verdict
                .structure_mismatch
                .as_deref()
                .is_some_and(|message| message.contains("no executed study")),
            "unexpected refusal: {verdict:?}"
        );
    }

    #[test]
    fn explain_json_escapes_hostile_artifact_kinds() {
        let ledger = mem();
        let op = ledger
            .begin_op_on(
                MAIN_BRANCH,
                ExecMode::Deterministic,
                None,
                "{\"op\":\"explain\"}",
                &FX,
                1,
            )
            .expect("begin");
        let artifact = ledger
            .put_artifact("kind\"\\\n\u{0001}", b"payload", None)
            .expect("artifact");
        ledger
            .link(op, &artifact.hash, EdgeRole::Out)
            .expect("link");
        ledger
            .finish_op(
                op,
                OpOutcome::Ok,
                Some("{\"message\":\"quoted \\\"diagnostic\\\"\"}"),
                2,
            )
            .expect("finish");

        let json = ledger
            .explain(&artifact.hash, 4)
            .expect("explain")
            .expect("node")
            .to_json();
        assert!(json.starts_with('{') && json.ends_with('}'), "{json}");
        assert!(json.contains("\\u0001"), "{json}");
        assert!(json.contains("\"outcome\":\"ok\""), "{json}");
        assert!(json.contains("\"diag\":{\"message\":"), "{json}");
        assert!(!json.contains('\u{0001}'), "{json:?}");
        assert!(!json.chars().any(|ch| u32::from(ch) < 0x20), "{json:?}");
        ledger
            .append_event(&crate::EventRow {
                session: None,
                t: 3,
                kind: "explain-json",
                payload: Some(&json),
            })
            .expect("SQLite accepts strict explain JSON");
    }

    #[test]
    fn explain_refuses_malformed_input_edge_identity() {
        let ledger = mem();
        let input = ledger
            .put_artifact("input", b"input", None)
            .expect("input artifact");
        let op = ledger
            .begin_op_on(
                MAIN_BRANCH,
                ExecMode::Deterministic,
                None,
                "{\"op\":\"corrupt-lineage\"}",
                &FX,
                1,
            )
            .expect("begin");
        ledger
            .link(op, &input.hash, EdgeRole::In)
            .expect("input link");
        let output = ledger
            .put_artifact("output", b"output", None)
            .expect("output artifact");
        ledger
            .link(op, &output.hash, EdgeRole::Out)
            .expect("output link");
        ledger
            .finish_op(op, OpOutcome::Ok, None, 2)
            .expect("finish");

        ledger
            .conn
            .query("PRAGMA foreign_keys=OFF")
            .expect("disable foreign keys for corruption fixture");
        ledger
            .conn
            .execute_with_params(
                "UPDATE edges SET artifact = ?1 WHERE op = ?2 AND role = 'in'",
                &[
                    SqliteValue::Blob(vec![0_u8].into()),
                    SqliteValue::Integer(op),
                ],
            )
            .expect("inject malformed input identity");

        let error = ledger
            .explain(&output.hash, 4)
            .expect_err("must fail closed");
        assert_eq!(error.code(), "LedgerCorruption", "{error}");
        assert!(error.to_string().contains("expected 32 bytes"), "{error}");
    }

    #[test]
    fn at_time_masks_in_flight_outcomes() {
        let l = mem();
        // op1: t 10..15; op2: t 20..35; op3 starts at 40.
        unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 1, 10);
        let op2 = l
            .begin_op_on(MAIN_BRANCH, ExecMode::Deterministic, None, "{}", &FX, 20)
            .unwrap();
        let r2 = l.put_artifact("fixture", b"op2 output", None).unwrap();
        l.link(op2, &r2.hash, EdgeRole::Out).unwrap();
        l.finish_op(op2, OpOutcome::Ok, None, 35).unwrap();
        let _op3 = l
            .begin_op_on(MAIN_BRANCH, ExecMode::Deterministic, None, "{}", &FX, 40)
            .unwrap();

        let view = l.at_time(MAIN_BRANCH, 25).unwrap();
        assert_eq!(view.ops.len(), 2, "op3 has not begun at t=25");
        assert_eq!(view.in_flight, 1, "op2 is mid-flight at t=25");
        let masked = view.ops.iter().find(|o| o.id == op2).unwrap();
        assert_eq!(
            masked.outcome, None,
            "outcome written at t=35 is the future"
        );
        // op2's output does not exist yet at t=25; only op1's artifact shows.
        assert_eq!(view.artifacts.len(), 1);
        let full = l.at_time(MAIN_BRANCH, 100).unwrap();
        assert_eq!(full.in_flight, 1, "op3 never finished");
        assert_eq!(full.artifacts.len(), 2);
    }

    #[test]
    fn gc_reclaims_only_unreferenced() {
        let l = mem();
        let kept = unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 1, 10);
        let loose = l.put_artifact("scratch", b"never linked", None).unwrap();
        let run = fs_exec::RunId(0x6c);
        let gate = fs_exec::CancelGate::new_clock_free();
        let tracker = fs_exec::DrainTracker::new(run, &gate);
        let worker = tracker.register_worker().expect("checkpoint worker");
        gate.request();
        drop(worker);
        let report = tracker.finalize().expect("checkpoint worker drained");
        let snapshot = fs_exec::solver::envelope::seal(7, 1, run.0, b"rooted checkpoint");
        let checkpoint = l
            .put_artifact(
                crate::session_registry::SOLVER_STATE_ARTIFACT_KIND,
                &snapshot,
                None,
            )
            .expect("checkpoint artifact");
        l.attest_solver_checkpoint(crate::session_registry::SolverCheckpointClaim {
            session: 1,
            pause_authority: hash_bytes(b"gc-checkpoint-authority"),
            gate_generation: 0,
            solver_state_artifact: checkpoint.hash,
            drain_report: &report,
        })
        .expect("checkpoint receipt");
        let dry = l.gc_unreferenced_artifacts(true).unwrap();
        assert_eq!(dry.candidates, vec![loose.hash.to_hex()]);
        assert_eq!(dry.deleted, 0);
        assert!(
            l.artifact_info(&loose.hash).unwrap().is_some(),
            "dry run deletes nothing"
        );
        let real = l.gc_unreferenced_artifacts(false).unwrap();
        assert_eq!(real.deleted, 1);
        assert!(l.artifact_info(&loose.hash).unwrap().is_none());
        assert!(
            l.artifact_info(&kept).unwrap().is_some(),
            "referenced artifacts are immortal"
        );
        assert!(
            l.artifact_info(&checkpoint.hash).unwrap().is_some(),
            "checkpoint receipts root solver-state artifacts"
        );
        assert!(l.lint().unwrap().is_clean());
    }
}
