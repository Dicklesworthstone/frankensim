//! TOLERANCE-AWARE INVALIDATION (Proposal 2, bead lmp4.7) — the part
//! no build system has. Bit-level invalidation reruns everything
//! downstream of a change; tolerance-level invalidation proves "this
//! change does not matter past stage four" — the Error Ledger run in
//! REVERSE.
//!
//! The semantics: each edge carries a CERTIFIED local sensitivity
//! bound `L` (interval-derived by the supplier). A perturbation `δ`
//! propagates as `Σ L_e · δ(parent)` at each node. A node whose
//! accumulated bound fits STRICTLY inside its recorded slack is
//! SKIPPED — its cached artifact is reused, the bound burns into its
//! slack, and because its artifact BYTES are unchanged, propagation
//! STOPS there. A recomputed node passes its bound downstream.
//!
//! Fail-closed hardening (review round 3): a tie (bound == slack)
//! recomputes, never skips; non-finite sensitivities force recompute;
//! negative-slack nodes are never skippable; `δ = 0` yields an empty
//! frontier. When bounds are loose the frontier balloons and the
//! system DEGRADES GRACEFULLY to ordinary hash memoization — still
//! correct, just less clever — and the SKIP YIELD metric says which
//! ops need tighter bounds (risk R4, owned here).

use crate::{Store, StoreError};
use fs_evidence::Color;
use fs_ledger::ContentHash;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// One DAG edge with its CERTIFIED sensitivity bound (`|∂out/∂in|`
/// over the perturbation box, interval-derived by the supplier).
#[derive(Debug, Clone, Copy)]
pub struct Edge {
    /// Upstream node.
    pub from: ContentHash,
    /// Downstream node.
    pub to: ContentHash,
    /// Certified local sensitivity bound.
    pub sensitivity: f64,
}

/// Why a node recomputes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecomputeReason {
    /// The perturbation entered HERE (the edit itself).
    SourcePerturbed,
    /// The accumulated bound exceeds the slack.
    BoundExceedsSlack,
    /// Bound EXACTLY equals slack: ties resolve to recompute, never
    /// skip (deterministic fail-closed boundary).
    TieOnBoundary,
    /// An incoming sensitivity was non-finite: FAIL CLOSED.
    NonFiniteSensitivity,
    /// The node's recorded slack is negative (over budget): never
    /// skippable.
    NegativeSlack,
}

/// The per-node verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum Verdict {
    /// Reuse the cached artifact: the perturbation is ABSORBED. The
    /// bound is the certificate; `slack_left` is what remains after
    /// burning it.
    Skip {
        /// Accumulated incoming bound (the certified staleness).
        bound: f64,
        /// Slack remaining after absorption.
        slack_left: f64,
    },
    /// Recompute, and why.
    Recompute {
        /// Accumulated incoming bound.
        bound: f64,
        /// The reason.
        reason: RecomputeReason,
    },
}

/// The invalidation plan (pure — apply separately).
#[derive(Debug, Clone)]
pub struct InvalidationPlan {
    /// Verdicts in store (topological) order, touched nodes only.
    pub verdicts: Vec<(ContentHash, Verdict)>,
    /// Fraction of TOUCHED nodes certifiably skipped (the R4 health
    /// metric).
    pub skip_yield: f64,
    /// Canonical JSON rows (verdicts with verified-color skip claims).
    pub rows: Vec<String>,
}

/// Planning errors.
#[derive(Debug, Clone, PartialEq)]
pub enum InvalidateError {
    /// An edge references a node the store does not hold.
    UnknownNode {
        /// The hash.
        node: ContentHash,
    },
    /// Edges must run from earlier to later store sequence (a DAG in
    /// insertion order).
    NotTopological {
        /// The offending edge's endpoints.
        from: ContentHash,
        /// Downstream end.
        to: ContentHash,
    },
}

impl core::fmt::Display for InvalidateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            InvalidateError::UnknownNode { node } => {
                write!(f, "edge references unknown node {}", node.to_hex())
            }
            InvalidateError::NotTopological { from, to } => write!(
                f,
                "edge {} -> {} runs against insertion order; register nodes \
                 topologically",
                from.to_hex(),
                to.to_hex()
            ),
        }
    }
}

impl std::error::Error for InvalidateError {}

/// Plan the recompute frontier for `perturbations` (node, |δ| bound).
/// Pure: the store is not mutated (see [`apply_plan`]).
///
/// # Errors
/// [`InvalidateError`] teaching errors.
pub fn plan(
    store: &Store,
    edges: &[Edge],
    perturbations: &[(ContentHash, f64)],
) -> Result<InvalidationPlan, InvalidateError> {
    // Validate + index edges by downstream node, topologically.
    let seq_of = |h: &ContentHash| -> Result<u64, InvalidateError> {
        store
            .get(h)
            .map(|n| n.seq)
            .ok_or(InvalidateError::UnknownNode { node: *h })
    };
    let mut incoming: BTreeMap<[u8; 32], Vec<Edge>> = BTreeMap::new();
    for e in edges {
        let (sf, st) = (seq_of(&e.from)?, seq_of(&e.to)?);
        if sf >= st {
            return Err(InvalidateError::NotTopological {
                from: e.from,
                to: e.to,
            });
        }
        incoming.entry(*e.to.as_bytes()).or_default().push(*e);
    }
    // Seed deltas; process nodes in store sequence order (topological
    // by validation above).
    let mut delta: BTreeMap<[u8; 32], f64> = BTreeMap::new();
    let mut sources: BTreeMap<[u8; 32], f64> = BTreeMap::new();
    for &(h, d) in perturbations {
        seq_of(&h)?;
        if d != 0.0 {
            sources.insert(*h.as_bytes(), d.abs());
        }
    }
    let mut order: Vec<(u64, [u8; 32], ContentHash)> = Vec::new();
    for (h, node) in store.iter() {
        order.push((node.seq, h, node.record.content_hash()));
    }
    order.sort_unstable_by_key(|&(seq, _, _)| seq);
    let mut verdicts = Vec::new();
    let mut rows = Vec::new();
    let mut skips = 0usize;
    for (_, key, hash) in order {
        let source_delta = sources.get(&key).copied();
        // Accumulate incoming bound from RECOMPUTED/perturbed parents.
        let mut bound = source_delta.unwrap_or(0.0);
        let mut nonfinite = false;
        if let Some(es) = incoming.get(&key) {
            for e in es {
                if let Some(&d) = delta.get(e.from.as_bytes()) {
                    if !e.sensitivity.is_finite() {
                        nonfinite = true;
                    }
                    bound += e.sensitivity * d;
                }
            }
        }
        if bound == 0.0 && !nonfinite {
            continue; // untouched: not on the frontier at all
        }
        let node = store.get(&hash).expect("iterated from store");
        let slack = node.effective_slack();
        let verdict = if source_delta.is_some() {
            Verdict::Recompute {
                bound,
                reason: RecomputeReason::SourcePerturbed,
            }
        } else if nonfinite || !bound.is_finite() {
            Verdict::Recompute {
                bound: f64::INFINITY,
                reason: RecomputeReason::NonFiniteSensitivity,
            }
        } else if slack < 0.0 {
            Verdict::Recompute {
                bound,
                reason: RecomputeReason::NegativeSlack,
            }
        } else if bound == slack {
            Verdict::Recompute {
                bound,
                reason: RecomputeReason::TieOnBoundary,
            }
        } else if bound < slack {
            Verdict::Skip {
                bound,
                slack_left: slack - bound,
            }
        } else {
            Verdict::Recompute {
                bound,
                reason: RecomputeReason::BoundExceedsSlack,
            }
        };
        match &verdict {
            Verdict::Skip { bound, slack_left } => {
                skips += 1;
                // The skip claim is VERIFIED-color: the bound is an
                // interval certificate over the perturbation.
                let claim = Color::Verified {
                    lo: 0.0,
                    hi: *bound,
                };
                let mut row = String::new();
                let _ = write!(
                    row,
                    "{{\"event\":\"invalidation\",\"node\":\"{}\",\"verdict\":\"skip\",\
                     \"claim\":\"skipped: perturbation absorbed, bound {bound:.3e} <= \
                     slack {:.3e}\",\"color\":\"{}\",\"payload\":{},\
                     \"slack_left\":{slack_left:.3e}}}",
                    hash.to_hex(),
                    slack,
                    claim.name(),
                    claim.payload_json()
                );
                rows.push(row);
                // Absorbed: downstream sees UNCHANGED bytes.
            }
            Verdict::Recompute { bound, reason } => {
                rows.push(format!(
                    "{{\"event\":\"invalidation\",\"node\":\"{}\",\
                     \"verdict\":\"recompute\",\"bound\":{bound:.3e},\
                     \"reason\":\"{reason:?}\"}}",
                    hash.to_hex()
                ));
                // A recomputed node passes its output delta downstream.
                delta.insert(key, *bound);
            }
        }
        verdicts.push((hash, verdict));
    }
    let touched = verdicts.len();
    Ok(InvalidationPlan {
        skip_yield: if touched == 0 {
            1.0
        } else {
            skips as f64 / touched as f64
        },
        verdicts,
        rows,
    })
}

/// Apply a plan to the store: skipped nodes BURN their absorbed bound
/// into `achieved_error` (repeat perturbations see the reduced slack —
/// slack is a real, spendable resource).
///
/// # Errors
/// [`StoreError::UnknownNode`] if the store changed under the plan.
pub fn apply_plan(store: &mut Store, plan: &InvalidationPlan) -> Result<(), StoreError> {
    for (hash, verdict) in &plan.verdicts {
        if let Verdict::Skip { bound, .. } = verdict {
            store.burn_slack(hash, *bound)?;
        }
    }
    Ok(())
}
