//! TOLERANCE-AWARE INVALIDATION (Proposal 2, bead lmp4.7) — the part
//! no build system has. Bit-level invalidation reruns everything
//! downstream of a change; tolerance-level invalidation proves "this
//! change does not matter past stage four" — the Error Ledger run in
//! REVERSE.
//!
//! The semantics: each edge carries a CERTIFIED local sensitivity
//! bound `L` (interval-derived by the supplier). A perturbation `δ`
//! propagates as `Σ L_e · δ(parent)` at each node, flowing through
//! skipped AND recomputed nodes alike (a skipped node's cached bytes
//! are STALE by its bound, and that staleness reaches its consumers
//! scaled by their sensitivities). A node whose accumulated bound fits
//! STRICTLY inside its slack ABSORBS it (skip). The recompute set is
//! then closed UPWARD along delta-carrying edges: recomputing a node
//! is only meaningful with FRESH input bytes, so a recomputing node
//! pulls its stale ancestors in (`PulledByDescendant`). After the
//! closure, recomputed nodes reproduce ground truth exactly and every
//! skip's staleness is bounded by its own slack — the soundness
//! invariant the battery checks against full recomputation.
//!
//! Fail-closed hardening (review round 3): a tie (bound == slack)
//! recomputes, never skips; non-finite sensitivities force recompute;
//! negative-slack nodes are never skippable; `δ = 0` yields an empty
//! frontier. When bounds are loose the frontier balloons and the
//! system DEGRADES GRACEFULLY to ordinary hash memoization — still
//! correct, just less clever — and the SKIP YIELD metric says which
//! ops need tighter bounds (risk R4, owned here).

use crate::{Store, StoreError, json_f64};
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
    /// A sensitivity bound was negative and therefore not a bound:
    /// FAIL CLOSED.
    NegativeSensitivity,
    /// The node's recorded slack is negative (over budget): never
    /// skippable.
    NegativeSlack,
    /// A descendant recomputes and needs FRESH bytes from here
    /// (upward closure along delta-carrying edges).
    PulledByDescendant,
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
    verdicts: Vec<(ContentHash, Verdict)>,
    /// Fraction of TOUCHED nodes certifiably skipped (the R4 health
    /// metric).
    skip_yield: f64,
    /// Canonical JSON rows (verdicts with verified-color skip claims).
    rows: Vec<String>,
    /// Store mutation revision against which every skip certificate was
    /// computed. Private so callers cannot forge plan freshness.
    pub(crate) store_revision: u64,
    /// Fingerprint of the certificate-relevant store state.
    pub(crate) store_state: ContentHash,
}

impl InvalidationPlan {
    /// Read-only semantic verdicts. Keeping the backing storage private
    /// prevents a caller from rewriting recomputes into certified skips.
    #[must_use]
    pub fn verdicts(&self) -> &[(ContentHash, Verdict)] {
        &self.verdicts
    }

    /// Fraction of touched nodes certifiably skipped.
    #[must_use]
    pub fn skip_yield(&self) -> f64 {
        self.skip_yield
    }

    /// Read-only canonical evidence rows.
    #[must_use]
    pub fn rows(&self) -> &[String] {
        &self.rows
    }
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

/// A frontier node's pass-1 state.
struct Touch {
    key: [u8; 32],
    hash: ContentHash,
    bound: f64,
    slack: f64,
    source: bool,
    nonfinite: bool,
    negative_sensitivity: bool,
}

/// Plan the recompute frontier for `perturbations` (node, |δ| bound).
/// Pure: the store is not mutated (see [`apply_plan`]).
///
/// # Errors
/// [`InvalidateError`] teaching errors.
#[allow(clippy::too_many_lines)] // four passes, one algorithm — splitting hides the closure
#[allow(clippy::float_cmp)] // the tie (bound == slack) is DELIBERATELY exact: fail closed
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
    let mut delta: BTreeMap<[u8; 32], f64> = BTreeMap::new();
    let mut sources: BTreeMap<[u8; 32], f64> = BTreeMap::new();
    for &(h, d) in perturbations {
        seq_of(&h)?;
        if d != 0.0 {
            let magnitude = if d.is_finite() {
                d.abs()
            } else {
                f64::INFINITY
            };
            let accumulated = sources.entry(*h.as_bytes()).or_insert(0.0);
            *accumulated += magnitude;
            if !accumulated.is_finite() {
                *accumulated = f64::INFINITY;
            }
        }
    }
    let mut order: Vec<(u64, [u8; 32], ContentHash)> = Vec::new();
    for (h, node) in store.iter() {
        order.push((node.seq, h, node.record.content_hash()));
    }
    order.sort_unstable_by_key(|&(seq, _, _)| seq);
    // Pass 1: deltas flow through EVERY frontier node (skipped nodes
    // are STALE by their bound; recomputed nodes change bytes by at
    // most their bound — either way the staleness reaches consumers).
    let mut touched_nodes: Vec<Touch> = Vec::new();
    for (_, key, hash) in order {
        let source_delta = sources.get(&key).copied();
        let mut bound = source_delta.unwrap_or(0.0);
        let mut nonfinite = false;
        let mut negative_sensitivity = false;
        if let Some(es) = incoming.get(&key) {
            for e in es {
                if let Some(&d) = delta.get(e.from.as_bytes()) {
                    if !e.sensitivity.is_finite() {
                        nonfinite = true;
                    } else if e.sensitivity < 0.0 {
                        negative_sensitivity = true;
                    } else {
                        bound += e.sensitivity * d;
                    }
                }
            }
        }
        if bound == 0.0 && !nonfinite && !negative_sensitivity {
            continue; // untouched: not on the frontier at all
        }
        delta.insert(
            key,
            if bound.is_finite() && !nonfinite && !negative_sensitivity {
                bound
            } else {
                f64::INFINITY
            },
        );
        let node = store.get(&hash).expect("iterated from store");
        touched_nodes.push(Touch {
            key,
            hash,
            bound,
            slack: node.effective_slack(),
            source: source_delta.is_some(),
            nonfinite,
            negative_sensitivity,
        });
    }
    // Pass 2: initial verdict per node (fail-closed hardening).
    let mut reasons: Vec<Option<RecomputeReason>> = touched_nodes
        .iter()
        .map(|t| {
            if t.source {
                Some(RecomputeReason::SourcePerturbed)
            } else if t.nonfinite || !t.bound.is_finite() {
                Some(RecomputeReason::NonFiniteSensitivity)
            } else if t.negative_sensitivity {
                Some(RecomputeReason::NegativeSensitivity)
            } else if t.slack < 0.0 {
                Some(RecomputeReason::NegativeSlack)
            } else if t.bound == t.slack {
                Some(RecomputeReason::TieOnBoundary)
            } else if t.bound < t.slack {
                None // absorbable
            } else {
                Some(RecomputeReason::BoundExceedsSlack)
            }
        })
        .collect();
    // Pass 3: UPWARD CLOSURE — recomputing a node needs FRESH input
    // bytes, so every frontier ancestor along delta-carrying edges
    // joins the recompute set. One reverse-topological sweep reaches
    // the transitive closure (parents have smaller indices).
    let index_of: BTreeMap<[u8; 32], usize> = touched_nodes
        .iter()
        .enumerate()
        .map(|(i, t)| (t.key, i))
        .collect();
    for i in (0..touched_nodes.len()).rev() {
        if reasons[i].is_none() {
            continue;
        }
        if let Some(es) = incoming.get(&touched_nodes[i].key) {
            for e in es {
                if let Some(&pi) = index_of.get(e.from.as_bytes())
                    && reasons[pi].is_none()
                {
                    reasons[pi] = Some(RecomputeReason::PulledByDescendant);
                }
            }
        }
    }
    // Pass 4: emit verdicts and rows.
    let mut verdicts = Vec::new();
    let mut rows = Vec::new();
    let mut skips = 0usize;
    for (t, reason) in touched_nodes.iter().zip(&reasons) {
        let verdict = match reason {
            None => {
                skips += 1;
                Verdict::Skip {
                    bound: t.bound,
                    slack_left: t.slack - t.bound,
                }
            }
            Some(r) => Verdict::Recompute {
                bound: if matches!(
                    r,
                    RecomputeReason::NonFiniteSensitivity | RecomputeReason::NegativeSensitivity
                ) {
                    f64::INFINITY
                } else {
                    t.bound
                },
                reason: *r,
            },
        };
        match &verdict {
            Verdict::Skip { bound, slack_left } => {
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
                    t.hash.to_hex(),
                    t.slack,
                    claim.name(),
                    claim.payload_json()
                );
                rows.push(row);
            }
            Verdict::Recompute { bound, reason } => {
                rows.push(format!(
                    "{{\"event\":\"invalidation\",\"node\":\"{}\",\
                     \"verdict\":\"recompute\",\"bound\":{},\
                     \"reason\":\"{reason:?}\"}}",
                    t.hash.to_hex(),
                    json_f64(*bound)
                ));
            }
        }
        verdicts.push((t.hash, verdict));
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
        store_revision: store.revision(),
        store_state: store.state_fingerprint(),
    })
}

/// Apply a plan to the store: skipped nodes BURN their absorbed bound
/// into runtime state (repeat perturbations see the reduced slack —
/// slack is a real, spendable resource). The plan revision is validated
/// before any burn, so a stale plan cannot partially mutate the store.
///
/// # Errors
/// [`StoreError::StalePlan`] or another fail-closed store error if the
/// state no longer matches the plan's certificate.
pub fn apply_plan(store: &mut Store, plan: &InvalidationPlan) -> Result<(), StoreError> {
    let burns = plan
        .verdicts
        .iter()
        .filter_map(|(hash, verdict)| match verdict {
            Verdict::Skip { bound, .. } => Some((*hash, *bound)),
            Verdict::Recompute { .. } => None,
        })
        .collect::<Vec<_>>();
    store.commit_burns(plan.store_revision, plan.store_state, &burns)
}
