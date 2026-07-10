//! Proposal 2's operator-facing API + cache policy (bead lmp4.8) —
//! deliberately headline-status, not plumbing: agent-driven design
//! exploration IS thousands of nearby variants, and the delta between
//! a thousand full re-simulations and a thousand certified delta-solves
//! is the whole economic difference of agentic operation.
//!
//! - [`RecomputeApi::perturb`] returns a FIRST-CLASS plan: the minimal
//!   recompute frontier, its estimated cost (Proposal 8's planner
//!   input), and the verified-color certificates for everything
//!   skipped. Inspect and cost it BEFORE committing.
//! - [`RecomputeApi::commit`] executes the bookkeeping: slack burns,
//!   telemetry updates.
//! - Cache policy: COST-WEIGHTED eviction (recompute-cost ×
//!   hit-probability, both measured) with pins untouchable; a cache
//!   full of pinned nodes surfaces a STRUCTURED error, never OOM.
//! - SKIP-YIELD telemetry per op — risk R4's early-warning dashboard:
//!   the ops with the worst yield are where bound-tightening effort
//!   goes.

use crate::invalidate::{Edge, InvalidateError, InvalidationPlan, Verdict, apply_plan, plan};
use crate::{Store, StoreError, json_string};
use fs_ledger::ContentHash;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// A first-class perturbation plan (inspect, cost, then commit).
#[derive(Debug, Clone)]
pub struct PerturbPlan {
    /// The underlying invalidation plan (frontier + verdicts + rows).
    inner: InvalidationPlan,
    /// Estimated cost of executing the recompute frontier (Σ measured
    /// per-node costs; unknown costs count as `default_cost`).
    estimated_cost: f64,
    /// The baseline a bit-level (hash) invalidator would pay: every
    /// touched node recomputes.
    hash_memo_cost: f64,
    /// The verified-color skip certificates (one row per absorption).
    certificates: Vec<String>,
}

impl PerturbPlan {
    /// Read-only invalidation plan carrying the certified verdicts.
    #[must_use]
    pub fn invalidation(&self) -> &InvalidationPlan {
        &self.inner
    }

    /// Estimated recompute-frontier cost.
    #[must_use]
    pub fn estimated_cost(&self) -> f64 {
        self.estimated_cost
    }

    /// Cost a bit-level hash invalidator would pay.
    #[must_use]
    pub fn hash_memo_cost(&self) -> f64 {
        self.hash_memo_cost
    }

    /// Read-only verified skip-certificate rows.
    #[must_use]
    pub fn certificates(&self) -> &[String] {
        &self.certificates
    }

    /// Nodes the plan recomputes.
    #[must_use]
    pub fn recompute_count(&self) -> usize {
        self.inner
            .verdicts()
            .iter()
            .filter(|(_, v)| matches!(v, Verdict::Recompute { .. }))
            .count()
    }

    /// Nodes the plan skips with certificates.
    #[must_use]
    pub fn skip_count(&self) -> usize {
        self.inner.verdicts().len() - self.recompute_count()
    }
}

/// Per-op skip-yield telemetry (the R4 dashboard).
#[derive(Debug, Clone, Default)]
pub struct SkipYield {
    per_op: BTreeMap<String, (u64, u64)>, // (skips, touches)
}

impl SkipYield {
    /// Yield for one op (`None` before any touches).
    #[must_use]
    pub fn of(&self, op: &str) -> Option<f64> {
        self.per_op
            .get(op)
            .map(|&(s, t)| s as f64 / t.max(1) as f64)
    }

    /// The dashboard rows (op, skips, touches, yield), canonically
    /// ordered.
    #[must_use]
    pub fn dashboard_json(&self) -> String {
        let mut s = String::from("[");
        for (i, (op, &(sk, t))) in self.per_op.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(
                s,
                "{{\"op\":{},\"skips\":{sk},\"touches\":{t},\
                 \"yield\":{:.4}}}",
                json_string(op),
                sk as f64 / t.max(1) as f64
            );
        }
        s.push(']');
        s
    }

    /// Ops sorted worst-yield-first (where tightening effort goes).
    #[must_use]
    pub fn worst_first(&self) -> Vec<(String, f64)> {
        let mut v: Vec<(String, f64)> = self
            .per_op
            .iter()
            .map(|(op, &(s, t))| (op.clone(), s as f64 / t.max(1) as f64))
            .collect();
        v.sort_by(|a, b| {
            a.1.partial_cmp(&b.1)
                .expect("finite yields")
                .then(a.0.cmp(&b.0))
        });
        v
    }
}

/// The operator-facing surface over a store + sensitivity DAG.
#[derive(Debug)]
pub struct RecomputeApi {
    /// The underlying store (public: lookups, pinning).
    pub store: Store,
    edges: Vec<Edge>,
    costs: BTreeMap<[u8; 32], f64>,
    hits: BTreeMap<[u8; 32], u64>,
    plans_seen: u64,
    default_cost: f64,
    yield_telemetry: SkipYield,
}

impl RecomputeApi {
    /// Wrap a store and its sensitivity edges.
    #[must_use]
    pub fn new(store: Store, edges: Vec<Edge>, default_cost: f64) -> Self {
        RecomputeApi {
            store,
            edges,
            costs: BTreeMap::new(),
            hits: BTreeMap::new(),
            plans_seen: 0,
            default_cost,
            yield_telemetry: SkipYield::default(),
        }
    }

    /// Record a node's MEASURED recompute cost (telemetry).
    pub fn record_cost(&mut self, node: &ContentHash, cost: f64) {
        self.costs.insert(*node.as_bytes(), cost);
    }

    /// The skip-yield dashboard.
    #[must_use]
    pub fn skip_yield(&self) -> &SkipYield {
        &self.yield_telemetry
    }

    fn cost_of(&self, node: &ContentHash) -> f64 {
        self.costs
            .get(node.as_bytes())
            .copied()
            .unwrap_or(self.default_cost)
    }

    /// Plan a perturbation: the minimal frontier, its estimated cost,
    /// and the certificates for everything skipped. Pure — nothing
    /// burns until [`RecomputeApi::commit`].
    ///
    /// # Errors
    /// [`InvalidateError`] teaching errors.
    pub fn perturb(
        &mut self,
        node: &ContentHash,
        delta: f64,
    ) -> Result<PerturbPlan, InvalidateError> {
        let inner = plan(&self.store, &self.edges, &[(*node, delta)])?;
        let mut estimated_cost = 0.0;
        let mut hash_memo_cost = 0.0;
        let mut certificates = Vec::new();
        for (h, v) in inner.verdicts() {
            hash_memo_cost += self.cost_of(h);
            match v {
                Verdict::Recompute { .. } => estimated_cost += self.cost_of(h),
                Verdict::Skip { .. } => {}
            }
        }
        for row in inner.rows() {
            if row.contains("\"verdict\":\"skip\"") {
                certificates.push(row.clone());
            }
        }
        Ok(PerturbPlan {
            inner,
            estimated_cost,
            hash_memo_cost,
            certificates,
        })
    }

    /// Commit a plan: burn absorbed slack, update hit and skip-yield
    /// telemetry.
    ///
    /// # Errors
    /// [`StoreError`] if the store changed under the plan.
    pub fn commit(&mut self, p: &PerturbPlan) -> Result<(), StoreError> {
        apply_plan(&mut self.store, &p.inner)?;
        self.plans_seen += 1;
        for (h, v) in p.inner.verdicts() {
            let op = self
                .store
                .get(h)
                .map(|n| n.record.op_id.clone())
                .unwrap_or_default();
            let e = self.yield_telemetry.per_op.entry(op).or_insert((0, 0));
            e.1 += 1;
            if matches!(v, Verdict::Skip { .. }) {
                e.0 += 1;
                *self.hits.entry(*h.as_bytes()).or_insert(0) += 1;
            }
        }
        Ok(())
    }

    /// COST-WEIGHTED eviction: score = recompute-cost ×
    /// hit-probability; evict lowest-score unpinned nodes until at
    /// most `keep_unpinned` remain. A cache whose pinned population
    /// alone exceeds `max_total` is a STRUCTURED refusal.
    ///
    /// # Errors
    /// [`StoreError::CacheFullOfPins`].
    pub fn ensure_capacity(&mut self, max_total: usize) -> Result<u32, StoreError> {
        let pinned = self
            .store
            .iter()
            .filter(|(_, n)| !n.pins.is_empty())
            .count();
        if pinned > max_total {
            return Err(StoreError::CacheFullOfPins {
                pinned,
                capacity: max_total,
            });
        }
        let keep_unpinned = max_total - pinned;
        let plans = self.plans_seen.max(1) as f64;
        let mut scored: Vec<([u8; 32], f64, u64)> = self
            .store
            .iter()
            .filter(|(_, n)| n.pins.is_empty())
            .map(|(k, n)| {
                let cost = self.costs.get(&k).copied().unwrap_or(self.default_cost);
                let hit_prob = self.hits.get(&k).copied().unwrap_or(0) as f64 / plans;
                (k, cost * hit_prob, n.seq)
            })
            .collect();
        // Lowest value first; deterministic tie-break by seq (oldest).
        scored.sort_by(|a, b| {
            // total_cmp is a total order: a non-finite recorded cost can make a
            // score NaN (inf * 0 hits), which would panic partial_cmp().expect().
            // For finite scores this matches the previous ordering exactly.
            a.1.total_cmp(&b.1).then(a.2.cmp(&b.2))
        });
        let excess = scored.len().saturating_sub(keep_unpinned);
        let mut evicted = 0;
        for &(k, _, _) in scored.iter().take(excess) {
            self.store.remove_by_key(k);
            evicted += 1;
        }
        Ok(evicted)
    }
}
