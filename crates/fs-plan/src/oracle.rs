//! Consumer wiring (plan §11.4 "consumers wired"): fitted cost models
//! behind fs-geom's `CostOracle` so the Rep Router plans with THIS
//! machine's measured history, and a fs-ledger `tune`-table loader so
//! models are rebuilt deterministically from ledger snapshots.

use std::collections::BTreeMap;

use crate::cost::{CostModel, CostObservation};
use fs_ledger::Ledger;

/// A [`fs_geom::CostOracle`] backed by per-edge quantile cost models.
/// Each edge is registered with the reference problem size its routing
/// requests are quoted at; recorded actuals feed the online refits.
#[derive(Debug, Default)]
pub struct PlanCostOracle {
    models: BTreeMap<String, (f64, CostModel)>,
    errors: BTreeMap<String, Vec<f64>>,
}

impl PlanCostOracle {
    /// An empty oracle.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register an edge with its reference size (idempotent; keeps
    /// existing observations).
    pub fn register_edge(&mut self, edge: &str, reference_size: f64) {
        self.models
            .entry(edge.to_string())
            .or_insert_with(|| (reference_size, CostModel::new()))
            .0 = reference_size;
    }

    /// The fitted model for an edge, if any.
    #[must_use]
    pub fn model(&self, edge: &str) -> Option<&CostModel> {
        self.models.get(edge).map(|(_, m)| m)
    }
}

impl fs_geom::CostOracle for PlanCostOracle {
    fn measured_cost_s(&self, edge: &str) -> Option<f64> {
        let (size, model) = self.models.get(edge)?;
        model.predict(*size).ok().map(|p| p.p50)
    }

    fn measured_error_abs(&self, edge: &str) -> Option<f64> {
        // Conservative: the p90 of observed absolute errors.
        let errs = self.errors.get(edge)?;
        if errs.is_empty() {
            return None;
        }
        let mut sorted = errs.clone();
        sorted.sort_by(f64::total_cmp);
        let idx = ((sorted.len() as f64 - 1.0) * 0.9).round() as usize;
        Some(sorted[idx.min(sorted.len() - 1)])
    }

    fn record(&mut self, edge: &str, cost_s: f64, error_abs: f64) {
        let entry =
            self.models.entry(edge.to_string()).or_insert_with(|| (1.0, CostModel::new()));
        let size = entry.0;
        // Nonpositive actuals cannot enter the log-log model; drop them
        // (the router still ran — a zero-cost record carries no signal).
        let _ = entry.1.observe(CostObservation { size, cost_s: cost_s.max(1e-12) });
        if error_abs.is_finite() && error_abs >= 0.0 {
            self.errors.entry(edge.to_string()).or_default().push(error_abs);
        }
    }
}

/// Extract a numeric field from one of our canonical flat JSON lines
/// (fs-roofline/fs-ledger emit `"key":value` pairs; this is a scanner, not
/// a JSON parser — documented boundary, sufficient for our own payloads).
#[must_use]
pub fn json_f64_field(json: &str, key: &str) -> Option<f64> {
    let needle = format!("\"{key}\":");
    let start = json.find(&needle)? + needle.len();
    let rest = &json[start..];
    let end = rest
        .char_indices()
        .find(|(_, c)| !matches!(c, '0'..='9' | '.' | '-' | '+' | 'e' | 'E'))
        .map_or(rest.len(), |(i, _)| i);
    rest[..end].parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Rebuild a per-kernel cost model from the ledger `tune` table: every
/// roofline row for `kernel` contributes one observation (size = the
/// recorded elements-per-run when present, else the reference size).
/// Deterministic from a ledger snapshot (P2).
///
/// # Errors
/// Ledger read errors propagate.
pub fn cost_model_from_tune(
    ledger: &Ledger,
    kernel: &str,
    reference_size: f64,
) -> Result<CostModel, fs_ledger::LedgerError> {
    let mut model = CostModel::new();
    for row in ledger.tune_rows(kernel)? {
        // Roofline rows carry elems_per_sec; cost per reference run is
        // size / rate.
        if let Some(rate) = json_f64_field(&row.measured, "elems_per_sec")
            && rate > 0.0
        {
            let size = json_f64_field(&row.measured, "elements").unwrap_or(reference_size);
            let _ = model.observe(CostObservation { size, cost_s: size / rate });
        }
    }
    Ok(model)
}
