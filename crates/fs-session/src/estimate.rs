//! `estimate()` — the DRY RUN: predicted wall, memory, and energy from
//! the learned cost models WITHOUT EXECUTING, so agents plan before they
//! spend. Every estimate can later be scored against actuals; the
//! calibration report is the cost models' own report card, ledgerable as
//! an artifact.

use fs_ir::{Node, NodeKind};
use fs_plan::CostModel;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::sync::Mutex;

/// Assumed compute power for the energy estimate (J/s per core).
const WATTS_PER_CORE: f64 = 45.0;

/// A dry-run prediction.
#[derive(Debug, Clone, PartialEq)]
pub struct Estimate {
    /// Optimistic wall (sum of per-op p10), seconds.
    pub wall_p10_s: f64,
    /// Median wall, seconds.
    pub wall_p50_s: f64,
    /// Conservative wall (sum of per-op p90), seconds.
    pub wall_p90_s: f64,
    /// Declared memory ask in bytes (from the study's clauses), if any.
    pub mem_ask_bytes: Option<u64>,
    /// Energy estimate in joules (p50 wall × cores × W/core).
    pub energy_j: f64,
    /// Ops that had no cost model (their wall is NOT included) — an
    /// honest coverage statement, never silent.
    pub unmodeled_ops: Vec<String>,
}

fn size_of_call(items: &[Node]) -> f64 {
    for pair in items.windows(2) {
        if let NodeKind::Keyword(k) = &pair[0].kind
            && (k == "dof" || k == "size" || k == "modes")
        {
            match &pair[1].kind {
                NodeKind::Int(i) => {
                    #[allow(clippy::cast_precision_loss)]
                    return *i as f64;
                }
                NodeKind::Float(f) => return *f,
                _ => {}
            }
        }
    }
    1.0
}

fn walk_calls(node: &Node, out: &mut Vec<(String, f64)>) {
    if let NodeKind::List(items) = &node.kind {
        if let Some(h) = node.head()
            && h.contains('.')
        {
            out.push((h.to_string(), size_of_call(items)));
        }
        for child in items {
            walk_calls(child, out);
        }
    }
}

fn mem_ask(node: &Node) -> Option<u64> {
    if let NodeKind::List(items) = &node.kind {
        if node.head() == Some("mem")
            && let Some(NodeKind::Count { value, unit }) = items.get(1).map(|n| &n.kind)
        {
            let factor: f64 = match unit {
                fs_ir::CountUnit::B => 1.0,
                fs_ir::CountUnit::KiB => 1024.0,
                fs_ir::CountUnit::MiB => 1024.0 * 1024.0,
                fs_ir::CountUnit::GiB => 1024.0 * 1024.0 * 1024.0,
                fs_ir::CountUnit::Cores => return None,
            };
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            return Some((value * factor) as u64);
        }
        for child in items {
            if let Some(m) = mem_ask(child) {
                return Some(m);
            }
        }
    }
    None
}

/// Predict a study's cost without executing it.
#[must_use]
pub fn estimate(study: &Node, models: &BTreeMap<String, CostModel>, cores: f64) -> Estimate {
    let mut calls = Vec::new();
    walk_calls(study, &mut calls);
    let (mut p10, mut p50, mut p90) = (0.0f64, 0.0f64, 0.0f64);
    let mut unmodeled = Vec::new();
    for (verb, size) in &calls {
        match models.get(verb).and_then(|m| m.predict(*size).ok()) {
            Some(p) => {
                p10 += p.p10;
                p50 += p.p50;
                p90 += p.p90;
            }
            None => unmodeled.push(verb.clone()),
        }
    }
    unmodeled.sort_unstable();
    unmodeled.dedup();
    Estimate {
        wall_p10_s: p10,
        wall_p50_s: p50,
        wall_p90_s: p90,
        mem_ask_bytes: mem_ask(study),
        energy_j: p50 * cores * WATTS_PER_CORE,
        unmodeled_ops: unmodeled,
    }
}

fn quantiles_of(rows: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    let mut ratios: Vec<f64> = rows
        .iter()
        .filter(|(p, _)| *p > 0.0)
        .map(|(p, a)| a / p)
        .collect();
    if ratios.is_empty() {
        return None;
    }
    ratios.sort_by(f64::total_cmp);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let q = |f: f64| ratios[((ratios.len() - 1) as f64 * f).round() as usize];
    Some((q(0.1), q(0.5), q(0.9)))
}

/// Estimate-vs-actual tracking: the calibration curve the acceptance
/// criteria demand (`actual / predicted-p50` ratio quantiles).
#[derive(Debug, Default)]
pub struct CalibrationReport {
    rows: Mutex<Vec<(f64, f64)>>, // (predicted_p50, actual)
}

impl CalibrationReport {
    /// An empty report.
    #[must_use]
    pub fn new() -> Self {
        CalibrationReport::default()
    }

    /// Record one completed study's actual wall against its estimate.
    pub fn record(&self, estimate: &Estimate, actual_wall_s: f64) {
        self.rows
            .lock()
            .expect("calibration lock")
            .push((estimate.wall_p50_s, actual_wall_s));
    }

    /// Ratio quantiles `(p10, p50, p90)` of actual/predicted; None until
    /// at least one row exists or predictions were all zero.
    #[must_use]
    pub fn ratio_quantiles(&self) -> Option<(f64, f64, f64)> {
        quantiles_of(&self.rows.lock().expect("calibration lock"))
    }

    /// Canonical JSON rendering (the ledger artifact payload).
    #[must_use]
    pub fn to_json(&self) -> String {
        // One lock scope for rows AND quantiles: std mutexes are not
        // reentrant (a nested ratio_quantiles() call here self-deadlocks —
        // caught by the hung conformance run).
        let rows = self.rows.lock().expect("calibration lock");
        let mut out = String::from("{\"kind\":\"estimate-calibration\",\"rows\":[");
        for (i, (p, a)) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "[{p},{a}]");
        }
        out.push_str("],\"ratio_quantiles\":");
        match quantiles_of(&rows) {
            Some((a, b, c)) => {
                let _ = write!(out, "[{a},{b},{c}]");
            }
            None => out.push_str("null"),
        }
        out.push('}');
        out
    }

    /// Persist the calibration table as a content-addressed artifact.
    ///
    /// # Errors
    /// [`crate::SessionError::Persistence`] wrapping the ledger error.
    pub fn flush_to_ledger(
        &self,
        ledger: &fs_ledger::Ledger,
    ) -> Result<fs_ledger::ContentHash, crate::SessionError> {
        let receipt = ledger
            .put_artifact("estimate-calibration", self.to_json().as_bytes(), None)
            .map_err(|e| crate::SessionError::Persistence {
                what: format!("calibration artifact: {e}"),
            })?;
        Ok(receipt.hash)
    }
}
