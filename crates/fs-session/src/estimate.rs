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

#[allow(clippy::cast_precision_loss)]
fn integer_size(value: i64) -> f64 {
    value as f64
}

fn size_of_call(items: &[Node], verb: &str) -> Result<f64, crate::SessionError> {
    let mut size = None;
    for (index, item) in items.iter().enumerate() {
        let NodeKind::Keyword(keyword) = &item.kind else {
            continue;
        };
        if keyword != "dof" && keyword != "size" && keyword != "modes" {
            continue;
        }
        if size.is_some() {
            return Err(crate::SessionError::Submission {
                what: format!(
                    "operation {verb:?} declares more than one :dof/:size/:modes feature"
                ),
            });
        }
        let value = items
            .get(index + 1)
            .ok_or_else(|| crate::SessionError::Submission {
                what: format!("operation {verb:?} has no value after :{keyword}"),
            })?;
        size = Some(match &value.kind {
            NodeKind::Int(value) => integer_size(*value),
            NodeKind::Float(value) => *value,
            _ => {
                return Err(crate::SessionError::Submission {
                    what: format!("operation {verb:?} requires a numeric value after :{keyword}"),
                });
            }
        });
    }
    Ok(size.unwrap_or(1.0))
}

fn walk_calls(
    node: &Node,
    models: &BTreeMap<String, CostModel>,
    out: &mut Vec<(String, f64)>,
) -> Result<(), crate::SessionError> {
    if let NodeKind::List(items) = &node.kind {
        if let Some(h) = node.head()
            && (h.contains('.') || models.contains_key(h))
        {
            out.push((h.to_string(), size_of_call(items, h)?));
        }
        for child in items {
            walk_calls(child, models, out)?;
        }
    }
    Ok(())
}

fn mem_ask(budget: Option<&Node>) -> Result<Option<u64>, crate::SessionError> {
    let Some(budget) = budget else {
        return Ok(None);
    };
    let Some(items) = budget.items() else {
        return Err(crate::SessionError::Submission {
            what: "recognized budget pillar is not a list".to_string(),
        });
    };
    let mut ask = None;
    for clause in &items[1..] {
        let Some(values) = clause.items() else {
            return Err(crate::SessionError::Submission {
                what: "budget entries must be parenthesized clauses such as (mem 96GiB)"
                    .to_string(),
            });
        };
        let Some(resource) = clause.head() else {
            return Err(crate::SessionError::Submission {
                what: "budget entries must have a symbolic name and at least one value".to_string(),
            });
        };
        if values.len() < 2 {
            return Err(crate::SessionError::Submission {
                what: format!("budget entry {resource:?} has no value"),
            });
        }
        if resource != "mem" {
            continue;
        }
        if ask.is_some() {
            return Err(crate::SessionError::Submission {
                what: "duplicate memory budgets are ambiguous; retain exactly one (mem ...) clause"
                    .to_string(),
            });
        }
        if values.len() != 2 {
            return Err(crate::SessionError::Submission {
                what: "memory budget must have exactly the shape (mem COUNT)".to_string(),
            });
        }
        let NodeKind::Count { value, unit } = &values[1].kind else {
            return Err(crate::SessionError::Submission {
                what: "memory budget operand must be a byte count such as 96GiB".to_string(),
            });
        };
        let factor: f64 = match unit {
            fs_ir::CountUnit::B => 1.0,
            fs_ir::CountUnit::KiB => 1024.0,
            fs_ir::CountUnit::MiB => 1024.0 * 1024.0,
            fs_ir::CountUnit::GiB => 1024.0 * 1024.0 * 1024.0,
            fs_ir::CountUnit::Cores => {
                return Err(crate::SessionError::InvalidResource {
                    resource: "declared memory ask",
                    value: *value,
                    requirement: "must use a byte unit (B, KiB, MiB, or GiB)",
                });
            }
        };
        let bytes = value * factor;
        // 2^64 is exactly representable as f64 and is the first value that
        // cannot fit u64. `u64::MAX as f64` rounds to that same value, so a
        // cast-based upper-bound check would accidentally accept overflow.
        const U64_EXCLUSIVE_UPPER_BOUND: f64 = 18_446_744_073_709_551_616.0;
        if !value.is_finite()
            || *value <= 0.0
            || !bytes.is_finite()
            || bytes >= U64_EXCLUSIVE_UPPER_BOUND
            || bytes.fract() != 0.0
        {
            return Err(crate::SessionError::InvalidResource {
                resource: "declared memory ask",
                value: bytes,
                requirement: "must be a finite, positive whole-byte quantity below 2^64 bytes",
            });
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bytes = bytes as u64;
        ask = Some(bytes);
    }
    Ok(ask)
}

/// Predict a study's cost without executing it.
#[must_use]
///
/// # Errors
/// [`crate::SessionError::InvalidResource`] when cores, a declared memory ask,
/// or a derived wall/energy estimate is outside its finite non-negative domain;
/// [`crate::SessionError::Submission`] for malformed studies or explicit size
/// features and for cost-model refusals.
pub fn estimate(
    study: &Node,
    models: &BTreeMap<String, CostModel>,
    cores: f64,
) -> Result<Estimate, crate::SessionError> {
    if !cores.is_finite() || cores < 0.0 {
        return Err(crate::SessionError::InvalidResource {
            resource: "estimate cores",
            value: cores,
            requirement: "must be finite and non-negative",
        });
    }
    let recognized =
        fs_ir::Study::from_node(study).map_err(|error| crate::SessionError::Submission {
            what: format!("cannot estimate malformed study: {error}"),
        })?;
    let mem_ask_bytes = mem_ask(recognized.budget)?;
    let mut calls = Vec::new();
    for (_, expression) in &recognized.lets {
        walk_calls(expression, models, &mut calls)?;
    }
    for clause in &recognized.body {
        walk_calls(clause, models, &mut calls)?;
    }
    let (mut p10, mut p50, mut p90) = (0.0f64, 0.0f64, 0.0f64);
    let mut unmodeled = Vec::new();
    for (verb, size) in &calls {
        if !size.is_finite() || *size < 0.0 {
            return Err(crate::SessionError::InvalidResource {
                resource: "estimate operation size",
                value: *size,
                requirement: "must be finite and non-negative",
            });
        }
        let Some(model) = models.get(verb) else {
            unmodeled.push(verb.clone());
            continue;
        };
        let prediction = model
            .predict(*size)
            .map_err(|error| crate::SessionError::Submission {
                what: format!("cost model for operation {verb:?} refused size {size}: {error}"),
            })?;
        for value in [prediction.p10, prediction.p50, prediction.p90] {
            if !value.is_finite() || value < 0.0 {
                return Err(crate::SessionError::InvalidResource {
                    resource: "predicted operation wall-seconds",
                    value,
                    requirement: "must be finite and non-negative",
                });
            }
        }
        if prediction.p10 > prediction.p50 || prediction.p50 > prediction.p90 {
            return Err(crate::SessionError::Submission {
                what: format!("cost model for operation {verb:?} returned reversed wall quantiles"),
            });
        }
        p10 += prediction.p10;
        p50 += prediction.p50;
        p90 += prediction.p90;
    }
    unmodeled.sort_unstable();
    unmodeled.dedup();
    for (resource, value) in [
        ("estimated p10 wall-seconds", p10),
        ("estimated p50 wall-seconds", p50),
        ("estimated p90 wall-seconds", p90),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(crate::SessionError::InvalidResource {
                resource,
                value,
                requirement: "must remain finite and non-negative after model aggregation",
            });
        }
    }
    if p10 > p50 || p50 > p90 {
        return Err(crate::SessionError::Submission {
            what: "aggregated cost-model wall quantiles are reversed".to_string(),
        });
    }
    let energy_j = p50 * cores * WATTS_PER_CORE;
    if !energy_j.is_finite() || energy_j < 0.0 {
        return Err(crate::SessionError::InvalidResource {
            resource: "estimated energy",
            value: energy_j,
            requirement: "must remain finite and non-negative after model aggregation",
        });
    }
    Ok(Estimate {
        wall_p10_s: p10,
        wall_p50_s: p50,
        wall_p90_s: p90,
        mem_ask_bytes,
        energy_j,
        unmodeled_ops: unmodeled,
    })
}

fn zero_summary_of(rows: &[CalRow]) -> ZeroPredictionSummary {
    let zero_rows: Vec<&CalRow> = rows.iter().filter(|r| r.predicted == 0.0).collect();
    ZeroPredictionSummary {
        true_zero: zero_rows.iter().filter(|r| r.fully_modeled).count(),
        unmodeled: zero_rows.iter().filter(|r| !r.fully_modeled).count(),
        actual_quantiles_s: sorted_quantiles(zero_rows.iter().map(|r| r.actual).collect()),
    }
}

/// One recorded estimate-vs-actual observation.
#[derive(Debug, Clone, Copy, PartialEq)]
struct CalRow {
    predicted: f64,
    actual: f64,
    /// False when the scored estimate carried unmodeled ops: its
    /// prediction systematically EXCLUDES wall, so a zero prediction
    /// from such a row is a coverage gap, not a zero-cost observation.
    fully_modeled: bool,
}

fn sorted_quantiles(mut values: Vec<f64>) -> Option<(f64, f64, f64)> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let q = |f: f64| values[((values.len() - 1) as f64 * f).round() as usize];
    Some((q(0.1), q(0.5), q(0.9)))
}

fn quantiles_of(rows: &[CalRow]) -> Option<(f64, f64, f64)> {
    sorted_quantiles(
        rows.iter()
            .filter(|r| r.predicted > 0.0)
            .map(|r| r.actual / r.predicted)
            .collect(),
    )
}

/// What the zero-prediction rows look like (bead gp3.21): the rows the
/// ratio quantiles CANNOT see. Reported instead of silently dropped so
/// repeated zero predictions can never make calibration look healthier
/// than it is. No ratio is invented for them — the actual-time
/// distribution travels raw.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZeroPredictionSummary {
    /// Zero-prediction rows from FULLY MODELED estimates: the model
    /// genuinely asserted zero cost.
    pub true_zero: usize,
    /// Zero-prediction rows from estimates with unmodeled ops: the
    /// prediction excluded wall it could not see (a coverage gap).
    pub unmodeled: usize,
    /// (p10, p50, p90) of the ACTUAL wall over all zero-prediction
    /// rows; `None` when there are none.
    pub actual_quantiles_s: Option<(f64, f64, f64)>,
}

/// Governance-configurable calibration thresholds (bead gp3.21).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrationPolicy {
    /// Largest tolerated fraction of zero-prediction rows before the
    /// report is DEGRADED (quantiles alone are no longer trustworthy).
    pub max_zero_prediction_fraction: f64,
}

impl Default for CalibrationPolicy {
    fn default() -> Self {
        CalibrationPolicy {
            max_zero_prediction_fraction: 0.25,
        }
    }
}

/// The policy verdict over a calibration report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CalibrationHealth {
    /// Zero-prediction mass is within the declared threshold.
    Healthy,
    /// Too much of the evidence is invisible to the ratio quantiles.
    Degraded {
        /// Observed zero-prediction fraction.
        zero_fraction: f64,
        /// The policy threshold it exceeded.
        limit: f64,
    },
}

/// Estimate-vs-actual tracking: the calibration curve the acceptance
/// criteria demand (`actual / predicted-p50` ratio quantiles), plus the
/// zero-prediction telemetry the quantiles cannot carry.
#[derive(Debug, Default)]
pub struct CalibrationReport {
    rows: Mutex<Vec<CalRow>>,
}

impl CalibrationReport {
    /// An empty report.
    #[must_use]
    pub fn new() -> Self {
        CalibrationReport::default()
    }

    /// Record one completed study's actual wall against its estimate.
    ///
    /// # Errors
    /// [`crate::SessionError::InvalidResource`] if either value or their
    /// positive-prediction ratio would poison the calibration JSON with a
    /// negative or non-finite number.
    pub fn record(
        &self,
        estimate: &Estimate,
        actual_wall_s: f64,
    ) -> Result<(), crate::SessionError> {
        let predicted = estimate.wall_p50_s;
        if !predicted.is_finite() || predicted < 0.0 {
            return Err(crate::SessionError::InvalidResource {
                resource: "calibration predicted wall-seconds",
                value: predicted,
                requirement: "must be finite and non-negative",
            });
        }
        if !actual_wall_s.is_finite() || actual_wall_s < 0.0 {
            return Err(crate::SessionError::InvalidResource {
                resource: "calibration actual wall-seconds",
                value: actual_wall_s,
                requirement: "must be finite and non-negative",
            });
        }
        let ratio = if predicted > 0.0 {
            Some(actual_wall_s / predicted)
        } else {
            None
        };
        if let Some(ratio) = ratio
            && !ratio.is_finite()
        {
            return Err(crate::SessionError::InvalidResource {
                resource: "calibration actual/predicted ratio",
                value: ratio,
                requirement: "must be finite for a positive prediction",
            });
        }
        self.rows.lock().expect("calibration lock").push(CalRow {
            predicted: estimate.wall_p50_s,
            actual: actual_wall_s,
            fully_modeled: estimate.unmodeled_ops.is_empty(),
        });
        Ok(())
    }

    /// The zero-prediction telemetry: counts split by whether the
    /// estimate was fully modeled, plus the raw actual-time quantiles.
    #[must_use]
    pub fn zero_prediction_summary(&self) -> ZeroPredictionSummary {
        let rows = self.rows.lock().expect("calibration lock");
        zero_summary_of(&rows)
    }

    /// Judge this report against a governance policy.
    ///
    /// # Errors
    /// [`crate::SessionError::InvalidResource`] for a non-finite or
    /// out-of-[0, 1] threshold — an unusable policy cannot certify.
    pub fn health(
        &self,
        policy: &CalibrationPolicy,
    ) -> Result<CalibrationHealth, crate::SessionError> {
        let limit = policy.max_zero_prediction_fraction;
        if !limit.is_finite() || !(0.0..=1.0).contains(&limit) {
            return Err(crate::SessionError::InvalidResource {
                resource: "calibration zero-prediction threshold",
                value: limit,
                requirement: "must be finite and within [0, 1]",
            });
        }
        let rows = self.rows.lock().expect("calibration lock");
        if rows.is_empty() {
            return Ok(CalibrationHealth::Healthy);
        }
        #[allow(clippy::cast_precision_loss)]
        let zero_fraction =
            rows.iter().filter(|r| r.predicted == 0.0).count() as f64 / rows.len() as f64;
        Ok(if zero_fraction <= limit {
            CalibrationHealth::Healthy
        } else {
            CalibrationHealth::Degraded {
                zero_fraction,
                limit,
            }
        })
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
        for (i, row) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let modeled = u8::from(row.fully_modeled);
            let _ = write!(out, "[{},{},{modeled}]", row.predicted, row.actual);
        }
        out.push_str("],\"ratio_quantiles\":");
        match quantiles_of(&rows) {
            Some((a, b, c)) => {
                let _ = write!(out, "[{a},{b},{c}]");
            }
            None => out.push_str("null"),
        }
        let zero = zero_summary_of(&rows);
        let _ = write!(
            out,
            ",\"zero_predictions\":{{\"true_zero\":{},\"unmodeled\":{},\"actual_quantiles_s\":",
            zero.true_zero, zero.unmodeled
        );
        match zero.actual_quantiles_s {
            Some((a, b, c)) => {
                let _ = write!(out, "[{a},{b},{c}]");
            }
            None => out.push_str("null"),
        }
        out.push_str("}}");
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
