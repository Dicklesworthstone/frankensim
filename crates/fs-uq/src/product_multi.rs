//! Paired physical responses from ONE model call per Monte Carlo input.
//!
//! The scalar executor owns admission, sampling, cancellation and paid-call
//! accounting. This adapter validates the whole response before accepting any
//! column. Marginals reuse its summary; joint compliance uses the paired rows,
//! never a product of marginal probabilities.
//! Bead: `frankensim-extreal-program-f85xj.6.7` (multi-response propagation).

use core::fmt::Display;
use fs_evidence::Color;
use crate::product_qmc::summarize_replicates;
use crate::{QmcEstimate, UqExecution, UqPlan, UqResult, UqStatus};

/// One physical output, in the immutable callback response order.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiDefinition {
    pub name: String,
    /// Physical output unit; dimensional correctness remains a caller obligation.
    pub unit: String,
    /// Optional inclusive upper limit. None means this output is unconstrained.
    pub upper_limit: Option<f64>,
}
impl QoiDefinition {
    #[must_use]
    pub fn new(name: impl Into<String>, unit: impl Into<String>) -> Self {
        Self { name: name.into(), unit: unit.into(), upper_limit: None }
    }
    #[must_use]
    pub fn with_upper_limit(mut self, limit: f64) -> Self {
        self.upper_limit = Some(limit);
        self
    }
}

/// Canonical scalar statistics retaining the output's identity and units.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiReport {
    pub quantity: QoiDefinition,
    pub statistics: UqResult,
}

/// Paired-output Monte Carlo diagnostics; no independence between outputs.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiQoiReport {
    pub study_name: String,
    /// Original number of physical calls, NOT calls multiplied by output count.
    pub samples_planned: usize,
    /// Completed or terminally rejected calls. Unfinished calls do not count.
    pub samples_evaluated: usize,
    /// Complete finite response vectors, including after a terminal refusal.
    pub samples_accepted: usize,
    /// Declaration order. Empty before observations or after ANY refusal.
    pub outputs: Vec<QoiReport>,
    /// Probability that ALL declared limits hold on the SAME physical response.
    /// Absent without limits/observations or after refusal. Error is descriptive,
    /// not a confidence bound; zero observed failures do not imply zero risk.
    pub joint_compliance: Option<QmcEstimate>,
    pub joint_compliant_samples: Option<usize>,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
    /// Estimated only. There is no scalar dispersion for mixed physical units.
    pub evidence_color: Color,
}

/// One physical trajectory can return temperature, stress, pressure, etc.
///
/// The input UqPlan names the STUDY and its input law/call budget. Its scalar
/// compliance_threshold must be None; limits belong to the output definitions.
/// Input laws, draw order and refusal semantics are exactly those of UqExecution.
/// Supports 1..=32 outputs and at most 4,194,304 retained scalar observations.
/// Observation capacity is reserved before physics. Callback-owned storage and
/// solver work are not bounded by that allowance.
///
/// The deterministic model must retain its meaning and output order. A malformed
/// row or failure refuses the entire experiment; there is no per-output missing
/// data, censoring or replacement. Completed rows survive chunks, cancellation
/// and cloning. Long callbacks own their cancellation. Reports reuse synchronous
/// bounded scalar summaries, including their sorting and numerical refusals.
/// There is no durable codec, multivariate quantile or statistical stopping rule.
#[derive(Debug, Clone)]
pub struct MultiQoiExecution {
    plan: UqPlan,
    quantities: Vec<QoiDefinition>,
    primary: UqExecution,
    additional: Vec<Vec<f64>>,
}
impl MultiQoiExecution {
    pub fn new(plan: &UqPlan, quantities: &[QoiDefinition]) -> Result<Self, &'static str> {
        if quantities.is_empty() || quantities.len() > 32 {
            return Err("multi-QoI execution requires 1..=32 output definitions");
        }
        if plan.target_qoi.is_empty() || plan.compliance_threshold.is_some() {
            return Err("name the study and declare compliance limits on its individual outputs");
        }
        if plan.budget_max_samples.checked_mul(quantities.len())
            .is_none_or(|n| n > 4_194_304)
        {
            return Err("multi-QoI retained observations exceed the storage envelope");
        }
        let mut names = std::collections::BTreeSet::new();
        for q in quantities {
            if q.name.is_empty() || q.unit.is_empty() || !names.insert(&q.name)
                || q.upper_limit.is_some_and(|limit| !limit.is_finite())
            {
                return Err("outputs require unique names, explicit units and finite optional limits");
            }
        }
        let mut scalar = plan.clone();
        scalar.target_qoi = quantities[0].name.clone();
        scalar.compliance_threshold = quantities[0].upper_limit;
        let mut primary = UqExecution::new(&scalar)?;
        primary.values.try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve primary output storage")?;
        let mut additional = Vec::with_capacity(quantities.len() - 1);
        for _ in 1..quantities.len() {
            let mut values = Vec::new();
            values.try_reserve_exact(plan.budget_max_samples)
                .map_err(|_| "cannot reserve paired output storage")?;
            additional.push(values);
        }
        Ok(Self { plan: plan.clone(), quantities: quantities.to_vec(), primary, additional })
    }

    #[must_use]
    pub const fn plan(&self) -> &UqPlan { &self.plan }
    #[must_use]
    pub fn quantities(&self) -> &[QoiDefinition] { &self.quantities }
    #[must_use]
    pub fn evaluations_attempted(&self) -> usize { self.primary.evaluations_attempted() }
    #[must_use]
    pub fn samples_accepted(&self) -> usize { self.primary.observations().len() }
    /// Raw paired column, in original sample order, even after refusal.
    #[must_use]
    pub fn observations(&self, output: usize) -> Option<&[f64]> {
        if output == 0 { Some(self.primary.observations()) }
        else { self.additional.get(output - 1).map(Vec::as_slice) }
    }
    /// Reconstruct one accepted physical response without evaluating the model.
    #[must_use]
    pub fn observation(&self, ordinal: usize) -> Option<Vec<f64>> {
        if ordinal >= self.samples_accepted() { return None; }
        Some((0..self.quantities.len()).map(|j| self.column(j)[ordinal]).collect())
    }

    fn column(&self, output: usize) -> &[f64] {
        if output == 0 { self.primary.observations() } else { &self.additional[output - 1] }
    }
    fn has_limits(&self) -> bool { self.quantities.iter().any(|q| q.upper_limit.is_some()) }
    fn complies(&self, ordinal: usize) -> bool {
        self.quantities.iter().enumerate().all(|(j, q)| {
            q.upper_limit.is_none_or(|limit| self.column(j)[ordinal] <= limit)
        })
    }

    /// Read-only marginal and joint summaries, with the scalar numerical policy.
    /// Temporary report-only views clone one additional column at a time; they
    /// never sample or factor the input law. Canonical quantiles also use scratch.
    #[must_use]
    pub fn report(&self) -> MultiQoiReport {
        let mut result = MultiQoiReport {
            study_name: self.plan.target_qoi.clone(), samples_planned: self.plan.budget_max_samples,
            samples_evaluated: self.evaluations_attempted(), samples_accepted: self.samples_accepted(),
            outputs: Vec::new(), joint_compliance: None, joint_compliant_samples: None,
            status: self.primary.status, rejection_reason: self.primary.failure.clone(),
            evidence_color: Color::Estimated {
                estimator: "paired physical Monte Carlo; no confidence or model-error bound".into(),
                dispersion: f64::INFINITY,
            },
        };
        if self.primary.status == UqStatus::Refused || result.samples_accepted == 0 { return result; }
        for (j, quantity) in self.quantities.iter().enumerate() {
            let statistics = if j == 0 { self.primary.report() } else {
                let mut plan = self.plan.clone();
                plan.target_qoi = quantity.name.clone();
                plan.compliance_threshold = quantity.upper_limit;
                // This internal view is used ONLY for the existing report algebra.
                // The primary owner alone samples and retains the admitted factor.
                UqExecution {
                    plan, factor: None, values: self.column(j).to_vec(),
                    attempted: self.evaluations_attempted(), status: self.primary.status,
                    failure: None,
                }.report()
            };
            if statistics.status == UqStatus::Refused {
                result.status = UqStatus::Refused;
                result.rejection_reason = Some(format!("output {}: {}", quantity.name,
                    statistics.rejection_reason.as_deref().unwrap_or("statistics refused")));
                result.outputs.clear();
                return result;
            }
            result.outputs.push(QoiReport { quantity: quantity.clone(), statistics });
        }
        if self.has_limits() {
            let indicators: Vec<_> = (0..result.samples_accepted)
                .map(|i| f64::from(u8::from(self.complies(i)))).collect();
            result.joint_compliant_samples = Some(indicators.iter().filter(|&&x| x == 1.0).count());
            result.joint_compliance = summarize_replicates(&indicators);
        }
        result
    }

    /// At most allowance additional physical calls, independently of output count.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> MultiQoiReport
    where F: FnMut(&[f64]) -> Result<Vec<f64>, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// Ok(None) is unfinished work and retries the exact input on resume. Only
    /// a complete finite vector is accepted. A bad component invalidates the
    /// entire draw and is terminal; it cannot silently shorten a marginal sample.
    /// Callback panics propagate before any columns change; side effects must
    /// tolerate unfinished-call retry. Zero allowance/terminal calls are no-ops.
    pub fn advance_interruptible<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F)
        -> MultiQoiReport
    where F: FnMut(&[f64]) -> Result<Option<Vec<f64>>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.primary.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let dimension = self.quantities.len();
        let additional = &mut self.additional;
        self.primary.advance_interruptible(allowance, cancelled, |x| -> Result<Option<f64>, String> {
            let Some(row) = evaluate(x).map_err(|e| e.to_string())? else { return Ok(None); };
            if row.len() != dimension {
                return Err(format!("expected {dimension} outputs, received {}", row.len()));
            }
            if let Some(j) = row.iter().position(|v| !v.is_finite()) {
                return Err(format!("output {j} is non-finite"));
            }
            // All columns were reserved at admission. Publish additional values
            // only after complete validation; the scalar owner accepts row[0].
            for (column, &value) in additional.iter_mut().zip(&row[1..]) { column.push(value); }
            Ok(Some(row[0]))
        });
        let report = self.report();
        if report.status == UqStatus::Refused {
            self.primary.status = UqStatus::Refused;
            self.primary.failure = report.rejection_reason.clone();
        }
        report
    }
}
