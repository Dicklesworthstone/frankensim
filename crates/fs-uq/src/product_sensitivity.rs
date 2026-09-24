//! Direct-model Sobol sensitivity with the Jansen pick-freeze estimators.
//!
//! A row evaluates A, B, then A with coordinate i replaced by B_i, for each i.
//! Thus N rows cost exactly N*(d+2) physical model calls. A and B reuse the
//! product Monte Carlo sampler at ordinals 2*r and 2*r+1; no surrogate is fit.
//! See Jansen (1999), doi:10.1016/S0010-4655(98)00154-4.

use core::fmt::Display;
use crate::product_execution::sample_parameters;
use crate::product_plan::admit_plan;
use crate::{CorrelationModel, UncertaintyKind, UqPlan, UqStatus};

/// Estimated contribution of one independently distributed physical input.
#[derive(Debug, Clone, PartialEq)]
pub struct SobolEffect {
    /// Exact plan parameter name, in declaration order.
    pub parameter: String,
    /// Declared physical input unit; dimensional applicability is caller-owned.
    pub unit: String,
    /// Main effect: 1 - mean((Y_B-Y_ABi)^2)/(2*V).
    pub first_order: f64,
    /// Main effect plus all interactions: mean((Y_A-Y_ABi)^2)/(2*V).
    pub total_order: f64,
}

/// A complete fixed-design variance decomposition estimate, not a certificate.
/// Finite-sample indices can lie outside [0,1]; they are NEVER clipped. These
/// are not derivatives, causal effects, confidence intervals or stopping rules.
#[derive(Debug, Clone, PartialEq)]
pub struct SobolEstimate {
    /// All input effects. Total-order effects generally do not sum to one.
    pub effects: Vec<SobolEffect>,
    /// Sample standard deviation of the 2*N base observations (not hybrids).
    /// None when the dimensional value exceeds f64, even if ratios are finite.
    pub base_output_std_dev: Option<f64>,
}

/// Accounting and optional estimates for a bounded sensitivity execution.
#[derive(Debug, Clone, PartialEq)]
pub struct SobolReport {
    /// Physical quantity of interest from the plan.
    pub qoi_name: String,
    /// Predeclared number of complete pick-freeze rows.
    pub base_samples: usize,
    /// Original lifetime budget, counting base AND hybrid evaluations.
    pub evaluations_planned: usize,
    /// Completed or terminally rejected calls; unfinished calls do not count.
    pub evaluations_attempted: usize,
    /// Finite accepted observations, including a partly completed row.
    pub evaluations_accepted: usize,
    /// Fully completed rows, retained for accounting even on refusal.
    pub completed_rows: usize,
    /// Only available after the ENTIRE fixed design completed successfully.
    pub estimate: Option<SobolEstimate>,
    /// Why normalization is unavailable after a successful completed design.
    /// Zero sampled variance does not establish that the physical model is constant.
    pub unavailable_reason: Option<&'static str>,
    /// Completion, per-call budget exhaustion, cancellation, or terminal refusal.
    pub status: UqStatus,
    /// Terminal sample/model failure. Failed samples are never silently skipped.
    pub rejection_reason: Option<String>,
}

/// Direct-model sensitivity under explicit independent Gaussian/uniform inputs.
///
/// `plan.budget_max_samples` is the TOTAL solver-call budget, not the row count:
/// it must be divisible by d+2, with at least two rows. The existing product
/// admission bounds it at one million calls and 256 inputs. Parameters listed
/// here must have nonzero uncertainty; fixed values belong in the base model.
/// Correlated, joint-Gaussian and unstated dependence refuse rather than
/// pretending that a coordinate swap preserves their probability law.
///
/// The same deterministic evaluator must be used on every call. Clone retains
/// an in-memory checkpoint including a partly evaluated row. There is no
/// durable format or statistical optional-stopping rule. Long model calls must
/// implement their own cancellation through `advance_interruptible`.
#[derive(Debug, Clone)]
pub struct SobolExecution {
    plan: UqPlan,
    rows: usize,
    width: usize,
    values: Vec<f64>,
    attempted: usize,
    status: UqStatus,
    failure: Option<String>,
    cached_row: Option<usize>,
    a: Vec<f64>,
    b: Vec<f64>,
    hybrid: Vec<f64>,
}

impl SobolExecution {
    /// Admit the probability law, exact work layout and observation storage.
    /// No model callback is invoked, and no alternate sampling method is used.
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let _ = admit_plan(plan)?;
        if !matches!(&plan.correlation, CorrelationModel::Independent) {
            return Err("Sobol physical-input sensitivity requires explicit independence");
        }
        for parameter in &plan.parameters {
            let varies = match parameter.kind {
                UncertaintyKind::AleatoryGaussian { std_dev, .. } => std_dev > 0.0,
                UncertaintyKind::AleatoryUniform { lo, hi } => lo < hi,
                _ => false,
            };
            if !varies {
                return Err("Sobol inputs must have nonzero Gaussian or uniform uncertainty; keep constants in the base model");
            }
        }
        let width = plan.parameters.len() + 2;
        if plan.budget_max_samples % width != 0 || plan.budget_max_samples / width < 2 {
            return Err("Sobol total evaluation budget must equal N*(parameters+2), with N >= 2");
        }
        let mut values = Vec::new();
        values.try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve the admitted Sobol observation budget")?;
        Ok(Self {
            plan: plan.clone(), rows: plan.budget_max_samples / width, width,
            values, attempted: 0, status: UqStatus::BudgetTruncated, failure: None,
            cached_row: None, a: Vec::new(), b: Vec::new(),
            hybrid: vec![0.0; plan.parameters.len()],
        })
    }

    /// The immutable probability model and original total evaluation budget.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan { &self.plan }

    /// Accepted QoIs in row-major A, B, AB_0, ..., AB_(d-1) order.
    #[must_use]
    pub fn observations(&self) -> &[f64] { &self.values }

    /// Read current accounting. Only full fixed designs produce estimates,
    /// so a value-dependent interruption cannot silently select the sample set.
    #[must_use]
    pub fn report(&self) -> SobolReport {
        let (estimate, unavailable_reason) = if self.status == UqStatus::Complete {
            match self.estimate() {
                Ok(estimate) => (Some(estimate), None),
                Err(reason) => (None, Some(reason)),
            }
        } else { (None, None) };
        SobolReport {
            qoi_name: self.plan.target_qoi.clone(), base_samples: self.rows,
            evaluations_planned: self.plan.budget_max_samples,
            evaluations_attempted: self.attempted, evaluations_accepted: self.values.len(),
            completed_rows: self.values.len() / self.width,
            estimate, unavailable_reason, status: self.status,
            rejection_reason: self.failure.clone(),
        }
    }

    /// Spend at most `allowance` additional physical evaluations. Poll before
    /// every call, not only before a row. Zero allowance and terminal runs are
    /// no-ops. Accepted observations survive cancellation and in-memory resume.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> SobolReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// `Ok(None)` means a model call interrupted WITHOUT an observed result.
    /// Its exact input is retried on resume. Errors and nonfinite QoIs refuse
    /// permanently; they cannot be used to filter the design. Panics propagate.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, mut cancelled: C, mut evaluate: F,
    ) -> SobolReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let count = allowance.min(self.plan.budget_max_samples - self.values.len());
        for _ in 0..count {
            if cancelled() { self.status = UqStatus::Cancelled; return self.report(); }
            let parameters = match self.next_parameters() {
                Ok(parameters) => parameters,
                Err(reason) => return self.fail(reason.into()),
            };
            let result = evaluate(parameters);
            if matches!(&result, Ok(None)) {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            self.attempted += 1;
            match result {
                Ok(Some(value)) if value.is_finite() => self.values.push(value),
                Ok(Some(_)) => return self.fail("Sobol model returned a non-finite QoI".into()),
                Err(error) => return self.fail(format!(
                    "Sobol model evaluation {} refused: {error}", self.values.len()
                )),
                Ok(None) => unreachable!("interruption does not consume an ordinal"),
            }
        }
        self.status = if self.values.len() == self.plan.budget_max_samples {
            UqStatus::Complete
        } else { UqStatus::BudgetTruncated };
        self.report()
    }

    fn fail(&mut self, reason: String) -> SobolReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn next_parameters(&mut self) -> Result<&[f64], &'static str> {
        let row = self.values.len() / self.width;
        if self.cached_row != Some(row) {
            // 2*row+1 < total budget <= 1_000_000, inside the sampler's u32 cap.
            self.a = sample_parameters(&self.plan, None, 2 * row)?;
            self.b = sample_parameters(&self.plan, None, 2 * row + 1)?;
            self.cached_row = Some(row);
        }
        match self.values.len() % self.width {
            0 => Ok(&self.a),
            1 => Ok(&self.b),
            slot => {
                self.hybrid.copy_from_slice(&self.a);
                self.hybrid[slot - 2] = self.b[slot - 2];
                Ok(&self.hybrid)
            }
        }
    }

    fn estimate(&self) -> Result<SobolEstimate, &'static str> {
        let normalize = Normalization::new(&self.values);
        let mut center = Sum::default();
        for row in self.values.chunks_exact(self.width) {
            center.add(normalize.value(row[0]));
            center.add(normalize.value(row[1]));
        }
        let mean = center.value / (2 * self.rows) as f64;
        let mut squares = Sum::default();
        for row in self.values.chunks_exact(self.width) {
            squares.add((normalize.value(row[0]) - mean).powi(2));
            squares.add((normalize.value(row[1]) - mean).powi(2));
        }
        // Only the independent A/B observations estimate the output variance.
        let variance = squares.value / (2 * self.rows - 1) as f64;
        if !(variance > 0.0 && variance.is_finite()) {
            return Err("sampled base-output variance is zero or unrepresentable; Sobol indices are undefined");
        }
        let mut effects = Vec::with_capacity(self.plan.parameters.len());
        for (index, parameter) in self.plan.parameters.iter().enumerate() {
            let mut total = Sum::default();
            let mut complement = Sum::default();
            for row in self.values.chunks_exact(self.width) {
                let hybrid = normalize.value(row[index + 2]);
                total.add((normalize.value(row[0]) - hybrid).powi(2));
                complement.add((normalize.value(row[1]) - hybrid).powi(2));
            }
            let total_order = (total.value / (2 * self.rows) as f64) / variance;
            let first_order = 1.0 - (complement.value / (2 * self.rows) as f64) / variance;
            if !total_order.is_finite() || !first_order.is_finite() {
                return Err("sampled Sobol variance ratios exceed finite f64 range");
            }
            effects.push(SobolEffect {
                parameter: parameter.name.clone(), unit: parameter.unit.clone(),
                first_order, total_order,
            });
        }
        let std_dev = fs_math::det::sqrt(variance) * normalize.scale;
        Ok(SobolEstimate { effects, base_output_std_dev: std_dev.is_finite().then_some(std_dev) })
    }
}

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
}

// Center BEFORE scaling when representable, retaining small variations on a
// large common offset. Opposite f64 extremes instead need scaling first. Both
// operations cancel from the indices; no dimensional square is ever required.
struct Normalization { origin: f64, scale: f64, scale_first: bool }
impl Normalization {
    fn new(values: &[f64]) -> Self {
        let origin = values[0];
        let centered_scale = values.iter().fold(0.0_f64, |s, &x| s.max((x - origin).abs()));
        let scale_first = !centered_scale.is_finite();
        let scale = if scale_first {
            values.iter().fold(0.0_f64, |s, x| s.max(x.abs()))
        } else { centered_scale };
        Self { origin, scale: if scale == 0.0 { 1.0 } else { scale }, scale_first }
    }
    fn value(&self, x: f64) -> f64 {
        if self.scale_first { x / self.scale - self.origin / self.scale }
        else { (x - self.origin) / self.scale }
    }
}
