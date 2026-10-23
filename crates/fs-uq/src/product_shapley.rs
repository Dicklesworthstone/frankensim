//! Physical-input Shapley variance attribution under a declared Gaussian law.
//!
//! For each coalition S, two samples share X_S and have independent conditional
//! complements. Half their mean squared QoI difference estimates
//! c(S) = E[Var(Y | X_S)]. Subset-weighted decreases of c give Shapley effects.
//! Dependence and interactions both contribute; these are not causal effects.
//!
//! Beads: `frankensim-ext-uq-reliability-upgrades-pwui`,
//! `frankensim-extreal-program-f85xj.6.7`.

mod checkpoint;

/// Concurrent conditional-pair evaluations with original-ordinal recovery.
pub mod batch;
pub use batch::{ShapleyBatchExecution, ShapleyBatchReport};

use core::fmt::Display;
use fs_evidence::Color;

use crate::product_execution::sample_parameters;
use crate::product_plan::admit_plan;
use crate::product_sensitivity::{Normalization, Sum};
use crate::{CorrelationModel, UncertaintyKind, UqPlan, UqStatus};

/// Estimated variance attribution to one declared physical input.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapleyEffect {
    /// Exact physical parameter name in the original declaration order.
    pub parameter: String,
    /// Original physical unit, not a latent-normal unit.
    pub unit: String,
    /// Subset-weighted conditional variance reduction divided by output variance.
    /// Finite-sample values can leave [0,1] and are never clipped.
    pub index: f64,
}

/// Complete fixed-design sensitivity estimate, without a confidence guarantee.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapleyEstimate {
    /// Physical-input effects. Their algebraic sum is one, up to rounding.
    pub effects: Vec<ShapleyEffect>,
    /// c(S)/c(empty), indexed by the physical-coordinate subset bitmask.
    /// Entry zero is 1; the full-subset entry is 0. Sampling noise can violate
    /// monotonicity between other subsets; no isotonic repair is applied.
    pub conditional_variance_fractions: Vec<f64>,
    /// sqrt(c(empty)) in physical output units, absent if unrepresentable.
    /// The dimensionless effects can remain available when this overflows.
    pub base_output_std_dev: Option<f64>,
}

/// Fixed work accounting and optional complete-design attribution.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapleyReport {
    /// Quantity of interest from the immutable plan.
    pub qoi_name: String,
    /// Independent base pairs; each row covers every non-full coalition.
    pub base_samples: usize,
    /// Original lifetime budget, N * 2 * (2^d - 1) physical calls.
    pub evaluations_planned: usize,
    /// Completed or terminally rejected calls; unfinished calls do not count.
    pub evaluations_attempted: usize,
    /// Finite accepted observations, including partially completed pairs/rows.
    pub evaluations_accepted: usize,
    /// Complete rows retained for accounting, even when the execution refuses.
    pub completed_rows: usize,
    /// Available only after the full fixed allocation and defined normalization.
    pub estimate: Option<ShapleyEstimate>,
    /// Why a completed design could not define finite variance attribution.
    pub unavailable_reason: Option<&'static str>,
    /// Always Estimated, with no finite total-error or confidence bound.
    pub evidence_color: Color,
    /// Complete, paused by allowance/cancellation, or terminally refused.
    pub status: UqStatus,
    /// Original terminal model/sampling diagnostic. Samples are never replaced.
    pub rejection_reason: Option<String>,
}

#[derive(Debug, Clone)]
struct Coalition {
    // Conditioned physical indices first, then the remaining indices; each
    // group stays in declaration order. The vector maps factor row -> physical.
    order: Vec<usize>,
    conditioned: usize,
    factor: Option<Vec<Vec<f64>>>,
}

/// Direct physical Shapley effects for 1..=8 nonconstant Gaussian parameters.
///
/// The plan must request Monte Carlo and declare Independent or JointGaussian
/// dependence. Every coalition's permuted law passes the existing numerical PSD
/// policy before any model call. Singular laws are supported without covariance
/// inversion. Zero-variance inputs must remain fixed in the caller's model:
/// conditioning on a constant physical input must not reveal a hidden germ.
///
/// Each row costs exactly 2*(2^d-1) physical calls: one conditional pair for
/// each subset except the full set. The empty subset supplies independent base
/// observations; c(full)=0 follows from the deterministic evaluator contract.
/// The original total budget must contain at least two whole rows and retains
/// the shared one-million-call cap. Effects require the complete allocation.
///
/// Physical inputs always use original parameter order and units. Long solvers
/// must implement their own cancellation; this owner polls between every call.
/// Clone retains paid work, including one completed half of a conditional pair.
/// `checkpoint`/`restore` retain that work across process lifetimes under the
/// same model and deterministic runtime profile. Reports use bounded synchronous
/// reductions. No adaptive subset selection, causal interpretation or anytime
/// stopping is supplied.
#[derive(Debug, Clone)]
pub struct ShapleyExecution {
    plan: UqPlan,
    germ_plan: UqPlan,
    coalitions: Vec<Coalition>,
    rows: usize,
    width: usize,
    values: Vec<f64>,
    attempted: usize,
    status: UqStatus,
    failure: Option<String>,
    cached_row: Option<usize>,
    a: Vec<f64>,
    b: Vec<f64>,
    parameters: Vec<f64>,
}

impl ShapleyExecution {
    /// Admit physical laws, all permuted PSD factors, exact call allocation and
    /// bounded observation storage before any callback. No uniform/unknown law
    /// is silently replaced, and a generic correlation matrix is insufficient.
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let _ = admit_plan(plan)?;
        let dim = plan.parameters.len();
        if dim > 8 {
            return Err("physical Shapley execution supports 1..=8 Gaussian inputs");
        }
        if !matches!(
            plan.correlation,
            CorrelationModel::Independent | CorrelationModel::JointGaussian { .. }
        ) {
            return Err(
                "physical Shapley sensitivity requires explicit Independent or JointGaussian dependence",
            );
        }
        if plan.parameters.iter().any(|p| {
            !matches!(p.kind,
            UncertaintyKind::AleatoryGaussian { std_dev, .. } if std_dev > 0.0)
        }) {
            return Err(
                "physical Shapley inputs need nonzero Gaussian uncertainty; keep constants in the physical model",
            );
        }
        let subset_count = 1usize << dim;
        let width = 2 * (subset_count - 1);
        if plan.budget_max_samples % width != 0 || plan.budget_max_samples / width < 2 {
            return Err("Shapley budget must equal N*2*(2^d-1), with N>=2");
        }
        let mut coalitions = Vec::with_capacity(subset_count - 1);
        for mask in 0..subset_count - 1 {
            let mut order: Vec<_> = (0..dim).filter(|&i| mask & (1 << i) != 0).collect();
            let conditioned = order.len();
            order.extend((0..dim).filter(|&i| mask & (1 << i) == 0));
            let mut permuted = plan.clone();
            permuted.parameters = order.iter().map(|&i| plan.parameters[i].clone()).collect();
            if let CorrelationModel::JointGaussian { matrix } = &plan.correlation {
                permuted.correlation = CorrelationModel::JointGaussian {
                    matrix: order
                        .iter()
                        .map(|&i| order.iter().map(|&j| matrix[i][j]).collect())
                        .collect(),
                };
            }
            // Factor EACH order under the existing numerical policy, including
            // its singular-pivot consistency checks. A rejected permutation is
            // not repaired, jittered or replaced by an independence assumption.
            let factor = admit_plan(&permuted)?;
            coalitions.push(Coalition {
                order,
                conditioned,
                factor,
            });
        }
        let mut germ_plan = plan.clone();
        germ_plan.correlation = CorrelationModel::Independent;
        for parameter in &mut germ_plan.parameters {
            parameter.kind = UncertaintyKind::AleatoryGaussian {
                mean: 0.0,
                std_dev: 1.0,
            };
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve admitted Shapley observation storage")?;
        Ok(Self {
            plan: plan.clone(),
            germ_plan,
            coalitions,
            rows: plan.budget_max_samples / width,
            width,
            values,
            attempted: 0,
            status: UqStatus::BudgetTruncated,
            failure: None,
            cached_row: None,
            a: Vec::new(),
            b: Vec::new(),
            parameters: vec![0.0; dim],
        })
    }

    /// Immutable physical law and original total solver-call budget.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan {
        &self.plan
    }

    /// Accepted QoIs, row-major, then subset bitmask, then first/second sample.
    /// All bitmasks refer to original physical parameter declaration order.
    #[must_use]
    pub fn observations(&self) -> &[f64] {
        &self.values
    }

    /// Current paid work, with effects only after a complete fixed design.
    #[must_use]
    pub fn report(&self) -> ShapleyReport {
        let (estimate, unavailable_reason) = if self.status == UqStatus::Complete {
            match self.estimate() {
                Ok(estimate) => (Some(estimate), None),
                Err(reason) => (None, Some(reason)),
            }
        } else {
            (None, None)
        };
        ShapleyReport {
            qoi_name: self.plan.target_qoi.clone(), base_samples: self.rows,
            evaluations_planned: self.plan.budget_max_samples,
            evaluations_attempted: self.attempted, evaluations_accepted: self.values.len(),
            completed_rows: self.values.len() / self.width, estimate, unavailable_reason,
            evidence_color: Color::Estimated {
                estimator: "fixed-design physical Gaussian Shapley attribution; no causal, confidence or model-error bound".into(),
                dispersion: f64::INFINITY,
            }, status: self.status, rejection_reason: self.failure.clone(),
        }
    }

    /// Spend at most this many additional physical calls. Zero allowance and
    /// completed/refused executions are no-ops; paid observations are retained.
    pub fn advance<F, E, C>(
        &mut self,
        allowance: usize,
        cancelled: C,
        mut evaluate: F,
    ) -> ShapleyReport
    where
        F: FnMut(&[f64]) -> Result<f64, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// `Ok(None)` pauses an unfinished call and retries identical inputs on
    /// resume. The solver owns interrupted costs. Errors and nonfinite outputs
    /// permanently refuse; no sample may be censored or replaced. Callback
    /// panics propagate, and side effects must tolerate unfinished-call retry.
    pub fn advance_interruptible<F, E, C>(
        &mut self,
        allowance: usize,
        mut cancelled: C,
        mut evaluate: F,
    ) -> ShapleyReport
    where
        F: FnMut(&[f64]) -> Result<Option<f64>, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let count = allowance.min(self.plan.budget_max_samples - self.values.len());
        for _ in 0..count {
            if cancelled() {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            let parameters = match self.parameters_at(self.values.len()) {
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
                Ok(Some(_)) => {
                    return self.fail("Shapley physical model returned a non-finite QoI".into());
                }
                Err(error) => {
                    return self.fail(format!(
                        "Shapley physical evaluation {} refused: {error}",
                        self.values.len()
                    ));
                }
                Ok(None) => unreachable!("interruption did not consume an ordinal"),
            }
        }
        self.status = if self.values.len() == self.plan.budget_max_samples {
            UqStatus::Complete
        } else {
            UqStatus::BudgetTruncated
        };
        self.report()
    }

    fn fail(&mut self, reason: String) -> ShapleyReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn parameters_at(&mut self, ordinal: usize) -> Result<&[f64], &'static str> {
        let row = ordinal / self.width;
        if self.cached_row != Some(row) {
            // 2*row+1 < original budget <=1_000_000, within u32 tile addressing.
            self.a = sample_parameters(&self.germ_plan, None, 2 * row)?;
            self.b = sample_parameters(&self.germ_plan, None, 2 * row + 1)?;
            self.cached_row = Some(row);
        }
        let slot = ordinal % self.width;
        let coalition = &self.coalitions[slot / 2];
        let second = slot % 2 == 1;
        let germ = |j: usize| {
            if second && j >= coalition.conditioned {
                self.b[j]
            } else {
                self.a[j]
            }
        };
        for (i, &physical) in coalition.order.iter().enumerate() {
            let z = coalition
                .factor
                .as_ref()
                .map_or_else(|| germ(i), |l| (0..=i).map(|j| l[i][j] * germ(j)).sum());
            let UncertaintyKind::AleatoryGaussian { mean, std_dev } =
                self.plan.parameters[physical].kind
            else {
                unreachable!("admission requires Gaussian physical marginals")
            };
            let value = mean + std_dev * z;
            if !value.is_finite() {
                return Err("Shapley conditional physical sample overflowed");
            }
            self.parameters[physical] = value;
        }
        Ok(&self.parameters)
    }

    fn estimate(&self) -> Result<ShapleyEstimate, &'static str> {
        let normalize = Normalization::new(&self.values);
        let mut fractions = Vec::with_capacity(self.coalitions.len() + 1);
        for subset in 0..self.coalitions.len() {
            let mut sum = Sum::default();
            for row in self.values.chunks_exact(self.width) {
                let difference =
                    normalize.value(row[2 * subset]) - normalize.value(row[2 * subset + 1]);
                sum.add(difference * difference);
            }
            fractions.push(sum.value / (2 * self.rows) as f64);
        }
        let base_variance = fractions[0];
        if !(base_variance.is_finite() && base_variance > 0.0) {
            return Err(
                "sampled base variance is zero or unrepresentable; physical Shapley effects are undefined",
            );
        }
        for fraction in &mut fractions {
            *fraction /= base_variance;
            if !fraction.is_finite() {
                return Err("Shapley conditional variance ratio is not representable");
            }
        }
        fractions.push(0.0); // A deterministic model has Var(Y | all physical inputs)=0.
        let dim = self.plan.parameters.len();
        let mut effects = Vec::with_capacity(dim);
        for (i, parameter) in self.plan.parameters.iter().enumerate() {
            let mut effect = Sum::default();
            for subset in 0..self.coalitions.len() {
                if subset & (1 << i) == 0 {
                    let size = subset.count_ones() as usize;
                    // |S|! (d-|S|-1)! / d! = 1 / (d * choose(d-1, |S|)).
                    // d<=8 bounds all integer products without factorial overflow.
                    let combinations = (1..=size).fold(1usize, |c, k| c * (dim - k) / k);
                    let weight = 1.0 / (dim * combinations) as f64;
                    effect.add(weight * (fractions[subset] - fractions[subset | (1 << i)]));
                }
            }
            if !effect.value.is_finite() {
                return Err("Shapley effect is not representable");
            }
            effects.push(ShapleyEffect {
                parameter: parameter.name.clone(),
                unit: parameter.unit.clone(),
                index: effect.value,
            });
        }
        let std_dev = fs_math::det::sqrt(base_variance) * normalize.scale;
        Ok(ShapleyEstimate {
            effects,
            conditional_variance_fractions: fractions,
            base_output_std_dev: std_dev.is_finite().then_some(std_dev),
        })
    }
}
