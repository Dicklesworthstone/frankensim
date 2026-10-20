//! Product-plan multilevel Monte Carlo with owned coarse/fine coupling.
//!
//! Each correction uses one physical parameter vector at both resolutions.
//! Different corrections use disjoint Philox sample ordinals. Fixed level
//! allocations target the finest declared discretization; the sampling error
//! does not bound its bias or the physical model error.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`.

use core::fmt::Display;
use fs_evidence::Color;

use crate::product_execution::sample_parameters;
use crate::product_plan::admit_plan;
use crate::product_qmc::summarize_replicates;
use crate::{PropagationMethod, UqPlan, UqStatus};

/// Mean estimate with descriptive sampling dispersion. Missing dispersion is
/// unavailable, including a singleton sample; it never means zero uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct MlmcEstimate {
    /// Sum of level-correction means, or one correction mean in a level report.
    pub mean: f64,
    /// Estimated variance of the mean, absent when unrepresentable or undefined.
    pub estimator_variance: Option<f64>,
    /// Square root of estimator variance, retained when its square overflows.
    pub standard_error: Option<f64>,
}

/// Work and observed correction statistics for one declared resolution level.
#[derive(Debug, Clone, PartialEq)]
pub struct MlmcLevelReport {
    /// Resolution passed to the evaluator; zero is the coarsest model.
    pub level: usize,
    /// Original number of complete corrections required on this level.
    pub samples_planned: usize,
    /// Accepted complete corrections; an unfinished pair is excluded.
    pub samples_accepted: usize,
    /// Statistics of P0, or of P(level) - P(level-1). Absent on failure.
    pub correction: Option<MlmcEstimate>,
    /// Statistics of the corresponding difference of compliance indicators.
    pub compliance_correction: Option<MlmcEstimate>,
}

/// Current multilevel work and estimates. Global estimates appear only after
/// the complete predeclared allocation; partial levels remain diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct MlmcExecutionReport {
    /// Physical quantity of interest from the immutable input plan.
    pub qoi_name: String,
    /// Original physical-call budget: N0 + 2 * sum(Nl, l > 0).
    pub model_calls_planned: usize,
    /// Completed or terminally rejected physical calls, including half-pairs.
    /// An interrupted call has no result and is retried without consuming a call.
    pub model_calls_evaluated: usize,
    /// Complete accepted level corrections across the ladder.
    pub samples_accepted: usize,
    /// A paid fine result is waiting for its coarse partner.
    pub pending_pair: bool,
    /// Level-major diagnostics in the original declaration order.
    pub levels: Vec<MlmcLevelReport>,
    /// Finest-declared-level expectation and sampling error estimate.
    pub estimate: Option<MlmcEstimate>,
    /// Finest-level compliance probability estimate. Finite-sample telescoping
    /// may leave [0,1]; the estimator is never clipped to manufacture a bound.
    pub compliance: Option<MlmcEstimate>,
    /// Always Estimated; unavailable dispersion is infinite, never zero.
    pub evidence_color: Color,
    /// Complete, paused by call allowance/cancellation, or terminally refused.
    pub status: UqStatus,
    /// Original numerical/model refusal; no failed sample is skipped.
    pub rejection_reason: Option<String>,
}

#[derive(Debug, Clone)]
struct LevelSamples {
    values: Vec<f64>,
    compliance: Vec<f64>,
}

#[derive(Debug, Clone)]
struct PendingPair {
    parameters: Vec<f64>,
    fine: f64,
}

/// Bounded, resumable MLMC on 1..=64 fixed resolution levels, with at least two
/// samples per level. The plan budget counts physical calls, not corrections.
///
/// The callback must keep the same deterministic physical meaning and units at
/// every level and throughout resumption. This owner couples input parameters;
/// it cannot couple hidden solver randomness or certify a discretization ladder.
/// Clone preserves accepted corrections and a paid half-pair. No durable
/// serialization, adaptive level selection, bias estimate or stopping bound is
/// supplied. Cancellation is polled between model calls; long solvers must also
/// poll their own scope, and bounded report reductions are synchronous.
#[derive(Debug, Clone)]
pub struct MlmcExecution {
    plan: UqPlan,
    allocation: Vec<usize>,
    factor: Option<Vec<Vec<f64>>>,
    levels: Vec<LevelSamples>,
    level: usize,
    accepted: usize,
    attempted: usize,
    pending: Option<PendingPair>,
    status: UqStatus,
    failure: Option<String>,
}

impl MlmcExecution {
    /// Admit the fixed allocation and full probability law before model calls.
    /// The exact call count must equal `plan.budget_max_samples`, bounded by
    /// shared product admission at one million. Marginal/dependence admission
    /// and parameter sampling are the same owners used by Monte Carlo.
    pub fn new(plan: &UqPlan, samples_per_level: &[usize]) -> Result<Self, &'static str> {
        if plan.method != PropagationMethod::MultilevelMonteCarlo {
            return Err("MLMC execution requires a MultilevelMonteCarlo plan");
        }
        if samples_per_level.is_empty()
            || samples_per_level.len() > 64
            || samples_per_level.iter().any(|&n| n < 2)
        {
            return Err("MLMC needs 1..=64 fixed levels with at least two samples each");
        }
        let calls = samples_per_level
            .iter()
            .enumerate()
            .try_fold(0usize, |sum, (l, &n)| {
                n.checked_mul(if l == 0 { 1 } else { 2 })
                    .and_then(|work| sum.checked_add(work))
            })
            .ok_or("MLMC physical-call allocation overflows")?;
        if calls != plan.budget_max_samples {
            return Err("MLMC allocation must match the original physical-call budget exactly");
        }
        let mut admission = plan.clone();
        admission.method = PropagationMethod::MonteCarlo;
        let factor = admit_plan(&admission)?;
        let mut levels = Vec::new();
        levels
            .try_reserve_exact(samples_per_level.len())
            .map_err(|_| "MLMC level storage allocation refused")?;
        for &count in samples_per_level {
            let mut values = Vec::new();
            values
                .try_reserve_exact(count)
                .map_err(|_| "MLMC correction storage allocation refused")?;
            let mut compliance = Vec::new();
            if plan.compliance_threshold.is_some() {
                compliance
                    .try_reserve_exact(count)
                    .map_err(|_| "MLMC compliance storage allocation refused")?;
            }
            levels.push(LevelSamples { values, compliance });
        }
        Ok(Self {
            plan: plan.clone(),
            allocation: samples_per_level.to_vec(),
            factor,
            levels,
            level: 0,
            accepted: 0,
            attempted: 0,
            pending: None,
            status: UqStatus::BudgetTruncated,
            failure: None,
        })
    }

    /// Immutable physical plan; the stored method remains MLMC.
    #[must_use]
    pub fn plan(&self) -> &UqPlan {
        &self.plan
    }

    /// Immutable level-major allocation of complete corrections.
    #[must_use]
    pub fn allocation(&self) -> &[usize] {
        &self.allocation
    }

    /// Completed or terminally rejected physical model calls.
    #[must_use]
    pub const fn evaluations_attempted(&self) -> usize {
        self.attempted
    }

    /// Accepted correction values, never a pooled physical-QoI distribution.
    #[must_use]
    pub fn corrections(&self, level: usize) -> Option<&[f64]> {
        self.levels
            .get(level)
            .map(|samples| samples.values.as_slice())
    }

    /// Accepted differences of physical compliance indicators, if requested.
    #[must_use]
    pub fn compliance_corrections(&self, level: usize) -> Option<&[f64]> {
        self.plan.compliance_threshold?;
        self.levels
            .get(level)
            .map(|samples| samples.compliance.as_slice())
    }

    /// Summarize retained work without invoking the physical model. A singleton
    /// level has a mean but undefined dispersion. Errors hide all statistics;
    /// paid-work counts remain available. Cost is O(accepted corrections).
    #[must_use]
    pub fn report(&self) -> MlmcExecutionReport {
        let mut levels: Vec<_> = self
            .levels
            .iter()
            .enumerate()
            .map(|(level, samples)| MlmcLevelReport {
                level,
                samples_planned: self.allocation[level],
                samples_accepted: samples.values.len(),
                correction: if self.failure.is_none() {
                    summarize(&samples.values)
                } else {
                    None
                },
                compliance_correction: if self.failure.is_none() {
                    summarize(&samples.compliance)
                } else {
                    None
                },
            })
            .collect();
        let mut status = self.status;
        let mut rejection_reason = self.failure.clone();
        let mut estimate = None;
        let mut compliance = None;
        if status == UqStatus::Complete && rejection_reason.is_none() {
            estimate = combine(&levels, false);
            if self.plan.compliance_threshold.is_some() {
                compliance = combine(&levels, true);
            }
            if estimate.is_none()
                || (self.plan.compliance_threshold.is_some() && compliance.is_none())
            {
                status = UqStatus::Refused;
                rejection_reason = Some("MLMC telescoping mean is not representable".into());
                estimate = None;
                compliance = None;
                for level in &mut levels {
                    level.correction = None;
                    level.compliance_correction = None;
                }
            }
        }
        let dispersion = estimate
            .as_ref()
            .and_then(|e| e.standard_error)
            .unwrap_or(f64::INFINITY);
        let label = if rejection_reason.is_some() {
            "refused; no multilevel estimate"
        } else if status != UqStatus::Complete {
            "partial MLMC level diagnostics; no finest-level estimate"
        } else {
            "fixed-level MLMC sampling estimate; no bias, continuum, model-error or stopping bound"
        };
        MlmcExecutionReport {
            qoi_name: self.plan.target_qoi.clone(),
            model_calls_planned: self.plan.budget_max_samples,
            model_calls_evaluated: self.attempted,
            samples_accepted: self.accepted,
            pending_pair: self.pending.is_some(),
            levels,
            estimate,
            compliance,
            evidence_color: Color::Estimated {
                estimator: label.into(),
                dispersion,
            },
            status,
            rejection_reason,
        }
    }

    /// Execute at most `allowance` additional physical calls. Fine and coarse
    /// calls count separately. A zero allowance is a no-op, even if cancelled.
    pub fn advance<F, E, C>(
        &mut self,
        allowance: usize,
        cancelled: C,
        mut evaluate: F,
    ) -> MlmcExecutionReport
    where
        F: FnMut(usize, &[f64]) -> Result<f64, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        self.advance_interruptible(allowance, cancelled, |level, x| {
            evaluate(level, x).map(Some)
        })
    }

    /// `Ok(None)` is unfinished work: stop and retry the exact physical call.
    /// A completed fine result survives interruption of its coarse partner.
    /// `Err` and nonfinite outputs permanently refuse the execution; they cannot
    /// censor an inconvenient sample. Panics are not caught, and callback side
    /// effects must tolerate retry. Interrupted work costs belong to the solver.
    pub fn advance_interruptible<F, E, C>(
        &mut self,
        allowance: usize,
        mut cancelled: C,
        mut evaluate: F,
    ) -> MlmcExecutionReport
    where
        F: FnMut(usize, &[f64]) -> Result<Option<f64>, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let work = allowance.min(self.plan.budget_max_samples - self.attempted);
        for _ in 0..work {
            if cancelled() {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            let fresh = if self.pending.is_none() {
                match sample_parameters(&self.plan, self.factor.as_deref(), self.accepted) {
                    Ok(parameters) => Some(parameters),
                    Err(reason) => return self.fail(reason.into()),
                }
            } else {
                None
            };
            let parameters = fresh.as_ref().unwrap_or_else(|| {
                &self
                    .pending
                    .as_ref()
                    .expect("pending pair owns its parameters")
                    .parameters
            });
            let model_level = if self.pending.is_some() {
                self.level - 1
            } else {
                self.level
            };
            let result = evaluate(model_level, parameters);
            if matches!(&result, Ok(None)) {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            self.attempted += 1;
            let value = match result {
                Ok(Some(value)) if value.is_finite() => value,
                Ok(Some(_)) => {
                    return self.fail("MLMC physical model returned a non-finite QoI".into());
                }
                Ok(None) => unreachable!("interruption returned without consuming a call"),
                Err(error) => {
                    return self.fail(format!(
                        "MLMC level {model_level}, correction {} refused: {error}",
                        self.accepted
                    ));
                }
            };
            if self.level > 0 && self.pending.is_none() {
                self.pending = Some(PendingPair {
                    parameters: fresh.expect("fine call has fresh parameters"),
                    fine: value,
                });
            } else {
                let (fine, coarse) = match self.pending.take() {
                    Some(pair) => (pair.fine, Some(value)),
                    None => (value, None),
                };
                let correction = fine - coarse.unwrap_or(0.0);
                if !correction.is_finite() {
                    return self.fail(
                        "MLMC finite physical outputs have an unrepresentable correction".into(),
                    );
                }
                self.levels[self.level].values.push(correction);
                if let Some(threshold) = self.plan.compliance_threshold {
                    let indicator = f64::from(u8::from(fine <= threshold))
                        - coarse.map_or(0.0, |v| f64::from(u8::from(v <= threshold)));
                    self.levels[self.level].compliance.push(indicator);
                }
                self.accepted += 1;
                if self.levels[self.level].values.len() == self.allocation[self.level] {
                    self.level += 1;
                }
            }
        }
        self.status = if self.level == self.levels.len() {
            UqStatus::Complete
        } else {
            UqStatus::BudgetTruncated
        };
        let report = self.report();
        if report.status == UqStatus::Refused {
            self.status = UqStatus::Refused;
            self.failure = report.rejection_reason.clone();
        }
        report
    }

    fn fail(&mut self, reason: String) -> MlmcExecutionReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }
}

// Reuse the ordered scaled moments of independent observations. The QMC owner
// applies this same algebra to independent replicate means; here the inputs
// are independent complete corrections within one resolution level.
fn summarize(values: &[f64]) -> Option<MlmcEstimate> {
    let moments = summarize_replicates(values)?;
    Some(MlmcEstimate {
        mean: moments.mean,
        standard_error: moments.standard_error,
        estimator_variance: moments.standard_error.and_then(representable_square),
    })
}

fn representable_square(value: f64) -> Option<f64> {
    let square = value * value;
    (square.is_finite() && (square > 0.0 || value == 0.0)).then_some(square)
}

fn combine(levels: &[MlmcLevelReport], compliance: bool) -> Option<MlmcEstimate> {
    let estimates: Option<Vec<_>> = levels
        .iter()
        .map(|level| {
            if compliance {
                level.compliance_correction.as_ref()
            } else {
                level.correction.as_ref()
            }
        })
        .collect();
    let estimates = estimates?;
    let scale = estimates.iter().fold(0.0_f64, |s, e| s.max(e.mean.abs()));
    let mean = if scale == 0.0 {
        0.0
    } else {
        estimates.iter().map(|e| e.mean / scale).sum::<f64>() * scale
    };
    if !mean.is_finite() {
        return None;
    }
    let errors: Option<Vec<_>> = estimates.iter().map(|e| e.standard_error).collect();
    let standard_error = errors.and_then(|errors| {
        let scale = errors.iter().copied().fold(0.0_f64, f64::max);
        let norm = if scale == 0.0 {
            0.0
        } else {
            fs_math::det::sqrt(errors.iter().map(|e| (e / scale).powi(2)).sum::<f64>()) * scale
        };
        norm.is_finite().then_some(norm)
    });
    Some(MlmcEstimate {
        mean,
        standard_error,
        estimator_variance: standard_error.and_then(representable_square),
    })
}
