//! Replicated Latin-hypercube propagation for independent declared marginals.
//!
//! Each complete replicate visits every probability stratum exactly once in
//! every parameter. Philox Fisher-Yates permutations and open-cell jitter are
//! keyed by seed, replicate and logical coordinate, never by chunk/worker.
//! Gaussian coordinates use fs-bo's existing inverse-normal approximation.
//!
//! Only COMPLETE replicates contribute estimates. Their means, not the
//! dependent points within a hypercube, supply descriptive standard errors.
//! Finite-grid bias, inverse-normal approximation, model/discretization error,
//! and optional stopping are NOT bounded by these reports. Same-profile replay
//! is deterministic; no cross-ISA evidence is claimed by this module.
//!
//! The existing product checkpoint encoder is reused with a distinct LHS
//! identity binding. Its internal MonteCarlo plan is an admission/storage
//! carrier only: this module never executes or reports Monte Carlo sampling.

pub mod mean_control;

use core::fmt::Display;
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use fs_rand::StreamKey;
use crate::{CorrelationModel, ParameterUncertainty, PropagationMethod, UncertaintyKind,
    UqCheckpointError, UqExecution, UqPlan, UqStatus};

const PERMUTATION_KERNEL: u32 = 0x4c48_0000;
const JITTER_KERNEL: u32 = 0x4c49_0000;
const MAX_DESIGN_CELLS: usize = 4_000_000;
const IDENTITY_DOMAIN: &str = "org.frankensim.uq.latin-hypercube.execution.v1";

/// An explicit probability model and immutable, replicated stratification budget.
#[derive(Debug, Clone, PartialEq)]
pub struct LhsPlan {
    /// Named output of the deterministic evaluator.
    pub target_qoi: String,
    /// Optional empirical compliance event, QoI <= threshold.
    pub compliance_threshold: Option<f64>,
    /// Named, unit-bearing Gaussian or uniform marginal declarations.
    pub parameters: Vec<ParameterUncertainty>,
    /// Must explicitly be Independent. Correlation/coplanar rank rearrangement
    /// is not silently substituted for the declared joint distribution.
    pub correlation: CorrelationModel,
    /// Independently randomized hypercubes, in 2..=256.
    pub replicates: usize,
    /// Strata per coordinate per replicate, at least two; not power-of-two only.
    pub samples_per_replicate: usize,
    /// Original Philox seed; immutable across interruption and resumption.
    pub seed: u64,
}

impl LhsPlan {
    fn admission_plan(&self) -> Result<UqPlan, &'static str> {
        if self.correlation != CorrelationModel::Independent {
            return Err("Latin hypercubes require explicitly independent marginals");
        }
        if !(2..=256).contains(&self.replicates) || self.samples_per_replicate < 2 {
            return Err("Latin hypercubes need 2..=256 replicates with at least two strata each");
        }
        if !(1..=256).contains(&self.parameters.len()) {
            return Err("Latin hypercubes support 1..=256 declared parameters");
        }
        let budget = self.replicates.checked_mul(self.samples_per_replicate)
            .ok_or("Latin-hypercube sample budget overflow")?;
        let cells = self.parameters.len().checked_mul(self.samples_per_replicate)
            .ok_or("Latin-hypercube design allocation overflow")?;
        if cells > MAX_DESIGN_CELLS {
            return Err("Latin-hypercube design exceeds four million u32 cells");
        }
        Ok(UqPlan {
            target_qoi: self.target_qoi.clone(), compliance_threshold: self.compliance_threshold,
            parameters: self.parameters.clone(), correlation: self.correlation.clone(),
            method: PropagationMethod::MonteCarlo, budget_max_samples: budget, seed: self.seed,
        })
    }
}

/// Mean of equally sized complete hypercubes and between-replicate dispersion.
#[derive(Debug, Clone, PartialEq)]
pub struct LhsEstimate {
    /// Average of complete replicate means/proportions.
    pub mean: f64,
    /// Descriptive standard error; absent until two complete replicates exist.
    /// Not an anytime-valid confidence interval or a model-error certificate.
    pub standard_error: Option<f64>,
}

/// A report never promotes empirical quadrature to verified evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct LhsReport {
    /// Physical output name from the admitted plan.
    pub qoi_name: String,
    /// Original lifetime evaluator budget, replicates times strata.
    pub samples_planned: usize,
    /// Completed or terminally failed calls; unfinished interrupted calls do not count.
    pub samples_evaluated: usize,
    /// Accepted observations, including an incomplete hypercube.
    pub samples_accepted: usize,
    /// Number of complete hypercubes retained, even on terminal failure.
    pub completed_replicates: usize,
    /// Complete replicate means in logical order; absent on refusal.
    pub replicate_means: Vec<f64>,
    /// QoI estimate using complete replicates only.
    pub estimate: Option<LhsEstimate>,
    /// Empirical compliance estimate, only when a threshold was declared.
    pub compliance: Option<LhsEstimate>,
    /// Always Estimated, with no within-hypercube iid standard-error claim.
    pub evidence_color: Color,
    /// Complete only after the entire predeclared lifetime budget is consumed.
    pub status: UqStatus,
    /// Original terminal diagnostic, without a substituted or skipped sample.
    pub rejection_reason: Option<String>,
}

/// Bounded resumable execution. The original plan and evaluator meaning must
/// remain fixed. At most four million permutation cells and one million
/// observations are retained. Cancellation polls occur before model calls and
/// every 1024 permutation swaps; expensive callbacks must also poll their Cx.
#[derive(Debug, Clone)]
pub struct LhsExecution {
    plan: LhsPlan,
    storage: UqExecution,
    permutations: Vec<Vec<u32>>,
    current_replicate: Option<usize>,
}

impl LhsExecution {
    /// Admit names, units, finite marginal parameters, explicit independence,
    /// total sample count and design memory before model evaluation/allocation.
    /// Reuses the established product probability admission, not a second policy.
    pub fn new(plan: &LhsPlan) -> Result<Self, &'static str> {
        let storage = UqExecution::new(&plan.admission_plan()?)?;
        Ok(Self { plan: plan.clone(), storage, permutations: Vec::new(), current_replicate: None })
    }

    /// Immutable admitted physical model and stratification layout.
    #[must_use]
    pub fn plan(&self) -> &LhsPlan { &self.plan }
    /// Accepted observations in original ordinal order, including a partial replicate.
    #[must_use]
    pub fn observations(&self) -> &[f64] { &self.storage.values }
    /// Completed/terminally failed evaluations. Interruptions do not consume ordinals.
    #[must_use]
    pub fn evaluations_attempted(&self) -> usize { self.storage.attempted }

    /// Recompute statistics without invoking the physical model. Incomplete
    /// hypercubes remain paid work, but are never presented as a stratified rule.
    #[must_use]
    pub fn report(&self) -> LhsReport {
        let n = self.plan.samples_per_replicate;
        let mut means = Vec::new();
        let mut proportions = Vec::new();
        if self.storage.failure.is_none() {
            for block in self.storage.values.chunks_exact(n) {
                means.push(scaled_mean(block));
                if let Some(limit) = self.plan.compliance_threshold {
                    proportions.push(block.iter().filter(|&&v| v <= limit).count() as f64 / n as f64);
                }
            }
        }
        let estimate = replicate_estimate(&means);
        let compliance = replicate_estimate(&proportions);
        let dispersion = estimate.as_ref().and_then(|e| e.standard_error).unwrap_or(0.0);
        LhsReport {
            qoi_name: self.plan.target_qoi.clone(), samples_planned: self.storage.plan.budget_max_samples,
            samples_evaluated: self.storage.attempted, samples_accepted: self.storage.values.len(),
            completed_replicates: self.storage.values.len() / n, replicate_means: means,
            estimate, compliance, evidence_color: Color::Estimated {
                estimator: "replicated Latin hypercube; complete-replicate descriptive error only; no iid point, stopping, grid-bias or model-error bound".into(), dispersion,
            }, status: self.storage.status, rejection_reason: self.storage.failure.clone(),
        }
    }

    /// Execute at most `allowance` more evaluations under the original lifetime
    /// budget. Completed and refused states are terminal; zero is a no-op.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> LhsReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// An unfinished `Ok(None)` stops and retries the SAME point on resume.
    /// Errors/nonfinite observations refuse permanently: no censoring or dummy
    /// zero is allowed. Callback side effects must tolerate retry; panics propagate.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, mut cancelled: C, mut evaluate: F,
    ) -> LhsReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.storage.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let count = allowance.min(self.storage.plan.budget_max_samples - self.storage.values.len());
        for _ in 0..count {
            if cancelled() { self.storage.status = UqStatus::Cancelled; return self.report(); }
            let ordinal = self.storage.values.len();
            let replicate = ordinal / self.plan.samples_per_replicate;
            if self.current_replicate != Some(replicate) {
                match self.build_design(replicate, &mut cancelled) {
                    Ok(true) => {}
                    Ok(false) => { self.storage.status = UqStatus::Cancelled; return self.report(); }
                    Err(reason) => return self.fail(reason.into()),
                }
            }
            if cancelled() { self.storage.status = UqStatus::Cancelled; return self.report(); }
            let parameters = match self.parameters(ordinal) {
                Ok(x) => x,
                Err(reason) => return self.fail(reason.into()),
            };
            let result = evaluate(&parameters);
            if matches!(&result, Ok(None)) {
                self.storage.status = UqStatus::Cancelled;
                return self.report();
            }
            self.storage.attempted += 1;
            match result {
                Ok(Some(value)) if value.is_finite() => self.storage.values.push(value),
                Ok(Some(_)) => return self.fail("model returned a non-finite QoI".into()),
                Err(error) => return self.fail(format!("model evaluation {ordinal} refused: {error}")),
                Ok(None) => unreachable!("interruption returned without consuming the ordinal"),
            }
        }
        self.storage.status = if self.storage.values.len() == self.storage.plan.budget_max_samples {
            UqStatus::Complete
        } else { UqStatus::BudgetTruncated };
        self.report()
    }

    fn fail(&mut self, reason: String) -> LhsReport {
        self.storage.status = UqStatus::Refused;
        self.storage.failure = Some(reason);
        self.report()
    }

    fn build_design<C: FnMut() -> bool>(&mut self, replicate: usize, cancelled: &mut C) -> Result<bool, &'static str> {
        self.current_replicate = None;
        self.permutations.clear();
        let n = self.plan.samples_per_replicate;
        for dimension in 0..self.plan.parameters.len() {
            if cancelled() { return Ok(false); }
            let mut permutation = Vec::new();
            permutation.try_reserve_exact(n).map_err(|_| "Latin-hypercube permutation allocation refused")?;
            // Shared admission caps the TOTAL budget at 1,000,000, below u32::MAX.
            permutation.extend((0..n).map(|i| i as u32));
            let mut stream = StreamKey { seed: self.plan.seed,
                kernel: PERMUTATION_KERNEL + replicate as u32, tile: dimension as u32 }.stream();
            for i in (1..n).rev() {
                if i % 1024 == 0 && cancelled() { return Ok(false); }
                let j = stream.next_below((i + 1) as u64) as usize;
                permutation.swap(i, j);
            }
            self.permutations.push(permutation);
        }
        self.current_replicate = Some(replicate);
        Ok(true)
    }

    fn parameters(&self, ordinal: usize) -> Result<Vec<f64>, &'static str> {
        let n = self.plan.samples_per_replicate;
        let replicate = ordinal / n;
        let index = ordinal % n;
        if self.current_replicate != Some(replicate) {
            return Err("Latin-hypercube design does not match the sample ordinal");
        }
        let mut jitter = StreamKey { seed: self.plan.seed,
            kernel: JITTER_KERNEL + replicate as u32, tile: index as u32 }.stream();
        let mut parameters = Vec::with_capacity(self.plan.parameters.len());
        for (dimension, parameter) in self.plan.parameters.iter().enumerate() {
            // 32-bit cell MIDPOINTS leave >1 ulp of room at stratum boundaries
            // throughout the admitted n<=500,000 range, including the tails.
            let offset = ((jitter.next_u64() >> 32) as f64 + 0.5) / 4_294_967_296.0;
            let u = (f64::from(self.permutations[dimension][index]) + offset) / n as f64;
            let value = match parameter.kind {
                UncertaintyKind::AleatoryUniform { lo, hi } => {
                    if lo == hi { lo } else { (1.0 - u) * lo + u * hi }
                }
                UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                    if std_dev == 0.0 { mean } else { mean + std_dev * fs_bo::phi_inv(u) }
                }
                _ => return Err("Latin-hypercube marginal has no admitted probability measure"),
            };
            if !value.is_finite() { return Err("Latin-hypercube parameter conversion overflowed"); }
            parameters.push(value);
        }
        Ok(parameters)
    }

    fn checkpoint_identity(&self, model_identity: ContentHash) -> ContentHash {
        let mut bytes = Vec::with_capacity(52);
        bytes.extend_from_slice(&model_identity.0);
        bytes.extend_from_slice(&(self.plan.replicates as u64).to_le_bytes());
        bytes.extend_from_slice(&(self.plan.samples_per_replicate as u64).to_le_bytes());
        bytes.extend_from_slice(&fs_rand::STREAM_SEMANTICS_VERSION.to_le_bytes());
        hash_domain(IDENTITY_DOMAIN, &bytes)
    }

    /// Persist the paid prefix, exact plan, layout and LHS sampler identity.
    /// `model_identity` must cover implementation, all fixed inputs and solver
    /// settings. Checksums detect corruption, not dishonest model observations.
    pub fn checkpoint(&self, model_identity: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        self.storage.checkpoint(self.checkpoint_identity(model_identity))
    }

    /// Restore without re-evaluating retained observations. Rebuild permutations
    /// deterministically at the next advance; changing seed/layout/model refuses.
    pub fn restore(plan: &LhsPlan, model_identity: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError> {
        let mut execution = Self::new(plan).map_err(UqCheckpointError::InvalidPlan)?;
        execution.storage = UqExecution::restore(&execution.storage.plan,
            execution.checkpoint_identity(model_identity), bytes)?;
        Ok(execution)
    }
}

fn scaled_mean(values: &[f64]) -> f64 {
    scaled_mean_interruptible(values, &mut || false).expect("uncancelled finite mean")
}

fn scaled_mean_interruptible(values: &[f64], cancelled: &mut impl FnMut() -> bool)
    -> Result<f64, crate::UqControlError>
{
    use crate::product_control_variate::poll;
    let mut scale = 0.0_f64;
    for (i, value) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 { return Ok(0.0); }
    let mut sum = 0.0;
    let mut correction = 0.0;
    for (i, &value) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        let y = value / scale - correction;
        let next = sum + y;
        correction = (next - sum) - y;
        sum = next;
    }
    // A finite convex average lies in [-scale, scale]; contain roundoff overshoot.
    Ok((sum / values.len() as f64).clamp(-1.0, 1.0) * scale)
}

fn replicate_estimate(replicate_means: &[f64]) -> Option<LhsEstimate> {
    if replicate_means.is_empty() { return None; }
    let mean = scaled_mean(replicate_means);
    let count = replicate_means.len();
    let standard_error = if count < 2 { None } else {
        let scale = replicate_means.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 { Some(0.0) } else {
            let square_sum: f64 = replicate_means.iter().map(|v| (v / scale - mean / scale).powi(2)).sum();
            let error = scale * (square_sum / (count as f64 * (count - 1) as f64)).sqrt();
            if error.is_finite() { Some(error) } else { None }
        }
    };
    Some(LhsEstimate { mean, standard_error })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plan(n: usize) -> LhsPlan {
        LhsPlan { target_qoi: "temperature".into(), compliance_threshold: Some(0.5),
            parameters: vec![ParameterUncertainty::uniform("a", 0.0, 1.0, "K"),
                ParameterUncertainty::uniform("b", 0.0, 1.0, "K")], correlation: CorrelationModel::Independent,
            replicates: 3, samples_per_replicate: n, seed: 42 }
    }
    fn model(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + 2.0 * x[1]) }
    fn model_id() -> ContentHash { hash_domain("lhs-test-model", b"a+2b") }

    #[test]
    fn g0_every_replicate_stratifies_every_coordinate_exactly_once() {
        for n in [2, 3, 17, 127] {
            let p = plan(n);
            let mut execution = LhsExecution::new(&p).unwrap();
            let mut coordinates = Vec::new();
            let report = execution.advance(usize::MAX, || false, |x| {
                coordinates.push(x.to_vec()); model(x)
            });
            assert_eq!(report.status, UqStatus::Complete);
            for replicate in coordinates.chunks_exact(n) {
                for dimension in 0..p.parameters.len() {
                    let mut bins: Vec<_> = replicate.iter().map(|x| {
                        assert!(x[dimension] > 0.0 && x[dimension] < 1.0);
                        (x[dimension] * n as f64).floor() as usize
                    }).collect();
                    bins.sort();
                    assert_eq!(bins, (0..n).collect::<Vec<_>>());
                }
            }
        }
    }

    #[test]
    fn g5_chunking_checkpoint_and_mid_evaluation_interruption_preserve_points() {
        let p = plan(17);
        let mut whole = LhsExecution::new(&p).unwrap();
        whole.advance(usize::MAX, || false, model);
        let mut split = LhsExecution::new(&p).unwrap();
        split.advance(9, || false, model);
        let mut interrupted = Vec::new();
        split.advance_interruptible(10, || false, |x| {
            interrupted = x.to_vec(); Ok::<Option<f64>, &'static str>(None)
        });
        assert_eq!(split.evaluations_attempted(), 9);
        let bytes = split.checkpoint(model_id()).unwrap();
        let mut resumed = LhsExecution::restore(&p, model_id(), &bytes).unwrap();
        resumed.advance(1, || false, |x| { assert_eq!(x, interrupted.as_slice()); model(x) });
        for allowance in [1, 7, 5, 200] { resumed.advance(allowance, || false, model); }
        assert_eq!(whole.observations(), resumed.observations());
        assert_eq!(whole.report(), resumed.report());
    }

    #[test]
    fn g0_only_complete_replicates_contribute_statistics() {
        let mut e = LhsExecution::new(&plan(8)).unwrap();
        let r = e.advance(7, || false, model);
        assert!(r.estimate.is_none());
        assert_eq!(r.samples_accepted, 7);
        let r = e.advance(1, || false, model);
        assert_eq!(r.completed_replicates, 1);
        assert!(r.estimate.unwrap().standard_error.is_none());
        let r = e.advance(8, || false, model);
        assert!(r.estimate.unwrap().standard_error.is_some());
        assert!(r.compliance.unwrap().standard_error.is_some());
    }

    #[test]
    fn g4_cancellation_during_design_rebuild_consumes_no_observation() {
        let p = plan(4096);
        let mut e = LhsExecution::new(&p).unwrap();
        let mut polls = 0;
        let r = e.advance(1, || { polls += 1; polls == 3 }, |_| -> Result<f64, &'static str> {
            panic!("must not evaluate an incomplete permutation")
        });
        assert_eq!(r.status, UqStatus::Cancelled);
        assert_eq!(r.samples_evaluated, 0);
        let mut fresh = LhsExecution::new(&p).unwrap();
        e.advance(1, || false, model); fresh.advance(1, || false, model);
        assert_eq!(e.observations(), fresh.observations());
    }

    #[test]
    fn g0_model_failure_is_terminal_and_cannot_resume_a_censored_prefix() {
        let mut e = LhsExecution::new(&plan(3)).unwrap();
        e.advance(3, || false, model);
        let r = e.advance(1, || false, |_| Err::<f64, _>("solver refused"));
        assert_eq!(r.status, UqStatus::Refused);
        assert_eq!(r.samples_evaluated, 4);
        assert!(r.estimate.is_none());
        assert!(r.rejection_reason.unwrap().contains("solver refused"));
        assert_eq!(e.checkpoint(model_id()).unwrap_err(), UqCheckpointError::NotResumable);
        let r = e.advance(100, || false, |_| -> Result<f64, &'static str> { panic!("terminal") });
        assert_eq!(r.samples_evaluated, 4);
    }

    #[test]
    fn g0_checkpoint_binds_layout_even_when_total_budget_is_unchanged() {
        let p = plan(12);
        let mut e = LhsExecution::new(&p).unwrap(); e.advance(5, || false, model);
        let bytes = e.checkpoint(model_id()).unwrap();
        let mut changed = p.clone(); changed.replicates = 6; changed.samples_per_replicate = 6;
        assert_eq!(LhsExecution::restore(&changed, model_id(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        assert_eq!(LhsExecution::restore(&p, hash_domain("other", b"model"), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        let mut changed = p.clone(); changed.seed += 1;
        assert!(LhsExecution::restore(&changed, model_id(), &bytes).is_err());
        let mut corrupt = bytes; corrupt[49] ^= 1;
        assert_eq!(LhsExecution::restore(&p, model_id(), &corrupt).unwrap_err(), UqCheckpointError::IntegrityMismatch);
    }

    #[test]
    fn g0_shared_admission_refuses_invalid_measures_and_unknown_dependence() {
        let mut p = plan(8); p.correlation = CorrelationModel::Unknown;
        assert!(LhsExecution::new(&p).is_err());
        let mut p = plan(8); p.parameters[0] = ParameterUncertainty::interval("a", 0.0, 1.0, "K");
        assert!(LhsExecution::new(&p).is_err());
        let mut p = plan(8); p.parameters[1].name = p.parameters[0].name.clone();
        assert!(LhsExecution::new(&p).is_err());
        let mut p = plan(8); p.samples_per_replicate = usize::MAX;
        assert!(LhsExecution::new(&p).is_err());
        let mut p = plan(400_000); p.replicates = 3;
        assert!(LhsExecution::new(&p).is_err());
    }

    #[test]
    fn g0_degenerate_uniform_and_gaussian_inputs_remain_exact() {
        let mut p = plan(17);
        p.parameters = vec![ParameterUncertainty::uniform("a", 3.0, 3.0, "K"),
            ParameterUncertainty::gaussian("b", 4.0, 0.0, "K")];
        let mut e = LhsExecution::new(&p).unwrap();
        let r = e.advance(usize::MAX, || false, |x| {
            assert_eq!(x, &[3.0, 4.0]); Ok::<_, &'static str>(f64::MAX)
        });
        assert_eq!(r.estimate.unwrap(), LhsEstimate { mean: f64::MAX, standard_error: Some(0.0) });
    }

    #[test]
    fn g1_additive_uniform_integral_has_the_stratification_error_bound() {
        let n = 127;
        let mut e = LhsExecution::new(&plan(n)).unwrap();
        let r = e.advance(usize::MAX, || false, model);
        // Each uniform marginal mean differs from 1/2 by less than 1/(2n).
        for mean in r.replicate_means { assert!((mean - 1.5).abs() < 1.5 / n as f64); }
    }

    #[test]
    fn g0_open_jitter_stays_inside_extreme_strata_at_maximum_admitted_size() {
        for n in [2, 17, 500_000] {
            for stratum in [0, 1, n - 1] {
                for bits in [0_u32, u32::MAX] {
                    let u = (stratum as f64 + (f64::from(bits) + 0.5) / 4_294_967_296.0) / n as f64;
                    assert!(u > 0.0 && u < 1.0);
                    assert_eq!((u * n as f64).floor() as usize, stratum);
                    assert!(fs_bo::phi_inv(u).is_finite());
                }
            }
        }
    }
}
