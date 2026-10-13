//! Greedy, cost-budgeted selection of a jointly informative measurement set.
//!
//! Recondition the target AND every remaining sensor after each selection.
//! Square-root residual rows retain cross-covariance, including the independent
//! noise of earlier selected sensors. No invented readings or mean updates are
//! used. This is a conditional variance-per-cost heuristic, not a globally
//! optimal subset, EVPI, or a certificate of the supplied Gaussian model.

use super::{SharedNoiseProbe, SharedNoiseStream, copy_string, norm, project, reserve, square, zeros};
use crate::{AssimError, Cx, WorkPlan, WorkProgress, validate_observation_numeric_parts};
use fs_blake3::ContentHash;

mod conditional;
use conditional::condition;

/// A candidate with a positive cost in one caller-declared integer unit.
#[derive(Debug, Clone, PartialEq)]
pub struct CostedProbe {
    probe: SharedNoiseProbe,
    cost: u32,
}

impl CostedProbe {
    /// Bind an acquisition cost, distinct from the compute budget in `Cx`.
    /// All candidates and the budget must use the SAME caller-selected unit.
    ///
    /// # Errors
    /// Refuses zero cost; free acquisitions need an explicit separate policy.
    pub fn new(probe: SharedNoiseProbe, cost: u32) -> Result<Self, AssimError> {
        if cost == 0 {
            return Err(declaration("probe acquisition cost must be positive"));
        }
        Ok(Self { probe, cost })
    }

    /// Immutable prospective measurement; contains no measured outcome.
    #[must_use]
    pub const fn probe(&self) -> &SharedNoiseProbe { &self.probe }

    /// Positive cost in the caller's common acquisition unit.
    #[must_use]
    pub const fn cost(&self) -> u32 { self.cost }
}

/// Explicit acquisition limits; computation is separately bounded by `Cx`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProbeBatchBudget {
    max_readings: usize,
    max_cost: u64,
}

impl ProbeBatchBudget {
    /// Construct acquisition limits. Either zero limit gives an empty plan.
    ///
    /// # Errors
    /// Refuses more than 256 requested measurements in a single plan.
    pub fn new(max_readings: usize, max_cost: u64) -> Result<Self, AssimError> {
        if max_readings > 256 {
            return Err(AssimError::ObservationCountLimit { count: max_readings, max: 256 });
        }
        Ok(Self { max_readings, max_cost })
    }

    /// Maximum measurements selected in this plan.
    #[must_use]
    pub const fn max_readings(&self) -> usize { self.max_readings }

    /// Maximum total acquisition cost, in the same unit as `CostedProbe`.
    #[must_use]
    pub const fn max_cost(&self) -> u64 { self.max_cost }
}

/// Reason selection ended; none implies global optimality or sufficiency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeBatchStop {
    /// The requested count was reached, including an explicit zero count.
    ReadingLimit,
    /// The live stream has no further replay-protected row slots available.
    RowCapacity,
    /// All candidates have been selected.
    CandidatesExhausted,
    /// No remaining candidate fits the unspent acquisition budget.
    CostLimit,
    /// Every affordable candidate has exactly zero computed target information.
    NoInformation,
}

/// One chosen measurement scored CONDITIONAL on the preceding plan steps.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeBatchStep {
    candidate_index: usize,
    candidate: CostedProbe,
    before: f64,
    after: f64,
    reduction: f64,
    innovation_variance: f64,
}

impl ProbeBatchStep {
    /// Index in the original candidate list, not its sorted position.
    #[must_use]
    pub const fn candidate_index(&self) -> usize { self.candidate_index }
    /// The exact response/noise/row/cost declaration selected.
    #[must_use]
    pub const fn candidate(&self) -> &CostedProbe { &self.candidate }
    /// Target variance before this conditional step.
    #[must_use]
    pub const fn prior_target_variance(&self) -> f64 { self.before }
    /// Target variance after this step, independent of acquired outcomes.
    #[must_use]
    pub const fn posterior_target_variance(&self) -> f64 { self.after }
    /// Computed conditional marginal reduction, not the initial one-step score.
    #[must_use]
    pub const fn variance_reduction(&self) -> f64 { self.reduction }
    /// Conditional measurement variance including its independent noise.
    #[must_use]
    pub const fn innovation_variance(&self) -> f64 { self.innovation_variance }
}

/// Immutable selection bound to a retained posterior history and physical goal.
/// The live stream, accepted rows, and history are not modified by planning.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeBatchPlan {
    history: ContentHash,
    model_id: String,
    physical_dim: usize,
    target: Vec<f64>,
    budget: ProbeBatchBudget,
    steps: Vec<ProbeBatchStep>,
    prior_variance: f64,
    posterior_variance: f64,
    total_cost: u64,
    stop: ProbeBatchStop,
}

impl ProbeBatchPlan {
    /// Posterior history on which every conditional forecast depends.
    #[must_use]
    pub const fn history_id(&self) -> ContentHash { self.history }
    /// Linear physical quantity whose uncertainty is being reduced.
    #[must_use]
    pub fn target(&self) -> &[f64] { &self.target }
    /// Explicit acquisition limits used for selection.
    #[must_use]
    pub const fn budget(&self) -> ProbeBatchBudget { self.budget }
    /// Selection order with successively conditioned scores.
    #[must_use]
    pub fn steps(&self) -> &[ProbeBatchStep] { &self.steps }
    /// Target variance at the live posterior, before any selected acquisition.
    #[must_use]
    pub const fn prior_target_variance(&self) -> f64 { self.prior_variance }
    /// Predicted target variance after acquiring the entire selected set.
    #[must_use]
    pub const fn posterior_target_variance(&self) -> f64 { self.posterior_variance }
    /// Sum of selected costs; never exceeds the requested cost limit.
    #[must_use]
    pub const fn total_cost(&self) -> u64 { self.total_cost }
    /// Explicit reason selection stopped.
    #[must_use]
    pub const fn stop(&self) -> ProbeBatchStop { self.stop }
}

impl SharedNoiseStream {
    /// Select a set by greedy conditional target-variance reduction per cost.
    /// Candidates may be redundant or share signed calibration uncertainty.
    /// Each selection changes the forecast covariance of every remaining row;
    /// taking the first few independent `rank_probes` results is NOT equivalent.
    /// Ties in the computed cost-normalized score use row identity.
    ///
    /// Project once through the retained UD factor, then append one square-root
    /// noise coordinate per selected reading and apply Joseph residual actions.
    /// Storage is O((m+1)(d+k)); work is O(m d² + k m (d+k)). Here d is joint
    /// state size, m candidates, k selections, all bounded by 256. No dense
    /// covariance inverse, mock observations, or candidate mean is required.
    /// One compute accountant covers admission, projection and ALL selections.
    /// The candidate noise must be mutually independent conditional on the
    /// retained nuisance, just as in `JacobianBatch` acquisition.
    ///
    /// # Errors
    /// Refuses bad target/model/dimensions, duplicate or previously acquired
    /// rows, oversized candidate lists, numerical range/contraction failures,
    /// allocation failure, cancellation, or exhausted cumulative compute work.
    /// Even unaffordable candidates undergo identity/dimension admission.
    pub fn plan_probe_batch(
        &self,
        target: &[f64],
        candidates: &[CostedProbe],
        budget: ProbeBatchBudget,
        cx: &Cx<'_>,
    ) -> Result<ProbeBatchPlan, AssimError> {
        if target.len() != self.physical_dim {
            return Err(AssimError::DimMismatch { state: self.physical_dim, operator: target.len() });
        }
        validate_observation_numeric_parts(target, 0.0, 1.0)?;
        if candidates.is_empty() { return Err(AssimError::EmptyObservations); }
        if candidates.len() > 256 {
            return Err(AssimError::ObservationCountLimit { count: candidates.len(), max: 256 });
        }
        let m = candidates.len();
        let d = self.joint.dim();
        let slots = self.row_capacity - self.seen.len();
        let k = budget.max_readings.min(slots).min(m);
        let (mu, du, ku) = (m as u128, d as u128, k as u128);
        let work = (mu + 1) * (8 * du * du + 64 * du + 16_384)
            + 256 * mu * mu + (ku + 1) * (mu + 1) * (256 * (du + ku + 1) + 1_024);
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, work, 0, 0)?)?;
        progress.checkpoint("batch probe admission")?;
        for (index, candidate) in candidates.iter().enumerate() {
            let probe = &candidate.probe;
            if probe.model_id != self.model_id {
                return Err(declaration("batch probe differs from retained shared model"));
            }
            if probe.physical_dim != self.physical_dim {
                return Err(AssimError::DimMismatch { state: self.physical_dim, operator: probe.physical_dim });
            }
            progress.comparison_bytes("batch probe row lookup", 13 * 256)?;
            if self.seen.binary_search_by(|seen| seen.row.cmp(&probe.row)).is_ok() {
                return Err(AssimError::DuplicateDatasetRow { id: copy_string(probe.row.as_str())? });
            }
            for earlier in &candidates[..index] {
                progress.comparison_bytes("batch probe uniqueness", 256)?;
                if earlier.probe.row == probe.row {
                    return Err(AssimError::DuplicateDatasetRow { id: copy_string(probe.row.as_str())? });
                }
            }
        }
        let mut target_copy = reserve(target.len())?;
        target_copy.extend_from_slice(target);
        let model_id = copy_string(&self.model_id)?;
        let mut augmented = zeros(d)?;
        augmented[..target.len()].copy_from_slice(target);
        let mut target_factor = project(&self.joint, &augmented, &mut progress)?;
        let prior_variance = square(norm(&target_factor, &mut progress)?)?;
        let mut posterior_variance = prior_variance;
        let mut features = reserve(m)?;
        for candidate in candidates {
            features.push(project(&self.joint, &candidate.probe.operator, &mut progress)?);
        }
        let mut selected = reserve(m)?;
        selected.resize(m, false);
        let mut steps = reserve(k)?;
        let mut total_cost = 0_u64;
        let stop = loop {
            progress.checkpoint("batch probe selection")?;
            if steps.len() == budget.max_readings { break ProbeBatchStop::ReadingLimit; }
            if steps.len() == slots { break ProbeBatchStop::RowCapacity; }
            if steps.len() == m { break ProbeBatchStop::CandidatesExhausted; }
            let mut best: Option<(usize, f64, conditional::Conditional)> = None;
            let mut affordable = false;
            for (index, candidate) in candidates.iter().enumerate() {
                progress.scalar("batch probe candidate", 1)?;
                if selected[index] || u64::from(candidate.cost) > budget.max_cost - total_cost { continue; }
                affordable = true;
                let forecast = condition(&target_factor, &features[index], candidate.probe.noise_variance, &mut progress)?;
                if forecast.reduction == 0.0 { continue; }
                let score = forecast.reduction / f64::from(candidate.cost);
                if score == 0.0 { return Err(AssimError::CovarianceCertificationUnresolved); }
                progress.comparison_bytes("batch probe tie", 256)?;
                let better = match &best {
                    None => true,
                    Some((previous, previous_score, _)) => score > *previous_score
                        || (score == *previous_score && candidate.probe.row < candidates[*previous].probe.row),
                };
                if better { best = Some((index, score, forecast)); }
            }
            let Some((index, _, forecast)) = best else {
                break if affordable { ProbeBatchStop::NoInformation } else { ProbeBatchStop::CostLimit };
            };
            // All rows use the same selected sensor, including the SAME sign
            // for its newly introduced independent-noise coordinate.
            let mut chosen = reserve(features[index].len())?;
            chosen.extend_from_slice(&features[index]);
            for other in 0..m {
                if other != index && !selected[other] {
                    features[other] = condition(&features[other], &chosen,
                        candidates[index].probe.noise_variance, &mut progress)?.residual;
                }
            }
            selected[index] = true;
            let candidate = copy_candidate(&candidates[index])?;
            let before = posterior_variance;
            posterior_variance = square(norm(&forecast.residual, &mut progress)?)?;
            if posterior_variance > before { return Err(AssimError::CovarianceCertificationUnresolved); }
            target_factor = forecast.residual;
            total_cost += u64::from(candidate.cost); // <= 256 * u32::MAX, admitted above
            steps.push(ProbeBatchStep { candidate_index: index, candidate, before,
                after: posterior_variance, reduction: forecast.reduction,
                innovation_variance: forecast.innovation_variance });
        };
        progress.checkpoint("batch probe publication")?;
        Ok(ProbeBatchPlan { history: self.history, model_id, physical_dim: self.physical_dim,
            target: target_copy, budget, steps, prior_variance, posterior_variance, total_cost, stop })
    }
}

fn declaration(reason: &'static str) -> AssimError {
    AssimError::SharedSourceDeclaration { id: "shared-probe-batch".to_owned(), reason }
}

fn copy_candidate(candidate: &CostedProbe) -> Result<CostedProbe, AssimError> {
    let probe = &candidate.probe;
    let mut operator = reserve(probe.operator.len())?;
    operator.extend_from_slice(&probe.operator);
    Ok(CostedProbe { cost: candidate.cost, probe: SharedNoiseProbe {
        row: crate::groups::RowId::try_new(copy_string(probe.row.as_str())?)?,
        model_id: copy_string(&probe.model_id)?, physical_dim: probe.physical_dim,
        operator, noise_variance: probe.noise_variance,
    } })
}
