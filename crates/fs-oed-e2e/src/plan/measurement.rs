//! Outcome-free, resumable SensorForge decisions followed by actual readings.
//!
//! This uses the same independent-Gaussian opportunity-loss quadrature and
//! scalar Kalman kernel as the worked campaign. There is no oracle truth in
//! an alternative, study, or request. Model/roundoff/quadrature uncertainty is
//! not certified; EVPI dispersion is infinity. Each public operation has its
//! own bounded compute scope, not a renewed claim about a lifetime quota.

use fs_assimilate::{Belief, assimilate_colored_with_shared_poll_quota, point_sensor};
use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::{Color, color_leaf_identity_reason};
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_voi::{Action, ActionKind, DesignEstimate, Recommendation, Uncertainty};

use crate::{
    AcquisitionCost, ByteShape, CampaignBytePlan, CampaignProgress, CampaignWorkPlan,
    CandidateError, CanonicalDesignMenu, MAX_CAMPAIGN_CANDIDATES, MAX_CAMPAIGN_SENSORS,
    MAX_CANDIDATE_NAME_BYTES, NORMAL_EXPECTATION_RULE, ObjectiveSpec, ObjectiveValue,
    OedError, canonicalize_zero, predicted_posterior_variance, recommend_with_cancellation,
};

/// One design alternative and its independently noisy scalar measurement.
/// The declaration contains a prior, NOT a known true value or future reading.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasurementAlternative {
    name: String,
    prior: ObjectiveValue,
    variance: f64,
    noise: f64,
    cost: AcquisitionCost,
}

impl MeasurementAlternative {
    /// Check name, finite prior, Q-squared variances, and positive sensor noise.
    /// Acquisition-cost role and quantity schema are checked across the study.
    pub fn new(name: impl Into<String>, prior: ObjectiveValue, variance: QtyAny,
        noise: QtyAny, cost: AcquisitionCost) -> Result<Self, MeasurementError> {
        let name = name.into();
        if name.len() > MAX_CANDIDATE_NAME_BYTES || color_leaf_identity_reason(&name).is_some() {
            return Err(MeasurementError::Invalid("invalid alternative identity"));
        }
        let dims = prior.spec().variance_dims()?;
        for (field, quantity) in [("prior_variance", variance), ("sensor_noise_variance", noise)] {
            if quantity.dims != dims {
                return Err(CandidateError::DimensionMismatch { field, actual: quantity.dims, expected: dims }.into());
            }
        }
        if !variance.value.is_finite() || variance.value < 0.0 || !noise.value.is_finite() || noise.value <= 0.0 {
            return Err(MeasurementError::Invalid("variance must be nonnegative and noise strictly positive and finite"));
        }
        Ok(Self { name, prior, variance: canonicalize_zero(variance.value), noise: noise.value, cost })
    }

    /// Stable alternative identity; canonical ordering uses this token.
    #[must_use]
    pub fn name(&self) -> &str { &self.name }
}

/// Failures never publish a partial measurement update.
#[derive(Debug, Clone, PartialEq)]
pub enum MeasurementError {
    /// Existing typed OED, budget, or numerical refusal.
    Oed(OedError),
    /// Invalid dimensional or semantic declaration.
    Candidate(CandidateError),
    /// Invalid measurement protocol input.
    Invalid(&'static str),
    /// The request was planned for a different retained study state.
    StaleRequest,
    /// The reading is bound to a different measurement request.
    RequestMismatch,
}
impl From<OedError> for MeasurementError { fn from(e: OedError) -> Self { Self::Oed(e) } }
impl From<CandidateError> for MeasurementError { fn from(e: CandidateError) -> Self { Self::Candidate(e) } }
impl core::fmt::Display for MeasurementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Oed(e) => write!(f, "{e}"), Self::Candidate(e) => write!(f, "{e}"),
            Self::Invalid(reason) => f.write_str(reason),
            Self::StaleRequest => f.write_str("measurement request is stale or belongs to another study"),
            Self::RequestMismatch => f.write_str("reading does not belong to this request"),
        }
    }
}
impl std::error::Error for MeasurementError {}

/// A stop above the EVPI threshold is deliberately not called robust.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasurementStop { ThresholdReached, ReadingLimit, NoUsefulAction }

/// Planning never performs an acquisition or installs a predicted posterior.
#[derive(Debug, Clone, PartialEq)]
pub enum MeasurementDecision { Acquire(MeasurementRequest), Stop(MeasurementStop) }

/// A selected acquisition bound to the exact prior and complete study settings.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasurementRequest {
    history: ContentHash,
    identity: ContentHash,
    index: usize,
    ordinal: usize,
    alternative: MeasurementAlternative,
    predicted_variance: f64,
    score: f64,
    score_spec: ObjectiveSpec,
}
impl MeasurementRequest {
    #[must_use]
    pub fn alternative(&self) -> &str { &self.alternative.name }
    /// Zero-based accepted-reading ordinal, not an external attempt counter.
    #[must_use]
    pub const fn ordinal(&self) -> usize { self.ordinal }
    #[must_use]
    pub const fn identity(&self) -> ContentHash { self.identity }
    #[must_use]
    pub fn noise_variance(&self) -> QtyAny { self.variance_quantity(self.alternative.noise) }
    #[must_use]
    pub fn predicted_posterior_variance(&self) -> QtyAny { self.variance_quantity(self.predicted_variance) }
    /// Estimated opportunity-loss reduction per declared acquisition cost.
    #[must_use]
    pub fn value_per_cost(&self) -> ObjectiveValue { ObjectiveValue::from_raw(self.score, self.score_spec) }
    #[must_use]
    pub const fn acquisition_cost(&self) -> AcquisitionCost { self.alternative.cost }
    fn variance_quantity(&self, value: f64) -> QtyAny {
        QtyAny::new(value, self.alternative.prior.spec().variance_dims().expect("checked alternative"))
    }
    /// Attach a real reading under the EXACT objective schema and instrument ID.
    /// No acquisition occurs here; this is also usable after external I/O.
    pub fn reading(&self, value: ObjectiveValue, instrument: impl Into<String>) -> Result<MeasurementOutcome, MeasurementError> {
        if value.spec() != self.alternative.prior.spec() {
            return Err(MeasurementError::Invalid("reading objective schema differs from requested measurement"));
        }
        let instrument = instrument.into();
        if instrument.len() > 256 || color_leaf_identity_reason(&instrument).is_some() {
            return Err(MeasurementError::Invalid("invalid measurement instrument identity"));
        }
        Ok(MeasurementOutcome { request: self.identity, value, instrument })
    }
}

/// Actual acquired value, separated from both planning and simulation truth.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasurementOutcome { request: ContentHash, value: ObjectiveValue, instrument: String }
impl MeasurementOutcome {
    #[must_use]
    pub const fn value(&self) -> ObjectiveValue { self.value }
    #[must_use]
    pub fn instrument(&self) -> &str { &self.instrument }
    #[must_use]
    pub const fn request_identity(&self) -> ContentHash { self.request }
}

/// A completed actual update with its instrument-bound assimilation evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasurementUpdate {
    request: MeasurementRequest,
    outcome: MeasurementOutcome,
    posterior_mean: ObjectiveValue,
    posterior_variance: QtyAny,
    evpi: ObjectiveValue,
    color: Color,
}
impl MeasurementUpdate {
    #[must_use]
    pub const fn request(&self) -> &MeasurementRequest { &self.request }
    #[must_use]
    pub const fn outcome(&self) -> &MeasurementOutcome { &self.outcome }
    #[must_use]
    pub const fn posterior_mean(&self) -> ObjectiveValue { self.posterior_mean }
    #[must_use]
    pub const fn posterior_variance(&self) -> QtyAny { self.posterior_variance }
    #[must_use]
    pub const fn evpi(&self) -> ObjectiveValue { self.evpi }
    #[must_use]
    pub const fn assimilation_color(&self) -> &Color { &self.color }
}

/// Resumable decision-relevant sensing for independent Gaussian alternatives.
/// A failed operation leaves all accepted measurements and posteriors intact.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasurementStudy {
    alternatives: Vec<MeasurementAlternative>,
    beliefs: Vec<Belief>,
    threshold: ObjectiveValue,
    score_spec: ObjectiveSpec,
    max_readings: usize,
    updates: Vec<MeasurementUpdate>,
    evpi: f64,
    history: ContentHash,
}
impl MeasurementStudy {
    /// Admit the complete menu, canonicalize once, and evaluate the initial EVPI.
    pub fn new(mut alternatives: Vec<MeasurementAlternative>, threshold: ObjectiveValue,
        max_readings: usize, cx: &Cx<'_>) -> Result<Self, MeasurementError> {
        if alternatives.is_empty() { return Err(OedError::NoCandidates.into()); }
        if alternatives.len() > MAX_CAMPAIGN_CANDIDATES {
            return Err(OedError::TooManyCandidates { count: alternatives.len(), max: MAX_CAMPAIGN_CANDIDATES }.into());
        }
        if max_readings > MAX_CAMPAIGN_SENSORS {
            return Err(OedError::TooManySensors { count: max_readings, max: MAX_CAMPAIGN_SENSORS }.into());
        }
        let (plan, shape, mut progress) = work(&alternatives, cx)?;
        progress.checkpoint(cx, plan, "measurement study admission")?;
        let spec = alternatives[0].prior.spec();
        let cost_spec = alternatives[0].cost.spec();
        if threshold.spec() != spec.decision_spec() {
            return Err(OedError::ThresholdSchemaMismatch { actual: threshold.spec(), expected: spec.decision_spec() }.into());
        }
        if threshold.value() < 0.0 { return Err(OedError::InvalidThreshold.into()); }
        let score_spec = cost_spec.score_spec(spec.decision_spec()).ok_or(OedError::ScoreDimensionOverflow {
            benefit: spec.decision_spec().dims(), cost: cost_spec.dims(),
        })?;
        for alternative in &alternatives {
            progress.checkpoint(cx, plan, "measurement schema")?;
            if alternative.prior.spec() != spec {
                return Err(OedError::ObjectiveSchemaMismatch { candidate: alternative.name.clone(), actual: alternative.prior.spec(), expected: spec }.into());
            }
            if alternative.cost.spec() != cost_spec {
                return Err(OedError::AcquisitionCostSchemaMismatch { candidate: alternative.name.clone(), actual: alternative.cost.spec(), expected: cost_spec }.into());
            }
        }
        alternatives.sort_unstable_by(|a, b| a.name.cmp(&b.name));
        if let Some(pair) = alternatives.windows(2).find(|p| p[0].name == p[1].name) {
            return Err(OedError::DuplicateCandidate { name: pair[0].name.clone() }.into());
        }
        let mut beliefs = reserve(alternatives.len())?;
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.measurement-study.v1");
        hash.update(&spec.canonical_bytes());
        hash.update(&cost_spec.canonical_bytes());
        hash.update(&score_spec.canonical_bytes());
        hash.update(&(max_readings as u64).to_le_bytes());
        hash.update(&(alternatives.len() as u64).to_le_bytes());
        hash.update(&fs_voi::EVPI_SEMANTICS_VERSION.to_le_bytes());
        number(&mut hash, threshold.value());
        hash.update(&crate::CAMPAIGN_PLANNING_POLICY_VERSION.to_le_bytes());
        hash.update(&(fs_voi::EOL_QUADRATURE_PANELS as u64).to_le_bytes());
        for (node, weight) in NORMAL_EXPECTATION_RULE { number(&mut hash, node); number(&mut hash, weight); }
        for alternative in &alternatives {
            beliefs.push(Belief::scalar(alternative.prior.value(), alternative.variance).map_err(OedError::BeliefInvariant)?);
            field(&mut hash, alternative.name.as_bytes());
            for value in [alternative.prior.value(), alternative.variance, alternative.noise, alternative.cost.value()] { number(&mut hash, value); }
            progress.advance(1);
        }
        progress.charge_bytes("measurement declarations", shape.records() + shape.sum_name + 256)?;
        let mut study = Self { alternatives, beliefs, threshold, score_spec, max_readings,
            updates: reserve(max_readings)?, evpi: 0.0, history: hash.finalize() };
        study.evpi = study.menu(None, &mut progress, plan, shape, cx)?.full_opportunity_loss_checked()?;
        progress.charge_bytes("initial measurement EVPI", shape.full_scan())?;
        progress.checkpoint(cx, plan, "measurement study publication")?;
        Ok(study)
    }

    #[must_use]
    pub const fn history_id(&self) -> ContentHash { self.history }
    #[must_use]
    pub fn updates(&self) -> &[MeasurementUpdate] { &self.updates }
    #[must_use]
    pub fn evpi(&self) -> ObjectiveValue { ObjectiveValue::from_raw(self.evpi, self.threshold.spec()) }
    /// Minimum posterior-mean alternative, with canonical-name tie breaking.
    #[must_use]
    pub fn best_alternative(&self) -> &str {
        let index = (0..self.beliefs.len()).min_by(|&a, &b| {
            self.beliefs[a].mean()[0].total_cmp(&self.beliefs[b].mean()[0]).then(a.cmp(&b))
        }).expect("admitted nonempty menu");
        &self.alternatives[index].name
    }
    /// Current mean under the original exact objective schema.
    #[must_use]
    pub fn posterior_mean(&self, name: &str) -> Option<ObjectiveValue> {
        let i = self.alternatives.binary_search_by(|a| a.name.as_str().cmp(name)).ok()?;
        Some(ObjectiveValue::from_raw(self.beliefs[i].mean()[0], self.alternatives[i].prior.spec()))
    }
    /// Current variance in Q-squared dimensions.
    #[must_use]
    pub fn posterior_variance(&self, name: &str) -> Option<QtyAny> {
        let i = self.alternatives.binary_search_by(|a| a.name.as_str().cmp(name)).ok()?;
        Some(QtyAny::new(self.beliefs[i].variance(0).ok()?, self.alternatives[i].prior.spec().variance_dims().ok()?))
    }
    /// Numerical model output, never a finite estimator-error bound or certificate.
    #[must_use]
    pub fn evpi_color(&self) -> Color {
        Color::Estimated { estimator: format!("sensorforge-live-evpi:v1:{}", self.history), dispersion: f64::INFINITY }
    }

    /// Choose the next measurement by outcome-integrated decision value/cost.
    /// This is the existing nine-node SensorForge rule, not variance-only ranking.
    pub fn recommend(&self, cx: &Cx<'_>) -> Result<MeasurementDecision, MeasurementError> {
        let (plan, shape, mut progress) = work(&self.alternatives, cx)?;
        progress.checkpoint(cx, plan, "measurement recommendation")?;
        if self.evpi <= self.threshold.value() { return Ok(MeasurementDecision::Stop(MeasurementStop::ThresholdReached)); }
        if self.updates.len() == self.max_readings { return Ok(MeasurementDecision::Stop(MeasurementStop::ReadingLimit)); }
        let menu = self.menu(None, &mut progress, plan, shape, cx)?;
        let mut actions = reserve(self.alternatives.len())?;
        for (alternative, belief) in self.alternatives.iter().zip(&self.beliefs) {
            let prior = belief.variance(0).map_err(OedError::BeliefInvariant)?;
            let posterior = predicted_posterior_variance(prior, alternative.noise)?;
            let reduction = if prior == 0.0 { 0.0 } else { (1.0 - (posterior / prior).clamp(0.0, 1.0).sqrt()).clamp(0.0, 1.0) };
            actions.push(Action { name: format!("measure-{}", alternative.name), kind: ActionKind::Sample,
                target_design: alternative.name.clone(), reduction, cost: alternative.cost.value() });
        }
        progress.advance(self.alternatives.len() as u128);
        progress.charge_bytes("measurement actions", shape.action_construction())?;
        let recommendation = recommend_with_cancellation(&menu, &actions, self.evpi, self.threshold.value(), plan, shape, &mut progress, cx)?;
        let Recommendation::Act { action, value_per_cost } = recommendation else {
            return Ok(MeasurementDecision::Stop(MeasurementStop::NoUsefulAction));
        };
        let index = actions.iter().position(|a| a.name == action).ok_or(OedError::UnknownRecommendation { action })?;
        let alternative = &self.alternatives[index];
        let predicted_variance = predicted_posterior_variance(self.beliefs[index].variance(0).map_err(OedError::BeliefInvariant)?, alternative.noise)?;
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.measurement-request.v1");
        hash.update(self.history.as_bytes());
        hash.update(&(index as u64).to_le_bytes());
        hash.update(&(self.updates.len() as u64).to_le_bytes());
        number(&mut hash, predicted_variance); number(&mut hash, value_per_cost);
        progress.checkpoint(cx, plan, "measurement request publication")?;
        Ok(MeasurementDecision::Acquire(MeasurementRequest { history: self.history, identity: hash.finalize(), index,
            ordinal: self.updates.len(), alternative: alternative.clone(), predicted_variance, score: value_per_cost, score_spec: self.score_spec }))
    }

    /// Apply one real outcome to its current request. Borrowed inputs survive a
    /// failed attempt. Publication occurs only after the next EVPI and evidence
    /// have been computed; no callback or fallible work follows the commit.
    pub fn assimilate(&mut self, request: &MeasurementRequest, outcome: &MeasurementOutcome,
        cx: &Cx<'_>) -> Result<&MeasurementUpdate, MeasurementError> {
        self.check_current(request)?;
        if outcome.request != request.identity { return Err(MeasurementError::RequestMismatch); }
        let (plan, shape, mut progress) = work(&self.alternatives, cx)?;
        progress.checkpoint(cx, plan, "actual measurement admission")?;
        let mut instrument_hash = DomainHasher::new("org.frankensim.sensorforge.measurement-instrument.v1");
        instrument_hash.update(request.identity.as_bytes()); field(&mut instrument_hash, outcome.instrument.as_bytes());
        let observation = point_sensor(0, 1, outcome.value.value(), request.alternative.noise,
            format!("measurement-{}", instrument_hash.finalize())).map_err(OedError::BeliefInvariant)?;
        let posterior = assimilate_colored_with_shared_poll_quota(&self.beliefs[request.index], &[observation],
            "sensor_count", 0.0, (self.updates.len() + 1) as f64, cx, &mut progress.polls_remaining)
            .map_err(|source| OedError::Assimilation { candidate: request.alternative.name.clone(), source })?;
        progress.advance(2);
        let next_evpi = self.menu(Some((request.index, posterior.belief())), &mut progress, plan, shape, cx)?.full_opportunity_loss_checked()?;
        progress.charge_bytes("actual measurement EVPI", shape.full_scan())?;
        let mean = posterior.belief().component_mean(0).map_err(OedError::BeliefInvariant)?;
        let variance = posterior.belief().variance(0).map_err(OedError::BeliefInvariant)?;
        let update = MeasurementUpdate { request: request.clone(), outcome: outcome.clone(),
            posterior_mean: ObjectiveValue::from_raw(mean, request.alternative.prior.spec()),
            posterior_variance: request.variance_quantity(variance), evpi: ObjectiveValue::from_raw(next_evpi, self.threshold.spec()),
            color: posterior.color().clone() };
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.measurement-history.v1");
        hash.update(self.history.as_bytes()); hash.update(request.identity.as_bytes());
        field(&mut hash, outcome.instrument.as_bytes());
        for value in [outcome.value.value(), mean, variance, next_evpi] { number(&mut hash, value); }
        let history = hash.finalize();
        // Cloning a study need not preserve spare Vec capacity. Reserve before
        // publication even though the original constructor reserved the cap.
        self.updates.try_reserve(1).map_err(|_| OedError::OutputAllocationRefused { what: "measurement update" })?;
        let next_belief = posterior.belief().clone();
        progress.checkpoint(cx, plan, "actual measurement commit")?;
        self.beliefs[request.index] = next_belief;
        self.evpi = next_evpi; self.history = history;
        self.updates.push(update); // capacity was checked before publication
        Ok(self.updates.last().expect("just committed one measurement"))
    }

    /// Validate a request before doing external I/O. This reserves no instrument.
    pub fn check_current(&self, request: &MeasurementRequest) -> Result<(), MeasurementError> {
        if request.history != self.history || request.ordinal != self.updates.len()
            || request.index >= self.alternatives.len() || self.updates.len() >= self.max_readings {
            return Err(MeasurementError::StaleRequest);
        }
        Ok(())
    }

    fn menu(&self, replacement: Option<(usize, &Belief)>, progress: &mut CampaignProgress<'_>,
        plan: CampaignWorkPlan, shape: ByteShape, cx: &Cx<'_>) -> Result<CanonicalDesignMenu, MeasurementError> {
        let mut estimates = reserve(self.alternatives.len())?;
        for (index, alternative) in self.alternatives.iter().enumerate() {
            progress.checkpoint(cx, plan, "measurement posterior menu")?;
            let belief = replacement.filter(|(i, _)| *i == index).map_or(&self.beliefs[index], |(_, b)| b);
            estimates.push(DesignEstimate::new(alternative.name.clone(), belief.component_mean(0).map_err(OedError::BeliefInvariant)?,
                Uncertainty { numerical: 0.0, statistical: belief.variance(0).map_err(OedError::BeliefInvariant)?.sqrt(), model: 0.0 }));
        }
        progress.advance(self.alternatives.len() as u128);
        progress.charge_bytes("measurement posterior menu", shape.menu_build())?;
        Ok(CanonicalDesignMenu::from_canonical(estimates)?)
    }
}

fn work<'s>(alternatives: &[MeasurementAlternative], cx: &Cx<'s>) -> Result<(CampaignWorkPlan, ByteShape, CampaignProgress<'s>), MeasurementError> {
    let n = alternatives.len();
    let sum_name = alternatives.iter().map(|a| a.name.len()).sum();
    let max_name = alternatives.iter().map(|a| a.name.len()).max().unwrap_or(0);
    // Each API call is bounded by one complete legacy campaign step, including
    // its setup/finalization reserve. No lifetime or exact realized-work claim.
    let plan = CampaignWorkPlan::checked(n, 1)?;
    let bytes = CampaignBytePlan::checked(n, 1, sum_name, max_name)?;
    let progress = CampaignProgress::admit(cx, n as u128, plan.admitted_work_units as u64, bytes)?;
    Ok((plan, ByteShape { n: n as u128, sum_name: sum_name as u128, max_name: max_name as u128 }, progress))
}
fn reserve<T>(count: usize) -> Result<Vec<T>, MeasurementError> {
    let mut result = Vec::new();
    result.try_reserve_exact(count).map_err(|_| OedError::OutputAllocationRefused { what: "measurement study storage" })?;
    Ok(result)
}
fn number(hash: &mut DomainHasher, value: f64) { hash.update(&canonicalize_zero(value).to_le_bytes()); }
fn field(hash: &mut DomainHasher, value: &[u8]) { hash.update(&(value.len() as u64).to_le_bytes()); hash.update(value); }
