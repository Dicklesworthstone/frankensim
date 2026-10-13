//! Explicit synthetic acquisitions, never a source of truth for the planner.
//!
//! Every draw uses the requested Gaussian instrument variance. Its counter-based
//! stream is a function of the explicit seed and complete request identity, not
//! call order, worker identity, wall time or the Cx's unrelated random stream.
//! Re-reading one request replays one acquisition; a new request gets a new draw.

use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::color_leaf_identity_reason;
use fs_exec::{AdmittedBudget, Cx};
use fs_qty::QtyAny;
use fs_qty::semantic::SemanticQty;
use fs_rand::StreamKey;

use super::MeasurementProvider;
use crate::plan::measurement::{MeasurementError, MeasurementOutcome, MeasurementRequest};
use crate::{MAX_CAMPAIGN_CANDIDATES, MAX_CANDIDATE_NAME_BYTES, ObjectiveValue, OedError};

/// Binds random-stream addressing, noise scaling and synthetic outcome identity.
pub const SIMULATION_VERSION: u32 = 1;

/// A caller-declared latent value confined to the simulation provider.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulationTruth {
    alternative: String,
    value: ObjectiveValue,
}

impl SimulationTruth {
    /// Admit a bounded alternative name and already-checked typed latent value.
    pub fn new(alternative: impl Into<String>, value: ObjectiveValue) -> Result<Self, MeasurementError> {
        let alternative = alternative.into();
        if alternative.len() > MAX_CANDIDATE_NAME_BYTES || color_leaf_identity_reason(&alternative).is_some() {
            return Err(MeasurementError::Invalid("invalid simulated alternative identity"));
        }
        Ok(Self { alternative, value })
    }
}

/// Stateless, versioned Gaussian measurement provider for synthetic experiments.
/// It does not claim that a supplied latent value or noise model is calibrated.
#[derive(Debug, Clone, PartialEq)]
pub struct GaussianMeasurementSimulator {
    seed: u64,
    truths: Vec<SimulationTruth>,
    identity: ContentHash,
}

/// Simulation-only diagnostics. Only `outcome()` is passed to production sensing.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulatedMeasurement {
    outcome: MeasurementOutcome,
    standard_normal: f64,
    perturbation: ObjectiveValue,
    key: StreamKey,
    source: ContentHash,
}

impl SimulatedMeasurement {
    /// The actual noisy value, explicitly identified as simulated acquisition.
    #[must_use]
    pub const fn outcome(&self) -> &MeasurementOutcome { &self.outcome }
    /// Dimensionless draw before scaling by the declared standard deviation.
    #[must_use]
    pub const fn standard_normal(&self) -> f64 { self.standard_normal }
    /// Additive noise in the objective's difference schema (delta for temperature).
    #[must_use]
    pub const fn perturbation(&self) -> ObjectiveValue { self.perturbation }
    /// Exact Philox key. This protocol consumes normal draws at indices 0 and 1.
    #[must_use]
    pub const fn stream_key(&self) -> StreamKey { self.key }
    /// Full input identity binding seed, named latent values and sampler versions.
    #[must_use]
    pub const fn simulator_identity(&self) -> ContentHash { self.source }
}

impl GaussianMeasurementSimulator {
    /// Admit a nonempty, unique, same-schema latent menu of at most 256 entries.
    /// Truth declaration order does not affect either draws or source identity.
    pub fn new(seed: u64, mut truths: Vec<SimulationTruth>, cx: &Cx<'_>) -> Result<Self, MeasurementError> {
        if truths.is_empty() || truths.len() > MAX_CAMPAIGN_CANDIDATES {
            return Err(MeasurementError::Invalid("simulation needs 1..=256 latent alternatives"));
        }
        let mut budget = AdmittedBudget::admit_ambient(cx, 32_768 * truths.len() as u64)
            .map_err(OedError::BudgetRefused)?;
        budget.checkpoint("simulation declarations", cx).map_err(OedError::BudgetRefused)?;
        truths.sort_unstable_by(|a, b| a.alternative.cmp(&b.alternative));
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.gaussian-simulator.v1");
        hash.update(&SIMULATION_VERSION.to_le_bytes());
        hash.update(&fs_rand::STREAM_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&seed.to_le_bytes());
        hash.update(&(truths.len() as u64).to_le_bytes());
        for (index, truth) in truths.iter().enumerate() {
            budget.checkpoint("simulation latent row", cx).map_err(OedError::BudgetRefused)?;
            if index > 0 && truths[index - 1].alternative == truth.alternative {
                return Err(MeasurementError::Invalid("duplicate simulated alternative"));
            }
            if truth.value.spec() != truths[0].value.spec() {
                return Err(MeasurementError::Invalid("simulated alternatives have different objective schemas"));
            }
            hash.update(&(truth.alternative.len() as u64).to_le_bytes());
            hash.update(truth.alternative.as_bytes());
            hash.update(&truth.value.spec().quantity_spec().canonical_bytes());
            hash.update(&truth.value.to_bits().to_le_bytes());
            budget.charge_cost("simulation latent row", 32_768).map_err(OedError::BudgetRefused)?;
        }
        budget.checkpoint("simulation publication", cx).map_err(OedError::BudgetRefused)?;
        Ok(Self { seed, truths, identity: hash.finalize() })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash { self.identity }

    /// Sample one request without changing either provider or numerical study.
    /// The full request hash and explicit seed are domain-hashed into a 128-bit
    /// Philox substream key; that derivation is versioned, not an injectivity or
    /// authentication claim. Latent values do NOT select the random substream,
    /// enabling paired comparisons with common noise for the same request.
    ///
    /// Nonfinite arithmetic, erased nonzero noise, invalid semantic outcomes,
    /// unknown alternatives, cancellation and exhausted budgets refuse. No
    /// clipping, rejection-resampling, jitter or implicit exact-truth fallback.
    pub fn sample(&self, request: &MeasurementRequest, cx: &Cx<'_>) -> Result<SimulatedMeasurement, MeasurementError> {
        let mut budget = AdmittedBudget::admit_ambient(cx, 16_384).map_err(OedError::BudgetRefused)?;
        budget.checkpoint("simulation before draw", cx).map_err(OedError::BudgetRefused)?;
        let index = self.truths.binary_search_by(|t| t.alternative.as_str().cmp(request.alternative()))
            .map_err(|_| MeasurementError::Invalid("requested alternative is absent from simulation"))?;
        let truth = self.truths[index].value;
        // Reuse exact schema admission before generating any outcome.
        request.reading(truth, "simulation-schema-check")?;
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.gaussian-draw.v1");
        hash.update(&SIMULATION_VERSION.to_le_bytes());
        hash.update(&fs_rand::STREAM_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&self.seed.to_le_bytes());
        hash.update(request.identity().as_bytes());
        let address = hash.finalize();
        let bytes = address.as_bytes();
        let key = StreamKey {
            seed: u64::from_le_bytes(bytes[0..8].try_into().expect("fixed hash slice")),
            kernel: u32::from_le_bytes(bytes[8..12].try_into().expect("fixed hash slice")),
            tile: u32::from_le_bytes(bytes[12..16].try_into().expect("fixed hash slice")),
        };
        let z = key.stream().next_normal();
        let sigma = request.noise_variance().value.sqrt();
        let noise = sigma * z;
        let measured = sigma.mul_add(z, truth.value());
        if !z.is_finite() || !noise.is_finite() || !measured.is_finite()
            || (z != 0.0 && noise == 0.0) || (noise != 0.0 && measured == truth.value()) {
            return Err(MeasurementError::Invalid("simulated Gaussian perturbation is not representable"));
        }
        let quantity = QtyAny::new(measured, truth.spec().dims());
        let value = match truth.spec().semantic_type() {
            Some(kind) => ObjectiveValue::semantic(SemanticQty::new(quantity, kind)
                .map_err(|_| MeasurementError::Invalid("Gaussian outcome violates objective semantics"))?),
            None => ObjectiveValue::dimensional(quantity)?,
        };
        let perturbation = truth.spec().decision_value(noise)?;
        let mut hash = DomainHasher::new("org.frankensim.sensorforge.simulated-outcome.v1");
        hash.update(self.identity.as_bytes());
        hash.update(request.identity().as_bytes());
        hash.update(address.as_bytes());
        for number in [z, noise, measured] { hash.update(&number.to_bits().to_le_bytes()); }
        let outcome = request.reading(value, format!("simulated-gaussian-v1-{}", hash.finalize()))?;
        budget.charge_cost("simulation Gaussian draw", 16_384).map_err(OedError::BudgetRefused)?;
        budget.checkpoint("simulation after draw", cx).map_err(OedError::BudgetRefused)?;
        Ok(SimulatedMeasurement { outcome, standard_normal: z, perturbation, key, source: self.identity })
    }
}

impl MeasurementProvider for GaussianMeasurementSimulator {
    type Error = MeasurementError;

    fn read(&mut self, request: &MeasurementRequest, cx: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
        self.sample(request, cx).map(|sample| sample.outcome)
    }
}

/// Resumable stochastic studies through the production ask/acquire/update path.
pub mod study;
