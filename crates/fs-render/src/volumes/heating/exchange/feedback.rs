//! Temperature-dependent volume emission and reabsorption on native enthalpy.
//!
//! The shared temperature-optical field supplies local absorption/scattering
//! packets. Every changed packet rebuilds the original reciprocal transport
//! operator; unchanged packets retain it and only reevaluate Planck exchange.
//! The fixed exchange owner remains responsible for signed heat, reciprocal
//! transfers, coefficient-error admission, native CFL and energy accounting.
//!
//! One replica retains common photon streams across its temperature trajectory.
//! Temperature, physical time and history count do not key those streams; H/2H
//! studies retain history prefixes. This defines a reproducible sampled optical
//! response as coefficients evolve. Its errors are correlated across time and
//! the empirical coefficient errors are not final-temperature uncertainty.
//! Frozen-coefficient secants provide the native explicit decay restriction;
//! they do not bound dκ/dT or replace feedback/time/photon refinement.

use super::super::feedback::OpticalFeedbackError;
use super::super::temperature_field::{TemperatureOpticalField, TemperatureOpticalState};
pub mod refinement;

use super::{
    EnergySums, RadiativeBath, RadiativeExchangeBudget, RadiativeExchangeEnergy,
    RadiativeExchangeError, RadiativeExchangePower, RadiativeExchangeTolerance,
    RadiativeThermalCell, VolumeRadiativeExchange,
};
use crate::volumes::transport::scattering::exchange::{
    RadiativeExchangeConfig, RadiativeExchangeOperator,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Aggregate work ceilings for one complete native advance, across every
/// source evaluation and optical rebuild. Field compilation/storage retains
/// the separate shared field budget; each operator also retains its own caps.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeFeedbackBudget {
    /// Includes cache hits: every actual native source query consumes one.
    pub max_source_evaluations: usize,
    /// Actual photon segments over every rebuilt optical operator together.
    pub max_traced_segments: u64,
    /// Maximum distinct native thermal-cell pairs in any admitted operator.
    pub max_thermal_pairs: usize,
}

/// Refusal of the original constitutive, radiative or native thermal owner.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeFeedbackError {
    /// Invalid physical inputs, clock or configuration.
    Invalid(&'static str),
    /// An explicit work/storage allowance was exhausted.
    Budget,
    /// The shared temperature coefficient field refused its state/domain.
    Coefficients(OpticalFeedbackError),
    /// Original reciprocal exchange, energy, arithmetic or thermal refusal.
    Exchange(RadiativeExchangeError),
}
impl From<OpticalFeedbackError> for RadiativeFeedbackError {
    fn from(error: OpticalFeedbackError) -> Self {
        Self::Coefficients(error)
    }
}
impl From<RadiativeExchangeError> for RadiativeFeedbackError {
    fn from(error: RadiativeExchangeError) -> Self {
        Self::Exchange(error)
    }
}
impl From<DepositionError> for RadiativeFeedbackError {
    fn from(error: DepositionError) -> Self {
        Self::Exchange(error.into())
    }
}
impl core::fmt::Display for RadiativeFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "temperature-dependent radiative exchange refused: {self:?}"
        )
    }
}
impl core::error::Error for RadiativeFeedbackError {}

struct Model<'g, const N: usize> {
    identity: ContentHash,
    field: Arc<TemperatureOpticalField<'g, N>>,
    weights_nm: [f64; N],
    bath: RadiativeBath,
    bindings: Vec<RadiativeThermalCell>,
    config: RadiativeExchangeConfig,
    tolerance: RadiativeExchangeTolerance,
}

/// Transactional temperature feedback with shared immutable constitutive data.
/// The sampled optical response is rebuilt only when its coefficient packet or
/// sampling policy changes. Every native temperature update still reevaluates
/// net Planck exchange, even if the optical packet is constant.
#[derive(Clone)]
pub struct TemperatureRadiativeExchange<'g, const N: usize> {
    model: Arc<Model<'g, N>>,
    state: TemperatureOpticalState<N>,
    current: Option<VolumeRadiativeExchange<N>>,
    budget: RadiativeFeedbackBudget,
    histories: u32,
    replica_seed: u64,
    time_s: f64,
    evaluations: usize,
    segments: u64,
    operator_builds: usize,
    cache_valid: bool,
    zero_energy: RadiativeExchangeEnergy<N>,
}
impl<const N: usize> core::fmt::Debug for TemperatureRadiativeExchange<'_, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("TemperatureRadiativeExchange")
            .field("model_identity", &self.model.identity)
            .field("time_s", &self.time_s)
            .field("histories_per_emitter_lane", &self.histories)
            .field("source_evaluations", &self.evaluations)
            .field("traced_segments", &self.segments)
            .field("operator_builds", &self.operator_builds)
            .field("energy", &self.energy())
            .finish_non_exhaustive()
    }
}

impl<'g, const N: usize> TemperatureRadiativeExchange<'g, N> {
    /// Bind one shared local constitutive field to volume radiation. Geometry,
    /// density, reciprocal phase laws, quadrature and bath remain prescribed;
    /// only supported local κa(T,λ), κs(T,λ) and native temperatures evolve.
    /// Admission takes no photon draws. The first actual native source query
    /// charges all of its transport work to the same advance budget as later
    /// queries. The field may also be shared with an incident-light observer.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        field: Arc<TemperatureOpticalField<'g, N>>,
        wavelength_weights_nm: [f64; N],
        bath: RadiativeBath,
        thermal: &DepositionSimulation<'_>,
        budget: RadiativeFeedbackBudget,
        config: RadiativeExchangeConfig,
        tolerance: RadiativeExchangeTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, RadiativeFeedbackError> {
        checkpoint()?;
        if thermal.cells().is_empty() {
            return Err(invalid(
                "temperature-dependent radiation requires active native thermal cells",
            ));
        }
        if N == 0
            || budget.max_source_evaluations == 0
            || config.histories_per_emitter_lane < 2
            || config.max_collisions_per_history == 0
            || field.material_bindings().len() > config.max_material_cells
        {
            return Err(RadiativeFeedbackError::Budget);
        }
        if !wavelength_weights_nm.iter().copied().all(super::positive)
            || matches!(bath, RadiativeBath::Blackbody { temperature_k } if !super::positive(temperature_k))
            || ![
                tolerance.coefficient_roundoff_m2_sr,
                tolerance.power_roundoff_w,
                tolerance.max_cell_standard_error_w,
                tolerance.max_total_standard_error_w,
                tolerance.max_unresolved_exchange_w,
                tolerance.max_unresolved_standard_error_w,
            ]
            .into_iter()
            .all(super::nonnegative)
        {
            return Err(invalid(
                "positive quadrature/bath and finite nonnegative source error allowances required",
            ));
        }
        let state = field.initial_state(thermal, &mut checkpoint)?;
        // Retain the transport owner's coefficient/phase/constituent-source
        // admission even when this shared field was compiled at another state.
        super::transport_checkpoint(&mut checkpoint, |check| {
            field.field(&state, check).map(|_| ())
        })?;
        let mut bindings = Vec::new();
        bindings
            .try_reserve_exact(field.material_bindings().len())
            .map_err(|_| RadiativeFeedbackError::Budget)?;
        for binding in field.material_bindings() {
            checkpoint()?;
            bindings.push(RadiativeThermalCell {
                grid_cell: binding.grid_cell,
                thermal_cell: binding.thermal_cell,
                minimum_temperature_k: binding.minimum_temperature_k,
                maximum_temperature_k: binding.maximum_temperature_k,
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-radiative-exchange.v1");
        hash.update(field.identity().as_bytes());
        for weight in wavelength_weights_nm {
            hash.update(&weight.to_bits().to_le_bytes());
        }
        match bath {
            RadiativeBath::Vacuum => hash.update(&[0]),
            RadiativeBath::Blackbody { temperature_k } => {
                hash.update(&[1]);
                hash.update(&temperature_k.to_bits().to_le_bytes());
            }
        }
        let model = Arc::new(Model {
            identity: hash.finalize(),
            field,
            weights_nm: wavelength_weights_nm,
            bath,
            bindings,
            config,
            tolerance,
        });
        let zero_energy = EnergySums::<N>::default().value()?;
        checkpoint()?;
        Ok(Self {
            model,
            state,
            current: None,
            budget,
            histories: config.histories_per_emitter_lane,
            replica_seed: 0,
            time_s: thermal.time_s(),
            evaluations: 0,
            segments: 0,
            operator_builds: 0,
            cache_valid: false,
            zero_energy,
        })
    }

    /// Shared immutable temperature laws/geometry, usable by incident heating.
    pub fn field(&self) -> &Arc<TemperatureOpticalField<'g, N>> {
        &self.model.field
    }
    /// Physical field, quadrature and bath identity, independent of replica/H.
    pub fn model_identity(&self) -> ContentHash {
        self.model.identity
    }
    /// Original accepted native physical time.
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Actual accepted radiation energy, including prior continuation history.
    pub fn energy(&self) -> &RadiativeExchangeEnergy<N> {
        self.current
            .as_ref()
            .map_or(&self.zero_energy, VolumeRadiativeExchange::energy)
    }
    /// Last source-start rates; absent before the first actual source query.
    pub fn power(&self) -> Option<&RadiativeExchangePower<N>> {
        self.current.as_ref().map(VolumeRadiativeExchange::power)
    }
    /// Last sampled coefficient operator; absent before the first source query.
    pub fn operator(&self) -> Option<&RadiativeExchangeOperator<N>> {
        self.current.as_ref().map(VolumeRadiativeExchange::operator)
    }
    /// Last prepared signed bath sources, before native internal transfers.
    pub fn sources_w(&self) -> Option<&BTreeMap<CellId, f64>> {
        self.current
            .as_ref()
            .map(VolumeRadiativeExchange::sources_w)
    }
    /// Source queries during the most recent native advance, including hits.
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    /// Actual sampled segments, summed across rebuilt operators in the advance.
    pub fn traced_segments(&self) -> u64 {
        self.segments
    }
    /// Actual operator estimates in the most recent native advance.
    pub fn operator_builds(&self) -> usize {
        self.operator_builds
    }
    /// Effective whole-advance work limits.
    pub fn budget(&self) -> RadiativeFeedbackBudget {
        self.budget
    }
    /// Histories for every emitting cell/lane on the next operator rebuild.
    pub fn histories_per_emitter_lane(&self) -> u32 {
        self.histories
    }
    /// Current sampling policy with the replica key, before remaining-work cuts.
    pub fn sampling_config(&self) -> RadiativeExchangeConfig {
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-radiative-stream.v1");
        for value in [
            self.model.config.seed,
            self.model.config.source_id,
            self.replica_seed,
        ] {
            hash.update(&value.to_le_bytes());
        }
        let identity = hash.finalize();
        let bytes = identity.as_bytes();
        RadiativeExchangeConfig {
            seed: u64::from_le_bytes(bytes[..8].try_into().expect("eight hash bytes")),
            source_id: u64::from_le_bytes(bytes[8..16].try_into().expect("eight hash bytes")),
            histories_per_emitter_lane: self.histories,
            ..self.model.config
        }
    }
    /// Resample only future response; already accepted energy is retained.
    pub fn with_replica_seed(&self, seed: u64) -> Self {
        let mut result = self.clone();
        if result.replica_seed != seed {
            result.cache_valid = false;
        }
        result.replica_seed = seed;
        result
    }
    /// History count is absent from stream keys, preserving H/2H prefixes.
    pub fn with_histories_per_emitter_lane(
        &self,
        histories: u32,
    ) -> Result<Self, RadiativeFeedbackError> {
        if histories < 2 {
            return Err(RadiativeFeedbackError::Budget);
        }
        let mut result = self.clone();
        if result.histories != histories {
            result.cache_valid = false;
        }
        result.histories = histories;
        Ok(result)
    }
    /// Tighten the next whole native advance without replenishing work.
    pub fn with_work_limits(&self, evaluations: usize, segments: u64) -> Self {
        let mut result = self.clone();
        result.budget.max_source_evaluations =
            evaluations.min(result.budget.max_source_evaluations);
        result.budget.max_traced_segments = segments.min(result.budget.max_traced_segments);
        result
    }

    // One clone per ensemble trial; copying mutable packets/rate maps several
    // times through chained public setters would waste O(cells+pairs) work.
    pub(crate) fn configured_trial(
        &self,
        seed: u64,
        histories: u32,
        evaluations: usize,
        segments: u64,
    ) -> Result<Self, RadiativeFeedbackError> {
        if histories < 2 {
            return Err(RadiativeFeedbackError::Budget);
        }
        let mut result = self.clone();
        result.replica_seed = seed;
        result.histories = histories;
        result.budget.max_source_evaluations =
            evaluations.min(result.budget.max_source_evaluations);
        result.budget.max_traced_segments = segments.min(result.budget.max_traced_segments);
        // Each independent trajectory samples its future response even if its
        // first state happens to match the original observer's cached packet.
        result.cache_valid = false;
        Ok(result)
    }

    /// Restore all effective continuation ceilings after a bounded study.
    pub(crate) fn restore_budget(&mut self, budget: RadiativeFeedbackBudget) {
        self.budget = budget;
        if let Some(current) = &mut self.current {
            current.evaluation_limit = budget.max_source_evaluations;
        }
    }

    /// Original compensated accounts for joint interval comparisons.
    pub(in crate::volumes::heating) fn energy_sum(
        &self,
        channel: usize,
        lane: usize,
    ) -> super::super::CheckedSum {
        let Some(current) = &self.current else {
            return super::super::CheckedSum::default();
        };
        match channel {
            0 => current.sums.bath[lane],
            1 => current.sums.internal[lane],
            2 => current.sums.external,
            3 => current.sums.unresolved,
            _ => unreachable!("internal radiative account channel"),
        }
    }
    fn sums(&self) -> EnergySums<N> {
        self.current
            .as_ref()
            .map_or_else(EnergySums::default, |value| value.sums.clone())
    }
    fn check_clock(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), RadiativeFeedbackError> {
        if thermal.time_s() != self.time_s {
            Err(invalid(
                "radiative feedback and native thermal clocks differ",
            ))
        } else {
            Ok(())
        }
    }
    fn current_mut(&mut self) -> Result<&mut VolumeRadiativeExchange<N>, RadiativeFeedbackError> {
        if !self.cache_valid {
            return Err(invalid(
                "radiation requires a successful source preparation after configuration/state changes",
            ));
        }
        self.current.as_mut().ok_or_else(|| {
            invalid("radiation requires a prepared source at the native start state")
        })
    }
}

impl<const N: usize> DepositionObserver for TemperatureRadiativeExchange<'_, N> {
    type Error = RadiativeFeedbackError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        self.segments = 0;
        self.operator_builds = 0;
        if let Some(current) = &mut self.current {
            current.begin()?;
        }
        Ok(())
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        checkpoint()?;
        if self.evaluations >= self.budget.max_source_evaluations {
            return Err(RadiativeFeedbackError::Budget);
        }
        // A failed packet preparation cannot leave a stale operator marked
        // reusable, including direct trait callers outside native transactions.
        let previous_cache_valid = self.cache_valid;
        self.cache_valid = false;
        let changed = self
            .model
            .field
            .prepare(thermal, &mut self.state, checkpoint)?;
        let rebuild = !previous_cache_valid || changed || self.current.is_none();
        let sources = if rebuild {
            let mut config = self.sampling_config();
            config.max_total_segments = config.max_total_segments.min(
                self.budget
                    .max_traced_segments
                    .saturating_sub(self.segments),
            );
            let report = super::transport_checkpoint(checkpoint, |check| {
                self.model
                    .field
                    .field(&self.state, check)?
                    .estimate_thermal_exchange_with_checkpoint(check, config)
            })?;
            let segments = self
                .segments
                .checked_add(report.traced_segments)
                .ok_or(RadiativeFeedbackError::Budget)?;
            if segments > self.budget.max_traced_segments {
                return Err(RadiativeFeedbackError::Budget);
            }
            let builds = self
                .operator_builds
                .checked_add(1)
                .ok_or(RadiativeFeedbackError::Budget)?;
            let mut next = VolumeRadiativeExchange::try_new(
                report,
                self.model.weights_nm,
                self.model.bath,
                &self.model.bindings,
                thermal,
                RadiativeExchangeBudget {
                    max_bindings: self.model.bindings.len(),
                    max_thermal_pairs: self.budget.max_thermal_pairs,
                    max_source_evaluations: self.budget.max_source_evaluations,
                },
                self.model.tolerance,
                &mut *checkpoint,
            )?;
            if let Some(previous) = &self.current {
                next.sums = previous.sums.clone();
                next.energy = previous.energy;
            }
            // Its constructor already evaluated this exact native state once;
            // use those admitted rates rather than duplicate Planck evaluation.
            next.evaluations = self.evaluations + 1;
            let sources = next.sources_w().clone();
            checkpoint()?;
            self.current = Some(next);
            self.segments = segments;
            self.operator_builds = builds;
            sources
        } else {
            self.current
                .as_mut()
                .expect("an admitted cache has a current operator")
                .sources(thermal, checkpoint)?
        };
        self.evaluations += 1;
        self.cache_valid = true;
        Ok(sources)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.check_clock(thermal)?;
        Ok(self.current_mut()?.heat_transfers(thermal, checkpoint)?)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.check_clock(thermal)?;
        Ok(self
            .current_mut()?
            .limit_timestep_s(thermal, proposed_s, checkpoint)?)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        self.model
            .field
            .validate_temperatures(thermal, &mut || Ok(()))?;
        if let Some(current) = &mut self.current {
            current.event(event, thermal)?;
        }
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s {
            return Err(invalid(
                "radiative feedback needs the matching accepted native interval",
            ));
        }
        self.model
            .field
            .validate_temperatures(thermal, checkpoint)?;
        self.current_mut()?
            .step(start_s, duration_s, old, thermal, checkpoint)?;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn invalid(message: &'static str) -> RadiativeFeedbackError {
    RadiativeFeedbackError::Invalid(message)
}
