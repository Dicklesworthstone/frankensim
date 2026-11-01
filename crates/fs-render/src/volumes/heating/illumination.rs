//! Several physical incident sources on one native temperature-dependent field.
//!
//! A finite source set may describe lamps or an explicitly weighted spatial/
//! angular quadrature. Every ray carries watts; its weights are already included
//! in that power. The caller owns quadrature accuracy. One coefficient packet is
//! prepared for the entire source set, and the original multiple-scattering and
//! heating owners retain absorption, escape and unresolved energy separately.
//!
//! Source IDs are unique and canonicalize reduction order. Independent source
//! and replica streams remain common across temperature, time and photon-count
//! comparisons. An unchanged optical packet reuses its immutable transport.
//! This sampled response field can have nonlinear finite-history feedback bias;
//! neither its source errors nor exact energy conservation certify that bias.

use super::feedback::OpticalFeedbackError;
use super::scattering::{ScatteringHeating, ScatteringHeatingTolerance};
use super::scattering_feedback::{
    EnergySums, IncidentScatteringBeam, ScatteringFeedbackEnergyBudget,
};
use super::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalState, transport_checkpoint,
};
use super::{CheckedSum, OpticalHeatingError, nonnegative};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::deposition::{
    ScatteringDepositionConfig, ScatteringPowerDeposition,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// One independently sampled ray or weighted source-quadrature node.
#[derive(Clone, Copy, Debug)]
pub struct IncidentIlluminationSource<const N: usize> {
    /// Stable physical source/node identity, unique within the illumination.
    /// Reordering sources does not alter streams or deterministic reductions.
    pub source_id: u64,
    pub beam: IncidentScatteringBeam<N>,
}

/// Whole-source-set admission and whole-native-advance work ceilings.
/// The admitted optical field already owns grid/table storage ceilings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IlluminationBudget {
    pub max_sources: usize,
    /// Counts each combined source query, including cached transport reuse.
    pub max_source_evaluations: usize,
    /// All source rays, wavelength lanes and rebuilds in one native advance.
    /// A cached or identically zero-power source set can use zero segments.
    pub max_traced_segments: u64,
    /// Sum of retained per-source cell records in ONE prepared illumination.
    /// A native candidate/replacement can coexist with its accepted predecessor;
    /// whole-trajectory memory admission must include both snapshots.
    pub max_retained_deposition_cells: usize,
}

struct Model<'g, const N: usize> {
    identity: ContentHash,
    field: Arc<TemperatureOpticalField<'g, N>>,
    sources: Vec<IncidentIlluminationSource<N>>,
    config: ScatteringDepositionConfig,
    tolerance: ScatteringHeatingTolerance,
}

/// Already-admitted per-source reports and their single combined heat field.
/// The original fixed heating owner integrates each physical report. Its own
/// elapsed-time state is unused here: the outer native observer owns intervals.
struct PreparedIllumination<const N: usize> {
    sources: Vec<ScatteringHeating<N>>,
    sources_w: BTreeMap<CellId, f64>,
    retained_cells: usize,
}

/// Transactional multiple-source absorption and scattering on native enthalpy.
#[derive(Clone)]
pub struct TemperatureIlluminationHeating<'g, const N: usize> {
    model: Arc<Model<'g, N>>,
    coefficients: TemperatureOpticalState<N>,
    current: Option<Arc<PreparedIllumination<N>>>,
    cache_valid: bool,
    budget: IlluminationBudget,
    histories: u32,
    replica_seed: u64,
    time_s: f64,
    evaluations: usize,
    segments: u64,
    builds: usize,
    sums: EnergySums<N>,
    energy: ScatteringFeedbackEnergyBudget<N>,
}

impl<const N: usize> core::fmt::Debug for TemperatureIlluminationHeating<'_, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("TemperatureIlluminationHeating")
            .field("model_identity", &self.model.identity)
            .field("sources", &self.source_count())
            .field("time_s", &self.time_s)
            .field("histories_per_lane", &self.histories)
            .field("source_evaluations", &self.evaluations)
            .field("traced_segments", &self.segments)
            .field("transport_builds", &self.builds)
            .field("energy", &self.energy)
            .finish_non_exhaustive()
    }
}

impl<'g, const N: usize> TemperatureIlluminationHeating<'g, N> {
    /// Admit a source set without tracing photons. Every source must use this
    /// field's supported spectral packet. The thermal-radiation observer may
    /// independently use a different packet; no infrared support is inferred.
    /// Empty illumination is legal and deposits exactly zero heat.
    #[allow(clippy::too_many_arguments)]
    pub fn from_field(
        field: Arc<TemperatureOpticalField<'g, N>>,
        sources: &[IncidentIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: IlluminationBudget,
        config: ScatteringDepositionConfig,
        tolerance: ScatteringHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        if N == 0
            || sources.len() > budget.max_sources
            || budget.max_source_evaluations == 0
            || config.histories_per_lane < 2
            || config.max_collisions_per_history == 0
        {
            return Err(OpticalFeedbackError::Budget);
        }
        if thermal.cells().is_empty()
            || ![
                tolerance.power_roundoff_w,
                tolerance.max_unresolved_w,
                tolerance.max_total_standard_error_w,
                tolerance.max_cell_standard_error_w,
            ]
            .into_iter()
            .all(nonnegative)
        {
            return Err(invalid(
                "active thermal cells and finite nonnegative illumination allowances required",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(sources.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        for source in sources {
            checkpoint()?;
            validate_source(&source.beam, &field)?;
            retained.push(*source);
        }
        retained.sort_unstable_by_key(|source| source.source_id);
        if retained
            .windows(2)
            .any(|pair| pair[0].source_id == pair[1].source_id)
        {
            return Err(invalid("illumination source identities must be unique"));
        }
        let mut identity = DomainHasher::new("org.frankensim.render.temperature-illumination.v1");
        identity.update(field.identity().as_bytes());
        identity.update(&(retained.len() as u64).to_le_bytes());
        for source in &retained {
            checkpoint()?;
            identity.update(&source.source_id.to_le_bytes());
            for value in [
                source.beam.origin.x,
                source.beam.origin.y,
                source.beam.origin.z,
                source.beam.unit_direction.x,
                source.beam.unit_direction.y,
                source.beam.unit_direction.z,
            ]
            .into_iter()
            .chain(source.beam.incident_power_w)
            {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        let coefficients = field.initial_state(thermal, &mut checkpoint)?;
        transport_checkpoint(&mut checkpoint, |check| {
            field.field(&coefficients, check).map(|_| ())
        })?;
        let sums = EnergySums::default();
        let energy = sums.value()?;
        Ok(Self {
            model: Arc::new(Model {
                identity: identity.finalize(),
                field,
                sources: retained,
                config,
                tolerance,
            }),
            coefficients,
            current: None,
            cache_valid: false,
            budget,
            histories: config.histories_per_lane,
            replica_seed: 0,
            time_s: thermal.time_s(),
            evaluations: 0,
            segments: 0,
            builds: 0,
            sums,
            energy,
        })
    }

    pub fn optical_field(&self) -> &Arc<TemperatureOpticalField<'g, N>> {
        &self.model.field
    }
    pub fn model_identity(&self) -> ContentHash {
        self.model.identity
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn energy(&self) -> &ScatteringFeedbackEnergyBudget<N> {
        &self.energy
    }
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn traced_segments(&self) -> u64 {
        self.segments
    }
    /// Number of complete source-set transports rebuilt in this advance.
    pub fn transport_builds(&self) -> usize {
        self.builds
    }
    pub fn histories_per_lane(&self) -> u32 {
        self.histories
    }
    pub fn budget(&self) -> IlluminationBudget {
        self.budget
    }
    pub fn source_count(&self) -> usize {
        self.model.sources.len()
    }
    pub fn has_scattering(&self) -> bool {
        self.model.field.has_scattering()
    }
    pub fn sampling_config(&self) -> ScatteringDepositionConfig {
        let mut config = self.model.config;
        config.histories_per_lane = self.histories;
        config
    }
    /// Canonical source-id order, independent of the caller's insertion order.
    pub fn incident_sources(&self) -> &[IncidentIlluminationSource<N>] {
        &self.model.sources
    }
    pub fn sources_w(&self) -> Option<&BTreeMap<CellId, f64>> {
        self.current
            .as_ref()
            .filter(|_| self.cache_valid)
            .map(|current| &current.sources_w)
    }
    /// Conditional report from the current prepared packet; no combined
    /// standard error is manufactured from independent scalar error bars.
    pub fn deposition(&self, source_id: u64) -> Option<&ScatteringPowerDeposition<N>> {
        if !self.cache_valid {
            return None;
        }
        let index = self
            .model
            .sources
            .binary_search_by_key(&source_id, |source| source.source_id)
            .ok()?;
        Some(self.current.as_ref()?.sources[index].deposition())
    }
    pub fn retained_deposition_cells(&self) -> usize {
        self.current
            .as_ref()
            .map_or(0, |current| current.retained_cells)
    }
    pub fn with_replica_seed(&self, replica_seed: u64) -> Self {
        let mut result = self.clone();
        result.cache_valid &= result.replica_seed == replica_seed;
        result.replica_seed = replica_seed;
        result
    }
    pub fn with_histories_per_lane(&self, histories: u32) -> Result<Self, OpticalFeedbackError> {
        if histories < 2 {
            return Err(OpticalFeedbackError::Budget);
        }
        let mut result = self.clone();
        result.cache_valid &= result.histories == histories;
        result.histories = histories;
        Ok(result)
    }
    /// Tighten future whole-advance work without invalidating reusable physics.
    pub fn with_work_limits(&self, evaluations: usize, segments: u64) -> Self {
        let mut result = self.clone();
        result.budget.max_source_evaluations = evaluations.min(self.budget.max_source_evaluations);
        result.budget.max_traced_segments = segments.min(self.budget.max_traced_segments);
        result
    }
    pub(super) fn configured_trial(
        &self,
        seed: u64,
        histories: u32,
        evaluations: usize,
        segments: u64,
    ) -> Result<Self, OpticalFeedbackError> {
        if histories < 2 {
            return Err(OpticalFeedbackError::Budget);
        }
        let mut result = self.clone();
        result.replica_seed = seed;
        result.histories = histories;
        result.budget.max_source_evaluations = evaluations.min(self.budget.max_source_evaluations);
        result.budget.max_traced_segments = segments.min(self.budget.max_traced_segments);
        result.cache_valid = false;
        Ok(result)
    }
    pub(super) fn restore_budget(&mut self, budget: IlluminationBudget) {
        self.budget = budget;
    }
    /// Compensated cumulative channels: incident, ballistic, diffuse, absorbed,
    /// unresolved, residual (each N lanes), followed by scalar thermal input.
    pub(super) fn energy_sum(&self, channel: usize, lane: usize) -> CheckedSum {
        self.sums.energy_sum(channel, lane)
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), OpticalFeedbackError> {
        if self.time_s != thermal.time_s() {
            Err(invalid("illumination and native thermal clocks differ"))
        } else {
            Ok(())
        }
    }
    fn source_config(
        &self,
        id: u64,
        remaining_segments: u64,
        remaining_cells: usize,
    ) -> ScatteringDepositionConfig {
        let mut hash =
            DomainHasher::new("org.frankensim.render.temperature-illumination-stream.v1");
        for value in [
            self.model.config.seed,
            self.model.config.source_id,
            self.replica_seed,
            id,
        ] {
            hash.update(&value.to_le_bytes());
        }
        let identity = hash.finalize();
        let bytes = identity.as_bytes();
        let mut config = self.model.config;
        config.seed = u64::from_le_bytes(bytes[..8].try_into().expect("eight hash bytes"));
        config.source_id = u64::from_le_bytes(bytes[8..16].try_into().expect("eight hash bytes"));
        config.histories_per_lane = self.histories;
        config.max_total_segments = config.max_total_segments.min(remaining_segments);
        config.max_deposition_cells = config.max_deposition_cells.min(remaining_cells);
        config
    }

    fn prepare_sources(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(PreparedIllumination<N>, u64), OpticalFeedbackError> {
        let available = self
            .budget
            .max_traced_segments
            .saturating_sub(self.segments);
        // Borrow/validate the shared optical packet ONCE for all source rays.
        let (reports, traced, retained_cells) = transport_checkpoint(checkpoint, |check| {
            let field = self.model.field.field(&self.coefficients, check)?;
            let mut reports = Vec::new();
            reports
                .try_reserve_exact(self.source_count())
                .map_err(|_| VolumeTransportError::Budget)?;
            let mut traced = 0_u64;
            let mut cells = 0_usize;
            for source in &self.model.sources {
                check()?;
                let config = self.source_config(
                    source.source_id,
                    available.saturating_sub(traced),
                    self.budget
                        .max_retained_deposition_cells
                        .saturating_sub(cells),
                );
                let report = field.deposit_multiple_scattering_power_with_checkpoint(
                    &mut || check(),
                    source.beam.origin,
                    source.beam.unit_direction,
                    &source.beam.incident_power_w,
                    config,
                )?;
                traced = traced
                    .checked_add(report.traced_segments)
                    .ok_or(VolumeTransportError::Budget)?;
                cells = cells
                    .checked_add(report.cells.len())
                    .ok_or(VolumeTransportError::Budget)?;
                if traced > available || cells > self.budget.max_retained_deposition_cells {
                    return Err(VolumeTransportError::Budget);
                }
                reports.push(report);
            }
            Ok((reports, traced, cells))
        })?;
        let mut accepted = Vec::new();
        accepted
            .try_reserve_exact(reports.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut errors = BTreeMap::<CellId, CheckedSum>::new();
        let mut total_absorbed_error = CheckedSum::default();
        let mut total_unresolved_error = CheckedSum::default();
        let mut unresolved = CheckedSum::default();
        let mut absorbed = CheckedSum::default();
        let mut residuals = [CheckedSum::default(); N];
        for report in reports {
            checkpoint()?;
            let mut bindings = Vec::new();
            bindings
                .try_reserve_exact(report.cells.len())
                .map_err(|_| OpticalFeedbackError::Budget)?;
            for cell in &report.cells {
                checkpoint()?;
                let receiver = self
                    .model
                    .field
                    .receiver(cell.cell_index)
                    .ok_or(OpticalHeatingError::UnmappedCell(cell.cell_index))?;
                bindings.push((cell.cell_index, receiver));
                for error in cell.absorbed.standard_error_w {
                    errors.entry(receiver).or_default().add(error)?;
                }
            }
            for lane in 0..N {
                total_absorbed_error.add(report.absorbed.standard_error_w[lane])?;
                total_unresolved_error.add(report.unresolved.standard_error_w[lane])?;
                unresolved.add(report.unresolved.mean_w[lane])?;
                absorbed.add(report.absorbed.mean_w[lane])?;
                residuals[lane].add(report.incident_power_w[lane])?;
                for power in [
                    report.ballistic_escape.mean_w[lane],
                    report.diffuse_escape.mean_w[lane],
                    report.absorbed.mean_w[lane],
                    report.unresolved.mean_w[lane],
                ] {
                    residuals[lane].add(-power)?;
                }
            }
            let source = ScatteringHeating::try_new(
                report,
                &bindings,
                thermal,
                self.model.tolerance,
                &mut *checkpoint,
            )?;
            for (&cell, &power) in source.sources_w() {
                checkpoint()?;
                sources.entry(cell).or_default().add(power)?;
            }
            accepted.push(source);
        }
        let tolerance = self.model.tolerance;
        if total_absorbed_error.value()? > tolerance.max_total_standard_error_w
            || total_unresolved_error.value()? > tolerance.max_total_standard_error_w
            || unresolved.value()? > tolerance.max_unresolved_w
        {
            return Err(invalid(
                "whole illumination sampling/truncation allowance exceeded",
            ));
        }
        for error in errors.values() {
            checkpoint()?;
            if error.value()? > tolerance.max_cell_standard_error_w {
                return Err(invalid(
                    "whole illumination receiving-cell sampling allowance exceeded",
                ));
            }
        }
        for residual in residuals {
            if residual.value()?.abs() > tolerance.power_roundoff_w {
                return Err(invalid(
                    "whole illumination power balance exceeds roundoff allowance",
                ));
            }
        }
        let mut sources_w = BTreeMap::new();
        let mut thermal_input = CheckedSum::default();
        for (cell, sum) in sources {
            checkpoint()?;
            let power = sum.value()?;
            thermal_input.add(power)?;
            sources_w.insert(cell, power);
        }
        if (thermal_input.value()? - absorbed.value()?).abs() > tolerance.power_roundoff_w {
            return Err(invalid(
                "combined illumination heat disagrees with original absorption",
            ));
        }
        checkpoint()?;
        Ok((
            PreparedIllumination {
                sources: accepted,
                sources_w,
                retained_cells,
            },
            traced,
        ))
    }
}

impl<const N: usize> DepositionObserver for TemperatureIlluminationHeating<'_, N> {
    type Error = OpticalFeedbackError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        self.segments = 0;
        self.builds = 0;
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
            return Err(OpticalFeedbackError::Budget);
        }
        let reusable = self.cache_valid;
        self.cache_valid = false;
        let changed = self
            .model
            .field
            .prepare(thermal, &mut self.coefficients, checkpoint)?;
        if !reusable || changed || self.current.is_none() {
            let (prepared, traced) = self.prepare_sources(thermal, checkpoint)?;
            let segments = self
                .segments
                .checked_add(traced)
                .ok_or(OpticalFeedbackError::Budget)?;
            let builds = self
                .builds
                .checked_add(1)
                .ok_or(OpticalFeedbackError::Budget)?;
            if segments > self.budget.max_traced_segments {
                return Err(OpticalFeedbackError::Budget);
            }
            self.current = Some(Arc::new(prepared));
            self.segments = segments;
            self.builds = builds;
        }
        checkpoint()?;
        self.cache_valid = true;
        self.evaluations += 1;
        Ok(self
            .current
            .as_ref()
            .expect("prepared illumination")
            .sources_w
            .clone())
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        self.model
            .field
            .validate_temperatures(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || thermal.time_s() <= start_s
            || !self.cache_valid
        {
            return Err(invalid(
                "matching prepared illumination and accepted native interval required",
            ));
        }
        self.model
            .field
            .validate_temperatures(thermal, checkpoint)?;
        let current = self
            .current
            .as_ref()
            .ok_or_else(|| invalid("native step needs prepared illumination"))?;
        let mut sums = self.sums.clone();
        sums.add_duration(duration_s)?;
        for source in &current.sources {
            checkpoint()?;
            sums.add_source(&source.interval_energy(duration_s)?)?;
        }
        let energy = sums.value()?;
        checkpoint()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn validate_source<const N: usize>(
    beam: &IncidentScatteringBeam<N>,
    field: &TemperatureOpticalField<'_, N>,
) -> Result<(), OpticalFeedbackError> {
    if beam.wavelengths_nm.map(f64::to_bits) != field.wavelengths_nm().map(f64::to_bits)
        || !beam.incident_power_w.iter().copied().all(nonnegative)
    {
        return Err(invalid(
            "illumination requires supported field wavelengths and finite nonnegative powers",
        ));
    }
    if [
        beam.origin.x,
        beam.origin.y,
        beam.origin.z,
        beam.unit_direction.x,
        beam.unit_direction.y,
        beam.unit_direction.z,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || (beam.unit_direction.dot(beam.unit_direction) - 1.0).abs() > 2.0e-10
    {
        return Err(VolumeTransportError::Ray.into());
    }
    Ok(())
}
fn invalid(message: &'static str) -> OpticalFeedbackError {
    OpticalFeedbackError::Invalid(message)
}
