//! Temperature-dependent multiple-scattering heat on the native thermal clock.
//!
//! Local absorption and scattering are evaluated from the original enthalpy
//! owner's temperatures, then the shared exact-cell photon transport supplies
//! the next explicit thermal substep. Geometry, density, source power, phase
//! laws and spectral quadrature stay prescribed. Every positive-density cell
//! has a material/thermal binding, including cells outside the original beam.
//!
//! Independent physical-start-time streams are used within each trajectory.
//! A source report's empirical error is conditional on that step's temperature
//! state; it is not an error bar on the final nonlinear thermal state. Cumulative
//! energies below account for the actual sampled heat, without incorrectly
//! applying a frozen-field standard-error formula to changing optical fields.
//!
//! Absorption estimation is an immutable trajectory choice. Collision capture
//! remains the default; continuous absorption scores every travelled segment
//! while sampling scattering alone. Both choices use these same native thermal
//! substeps, material tables, work bounds, and heat-admission policy.

use super::feedback::{OpticalFeedbackError, TemperatureAbsorptionTable};
use super::scattering::{ScatteringEnergyBudget, ScatteringHeating, ScatteringHeatingTolerance};
use super::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget, TemperatureOpticalState,
    transport_checkpoint_with_error,
};
use super::{CheckedSum, OpticalHeatingError, nonnegative};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::volumes::transport::scattering::PhaseFunction;
use crate::volumes::transport::scattering::deposition::{
    AbsorptionEstimator, PolarizedDepositionConfig, ScatteringDepositionConfig,
    ScatteringPowerDeposition,
};
use crate::volumes::transport::scattering::polarized::PolarizedScatteringError;
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};
use fs_blake3::DomainHasher;
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;
use std::sync::Arc;

mod refinement;
pub use refinement::{
    ScatteringCellEstimate, ScatteringEnsembleAssessment, ScatteringEnsembleEnergy,
    ScatteringEnsembleTolerance, ScatteringEnsembleWork, ScatteringFeedbackEnsemble,
    ScatteringFeedbackRefinement, ScatteringFeedbackTrajectory, ScatteringRefinementDiagnostics,
    ScatteringRefinementError, ScatteringScalarEstimate,
};

pub use super::temperature_field::{TemperatureScatteringTable, ThermalScatteringCell};

/// A collimated source followed until every represented packet escapes the
/// open grid or is retained unresolved. Powers are quadrature-weighted watts.
/// This deliberately has no finite directed-segment cutoff: scattered photons
/// may reach any cell of the original geometry.
#[derive(Clone, Copy, Debug)]
pub struct IncidentScatteringBeam<const N: usize> {
    pub origin: Point3,
    pub unit_direction: Vec3,
    pub wavelengths_nm: [f64; N],
    pub incident_power_w: [f64; N],
}

/// Full-Stokes collimated source with spectral-quadrature-weighted watts.
/// Polarization is referred to the common physical frame. Absorption remains
/// isotropic; the admitted molecular/spherical phase law supplies the angular
/// response, including polarization developed during later scattering events.
#[derive(Clone, Copy, Debug)]
pub struct IncidentPolarizedScatteringBeam<const N: usize> {
    pub origin: Point3,
    pub frame: PolarizationFrame,
    pub wavelengths_nm: [f64; N],
    pub incident_stokes_w: [Stokes; N],
}

#[derive(Clone, Copy)]
struct PolarizedSource<const N: usize> {
    beam: IncidentPolarizedScatteringBeam<N>,
    max_phase_proposals: u32,
}

/// Storage bounds at admission and aggregate work ceilings for ONE native
/// advance, including all its substeps. The transport config also limits every
/// individual source evaluation; neither limit is silently renewed per step.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringFeedbackBudget {
    pub max_grid_cells: usize,
    /// Combined count of absorption and scattering temperature tables.
    pub max_tables: usize,
    pub max_total_temperature_knots: usize,
    pub max_source_evaluations: usize,
    pub max_traced_segments: u64,
}

/// Actual radiation and heat accumulated over accepted native substeps. These
/// are sampled trajectory accounts, not ensemble means or final-state error
/// bars. Optical heat is already included once in the native external-heat audit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringFeedbackEnergyBudget<const N: usize> {
    pub duration_s: f64,
    pub incident_j: [f64; N],
    pub ballistic_escape_j: [f64; N],
    pub diffuse_escape_j: [f64; N],
    pub absorbed_j: [f64; N],
    pub unresolved_j: [f64; N],
    pub residual_j: [f64; N],
    pub thermal_input_j: f64,
}

/// Largest conditional source standard errors encountered on this trajectory.
/// Units stay watts. These diagnostics do NOT propagate nonlinear temperature
/// feedback, temporal error or unresolved radiation into final uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConditionalSourcePrecision<const N: usize> {
    pub maximum_absorbed_standard_error_w: [f64; N],
    pub maximum_unresolved_standard_error_w: [f64; N],
}

struct FeedbackModel<'g, const N: usize> {
    field: Arc<TemperatureOpticalField<'g, N>>,
    beam: IncidentScatteringBeam<N>,
    config: ScatteringDepositionConfig,
    absorption_estimator: AbsorptionEstimator,
    polarized: Option<PolarizedSource<N>>,
    tolerance: ScatteringHeatingTolerance,
}

/// Cloneable transactional observer sharing immutable geometry, original source
/// data, compiled spectral samples and bindings. Mutable coefficient packets and
/// physical accounts are candidate-local, so native rollback includes photons.
#[derive(Clone)]
pub struct TemperatureScatteringHeating<'g, const N: usize> {
    model: Arc<FeedbackModel<'g, N>>,
    coefficients: TemperatureOpticalState<N>,
    current: Option<ScatteringHeating<N>>,
    budget: ScatteringFeedbackBudget,
    histories: u32,
    replica_seed: u64,
    time_s: f64,
    evaluations: usize,
    segments: u64,
    sums: EnergySums<N>,
    energy: ScatteringFeedbackEnergyBudget<N>,
    precision: ConditionalSourcePrecision<N>,
}

impl<const N: usize> core::fmt::Debug for TemperatureScatteringHeating<'_, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("TemperatureScatteringHeating")
            .field("model_identity", &self.model.field.identity())
            .field("time_s", &self.time_s)
            .field("histories_per_lane", &self.histories)
            .field("absorption_estimator", &self.model.absorption_estimator)
            .field("polarized", &self.is_polarized())
            .field("source_evaluations", &self.evaluations)
            .field("traced_segments", &self.segments)
            .field("energy", &self.energy)
            .finish_non_exhaustive()
    }
}

impl<'g, const N: usize> TemperatureScatteringHeating<'g, N> {
    /// Compile source spectra once and bind every nonvacuum original grid cell.
    /// The geometry is moved without copying its borrowed density. Its frozen
    /// absorption/LTE temperature fields do not replace the supplied local laws
    /// or native thermal temperatures. No initial photon draw is taken here.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        geometry: AbsorbingGrid<'g>,
        beam: IncidentScatteringBeam<N>,
        absorption: &[TemperatureAbsorptionTable],
        scattering: &[TemperatureScatteringTable],
        phases: &[PhaseFunction],
        bindings: &[ThermalScatteringCell],
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        tolerance: ScatteringHeatingTolerance,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        Self::try_new_with_estimator(
            geometry,
            beam,
            absorption,
            scattering,
            phases,
            bindings,
            thermal,
            budget,
            config,
            AbsorptionEstimator::Collision,
            tolerance,
            checkpoint,
        )
    }

    /// Compile the native temperature feedback with an explicit absorption
    /// estimator. Continuous absorption retains Beer heat on each traversed cell,
    /// including histories without a scattering event. It can reduce variance
    /// for weak absorption; neither choice promises lower variance in all media.
    /// The choice is retained by clones, replica streams and trajectory refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_estimator(
        geometry: AbsorbingGrid<'g>,
        beam: IncidentScatteringBeam<N>,
        absorption: &[TemperatureAbsorptionTable],
        scattering: &[TemperatureScatteringTable],
        phases: &[PhaseFunction],
        bindings: &[ThermalScatteringCell],
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        absorption_estimator: AbsorptionEstimator,
        tolerance: ScatteringHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        let count = geometry.cell_count();
        let table_count = absorption
            .len()
            .checked_add(scattering.len())
            .ok_or(OpticalFeedbackError::Budget)?;
        if N == 0
            || count > budget.max_grid_cells
            || table_count > budget.max_tables
            || bindings.len() > count
            || phases.len() > count
        {
            return Err(OpticalFeedbackError::Budget);
        }
        validate_source(&beam, thermal, budget, config, tolerance)?;
        let components = [
            beam.origin.x,
            beam.origin.y,
            beam.origin.z,
            beam.unit_direction.x,
            beam.unit_direction.y,
            beam.unit_direction.z,
        ];
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-scattering-model.v1");
        hash.update(geometry.content_identity().as_bytes());
        hash.update(&(N as u64).to_le_bytes());
        for value in components {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for lane in 0..N {
            hash.update(&beam.wavelengths_nm[lane].to_bits().to_le_bytes());
            hash.update(&beam.incident_power_w[lane].to_bits().to_le_bytes());
        }
        let (field, coefficients) = TemperatureOpticalField::compile_with_hash(
            geometry,
            beam.wavelengths_nm,
            absorption,
            scattering,
            phases,
            bindings,
            thermal,
            field_budget(budget),
            hash,
            &mut checkpoint,
        )?;
        Self::assembled(
            Arc::new(field),
            coefficients,
            beam,
            thermal,
            budget,
            config,
            absorption_estimator,
            None,
            tolerance,
        )
    }

    /// Use the SAME compiled constitutive field as other optical/thermal
    /// observers. Candidate coefficient packets and sampled accounts remain
    /// private; immutable source tables, geometry, and phase data share one Arc.
    #[allow(clippy::too_many_arguments)]
    pub fn from_field(
        field: Arc<TemperatureOpticalField<'g, N>>,
        beam: IncidentScatteringBeam<N>,
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        tolerance: ScatteringHeatingTolerance,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        Self::from_field_with_estimator(
            field,
            beam,
            thermal,
            budget,
            config,
            AbsorptionEstimator::Collision,
            tolerance,
            checkpoint,
        )
    }

    /// Share an admitted constitutive field while choosing the immutable
    /// trajectory estimator. The field identity continues to describe physical
    /// source data; the deposition report separately records this estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn from_field_with_estimator(
        field: Arc<TemperatureOpticalField<'g, N>>,
        beam: IncidentScatteringBeam<N>,
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        absorption_estimator: AbsorptionEstimator,
        tolerance: ScatteringHeatingTolerance,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        Self::from_field_source(
            field,
            beam,
            None,
            thermal,
            budget,
            config,
            absorption_estimator,
            tolerance,
            checkpoint,
        )
    }

    /// Couple an admitted full-Stokes source to the SAME changing native field.
    ///
    /// Original Stokes/frame admission and active phase matrices are checked
    /// before any history is sampled. Every native source query then prepares
    /// the accepted temperatures and carries polarization through each actual
    /// scattering event. Isotropic absorption changes watts, not normalized
    /// Stokes state. The source, proposal bound and estimator are immutable
    /// across clones, replica streams, history changes and trajectory refinement.
    ///
    /// A scalar phase without a retained polarized response refuses whenever
    /// its prepared cell has scattering; pure absorbers need no such response.
    #[allow(clippy::too_many_arguments)]
    pub fn from_polarized_field(
        field: Arc<TemperatureOpticalField<'g, N>>,
        beam: IncidentPolarizedScatteringBeam<N>,
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: PolarizedDepositionConfig,
        absorption_estimator: AbsorptionEstimator,
        tolerance: ScatteringHeatingTolerance,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        let scalar = IncidentScatteringBeam {
            origin: beam.origin,
            unit_direction: beam.frame.direction(),
            wavelengths_nm: beam.wavelengths_nm,
            incident_power_w: beam.incident_stokes_w.map(|s| s.i),
        };
        Self::from_field_source(
            field,
            scalar,
            Some(PolarizedSource {
                beam,
                max_phase_proposals: config.max_phase_proposals,
            }),
            thermal,
            budget,
            config.transport,
            absorption_estimator,
            tolerance,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_field_source(
        field: Arc<TemperatureOpticalField<'g, N>>,
        beam: IncidentScatteringBeam<N>,
        polarized: Option<PolarizedSource<N>>,
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        absorption_estimator: AbsorptionEstimator,
        tolerance: ScatteringHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        if field.geometry().cell_count() > budget.max_grid_cells
            || field.table_count() > budget.max_tables
            || field.temperature_knot_count() > budget.max_total_temperature_knots
        {
            return Err(OpticalFeedbackError::Budget);
        }
        validate_source(&beam, thermal, budget, config, tolerance)?;
        if beam.wavelengths_nm.map(f64::to_bits) != field.wavelengths_nm().map(f64::to_bits) {
            return Err(invalid(
                "incident wavelengths must match the compiled optical field",
            ));
        }
        let coefficients = field.initial_state(thermal, &mut checkpoint)?;
        transport_checkpoint_with_error(
            &mut checkpoint,
            |check| -> Result<(), PolarizedScatteringError> {
                let prepared = field.field(&coefficients, check)?;
                if let Some(source) = polarized {
                    prepared.admit_polarized_source(
                        check,
                        source.beam.origin,
                        source.beam.frame,
                        &source.beam.incident_stokes_w,
                        source.max_phase_proposals,
                    )?;
                }
                Ok(())
            },
        )?;
        Self::assembled(
            field,
            coefficients,
            beam,
            thermal,
            budget,
            config,
            absorption_estimator,
            polarized,
            tolerance,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn assembled(
        field: Arc<TemperatureOpticalField<'g, N>>,
        coefficients: TemperatureOpticalState<N>,
        beam: IncidentScatteringBeam<N>,
        thermal: &DepositionSimulation<'_>,
        budget: ScatteringFeedbackBudget,
        config: ScatteringDepositionConfig,
        absorption_estimator: AbsorptionEstimator,
        polarized: Option<PolarizedSource<N>>,
        tolerance: ScatteringHeatingTolerance,
    ) -> Result<Self, OpticalFeedbackError> {
        let sums = EnergySums::default();
        let energy = sums.value()?;
        Ok(Self {
            model: Arc::new(FeedbackModel {
                field,
                beam,
                config,
                absorption_estimator,
                polarized,
                tolerance,
            }),
            coefficients,
            current: None,
            budget,
            histories: config.histories_per_lane,
            replica_seed: 0,
            time_s: thermal.time_s(),
            evaluations: 0,
            segments: 0,
            sums,
            energy,
            precision: ConditionalSourcePrecision {
                maximum_absorbed_standard_error_w: [0.0; N],
                maximum_unresolved_standard_error_w: [0.0; N],
            },
        })
    }

    pub fn optical_field(&self) -> &Arc<TemperatureOpticalField<'g, N>> {
        &self.model.field
    }

    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn energy(&self) -> &ScatteringFeedbackEnergyBudget<N> {
        &self.energy
    }
    /// Most recent conditional field, used over the last accepted interval.
    /// Before the first actual native source query this is absent.
    pub fn deposition(&self) -> Option<&ScatteringPowerDeposition<N>> {
        self.current.as_ref().map(ScatteringHeating::deposition)
    }
    pub fn sources_w(&self) -> Option<&BTreeMap<CellId, f64>> {
        self.current.as_ref().map(ScatteringHeating::sources_w)
    }
    pub fn conditional_source_precision(&self) -> &ConditionalSourcePrecision<N> {
        &self.precision
    }
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn traced_segments(&self) -> u64 {
        self.segments
    }
    pub fn budget(&self) -> ScatteringFeedbackBudget {
        self.budget
    }
    pub fn histories_per_lane(&self) -> u32 {
        self.histories
    }
    /// Estimator selected for every source draw in this thermal trajectory.
    pub fn absorption_estimator(&self) -> AbsorptionEstimator {
        self.model.absorption_estimator
    }
    /// Whether native histories carry the admitted full-Stokes source state.
    pub fn is_polarized(&self) -> bool {
        self.model.polarized.is_some()
    }
    /// Conservative over all admitted temperatures; an initially transparent
    /// scattering law may become active later along the thermal trajectory.
    pub fn has_scattering(&self) -> bool {
        self.model.field.has_scattering()
    }

    /// Start a future conditional replica from this accepted physical state.
    /// Past heat is retained; past uncertainty is not resampled by this clone.
    pub fn with_replica_seed(&self, replica_seed: u64) -> Self {
        let mut result = self.clone();
        result.replica_seed = replica_seed;
        result
    }
    /// Change fixed histories for future source evaluations. History count is
    /// absent from stream keys, so resolution studies retain common prefixes.
    pub fn with_histories_per_lane(&self, histories: u32) -> Result<Self, OpticalFeedbackError> {
        if histories < 2 {
            return Err(OpticalFeedbackError::Budget);
        }
        let mut result = self.clone();
        result.histories = histories;
        Ok(result)
    }
    /// Tighten, never replenish, the caller's per-advance optical work ceilings.
    pub fn with_work_limits(&self, evaluations: usize, segments: u64) -> Self {
        let mut result = self.clone();
        result.budget.max_source_evaluations = evaluations.min(self.budget.max_source_evaluations);
        result.budget.max_traced_segments = segments.min(self.budget.max_traced_segments);
        result
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), OpticalFeedbackError> {
        if self.time_s != thermal.time_s() {
            return Err(invalid(
                "scattering-feedback and native thermal clocks differ",
            ));
        }
        Ok(())
    }

    fn step_config(&self) -> ScatteringDepositionConfig {
        let mut hash =
            DomainHasher::new("org.frankensim.render.temperature-scattering-step-stream.v1");
        // Temperature/coefficient state and history count must not key streams:
        // conditional fresh draws and paired studies depend on this separation.
        for value in [
            self.model.config.seed,
            self.model.config.source_id,
            self.replica_seed,
            if self.time_s == 0.0 {
                0
            } else {
                self.time_s.to_bits()
            },
        ] {
            hash.update(&value.to_le_bytes());
        }
        let identity = hash.finalize();
        let bytes = identity.as_bytes();
        let mut config = self.model.config;
        config.seed = u64::from_le_bytes(bytes[..8].try_into().expect("eight hash bytes"));
        config.source_id = u64::from_le_bytes(bytes[8..16].try_into().expect("eight hash bytes"));
        config.histories_per_lane = self.histories;
        config.max_total_segments = config.max_total_segments.min(
            self.budget
                .max_traced_segments
                .saturating_sub(self.segments),
        );
        config
    }
}

impl<const N: usize> DepositionObserver for TemperatureScatteringHeating<'_, N> {
    type Error = OpticalFeedbackError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        self.segments = 0;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        checkpoint()?;
        if self.evaluations >= self.budget.max_source_evaluations
            || self.segments >= self.budget.max_traced_segments
        {
            return Err(OpticalFeedbackError::Budget);
        }
        self.model
            .field
            .prepare(thermal, &mut self.coefficients, checkpoint)?;
        let config = self.step_config();
        let report = transport_checkpoint_with_error(
            checkpoint,
            |check| -> Result<_, PolarizedScatteringError> {
                let field = self.model.field.field(&self.coefficients, check)?;
                if let Some(source) = self.model.polarized {
                    return field.deposit_polarized_power_with_checkpoint(
                        check,
                        source.beam.origin,
                        source.beam.frame,
                        &source.beam.incident_stokes_w,
                        PolarizedDepositionConfig {
                            transport: config,
                            max_phase_proposals: source.max_phase_proposals,
                        },
                        self.model.absorption_estimator,
                    );
                }
                match self.model.absorption_estimator {
                    AbsorptionEstimator::Collision => field
                        .deposit_multiple_scattering_power_with_checkpoint(
                            &mut || check(),
                            self.model.beam.origin,
                            self.model.beam.unit_direction,
                            &self.model.beam.incident_power_w,
                            config,
                        ),
                    AbsorptionEstimator::Continuous => field
                        .deposit_continuous_absorption_power_with_checkpoint(
                            &mut || check(),
                            self.model.beam.origin,
                            self.model.beam.unit_direction,
                            &self.model.beam.incident_power_w,
                            config,
                        ),
                }
                .map_err(Into::into)
            },
        )?;
        let segments = self
            .segments
            .checked_add(report.traced_segments)
            .ok_or(OpticalFeedbackError::Budget)?;
        if segments > self.budget.max_traced_segments {
            return Err(OpticalFeedbackError::Budget);
        }
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
        }
        // Preserve one owner for conservation, uncertainty and many-to-one
        // thermal admission; no second report-validation policy is introduced.
        let accepted = ScatteringHeating::try_new(
            report,
            &bindings,
            thermal,
            self.model.tolerance,
            &mut *checkpoint,
        )?;
        let sources = accepted.sources_w().clone();
        self.current = Some(accepted);
        self.evaluations += 1;
        self.segments = segments;
        Ok(sources)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        // The native owner brackets every event with its cancellation callback.
        // Endpoint events still validate material support with no next step.
        self.model
            .field
            .validate_temperatures(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "matching accepted interval required for temperature-dependent scattered heat",
            ));
        }
        self.model
            .field
            .validate_temperatures(thermal, checkpoint)?;
        let current = self
            .current
            .as_mut()
            .ok_or_else(|| invalid("native step needs a prepared scattering source"))?;
        current.step(start_s, duration_s, old, thermal, checkpoint)?;
        let mut sums = self.sums.clone();
        sums.add(current.energy())?;
        let energy = sums.value()?;
        checkpoint()?;
        for lane in 0..N {
            self.precision.maximum_absorbed_standard_error_w[lane] =
                self.precision.maximum_absorbed_standard_error_w[lane]
                    .max(current.deposition().absorbed.standard_error_w[lane]);
            self.precision.maximum_unresolved_standard_error_w[lane] =
                self.precision.maximum_unresolved_standard_error_w[lane]
                    .max(current.deposition().unresolved.standard_error_w[lane]);
        }
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub(super) struct EnergySums<const N: usize> {
    duration: CheckedSum,
    incident: [CheckedSum; N],
    ballistic_escape: [CheckedSum; N],
    diffuse_escape: [CheckedSum; N],
    absorbed: [CheckedSum; N],
    unresolved: [CheckedSum; N],
    residual: [CheckedSum; N],
    thermal: CheckedSum,
}
impl<const N: usize> Default for EnergySums<N> {
    fn default() -> Self {
        Self {
            duration: CheckedSum::default(),
            incident: [CheckedSum::default(); N],
            ballistic_escape: [CheckedSum::default(); N],
            diffuse_escape: [CheckedSum::default(); N],
            absorbed: [CheckedSum::default(); N],
            unresolved: [CheckedSum::default(); N],
            residual: [CheckedSum::default(); N],
            thermal: CheckedSum::default(),
        }
    }
}
impl<const N: usize> EnergySums<N> {
    fn add(&mut self, interval: &ScatteringEnergyBudget<N>) -> Result<(), OpticalFeedbackError> {
        self.add_duration(interval.duration_s)?;
        self.add_source(interval)
    }
    pub(super) fn add_duration(&mut self, duration_s: f64) -> Result<(), OpticalFeedbackError> {
        self.duration.add(duration_s)?;
        Ok(())
    }
    /// Several simultaneously active sources share one accepted duration.
    pub(super) fn add_source(
        &mut self,
        interval: &ScatteringEnergyBudget<N>,
    ) -> Result<(), OpticalFeedbackError> {
        self.thermal.add(interval.thermal_input_j)?;
        self.add_source_radiation(interval)
    }
    /// Scheduled sources reduce their actual receiving-cell heat once, after
    /// combining all independently scaled spectral source reports.
    pub(super) fn add_thermal_input(&mut self, joules: f64) -> Result<(), OpticalFeedbackError> {
        self.thermal.add(joules)?;
        Ok(())
    }
    pub(super) fn add_source_radiation(
        &mut self,
        interval: &ScatteringEnergyBudget<N>,
    ) -> Result<(), OpticalFeedbackError> {
        for lane in 0..N {
            self.incident[lane].add(interval.incident_j[lane])?;
            self.ballistic_escape[lane].add(interval.ballistic_escape_j[lane])?;
            self.diffuse_escape[lane].add(interval.diffuse_escape_j[lane])?;
            self.absorbed[lane].add(interval.absorbed_j[lane])?;
            self.unresolved[lane].add(interval.unresolved_j[lane])?;
            self.residual[lane].add(interval.residual_j[lane])?;
        }
        Ok(())
    }
    pub(super) fn value(&self) -> Result<ScatteringFeedbackEnergyBudget<N>, OpticalFeedbackError> {
        let mut result = ScatteringFeedbackEnergyBudget {
            duration_s: self.duration.value()?,
            incident_j: [0.0; N],
            ballistic_escape_j: [0.0; N],
            diffuse_escape_j: [0.0; N],
            absorbed_j: [0.0; N],
            unresolved_j: [0.0; N],
            residual_j: [0.0; N],
            thermal_input_j: self.thermal.value()?,
        };
        for lane in 0..N {
            result.incident_j[lane] = self.incident[lane].value()?;
            result.ballistic_escape_j[lane] = self.ballistic_escape[lane].value()?;
            result.diffuse_escape_j[lane] = self.diffuse_escape[lane].value()?;
            result.absorbed_j[lane] = self.absorbed[lane].value()?;
            result.unresolved_j[lane] = self.unresolved[lane].value()?;
            result.residual_j[lane] = self.residual[lane].value()?;
        }
        Ok(result)
    }
    pub(super) fn energy_sum(&self, channel: usize, lane: usize) -> CheckedSum {
        match channel {
            0 => self.incident[lane],
            1 => self.ballistic_escape[lane],
            2 => self.diffuse_escape[lane],
            3 => self.absorbed[lane],
            4 => self.unresolved[lane],
            5 => self.residual[lane],
            6 => self.thermal,
            _ => unreachable!("private optical energy channel"),
        }
    }
}

fn field_budget(budget: ScatteringFeedbackBudget) -> TemperatureOpticalFieldBudget {
    TemperatureOpticalFieldBudget {
        max_grid_cells: budget.max_grid_cells,
        max_tables: budget.max_tables,
        max_total_temperature_knots: budget.max_total_temperature_knots,
    }
}

fn validate_source<const N: usize>(
    beam: &IncidentScatteringBeam<N>,
    thermal: &DepositionSimulation<'_>,
    budget: ScatteringFeedbackBudget,
    config: ScatteringDepositionConfig,
    tolerance: ScatteringHeatingTolerance,
) -> Result<(), OpticalFeedbackError> {
    if N == 0
        || budget.max_source_evaluations == 0
        || budget.max_traced_segments == 0
        || config.histories_per_lane < 2
        || config.max_collisions_per_history == 0
        || config.max_total_segments == 0
    {
        return Err(OpticalFeedbackError::Budget);
    }
    if thermal.cells().is_empty()
        || !beam.incident_power_w.iter().copied().all(nonnegative)
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
            "active thermal cells and finite power/error allowances required",
        ));
    }
    let components = [
        beam.origin.x,
        beam.origin.y,
        beam.origin.z,
        beam.unit_direction.x,
        beam.unit_direction.y,
        beam.unit_direction.z,
    ];
    if components.iter().any(|v| !v.is_finite())
        || (beam.unit_direction.dot(beam.unit_direction) - 1.0).abs() > 2.0e-10
    {
        return Err(VolumeTransportError::Ray.into());
    }
    Ok(())
}

fn invalid(reason: &'static str) -> OpticalFeedbackError {
    OpticalFeedbackError::Invalid(reason)
}
