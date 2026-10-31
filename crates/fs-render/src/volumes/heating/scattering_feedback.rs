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

use super::feedback::{
    CompiledTable, OpticalFeedbackError, TemperatureAbsorptionTable, TemperatureSpectrum,
};
use super::scattering::{ScatteringEnergyBudget, ScatteringHeating, ScatteringHeatingTolerance};
use super::{CheckedSum, OpticalHeatingError, nonnegative};
use crate::volumes::transport::scattering::deposition::{
    ScatteringDepositionConfig, ScatteringPowerDeposition,
};
use crate::volumes::transport::scattering::{
    CellPhaseFunction, CellScatteringGrid, MassScatteringSpectrum, PhaseFunction,
};
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Source-bound mass scattering in m²/kg versus kelvin. The retained rows are
/// actual scattering spectra, including their original source support and any
/// analytic wavelength law. Interpolation in kelvin reuses the absorption
/// feedback owner's compiled coefficient interpolation and never extrapolates.
#[derive(Clone, Debug)]
pub struct TemperatureScatteringTable {
    source_identity: ContentHash,
    rows: Vec<(f64, MassScatteringSpectrum)>,
    wavelength_bounds_nm: (f64, f64),
}

impl TemperatureScatteringTable {
    pub fn try_new(
        source_identity: ContentHash,
        rows: &[(f64, MassScatteringSpectrum)],
        max_temperature_knots: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        if rows.len() < 2
            || rows.len() > max_temperature_knots
            || source_identity == ContentHash([0; 32])
            || rows.iter().any(|(t, _)| !t.is_finite() || *t <= 0.0)
            || rows.windows(2).any(|p| p[0].0 >= p[1].0)
        {
            return Err(invalid(
                "ordered positive scattering temperatures, bounded rows and source identity required",
            ));
        }
        let bounds = rows[0].1.wavelength_bounds_nm();
        if rows.iter().any(|(_, s)| s.wavelength_bounds_nm() != bounds) {
            return Err(invalid(
                "scattering temperature rows must share a wavelength domain",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        retained.extend_from_slice(rows);
        Ok(Self {
            source_identity,
            rows: retained,
            wavelength_bounds_nm: bounds,
        })
    }

    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        (self.rows[0].0, self.rows[self.rows.len() - 1].0)
    }
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }
    pub fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
}

impl TemperatureSpectrum for MassScatteringSpectrum {
    fn source_identity(&self) -> Option<ContentHash> {
        MassScatteringSpectrum::source_identity(self)
    }
    fn at(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        MassScatteringSpectrum::at(self, wavelength_nm)
    }
}

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

/// Local material laws, heat receiver and phase-palette selection. Multiple
/// optical cells may share one thermal cell; all their true heat is summed.
#[derive(Clone, Copy, Debug)]
pub struct ThermalScatteringCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub absorption_table: usize,
    pub scattering_table: usize,
    pub phase: usize,
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

#[derive(Clone, Copy)]
struct BoundMaterial {
    thermal_cell: CellId,
    absorption_table: usize,
    scattering_table: usize,
}

struct FeedbackModel<'g, const N: usize> {
    identity: ContentHash,
    geometry: AbsorbingGrid<'g>,
    beam: IncidentScatteringBeam<N>,
    bounds_nm: (f64, f64),
    absorption: Vec<CompiledTable<N>>,
    scattering: Vec<CompiledTable<N>>,
    materials: Vec<BoundMaterial>,
    cell_materials: Vec<Option<usize>>,
    phases: Vec<PhaseFunction>,
    cell_phases: Vec<usize>,
    config: ScatteringDepositionConfig,
    tolerance: ScatteringHeatingTolerance,
    has_scattering: bool,
}

#[derive(Clone, Debug)]
struct PreparedCoefficients<const N: usize> {
    absorption: Vec<[f64; N]>,
    scattering: Vec<[f64; N]>,
    materials: Vec<([f64; N], [f64; N])>,
}

/// Cloneable transactional observer sharing immutable geometry, original source
/// data, compiled spectral samples and bindings. Mutable coefficient packets and
/// physical accounts are candidate-local, so native rollback includes photons.
#[derive(Clone)]
pub struct TemperatureScatteringHeating<'g, const N: usize> {
    model: Arc<FeedbackModel<'g, N>>,
    coefficients: PreparedCoefficients<N>,
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
            .field("model_identity", &self.model.identity)
            .field("time_s", &self.time_s)
            .field("histories_per_lane", &self.histories)
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
            || budget.max_source_evaluations == 0
            || budget.max_traced_segments == 0
            || config.histories_per_lane < 2
            || config.max_collisions_per_history == 0
            || config.max_total_segments == 0
        {
            return Err(OpticalFeedbackError::Budget);
        }
        if absorption.is_empty()
            || scattering.is_empty()
            || phases.is_empty()
            || thermal.cells().is_empty()
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
                "nonempty optical laws, active thermal cells and finite power/error allowances required",
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
        let mut bounds = (0.0_f64, f64::MAX);
        let mut knots = 0usize;
        let mut compiled_absorption = Vec::new();
        let mut compiled_scattering = Vec::new();
        compiled_absorption
            .try_reserve_exact(absorption.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        compiled_scattering
            .try_reserve_exact(scattering.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        hash.update(&(absorption.len() as u64).to_le_bytes());
        for table in absorption {
            checkpoint()?;
            knots = knots
                .checked_add(table.knot_count())
                .ok_or(OpticalFeedbackError::Budget)?;
            if knots > budget.max_total_temperature_knots {
                return Err(OpticalFeedbackError::Budget);
            }
            let domain = table.wavelength_bounds_nm();
            bounds.0 = bounds.0.max(domain.0);
            bounds.1 = bounds.1.min(domain.1);
            compiled_absorption.push(table.compile(
                &beam.wavelengths_nm,
                &mut hash,
                &mut checkpoint,
            )?);
        }
        hash.update(&(scattering.len() as u64).to_le_bytes());
        for table in scattering {
            checkpoint()?;
            knots = knots
                .checked_add(table.rows.len())
                .ok_or(OpticalFeedbackError::Budget)?;
            if knots > budget.max_total_temperature_knots {
                return Err(OpticalFeedbackError::Budget);
            }
            bounds.0 = bounds.0.max(table.wavelength_bounds_nm.0);
            bounds.1 = bounds.1.min(table.wavelength_bounds_nm.1);
            compiled_scattering.push(CompiledTable::compile(
                table.source_identity,
                table.wavelength_bounds_nm,
                &table.rows,
                &beam.wavelengths_nm,
                &mut hash,
                &mut checkpoint,
            )?);
        }
        let mut selected = allocate(count, None)?;
        for binding in bindings {
            checkpoint()?;
            if binding.grid_cell >= count
                || binding.absorption_table >= absorption.len()
                || binding.scattering_table >= scattering.len()
                || binding.phase >= phases.len()
            {
                return Err(invalid(
                    "optical binding references an unknown cell, table or phase",
                ));
            }
            if selected[binding.grid_cell].replace(*binding).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell).into());
            }
            if !thermal.cells().contains_key(&binding.thermal_cell) {
                return Err(OpticalHeatingError::InactiveCell(binding.thermal_cell).into());
            }
        }
        let mut materials = Vec::new();
        materials
            .try_reserve_exact(bindings.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        let mut unique_materials = BTreeMap::new();
        let mut cell_materials = allocate(count, None)?;
        let mut cell_phases = allocate(count, 0usize)?;
        for (index, binding) in selected.into_iter().enumerate() {
            checkpoint()?;
            let Some(binding) = binding else {
                if geometry
                    .cell_density_kg_m3(index)
                    .is_some_and(|rho| rho > 0.0)
                {
                    return Err(OpticalHeatingError::UnmappedCell(index).into());
                }
                continue;
            };
            let key = (
                binding.thermal_cell,
                binding.absorption_table,
                binding.scattering_table,
            );
            let slot = *unique_materials.entry(key).or_insert_with(|| {
                let slot = materials.len();
                materials.push(BoundMaterial {
                    thermal_cell: binding.thermal_cell,
                    absorption_table: binding.absorption_table,
                    scattering_table: binding.scattering_table,
                });
                slot
            });
            cell_materials[index] = Some(slot);
            cell_phases[index] = binding.phase;
            for value in [
                index as u64,
                binding.thermal_cell.0,
                binding.absorption_table as u64,
                binding.scattering_table as u64,
                binding.phase as u64,
            ] {
                hash.update(&value.to_le_bytes());
            }
        }
        let mut retained_phases = Vec::new();
        retained_phases
            .try_reserve_exact(phases.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        hash.update(&(phases.len() as u64).to_le_bytes());
        for phase in phases {
            checkpoint()?;
            phase.hash_into(&mut hash);
            retained_phases.push(phase.clone());
        }
        let has_scattering = compiled_scattering
            .iter()
            .any(CompiledTable::has_nonzero_coefficient);
        let mut coefficients = PreparedCoefficients {
            absorption: allocate(count, [0.0; N])?,
            scattering: allocate(count, [0.0; N])?,
            materials: allocate(materials.len(), ([0.0; N], [0.0; N]))?,
        };
        let model = Arc::new(FeedbackModel {
            identity: hash.finalize(),
            geometry,
            beam,
            bounds_nm: bounds,
            absorption: compiled_absorption,
            scattering: compiled_scattering,
            materials,
            cell_materials,
            phases: retained_phases,
            cell_phases,
            config,
            tolerance,
            has_scattering,
        });
        model.prepare(thermal, &mut coefficients, &mut checkpoint)?;
        // Validate full packet/phase support before exposing the observer,
        // without consuming any of the future trajectory's photon streams.
        transport_checkpoint(&mut checkpoint, |check| {
            model.field(&coefficients, check).map(|_| ())
        })?;
        let sums = EnergySums::default();
        let energy = sums.value()?;
        Ok(Self {
            model,
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
    /// Conservative over all admitted temperatures; an initially transparent
    /// scattering law may become active later along the thermal trajectory.
    pub fn has_scattering(&self) -> bool {
        self.model.has_scattering
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

impl<const N: usize> FeedbackModel<'_, N> {
    fn prepare(
        &self,
        thermal: &DepositionSimulation<'_>,
        fields: &mut PreparedCoefficients<N>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalFeedbackError> {
        for (material, values) in self.materials.iter().zip(&mut fields.materials) {
            checkpoint()?;
            let temperature = thermal.temperature_k(material.thermal_cell)?;
            values.0 = self.absorption[material.absorption_table]
                .at(material.thermal_cell, temperature)?;
            values.1 = self.scattering[material.scattering_table]
                .at(material.thermal_cell, temperature)?;
        }
        for (index, slot) in self.cell_materials.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Some(slot) = slot {
                (fields.absorption[index], fields.scattering[index]) = fields.materials[*slot];
            }
        }
        Ok(())
    }

    fn validate_temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalFeedbackError> {
        for material in &self.materials {
            checkpoint()?;
            let t = thermal.temperature_k(material.thermal_cell)?;
            self.absorption[material.absorption_table]
                .check_temperature(material.thermal_cell, t)?;
            self.scattering[material.scattering_table]
                .check_temperature(material.thermal_cell, t)?;
        }
        Ok(())
    }
}

impl<'g, const N: usize> FeedbackModel<'g, N> {
    fn field<'s>(
        &'s self,
        values: &'s PreparedCoefficients<N>,
        check: &mut dyn FnMut() -> Result<(), VolumeTransportError>,
    ) -> Result<CellScatteringGrid<'s, 'g, N>, VolumeTransportError> {
        CellScatteringGrid::try_new_with_checkpoint(
            &mut || check(),
            &self.geometry,
            self.beam.wavelengths_nm,
            &values.absorption,
            &values.scattering,
            CellPhaseFunction::Indexed {
                palette: &self.phases,
                cell_indices: &self.cell_phases,
            },
            self.bounds_nm,
            self.identity,
        )
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
            .prepare(thermal, &mut self.coefficients, checkpoint)?;
        let config = self.step_config();
        let report = transport_checkpoint(checkpoint, |check| {
            let field = self.model.field(&self.coefficients, check)?;
            field.deposit_multiple_scattering_power_with_checkpoint(
                &mut || check(),
                self.model.beam.origin,
                self.model.beam.unit_direction,
                &self.model.beam.incident_power_w,
                config,
            )
        })?;
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
            let material = self
                .model
                .cell_materials
                .get(cell.cell_index)
                .copied()
                .flatten()
                .ok_or(OpticalHeatingError::UnmappedCell(cell.cell_index))?;
            bindings.push((cell.cell_index, self.model.materials[material].thermal_cell));
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
        self.model.validate_temperatures(thermal, &mut || Ok(()))
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
        self.model.validate_temperatures(thermal, checkpoint)?;
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
struct EnergySums<const N: usize> {
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
        self.duration.add(interval.duration_s)?;
        self.thermal.add(interval.thermal_input_j)?;
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
    fn value(&self) -> Result<ScatteringFeedbackEnergyBudget<N>, OpticalFeedbackError> {
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
}

fn allocate<T: Clone>(count: usize, value: T) -> Result<Vec<T>, OpticalFeedbackError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| OpticalFeedbackError::Budget)?;
    result.resize(count, value);
    Ok(result)
}

fn transport_checkpoint<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(
        &mut dyn FnMut() -> Result<(), VolumeTransportError>,
    ) -> Result<T, VolumeTransportError>,
) -> Result<T, OpticalFeedbackError> {
    let mut refusal = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => Ok(()),
        Err(error) => {
            refusal = Some(error);
            Err(VolumeTransportError::Cancelled)
        }
    });
    if let Some(error) = refusal {
        return Err(error.into());
    }
    result.map_err(Into::into)
}

fn invalid(reason: &'static str) -> OpticalFeedbackError {
    OpticalFeedbackError::Invalid(reason)
}
