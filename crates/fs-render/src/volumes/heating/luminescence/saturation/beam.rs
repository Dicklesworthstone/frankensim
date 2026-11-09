//! Depleting spectral pumps on the original cell and thermal timeline.
//!
//! One collimated, quadrature-weighted beam is transported in native DDA order.
//! Its molecular optical depth is cross section times the complete cell's
//! ground inventory times ray length / cell volume. The coupled thermal owner
//! solves that inventory and the transmitted photons together. This is an
//! explicitly cell-uniform finite-volume closure: strong longitudinal or
//! transverse saturation requires spatial refinement. No beam area is invented.
//! Molecular cross sections and kinetic rates refresh at the native thermal
//! interval; passive host coefficients retain the original optical snapshot.
//! Inner refinement controls population/opacity coupling, not temperature
//! feedback error.

use super::{
    CheckedSum, ContentHash, DepositionError, DepositionSimulation, DomainHasher,
    OpticalHeatingError, SaturableLuminescentHeating, finite, invalid, nonnegative,
    preserve_checkpoint, resolve_model,
};
use crate::volumes::heating::feedback::IncidentBeam;
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError, positive_product};
use fs_exec::Cx;
use fs_thermal::enthalpy::CellId;
use fs_thermal::luminescence::transfer::pump_column::{
    PumpColumnStep, SaturablePumpCell, SaturablePumpColumn,
};
use fs_thermal::luminescence::transfer::saturation::SaturableSingletTripletStep;
use std::{collections::BTreeMap, sync::Arc};

pub use fs_thermal::luminescence::transfer::pump_column::PumpColumnLimits;

#[derive(Clone, Debug)]
enum CrossSections {
    Authored {
        identity: ContentHash,
        values: [[f64; 2]; 2],
    },
    #[cfg(feature = "material-state")]
    Material {
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    },
}

/// Molecular cross sections for one original emitter cell. Rows are the two
/// incident wavelength lanes; columns are transitions [ground->singlet,
/// ground->triplet]. Both transitions compete for the SAME incident photons.
/// Bulk mass absorption is a separate passive host channel.
#[derive(Clone, Debug)]
pub struct BeamEmitterCrossSections {
    grid_cell: usize,
    source: CrossSections,
}

impl BeamEmitterCrossSections {
    /// Bind explicit, nonnegative square metres per ground-state emitter.
    /// The source identity attributes authored data; it does not certify it.
    pub fn authored(
        grid_cell: usize,
        source_identity: ContentHash,
        cross_sections_m2: [[f64; 2]; 2],
    ) -> Result<Self, OpticalHeatingError> {
        if source_identity == ContentHash([0; 32])
            || !cross_sections_m2.into_iter().flatten().all(nonnegative)
        {
            return Err(invalid("invalid molecular pump cross sections"));
        }
        Ok(Self {
            grid_cell,
            source: CrossSections::Authored {
                identity: source_identity,
                values: cross_sections_m2,
            },
        })
    }

    /// Bind the original selected molecular cross-section claims. Both vacuum
    /// wavelength and absolute temperature must be explicit owned coordinates;
    /// every other physical condition remains fixed by the native query.
    #[cfg(feature = "material-state")]
    pub fn from_material(
        grid_cell: usize,
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    ) -> Result<Self, OpticalHeatingError> {
        if source.temperature_axis().is_none() {
            return Err(invalid(
                "thermal pump source requires an owned absolute-temperature axis",
            ));
        }
        // Re-query at the exact initial point to reject a different card now.
        source
            .at_query_point(&card, &source.resolved().query())
            .map_err(|_| invalid("pump cross sections differ from original material card"))?;
        Ok(Self {
            grid_cell,
            source: CrossSections::Material { card, source },
        })
    }

    pub const fn grid_cell(&self) -> usize {
        self.grid_cell
    }

    fn identity(&self) -> ContentHash {
        match &self.source {
            CrossSections::Authored { identity, .. } => *identity,
            #[cfg(feature = "material-state")]
            CrossSections::Material { source, .. } => source.identity(),
        }
    }

    fn evaluate(
        &self,
        wavelengths_nm: [f64; 2],
        temperature_k: f64,
    ) -> Result<[[f64; 2]; 2], OpticalHeatingError> {
        match &self.source {
            CrossSections::Authored { values, .. } => {
                let _ = (wavelengths_nm, temperature_k);
                Ok(*values)
            }
            #[cfg(feature = "material-state")]
            CrossSections::Material { card, source } => {
                let mut values = [[0.0; 2]; 2];
                for lane in 0..2 {
                    values[lane] = source.at_wavelength_temperature(
                        card,
                        finite(wavelengths_nm[lane] / 1e9)?,
                        Some(temperature_k),
                    ).map_err(|_| invalid("molecular pump wavelength/temperature is outside original material support"))?
                        .cross_sections_m2();
                }
                Ok(values)
            }
        }
    }
}

/// Accepted cumulative pump ledger, in incident wavelength order. Molecular
/// absorption supplies excitation plus relaxation heat; passive host absorption
/// supplies heat directly. Neither emission nor a detector alters this ledger.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DepletingPumpEnergy {
    pub duration_s: f64,
    pub incident_j: [f64; 2],
    pub transmitted_j: [f64; 2],
    pub absorbed_j: [f64; 2],
    pub passive_absorbed_j: [f64; 2],
    /// Incident minus transmitted, molecular absorption and passive absorption.
    pub residual_j: [f64; 2],
    pub accepted_kinetic_substeps: u64,
    pub cell_evaluations: u64,
    /// Largest accepted whole-column time-refinement discrepancy, joules.
    /// This is an estimator, not a certified continuum or experimental bound.
    pub maximum_energy_discrepancy_j: f64,
}

#[derive(Clone, Debug, Default)]
pub(super) struct BeamAccounts {
    duration: CheckedSum,
    incident: [CheckedSum; 2],
    transmitted: [CheckedSum; 2],
    absorbed: [CheckedSum; 2],
    passive: [CheckedSum; 2],
    report: DepletingPumpEnergy,
}

impl BeamAccounts {
    fn add(
        &mut self,
        step: &PumpColumnStep<2>,
        duration_s: f64,
        energy_tolerance_j: f64,
    ) -> Result<(), OpticalHeatingError> {
        self.duration.add(duration_s)?;
        self.report.duration_s = self.duration.value()?;
        for lane in 0..2 {
            self.incident[lane].add(step.incident_energy_j[lane])?;
            self.transmitted[lane].add(step.transmitted_energy_j[lane])?;
            self.absorbed[lane].add(step.absorbed_energy_j[lane])?;
            self.passive[lane].add(step.passive_absorbed_energy_j[lane])?;
            self.report.incident_j[lane] = self.incident[lane].value()?;
            self.report.transmitted_j[lane] = self.transmitted[lane].value()?;
            self.report.absorbed_j[lane] = self.absorbed[lane].value()?;
            self.report.passive_absorbed_j[lane] = self.passive[lane].value()?;
            let mut balance = CheckedSum::default();
            for value in [
                self.report.incident_j[lane],
                -self.report.transmitted_j[lane],
                -self.report.absorbed_j[lane],
                -self.report.passive_absorbed_j[lane],
            ] {
                balance.add(value)?;
            }
            self.report.residual_j[lane] = balance.value()?;
            if self.report.residual_j[lane].abs() > energy_tolerance_j {
                return Err(invalid(
                    "depleting pump incident, transmitted and absorbed energies do not balance",
                ));
            }
        }
        self.report.accepted_kinetic_substeps = self
            .report
            .accepted_kinetic_substeps
            .checked_add(step.substeps as u64)
            .ok_or_else(|| invalid("depleting pump work counter overflow"))?;
        self.report.cell_evaluations = self
            .report
            .cell_evaluations
            .checked_add(step.cell_evaluations as u64)
            .ok_or_else(|| invalid("depleting pump work counter overflow"))?;
        self.report.maximum_energy_discrepancy_j = self
            .report
            .maximum_energy_discrepancy_j
            .max(step.max_energy_discrepancy_j);
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct BeamCell {
    emitter_index: usize,
    length_m: f64,
    volume_m3: f64,
    passive_optical_depth: [f64; 2],
    cross_sections: BeamEmitterCrossSections,
}

#[derive(Debug)]
struct BeamModel {
    cells: Vec<BeamCell>,
    wavelengths_nm: [f64; 2],
    photon_energies_j: [f64; 2],
    incident_photon_rates_s: [f64; 2],
    limits: PumpColumnLimits,
}

#[derive(Clone, Debug)]
pub(super) struct DepletingBeam {
    model: Arc<BeamModel>,
    prepared: Option<Arc<PreparedColumn>>,
    accounts: BeamAccounts,
}

/// Constitutive values only: populations and photon histories are never cached.
/// Accepted snapshots share this bounded immutable column across observer clones.
#[derive(Debug)]
pub(super) struct PreparedColumn {
    temperatures: Vec<u64>,
    column: SaturablePumpColumn<2>,
}

pub(super) struct BeamCandidate {
    pub steps: Vec<Option<SaturableSingletTripletStep>>,
    pub passive_heat_j: Vec<(CellId, f64)>,
    pub absolute_roundoff_j: f64,
    pub accounts: BeamAccounts,
}

impl<const M: usize> SaturableLuminescentHeating<M> {
    /// Couple one two-wavelength beam to the original finite emitter inventory.
    /// Bind before advancing; existing prescribed pump hazards must be zero.
    /// Each wavelength can excite both states through its own cross sections.
    /// Multiple overlapping beams require a simultaneous shared-inventory solve
    /// and are not represented by independently summing this observer.
    ///
    /// Original mass absorption is declared PASSIVE host absorption, in addition
    /// to molecular absorption. Its coefficients must exclude these same active
    /// emitters to avoid double counting. Every traversed emitter needs explicit
    /// cross sections. A passive absorbing cell needs an original thermal
    /// receiver, which can be bound as a zero-population emitter. Transparent
    /// cells with no inventory need no binding. Emission-band propagation keeps
    /// the existing explicit passive absorption/scattering environment.
    #[allow(clippy::too_many_arguments)]
    pub fn with_depleting_beam(
        mut self,
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        source_identity: ContentHash,
        beam: IncidentBeam<2>,
        cross_sections: &[BeamEmitterCrossSections],
        thermal: &DepositionSimulation<'_>,
        max_path_cells: usize,
        limits: PumpColumnLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.check_state(thermal, true, &mut checkpoint)?;
        if self.energy.duration_s != 0.0
            || self.pending.is_some()
            || self.beam.is_some()
            || source_identity == ContentHash([0; 32])
            || medium.content_identity() != self.medium_identity
            || !beam
                .wavelengths_nm
                .into_iter()
                .all(|v| v.is_finite() && v > 0.0)
            || !beam.incident_power_w.into_iter().all(nonnegative)
            || cross_sections.len() > max_path_cells
            || self
                .emitters
                .iter()
                .any(|e| e.binding.pump.rates_s != [0.0; 2])
        {
            return Err(invalid(
                "invalid depleting beam, source snapshot or nonzero prescribed pump",
            ));
        }
        let path = medium
            .population_beam_cell_path(
                cx,
                beam.origin,
                beam.unit_direction,
                beam.max_distance_m,
                max_path_cells,
            )
            .map_err(|error| match error {
                VolumeTransportError::Cancelled => {
                    OpticalHeatingError::Thermal(DepositionError::Cancelled)
                }
                _ => invalid("depleting pump geometry or traversal budget"),
            })?;
        let mut bindings = BTreeMap::new();
        for binding in cross_sections {
            checkpoint()?;
            if bindings
                .insert(binding.grid_cell, binding.clone())
                .is_some()
            {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
            }
        }
        let mut photon_energies_j = [0.0; 2];
        let mut incident_photon_rates_s = [0.0; 2];
        let mut hash = DomainHasher::new("org.frankensim.render.depleting-luminescence-pump.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(source_identity.as_bytes());
        for lane in 0..2 {
            photon_energies_j[lane] = finite(super::HC_J_NM / beam.wavelengths_nm[lane])?;
            if photon_energies_j[lane] <= 0.0 {
                return Err(invalid("unrepresentable pump photon energy"));
            }
            incident_photon_rates_s[lane] =
                finite(beam.incident_power_w[lane] / photon_energies_j[lane])?;
            if beam.incident_power_w[lane] > 0.0 && incident_photon_rates_s[lane] == 0.0 {
                return Err(invalid("positive pump photon rate is below representation"));
            }
            hash.update(&beam.wavelengths_nm[lane].to_bits().to_le_bytes());
            hash.update(&beam.incident_power_w[lane].to_bits().to_le_bytes());
        }
        let mut cells = Vec::with_capacity(path.len());
        for segment in path {
            checkpoint()?;
            let coefficients = medium
                .coefficients_at_cell(segment.cell_index, &beam.wavelengths_nm)
                .map_err(|_| invalid("passive host does not support pump wavelength"))?;
            let density = medium
                .cell_density_kg_m3(segment.cell_index)
                .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
            let mut passive_optical_depth = [0.0; 2];
            for lane in 0..2 {
                if coefficients[lane] > 0.0 && density > 0.0 {
                    passive_optical_depth[lane] = finite(positive_product(
                        coefficients[lane],
                        density,
                        segment.length_m,
                    ))?;
                    if passive_optical_depth[lane] == 0.0 {
                        return Err(invalid(
                            "positive passive pump optical depth is below representation",
                        ));
                    }
                }
            }
            let index = self
                .emitters
                .binary_search_by_key(&segment.cell_index, |e| e.binding.grid_cell);
            let emitter_index = match index {
                Ok(index) => index,
                Err(_) if passive_optical_depth == [0.0; 2] => continue,
                Err(_) => return Err(OpticalHeatingError::UnmappedCell(segment.cell_index)),
            };
            let source = bindings
                .remove(&segment.cell_index)
                .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
            hash.update(&(segment.cell_index as u64).to_le_bytes());
            hash.update(&segment.length_m.to_bits().to_le_bytes());
            hash.update(&segment.cell_volume_m3.to_bits().to_le_bytes());
            hash.update(source.identity().as_bytes());
            for depth in passive_optical_depth {
                hash.update(&depth.to_bits().to_le_bytes());
            }
            match &source.source {
                CrossSections::Authored { values, .. } => {
                    for value in values.iter().flatten() {
                        hash.update(&value.to_bits().to_le_bytes());
                    }
                }
                #[cfg(feature = "material-state")]
                CrossSections::Material { .. } => {}
            }
            cells.push(BeamCell {
                emitter_index,
                length_m: segment.length_m,
                volume_m3: segment.cell_volume_m3,
                passive_optical_depth,
                cross_sections: source,
            });
        }
        if !bindings.is_empty() {
            return Err(invalid(
                "pump cross-section binding is outside the traversed inventory",
            ));
        }
        let mut pump = DepletingBeam {
            model: Arc::new(BeamModel {
                cells,
                wavelengths_nm: beam.wavelengths_nm,
                photon_energies_j,
                incident_photon_rates_s,
                limits,
            }),
            prepared: None,
            accounts: BeamAccounts::default(),
        };
        pump.prepared = Some(pump.prepare(&self, thermal, &mut checkpoint)?);
        checkpoint()?;
        self.source_identity = hash.finalize();
        self.beam = Some(pump);
        Ok(self)
    }

    /// Accepted incident/transmitted/molecular/passive pump energies by actual
    /// wavelength. Absent for the original externally maintained local field.
    pub fn depleting_pump_energy(&self) -> Option<&DepletingPumpEnergy> {
        self.beam.as_ref().map(|beam| &beam.accounts.report)
    }
}

impl DepletingBeam {
    fn prepare<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedColumn>, OpticalHeatingError> {
        // Native lifecycle hooks can inspect an unchanged accepted state several
        // times. Reuse its exact constitutive snapshot, while checking every
        // original cell and cancellation without allocating or re-querying data.
        if let Some(prepared) = &self.prepared {
            let mut matches = true;
            for (cell, expected) in self.model.cells.iter().zip(&prepared.temperatures) {
                checkpoint()?;
                let thermal_cell = owner.emitters[cell.emitter_index].binding.thermal_cell;
                if thermal.temperature_k(thermal_cell)?.to_bits() != *expected {
                    matches = false;
                    break;
                }
            }
            if matches {
                return Ok(Arc::clone(prepared));
            }
        }
        let mut cells = Vec::with_capacity(self.model.cells.len());
        let mut temperatures = Vec::with_capacity(self.model.cells.len());
        let mut cache = None;
        for cell in &self.model.cells {
            checkpoint()?;
            let binding = owner.emitters[cell.emitter_index].binding;
            let temperature = thermal.temperature_k(binding.thermal_cell)?;
            temperatures.push(temperature.to_bits());
            let model = resolve_model(&owner.singlet, &owner.triplet, temperature, &mut cache)?;
            let sections = cell
                .cross_sections
                .evaluate(self.model.wavelengths_nm, temperature)?;
            let mut depth = [[0.0; 2]; 2];
            for lane in 0..2 {
                for target in 0..2 {
                    let source = if target == 0 {
                        &owner.singlet
                    } else {
                        owner.triplet.triplet()
                    };
                    let (lo, hi) = source.fluorescence().excitation_bounds_nm();
                    if sections[lane][target] > 0.0
                        && !(lo..=hi).contains(&self.model.wavelengths_nm[lane])
                    {
                        return Err(invalid(
                            "molecular pump cross section exceeds original excitation support",
                        ));
                    }
                    depth[lane][target] = super::super::rate_power(
                        sections[lane][target],
                        cell.length_m,
                        cell.volume_m3,
                    )?;
                    if sections[lane][target] > 0.0 && depth[lane][target] == 0.0 {
                        return Err(invalid(
                            "positive molecular pump optical depth is below representation",
                        ));
                    }
                }
            }
            cells.push(SaturablePumpCell {
                model,
                optical_depth_per_ground: depth,
                passive_optical_depth: cell.passive_optical_depth,
            });
        }
        let column = preserve_checkpoint(checkpoint, |cancelled| {
            SaturablePumpColumn::try_new(
                &cells,
                self.model.photon_energies_j,
                self.model.cells.len().max(1),
                cancelled,
            )
        })?;
        Ok(Arc::new(PreparedColumn {
            temperatures,
            column,
        }))
    }

    pub(super) fn check_state<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedColumn>, OpticalHeatingError> {
        self.prepare(owner, thermal, checkpoint)
    }

    pub(super) fn advance<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        duration_s: f64,
        pumped: bool,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BeamCandidate, OpticalHeatingError> {
        let prepared = self.prepare(owner, thermal, checkpoint)?;
        let populations: Vec<_> = self
            .model
            .cells
            .iter()
            .map(|cell| owner.emitters[cell.emitter_index].populations)
            .collect();
        let rates = if pumped {
            self.model.incident_photon_rates_s
        } else {
            [0.0; 2]
        };
        let step = preserve_checkpoint(checkpoint, |cancelled| {
            prepared.column.advance(
                &populations,
                rates,
                duration_s,
                self.model.limits,
                cancelled,
            )
        })?;
        let mut accounts = self.accounts.clone();
        accounts.add(&step, duration_s, owner.energy_tolerance_j)?;
        let mut roundoff = CheckedSum::default();
        roundoff.add(step.residual_j.abs())?;
        for lane in 0..2 {
            roundoff.add(finite(
                step.photon_residual[lane].abs() * self.model.photon_energies_j[lane],
            )?)?;
        }
        let mut steps = vec![None; owner.emitters.len()];
        let mut passive_heat_j = Vec::with_capacity(self.model.cells.len());
        for (cell, result) in self.model.cells.iter().zip(step.cells) {
            checkpoint()?;
            let binding = owner.emitters[cell.emitter_index].binding;
            passive_heat_j.push((binding.thermal_cell, result.passive_heat_j));
            steps[cell.emitter_index] = Some(result.kinetics);
        }
        checkpoint()?;
        Ok(BeamCandidate {
            steps,
            passive_heat_j,
            absolute_roundoff_j: roundoff.value()?,
            accounts,
        })
    }

    pub(super) fn accept(&mut self, accounts: BeamAccounts, prepared: Arc<PreparedColumn>) {
        self.accounts = accounts;
        self.prepared = Some(prepared);
    }
}
