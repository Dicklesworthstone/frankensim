//! Shared temperature-dependent optical constitutive state for native heating.
//!
//! One immutable field owns source-bound absorption/scattering tables, fixed
//! geometry and phase laws, and mappings onto original enthalpy cells. Incident
//! illumination and thermal emission can share it in Arc without inventing a
//! beam, a second interpolation law, or an independently integrated temperature.
//! Spectral samples are compiled once per table knot. Each state preparation
//! reads every distinct native temperature once and evaluates each distinct
//! receiver/table-pair once; unchanged coefficient bits retain the old packets.

use super::OpticalHeatingError;
use super::feedback::{
    CompiledTable, OpticalFeedbackError, TemperatureAbsorptionTable, TemperatureSpectrum,
};
use crate::volumes::transport::scattering::{
    CellPhaseFunction, CellScatteringGrid, MassScatteringSpectrum, PhaseFunction,
};
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// Source-bound mass scattering in m²/kg versus kelvin. The retained rows keep
/// their original spectral owners, including source support and analytic laws.
/// Temperature interpolation reuses the original absorption table compiler.
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

/// Local material laws, native thermal receiver and fixed phase-palette entry.
/// Multiple original optical cells may share one thermal receiver.
#[derive(Clone, Copy, Debug)]
pub struct ThermalScatteringCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub absorption_table: usize,
    pub scattering_table: usize,
    pub phase: usize,
}

/// Immutable source/geometry storage ceilings, independent of photon or thermal
/// trajectory work. Consumers retain their original per-advance work budgets.
#[derive(Clone, Copy, Debug)]
pub struct TemperatureOpticalFieldBudget {
    pub max_grid_cells: usize,
    /// Combined absorption and scattering table count.
    pub max_tables: usize,
    pub max_total_temperature_knots: usize,
}

/// Every positive-density material cell, including transparent material and
/// pure scatterers, with the intersection of its two coefficient-table domains.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureOpticalBinding {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
}

#[derive(Clone, Copy)]
struct BoundMaterial {
    thermal_cell: CellId,
    temperature_slot: usize,
    absorption_table: usize,
    scattering_table: usize,
    representative_cell: usize,
}
#[derive(Clone, Copy)]
struct BoundTemperature {
    thermal_cell: CellId,
    minimum_k: f64,
    maximum_k: f64,
}

/// Immutable source-bound constitutive data, independent of incident sources,
/// bath temperatures, photon histories, or accumulated heat. Share in Arc when
/// multiple observers read the same optical specimen. The density grid remains
/// borrowed; the source tables are compiled into wavelength packets once.
pub struct TemperatureOpticalField<'g, const N: usize> {
    identity: ContentHash,
    geometry: AbsorbingGrid<'g>,
    wavelengths_nm: [f64; N],
    bounds_nm: (f64, f64),
    absorption: Vec<CompiledTable<N>>,
    scattering: Vec<CompiledTable<N>>,
    materials: Vec<BoundMaterial>,
    temperatures: Vec<BoundTemperature>,
    cell_materials: Vec<Option<usize>>,
    phases: Vec<PhaseFunction>,
    cell_phases: Vec<usize>,
    material_bindings: Vec<TemperatureOpticalBinding>,
    temperature_knots: usize,
    has_scattering: bool,
}

impl<const N: usize> core::fmt::Debug for TemperatureOpticalField<'_, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("TemperatureOpticalField")
            .field("identity", &self.identity)
            .field("wavelengths_nm", &self.wavelengths_nm)
            .field("grid_cells", &self.geometry.cell_count())
            .field("material_cells", &self.material_bindings.len())
            .field("unique_native_temperatures", &self.temperatures.len())
            .field("unique_table_pairs", &self.materials.len())
            .finish_non_exhaustive()
    }
}

/// Candidate-local coefficient storage. Public access is immutable; native
/// observers refresh their private candidate before requesting transport. A
/// state from a different source field cannot be used as that field's packet.
#[derive(Clone, Debug)]
pub struct TemperatureOpticalState<const N: usize> {
    identity: ContentHash,
    absorption: Vec<[f64; N]>,
    scattering: Vec<[f64; N]>,
    materials: Vec<([f64; N], [f64; N])>,
    temperatures: Vec<f64>,
    initialized: bool,
}
impl<const N: usize> TemperatureOpticalState<N> {
    pub fn absorption_m2_kg(&self) -> &[[f64; N]] {
        &self.absorption
    }
    pub fn scattering_m2_kg(&self) -> &[[f64; N]] {
        &self.scattering
    }
}

impl<'g, const N: usize> TemperatureOpticalField<'g, N> {
    /// Compile and admit one shared material field at an active native state.
    /// Every positive-density geometry cell must be mapped, even if its current
    /// absorption/scattering is zero. No wavelength or temperature extrapolates.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        geometry: AbsorbingGrid<'g>,
        wavelengths_nm: [f64; N],
        absorption: &[TemperatureAbsorptionTable],
        scattering: &[TemperatureScatteringTable],
        phases: &[PhaseFunction],
        bindings: &[ThermalScatteringCell],
        thermal: &DepositionSimulation<'_>,
        budget: TemperatureOpticalFieldBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-optical-field.v1");
        hash.update(geometry.content_identity().as_bytes());
        hash.update(&(N as u64).to_le_bytes());
        for wavelength in wavelengths_nm {
            hash.update(&wavelength.to_bits().to_le_bytes());
        }
        Self::compile_with_hash(
            geometry,
            wavelengths_nm,
            absorption,
            scattering,
            phases,
            bindings,
            thermal,
            budget,
            hash,
            &mut checkpoint,
        )
        .map(|(field, _)| field)
    }

    // Existing incident heating seeds this with its original source prefix.
    // The remaining source/table/binding/phase hash order is unchanged, keeping
    // established medium identities and incident trajectories byte-for-byte.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn compile_with_hash(
        geometry: AbsorbingGrid<'g>,
        wavelengths_nm: [f64; N],
        absorption: &[TemperatureAbsorptionTable],
        scattering: &[TemperatureScatteringTable],
        phases: &[PhaseFunction],
        bindings: &[ThermalScatteringCell],
        thermal: &DepositionSimulation<'_>,
        budget: TemperatureOpticalFieldBudget,
        mut hash: DomainHasher,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, TemperatureOpticalState<N>), OpticalFeedbackError> {
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
        if absorption.is_empty()
            || scattering.is_empty()
            || phases.is_empty()
            || thermal.cells().is_empty()
        {
            return Err(invalid(
                "nonempty optical laws and active thermal cells required",
            ));
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
            compiled_absorption.push(table.compile(&wavelengths_nm, &mut hash, checkpoint)?);
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
                &wavelengths_nm,
                &mut hash,
                checkpoint,
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
        let mut temperatures: Vec<BoundTemperature> = Vec::new();
        let mut material_bindings = Vec::new();
        materials
            .try_reserve_exact(bindings.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        temperatures
            .try_reserve_exact(bindings.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        material_bindings
            .try_reserve_exact(bindings.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        let mut unique_materials = BTreeMap::new();
        let mut unique_temperatures = BTreeMap::new();
        let mut cell_materials = allocate(count, None)?;
        let mut cell_phases = allocate(count, 0usize)?;
        for (index, binding) in selected.into_iter().enumerate() {
            checkpoint()?;
            let positive_density = geometry
                .cell_density_kg_m3(index)
                .is_some_and(|rho| rho > 0.0);
            let Some(binding) = binding else {
                if positive_density {
                    return Err(OpticalHeatingError::UnmappedCell(index).into());
                }
                continue;
            };
            let a_domain = absorption[binding.absorption_table].temperature_bounds_k();
            let s_domain = scattering[binding.scattering_table].temperature_bounds_k();
            let minimum_k = a_domain.0.max(s_domain.0);
            let maximum_k = a_domain.1.min(s_domain.1);
            let temperature_slot = *unique_temperatures
                .entry(binding.thermal_cell)
                .or_insert_with(|| {
                    let slot = temperatures.len();
                    temperatures.push(BoundTemperature {
                        thermal_cell: binding.thermal_cell,
                        minimum_k,
                        maximum_k,
                    });
                    slot
                });
            temperatures[temperature_slot].minimum_k =
                temperatures[temperature_slot].minimum_k.max(minimum_k);
            temperatures[temperature_slot].maximum_k =
                temperatures[temperature_slot].maximum_k.min(maximum_k);
            let key = (
                binding.thermal_cell,
                binding.absorption_table,
                binding.scattering_table,
            );
            let slot = *unique_materials.entry(key).or_insert_with(|| {
                let slot = materials.len();
                materials.push(BoundMaterial {
                    thermal_cell: binding.thermal_cell,
                    temperature_slot,
                    absorption_table: binding.absorption_table,
                    scattering_table: binding.scattering_table,
                    representative_cell: index,
                });
                slot
            });
            cell_materials[index] = Some(slot);
            cell_phases[index] = binding.phase;
            if positive_density {
                material_bindings.push(TemperatureOpticalBinding {
                    grid_cell: index,
                    thermal_cell: binding.thermal_cell,
                    minimum_temperature_k: minimum_k,
                    maximum_temperature_k: maximum_k,
                });
            }
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
        let model = Self {
            identity: hash.finalize(),
            geometry,
            wavelengths_nm,
            bounds_nm: bounds,
            absorption: compiled_absorption,
            scattering: compiled_scattering,
            materials,
            temperatures,
            cell_materials,
            phases: retained_phases,
            cell_phases,
            material_bindings,
            temperature_knots: knots,
            has_scattering,
        };
        let state = model.initial_state(thermal, &mut *checkpoint)?;
        transport_checkpoint(checkpoint, |check| model.field(&state, check).map(|_| ()))?;
        Ok((model, state))
    }

    pub fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn geometry(&self) -> &AbsorbingGrid<'g> {
        &self.geometry
    }
    pub fn wavelengths_nm(&self) -> &[f64; N] {
        &self.wavelengths_nm
    }
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.bounds_nm
    }
    pub fn material_bindings(&self) -> &[TemperatureOpticalBinding] {
        &self.material_bindings
    }
    /// Conservative over all admitted temperatures, including initially zero laws.
    pub fn has_scattering(&self) -> bool {
        self.has_scattering
    }
    pub fn unique_temperature_count(&self) -> usize {
        self.temperatures.len()
    }
    pub fn unique_material_count(&self) -> usize {
        self.materials.len()
    }
    pub fn table_count(&self) -> usize {
        self.absorption.len() + self.scattering.len()
    }
    pub fn temperature_knot_count(&self) -> usize {
        self.temperature_knots
    }

    /// Allocate one candidate's packet/scratch buffers and read the actual
    /// native temperature state. The field's source identity remains immutable.
    pub fn initial_state(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TemperatureOpticalState<N>, OpticalFeedbackError> {
        checkpoint()?;
        let mut state = TemperatureOpticalState {
            identity: self.identity,
            absorption: allocate(self.geometry.cell_count(), [0.0; N])?,
            scattering: allocate(self.geometry.cell_count(), [0.0; N])?,
            materials: allocate(self.materials.len(), ([0.0; N], [0.0; N]))?,
            temperatures: allocate(self.temperatures.len(), 0.0)?,
            initialized: false,
        };
        self.prepare(thermal, &mut state, &mut checkpoint)?;
        Ok(state)
    }

    /// Refresh private native-candidate buffers without allocating. A refusal
    /// discards the enclosing native candidate; scratch state is not published.
    /// Returns true only if optical coefficient bits changed (or at first use).
    /// A temperature change along a constant table therefore preserves transport.
    pub(crate) fn prepare(
        &self,
        thermal: &DepositionSimulation<'_>,
        state: &mut TemperatureOpticalState<N>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<bool, OpticalFeedbackError> {
        checkpoint()?;
        if state.identity != self.identity {
            return Err(VolumeTransportError::SourceIdentity.into());
        }
        for (binding, value) in self.temperatures.iter().zip(&mut state.temperatures) {
            checkpoint()?;
            *value = thermal.temperature_k(binding.thermal_cell)?;
            check_temperature(*binding, *value)?;
        }
        let mut changed = !state.initialized;
        for (material, values) in self.materials.iter().zip(&mut state.materials) {
            checkpoint()?;
            let temperature = state.temperatures[material.temperature_slot];
            values.0 = self.absorption[material.absorption_table]
                .at(material.thermal_cell, temperature)?;
            values.1 = self.scattering[material.scattering_table]
                .at(material.thermal_cell, temperature)?;
            for lane in 0..N {
                changed |= values.0[lane].to_bits()
                    != state.absorption[material.representative_cell][lane].to_bits()
                    || values.1[lane].to_bits()
                        != state.scattering[material.representative_cell][lane].to_bits();
            }
        }
        if changed {
            for (index, slot) in self.cell_materials.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint()?;
                }
                if let Some(slot) = slot {
                    (state.absorption[index], state.scattering[index]) = state.materials[*slot];
                }
            }
        }
        checkpoint()?;
        state.initialized = true;
        Ok(changed)
    }

    /// Validate all source domains at an accepted endpoint/event, including
    /// cells that presently emit or absorb no radiation.
    pub fn validate_temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalFeedbackError> {
        for binding in &self.temperatures {
            checkpoint()?;
            check_temperature(*binding, thermal.temperature_k(binding.thermal_cell)?)?;
        }
        checkpoint()?;
        Ok(())
    }

    /// Borrow the original exact-cell transport owner over this prepared packet.
    /// Its native coefficient/phase/Kirchhoff admission remains authoritative.
    pub fn field<'s>(
        &'s self,
        values: &'s TemperatureOpticalState<N>,
        check: &mut dyn FnMut() -> Result<(), VolumeTransportError>,
    ) -> Result<CellScatteringGrid<'s, 'g, N>, VolumeTransportError> {
        check()?;
        if !values.initialized || values.identity != self.identity {
            return Err(VolumeTransportError::SourceIdentity);
        }
        CellScatteringGrid::try_new_with_checkpoint(
            &mut || check(),
            &self.geometry,
            self.wavelengths_nm,
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

    pub(crate) fn receiver(&self, grid_cell: usize) -> Option<CellId> {
        self.cell_materials
            .get(grid_cell)
            .copied()
            .flatten()
            .map(|slot| self.materials[slot].thermal_cell)
    }
}

fn check_temperature(
    binding: BoundTemperature,
    temperature_k: f64,
) -> Result<(), OpticalFeedbackError> {
    if !temperature_k.is_finite()
        || !(binding.minimum_k..=binding.maximum_k).contains(&temperature_k)
    {
        return Err(OpticalFeedbackError::TemperatureDomain {
            thermal_cell: binding.thermal_cell,
            temperature_k,
            minimum_k: binding.minimum_k,
            maximum_k: binding.maximum_k,
        });
    }
    Ok(())
}
fn allocate<T: Clone>(count: usize, value: T) -> Result<Vec<T>, OpticalFeedbackError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| OpticalFeedbackError::Budget)?;
    result.resize(count, value);
    Ok(result)
}
pub(super) fn transport_checkpoint<T>(
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

#[cfg(test)]
mod tests;
