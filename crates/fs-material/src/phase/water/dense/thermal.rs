//! Original dense-water internal energy on the shared thermal timeline.
//!
//! Each sealed volume stores m*u, has fixed mass and volume, and obtains its
//! temperature and pressure from the original IF97 isochore. The timeline's
//! legacy `specific_h_j_kg` field therefore carries u for this explicitly
//! InternalEnergy law. It must not enter open-stream enthalpy advection.
//!
//! Contacts and bath conductances are authored network inputs, not inferred
//! dense-water conductivity. Each node is a spatially uniform equilibrium
//! volume; wall storage, compliance, internal convection and resolved spatial
//! gradients require their own owners. The original timeline owns every heat
//! flux, event, timestep, entropy account and all-or-nothing state publication.

use super::{DenseIsochore, DenseWaterError, DenseWaterState};
use fs_thermal::enthalpy::entropy::EntropyAudit;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis, TimedEvent,
};

/// The original prepared isochore, expressed as internal-energy storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsochoricThermalLaw {
    chart: DenseIsochore,
}
impl DenseIsochoricThermalLaw {
    pub const fn new(chart: DenseIsochore) -> Self {
        Self { chart }
    }
    pub const fn chart(&self) -> &DenseIsochore {
        &self.chart
    }
    pub fn state_at_specific_internal_energy(
        &self,
        specific_u_j_kg: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DenseWaterState, DenseWaterError> {
        self.chart
            .state_at_specific_internal_energy(specific_u_j_kg, cancelled)
            .map(|solution| solution.into_state())
    }
    /// Local analytic dT/du [kg K/J] from the SAME original source Cv.
    pub fn temperature_internal_energy_derivative_kg_k_j(
        &self,
        specific_u_j_kg: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, DenseWaterError> {
        let state = self.state_at_specific_internal_energy(specific_u_j_kg, cancelled)?;
        Ok(1.0 / state.source().cv().value())
    }
}
impl ThermalLaw for DenseIsochoricThermalLaw {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        let expected_mass = self.chart.density().value() * cell.geometry.volume_m3;
        // This is an arithmetic admission, not a density uncertainty allowance.
        if !cell.geometry.volume_m3.is_finite()
            || cell.geometry.volume_m3 <= 0.0
            || !expected_mass.is_finite()
            || expected_mass <= 0.0
            || !cell.mass_kg.is_finite()
            || cell.mass_kg <= 0.0
            || (cell.mass_kg - expected_mass).abs()
                > 64.0 * f64::EPSILON * cell.mass_kg.max(expected_mass)
        {
            return Err(DepositionError::InvalidInput(
                "rigid-water cell mass/volume must match its original isochore density",
            ));
        }
        let [lo, hi] = self.chart.internal_energy_bounds_j_kg();
        if !cell.specific_h_j_kg.is_finite() || !(lo..=hi).contains(&cell.specific_h_j_kg) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
    fn temperature_k(&self, specific_u_j_kg: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_internal_energy(specific_u_j_kg, || false)
            .map(|state| state.temperature_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.chart.maximum_temperature_slope()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_u_j_kg: f64,
        final_u_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let a = self
            .state_at_specific_internal_energy(initial_u_j_kg, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if initial_u_j_kg == final_u_j_kg {
            return Ok(Some(0.0));
        }
        let b = self
            .state_at_specific_internal_energy(final_u_j_kg, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let delta_u = final_u_j_kg - initial_u_j_kg;
        // Use the source integral's mean reciprocal temperature with the
        // ACTUAL stored increment. Distinct u coordinates can round to the
        // same inverse T. Never erase that heat or its du/T contribution.
        // Source endpoint entropy is separately reported by the network.
        let reciprocal_temperature = if a.temperature_k() == b.temperature_k() {
            1.0 / a.temperature_k()
        } else {
            let difference = self
                .chart
                .source_interval()
                .difference(a.source().temperature(), b.source().temperature())
                .map_err(|_| DepositionError::MaterialDomain)?;
            difference.specific_entropy_change().value()
                / difference.specific_internal_energy_change().value()
        };
        let entropy = delta_u * reciprocal_temperature;
        if !reciprocal_temperature.is_finite()
            || reciprocal_temperature <= 0.0
            || !entropy.is_finite()
            || entropy == 0.0
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Some(entropy))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DenseWaterThermalError {
    Thermal(DepositionError),
    Material(DenseWaterError),
    Invalid(&'static str),
}
impl From<DepositionError> for DenseWaterThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<DenseWaterError> for DenseWaterThermalError {
    fn from(error: DenseWaterError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for DenseWaterThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid dense-water thermal network: {self:?}")
    }
}
impl std::error::Error for DenseWaterThermalError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterThermalCell {
    pub cell_id: CellId,
    pub mass_kg: f64,
    pub volume_m3: f64,
    pub stored_specific_internal_energy_j_kg: f64,
    pub material: DenseWaterState,
}

/// Source differences and the existing timeline's physical heat/entropy account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterThermalAudit {
    pub thermal: EnergyAudit,
    pub entropy: EntropyAudit,
    pub internal_energy_change_j: f64,
    pub enthalpy_change_j: f64,
    /// Sum V*Delta p. This is part of Delta H, not moving-boundary work.
    pub pressure_enthalpy_change_j: f64,
    /// Exactly zero for these fixed-volume, sealed control volumes.
    pub boundary_work_j: f64,
    /// Original-source entropy increments; no coordinate correction applied.
    pub source_entropy_change_j_k: f64,
    /// Thermal-coordinate entropy minus original-source endpoint entropy.
    /// Exposes inversion/representation effects; not physical entropy production.
    pub entropy_storage_source_defect_j_k: f64,
    pub internal_energy_residual_j: f64,
    pub enthalpy_residual_j: f64,
    pub represented_internal_energy_residual_j: f64,
    pub represented_enthalpy_residual_j: f64,
}

/// Sealed rigid volumes on one original fixed-density caloric chart.
///
/// Use MaterialId(0) and put every birth at time zero; each birth's stored
/// coordinate is original specific INTERNAL energy. Existing contact and bath
/// events retain their native semantics. Later mass insertion is not this model.
#[derive(Clone)]
pub struct DenseWaterThermalNetwork<'a> {
    law: &'a DenseIsochoricThermalLaw,
    thermal: DepositionSimulation<'a>,
}
impl<'a> DenseWaterThermalNetwork<'a> {
    pub fn new(
        law: &'a DenseIsochoricThermalLaw,
        events: Vec<TimedEvent>,
        budget: DepositionBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DenseWaterThermalError> {
        checkpoint()?;
        let mut births = 0;
        for event in &events {
            checkpoint()?;
            if let Event::Birth(cell) = &event.event {
                if event.time_s != 0.0 || cell.material_id != MaterialId(0) {
                    return Err(DenseWaterThermalError::Invalid(
                        "sealed water volumes require time-zero births under MaterialId(0)",
                    ));
                }
                births += 1;
            }
        }
        if births == 0 {
            return Err(DenseWaterThermalError::Invalid(
                "at least one sealed water volume required",
            ));
        }
        let thermal =
            DepositionSimulation::new(events, &[(MaterialId(0), law as &dyn ThermalLaw)], budget)?
                .advance_to(
                    0.0,
                    StepConfig {
                        max_dt_s: 1.0,
                        energy_tolerance_j: 0.0,
                    },
                    &mut checkpoint,
                )?
                .0;
        checkpoint()?;
        Ok(Self { law, thermal })
    }
    pub const fn law(&self) -> &'a DenseIsochoricThermalLaw {
        self.law
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'a> {
        &self.thermal
    }
    pub fn time_s(&self) -> f64 {
        self.thermal.time_s()
    }
    pub fn profile(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<DenseWaterThermalCell>, DenseWaterThermalError> {
        let mut profile = Vec::with_capacity(self.thermal.cells().len());
        for cell in self.thermal.cells().values() {
            checkpoint()?;
            self.law.validate_cell(cell)?;
            // The source inverse is explicitly bounded; cancellation is checked
            // before and after it, with no external physical history to mutate.
            let material = self
                .law
                .state_at_specific_internal_energy(cell.specific_h_j_kg, || false)?;
            profile.push(DenseWaterThermalCell {
                cell_id: cell.id,
                mass_kg: cell.mass_kg,
                volume_m3: cell.geometry.volume_m3,
                stored_specific_internal_energy_j_kg: cell.specific_h_j_kg,
                material,
            });
        }
        checkpoint()?;
        Ok(profile)
    }
    pub fn advance_to(
        &self,
        target_s: f64,
        controls: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, DenseWaterThermalAudit), DenseWaterThermalError> {
        let (thermal, energy, entropy) =
            self.thermal
                .advance_to_with_entropy(target_s, controls, &mut checkpoint)?;
        let next = Self {
            law: self.law,
            thermal,
        };
        let before = self.profile(&mut checkpoint)?;
        let after = next.profile(&mut checkpoint)?;
        // Births were all consumed at construction. Every subsequent increment
        // is heat on the same closed masses and rigid geometries.
        let mut sums = [Sum::default(); 7];
        for (a, b) in before.iter().zip(&after) {
            checkpoint()?;
            let difference = self
                .law
                .chart
                .source_interval()
                .difference(
                    a.material.source().temperature(),
                    b.material.source().temperature(),
                )
                .map_err(DenseWaterError::from)?;
            for (sum, value) in sums.iter_mut().zip([
                a.mass_kg * difference.specific_internal_energy_change().value(),
                a.mass_kg * difference.specific_enthalpy_change().value(),
                a.volume_m3 * difference.pressure_change().value(),
                a.mass_kg * difference.specific_entropy_change().value(),
                a.mass_kg
                    * (b.material.specific_internal_energy_j_kg()
                        - a.material.specific_internal_energy_j_kg()),
                a.mass_kg
                    * (b.material.specific_enthalpy_j_kg() - a.material.specific_enthalpy_j_kg()),
                a.volume_m3 * (b.material.pressure_pa() - a.material.pressure_pa()),
            ]) {
                sum.add(value)?;
            }
        }
        let [
            du,
            dh,
            vdp,
            ds,
            represented_u,
            represented_h,
            represented_vdp,
        ] = sums.map(Sum::total);
        let residuals = [
            du - energy.external_heat_j,
            dh - vdp - energy.external_heat_j,
            represented_u - energy.external_heat_j,
            represented_h - represented_vdp - energy.external_heat_j,
            energy.residual_j,
        ];
        if residuals.iter().any(|value| !value.is_finite()) || !ds.is_finite() {
            return Err(DenseWaterError::Unrepresentable("thermal-network source audit").into());
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > controls.energy_tolerance_j {
            return Err(DenseWaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: controls.energy_tolerance_j,
            }
            .into());
        }
        let entropy_defect = entropy.thermal_storage_change_j_k - ds;
        if !entropy_defect.is_finite() {
            return Err(DenseWaterError::Unrepresentable("thermal-network entropy audit").into());
        }
        checkpoint()?;
        Ok((
            next,
            DenseWaterThermalAudit {
                thermal: energy,
                entropy,
                internal_energy_change_j: du,
                enthalpy_change_j: dh,
                pressure_enthalpy_change_j: vdp,
                boundary_work_j: 0.0,
                source_entropy_change_j_k: ds,
                entropy_storage_source_defect_j_k: entropy_defect,
                internal_energy_residual_j: residuals[0],
                enthalpy_residual_j: residuals[1],
                represented_internal_energy_residual_j: residuals[2],
                represented_enthalpy_residual_j: residuals[3],
            },
        ))
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), DepositionError> {
        let next = self.value + term;
        if !term.is_finite() || !next.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        self.correction += if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.value = next;
        if !self.correction.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok(())
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}

#[cfg(test)]
mod tests;
