//! Original bulk acoustic viscous heat on the finite thermal timeline.
//!
//! A loaded cylinder supplies its coherent analytic volume integrals. Only
//! `eta |du/dx|² / 2` enters the original fixed-volume internal-energy cells.
//! Conductive entropy production and its acoustic availability loss remain a
//! separate account: neither is an additional caloric heat source. This is not
//! a reconstruction of mean compressional work, acoustic streaming, thermal
//! end layers, molecular relaxation, or wall-layer heat.
//!
//! The original harmonic field is held over each requested thermal advance.
//! Its cycle-mean watts are not an exact instantaneous short-window integral.
//! The caller supplies the original fluid caloric law and justifies uniform
//! thermal cells, fixed volume, and a sufficiently small background change.
//! A declared temperature excursion is checked after every native step; it
//! does not certify material validity or unchanged phase/composition. Refresh
//! the original material/acoustic solve explicitly when that model ceases to
//! apply. No caloric law, conduction solver, or pressure work is invented here.

use fs_duct::fluid::{FluidDuctError, bulk::BulkTractionAbsorptionResponse};
use fs_material::fluid_acoustics::thermoviscous::absorption::BulkCoherentCellAbsorption;
use fs_math::{dd::Dd, det};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalStorageBasis,
};
use std::collections::BTreeMap;

/// Original coherent spectral heat integrated on each actual thermal interval.
#[path = "bulk_viscous_heating/finite_window.rs"]
pub mod finite_window;

/// Hard bound on a complete, explicitly authored thermal partition.
pub const MAX_BULK_HEATING_CELLS: usize = 4096;

/// Existing cell in inlet-to-outlet order; its near face is the previous end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkHeatCell {
    pub cell: CellId,
    /// Far face distance from the original cylinder inlet [m].
    pub end_distance_m: f64,
}

/// Cached physical integral and retained original caloric/geometry binding.
#[derive(Clone, Debug, PartialEq)]
pub struct BulkViscousCellHeating {
    pub cell: CellId,
    pub start_distance_m: f64,
    pub end_distance_m: f64,
    /// Already cycle averaged and integrated over the actual full cell volume.
    /// Only `viscous_heat_w` is deposited in thermal storage.
    pub absorption: BulkCoherentCellAbsorption,
    geometry: CellGeometry,
    material: MaterialId,
    mass_kg: f64,
    body_id: u64,
    segment_id: u64,
    caloric_identity: [u8; 32],
}

/// Refusal leaves the original acoustic response and thermal state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum BulkViscousHeatingError {
    Acoustic(FluidDuctError),
    Thermal(DepositionError),
    Invalid(&'static str),
    TemperatureExcursion {
        cell: CellId,
        actual_k: f64,
        maximum_k: f64,
    },
    PowerBalance {
        residual_w: f64,
        tolerance_w: f64,
    },
}
impl core::fmt::Display for BulkViscousHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bulk viscous heating refused: {self:?}")
    }
}
impl std::error::Error for BulkViscousHeatingError {}
impl From<FluidDuctError> for BulkViscousHeatingError {
    fn from(error: FluidDuctError) -> Self {
        Self::Acoustic(error)
    }
}
impl From<DepositionError> for BulkViscousHeatingError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Distinct acoustic quantities integrated on the actual accepted thermal clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkAcousticEnergy {
    /// Full acoustic work through the inlet face, signed in the +x direction.
    pub left_acoustic_work_j: f64,
    /// Full acoustic work through the outlet face, signed in the +x direction.
    pub right_acoustic_work_j: f64,
    /// Direct viscous heat, already included ONCE in the native heat report.
    pub viscous_heat_j: f64,
    /// Conductive entropy production [J/K], not extra caloric energy.
    pub conductive_entropy_j_k: f64,
    /// T0 times conductive entropy production; never deposited as extra heat.
    pub thermal_availability_loss_j: f64,
    /// Left minus right minus viscous heat minus thermal availability loss.
    pub balance_residual_j: f64,
}

/// Native thermal candidate and independent acoustic diagnostic account.
pub struct BulkViscousHeatingStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// Includes viscous heat AND original sources, reservoirs and exchanges.
    pub heat_report: EnergyAudit,
    pub acoustic: BulkAcousticEnergy,
    /// Viscous joules only, in original physical partition order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

/// One prepared source on original, fixed-volume internal-energy storage.
/// Preparation integrates acoustics once per cell. Native substeps copy cached
/// watts and check original bindings; no acoustic solve or spatial quadrature
/// is performed on the thermal hot path. Storage is O(cells).
pub struct BulkViscousHeating<'a> {
    response: &'a BulkTractionAbsorptionResponse,
    cells: Vec<BulkViscousCellHeating>,
    sources_w: BTreeMap<CellId, f64>,
    whole: BulkCoherentCellAbsorption,
    maximum_temperature_change_k: f64,
}

fn finite(value: f64) -> Result<f64, BulkViscousHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn norm(value: [f64; 3]) -> f64 {
    det::hypot(det::hypot(value[0], value[1]), value[2])
}
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 256.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn power_balance(residual: f64, scale: f64, cells: usize) -> Result<(), BulkViscousHeatingError> {
    let residual_w = finite(residual)?;
    let tolerance_w = finite(512.0 * f64::EPSILON * (cells + 1) as f64 * scale)?;
    if residual_w.abs() > tolerance_w {
        Err(BulkViscousHeatingError::PowerBalance {
            residual_w,
            tolerance_w,
        })
    } else {
        Ok(())
    }
}

impl<'a> BulkViscousHeating<'a> {
    /// Bind a complete partition of the ORIGINAL physical cylinder.
    ///
    /// Cell centreline endpoints are the centres of its opposite full-area
    /// faces. Cells must adjoin on one straight axis and physical body, with
    /// original source density, temperature, and an identified internal-energy
    /// law. Geometry, mass and caloric identity are retained and rechecked.
    /// Numeric rho/T agreement does not establish source provenance: the caller
    /// must supply the original caloric law, including its actual Cv and datum.
    /// An enthalpy/Cp law cannot stand in for fixed-volume gas storage.
    pub fn try_new(
        response: &'a BulkTractionAbsorptionResponse,
        thermal: &DepositionSimulation<'_>,
        partition: &[BulkHeatCell],
        maximum_cells: usize,
        maximum_temperature_change_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, BulkViscousHeatingError> {
        checkpoint()?;
        if maximum_cells == 0
            || maximum_cells > MAX_BULK_HEATING_CELLS
            || partition.is_empty()
            || partition.len() > maximum_cells
        {
            return Err(DepositionError::Budget("bulk viscous heating cells").into());
        }
        if !maximum_temperature_change_k.is_finite() || maximum_temperature_change_k <= 0.0 {
            return Err(BulkViscousHeatingError::Invalid(
                "positive finite temperature excursion required",
            ));
        }
        if partition.last().expect("nonempty partition").end_distance_m != response.length_m() {
            return Err(BulkViscousHeatingError::Invalid(
                "bulk heating must cover the complete cylinder",
            ));
        }
        let mut cells: Vec<BulkViscousCellHeating> = Vec::with_capacity(partition.len());
        let mut sources_w = BTreeMap::new();
        let mut sums = [Dd::ZERO; 3];
        let mut start = 0.0;
        let mut axis = [0.0; 3];
        let density = response
            .equilibrium()
            .absorption()
            .state()
            .acoustic()
            .density_kg_m3();
        for route in partition {
            checkpoint()?;
            let end = route.end_distance_m;
            if !end.is_finite()
                || end <= start
                || end > response.length_m()
                || sources_w.contains_key(&route.cell)
            {
                return Err(BulkViscousHeatingError::Invalid(
                    "ordered nonoverlapping depths and unique cells required",
                ));
            }
            let cell = thermal
                .cells()
                .get(&route.cell)
                .ok_or(BulkViscousHeatingError::Invalid(
                    "bulk heat destination is not active",
                ))?;
            let width = end - start;
            let delta = std::array::from_fn(|i| cell.geometry.end_m[i] - cell.geometry.start_m[i]);
            let length = finite(norm(delta))?;
            let volume = finite(response.area_m2() * width)?;
            if length <= 0.0
                || volume <= 0.0
                || !close(length, width)
                || !close(cell.geometry.volume_m3, volume)
                || !close(cell.mass_kg, finite(density * volume)?)
                || cell.geometry.surface_area_m2
                    < 2.0 * response.area_m2() * (1.0 - 256.0 * f64::EPSILON)
            {
                return Err(BulkViscousHeatingError::Invalid(
                    "original full-area volume, length and fluid mass required",
                ));
            }
            if let Some(previous) = cells.last() {
                if previous.geometry.end_m != cell.geometry.start_m
                    || previous.body_id != cell.body_id
                    || norm(std::array::from_fn(|i| delta[i] - axis[i] * width))
                        > 256.0 * f64::EPSILON * width
                {
                    return Err(BulkViscousHeatingError::Invalid(
                        "bulk heat cells must adjoin on one straight physical body",
                    ));
                }
            } else {
                axis = delta.map(|v| v / length);
            }
            if thermal.material_storage_basis(cell.material_id)?
                != ThermalStorageBasis::InternalEnergy
                || !close(
                    thermal.temperature_k(cell.id)?,
                    response.equilibrium().temperature_k(),
                )
            {
                return Err(BulkViscousHeatingError::Invalid(
                    "original equilibrium temperature and internal-energy storage required",
                ));
            }
            let caloric_identity = thermal.material_caloric_identity(cell.material_id)?.ok_or(
                BulkViscousHeatingError::Invalid("original caloric identity required"),
            )?;
            let absorption = response.absorption_between(start, end)?;
            for (sum, power) in sums.iter_mut().zip([
                absorption.viscous_heat_w,
                absorption.conductive_entropy_w_k,
                absorption.thermal_availability_loss_w,
            ]) {
                *sum = *sum + Dd::from_f64(power);
            }
            sources_w.insert(cell.id, absorption.viscous_heat_w);
            cells.push(BulkViscousCellHeating {
                cell: cell.id,
                start_distance_m: start,
                end_distance_m: end,
                absorption,
                geometry: cell.geometry,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                caloric_identity,
            });
            start = end;
        }
        let whole = response.absorption_between(0.0, response.length_m())?;
        for (sum, power) in sums.into_iter().zip([
            whole.viscous_heat_w,
            whole.conductive_entropy_w_k,
            whole.thermal_availability_loss_w,
        ]) {
            power_balance(
                (sum - Dd::from_f64(power)).to_f64(),
                sum.to_f64().abs().max(power),
                cells.len(),
            )?;
        }
        let residual = Dd::from_f64(whole.left_acoustic_power_w)
            - Dd::from_f64(whole.right_acoustic_power_w)
            - sums[0]
            - sums[2];
        power_balance(
            residual.to_f64(),
            whole
                .left_acoustic_power_w
                .abs()
                .max(whole.right_acoustic_power_w.abs())
                .max(whole.acoustic_power_loss_w),
            cells.len(),
        )?;
        checkpoint()?;
        Ok(Self {
            response,
            cells,
            sources_w,
            whole,
            maximum_temperature_change_k,
        })
    }

    pub const fn response(&self) -> &'a BulkTractionAbsorptionResponse {
        self.response
    }
    pub fn cells(&self) -> &[BulkViscousCellHeating] {
        &self.cells
    }
    pub const fn maximum_temperature_change_k(&self) -> f64 {
        self.maximum_temperature_change_k
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), BulkViscousHeatingError> {
        for bound in &self.cells {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&bound.cell)
                .ok_or(BulkViscousHeatingError::Invalid(
                    "bulk heat destination is no longer active",
                ))?;
            if cell.geometry != bound.geometry
                || cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.body_id != bound.body_id
                || cell.segment_id != bound.segment_id
                || thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::InternalEnergy
                || thermal.material_caloric_identity(cell.material_id)?
                    != Some(bound.caloric_identity)
            {
                return Err(BulkViscousHeatingError::Invalid(
                    "original geometry, inventory or caloric binding changed",
                ));
            }
            let actual_k = (thermal.temperature_k(cell.id)?
                - self.response.equilibrium().temperature_k())
            .abs();
            if actual_k > self.maximum_temperature_change_k {
                return Err(BulkViscousHeatingError::TemperatureExcursion {
                    cell: cell.id,
                    actual_k,
                    maximum_k: self.maximum_temperature_change_k,
                });
            }
        }
        Ok(())
    }

    /// Apply held cycle-mean viscous watts on actual native accepted substeps.
    /// Existing conduction, caloric nonlinearity, boundaries and events retain
    /// their original owners. Additional sources never become persistent
    /// boundaries. Do not author these same watts separately on the timeline.
    /// Every candidate, including one refused for temperature excursion, is
    /// transactional; the accepted input and acoustic field remain unchanged.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BulkViscousHeatingStep<'t>, BulkViscousHeatingError> {
        checkpoint()?;
        self.validate(initial, &mut checkpoint)?;
        let observer = HeldHeat {
            owner: self,
            cell_energy: vec![Dd::ZERO; self.cells.len()],
            account: [Dd::ZERO; 5],
        };
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.validate(&thermal, &mut checkpoint)?;
        let [left, right, viscous, entropy, availability] = observer.account;
        let acoustic = BulkAcousticEnergy {
            left_acoustic_work_j: finite(left.to_f64())?,
            right_acoustic_work_j: finite(right.to_f64())?,
            viscous_heat_j: finite(viscous.to_f64())?,
            conductive_entropy_j_k: finite(entropy.to_f64())?,
            thermal_availability_loss_j: finite(availability.to_f64())?,
            balance_residual_j: finite((left - right - viscous - availability).to_f64())?,
        };
        let cell_energy_j = self
            .cells
            .iter()
            .zip(observer.cell_energy)
            .map(|(cell, energy)| Ok((cell.cell, finite(energy.to_f64())?)))
            .collect::<Result<_, BulkViscousHeatingError>>()?;
        checkpoint()?;
        Ok(BulkViscousHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

#[derive(Clone)]
struct HeldHeat<'h, 'a> {
    owner: &'h BulkViscousHeating<'a>,
    cell_energy: Vec<Dd>,
    account: [Dd; 5],
}
fn energy(power: f64, dt: f64) -> Result<Dd, BulkViscousHeatingError> {
    let joules = finite(power * dt)?;
    if power != 0.0 && joules == 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    Ok(Dd::from_f64(joules))
}
impl DepositionObserver for HeldHeat<'_, '_> {
    type Error = BulkViscousHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        Ok(self.owner.sources_w.clone())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        for (cell, total) in self.owner.cells.iter().zip(&mut self.cell_energy) {
            checkpoint()?;
            let heat = energy(cell.absorption.viscous_heat_w, dt)?;
            *total = *total + heat;
            self.account[2] = self.account[2] + heat;
            self.account[3] = self.account[3] + energy(cell.absorption.conductive_entropy_w_k, dt)?;
            self.account[4] =
                self.account[4] + energy(cell.absorption.thermal_availability_loss_w, dt)?;
        }
        self.account[0] = self.account[0] + energy(self.owner.whole.left_acoustic_power_w, dt)?;
        self.account[1] = self.account[1] + energy(self.owner.whole.right_acoustic_power_w, dt)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "bulk_viscous_heating/tests.rs"]
mod tests;
