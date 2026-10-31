//! Closed real-gas storage with one source-bound equation of state.
//!
//! Mass is fixed, and heating holds the current volume. The original NASA
//! reference calorics and Peng--Robinson departures determine temperature,
//! pressure and every published material view together. This is a homogeneous
//! frozen-composition branch, not an equilibrium flash or a vessel-wall model.
//! An opt-in immutable isentropic volume candidate uses this same storage and
//! source owner; heating the returned parcel then holds its new volume fixed.

mod volume;

/// Simultaneous heat and volume updates with source-backed midpoint pressure work.
pub mod midpoint;

/// Conserved-mass volume changes through the original real-gas energy/entropy laws.
pub mod closed;
pub use closed::ClosedRealGasParcel;

use super::super::{ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, GasHeatStep, poll, updated_point,
};
use crate::thermochemical_gas::positive;
use fs_matdb::QueryPoint;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

/// A fixed-mass parcel retaining its real-gas material source and held geometry.
///
/// The initial PT phase choice selects the starting density. Subsequent heat
/// updates keep that density and follow the mechanically and calorically stable
/// homogeneous branch admitted by the caller's temperature bracket. A locally
/// stable state may be metastable; no phase amounts or reaction are inferred.
/// `at_volume_isentropic` constructs a separate admitted parcel at a new volume;
/// it does not mutate this parcel or change the fixed-volume heating contract.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricRealGasParcel {
    model: ThermochemicalPengRobinsonModel,
    point: QueryPoint,
    state: ResolvedThermochemicalRealGasState,
    density_kg_m3: f64,
    mass_kg: f64,
    volume_m3: f64,
}

impl IsochoricRealGasParcel {
    /// Resolve the initial PT state and derive mass from its real density.
    /// The retained ideal-reference source density never determines this mass.
    pub fn new(
        model: ThermochemicalPengRobinsonModel,
        point: QueryPoint,
        volume_m3: f64,
        phase: PengRobinsonPhaseV1,
    ) -> Result<Self, GasEnergyError> {
        if !positive(volume_m3) {
            return Err(GasEnergyError::Invalid(
                "positive finite real-gas parcel volume required",
            ));
        }
        let state = model.evaluate_real_gas(&point, phase)?;
        let density_kg_m3 = state.density_kg_m3();
        let mass_kg = density_kg_m3 * volume_m3;
        if !positive(mass_kg) {
            return Err(GasEnergyError::Invalid(
                "real-gas parcel mass is unrepresentable",
            ));
        }
        Ok(Self {
            model,
            point,
            state,
            density_kg_m3,
            mass_kg,
            volume_m3,
        })
    }

    /// Original immutable caloric, interaction and transport sources.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    /// Complete accepted thermal, acoustic and fluid-transport state.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalRealGasState {
        &self.state
    }
    /// Current held density [kg/m3], retained without repeated mass/volume rounding.
    /// Heating conserves it; an isentropic volume candidate derives its new value
    /// from the same conserved mass and its explicitly requested volume.
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Conserved mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Current held volume [m3], fixed during every `add_heat` operation.
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Accepted temperature and pressure plus the original extra source axes.
    #[must_use]
    pub const fn point(&self) -> &QueryPoint {
        &self.point
    }

    /// Add signed heat [J], atomically accepting a complete real-gas state.
    ///
    /// The source-law energy increment cancels common NASA formation references
    /// before evaluation and includes the PR residual-energy change at the same
    /// fixed density. No Cp-based temperature update or ideal pressure closure
    /// is substituted. Cancellation, source-domain failure, unstable trials and
    /// insufficient nonlinear work leave the accepted parcel unchanged.
    pub fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(GasEnergyError::Invalid("finite signed heat required"));
        }
        if heat_j == 0.0 {
            return Ok(GasHeatStep {
                supplied_heat_j: heat_j,
                stored_energy_change_j: 0.0,
                residual_j: 0.0,
                iterations: 0,
            });
        }
        let change_j_kg = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.energy_tolerance_j_kg;
        if !change_j_kg.is_finite() || change_j_kg == 0.0 || !positive(tolerance_j) {
            return Err(GasEnergyError::Invalid(
                "specific heat input or total energy tolerance is unrepresentable",
            ));
        }
        let solved = self.model.reconstruct_density_internal_energy_increment(
            &self.point,
            self.density_kg_m3,
            change_j_kg,
            controls,
            &mut cancelled,
        )?;
        let stored_energy_change_j = self.mass_kg
            * self.model.specific_internal_energy_change_at_density(
                self.density_kg_m3,
                self.state.temperature_k(),
                solved.state().temperature_k(),
            )?;
        let residual_j = stored_energy_change_j - heat_j;
        if !stored_energy_change_j.is_finite()
            || !residual_j.is_finite()
            || residual_j.abs() > tolerance_j
        {
            return Err(GasEnergyError::Invalid(
                "accepted real-gas parcel did not close the requested joule balance",
            ));
        }
        let point = updated_point(
            &self.point,
            solved.state().temperature_k(),
            solved.state().pressure_pa(),
        )?;
        poll(&mut cancelled)?;
        let report = GasHeatStep {
            supplied_heat_j: heat_j,
            stored_energy_change_j,
            residual_j,
            iterations: solved.iterations(),
        };
        self.point = point;
        self.state = solved.into_state();
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
