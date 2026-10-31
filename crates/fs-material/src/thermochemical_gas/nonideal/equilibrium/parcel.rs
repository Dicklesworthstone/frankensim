//! Closed, rigid-vessel heating through the original pure-fluid equilibrium law.
//!
//! Mass and volume are conserved. Supplied heat changes internal energy, while
//! the retained NASA/PR owner determines temperature, pressure and phase amount
//! together. Phase equilibration is instantaneous in this caloric model; it does
//! not supply evaporation kinetics, wall transport, slip or an acoustic closure.

use super::{
    PureFluidEquilibriumControls, PureFluidEquilibriumError, PureFluidEquilibriumState,
    ThermochemicalPengRobinsonModel, poll,
};
use crate::thermochemical_gas::positive;
use fs_qty::Temperature;

/// Accepted heat balance and actual bounded equilibrium work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumHeatStep {
    /// Signed heat supplied to the closed vessel [J].
    pub supplied_heat_j: f64,
    /// Independently evaluated source-law internal-energy increment [J].
    pub stored_energy_change_j: f64,
    /// Stored increment minus supplied heat [J].
    pub residual_j: f64,
    /// Temperature samples, including the original state and bracket endpoints.
    pub temperature_evaluations: usize,
    /// Original EOS evaluations across all nested coexistence solves.
    pub eos_evaluations: usize,
}

/// Conserved mass and volume with one accepted, source-bound equilibrium state.
///
/// Heating and cooling may move between admitted liquid, coexistence and vapor
/// states. The original energy reference and complete phase departure travel
/// together; latent energy is never added as a second heat source.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricEquilibriumParcel {
    model: ThermochemicalPengRobinsonModel,
    state: PureFluidEquilibriumState,
    mass_kg: f64,
    volume_m3: f64,
    density_kg_m3: f64,
}

impl IsochoricEquilibriumParcel {
    /// Resolve the initial equilibrium at the specified conserved mass/volume.
    ///
    /// The controls declare the original source/phase admission interval and
    /// bounded numerical work. No transport card is required to store calorics.
    pub fn new(
        model: ThermochemicalPengRobinsonModel,
        mass_kg: f64,
        volume_m3: f64,
        temperature: Temperature,
        controls: PureFluidEquilibriumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PureFluidEquilibriumError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !positive(mass_kg) || !positive(volume_m3) {
            return Err(PureFluidEquilibriumError::Invalid(
                "positive finite equilibrium parcel mass and volume required",
            ));
        }
        let density_kg_m3 = mass_kg / volume_m3;
        if !positive(density_kg_m3) {
            return Err(PureFluidEquilibriumError::Invalid(
                "equilibrium parcel density is unrepresentable",
            ));
        }
        let state = model.evaluate_pure_equilibrium_density_temperature(
            density_kg_m3,
            temperature,
            controls,
            &mut cancelled,
        )?;
        if !(mass_kg * state.specific_internal_energy_j_kg()).is_finite() {
            return Err(PureFluidEquilibriumError::Invalid(
                "equilibrium parcel total internal energy is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            model,
            state,
            mass_kg,
            volume_m3,
            density_kg_m3,
        })
    }

    /// Original immutable thermochemical and interaction sources.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    /// Complete accepted equilibrium calorics and phase amount.
    #[must_use]
    pub const fn state(&self) -> &PureFluidEquilibriumState {
        &self.state
    }
    /// Conserved mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Fixed vessel volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Conserved density [kg/m3], retained without repeated division.
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Accepted total internal energy on the original source reference [J].
    #[must_use]
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_internal_energy_j_kg()
    }

    /// Add signed heat [J] and atomically publish the complete new equilibrium.
    ///
    /// The source integrates the ideal-reference increment before combining it
    /// with the original PR departure change. This avoids subtracting large
    /// formation energies. The returned physical increment independently checks
    /// the requested joule balance before any stored field changes. Cancellation,
    /// insufficient work, invalid source/phase trials and failed balance leave
    /// the parcel exactly unchanged.
    pub fn add_heat(
        &mut self,
        heat_j: f64,
        controls: PureFluidEquilibriumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumHeatStep, PureFluidEquilibriumError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(PureFluidEquilibriumError::Invalid(
                "finite signed heat required",
            ));
        }
        if heat_j == 0.0 {
            return Ok(EquilibriumHeatStep {
                supplied_heat_j: heat_j,
                stored_energy_change_j: 0.0,
                residual_j: 0.0,
                temperature_evaluations: 0,
                eos_evaluations: 0,
            });
        }
        let change_j_kg = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.energy_tolerance_j_kg;
        if !change_j_kg.is_finite() || change_j_kg == 0.0 || !positive(tolerance_j) {
            return Err(PureFluidEquilibriumError::Invalid(
                "specific heat input or total energy tolerance is unrepresentable",
            ));
        }
        let solved = self
            .model
            .reconstruct_pure_density_internal_energy_increment(
                self.density_kg_m3,
                Temperature::new(self.state.temperature_k()),
                change_j_kg,
                controls,
                &mut cancelled,
            )?;
        let increment =
            solved
                .accepted_increment_j_kg()
                .ok_or(PureFluidEquilibriumError::Invalid(
                    "equilibrium increment receipt is missing",
                ))?;
        let stored_energy_change_j = self.mass_kg * increment;
        let residual_j = stored_energy_change_j - heat_j;
        if !stored_energy_change_j.is_finite()
            || !residual_j.is_finite()
            || stored_energy_change_j == 0.0
            || stored_energy_change_j.signum() != heat_j.signum()
            || residual_j.abs() > tolerance_j
            || !(self.mass_kg * solved.state().specific_internal_energy_j_kg()).is_finite()
            || solved.state().density_kg_m3() != self.density_kg_m3
        {
            return Err(PureFluidEquilibriumError::Invalid(
                "accepted equilibrium parcel did not close its mass/volume/joule balance",
            ));
        }
        let report = EquilibriumHeatStep {
            supplied_heat_j: heat_j,
            stored_energy_change_j,
            residual_j,
            temperature_evaluations: solved.temperature_evaluations(),
            eos_evaluations: solved.eos_evaluations(),
        };
        poll(&mut cancelled)?;
        self.state = solved.into_state();
        Ok(report)
    }
}
