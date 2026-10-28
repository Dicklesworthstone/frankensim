//! Closed gas storage with material-driven heating and reversible boundary work.
//!
//! One existing isochoric parcel owns mass, geometry and the resolved state.
//! Heating holds its current volume; a separate isentropic operation changes
//! volume without heat. Composition, source selection and extra query axes
//! survive both operations. The separate [`transfer`] operation exchanges mass
//! between volumes; no piston inertia or combustion is implied. The opt-in [`midpoint`] update couples heat and volume work in
//! conservative coordinates; it is not an exact finite isentropic operation.

/// Conservative composition/model replacement at fixed mass, volume and total internal energy.
pub mod composition;
/// Work-conjugate mechanical loads using the common acoustic compressibility.
pub mod mechanics;
/// Simultaneous heating and volume work through the original material laws.
pub mod midpoint;
/// Conservative paired gas mass, species and enthalpy transfer at fixed volumes.
pub mod transfer;

use super::entropy::{GasEntropyError, GasEntropySolve, map_context};
use super::{GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
    ResolvedThermochemicalGasState, ThermochemicalGasModel, poll, positive, updated_point};
use fs_matdb::QueryPoint;

/// Closed, variable-volume gas. Storage and heating retain the existing owner;
/// no second species registry, caloric law, or acoustic medium is introduced.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedGasParcel {
    parcel: IsochoricGasParcel,
}

/// Accepted reversible adiabatic volume change. Positive work leaves the gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVolumeStep {
    /// Initial volume [m3].
    pub initial_volume_m3: f64,
    /// Accepted final volume [m3].
    pub final_volume_m3: f64,
    /// Actual `m*(u_final-u_initial)` [J].
    pub stored_energy_change_j: f64,
    /// Work delivered by the gas [J], equal to minus its stored-energy change.
    /// This is the thermodynamic work on the admitted isentrope, NOT an
    /// independent pressure-quadrature or measured mechanical-work certificate.
    pub work_by_gas_j: f64,
    /// Actual `m*(s_final-s_initial)` [J/K], bounded by the caller's tolerance.
    pub entropy_change_j_k: f64,
    /// Interior entropy-inversion iterations, excluding bounded endpoint checks.
    pub iterations: usize,
}
impl ClosedGasParcel {
    /// Resolve initial T/p and obtain mass from that same material density.
    pub fn new(model: ThermochemicalGasModel, point: QueryPoint, volume_m3: f64)
        -> Result<Self, GasEnergyError>
    {
        IsochoricGasParcel::new(model, point, volume_m3).map(Self::from_isochoric)
    }
    /// Preserve an already resolved parcel without resampling or losing sources.
    #[must_use]
    pub fn from_isochoric(parcel: IsochoricGasParcel) -> Self { Self { parcel } }
    /// Freeze the current geometry without re-evaluating any material property.
    #[must_use]
    pub fn into_isochoric(self) -> IsochoricGasParcel { self.parcel }
    /// Common gas state for momentum, thermal and acoustic consumers.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalGasState { self.parcel.state() }
    /// Conserved mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.parcel.mass_kg() }
    /// Current volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.parcel.volume_m3() }
    /// Updated T/p and all retained external material conditions.
    #[must_use]
    pub const fn point(&self) -> &QueryPoint { self.parcel.point() }
    /// Original active species, caloric/transport sources and composition.
    /// Immutable access lets reaction/flow consumers validate common laws
    /// without introducing a second gas state or changing accepted storage.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalGasModel { &self.parcel.model }
    /// Prescribed signed heat at the CURRENT held volume. The existing energy
    /// closure, source-domain checks and atomic publication are unchanged.
    pub fn add_heat(&mut self, heat_j: f64, controls: GasEnergySolve, cancelled: impl FnMut() -> bool)
        -> Result<GasHeatStep, GasEnergyError>
    {
        self.parcel.add_heat(heat_j, controls, cancelled)
    }

    /// Reversible adiabatic compression/expansion at conserved mass and entropy.
    /// The caller must admit one continuous, stable gas branch over the supplied
    /// temperature bracket, INCLUDING the initial state. The endpoint is fully
    /// re-resolved through the same material law; no constant-gamma rescaling
    /// or extrapolated source property is substituted.
    ///
    /// Heat and volume changes are separate physical operations. Sequentially
    /// composing them does not establish a time integrator for simultaneous
    /// heating and motion. Work here excludes any external-pressure load, piston
    /// inertia or friction; the mechanical consumer must account for those.
    ///
    /// Failed source queries, unrepresentable work, exhausted iteration budgets
    /// and cancellation preserve the entire parcel, including its geometry.
    pub fn set_volume_isentropic(
        &mut self, volume_m3: f64, controls: GasEntropySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasVolumeStep, GasEntropyError> {
        controls.validate()?;
        poll(&mut cancelled).map_err(map_context)?;
        let old = &self.parcel;
        let temperature = old.state.state().temperature;
        if !positive(volume_m3) || temperature < controls.minimum_temperature_k
            || temperature > controls.maximum_temperature_k
        {
            return Err(GasEntropyError::Invalid("positive finite volume and a bracket containing the initial temperature required"));
        }
        let tolerance = old.mass_kg * controls.entropy_tolerance_j_kg_k;
        if !positive(tolerance) {
            return Err(GasEntropyError::Invalid("total entropy tolerance is unrepresentable"));
        }
        let mut report = GasVolumeStep {
            initial_volume_m3: old.volume_m3, final_volume_m3: volume_m3,
            stored_energy_change_j: 0.0, work_by_gas_j: 0.0,
            entropy_change_j_k: 0.0, iterations: 0,
        };
        if volume_m3 == old.volume_m3 {
            poll(&mut cancelled).map_err(map_context)?;
            return Ok(report);
        }
        let density = old.mass_kg / volume_m3;
        if !positive(density) || density == old.mass_kg / old.volume_m3 {
            return Err(GasEntropyError::Invalid("volume change does not give a representable positive density change"));
        }
        let solved = old.model.reconstruct_density_entropy_change(&old.point, density,
            0.0, controls, &mut cancelled)?;
        let stored = old.mass_kg * (solved.state().specific_internal_energy_j_kg()
            - old.state.specific_internal_energy_j_kg());
        let entropy_change = old.mass_kg * solved.residual_j_kg_k();
        // Positive cv and pressure imply heating on compression and cooling on
        // expansion. Do not accept vanished work under an overly loose tolerance.
        if !stored.is_finite() || !entropy_change.is_finite() || entropy_change.abs() > tolerance
            || (volume_m3 < old.volume_m3 && stored <= 0.0)
            || (volume_m3 > old.volume_m3 && stored >= 0.0)
        {
            return Err(GasEntropyError::Invalid("isentropic volume work or entropy balance is unrepresentable/inconsistent"));
        }
        let point = updated_point(&old.point, solved.state().state().temperature,
            solved.state().state().pressure).map_err(map_context)?;
        report.stored_energy_change_j = stored;
        report.work_by_gas_j = -stored;
        report.entropy_change_j_k = entropy_change;
        report.iterations = solved.iterations();
        poll(&mut cancelled).map_err(map_context)?;
        self.parcel.state = solved.into_state();
        self.parcel.point = point;
        self.parcel.volume_m3 = volume_m3;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;

