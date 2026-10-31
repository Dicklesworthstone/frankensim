//! Variable-volume real-gas storage with original-law work and caloric changes.
//!
//! Reversible adiabatic volume changes and fixed-volume heating are separate
//! operations. The conservative volume/energy endpoint is also available to a
//! mechanical integrator that owns its actual heat and boundary-work balance.

use super::IsochoricRealGasParcel;
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, GasHeatStep,
    entropy::{GasEntropyError, GasEntropySolve, map_context},
    parcel::GasVolumeStep,
    poll, updated_point,
};
use crate::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::positive;
use fs_matdb::QueryPoint;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

/// One closed inventory whose density follows its changing physical volume.
///
/// Initial density comes from the selected PR state. Mass, composition, source
/// cards and auxiliary query axes survive every accepted update. The caller
/// admits a stable continuous homogeneous branch; no equilibrium phase split,
/// piston dynamics, shocks or heat/work time integration is inferred here.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedRealGasParcel {
    parcel: IsochoricRealGasParcel,
}

impl ClosedRealGasParcel {
    /// Resolve initial temperature/pressure and derive mass from real density.
    pub fn new(
        model: ThermochemicalPengRobinsonModel,
        point: QueryPoint,
        volume_m3: f64,
        phase: PengRobinsonPhaseV1,
    ) -> Result<Self, GasEnergyError> {
        IsochoricRealGasParcel::new(model, point, volume_m3, phase).map(Self::from_isochoric)
    }

    /// Retain an existing accepted state and inventory without re-querying sources.
    #[must_use]
    pub fn from_isochoric(parcel: IsochoricRealGasParcel) -> Self {
        Self { parcel }
    }
    /// Hold the current volume for the original isochoric thermal consumer.
    #[must_use]
    pub fn into_isochoric(self) -> IsochoricRealGasParcel {
        self.parcel
    }
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalRealGasState {
        self.parcel.state()
    }
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        self.parcel.model()
    }
    /// Conserved specimen mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.parcel.mass_kg()
    }
    /// Current physical volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.parcel.volume_m3()
    }
    /// Current real-fluid density [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.parcel.density_kg_m3()
    }
    /// Updated temperature/pressure and unchanged external material conditions.
    #[must_use]
    pub const fn point(&self) -> &QueryPoint {
        self.parcel.point()
    }

    /// Signed heat [J] at the current held volume, through the existing owner.
    pub fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        self.parcel.add_heat(heat_j, controls, cancelled)
    }

    /// Propose an endpoint from volume and total internal-energy change [J].
    ///
    /// Internal energy depends on BOTH density and temperature for this fluid.
    /// The reference departure is therefore evaluated at the original density,
    /// not the candidate density. Zero Delta U under expansion can change
    /// temperature. The enclosing consumer must account for heat and work;
    /// this method supplies neither their split nor a mechanical trajectory.
    /// All failures leave the original parcel immutable.
    pub fn at_volume_internal_energy_increment(
        &self,
        volume_m3: f64,
        energy_change_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, f64), GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let old = &self.parcel;
        let density = if volume_m3 == old.volume_m3 {
            old.density_kg_m3
        } else {
            old.mass_kg / volume_m3
        };
        let requested = energy_change_j / old.mass_kg;
        let tolerance = old.mass_kg * controls.energy_tolerance_j_kg;
        if !positive(volume_m3)
            || !positive(density)
            || !positive(tolerance)
            || !energy_change_j.is_finite()
            || !requested.is_finite()
            || (energy_change_j != 0.0 && requested == 0.0)
            || (volume_m3 != old.volume_m3 && density == old.density_kg_m3)
        {
            return Err(GasEnergyError::Invalid(
                "positive represented real-gas volume/density and finite energy change required",
            ));
        }
        let solved = old
            .model
            .reconstruct_density_internal_energy_increment_between_densities(
                &old.point,
                old.density_kg_m3,
                density,
                requested,
                controls,
                &mut cancelled,
            )?;
        let state = solved.state();
        let stored = old.mass_kg
            * old
                .model
                .specific_internal_energy_change_between_densities(
                    old.density_kg_m3,
                    old.state.temperature_k(),
                    density,
                    state.temperature_k(),
                )?;
        if !stored.is_finite() || (stored - energy_change_j).abs() > tolerance {
            return Err(GasEnergyError::Invalid(
                "real-gas volume endpoint did not retain its energy change",
            ));
        }
        let point = updated_point(&old.point, state.temperature_k(), state.pressure_pa())?;
        poll(&mut cancelled)?;
        let mut parcel = old.clone();
        parcel.point = point;
        parcel.state = solved.into_state();
        parcel.density_kg_m3 = density;
        parcel.volume_m3 = volume_m3;
        Ok((Self::from_isochoric(parcel), stored))
    }

    /// Reversible adiabatic compression/expansion at conserved mass and entropy.
    ///
    /// Work delivered by the gas is `-m Delta u` on the admitted isentrope.
    /// It is not an independent pressure-quadrature certificate or a simultaneous
    /// heating/motion update. Positive pressure determines the work sign; no
    /// ideal-gas or constant-gamma temperature rescaling is used.
    pub fn set_volume_isentropic(
        &mut self,
        volume_m3: f64,
        controls: GasEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasVolumeStep, GasEntropyError> {
        let (candidate, report) =
            self.parcel
                .at_volume_isentropic(volume_m3, controls, &mut cancelled)?;
        poll(&mut cancelled).map_err(map_context)?;
        self.parcel = candidate;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
