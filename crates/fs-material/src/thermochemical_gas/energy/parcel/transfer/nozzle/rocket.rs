//! Frozen-composition, thermally perfect converging-diverging rocket nozzles.
//!
//! One source-backed caloric owner supplies isentropic acceleration, choking,
//! normal shocks and subsonic pressure recovery. The model is quasi-steady,
//! quasi-one-dimensional, adiabatic and inviscid with one ordinary sonic
//! crossing. Back pressure selects the unchoked branch, an internal normal
//! shock, or supersonic discharge with external expansion/compression.
//!
//! External compression is an attached ideal-flow solution, NOT a prediction
//! of viscous separation, side loads, shock trains or plume geometry. Frozen
//! composition does not imply combustion equilibrium or real-gas validity.
//! References: NASA Glenn, Rocket Thrust and Normal Shock Wave Equations:
//! <https://www.grc.nasa.gov/www/BGH/rktthsum.html>
//! <https://www.grc.nasa.gov/www/BGH/normal.html>

mod offdesign;
/// Smooth choked expansion with composition in chemical equilibrium throughout.
pub mod equilibrium;
/// Source-backed conservation of mass, momentum and enthalpy across a normal shock.
pub mod shock;

use super::{
    ClosedGasParcel, GasEnergyError, GasEnergySolve, ResolvedThermochemicalGasState, det, invalid,
    poll, positive, source_error, updated_point,
};
use fs_qty::{Pressure, Temperature};
use shock::NormalShockFlow;

/// Standard gravity used only to express specific impulse in seconds [m/s²].
/// It is not an assumed gravitational acceleration acting on a vehicle.
pub const STANDARD_GRAVITY_M_S2: f64 = 9.80665;

/// Physical flow areas of an ideal rocket nozzle, in square metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketNozzle {
    /// Minimum physical area, strictly positive [m²].
    pub throat_area_m2: f64,
    /// Exit area, at least the throat area [m²].
    pub exit_area_m2: f64,
}
impl RocketNozzle {
    fn validate(self) -> Result<(), GasEnergyError> {
        if !positive(self.throat_area_m2)
            || !positive(self.exit_area_m2)
            || self.exit_area_m2 < self.throat_area_m2
            || !(self.exit_area_m2 / self.throat_area_m2).is_finite()
        {
            return Err(invalid(
                "rocket nozzle requires finite positive areas with exit area at least throat area",
            ));
        }
        Ok(())
    }
}

/// Explicit source-temperature, residual and root-work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketNozzleSolve {
    /// Source bracket and absolute energy/entropy-energy tolerance [J/kg].
    /// Each individual root has this iteration cap; the reported work sums
    /// all nested roots. Internal shocks use at most one outer pressure root,
    /// one shock root and one downstream exit root per outer iteration.
    pub nozzle: GasEnergySolve,
    /// Relative mass-flow closure and reservoir-normalized back-pressure tolerance.
    pub relative_mass_flow_tolerance: f64,
}
impl RocketNozzleSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.nozzle.validate()?;
        if !positive(self.relative_mass_flow_tolerance) || self.relative_mass_flow_tolerance >= 1.0
        {
            return Err(invalid(
                "rocket nozzle relative mass-flow tolerance must lie in (0,1)",
            ));
        }
        Ok(())
    }
}

/// Back-pressure-selected regime of the ideal one-dimensional nozzle.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RocketNozzleRegime {
    /// Subsonic throughout, with an exit matching the ambient pressure.
    UnchokedSubsonic,
    /// Sonic throat and supersonic/sonic exit, matched or underexpanded.
    ChokedSupersonic,
    /// Supersonic exit with external compression; no separation prediction.
    ChokedExternalCompression,
    /// Sonic throat, normal shock at/before the exit and subsonic recovery.
    ChokedNormalShock,
}

/// Accepted internal shock location expressed as physical cross-sectional area.
/// Axial location requires the caller's actual nozzle area/length geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct RocketNormalShock {
    /// Cross-sectional area at the discontinuity [m²].
    pub area_m2: f64,
    /// Actual source-resolved jump and conservation residuals.
    pub flow: NormalShockFlow,
}

/// One accepted source-backed rocket flow at the unchanged reservoir state.
#[derive(Clone, Debug, PartialEq)]
pub struct RocketNozzleFlow {
    throat: ResolvedThermochemicalGasState,
    exit: ResolvedThermochemicalGasState,
    /// Selected ideal-flow regime; this is not a viscous-separation assessment.
    pub regime: RocketNozzleRegime,
    /// Actual throat speed divided by its local frozen sound speed.
    pub throat_mach_number: f64,
    /// Internal normal shock, when that regime is selected.
    pub normal_shock: Option<RocketNormalShock>,
    /// Positive propellant mass rate from the current chamber/nozzle state [kg/s].
    pub mass_flow_kg_s: f64,
    /// Exit speed relative to the nozzle [m/s].
    pub exit_velocity_m_s: f64,
    /// Exit velocity divided by the local frozen sound speed.
    pub exit_mach_number: f64,
    /// Axial exhaust momentum contribution, mdot*v_exit [N].
    pub momentum_thrust_n: f64,
    /// Axial pressure contribution, (p_exit-p_ambient)*A_exit [N].
    pub pressure_thrust_n: f64,
    /// Total axial thrust on the vehicle, opposite exhaust direction [N].
    pub thrust_n: f64,
    /// Total thrust divided by actual propellant mass rate [m/s].
    pub effective_exhaust_velocity_m_s: f64,
    /// Characteristic velocity p_chamber*A_throat/mdot [m/s].
    pub characteristic_velocity_m_s: f64,
    /// Effective exhaust velocity divided by standard gravity [s].
    pub specific_impulse_s: f64,
    /// Accepted exit mass rate divided by throat mass rate, minus one.
    pub relative_mass_flow_residual: f64,
    /// Actual h_exit-h_chamber + v_exit²/2 [J/kg], reference independent.
    pub energy_residual_j_kg: f64,
    /// T_chamber times (actual entropy change minus the shock entropy rise) [J/kg].
    pub entropy_energy_residual_j_kg: f64,
    /// Physical entropy production, zero for the smooth isentropic regimes [J/(kg K)].
    pub entropy_increase_j_kg_k: f64,
    /// Actual exit-minus-ambient pressure for subsonic discharge [Pa].
    /// Supersonic discharge has no exit pressure-matching condition.
    pub back_pressure_residual_pa: Option<f64>,
    /// Total interior iterations over every bounded physical root.
    pub iterations: usize,
}
impl RocketNozzleFlow {
    /// Source-resolved throat, sonic only when the flow is choked.
    #[must_use]
    pub const fn throat(&self) -> &ResolvedThermochemicalGasState {
        &self.throat
    }
    /// Source-resolved exit, subsonic or supersonic according to the regime.
    #[must_use]
    pub const fn exit(&self) -> &ResolvedThermochemicalGasState {
        &self.exit
    }
}

impl ClosedGasParcel {
    /// Resolve the ideal rocket nozzle at an explicit lower ambient pressure.
    ///
    /// Bounded roots poll cancellation and use allocation-free caloric samples.
    /// Only accepted physical states resolve full transport/source data. Internal
    /// shocks conserve stagnation enthalpy while increasing entropy; downstream
    /// subsonic recovery matches back pressure. External compression does not
    /// resolve a plume or assess boundary-layer separation. Equal/higher ambient
    /// pressure has no forward-flow solution in this upwind propulsion API.
    /// No inventory, momentum or thermal state is mutated on success or failure.
    pub fn rocket_nozzle_flow(
        &self,
        ambient_pressure_pa: f64,
        nozzle: RocketNozzle,
        controls: RocketNozzleSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RocketNozzleFlow, GasEnergyError> {
        nozzle.validate()?;
        controls.validate()?;
        if !ambient_pressure_pa.is_finite()
            || ambient_pressure_pa < 0.0
            || ambient_pressure_pa >= self.state().state().pressure
        {
            return Err(invalid(
                "rocket ambient pressure must be finite, nonnegative and below chamber pressure",
            ));
        }
        poll(&mut cancelled)?;
        offdesign::flow(self, ambient_pressure_pa, nozzle, controls, &mut cancelled)
    }
}

#[cfg(test)]
mod tests;
