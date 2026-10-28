//! Transactional finite-gas propulsion and variable-mass rocket translation.
//!
//! The material owner advances the tank and nozzle first. The dynamics owner
//! then consumes exactly that accepted expelled mass and thrust impulse over
//! the same duration: mdot=dm/dt and ve=impulse/dm. Independently resampling the
//! endpoint nozzle rate would invent a second propellant inventory and impulse.
//!
//! The tank is fixed-volume, spatially mixed and comoving with a point vehicle.
//! Frozen-gas adiabatic depletion follows the exact thermodynamic path, while
//! smooth propulsion quadrature is second order in time. The dynamics step
//! integrates its piecewise-constant effective rate, exhaust vector, gravity
//! and external force analytically. Refine changing
//! inputs. This adapter has no attitude, rotating-frame transformation, slosh,
//! combustion, wall heating, or full exhaust/world-frame energy audit.

/// Frame-checked atmosphere, relative-flow loads and central gravity.
pub mod environment;

use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::blowdown::{
    GasBlowdownError, GasBlowdownStep, RocketBlowdownSolve,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, RocketNozzleFlow,
};
use fs_mbd::Vec3;
use fs_mbd::rocket::{RocketError, RocketForcing, RocketState, RocketStep};

/// Explicit common inertial-frame forcing, held fixed over one accepted step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasRocketEnvironment {
    /// Caller-owned frame identity, exactly matching the assembled state.
    /// This is an identity check, not a transformation or inertial-frame proof.
    pub inertial_frame: u64,
    /// Unit vector in the common inertial frame pointing along vehicle thrust.
    /// A nonunit vector is refused rather than silently changing total impulse.
    pub thrust_direction_world: Vec3,
    /// Explicit nonnegative ambient pressure [Pa]; zero means vacuum.
    pub ambient_pressure_pa: f64,
    /// Gravity acceleration in the same inertial frame [m/s²].
    pub gravity_m_s2: Vec3,
    /// Other force [N], excluding gravity and the nozzle's total thrust.
    pub external_force_n: Vec3,
}

/// One assembled reservoir and vehicle. Their inventories cannot be authored
/// independently or modified through this interface.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRocketState {
    reservoir: ClosedGasParcel,
    vehicle: RocketState,
    inertial_frame: u64,
}

/// Failure before either new physical owner is published.
#[derive(Clone, Debug, PartialEq)]
pub enum GasRocketError {
    /// Invalid coupling input or unsupported direction/frame semantics.
    Invalid(&'static str),
    /// The reservoir/nozzle owner refused the same physical interval.
    Gas(GasBlowdownError),
    /// The original variable-mass dynamics owner refused its derived forcing.
    Dynamics(RocketError),
    /// Actual accepted owner inventories disagree beyond the declared tolerance.
    InventoryBalance {
        remaining_residual_kg: f64,
        consumed_residual_kg: f64,
    },
    /// Actual dynamical thrust impulse disagrees with the accepted nozzle impulse.
    ImpulseBalance { residual_n_s: Vec3 },
    /// Cancellation at the coupling boundary; owner cancellations retain theirs.
    Cancelled,
}
impl core::fmt::Display for GasRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled finite-gas rocket: {self:?}")
    }
}
impl std::error::Error for GasRocketError {}
impl From<GasBlowdownError> for GasRocketError {
    fn from(error: GasBlowdownError) -> Self {
        Self::Gas(error)
    }
}
impl From<RocketError> for GasRocketError {
    fn from(error: RocketError) -> Self {
        Self::Dynamics(error)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasRocketError> {
    if cancelled() {
        Err(GasRocketError::Cancelled)
    } else {
        Ok(())
    }
}

/// The original owners' accepted reports and the coupling's actual residuals.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRocketStep {
    /// Physical duration shared by gas, thrust and translation [s].
    pub duration_s: f64,
    /// Actual exhausted mass, species and stagnation enthalpy from the tank.
    /// Nozzle kinetic energy is already part of this enthalpy, not extra heat.
    pub gas: GasBlowdownStep,
    /// Original source-resolved midpoint expansion and performance.
    pub midpoint_nozzle: RocketNozzleFlow,
    /// Original variable-mass translation report, including actual engine impulse.
    pub mechanics: RocketStep,
    /// Vehicle remaining propellant minus reservoir remaining gas [kg].
    pub remaining_mass_residual_kg: f64,
    /// Vehicle consumed propellant minus actual expelled gas [kg].
    pub consumed_mass_residual_kg: f64,
    /// Actual mechanical thrust impulse minus the directed nozzle impulse [N s].
    pub thrust_impulse_residual_n_s: Vec3,
}

/// Complete tentative gas and mechanical endpoints, published together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRocketCandidate {
    endpoint: GasRocketState,
    step: GasRocketStep,
}
impl GasRocketCandidate {
    /// Joint next state; no independent gas or vehicle mutation is exposed.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    /// Source, mass, enthalpy and thrust accounting for this same interval.
    #[must_use]
    pub const fn step(&self) -> &GasRocketStep {
        &self.step
    }
    /// Accept both owners and their accounting without another physical solve.
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, GasRocketStep) {
        (self.endpoint, self.step)
    }
}

impl GasRocketState {
    /// Assemble a point vehicle with ALL its consumable propellant in this tank.
    /// Vehicle propellant mass is derived from the already resolved gas, never
    /// accepted as an independent input. Dry mass includes the tank/container.
    pub fn new(
        reservoir: ClosedGasParcel,
        dry_mass_kg: f64,
        position_m: Vec3,
        velocity_m_s: Vec3,
        inertial_frame: u64,
    ) -> Result<Self, GasRocketError> {
        let vehicle = RocketState {
            position_m,
            velocity_m_s,
            dry_mass_kg,
            propellant_mass_kg: reservoir.mass_kg(),
        };
        // Reuse the dynamics owner's original scalar/vector admission rules.
        vehicle.advance(
            RocketForcing {
                mass_flow_kg_s: 0.0,
                effective_exhaust_velocity_m_s: Vec3::ZERO,
                gravity_m_s2: Vec3::ZERO,
                external_force_n: Vec3::ZERO,
            },
            0.0,
            || false,
        )?;
        Ok(Self {
            reservoir,
            vehicle,
            inertial_frame,
        })
    }
    /// The original evolving thermodynamic/source state of the propellant.
    #[must_use]
    pub const fn reservoir(&self) -> &ClosedGasParcel {
        &self.reservoir
    }
    /// The original variable-mass point-vehicle state in the declared frame.
    #[must_use]
    pub const fn vehicle(&self) -> &RocketState {
        &self.vehicle
    }
    /// Caller-owned identity of the common inertial frame.
    #[must_use]
    pub const fn inertial_frame(&self) -> u64 {
        self.inertial_frame
    }

    /// Tentatively advance gas and vehicle over one shared physical interval.
    /// Uses the accepted finite-step mass and impulse, including pressure thrust
    /// once. The existing reservoir mass tolerance also bounds cross-owner mass
    /// discrepancies; floating-point impulse conversion is checked separately.
    /// Invalid directions/frames, source limits, unsupported nozzle regimes,
    /// failed balance and late cancellation return neither new owner.
    pub fn advance(
        &self,
        environment: GasRocketEnvironment,
        nozzle: RocketNozzle,
        duration_s: f64,
        controls: RocketBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasRocketCandidate, GasRocketError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame {
            return Err(GasRocketError::Invalid(
                "rocket forcing and state must share the declared inertial frame",
            ));
        }
        let direction = environment.thrust_direction_world;
        let norm_squared = direction.norm_squared();
        if !direction.is_finite()
            || !norm_squared.is_finite()
            || (norm_squared - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(GasRocketError::Invalid(
                "rocket thrust direction must be a finite unit vector",
            ));
        }
        if !environment.gravity_m_s2.is_finite() || !environment.external_force_n.is_finite() {
            return Err(GasRocketError::Invalid(
                "rocket gravity and external force must be finite inertial-frame vectors",
            ));
        }
        let gas = self.reservoir.rocket_blowdown(
            environment.ambient_pressure_pa,
            nozzle,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        let expelled = gas.reservoir().step().expelled_mass_kg;
        let impulse = gas.axial_impulse_n_s;
        let mass_flow = expelled / duration_s;
        let exhaust_velocity = impulse / expelled;
        if !mass_flow.is_finite()
            || mass_flow <= 0.0
            || !exhaust_velocity.is_finite()
            || exhaust_velocity <= 0.0
        {
            return Err(GasRocketError::Invalid(
                "accepted rocket mass or impulse cannot be represented as step forcing",
            ));
        }
        let forcing = RocketForcing {
            mass_flow_kg_s: mass_flow,
            effective_exhaust_velocity_m_s: direction.scale(exhaust_velocity),
            gravity_m_s2: environment.gravity_m_s2,
            external_force_n: environment.external_force_n,
        };
        let mechanics = self.vehicle.advance(forcing, duration_s, &mut cancelled)?;
        let remaining_residual =
            mechanics.state_after.propellant_mass_kg - gas.reservoir().endpoint().mass_kg();
        let consumed_residual = mechanics.burned_propellant_mass_kg - expelled;
        if !remaining_residual.is_finite()
            || !consumed_residual.is_finite()
            || remaining_residual.abs() > controls.reservoir.mass_tolerance_kg
            || consumed_residual.abs() > controls.reservoir.mass_tolerance_kg
            || mechanics.coast_duration_s != 0.0
            || mechanics.state_after.dry_mass_kg != self.vehicle.dry_mass_kg
        {
            return Err(GasRocketError::InventoryBalance {
                remaining_residual_kg: remaining_residual,
                consumed_residual_kg: consumed_residual,
            });
        }
        let residual = mechanics.thrust_impulse_n_s.sub(direction.scale(impulse));
        let impulse_tolerance = 128.0 * f64::EPSILON * impulse;
        if !residual.is_finite()
            || residual.x.abs() > impulse_tolerance
            || residual.y.abs() > impulse_tolerance
            || residual.z.abs() > impulse_tolerance
        {
            return Err(GasRocketError::ImpulseBalance {
                residual_n_s: residual,
            });
        }
        poll(&mut cancelled)?;
        let (reservoir, gas, midpoint_nozzle) = gas.into_parts();
        let endpoint = Self {
            reservoir,
            vehicle: mechanics.state_after,
            inertial_frame: self.inertial_frame,
        };
        Ok(GasRocketCandidate {
            endpoint,
            step: GasRocketStep {
                duration_s,
                gas,
                midpoint_nozzle,
                mechanics,
                remaining_mass_residual_kg: remaining_residual,
                consumed_mass_residual_kg: consumed_residual,
                thrust_impulse_residual_n_s: residual,
            },
        })
    }
}

#[cfg(test)]
mod tests;
