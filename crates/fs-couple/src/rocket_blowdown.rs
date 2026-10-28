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
//! inputs. This point-vehicle adapter has no attitude or rotating-frame
//! transformation. The separate `rigid` module connects the same source owner
//! to native attitude and exhaust spin. Neither adapter supplies slosh,
//! combustion, wall heating, or a full exhaust/world-frame energy audit.
//! `equilibrium` instead allows instantaneous composition equilibration in the
//! tank and nozzle, while sharing the same native translation and impulse audit.

/// Frame-checked atmosphere, relative-flow loads and central gravity.
pub mod environment;

/// Source-bearing gas coupled to native attitude and explicit closed-valve coast.
#[path = "rocket_blowdown/rigid.rs"]
pub mod rigid;

/// Midpoint source atmosphere, aerodynamic loads and gravity on the gas/rigid vehicle.
#[path = "rocket_blowdown/flight.rs"]
pub mod flight;

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
    fn from(error: GasBlowdownError) -> Self { Self::Gas(error) }
}
impl From<RocketError> for GasRocketError {
    fn from(error: RocketError) -> Self { Self::Dynamics(error) }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasRocketError> {
    if cancelled() { Err(GasRocketError::Cancelled) } else { Ok(()) }
}

fn admit_vehicle(mass_kg: f64, dry_mass_kg: f64, position_m: Vec3, velocity_m_s: Vec3)
    -> Result<RocketState, GasRocketError> {
    let vehicle = RocketState { position_m, velocity_m_s, dry_mass_kg, propellant_mass_kg: mass_kg };
    vehicle.advance(RocketForcing { mass_flow_kg_s: 0.0,
        effective_exhaust_velocity_m_s: Vec3::ZERO, gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::ZERO }, 0.0, || false)?;
    Ok(vehicle)
}
fn admit_environment(frame: u64, environment: GasRocketEnvironment) -> Result<(), GasRocketError> {
    if environment.inertial_frame != frame {
        return Err(GasRocketError::Invalid("rocket forcing and state must share the declared inertial frame"));
    }
    let direction = environment.thrust_direction_world;
    let norm_squared = direction.norm_squared();
    if !direction.is_finite() || !norm_squared.is_finite()
        || (norm_squared - 1.0).abs() > 64.0*f64::EPSILON {
        return Err(GasRocketError::Invalid("rocket thrust direction must be a finite unit vector"));
    }
    if !environment.gravity_m_s2.is_finite() || !environment.external_force_n.is_finite() {
        return Err(GasRocketError::Invalid("rocket gravity and external force must be finite inertial-frame vectors"));
    }
    Ok(())
}
struct PropulsionTranslation {
    mechanics: RocketStep,
    remaining_residual: f64,
    consumed_residual: f64,
    impulse_residual: Vec3,
}
// Both material owners cross this one accepted-mass/impulse boundary. It does
// not know a gas law or resample the nozzle. All dynamics stay in fs-mbd.
fn translate_propulsion(vehicle: &RocketState, remaining_mass: f64, expelled: f64, impulse: f64,
    environment: GasRocketEnvironment, duration_s: f64, mass_tolerance_kg: f64,
    cancelled: &mut impl FnMut() -> bool) -> Result<PropulsionTranslation, GasRocketError> {
    let mass_flow = expelled/duration_s;
    let exhaust_velocity = impulse/expelled;
    if !mass_flow.is_finite() || mass_flow <= 0.0
        || !exhaust_velocity.is_finite() || exhaust_velocity <= 0.0 {
        return Err(GasRocketError::Invalid("accepted rocket mass or impulse cannot be represented as step forcing"));
    }
    let direction = environment.thrust_direction_world;
    let forcing = RocketForcing { mass_flow_kg_s: mass_flow,
        effective_exhaust_velocity_m_s: direction.scale(exhaust_velocity),
        gravity_m_s2: environment.gravity_m_s2, external_force_n: environment.external_force_n };
    let mechanics = vehicle.advance(forcing, duration_s, &mut *cancelled)?;
    let remaining_residual = mechanics.state_after.propellant_mass_kg-remaining_mass;
    let consumed_residual = mechanics.burned_propellant_mass_kg-expelled;
    if !remaining_residual.is_finite() || !consumed_residual.is_finite()
        || remaining_residual.abs() > mass_tolerance_kg || consumed_residual.abs() > mass_tolerance_kg
        || mechanics.coast_duration_s != 0.0 || mechanics.state_after.dry_mass_kg != vehicle.dry_mass_kg {
        return Err(GasRocketError::InventoryBalance {
            remaining_residual_kg: remaining_residual, consumed_residual_kg: consumed_residual,
        });
    }
    let impulse_residual = mechanics.thrust_impulse_n_s.sub(direction.scale(impulse));
    let impulse_tolerance = 128.0*f64::EPSILON*impulse;
    if !impulse_residual.is_finite() || impulse_residual.x.abs() > impulse_tolerance
        || impulse_residual.y.abs() > impulse_tolerance || impulse_residual.z.abs() > impulse_tolerance {
        return Err(GasRocketError::ImpulseBalance { residual_n_s: impulse_residual });
    }
    poll(cancelled)?;
    Ok(PropulsionTranslation { mechanics, remaining_residual, consumed_residual, impulse_residual })
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
pub struct GasRocketCandidate { endpoint: GasRocketState, step: GasRocketStep }
impl GasRocketCandidate {
    /// Joint next state; no independent gas or vehicle mutation is exposed.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState { &self.endpoint }
    /// Source, mass, enthalpy and thrust accounting for this same interval.
    #[must_use]
    pub const fn step(&self) -> &GasRocketStep { &self.step }
    /// Accept both owners and their accounting without another physical solve.
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, GasRocketStep) { (self.endpoint, self.step) }
}

impl GasRocketState {
    /// Assemble a point vehicle with ALL its consumable propellant in this tank.
    /// Vehicle propellant mass is derived from the already resolved gas, never
    /// accepted as an independent input. Dry mass includes the tank/container.
    pub fn new(reservoir: ClosedGasParcel, dry_mass_kg: f64, position_m: Vec3,
        velocity_m_s: Vec3, inertial_frame: u64) -> Result<Self, GasRocketError> {
        let vehicle = admit_vehicle(reservoir.mass_kg(), dry_mass_kg, position_m, velocity_m_s)?;
        Ok(Self { reservoir, vehicle, inertial_frame })
    }
    /// Original evolving thermodynamic/source state of the propellant.
    #[must_use]
    pub const fn reservoir(&self) -> &ClosedGasParcel { &self.reservoir }
    /// Original variable-mass point-vehicle state in the declared frame.
    #[must_use]
    pub const fn vehicle(&self) -> &RocketState { &self.vehicle }
    /// Caller-owned identity of the common inertial frame.
    #[must_use]
    pub const fn inertial_frame(&self) -> u64 { self.inertial_frame }

    /// Tentatively advance gas and vehicle over one shared physical interval.
    /// Uses accepted finite-step mass and impulse, including pressure thrust
    /// once. The reservoir mass tolerance also bounds cross-owner discrepancies;
    /// floating-point impulse conversion is checked separately. A late refusal
    /// returns neither new owner.
    pub fn advance(&self, environment: GasRocketEnvironment, nozzle: RocketNozzle,
        duration_s: f64, controls: RocketBlowdownSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<GasRocketCandidate, GasRocketError> {
        poll(&mut cancelled)?;
        admit_environment(self.inertial_frame, environment)?;
        let gas = self.reservoir.rocket_blowdown(environment.ambient_pressure_pa,
            nozzle, duration_s, controls, &mut cancelled)?;
        let translated = translate_propulsion(&self.vehicle, gas.reservoir().endpoint().mass_kg(),
            gas.reservoir().step().expelled_mass_kg, gas.axial_impulse_n_s,
            environment, duration_s, controls.reservoir.mass_tolerance_kg, &mut cancelled)?;
        poll(&mut cancelled)?;
        let (reservoir, gas, midpoint_nozzle) = gas.into_parts();
        let endpoint = Self { reservoir, vehicle: translated.mechanics.state_after, inertial_frame: self.inertial_frame };
        Ok(GasRocketCandidate { endpoint, step: GasRocketStep { duration_s, gas, midpoint_nozzle,
            mechanics: translated.mechanics, remaining_mass_residual_kg: translated.remaining_residual,
            consumed_mass_residual_kg: translated.consumed_residual,
            thrust_impulse_residual_n_s: translated.impulse_residual } })
    }
}

/// Reacting equilibrium propellant coupled to the same native point dynamics.
pub mod equilibrium {
//! One physical propellant inventory and one accepted thrust impulse. Neither
//! individual species nor chamber composition are artificially held constant.
//! This is the well-mixed, adiabatic, infinitely fast neutral chemistry limit.
//! Held inertial forcing and smooth step refinement remain explicit. There is
//! no inferred attitude, atmosphere, finite-rate chemistry or exhaust energy audit.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::{
    EquilibriumGasReservoir, time::{EquilibriumRocketBlowdownError,
        EquilibriumRocketBlowdownSolve, EquilibriumRocketBlowdownStep},
};

/// Source-bearing gas and the native vehicle, admitted and accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasRocketState {
    reservoir: EquilibriumGasReservoir,
    vehicle: RocketState,
    inertial_frame: u64,
}
/// No failed operation publishes either physical endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumGasRocketError {
    /// Original reacting reservoir/nozzle/time refusal.
    Gas(EquilibriumRocketBlowdownError),
    /// Original shared frame, dynamics, inventory or impulse refusal.
    Coupling(GasRocketError),
}
impl From<EquilibriumRocketBlowdownError> for EquilibriumGasRocketError {
    fn from(error: EquilibriumRocketBlowdownError) -> Self { Self::Gas(error) }
}
impl From<GasRocketError> for EquilibriumGasRocketError {
    fn from(error: GasRocketError) -> Self { Self::Coupling(error) }
}
impl core::fmt::Display for EquilibriumGasRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled equilibrium rocket: {self:?}")
    }
}
impl std::error::Error for EquilibriumGasRocketError {}

/// Complete source and mechanical accounting for one physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasRocketStep {
    /// Duration shared by the original owners [s].
    pub duration_s: f64,
    /// Actual reacting mass, energy and nozzle impulse; None means an explicitly
    /// CLOSED valve, not an attempted nozzle solve that happened to fail.
    pub propulsion: Option<EquilibriumRocketBlowdownStep>,
    /// Original analytic variable-mass/closed-valve translation report.
    pub mechanics: RocketStep,
    /// Mechanical remaining propellant minus physical tank inventory [kg].
    pub remaining_mass_residual_kg: f64,
    /// Mechanical consumed propellant minus actual tank withdrawal [kg].
    pub consumed_mass_residual_kg: f64,
    /// Mechanical thrust impulse minus actual directed nozzle impulse [N s].
    pub thrust_impulse_residual_n_s: Vec3,
}
/// Atomic next reservoir/vehicle and their matching reports.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasRocketCandidate {
    endpoint: EquilibriumGasRocketState,
    step: EquilibriumGasRocketStep,
}
impl EquilibriumGasRocketCandidate {
    /// Joint physical endpoint, with original immutable material sources retained.
    pub const fn endpoint(&self) -> &EquilibriumGasRocketState { &self.endpoint }
    /// Accepted interval's actual accounts.
    pub const fn step(&self) -> &EquilibriumGasRocketStep { &self.step }
    /// Accept without reevaluating either physical owner.
    pub fn into_parts(self) -> (EquilibriumGasRocketState, EquilibriumGasRocketStep) { (self.endpoint, self.step) }
}
impl EquilibriumGasRocketState {
    /// Assemble ALL consumable propellant from this source-bearing finite tank.
    /// No separate vehicle propellant mass is accepted. Dry mass includes its tank.
    pub fn new(reservoir: EquilibriumGasReservoir, dry_mass_kg: f64, position_m: Vec3,
        velocity_m_s: Vec3, inertial_frame: u64) -> Result<Self, EquilibriumGasRocketError> {
        let vehicle = admit_vehicle(reservoir.mass_kg(), dry_mass_kg, position_m, velocity_m_s)?;
        Ok(Self { reservoir, vehicle, inertial_frame })
    }
    /// Original physical tank inventory and its source-bearing thermodynamics.
    pub const fn reservoir(&self) -> &EquilibriumGasReservoir { &self.reservoir }
    /// Original point-vehicle state in the declared inertial frame.
    pub const fn vehicle(&self) -> &RocketState { &self.vehicle }
    /// Caller-declared common inertial-frame identity.
    pub const fn inertial_frame(&self) -> u64 { self.inertial_frame }

    /// Reacting finite-tank burn and native variable-mass motion in one candidate.
    /// Lower only actual dm/dt and impulse/dm; chamber/exit endpoint resampling
    /// cannot replace those accepted extensive exchanges. Gravity and other force
    /// stay separate from propulsion. The original gas energy is exported once;
    /// this method does not add mechanical kinetic energy as a second tank debit.
    pub fn advance(&self, environment: GasRocketEnvironment, nozzle: RocketNozzle,
        duration_s: f64, controls: EquilibriumRocketBlowdownSolve,
        mut cancelled: impl FnMut() -> bool) -> Result<EquilibriumGasRocketCandidate, EquilibriumGasRocketError> {
        poll(&mut cancelled)?;
        admit_environment(self.inertial_frame, environment)?;
        let gas = self.reservoir.rocket_blowdown(environment.ambient_pressure_pa, nozzle,
            duration_s, controls, &mut cancelled)?;
        let translated = translate_propulsion(&self.vehicle, gas.endpoint().mass_kg(),
            gas.step().gas.expelled_mass_kg, gas.step().axial_impulse_n_s,
            environment, duration_s, controls.reservoir.mass_tolerance_kg, &mut cancelled)?;
        poll(&mut cancelled)?;
        let (reservoir, propulsion) = gas.into_parts();
        let endpoint = Self { reservoir, vehicle: translated.mechanics.state_after, inertial_frame: self.inertial_frame };
        Ok(EquilibriumGasRocketCandidate { endpoint, step: EquilibriumGasRocketStep {
            duration_s, propulsion: Some(propulsion), mechanics: translated.mechanics,
            remaining_mass_residual_kg: translated.remaining_residual,
            consumed_mass_residual_kg: translated.consumed_residual,
            thrust_impulse_residual_n_s: translated.impulse_residual } })
    }

    /// Explicit closed-valve coast with the original tank retained bit-for-bit.
    /// No nozzle or equilibrium solve is performed and ambient pressure does not
    /// drive leakage through a closed valve. The common frame, finite force and
    /// unit-direction declaration remain admitted, although direction is unused.
    /// Existing accepted mass residual is retained, not reset or accumulated.
    pub fn coast(&self, environment: GasRocketEnvironment, duration_s: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<EquilibriumGasRocketCandidate, EquilibriumGasRocketError> {
        poll(&mut cancelled)?;
        admit_environment(self.inertial_frame, environment)?;
        let mechanics = self.vehicle.advance(RocketForcing { mass_flow_kg_s: 0.0,
            effective_exhaust_velocity_m_s: Vec3::ZERO, gravity_m_s2: environment.gravity_m_s2,
            external_force_n: environment.external_force_n }, duration_s, &mut cancelled)
            .map_err(GasRocketError::from)?;
        if mechanics.burned_propellant_mass_kg != 0.0
            || mechanics.state_after.propellant_mass_kg != self.vehicle.propellant_mass_kg
            || mechanics.state_after.dry_mass_kg != self.vehicle.dry_mass_kg
            || mechanics.thrust_impulse_n_s != Vec3::ZERO {
            return Err(GasRocketError::Invalid("closed valve changed inventory or engine impulse").into());
        }
        poll(&mut cancelled)?;
        let remaining = mechanics.state_after.propellant_mass_kg-self.reservoir.mass_kg();
        let endpoint = Self { reservoir: self.reservoir.clone(), vehicle: mechanics.state_after,
            inertial_frame: self.inertial_frame };
        Ok(EquilibriumGasRocketCandidate { endpoint, step: EquilibriumGasRocketStep { duration_s,
            propulsion: None, mechanics, remaining_mass_residual_kg: remaining,
            consumed_mass_residual_kg: 0.0, thrust_impulse_residual_n_s: Vec3::ZERO } })
    }
}
}

#[cfg(test)]
mod tests;
