//! Smooth, choked, equilibrium-chemistry rocket expansion.
//!
//! The existing nozzle geometry and thermochemical Gibbs owner are composed,
//! not replaced. Every pressure trial re-equilibrates composition at the same
//! entropy and elemental inventory. Choking uses the reacting sound speed;
//! substituting frozen gamma would select the wrong throat. Formation energy
//! is already in the source enthalpy and is never added as a second heat term.
//!
//! This is the infinitely fast chemistry, neutral ideal-gas, quasi-steady,
//! quasi-one-dimensional, adiabatic inviscid limit. The supplied bracket must
//! describe an ordinary single-sonic branch; uniqueness is not certified for
//! arbitrary authored caloric laws. Matched/underexpanded operation and vacuum
//! are admitted. Reacting shocks, separation, losses, finite-rate chemistry,
//! transport coefficients and vehicle/tank evolution are not inferred. The
//! original frozen off-design nozzle remains unchanged.

use super::{RocketNozzle, STANDARD_GRAVITY_M_S2};
use fs_qty::{Amount, Pressure, Temperature};
use fs_thermochem::{SpeciesId, equilibrium::composition::{
    EquilibriumEntropySolveV1, EquilibriumExpansionErrorV1, EquilibriumThermodynamicStateV1,
    IdealGasEquilibriumV1, isentrope::{EquilibriumIsentropeSampleV1, EquilibriumIsentropeV1},
}};

/// Explicit source, numerical-residual and total nested-work ceilings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumRocketNozzleSolve {
    /// Original composition and entropy inversion controls. The entropy
    /// tolerance is extensive [J/K] for the supplied feed amount, not J/(kg K).
    pub thermochemistry: EquilibriumEntropySolveV1,
    /// Positive absolute pressure floor for BOTH throat and supersonic roots [Pa].
    /// Vacuum ambient does not ask the EOS to evaluate at zero pressure.
    pub minimum_pressure_pa: f64,
    /// Allowed absolute M_throat^2-1 residual, in (0,0.01).
    pub sonic_tolerance: f64,
    /// Allowed absolute exit/throat mass-rate ratio minus one, in (0,0.01).
    pub relative_mass_flow_tolerance: f64,
    /// Absolute h-h_chamber+v^2/2 arithmetic residual [J/kg].
    pub energy_tolerance_j_kg: f64,
    /// At most 256 bisections for each of the two pressure roots. Each sample
    /// inherits the original bounded temperature/composition/backtrack limits.
    pub maximum_pressure_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
/// Typed refusal; no partially resolved flow accompanies an error.
pub enum EquilibriumRocketNozzleError {
    /// Original thermochemical source, convergence or cancellation refusal.
    Thermochemistry(EquilibriumExpansionErrorV1),
    /// Invalid geometry, pressure or numerical controls.
    Invalid(&'static str),
    /// The supplied supported pressure floor does not bracket the needed root.
    NotBracketed(&'static str),
    /// A bounded root failed to close its physical residual.
    NoConvergence { quantity: &'static str, residual: f64 },
    /// A positive finite output or specified balance cannot be represented.
    Balance(&'static str),
    /// This smooth reacting lane does not resolve overexpanded flow.
    UnsupportedBackPressure { ambient_pressure_pa: f64, exit_pressure_pa: f64 },
    /// Cancellation outside a nested thermochemical solve.
    Cancelled,
}
impl From<EquilibriumExpansionErrorV1> for EquilibriumRocketNozzleError {
    fn from(value: EquilibriumExpansionErrorV1) -> Self { Self::Thermochemistry(value) }
}
impl core::fmt::Display for EquilibriumRocketNozzleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium rocket nozzle: {self:?}")
    }
}
impl std::error::Error for EquilibriumRocketNozzleError {}
type Error = EquilibriumRocketNozzleError;

/// One accepted flow station with its actual reacting composition and sources.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumRocketStation {
    state: EquilibriumThermodynamicStateV1,
    sample: EquilibriumIsentropeSampleV1,
    /// Speed relative to the nozzle [m/s].
    pub velocity_m_s: f64,
    /// Velocity divided by the chemical-equilibrium sound speed.
    pub mach_number: f64,
    /// Actual density times velocity [kg/(m2 s)].
    pub mass_flux_kg_m2_s: f64,
    /// Source enthalpy change plus v^2/2 [J/kg].
    pub energy_residual_j_kg: f64,
}
impl EquilibriumRocketStation {
    /// Complete original thermochemical endpoint and accepted source receipts.
    pub const fn state(&self) -> &EquilibriumThermodynamicStateV1 { &self.state }
    /// Includes the directly evaluated source-relative entropy residual [J/K].
    pub const fn sample(&self) -> &EquilibriumIsentropeSampleV1 { &self.sample }
}

/// Atomic equilibrium nozzle result. No feed or physical reservoir is consumed.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumRocketNozzleFlow {
    chamber: EquilibriumThermodynamicStateV1,
    throat: EquilibriumRocketStation,
    exit: EquilibriumRocketStation,
    /// Physical throat mass rate [kg/s], independent of arbitrary feed scale.
    pub mass_flow_kg_s: f64,
    /// Axial exhaust momentum rate [N].
    pub momentum_thrust_n: f64,
    /// Exit-minus-ambient pressure times actual exit area [N].
    pub pressure_thrust_n: f64,
    /// Sum of momentum and pressure thrust [N].
    pub thrust_n: f64,
    /// Total thrust divided by mass flow [m/s].
    pub effective_exhaust_velocity_m_s: f64,
    /// Chamber pressure times throat area divided by mass flow [m/s].
    pub characteristic_velocity_m_s: f64,
    /// Effective exhaust speed divided by standard gravity [s].
    pub specific_impulse_s: f64,
    /// Actual exit mass rate divided by throat mass rate minus one.
    pub relative_mass_flow_residual: f64,
    /// Total actual composition iterations over all successful trial queries.
    pub composition_iterations: usize,
    /// Total temperature iterations across successful pressure trials.
    pub temperature_iterations: usize,
    /// Number of SP queries, including final publication checks.
    pub pressure_samples: usize,
}
impl EquilibriumRocketNozzleFlow {
    /// Source-bearing equilibrium chamber at the declared stagnation state.
    pub const fn chamber(&self) -> &EquilibriumThermodynamicStateV1 { &self.chamber }
    /// Sonic station with reacting, rather than frozen, sound speed.
    pub const fn throat(&self) -> &EquilibriumRocketStation { &self.throat }
    /// Mass-conserving supersonic exit, or the throat for unit area ratio.
    pub const fn exit(&self) -> &EquilibriumRocketStation { &self.exit }
}

#[derive(Clone, Copy)]
struct Trial { sample: EquilibriumIsentropeSampleV1, velocity: f64, mach: f64, flux: f64, energy: f64 }
#[derive(Default)]
struct Work { composition: usize, temperature: usize, samples: usize }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }
fn sample(path: &mut EquilibriumIsentropeV1<'_>, p: f64, controls: EquilibriumRocketNozzleSolve,
    work: &mut Work, cancelled: &mut impl FnMut() -> bool) -> Result<Trial, Error> {
    poll(cancelled)?;
    let sample = path.sample(Pressure::new(p), &mut *cancelled)?;
    work.samples += 1;
    work.composition += sample.composition_iterations;
    work.temperature += sample.temperature_iterations;
    let speed_squared = -2.0 * sample.specific_enthalpy_change_j_kg;
    if !positive(speed_squared) { return Err(Error::Balance("nonpositive available nozzle energy")); }
    let velocity = fs_math::det::sqrt(speed_squared);
    let mach = velocity / sample.response.sound_speed_m_s;
    let flux = sample.response.density_kg_m3 * velocity;
    let energy = sample.specific_enthalpy_change_j_kg + 0.5 * velocity * velocity;
    if !positive(velocity) || !positive(mach) || !positive(flux)
        || !energy.is_finite() || energy.abs() > controls.energy_tolerance_j_kg {
        return Err(Error::Balance("flow arithmetic or stagnation enthalpy balance"));
    }
    Ok(Trial { sample, velocity, mach, flux, energy })
}
fn root(path: &mut EquilibriumIsentropeV1<'_>, low: f64, high: f64, mut f_low: f64, f_high: f64,
    tolerance: f64, quantity: &'static str, residual: impl Fn(Trial) -> f64,
    controls: EquilibriumRocketNozzleSolve, work: &mut Work, cancelled: &mut impl FnMut() -> bool)
    -> Result<Trial, Error> {
    if !f_low.is_finite() || !f_high.is_finite() || f_low.signum() == f_high.signum() {
        return Err(Error::NotBracketed(quantity));
    }
    let (mut lo, mut hi) = (fs_math::det::ln(low), fs_math::det::ln(high));
    let mut last = f_low;
    for _ in 0..controls.maximum_pressure_iterations {
        let mid = lo + 0.5 * (hi - lo);
        if mid == lo || mid == hi { break; }
        let trial = sample(path, fs_math::det::exp(mid), controls, work, cancelled)?;
        last = residual(trial);
        if !last.is_finite() { return Err(Error::Balance("nonfinite pressure-root residual")); }
        if last.abs() <= tolerance { return Ok(trial); }
        if last.signum() == f_low.signum() { lo = mid; f_low = last; } else { hi = mid; }
    }
    Err(Error::NoConvergence { quantity, residual: last })
}
fn publish(path: &mut EquilibriumIsentropeV1<'_>, trial: Trial, cancelled: &mut impl FnMut() -> bool)
    -> Result<EquilibriumRocketStation, Error> {
    let state = path.publish_current(&mut *cancelled)?;
    if state.equilibrium().temperature() != trial.sample.temperature
        || state.equilibrium().pressure() != trial.sample.pressure
        || state.response() != trial.sample.response {
        return Err(Error::Balance("accepted source state differs from the flow trial"));
    }
    Ok(EquilibriumRocketStation { state, sample: trial.sample, velocity_m_s: trial.velocity,
        mach_number: trial.mach, mass_flux_kg_m2_s: trial.flux, energy_residual_j_kg: trial.energy })
}

/// Solve smooth equilibrium expansion from an explicitly declared equilibrium
/// chamber. A previous HP chamber calculation can supply its accepted T,p and
/// original feed; no combustion temperature or reactant mixture is invented.
///
/// The feed determines elemental ratios, not mass flow. The sonic throat sets
/// mdot, the supersonic exit closes that same mass flux, and total thrust is
/// mdot*v_exit+(p_exit-p_ambient)*A_exit. The same source equilibrium workspace
/// is reused across all scalar roots; only chamber/throat/exit retain receipts.
/// Any source, chemistry, pressure, balance or cancellation refusal returns no
/// partial flow. Composition may change throughout the nozzle. No transport
/// data is fabricated for product species that lack a transport source.
pub fn equilibrium_rocket_nozzle_flow(
    chemistry: &IdealGasEquilibriumV1, chamber_temperature: Temperature, chamber_pressure: Pressure,
    feed: &[(SpeciesId, Amount)], ambient_pressure_pa: f64, nozzle: RocketNozzle,
    controls: EquilibriumRocketNozzleSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumRocketNozzleFlow, Error> {
    nozzle.validate().map_err(|_| Error::Invalid("invalid native rocket nozzle geometry"))?;
    let pc = chamber_pressure.value();
    if !positive(pc) || !positive(controls.minimum_pressure_pa) || controls.minimum_pressure_pa >= pc
        || !ambient_pressure_pa.is_finite() || ambient_pressure_pa < 0.0 || ambient_pressure_pa >= pc
        || !positive(controls.energy_tolerance_j_kg)
        || !positive(controls.sonic_tolerance) || controls.sonic_tolerance >= 0.01
        || !positive(controls.relative_mass_flow_tolerance) || controls.relative_mass_flow_tolerance >= 0.01
        || !(1..=256).contains(&controls.maximum_pressure_iterations) {
        return Err(Error::Invalid("invalid equilibrium nozzle pressures, residual tolerances or work ceiling"));
    }
    poll(&mut cancelled)?;
    let mut path = EquilibriumIsentropeV1::new(chemistry, chamber_temperature, chamber_pressure,
        feed, controls.thermochemistry, &mut cancelled)?;
    let chamber = path.reference_state().clone();
    let mut work = Work { composition: chamber.equilibrium().composition_iterations(), ..Work::default() };
    let minimum = sample(&mut path, controls.minimum_pressure_pa, controls, &mut work, &mut cancelled)?;
    let sonic = |v: Trial| v.mach * v.mach - 1.0;
    if sonic(minimum) <= 0.0 { return Err(Error::NotBracketed("sonic throat")); }
    let throat_trial = root(&mut path, controls.minimum_pressure_pa, pc, sonic(minimum), -1.0,
        controls.sonic_tolerance, "sonic throat", sonic, controls, &mut work, &mut cancelled)?;
    let throat = publish(&mut path, throat_trial, &mut cancelled)?;
    let rate = throat.mass_flux_kg_m2_s * nozzle.throat_area_m2;
    if !positive(rate) { return Err(Error::Balance("choked mass rate")); }
    let exit = if nozzle.exit_area_m2 == nozzle.throat_area_m2 {
        throat.clone()
    } else {
        let target = rate / nozzle.exit_area_m2;
        if !positive(target) { return Err(Error::Balance("exit mass flux target")); }
        let flux_residual = |v: Trial| v.flux / target - 1.0;
        if flux_residual(minimum) >= 0.0 { return Err(Error::NotBracketed("supersonic exit mass flux")); }
        let exit_trial = root(&mut path, controls.minimum_pressure_pa, throat_trial.sample.pressure.value(),
            flux_residual(minimum), flux_residual(throat_trial), controls.relative_mass_flow_tolerance,
            "supersonic exit mass flux", flux_residual, controls, &mut work, &mut cancelled)?;
        if exit_trial.mach <= 1.0 { return Err(Error::Balance("exit is not on the supersonic branch")); }
        publish(&mut path, exit_trial, &mut cancelled)?
    };
    let pe = exit.sample.pressure.value();
    if ambient_pressure_pa > pe {
        return Err(Error::UnsupportedBackPressure { ambient_pressure_pa, exit_pressure_pa: pe });
    }
    let mass_residual = exit.mass_flux_kg_m2_s * nozzle.exit_area_m2 / rate - 1.0;
    if !mass_residual.is_finite() || mass_residual.abs() > controls.relative_mass_flow_tolerance {
        return Err(Error::Balance("accepted exit mass flow"));
    }
    let momentum_thrust = rate * exit.velocity_m_s;
    let pressure_thrust = (pe - ambient_pressure_pa) * nozzle.exit_area_m2;
    let thrust = momentum_thrust + pressure_thrust;
    let effective_velocity = thrust / rate;
    let characteristic_velocity = pc * (nozzle.throat_area_m2 / rate);
    let specific_impulse = effective_velocity / STANDARD_GRAVITY_M_S2;
    if [momentum_thrust, thrust, effective_velocity, characteristic_velocity, specific_impulse]
        .into_iter().any(|v| !positive(v)) || !pressure_thrust.is_finite()
        || pressure_thrust < 0.0 || (pe != ambient_pressure_pa && pressure_thrust == 0.0) {
        return Err(Error::Balance("unrepresentable thrust or performance"));
    }
    poll(&mut cancelled)?;
    Ok(EquilibriumRocketNozzleFlow { chamber, throat, exit, mass_flow_kg_s: rate,
        momentum_thrust_n: momentum_thrust, pressure_thrust_n: pressure_thrust, thrust_n: thrust,
        effective_exhaust_velocity_m_s: effective_velocity, characteristic_velocity_m_s: characteristic_velocity,
        specific_impulse_s: specific_impulse, relative_mass_flow_residual: mass_residual,
        composition_iterations: work.composition, temperature_iterations: work.temperature,
        pressure_samples: work.samples })
}

#[cfg(test)]
mod tests;
