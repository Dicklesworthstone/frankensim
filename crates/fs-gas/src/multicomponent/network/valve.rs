//! Conservative, controllable valve links between finite rigid plenums.
//!
//! A converging-nozzle closure uses each donor's original temperature-dependent
//! calorics, species and stagnation enthalpy. An explicit discharge coefficient
//! and opening multiply geometric area once. The throat can choke; no constant
//! gamma, incompressible square-root surrogate or hidden leakage is introduced.
//!
//! Finite transfers follow the donor's adiabatic blowdown path and deposit the
//! identical species and energy into the receiver. A backward endpoint mass-rate
//! equation bounds the transferred inventory. This is a first-order, ordered
//! pair-split network model, not a simultaneous implicit junction solve. The
//! estimated closure has no valve inertia, lift/contact dynamics, moving-volume
//! work, resolved jet/mixing region, empirical calibration or entropy theorem.
//!
//! Reference equations: NASA Glenn, "Mass Flow Rate Equations" and
//! "Conservation of Energy"; <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/>.
//! Constant-Cp relations are test oracles only, not the production caloric law.

use super::coupled::{checked_totals, duplicate, initial_duration, reset};
use super::{MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum, PortExchange};
use crate::multicomponent::{
    Cell, Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, MixtureGas,
    Primitive, StepControl, StepReport, allocated, cell, positive, sum,
};
use fs_qty::{Area, Time};

mod flow;
mod pair;
mod network;
/// Reusable clock-preserving valve source with transactional accepted receipts.
pub mod operator;
pub use network::{ValvedNetwork1d, ValvedNetworkReport};
pub use operator::{NetworkValveOperator, SplitNetworkValveReport, ValveSolve, ValveSourceReport};

/// Maximum explicit inter-plenum links in this bounded source owner.
pub const MAX_VALVE_LINKS: usize = 4096;

/// Allowed flow direction in the link's retained first/second endpoint order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValveDirection {
    /// Either endpoint may become the donor after a pressure reversal.
    Bidirectional,
    /// Ideal check valve: reverse flow is exactly zero.
    FirstToSecond,
    /// Ideal check valve in the opposite orientation.
    SecondToFirst,
}

/// Caller-authored effective opening of a stationary converging orifice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orifice {
    area: Area,
    discharge_coefficient: f64,
    opening: f64,
    direction: ValveDirection,
}
impl Orifice {
    /// Area is nonnegative; Cd and fractional opening are in [0,1]. Zero is an
    /// exact closure. These are explicit engineering inputs, not measured data
    /// inferred from gas names. The caller owns their physical applicability.
    pub fn new(
        area: Area,
        discharge_coefficient: f64,
        opening: f64,
        direction: ValveDirection,
    ) -> Result<Self, ValveError> {
        if !area.value().is_finite() || area.value() < 0.0
            || !discharge_coefficient.is_finite() || !(0.0..=1.0).contains(&discharge_coefficient)
            || !opening.is_finite() || !(0.0..=1.0).contains(&opening)
        {
            return Err(ValveError::Invalid("orifice area, discharge coefficient or opening"));
        }
        let out = Self { area, discharge_coefficient, opening, direction };
        let effective = out.effective_area().value();
        if area.value() > 0.0 && discharge_coefficient > 0.0 && opening > 0.0 && effective == 0.0 {
            return Err(ValveError::Invalid("nonzero effective orifice area underflowed"));
        }
        Ok(out)
    }
    /// Geometric throat area before discharge and opening factors.
    pub const fn area(self) -> Area { self.area }
    /// Explicit discharge coefficient, with no hidden loss or area multiplier.
    pub const fn discharge_coefficient(self) -> f64 { self.discharge_coefficient }
    /// Retained fractional opening.
    pub const fn opening(self) -> f64 { self.opening }
    /// Retained ideal check-valve direction.
    pub const fn direction(self) -> ValveDirection { self.direction }
    /// Cd times opening times geometric area, applied once.
    pub fn effective_area(self) -> Area {
        Area::new((self.area.value() * self.opening) * self.discharge_coefficient)
    }
    /// Evaluate instantaneous flow between two source-matched stationary states.
    /// Rates are signed first-to-second. No numerical heap allocation occurs.
    pub fn flow(
        self,
        gas: &MixtureGas,
        first: Cell<'_>,
        second: Cell<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ValveFlow, ValveError> {
        flow::evaluate(gas, [first, second], self, cancelled)
    }
}

/// A stationary valve connecting two different finite plenums, not branch faces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveLink {
    /// First endpoint in the original network's plenum order.
    pub first: usize,
    /// Second endpoint in that same order.
    pub second: usize,
    /// Explicit effective opening and check-valve direction.
    pub orifice: Orifice,
}

/// Instantaneous quasi-steady nozzle exchange; positive means first-to-second.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveFlow {
    /// Canonical species mass rates [kg/s]; inactive entries are zero.
    pub species_kg_s: [f64; MAX_SPECIES],
    /// Active species count.
    pub species_count: usize,
    /// Sum of the species rates, not an independently evolved mass rate [kg/s].
    pub mass_kg_s: f64,
    /// Species-referenced stagnation enthalpy rate [W].
    pub energy_w: f64,
    /// 0 or 1 for the actual donor; None for an exactly closed/equilibrated link.
    pub donor: Option<usize>,
    /// Source-admitted exit/throat state. Speed is positive toward the receiver.
    pub exit: Option<Primitive>,
    /// True only when the admitted exit reaches its source-owned sonic limit.
    pub choked: bool,
}
impl ValveFlow {
    fn zero(n: usize) -> Self {
        Self { species_kg_s: [0.0; MAX_SPECIES], species_count: n,
            mass_kg_s: 0.0, energy_w: 0.0, donor: None, exit: None, choked: false }
    }
}

/// Actual integrated transfer, summed over the accepted valve source intervals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveExchange {
    /// Canonical species masses transferred from first to second [kg].
    pub species_to_second_kg: [f64; MAX_SPECIES],
    /// Active species count.
    pub species_count: usize,
    /// Referenced energy transferred first-to-second [J], including jet energy.
    pub energy_to_second_j: f64,
    /// Total accepted scalar mass-root iterations across the source intervals.
    pub iterations: usize,
    /// Number of accepted source endpoints with a sonic throat.
    pub choked_endpoints: usize,
    /// Maximum absolute endpoint mass-equation residual [kg].
    pub maximum_mass_residual_kg: f64,
}
impl ValveExchange {
    fn zero(n: usize) -> Self {
        Self { species_to_second_kg: [0.0; MAX_SPECIES], species_count: n,
            energy_to_second_j: 0.0, iterations: 0, choked_endpoints: 0,
            maximum_mass_residual_kg: 0.0 }
    }
    /// Total signed mass obtained from the actually transferred species.
    pub fn mass_to_second_kg(&self) -> f64 {
        sum(self.species_to_second_kg[..self.species_count].iter().copied())
    }
    fn add(&mut self, other: Self) -> Result<(), ValveError> {
        for k in 0..self.species_count {
            self.species_to_second_kg[k] += other.species_to_second_kg[k];
        }
        self.energy_to_second_j += other.energy_to_second_j;
        self.iterations += other.iterations;
        self.choked_endpoints += other.choked_endpoints;
        self.maximum_mass_residual_kg = self.maximum_mass_residual_kg.max(other.maximum_mass_residual_kg);
        if !self.energy_to_second_j.is_finite()
            || !self.species_to_second_kg[..self.species_count].iter().all(|x| x.is_finite())
        { return Err(ValveError::Invalid("composed valve exchange overflow")); }
        Ok(())
    }
}

/// Independent transport, valve root and complete-composition budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveControl {
    /// Existing conservative network wave controls.
    pub transport: StepControl,
    /// Absolute tolerance on transferred mass minus dt times endpoint flow [kg].
    pub absolute_mass_tolerance_kg: f64,
    /// Relative tolerance against the mass actually exchanged.
    pub relative_mass_tolerance: f64,
    /// Maximum safeguarded scalar mass iterations per link and half sweep.
    pub maximum_iterations: usize,
    /// Maximum whole-composition restarts after the first trial, at most 32.
    pub maximum_coupling_retries: usize,
}
impl Default for ValveControl {
    fn default() -> Self {
        Self { transport: StepControl::default(), absolute_mass_tolerance_kg: 1e-12,
            relative_mass_tolerance: 1e-9, maximum_iterations: 96, maximum_coupling_retries: 20 }
    }
}
impl ValveControl {
    fn validate(self) -> Result<(), ValveError> {
        ValveSolve::from(self).validate()?;
        if self.maximum_coupling_retries > 32 || !positive(self.transport.cfl)
            || self.transport.cfl > 0.45 || self.transport.maximum_halvings > 32
        { return Err(ValveError::Invalid("valve transport or retry controls")); }
        Ok(())
    }
}

/// A failed call preserves accepted fields, clocks and every accepted receipt.
#[derive(Clone, Debug, PartialEq)]
pub enum ValveError {
    /// Incompatible geometry, controls, source state or arithmetic.
    Invalid(&'static str),
    /// Original source/transport owner refused.
    Flow(FlowError),
    /// No representable transfer on the original caloric path.
    CaloricDomain,
    /// The bounded implicit mass relation did not converge.
    IterationLimit,
    /// All complete valve/wave compositions refused.
    RetryLimit { attempts: usize, last: Box<ValveError> },
    /// Cancellation occurred before publication.
    Cancelled,
}
impl From<FlowError> for ValveError {
    fn from(value: FlowError) -> Self {
        match value { FlowError::Cancelled => Self::Cancelled, other => Self::Flow(other) }
    }
}
impl core::fmt::Display for ValveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-plenum valve: {self:?}")
    }
}
impl std::error::Error for ValveError {}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ValveError> {
    if cancelled() { Err(ValveError::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
