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

/// Actual integrated transfer, summed over both accepted valve half sweeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveExchange {
    /// Canonical species masses transferred from first to second [kg].
    pub species_to_second_kg: [f64; MAX_SPECIES],
    /// Active species count.
    pub species_count: usize,
    /// Referenced energy transferred first-to-second [J], including jet energy.
    pub energy_to_second_j: f64,
    /// Total accepted scalar mass-root iterations across the two half sweeps.
    pub iterations: usize,
    /// Number of accepted half-sweep endpoints with a sonic throat (0, 1 or 2).
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
        if !positive(self.absolute_mass_tolerance_kg) || !positive(self.relative_mass_tolerance)
            || self.relative_mass_tolerance > 1e-3 || !(1..=128).contains(&self.maximum_iterations)
            || self.maximum_coupling_retries > 32 || !positive(self.transport.cfl)
            || self.transport.cfl > 0.45 || self.transport.maximum_halvings > 32
        { return Err(ValveError::Invalid("valve root, transport or retry controls")); }
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

/// Complete accepted valve/wave composition. Valve exchange is internal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValvedNetworkReport {
    /// Duration shared by transport and both valve half sweeps together.
    pub elapsed: Time,
    /// Complete discarded compositions before this accepted step.
    pub coupling_retries: usize,
    /// Intermediate branch/plenum advection balance, not the complete endpoint.
    pub transport: NetworkStepReport,
    /// Initial complete gas storage.
    pub before: NetworkTotals,
    /// Final complete gas storage after all valve transfers.
    pub after: NetworkTotals,
    /// After minus before minus transport's external boundary change.
    pub conservation_residual: NetworkTotals,
    /// Physical-energy residual including each original species reference [J].
    pub absolute_energy_residual_j: f64,
}

/// Prepared valved shockwave network, with one atomic publication point.
///
/// Each link connects two finite plenums. Resolved branches retain their original
/// characteristic ports and wave propagation. The half sweeps run in opposite
/// retained link orders; individual endpoint mass solves are still first order.
/// Construction allocates all network/report storage. Successful steps reuse it.
pub struct ValvedNetwork1d {
    accepted: MixtureNetwork1d,
    trial: MixtureNetwork1d,
    links: Vec<ValveLink>,
    exchanges: Vec<ValveExchange>,
    trial_exchanges: Vec<ValveExchange>,
    last_report: Option<ValvedNetworkReport>,
}
impl ValvedNetwork1d {
    /// Bind explicit inter-plenum topology. Parallel links are permitted and are
    /// split in caller order. Self-links and absent plenum indices refuse.
    pub fn new(
        network: MixtureNetwork1d,
        links: Vec<ValveLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ValveError> {
        poll(&mut cancelled)?;
        if links.len() > MAX_VALVE_LINKS { return Err(ValveError::Invalid("valve link count")); }
        for link in &links {
            poll(&mut cancelled)?;
            if link.first == link.second || link.first >= network.plenums.len()
                || link.second >= network.plenums.len()
            { return Err(ValveError::Invalid("valve plenum topology")); }
        }
        let trial = duplicate(&network, &mut cancelled)?;
        let empty = ValveExchange::zero(network.gas().species_count());
        let exchanges = allocated(links.len(), empty)?;
        let trial_exchanges = allocated(links.len(), empty)?;
        poll(&mut cancelled)?;
        Ok(Self { accepted: network, trial, links, exchanges, trial_exchanges, last_report: None })
    }
    /// Exact original common source bank and energy references.
    pub fn gas(&self) -> &MixtureGas { self.accepted.gas() }
    /// Accepted original branches; no independent mutable handle is exposed.
    pub fn branches(&self) -> &[MixtureEuler1d] { self.accepted.branches() }
    /// Accepted original compact volumes, including valve-connected plenums.
    pub fn plenums(&self) -> &[Plenum] { self.accepted.plenums() }
    /// Retained deterministic valve topology and current settings.
    pub fn links(&self) -> &[ValveLink] { &self.links }
    /// Common accepted physical clock.
    pub fn time(&self) -> Time { self.accepted.time() }
    /// Complete accepted species and energy storage.
    pub fn totals(&self) -> NetworkTotals { self.accepted.totals() }
    /// Last complete accepted report, unchanged after a refused trial.
    pub const fn last_report(&self) -> Option<ValvedNetworkReport> { self.last_report }
    /// Integrated inter-plenum transfers from both accepted half sweeps.
    pub fn accepted_valve_exchanges(&self) -> Option<&[ValveExchange]> {
        self.last_report.map(|_| self.exchanges.as_slice())
    }
    /// Intermediate advective branch/plenum exchanges, excluding valve links.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Intermediate branch advection balances, before the final valve half.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Change one admitted setting exactly between accepted steps. No field or
    /// clock changes; prior reports continue to describe the previous interval.
    pub fn set_orifice(&mut self, index: usize, orifice: Orifice) -> Result<(), ValveError> {
        let link = self.links.get_mut(index).ok_or(ValveError::Invalid("valve index"))?;
        link.orifice = orifice;
        Ok(())
    }
    /// Release the physical network, clearing its intermediate advection report.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }
    /// Advance V(dt/2), A(dt), V(dt/2). Source-domain failures and a shorter
    /// acoustic duration restart the entire composition from accepted fields.
    /// Only A advances the clock; no species clipping or separate jet heat is used.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: ValveControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ValvedNetworkReport, ValveError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(maximum_duration.value()) { return Err(ValveError::Invalid("valve duration")); }
        let before = checked_totals(&self.accepted, &mut cancelled)?;
        let mut dt = initial_duration(&self.accepted, &mut self.trial,
            maximum_duration.value(), controls.transport.cfl, &mut cancelled)?;
        let mut last = ValveError::Invalid("no admitted valve composition");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.time().value();
            if !positive(dt) || !positive(time + dt - time) || !(time + dt).is_finite()
                || !positive(0.5 * dt)
            { return Err(ValveError::Invalid("valve increment cannot advance represented clock")); }
            reset(&mut self.trial, &self.accepted, &mut cancelled)?;
            self.trial_exchanges.fill(ValveExchange::zero(self.accepted.gas().species_count()));
            let attempt = (|| {
                self.sweep(0.5 * dt, false, controls, &mut cancelled)?;
                let transport = self.trial.step(Time::new(dt), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != dt { return Ok((transport, false)); }
                self.sweep(0.5 * dt, true, controls, &mut cancelled)?;
                Ok::<_, ValveError>((transport, true))
            })();
            match attempt {
                Ok((transport, false)) => {
                    dt = transport.elapsed.value();
                    last = ValveError::Invalid("advection shortened valve composition");
                }
                Ok((transport, true)) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut residual = NetworkTotals { species_masses: [0.0; MAX_SPECIES],
                        species_count: n, energy: after.energy - before.energy - transport.boundary_change.energy };
                    for k in 0..n {
                        residual.species_masses[k] = after.species_masses[k] - before.species_masses[k]
                            - transport.boundary_change.species_masses[k];
                    }
                    let physical = residual.absolute_energy(self.gas())?;
                    if !residual.energy.is_finite()
                        || !residual.species_masses[..n].iter().all(|x| x.is_finite())
                    { return Err(ValveError::Invalid("complete valve balance overflow")); }
                    let report = ValvedNetworkReport { elapsed: Time::new(dt), coupling_retries: retry,
                        transport, before, after, conservation_residual: residual,
                        absolute_energy_residual_j: physical };
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    core::mem::swap(&mut self.exchanges, &mut self.trial_exchanges);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(ValveError::Cancelled) => return Err(ValveError::Cancelled),
                Err(error) => { last = error; dt *= 0.5; }
            }
        }
        Err(ValveError::RetryLimit { attempts: controls.maximum_coupling_retries + 1, last: Box::new(last) })
    }
    fn sweep(
        &mut self, dt: f64, reverse: bool, controls: ValveControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ValveError> {
        for j in 0..self.links.len() {
            poll(cancelled)?;
            let index = if reverse { self.links.len() - 1 - j } else { j };
            let link = self.links[index];
            let update = pair::solve(self.trial.gas(),
                [&self.trial.plenums[link.first], &self.trial.plenums[link.second]],
                link.orifice, dt, controls, cancelled)?;
            self.trial_exchanges[index].add(update.exchange)?;
            let n = self.trial.gas().species_count();
            for (side, plenum) in [link.first, link.second].into_iter().enumerate() {
                self.trial.plenums[plenum].state.partial_densities.copy_from_slice(&update.states[side][..n]);
                self.trial.plenums[plenum].state.energy = update.states[side][n + 1];
            }
        }
        self.trial.last_report = None;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
