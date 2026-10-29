//! Atomic molecular/advection/molecular stepping with optional walls and valves.
//! Every source half uses the actual accepted wave duration. Whole-step retries
//! preserve accepted gas, finite walls and every separately owned source report.

use super::super::dissipation::{DissipationControl, DissipationError, TransportLaw};
use super::super::wall_exchange::{WallControl, WallExchangeError};
use super::super::{
    Error as FlowError, MAX_SPECIES, MixtureEuler1d, MixtureGas, StepControl, StepReport,
    allocated, positive, sum,
};
use super::coupled::{SourceError, checked_totals, duplicate, initial_duration, reset};
pub use super::coupled::NetworkSources;
use super::dissipation::{NetworkDissipationOperator, SplitNetworkDissipationReport};
use super::wall_exchange::{NetworkWallOperator, SplitNetworkWallReport};
use super::valve::{SplitNetworkValveReport, ValveError};
use super::{MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum, PlenumConnection, PortExchange};
use fs_qty::{Length, Time};

mod step;

/// Explicit hyperbolic, molecular and complete-composition work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipatingStepControl {
    /// Conservative network transport controls.
    pub transport: StepControl,
    /// Each molecular half must cover its complete requested duration.
    pub molecular: DissipationControl,
    /// Complete composition restarts, at most 32.
    pub maximum_coupling_retries: usize,
}
impl Default for DissipatingStepControl {
    fn default() -> Self {
        Self { transport: StepControl::default(), molecular: DissipationControl::default(), maximum_coupling_retries: 20 }
    }
}

/// Refusal preserves the accepted network, clock and all accepted reports.
#[derive(Clone, Debug, PartialEq)]
pub enum DissipatingNetworkError {
    /// Source, duration, work limit or complete physical balance is inadmissible.
    Invalid(&'static str),
    /// Original conservative network or source reconstruction refused.
    Flow(FlowError),
    /// Original molecular, Maxwell--Stefan or transport-property owner refused.
    Molecular(DissipationError),
    /// Original-source wall drag, heat or finite-wall storage update refused.
    Wall(WallExchangeError),
    /// Original finite-plenum valve source refused.
    Valve(ValveError),
    /// No complete composition remained within the retry budget.
    RetryLimit { attempts: usize, last: Box<DissipatingNetworkError> },
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for DissipatingNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dissipating gas network: {self:?}")
    }
}
impl std::error::Error for DissipatingNetworkError {}
impl From<FlowError> for DissipatingNetworkError {
    fn from(value: FlowError) -> Self {
        match value { FlowError::Cancelled => Self::Cancelled, other => Self::Flow(other) }
    }
}
impl From<DissipationError> for DissipatingNetworkError {
    fn from(value: DissipationError) -> Self {
        match value { DissipationError::Cancelled => Self::Cancelled, other => Self::Molecular(other) }
    }
}
impl From<WallExchangeError> for DissipatingNetworkError {
    fn from(value: WallExchangeError) -> Self {
        match value { WallExchangeError::Cancelled => Self::Cancelled, other => Self::Wall(other) }
    }
}
impl From<ValveError> for DissipatingNetworkError {
    fn from(value: ValveError) -> Self {
        match value { ValveError::Cancelled => Self::Cancelled, other => Self::Valve(other) }
    }
}
impl From<SourceError> for DissipatingNetworkError {
    fn from(value: SourceError) -> Self {
        match value { SourceError::Wall(e) => e.into(), SourceError::Valve(e) => e.into() }
    }
}

/// Accepted D/2, A, D/2, optionally enclosed by conservative wall/valve halves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipatingNetworkReport {
    /// Duration covered by transport and the two molecular halves together.
    pub elapsed: Time,
    /// Complete compositions restarted from accepted fields.
    pub coupling_retries: usize,
    /// Intermediate advection-only report between the molecular half steps.
    pub transport: NetworkStepReport,
    /// Accepted original-source molecular intervals before and after transport.
    pub molecular: SplitNetworkDissipationReport,
    /// Both wall halves, including finite-wall storage and reservoir heat.
    pub wall: Option<SplitNetworkWallReport>,
    /// Both valve halves. Their transfers are internal, not boundary heat.
    pub valves: Option<SplitNetworkValveReport>,
    /// Complete initial storage across branches and plenums.
    pub before: NetworkTotals,
    /// Complete final storage across branches and plenums.
    pub after: NetworkTotals,
    /// External molecular/advection exchange and heat into gas from walls.
    /// Add finite-wall storage to close gas-plus-wall energy. Internal valve
    /// transfers are excluded; they cannot create an external energy supply.
    pub boundary_change: NetworkTotals,
    /// After minus before minus external exchange, per species and energy gauge.
    pub conservation_residual: NetworkTotals,
    /// Total mass residual from the retained species balances [kg].
    pub mass_residual_kg: f64,
    /// Physical total-energy residual with every species reference included [J].
    pub absolute_energy_residual_j: f64,
}

/// A source-bound molecular law composed with a complete shockwave network.
/// Accepted/trial networks and numerical/receipt workspaces are prepared once.
/// No independently mutable handle can bypass the complete gas transaction.
pub struct DissipatingNetwork1d<L: TransportLaw> {
    accepted: MixtureNetwork1d,
    trial: MixtureNetwork1d,
    molecular: NetworkDissipationOperator<L>,
    molecular_ports: Vec<PortExchange>,
    trial_molecular_ports: Vec<PortExchange>,
    last_report: Option<DissipatingNetworkReport>,
}
impl<L: TransportLaw> DissipatingNetwork1d<L> {
    /// Bind the original source and stationary constant-area branches. Positive
    /// plenum-side mixing lengths retain the original connection order.
    pub fn new(
        network: MixtureNetwork1d, law: L, plenum_lengths: Vec<Length>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        let molecular = NetworkDissipationOperator::new(&network, law, plenum_lengths, &mut cancelled)?;
        let trial = duplicate(&network, &mut cancelled)?;
        let zero = PortExchange {
            species_masses_to_plenum: [0.0; MAX_SPECIES], species_count: network.gas().species_count(),
            momentum_to_plenum: 0.0, energy_to_plenum: 0.0,
        };
        let molecular_ports = allocated(network.connections.len(), zero)?;
        let trial_molecular_ports = allocated(network.connections.len(), zero)?;
        poll(&mut cancelled)?;
        Ok(Self { accepted: network, trial, molecular, molecular_ports, trial_molecular_ports, last_report: None })
    }
    /// Original common source bank, caloric interval and species references.
    pub fn gas(&self) -> &MixtureGas { self.accepted.gas() }
    /// Common accepted clock after the complete composition.
    pub fn time(&self) -> Time { self.accepted.time() }
    /// Accepted branch states in original mesh order.
    pub fn branches(&self) -> &[MixtureEuler1d] { self.accepted.branches() }
    /// Accepted finite rigid volumes in original topology order.
    pub fn plenums(&self) -> &[Plenum] { self.accepted.plenums() }
    /// Retained deterministic connection order.
    pub fn connections(&self) -> &[PlenumConnection] { self.accepted.connections() }
    /// Complete accepted physical inventories.
    pub fn totals(&self) -> NetworkTotals { self.accepted.totals() }
    /// Most recent complete step; refused steps leave it unchanged.
    pub const fn last_report(&self) -> Option<DissipatingNetworkReport> { self.last_report }
    /// Intermediate advection exchange. Add molecular and separately supplied
    /// valve receipts when reconstructing each plenum's complete storage change.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Sum of both accepted molecular half-step exchanges, in connection order.
    pub fn accepted_molecular_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.molecular_ports.as_slice())
    }
    /// Branch reports for intermediate conservative advection only.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Release final states without a misleading intermediate transport report.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }
    /// Advance D/2, A, D/2, without optional wall or valve sources.
    pub fn step(
        &mut self, maximum_duration: Time, controls: DissipatingStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        self.step_with_sources(maximum_duration, controls, NetworkSources::default(), cancelled)
    }
    /// Advance W/2, V/2, D/2, A, D/2, V/2, W/2, omitting absent sources.
    /// Prepare the optional operators before moving the raw network into this
    /// owner. All source receipts, finite-wall storage and gas fields publish
    /// only after the complete accepted duration; all restart together on retry.
    pub fn step_with_sources(
        &mut self, maximum_duration: Time, controls: DissipatingStepControl,
        sources: NetworkSources<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        self.step_impl(maximum_duration, controls, sources, cancelled)
    }
    /// Existing wall/molecular API, using the common optional-source transaction.
    pub fn step_with_wall(
        &mut self, wall: &mut NetworkWallOperator, maximum_duration: Time,
        controls: DissipatingStepControl, wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        self.step_with_sources(maximum_duration, controls,
            NetworkSources { wall: Some((wall, wall_controls)), valves: None }, cancelled)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipatingNetworkError> {
    if cancelled() { Err(DissipatingNetworkError::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
