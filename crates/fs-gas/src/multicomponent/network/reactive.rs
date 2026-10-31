//! Finite-rate chemistry in every shock-carrying branch and compact plenum.
//!
//! Chemistry retains original fixed-U,V calorics and formation energy exactly
//! once. Optional molecular transport, valve transfers and finite-wall sources
//! participate in the same complete-duration retry and atomic publication.
//! Symmetric ordering is not a global order or entropy-stability certificate.

use super::super::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, TransportLaw,
};
pub use super::super::reactive::{ChemistryReport, ReactiveStepControl};
use super::super::reactive::{ReactiveError, react_cell};
use super::super::wall_exchange::{WallControl, WallExchangeError};
use super::super::{
    Cell, Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, MixtureGas, StepReport,
    allocated, positive, sum,
};
pub use super::coupled::NetworkSources;
use super::coupled::{SourceError, checked_totals, duplicate, initial_duration, reset};
use super::dissipation::{NetworkDissipationOperator, SplitNetworkDissipationReport};
use super::valve::{SplitNetworkValveReport, ValveError};
use super::wall_exchange::{NetworkWallOperator, SplitNetworkWallReport};
use super::{
    MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum, PlenumBoundaryWork,
    PlenumConnection, PortExchange,
};
use crate::chemistry::{IsochoricKineticsError, IsochoricKineticsSolve, IsochoricReactor};
use fs_qty::{Amount, Time};
use fs_thermochem::kinetics::IdealGasKineticsV1;

mod step;

/// Physical owner of one local chemistry or caloric reconstruction refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChemistryLocation {
    /// A zero-based branch and cell in its original mesh order.
    Branch { branch: usize, cell: usize },
    /// A zero-based compact plenum in retained topology order.
    Plenum { plenum: usize },
}

/// Refusal publishes no network field, clock, port exchange or coupled report.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactiveNetworkError {
    /// Source identity, controls or state is incompatible with this model.
    Invalid(&'static str),
    /// The original source-bound reactor could not be prepared.
    Preparation(IsochoricKineticsError),
    /// Conservative network transport or aggregate storage admission refused.
    Flow(FlowError),
    /// Original molecular constitutive law or complete diffusive interval refused.
    Dissipation(DissipationError),
    /// Stationary-wall drag, heat exchange or finite-wall storage refused.
    Wall(WallExchangeError),
    /// Original finite-volume valve transfer refused.
    Valve(ValveError),
    /// One original per-cell chemistry update or energy reconstruction refused.
    Cell { location: ChemistryLocation, source: ReactiveError },
    /// All complete composition attempts failed.
    RetryLimit { attempts: usize, last: Box<ReactiveNetworkError> },
    /// Cancellation occurred before publication.
    Cancelled,
}
impl core::fmt::Display for ReactiveNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting gas network: {self:?}")
    }
}
impl std::error::Error for ReactiveNetworkError {}
impl From<FlowError> for ReactiveNetworkError {
    fn from(value: FlowError) -> Self {
        match value { FlowError::Cancelled => Self::Cancelled, other => Self::Flow(other) }
    }
}
impl From<DissipationError> for ReactiveNetworkError {
    fn from(value: DissipationError) -> Self {
        match value { DissipationError::Cancelled => Self::Cancelled, other => Self::Dissipation(other) }
    }
}
impl From<WallExchangeError> for ReactiveNetworkError {
    fn from(value: WallExchangeError) -> Self {
        match value { WallExchangeError::Cancelled => Self::Cancelled, other => Self::Wall(other) }
    }
}
impl From<ValveError> for ReactiveNetworkError {
    fn from(value: ValveError) -> Self {
        match value { ValveError::Cancelled => Self::Cancelled, other => Self::Valve(other) }
    }
}
impl From<SourceError> for ReactiveNetworkError {
    fn from(value: SourceError) -> Self {
        match value { SourceError::Wall(e) => e.into(), SourceError::Valve(e) => e.into() }
    }
}

/// One accepted reaction/transport step with any requested source halves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveNetworkReport {
    /// Exact duration used by transport and the two chemistry half steps.
    pub elapsed: Time,
    /// Complete compositions restarted from the accepted network.
    pub coupling_retries: usize,
    /// Intermediate transport-only balance, not the final post-source state.
    pub transport: NetworkStepReport,
    /// Accepted chemistry in all branch cells and finite plenums.
    pub chemistry: ChemistryReport,
    /// Both accepted molecular halves, each covering elapsed/2.
    pub dissipation: Option<SplitNetworkDissipationReport>,
    /// Both accepted wall halves and their actual storage changes.
    pub wall: Option<SplitNetworkWallReport>,
    /// Both accepted valve halves. These transfers are internal, never external heat.
    pub valves: Option<SplitNetworkValveReport>,
    /// Original complete branch and plenum inventories.
    pub before: NetworkTotals,
    /// Final complete inventories after all source stages.
    pub after: NetworkTotals,
    /// External transport, molecular transfers and heat into gas from walls.
    /// Add the reported finite-wall storage to close gas-plus-wall energy.
    /// Inter-plenum valve transfers are deliberately excluded from this account.
    pub boundary_change: NetworkTotals,
    /// Mechanical work on the gas [J], distinct from heat and chemical reference changes.
    pub boundary_work_on_gas_j: f64,
    /// Species/gauge-energy residual after external transport, pressure work and chemistry changes.
    pub conservation_residual: NetworkTotals,
    /// Total mass change minus external transport [kg].
    pub mass_residual_kg: f64,
    /// Physical chemical-total-energy change minus external exchange and pressure work [J],
    /// evaluated from differences rather than large absolute energy totals.
    pub absolute_energy_residual_j: f64,
}

/// One exact mechanism coupled to an existing species-bearing shock network.
/// Two complete network workspaces, one reactor and amount/receipt buffers are
/// prepared once. No mutable network handle bypasses the coupled transaction.
pub struct ReactingNetwork1d {
    accepted: MixtureNetwork1d,
    trial: MixtureNetwork1d,
    reactor: IsochoricReactor,
    amounts: Vec<Amount>,
    molecular_ports: Vec<PortExchange>,
    trial_molecular_ports: Vec<PortExchange>,
    last_report: Option<ReactiveNetworkReport>,
}
impl ReactingNetwork1d {
    /// Bind the mechanism to the identical original species bank and energy gauges.
    /// Every branch cell and plenum must lie in the declared chemistry interval.
    pub fn new(
        network: MixtureNetwork1d,
        mechanism: IdealGasKineticsV1,
        chemistry: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        if network.gas().source() != mechanism.thermodynamics() {
            return Err(ReactiveNetworkError::Invalid("network and mechanism have different original species banks"));
        }
        let domain = network.gas().domain();
        if chemistry.minimum_temperature_k < domain.minimum_temperature
            || chemistry.maximum_temperature_k > domain.maximum_temperature
        { return Err(ReactiveNetworkError::Invalid("chemistry interval exceeds the network caloric path")); }
        let reactor = IsochoricReactor::new(mechanism, chemistry, &mut cancelled).map_err(|error| match error {
            IsochoricKineticsError::Cancelled => ReactiveNetworkError::Cancelled,
            other => ReactiveNetworkError::Preparation(other),
        })?;
        for (branch, flow) in network.branches.iter().enumerate() {
            for cell in 0..flow.cell_count() {
                admit_temperature(network.gas(), flow.cell(cell).expect("admitted cell index"),
                    chemistry, ChemistryLocation::Branch { branch, cell }, &mut cancelled)?;
            }
        }
        for (plenum, volume) in network.plenums.iter().enumerate() {
            admit_temperature(network.gas(), volume.state(), chemistry,
                ChemistryLocation::Plenum { plenum }, &mut cancelled)?;
        }
        let trial = duplicate(&network, &mut cancelled)?;
        let amounts = allocated(network.gas().species_count(), Amount::new(0.0))?;
        let empty = PortExchange {
            species_masses_to_plenum: [0.0; MAX_SPECIES], species_count: network.gas().species_count(),
            momentum_to_plenum: 0.0, energy_to_plenum: 0.0,
        };
        let molecular_ports = allocated(network.connections().len(), empty)?;
        let trial_molecular_ports = allocated(network.connections().len(), empty)?;
        poll(&mut cancelled)?;
        Ok(Self { accepted: network, trial, reactor, amounts, molecular_ports,
            trial_molecular_ports, last_report: None })
    }
    /// Exact common source bank, caloric path and original species references.
    pub fn gas(&self) -> &MixtureGas { self.accepted.gas() }
    /// Common accepted clock for every branch and finite volume.
    pub fn time(&self) -> Time { self.accepted.time() }
    /// Accepted branches after the complete composition, in original order.
    pub fn branches(&self) -> &[MixtureEuler1d] { self.accepted.branches() }
    /// Accepted finite volumes after the complete composition, in original order.
    pub fn plenums(&self) -> &[Plenum] { self.accepted.plenums() }
    /// Original deterministic connection order.
    pub fn connections(&self) -> &[PlenumConnection] { self.accepted.connections() }
    /// Complete accepted species and energy storage.
    pub fn totals(&self) -> NetworkTotals { self.accepted.totals() }
    /// Original source-bound finite-rate mechanism.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 { self.reactor.mechanism() }
    /// Declared chemistry tolerances and work limits.
    pub fn chemistry_controls(&self) -> IsochoricKineticsSolve { self.reactor.controls() }
    /// Most recent complete report; refused steps leave it unchanged.
    pub const fn last_report(&self) -> Option<ReactiveNetworkReport> { self.last_report }
    /// Intermediate transport port exchanges. Add molecular and valve receipts
    /// separately when reconstructing a finite plenum's complete storage change.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Sum of both accepted molecular halves, in retained connection order.
    pub fn accepted_molecular_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.and_then(|report| report.dissipation).map(|_| self.molecular_ports.as_slice())
    }
    /// Branch balances during transport only, between the chemical half steps.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Actual accepted chamber volumes and mechanical work during transport.
    /// The surrounding chemistry, molecular and wall stages do not change volume.
    pub fn accepted_plenum_boundary_work(&self) -> Option<&[PlenumBoundaryWork]> {
        self.last_report
            .map(|_| self.accepted.plenum_work.as_slice())
    }
    /// Release complete physical fields without a stale intermediate report.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }
    /// Advance C(dt/2), A(dt), C(dt/2) atomically at the accepted duration.
    pub fn step(
        &mut self, maximum_duration: Time, controls: ReactiveStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_sources(maximum_duration, controls, NetworkSources::default(), cancelled)
    }
    /// Enclose chemistry/advection in optional valve and finite-wall half steps.
    /// All supplied source operators are borrowed for this one complete step;
    /// accepted reports and wall storage publish only with the final gas state.
    pub fn step_with_sources(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        sources: NetworkSources<'_, '_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl::<ConstantTransport>(maximum_duration, controls, None,
            DissipationControl::default(), sources, None, cancelled)
    }
    /// Compose chemistry with prescribed finite-time chamber motion and port flow.
    /// Negative rates compress; one finite dV/dt [m³/s] is required per plenum.
    /// Every chemistry half uses its actual volume on that side of transport.
    pub fn step_with_volume_rates(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        volume_rates_m3_s: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_volume_rates_and_sources(
            maximum_duration,
            controls,
            volume_rates_m3_s,
            NetworkSources::default(),
            cancelled,
        )
    }
    /// Add finite-wall heat exchange to the moving reacting network transaction.
    /// Authored wall conductances and compact mixing lengths remain fixed model
    /// inputs; volume motion alone does not infer a changing contact area.
    /// Nonzero motion with valve sources explicitly refuses. Zero rates retain
    /// the original rigid composition, including its valve support, bit for bit.
    pub fn step_with_volume_rates_and_sources(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        volume_rates_m3_s: &[f64],
        sources: NetworkSources<'_, '_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            sources,
            Some(volume_rates_m3_s),
            cancelled,
        )
    }
    /// Advance D/2, C/2, A, C/2, D/2 with source-bound molecular subcycling.
    pub fn step_with_dissipation<L: TransportLaw>(
        &mut self, operator: &mut NetworkDissipationOperator<L>, maximum_duration: Time,
        controls: ReactiveStepControl, dissipation_controls: DissipationControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_dissipation_and_sources(operator, maximum_duration, controls,
            dissipation_controls, NetworkSources::default(), cancelled)
    }
    /// Advance W/2, V/2, D/2, C/2, A, C/2, D/2, V/2, W/2, omitting absent
    /// sources. Prepare all operators from the same raw network. Every half
    /// restarts from accepted storage when any stage shortens/refuses the trial.
    /// Final valve/wall temperatures must still fit the chemistry source path.
    pub fn step_with_dissipation_and_sources<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        sources: NetworkSources<'_, '_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            sources,
            None,
            cancelled,
        )
    }
    /// Compose molecular transport, chemistry, moving-volume advection and walls.
    ///
    /// All sources use the common accepted duration. Complete retries restore
    /// original fields and volumes; operator geometry and finite-wall storage
    /// publish only with the final gas state. Pressure work enters the energy
    /// balance exactly once and remains separate from wall heat and chemistry.
    /// The original compact-volume, transport-law and splitting limits remain;
    /// this API adds no temporal-order or spatial piston-shock claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_volume_rates_and_sources<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        volume_rates_m3_s: &[f64],
        sources: NetworkSources<'_, '_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            sources,
            Some(volume_rates_m3_s),
            cancelled,
        )
    }
    /// Existing wall/chemistry/advection path, with no molecular operator.
    pub fn step_with_wall(
        &mut self, wall: &mut NetworkWallOperator, maximum_duration: Time,
        controls: ReactiveStepControl, wall_controls: WallControl, cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_sources(maximum_duration, controls,
            NetworkSources { wall: Some((wall, wall_controls)), valves: None }, cancelled)
    }
    /// Existing molecular/reacting/wall API delegates to the same source composition.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_wall<L: TransportLaw>(
        &mut self, operator: &mut NetworkDissipationOperator<L>, wall: &mut NetworkWallOperator,
        maximum_duration: Time, controls: ReactiveStepControl,
        dissipation_controls: DissipationControl, wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_with_dissipation_and_sources(operator, maximum_duration, controls, dissipation_controls,
            NetworkSources { wall: Some((wall, wall_controls)), valves: None }, cancelled)
    }
}

fn react_network(
    network: &mut MixtureNetwork1d, reactor: &mut IsochoricReactor, amounts: &mut [Amount],
    duration: f64, report: &mut ChemistryReport, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveNetworkError> {
    let width = network.gas().species_count() + 2;
    for (branch, flow) in network.branches.iter_mut().enumerate() {
        for (cell, q) in flow.states.chunks_exact_mut(width).enumerate() {
            react_cell(&flow.gas, q, flow.mesh.volumes()[cell], cell, reactor, amounts,
                duration, report, cancelled)
                .map_err(|source| cell_error(ChemistryLocation::Branch { branch, cell }, source))?;
        }
    }
    let gas = &network.branches[0].gas;
    let n = gas.species_count();
    let mut q = [0.0; MAX_COMPONENTS];
    for (plenum, volume) in network.plenums.iter_mut().enumerate() {
        q[..n].copy_from_slice(&volume.state.partial_densities);
        q[n] = volume.state.momentum;
        q[n + 1] = volume.state.energy;
        react_cell(gas, &mut q[..width], volume.volume_m3, plenum, reactor, amounts,
            duration, report, cancelled)
            .map_err(|source| cell_error(ChemistryLocation::Plenum { plenum }, source))?;
        volume.state.partial_densities.copy_from_slice(&q[..n]);
        volume.state.momentum = q[n];
        volume.state.energy = q[n + 1];
    }
    checked_totals(network, cancelled)?;
    poll(cancelled)
}
fn admit_temperature(
    gas: &MixtureGas, state: Cell<'_>, controls: IsochoricKineticsSolve,
    location: ChemistryLocation, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveNetworkError> {
    let p = gas.primitive(state, cancelled).map_err(|source| cell_error(location, source.into()))?;
    if !(controls.minimum_temperature_k..=controls.maximum_temperature_k).contains(&p.temperature) {
        return Err(cell_error(location, ReactiveError::Invalid("network temperature is outside chemistry interval")));
    }
    Ok(())
}
fn cell_error(location: ChemistryLocation, source: ReactiveError) -> ReactiveNetworkError {
    match source {
        ReactiveError::Cancelled => ReactiveNetworkError::Cancelled,
        source => ReactiveNetworkError::Cell { location, source },
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReactiveNetworkError> {
    if cancelled() { Err(ReactiveNetworkError::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
