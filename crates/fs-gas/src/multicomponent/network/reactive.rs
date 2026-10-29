//! Finite-rate chemistry in every shock-carrying branch and compact plenum.
//!
//! A complete C(dt/2), network(dt), C(dt/2) composition publishes atomically.
//! Chemistry reuses the original single-cell fixed-U,V update at each actual
//! physical volume. Changing species references changes stored sensible energy,
//! never adds a second heat-release source to chemical total energy. This is
//! an estimated composition, not a global entropy or order certificate. Optional
//! molecular halves add original-source diffusion, conduction and viscosity.

use super::super::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, TransportLaw,
};
pub use super::super::reactive::{ChemistryReport, ReactiveStepControl};
use super::super::reactive::{ReactiveError, react_cell};
use super::super::{
    Cell, Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, MixtureGas, StepReport,
    allocated, positive, sum,
};
use super::coupled::{checked_totals, duplicate, initial_duration, reset};
use super::dissipation::{NetworkDissipationOperator, SplitNetworkDissipationReport};
use super::{
    MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum, PlenumConnection, PortExchange,
};
use crate::chemistry::{IsochoricKineticsError, IsochoricKineticsSolve, IsochoricReactor};
use fs_qty::{Amount, Time};
use fs_thermochem::kinetics::IdealGasKineticsV1;

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
    /// One original per-cell chemistry update or energy reconstruction refused.
    Cell {
        location: ChemistryLocation,
        source: ReactiveError,
    },
    /// All complete composition attempts failed.
    RetryLimit {
        attempts: usize,
        last: Box<ReactiveNetworkError>,
    },
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
        match value {
            FlowError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
impl From<DissipationError> for ReactiveNetworkError {
    fn from(value: DissipationError) -> Self {
        match value {
            DissipationError::Cancelled => Self::Cancelled,
            other => Self::Dissipation(other),
        }
    }
}

/// One accepted complete reaction/transport/reaction step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveNetworkReport {
    /// Exact duration used by transport and the two chemistry half steps.
    pub elapsed: Time,
    /// Number of complete compositions restarted from the accepted network.
    pub coupling_retries: usize,
    /// Transport-only balance between the chemistry halves. Its before/after
    /// states intentionally differ from this complete report's before/after.
    pub transport: NetworkStepReport,
    /// Accepted chemistry in all branch cells and all finite plenums.
    pub chemistry: ChemistryReport,
    /// Both accepted molecular halves, when requested; each spans elapsed/2.
    pub dissipation: Option<SplitNetworkDissipationReport>,
    /// Original complete branch and plenum inventories.
    pub before: NetworkTotals,
    /// Final complete branch and plenum inventories after the second chemistry half.
    pub after: NetworkTotals,
    /// External transport plus both molecular boundary transfers.
    pub boundary_change: NetworkTotals,
    /// Species/gauge-energy residual after external transport and chemistry changes.
    pub conservation_residual: NetworkTotals,
    /// Total mass change minus external transport [kg].
    pub mass_residual_kg: f64,
    /// Physical chemical-total-energy change minus external transport [J],
    /// computed from differences without subtracting large absolute totals.
    pub absolute_energy_residual_j: f64,
}

/// One exact mechanism coupled to an existing species-bearing shock network.
///
/// Two complete network workspaces and one reactor are prepared at construction.
/// Numerical steps reuse them and one dense amount buffer. All branch cells and
/// plenums react at fixed volume and momentum; chemical total energy is retained.
/// No independently mutable network handle can bypass the coupled transaction.
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
    /// Its admitted temperature window must fit inside the flow's caloric path,
    /// and every initial branch cell and compact volume must lie in that window.
    pub fn new(
        network: MixtureNetwork1d,
        mechanism: IdealGasKineticsV1,
        chemistry: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        if network.gas().source() != mechanism.thermodynamics() {
            return Err(ReactiveNetworkError::Invalid(
                "network and mechanism have different original species banks",
            ));
        }
        let domain = network.gas().domain();
        if chemistry.minimum_temperature_k < domain.minimum_temperature
            || chemistry.maximum_temperature_k > domain.maximum_temperature
        {
            return Err(ReactiveNetworkError::Invalid(
                "chemistry interval exceeds the network caloric path",
            ));
        }
        let reactor = IsochoricReactor::new(mechanism, chemistry, &mut cancelled).map_err(
            |error| match error {
                IsochoricKineticsError::Cancelled => ReactiveNetworkError::Cancelled,
                other => ReactiveNetworkError::Preparation(other),
            },
        )?;
        for (branch, flow) in network.branches.iter().enumerate() {
            for cell in 0..flow.cell_count() {
                admit_temperature(
                    network.gas(),
                    flow.cell(cell).expect("admitted cell index"),
                    chemistry,
                    ChemistryLocation::Branch { branch, cell },
                    &mut cancelled,
                )?;
            }
        }
        for (plenum, volume) in network.plenums.iter().enumerate() {
            admit_temperature(
                network.gas(),
                volume.state(),
                chemistry,
                ChemistryLocation::Plenum { plenum },
                &mut cancelled,
            )?;
        }
        let trial = duplicate(&network, &mut cancelled)?;
        let amounts = allocated(network.gas().species_count(), Amount::new(0.0))?;
        let empty = PortExchange {
            species_masses_to_plenum: [0.0; MAX_SPECIES],
            species_count: network.gas().species_count(),
            momentum_to_plenum: 0.0,
            energy_to_plenum: 0.0,
        };
        let molecular_ports = allocated(network.connections().len(), empty)?;
        let trial_molecular_ports = allocated(network.connections().len(), empty)?;
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: network,
            trial,
            reactor,
            amounts,
            molecular_ports,
            trial_molecular_ports,
            last_report: None,
        })
    }
    /// Exact common source bank, caloric path and original species references.
    pub fn gas(&self) -> &MixtureGas {
        self.accepted.gas()
    }
    /// Common accepted clock for every branch and finite volume.
    pub fn time(&self) -> Time {
        self.accepted.time()
    }
    /// Accepted branches after the complete composition, in original order.
    pub fn branches(&self) -> &[MixtureEuler1d] {
        self.accepted.branches()
    }
    /// Accepted finite volumes after the complete composition, in original order.
    pub fn plenums(&self) -> &[Plenum] {
        self.accepted.plenums()
    }
    /// Original deterministic connection order.
    pub fn connections(&self) -> &[PlenumConnection] {
        self.accepted.connections()
    }
    /// Complete accepted species and energy storage.
    pub fn totals(&self) -> NetworkTotals {
        self.accepted.totals()
    }
    /// Original source-bound finite-rate mechanism.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }
    /// Declared chemistry tolerances and work limits.
    pub fn chemistry_controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls()
    }
    /// Most recent complete report; refused steps leave it unchanged.
    pub const fn last_report(&self) -> Option<ReactiveNetworkReport> {
        self.last_report
    }
    /// Accepted internal port exchange during the transport portion of this step.
    /// Chemistry changes storage locally and introduces no extra port exchange.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report
            .map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Sum of the two accepted molecular halves, in retained connection order.
    /// This is separate from the advective exchange returned above. Refused
    /// coupled steps preserve the previously accepted molecular exchange.
    pub fn accepted_molecular_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report
            .and_then(|report| report.dissipation)
            .map(|_| self.molecular_ports.as_slice())
    }
    /// Branch balances during transport only, between the two chemistry halves.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report
            .map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Release the final physical network. Clear the intermediate transport
    /// report so it cannot be mistaken for the post-chemistry storage balance;
    /// the returned owner's accepted port/branch report accessors return None.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }

    /// Advance all owners in one C(dt/2), network(dt), C(dt/2) transaction.
    /// If heat release changes the admissible transport duration, restart both
    /// chemistry halves and transport from the original network at that duration.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            cancelled,
        )
    }

    /// Advance D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2) atomically.
    /// Molecular transport includes every connected plenum and may subcycle;
    /// every stage uses the actual duration accepted by network advection.
    /// Prepare the molecular operator from the network before moving that
    /// network into this owner. Its accepted reports also remain unchanged
    /// unless the entire five-stage composition succeeds.
    pub fn step_with_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            cancelled,
        )
    }

    fn step_impl<L: TransportLaw>(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut operator: Option<&mut NetworkDissipationOperator<L>>,
        dissipation_controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(ReactiveNetworkError::Invalid(
                "reacting network duration, CFL or retry controls",
            ));
        }
        let before = checked_totals(&self.accepted, &mut cancelled)?;
        let mut duration = initial_duration(
            &self.accepted,
            &mut self.trial,
            maximum_duration.value(),
            controls.transport.cfl,
            &mut cancelled,
        )?;
        let mut last = ReactiveNetworkError::Invalid("no representable reacting network increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.time().value();
            if !positive(duration)
                || !positive(time + duration - time)
                || !(time + duration).is_finite()
            {
                return Err(ReactiveNetworkError::Invalid(
                    "reacting network increment cannot advance the clock",
                ));
            }
            reset(&mut self.trial, &self.accepted, &mut cancelled)?;
            let mut chemistry = ChemistryReport::default();
            let attempt = (|| {
                let first = match operator.as_deref_mut() {
                    Some(operator) => {
                        let report = operator.advance_for_trial(
                            &mut self.trial,
                            Time::new(0.5 * duration),
                            dissipation_controls,
                            &mut cancelled,
                        )?;
                        for (destination, source) in self
                            .trial_molecular_ports
                            .iter_mut()
                            .zip(operator.trial_port_exchanges())
                        {
                            poll(&mut cancelled)?;
                            *destination = *source;
                        }
                        Some(report)
                    }
                    None => None,
                };
                react_network(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let transport =
                    self.trial
                        .step(Time::new(duration), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != duration {
                    return Ok((transport, false, None));
                }
                react_network(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let second = match operator.as_deref_mut() {
                    Some(operator) => {
                        let report = operator.advance_for_trial(
                            &mut self.trial,
                            Time::new(0.5 * duration),
                            dissipation_controls,
                            &mut cancelled,
                        )?;
                        for (destination, source) in self
                            .trial_molecular_ports
                            .iter_mut()
                            .zip(operator.trial_port_exchanges())
                        {
                            poll(&mut cancelled)?;
                            for k in 0..destination.species_count {
                                destination.species_masses_to_plenum[k] +=
                                    source.species_masses_to_plenum[k];
                            }
                            destination.momentum_to_plenum += source.momentum_to_plenum;
                            destination.energy_to_plenum += source.energy_to_plenum;
                            if !destination.species_masses_to_plenum[..destination.species_count]
                                .iter()
                                .all(|x| x.is_finite())
                                || !destination.momentum_to_plenum.is_finite()
                                || !destination.energy_to_plenum.is_finite()
                            {
                                return Err(ReactiveNetworkError::Invalid(
                                    "unrepresentable combined molecular port exchange",
                                ));
                            }
                        }
                        // Viscous heating after the last chemistry solve must
                        // still fit the mechanism's possibly narrower interval.
                        let chemistry_controls = self.reactor.controls();
                        for (branch, flow) in self.trial.branches.iter().enumerate() {
                            for cell in 0..flow.cell_count() {
                                admit_temperature(
                                    self.trial.gas(),
                                    flow.cell(cell).expect("admitted cell index"),
                                    chemistry_controls,
                                    ChemistryLocation::Branch { branch, cell },
                                    &mut cancelled,
                                )?;
                            }
                        }
                        for (plenum, volume) in self.trial.plenums.iter().enumerate() {
                            admit_temperature(
                                self.trial.gas(),
                                volume.state(),
                                chemistry_controls,
                                ChemistryLocation::Plenum { plenum },
                                &mut cancelled,
                            )?;
                        }
                        Some(report)
                    }
                    None => None,
                };
                let dissipation = first
                    .zip(second)
                    .map(|(first, second)| SplitNetworkDissipationReport { first, second });
                Ok::<_, ReactiveNetworkError>((transport, true, dissipation))
            })();
            match attempt {
                Ok((transport, false, _)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveNetworkError::Invalid(
                        "network transport shortened the coupled duration",
                    );
                }
                Ok((transport, true, dissipation)) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut boundary_change = transport.boundary_change;
                    if let Some(d) = dissipation {
                        boundary_change.energy +=
                            d.first.boundary_change.energy + d.second.boundary_change.energy;
                        for k in 0..n {
                            boundary_change.species_masses[k] +=
                                d.first.boundary_change.species_masses[k]
                                    + d.second.boundary_change.species_masses[k];
                        }
                    }
                    let mut residual = NetworkTotals {
                        species_masses: [0.0; MAX_SPECIES],
                        species_count: n,
                        energy: after.energy
                            - before.energy
                            - boundary_change.energy
                            - chemistry.referenced_energy_change_j,
                    };
                    for k in 0..n {
                        residual.species_masses[k] = after.species_masses[k]
                            - before.species_masses[k]
                            - boundary_change.species_masses[k]
                            - chemistry.species_change_kg[k];
                    }
                    let mass_residual_kg = sum((0..n).map(|k| {
                        after.species_masses[k]
                            - before.species_masses[k]
                            - boundary_change.species_masses[k]
                    }));
                    let absolute_energy_residual_j = sum(core::iter::once(
                        after.energy - before.energy - boundary_change.energy,
                    )
                    .chain((0..n).map(|k| {
                        (after.species_masses[k]
                            - before.species_masses[k]
                            - boundary_change.species_masses[k])
                            * self.gas().energy_references()[k]
                    })));
                    if !residual.energy.is_finite()
                        || !residual.species_masses[..n].iter().all(|x| x.is_finite())
                        || !mass_residual_kg.is_finite()
                        || !absolute_energy_residual_j.is_finite()
                        || !chemistry.entropy_change_j_k.is_finite()
                    {
                        return Err(ReactiveNetworkError::Invalid(
                            "non-finite reacting network balance",
                        ));
                    }
                    let report = ReactiveNetworkReport {
                        elapsed: Time::new(duration),
                        coupling_retries: retry,
                        transport,
                        chemistry,
                        dissipation,
                        before,
                        after,
                        boundary_change,
                        conservation_residual: residual,
                        mass_residual_kg,
                        absolute_energy_residual_j,
                    };
                    poll(&mut cancelled)?;
                    if let Some(operator) = operator.as_deref_mut() {
                        operator.publish_trial_report();
                        core::mem::swap(&mut self.molecular_ports, &mut self.trial_molecular_ports);
                    }
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(ReactiveNetworkError::Cancelled) => {
                    return Err(ReactiveNetworkError::Cancelled);
                }
                Err(error) => {
                    last = error;
                    duration *= 0.5;
                }
            }
        }
        Err(ReactiveNetworkError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn react_network(
    network: &mut MixtureNetwork1d,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveNetworkError> {
    let width = network.gas().species_count() + 2;
    for (branch, flow) in network.branches.iter_mut().enumerate() {
        for (cell, q) in flow.states.chunks_exact_mut(width).enumerate() {
            react_cell(
                &flow.gas,
                q,
                flow.mesh.volumes()[cell],
                cell,
                reactor,
                amounts,
                duration,
                report,
                cancelled,
            )
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
        react_cell(
            gas,
            &mut q[..width],
            volume.volume_m3,
            plenum,
            reactor,
            amounts,
            duration,
            report,
            cancelled,
        )
        .map_err(|source| cell_error(ChemistryLocation::Plenum { plenum }, source))?;
        volume.state.partial_densities.copy_from_slice(&q[..n]);
        volume.state.momentum = q[n];
        volume.state.energy = q[n + 1];
    }
    checked_totals(network, cancelled)?;
    poll(cancelled)
}

fn admit_temperature(
    gas: &MixtureGas,
    state: Cell<'_>,
    controls: IsochoricKineticsSolve,
    location: ChemistryLocation,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveNetworkError> {
    let p = gas
        .primitive(state, cancelled)
        .map_err(|source| cell_error(location, source.into()))?;
    if !(controls.minimum_temperature_k..=controls.maximum_temperature_k).contains(&p.temperature) {
        return Err(cell_error(
            location,
            ReactiveError::Invalid("network temperature is outside chemistry interval"),
        ));
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
    if cancelled() {
        Err(ReactiveNetworkError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
