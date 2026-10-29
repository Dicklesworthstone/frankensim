//! Atomic molecular/advection/molecular stepping of a complete gas network.
//!
//! Molecular diffusion, heat conduction and longitudinal viscous transport use
//! the same source-bound operator in branches and finite-plenum connections.
//! Every half step covers the actual accepted advection duration. A shortened
//! transport step restarts the entire composition from the accepted network.

use super::super::dissipation::{DissipationControl, DissipationError, TransportLaw};
use super::super::{
    Error as FlowError, MAX_SPECIES, MixtureEuler1d, MixtureGas, StepControl, StepReport,
    allocated, positive, sum,
};
use super::coupled::{checked_totals, duplicate, initial_duration, reset};
use super::dissipation::{NetworkDissipationOperator, SplitNetworkDissipationReport};
use super::{
    MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum, PlenumConnection, PortExchange,
};
use fs_qty::{Length, Time};

/// Explicit hyperbolic, molecular and complete-composition work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipatingStepControl {
    /// Conservative network transport controls.
    pub transport: StepControl,
    /// Each molecular half must cover its complete requested duration.
    pub molecular: DissipationControl,
    /// Complete D/A/D restarts, at most 32.
    pub maximum_coupling_retries: usize,
}
impl Default for DissipatingStepControl {
    fn default() -> Self {
        Self {
            transport: StepControl::default(),
            molecular: DissipationControl::default(),
            maximum_coupling_retries: 20,
        }
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
    /// No complete composition remained within the retry budget.
    RetryLimit {
        attempts: usize,
        last: Box<DissipatingNetworkError>,
    },
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
        match value {
            FlowError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
impl From<DissipationError> for DissipatingNetworkError {
    fn from(value: DissipationError) -> Self {
        match value {
            DissipationError::Cancelled => Self::Cancelled,
            other => Self::Molecular(other),
        }
    }
}

/// Accepted D(dt/2), network(dt), D(dt/2) physical storage and exchange.
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
    /// Complete initial storage across branches and plenums.
    pub before: NetworkTotals,
    /// Complete final storage across branches and plenums.
    pub after: NetworkTotals,
    /// All external molecular and advective exchange, signed into the network.
    pub boundary_change: NetworkTotals,
    /// `after - before - boundary_change`, per species and referenced energy.
    pub conservation_residual: NetworkTotals,
    /// Total mass residual, reconstructed from the retained species balances [kg].
    pub mass_residual_kg: f64,
    /// Physical total-energy residual with every species reference included [J].
    pub absolute_energy_residual_j: f64,
}

/// One source-bound molecular law composed with a complete shockwave network.
///
/// Accepted/trial networks and molecular/port workspaces are prepared once.
/// Steps neither clone branches nor allocate numerical storage. The compact
/// mixing distances belong to the molecular closure and remain explicit.
/// There is no independently mutable handle to the accepted physical network.
pub struct DissipatingNetwork1d<L: TransportLaw> {
    accepted: MixtureNetwork1d,
    trial: MixtureNetwork1d,
    molecular: NetworkDissipationOperator<L>,
    molecular_ports: Vec<PortExchange>,
    trial_molecular_ports: Vec<PortExchange>,
    last_report: Option<DissipatingNetworkReport>,
}
impl<L: TransportLaw> DissipatingNetwork1d<L> {
    /// Bind the exact original source and stationary constant-area branches.
    /// One positive plenum-side mixing distance is required per connection,
    /// in the same retained order as `network.connections()`.
    pub fn new(
        network: MixtureNetwork1d,
        law: L,
        plenum_lengths: Vec<Length>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        let molecular =
            NetworkDissipationOperator::new(&network, law, plenum_lengths, &mut cancelled)?;
        let trial = duplicate(&network, &mut cancelled)?;
        let zero = PortExchange {
            species_masses_to_plenum: [0.0; MAX_SPECIES],
            species_count: network.gas().species_count(),
            momentum_to_plenum: 0.0,
            energy_to_plenum: 0.0,
        };
        let molecular_ports = allocated(network.connections.len(), zero)?;
        let trial_molecular_ports = allocated(network.connections.len(), zero)?;
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: network,
            trial,
            molecular,
            molecular_ports,
            trial_molecular_ports,
            last_report: None,
        })
    }
    /// Original common source bank, caloric interval and species references.
    pub fn gas(&self) -> &MixtureGas {
        self.accepted.gas()
    }
    /// Common accepted clock after the complete composition.
    pub fn time(&self) -> Time {
        self.accepted.time()
    }
    /// Accepted branch states in original mesh order.
    pub fn branches(&self) -> &[MixtureEuler1d] {
        self.accepted.branches()
    }
    /// Accepted finite rigid volumes in original topology order.
    pub fn plenums(&self) -> &[Plenum] {
        self.accepted.plenums()
    }
    /// Retained deterministic connection order.
    pub fn connections(&self) -> &[PlenumConnection] {
        self.accepted.connections()
    }
    /// Complete accepted physical inventories.
    pub fn totals(&self) -> NetworkTotals {
        self.accepted.totals()
    }
    /// Most recent complete step; refused steps leave it unchanged.
    pub const fn last_report(&self) -> Option<DissipatingNetworkReport> {
        self.last_report
    }
    /// Accepted advective port exchange during the intermediate A stage only.
    /// Add `accepted_molecular_port_exchanges()` for the complete port transfer.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report
            .map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Sum of both accepted molecular half-step exchanges, in connection order.
    /// Each impulse belongs to its branch axis; unrelated axes are not summed.
    pub fn accepted_molecular_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.molecular_ports.as_slice())
    }
    /// Branch reports for intermediate conservative advection only.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report
            .map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Release final states without a misleading intermediate transport report.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }

    /// Advance every branch and plenum for one complete D/A/D transaction.
    ///
    /// Molecular stages internally subcycle their exact half duration. If
    /// heating or momentum diffusion shortens advection's accepted duration,
    /// all stages restart from accepted states at that new duration. No source
    /// term independently adds viscous heating to conservative total energy.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: DissipatingStepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(DissipatingNetworkError::Invalid(
                "molecular network duration, CFL or retry controls",
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
        let mut last = DissipatingNetworkError::Invalid("no representable network increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.time().value();
            if !positive(duration)
                || !positive(time + duration - time)
                || !(time + duration).is_finite()
            {
                return Err(DissipatingNetworkError::Invalid(
                    "molecular network increment cannot advance the clock",
                ));
            }
            reset(&mut self.trial, &self.accepted, &mut cancelled)?;
            let attempt = (|| {
                let first = self.molecular.advance_for_trial(
                    &mut self.trial,
                    Time::new(0.5 * duration),
                    controls.molecular,
                    &mut cancelled,
                )?;
                for (destination, source) in self
                    .trial_molecular_ports
                    .iter_mut()
                    .zip(self.molecular.trial_port_exchanges())
                {
                    poll(&mut cancelled)?;
                    *destination = *source;
                }
                let transport =
                    self.trial
                        .step(Time::new(duration), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != duration {
                    return Ok((transport, None));
                }
                let second = self.molecular.advance_for_trial(
                    &mut self.trial,
                    Time::new(0.5 * duration),
                    controls.molecular,
                    &mut cancelled,
                )?;
                for (total, second) in self
                    .trial_molecular_ports
                    .iter_mut()
                    .zip(self.molecular.trial_port_exchanges())
                {
                    poll(&mut cancelled)?;
                    total.momentum_to_plenum += second.momentum_to_plenum;
                    total.energy_to_plenum += second.energy_to_plenum;
                    for k in 0..total.species_count {
                        total.species_masses_to_plenum[k] += second.species_masses_to_plenum[k];
                    }
                    if !total.momentum_to_plenum.is_finite()
                        || !total.energy_to_plenum.is_finite()
                        || !total.species_masses_to_plenum[..total.species_count]
                            .iter()
                            .all(|x| x.is_finite())
                    {
                        return Err(DissipatingNetworkError::Invalid(
                            "non-finite composed molecular port exchange",
                        ));
                    }
                }
                Ok::<_, DissipatingNetworkError>((
                    transport,
                    Some(SplitNetworkDissipationReport { first, second }),
                ))
            })();
            match attempt {
                Ok((transport, None)) => {
                    duration = transport.elapsed.value();
                    last = DissipatingNetworkError::Invalid(
                        "network transport shortened the molecular composition",
                    );
                }
                Ok((transport, Some(molecular))) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut boundary_change = NetworkTotals {
                        species_masses: [0.0; MAX_SPECIES],
                        species_count: n,
                        energy: sum([
                            molecular.first.boundary_change.energy,
                            transport.boundary_change.energy,
                            molecular.second.boundary_change.energy,
                        ]),
                    };
                    let mut residual = NetworkTotals {
                        species_masses: [0.0; MAX_SPECIES],
                        species_count: n,
                        energy: after.energy - before.energy - boundary_change.energy,
                    };
                    for k in 0..n {
                        boundary_change.species_masses[k] = sum([
                            molecular.first.boundary_change.species_masses[k],
                            transport.boundary_change.species_masses[k],
                            molecular.second.boundary_change.species_masses[k],
                        ]);
                        residual.species_masses[k] = after.species_masses[k]
                            - before.species_masses[k]
                            - boundary_change.species_masses[k];
                    }
                    let mass_residual_kg = residual.mass();
                    let absolute_energy_residual_j =
                        sum(core::iter::once(residual.energy).chain((0..n).map(|k| {
                            residual.species_masses[k] * self.gas().energy_references()[k]
                        })));
                    if !boundary_change.energy.is_finite()
                        || !boundary_change.species_masses[..n]
                            .iter()
                            .all(|x| x.is_finite())
                        || !residual.energy.is_finite()
                        || !residual.species_masses[..n].iter().all(|x| x.is_finite())
                        || !mass_residual_kg.is_finite()
                        || !absolute_energy_residual_j.is_finite()
                    {
                        return Err(DissipatingNetworkError::Invalid(
                            "non-finite complete molecular network balance",
                        ));
                    }
                    let report = DissipatingNetworkReport {
                        elapsed: Time::new(duration),
                        coupling_retries: retry,
                        transport,
                        molecular,
                        before,
                        after,
                        boundary_change,
                        conservation_residual: residual,
                        mass_residual_kg,
                        absolute_energy_residual_j,
                    };
                    poll(&mut cancelled)?;
                    self.molecular.publish_trial_report();
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    core::mem::swap(&mut self.molecular_ports, &mut self.trial_molecular_ports);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(DissipatingNetworkError::Cancelled) => {
                    return Err(DissipatingNetworkError::Cancelled);
                }
                Err(error) => {
                    last = error;
                    duration *= 0.5;
                }
            }
        }
        Err(DissipatingNetworkError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipatingNetworkError> {
    if cancelled() {
        Err(DissipatingNetworkError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
