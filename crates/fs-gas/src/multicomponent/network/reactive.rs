//! Finite-rate chemistry in every shock-carrying branch and compact plenum.
//!
//! A complete C(dt/2), network(dt), C(dt/2) composition publishes atomically.
//! Chemistry reuses the original single-cell fixed-U,V update at each actual
//! physical volume. Changing species references changes stored sensible energy,
//! never adds a second heat-release source to chemical total energy. This is
//! an inviscid estimated composition, not a global entropy or order certificate.

use super::super::evolve::totals;
pub use super::super::reactive::{ChemistryReport, ReactiveStepControl};
use super::super::reactive::{ReactiveError, react_cell};
use super::super::{
    Cell, Conserved, Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, MixtureGas,
    StepReport, allocated, positive, sum,
};
use super::{
    Accumulator, BranchEnd, MixtureNetwork1d, NetworkStepReport, NetworkTotals, Plenum,
    PlenumConnection, PortExchange,
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
    /// Original complete branch and plenum inventories.
    pub before: NetworkTotals,
    /// Final complete branch and plenum inventories after the second chemistry half.
    pub after: NetworkTotals,
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
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: network,
            trial,
            reactor,
            amounts,
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
        let mut duration = self.initial_duration(
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
                    return Ok((transport, false));
                }
                react_network(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                Ok::<_, ReactiveNetworkError>((transport, true))
            })();
            match attempt {
                Ok((transport, false)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveNetworkError::Invalid(
                        "network transport shortened the coupled duration",
                    );
                }
                Ok((transport, true)) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut residual = NetworkTotals {
                        species_masses: [0.0; MAX_SPECIES],
                        species_count: n,
                        energy: after.energy
                            - before.energy
                            - transport.boundary_change.energy
                            - chemistry.referenced_energy_change_j,
                    };
                    for k in 0..n {
                        residual.species_masses[k] = after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                            - chemistry.species_change_kg[k];
                    }
                    let mass_residual_kg = sum((0..n).map(|k| {
                        after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                    }));
                    let absolute_energy_residual_j = sum(core::iter::once(
                        after.energy - before.energy - transport.boundary_change.energy,
                    )
                    .chain((0..n).map(|k| {
                        (after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k])
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
                        before,
                        after,
                        conservation_residual: residual,
                        mass_residual_kg,
                        absolute_energy_residual_j,
                    };
                    poll(&mut cancelled)?;
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

    // Bound the first chemical attempt by both original branch acoustics and
    // compact-volume exchange scales. Actual post-reaction port CFL remains the
    // network solver's responsibility and can restart this complete composition.
    fn initial_duration(
        &mut self,
        maximum: f64,
        cfl: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ReactiveNetworkError> {
        let mut rate = 0.0_f64;
        for branch in &self.accepted.branches {
            let width = branch.gas.species_count() + 2;
            for (i, q) in branch.states.chunks_exact(width).enumerate() {
                let p = branch.gas.decode(q, None, cancelled)?;
                let area =
                    0.5 * branch.mesh.face_areas()[i] + 0.5 * branch.mesh.face_areas()[i + 1];
                let local = (p.velocity.abs() + p.sound_speed) * (area / branch.mesh.volumes()[i]);
                if !positive(local) {
                    return Err(ReactiveNetworkError::Invalid(
                        "unrepresentable reacting branch CFL",
                    ));
                }
                rate = rate.max(local);
            }
        }
        self.trial.work.plenum_signal_rates.fill(0.0);
        for (i, volume) in self.accepted.plenums.iter().enumerate() {
            self.trial.work.plenum_primitive[i] = self
                .accepted
                .gas()
                .primitive(volume.state(), &mut *cancelled)?;
        }
        for connection in &self.accepted.connections {
            poll(cancelled)?;
            let branch = &self.accepted.branches[connection.branch];
            let face = match connection.end {
                BranchEnd::Left => 0,
                BranchEnd::Right => branch.cell_count(),
            };
            let local = 0.5
                * branch.mesh.face_areas()[face]
                * (self.trial.work.plenum_primitive[connection.plenum].sound_speed
                    / self.accepted.plenums[connection.plenum].volume_m3);
            let total = self.trial.work.plenum_signal_rates[connection.plenum] + local;
            if !positive(local) || !positive(total) {
                return Err(ReactiveNetworkError::Invalid(
                    "unrepresentable reacting plenum CFL",
                ));
            }
            self.trial.work.plenum_signal_rates[connection.plenum] = total;
            rate = rate.max(total);
        }
        Ok(maximum.min(0.8 * cfl / rate))
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

fn checked_totals(
    network: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NetworkTotals, ReactiveNetworkError> {
    let n = network.gas().species_count();
    let mut accumulated = Accumulator::new(n);
    for branch in &network.branches {
        let value = totals(&branch.states, branch.mesh.volumes(), n, cancelled)?;
        accumulated.add(&value.species_masses[..n], value.energy)?;
    }
    for volume in &network.plenums {
        poll(cancelled)?;
        accumulated.add_cell(volume.state(), volume.volume_m3)?;
    }
    Ok(accumulated.totals())
}

fn duplicate(
    network: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MixtureNetwork1d, ReactiveNetworkError> {
    let mut branches = Vec::new();
    branches
        .try_reserve_exact(network.branches.len())
        .map_err(|_| FlowError::Allocation)?;
    for branch in &network.branches {
        branches.push(branch.duplicate(cancelled)?);
    }
    let mut plenums = Vec::new();
    plenums
        .try_reserve_exact(network.plenums.len())
        .map_err(|_| FlowError::Allocation)?;
    for plenum in &network.plenums {
        poll(cancelled)?;
        let mut densities = allocated(plenum.state.partial_densities.len(), 0.0)?;
        densities.copy_from_slice(&plenum.state.partial_densities);
        plenums.push(Plenum::new(
            plenum.volume_m3,
            Conserved {
                partial_densities: densities,
                momentum: plenum.state.momentum,
                energy: plenum.state.energy,
            },
        )?);
    }
    let mut connections = Vec::new();
    connections
        .try_reserve_exact(network.connections.len())
        .map_err(|_| FlowError::Allocation)?;
    connections.extend_from_slice(&network.connections);
    Ok(MixtureNetwork1d::new(
        branches,
        plenums,
        connections,
        cancelled,
    )?)
}

// These networks were constructed from identical immutable owners and topology.
// Only their fields and accepted clocks can differ. Reset numerical candidates
// in bounded copy tiles without allocations or repeated source-bank comparisons.
fn reset(
    trial: &mut MixtureNetwork1d,
    accepted: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveNetworkError> {
    for (trial, accepted) in trial.branches.iter_mut().zip(&accepted.branches) {
        for (destination, source) in trial
            .states
            .chunks_mut(4096)
            .zip(accepted.states.chunks(4096))
        {
            poll(cancelled)?;
            destination.copy_from_slice(source);
        }
        trial.time = accepted.time;
    }
    for (trial, accepted) in trial.plenums.iter_mut().zip(&accepted.plenums) {
        poll(cancelled)?;
        trial
            .state
            .partial_densities
            .copy_from_slice(&accepted.state.partial_densities);
        trial.state.momentum = accepted.state.momentum;
        trial.state.energy = accepted.state.energy;
    }
    trial.last_report = None;
    Ok(())
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
            ReactiveError::Invalid("initial network temperature is outside chemistry interval"),
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
