//! Valve/wave composition using the same reusable source as reacting owners.
use super::*;

/// Complete accepted valve/wave composition. Valve exchange is internal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValvedNetworkReport {
    /// Duration shared by transport and both valve half sweeps together.
    pub elapsed: Time,
    /// Complete discarded compositions before this accepted step.
    pub coupling_retries: usize,
    /// Intermediate branch/plenum advection balance, not the complete endpoint.
    pub transport: NetworkStepReport,
    /// The accepted forward and reverse valve half sweeps.
    pub valves: SplitNetworkValveReport,
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
    valves: NetworkValveOperator,
    last_report: Option<ValvedNetworkReport>,
}
impl ValvedNetwork1d {
    /// Bind explicit inter-plenum topology using the original source operator.
    pub fn new(
        network: MixtureNetwork1d,
        links: Vec<ValveLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ValveError> {
        let valves = NetworkValveOperator::new(&network, links, &mut cancelled)?;
        let trial = duplicate(&network, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self { accepted: network, trial, valves, last_report: None })
    }
    /// Exact original common source bank and energy references.
    pub fn gas(&self) -> &MixtureGas { self.accepted.gas() }
    /// Accepted original branches; no independent mutable handle is exposed.
    pub fn branches(&self) -> &[MixtureEuler1d] { self.accepted.branches() }
    /// Accepted original compact volumes, including valve-connected plenums.
    pub fn plenums(&self) -> &[Plenum] { self.accepted.plenums() }
    /// Retained deterministic valve topology and current settings.
    pub fn links(&self) -> &[ValveLink] { self.valves.links() }
    /// Common accepted physical clock.
    pub fn time(&self) -> Time { self.accepted.time() }
    /// Complete accepted species and energy storage.
    pub fn totals(&self) -> NetworkTotals { self.accepted.totals() }
    /// Last complete accepted report, unchanged after a refused trial.
    pub const fn last_report(&self) -> Option<ValvedNetworkReport> { self.last_report }
    /// Integrated inter-plenum transfers from both accepted half sweeps.
    pub fn accepted_valve_exchanges(&self) -> Option<&[ValveExchange]> { self.valves.accepted_exchanges() }
    /// Intermediate advective branch/plenum exchanges, excluding valve links.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.accepted.port_exchanges.as_slice())
    }
    /// Intermediate branch advection balances, before the final valve half.
    pub fn accepted_branch_transport_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.accepted.branch_reports.as_slice())
    }
    /// Change an opening exactly between accepted steps without changing fields.
    pub fn set_orifice(&mut self, index: usize, orifice: Orifice) -> Result<(), ValveError> {
        self.valves.set_orifice(index, orifice)
    }
    /// Release the physical network, clearing its intermediate advection report.
    pub fn into_network(mut self) -> MixtureNetwork1d {
        self.accepted.last_report = None;
        self.accepted
    }
    /// Advance V(dt/2), A(dt), V(dt/2). A shorter acoustic step restarts both
    /// valve halves. The final gas and valve receipts publish together.
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
            self.valves.begin_trial(&mut cancelled)?;
            let attempt = (|| {
                let first = self.valves.advance_for_trial(&mut self.trial,
                    Time::new(0.5 * dt), controls.into(), false, &mut cancelled)?;
                let transport = self.trial.step(Time::new(dt), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != dt { return Ok((transport, None)); }
                let second = self.valves.advance_for_trial(&mut self.trial,
                    Time::new(0.5 * dt), controls.into(), true, &mut cancelled)?;
                Ok::<_, ValveError>((transport, Some(SplitNetworkValveReport { first, second })))
            })();
            match attempt {
                Ok((transport, None)) => {
                    dt = transport.elapsed.value();
                    last = ValveError::Invalid("advection shortened valve composition");
                }
                Ok((transport, Some(valves))) => {
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
                        transport, valves, before, after, conservation_residual: residual,
                        absolute_energy_residual_j: physical };
                    poll(&mut cancelled)?;
                    self.valves.publish_trial_report();
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(ValveError::Cancelled) => return Err(ValveError::Cancelled),
                Err(error) => { last = error; dt *= 0.5; }
            }
        }
        Err(ValveError::RetryLimit { attempts: controls.maximum_coupling_retries + 1, last: Box::new(last) })
    }
}
