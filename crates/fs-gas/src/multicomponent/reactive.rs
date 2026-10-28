//! Atomic coupling of conservative species transport and finite-rate chemistry.
//!
//! Chemistry advances each physical cell at fixed volume, momentum and total
//! chemical energy. A symmetric chemistry/transport/chemistry composition uses
//! the actual accepted transport duration. The chemistry solver's implicit
//! first-order base still limits temporal accuracy; this composition alone is
//! not a second-order or entropy-stability certificate.

use super::{
    Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, StepControl,
    StepReport as TransportStepReport, Totals, allocated, positive, sum,
};
use crate::chemistry::{IsochoricKineticsError, IsochoricKineticsSolve, IsochoricReactor};
use fs_qty::{Amount, Energy, Temperature, Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

/// Transport stability and bounded retries of the complete coupled step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepControl {
    /// Controls for the existing conservative transport owner.
    pub transport: StepControl,
    /// Maximum restarts after the first complete coupled attempt, at most 32.
    pub maximum_coupling_retries: usize,
}
impl Default for ReactiveStepControl {
    fn default() -> Self {
        Self {
            transport: StepControl::default(),
            maximum_coupling_retries: 16,
        }
    }
}

/// Work and source changes from the accepted chemistry half steps only.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChemistryReport {
    /// Number of completed physical-cell solves, normally twice the cell count.
    pub cell_solves: usize,
    /// Accepted adaptive chemistry substeps summed over physical cells.
    pub accepted_steps: usize,
    /// Rejected adaptive chemistry substeps summed over physical cells.
    pub rejected_steps: usize,
    /// Newton iterations summed over physical cells.
    pub nonlinear_iterations: usize,
    /// Largest absolute reactor internal-energy residual [J].
    pub maximum_energy_residual_j: f64,
    /// Largest difference between reactor and conservatively reconstructed U [J].
    pub maximum_reconstruction_residual_j: f64,
    /// Sum of chemistry-only physical entropy changes [J/K].
    /// This excludes transport and boundary entropy and is not a global theorem.
    pub entropy_change_j_k: f64,
    /// Species mass change caused by chemistry [kg], in canonical source order.
    pub species_change_kg: [f64; MAX_SPECIES],
    /// Change in the stored energy gauge due to chemistry [J].
    /// Physical total energy includes the opposite species-reference change.
    pub referenced_energy_change_j: f64,
}
impl Default for ChemistryReport {
    fn default() -> Self {
        Self {
            cell_solves: 0,
            accepted_steps: 0,
            rejected_steps: 0,
            nonlinear_iterations: 0,
            maximum_energy_residual_j: 0.0,
            maximum_reconstruction_residual_j: 0.0,
            entropy_change_j_k: 0.0,
            species_change_kg: [0.0; MAX_SPECIES],
            referenced_energy_change_j: 0.0,
        }
    }
}

/// One published reacting-flow step, with actual source/boundary balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepReport {
    /// Exact duration used by both chemistry half steps and transport.
    pub elapsed: Time,
    /// Whole-step restarts; failed attempts never change the published field.
    pub coupling_retries: usize,
    /// Accepted transport step, whose before/after states lie between chemistry.
    pub transport: TransportStepReport,
    /// Accepted chemistry work and changes.
    pub chemistry: ChemistryReport,
    /// Original complete field.
    pub before: Totals,
    /// Final complete field.
    pub after: Totals,
    /// Species/momentum/gauged-energy residual after boundary and source changes.
    pub conservation_residual: Totals,
    /// Total mass change minus transport through the boundaries [kg].
    pub mass_residual_kg: f64,
    /// Physical total-energy change minus transport through boundaries [J].
    /// Evaluated from differences, avoiding subtraction of large absolute totals.
    pub absolute_energy_residual_j: f64,
}

/// Refusal leaves the public flow and its time unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactiveError {
    /// Incompatible source identity, invalid controls or unrepresentable state.
    Invalid(&'static str),
    /// Conservative transport or caloric reconstruction refused.
    Flow(FlowError),
    /// One original-source finite-rate cell solve refused.
    Chemistry {
        /// Zero-based physical cell index.
        cell: usize,
        /// Original chemistry refusal.
        source: IsochoricKineticsError,
    },
    /// Conservative reconstruction exceeded the declared caloric tolerance.
    EnergyConsistency {
        /// Zero-based physical cell index.
        cell: usize,
        /// Reconstructed minus reactor internal energy [J].
        residual_j: f64,
        /// Admitted energy discrepancy [J].
        tolerance_j: f64,
    },
    /// All bounded complete-step attempts failed.
    RetryLimit {
        /// Number of complete-step attempts.
        attempts: usize,
        /// Last concrete refusal.
        last: Box<ReactiveError>,
    },
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for ReactiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting one-dimensional gas: {self:?}")
    }
}
impl std::error::Error for ReactiveError {}
impl From<FlowError> for ReactiveError {
    fn from(value: FlowError) -> Self {
        match value {
            FlowError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReactiveError> {
    if cancelled() {
        Err(ReactiveError::Cancelled)
    } else {
        Ok(())
    }
}
fn chemistry_error(cell: usize, source: IsochoricKineticsError) -> ReactiveError {
    match source {
        IsochoricKineticsError::Cancelled => ReactiveError::Cancelled,
        source => ReactiveError::Chemistry { cell, source },
    }
}

/// A conservative multicomponent flow coupled to its exact finite-rate mechanism.
///
/// Construction allocates two complete flow workspaces and one reusable reactor.
/// Numerical steps reuse them: the public field is swapped with a fully admitted
/// candidate only after both chemistry half steps and transport finish. Failed
/// attempts and cancellation cannot publish partially reacted cells.
pub struct ReactingEuler1d {
    accepted: MixtureEuler1d,
    trial: MixtureEuler1d,
    reactor: IsochoricReactor,
    amounts: Vec<Amount>,
}
impl ReactingEuler1d {
    /// Bind a complete mechanism to the flow's identical original species bank.
    /// The chemistry interval must fit inside the flow's admitted caloric path.
    pub fn new(
        flow: MixtureEuler1d,
        mechanism: IdealGasKineticsV1,
        chemistry: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(&mut cancelled)?;
        if flow.gas().source() != mechanism.thermodynamics() {
            return Err(ReactiveError::Invalid(
                "flow and mechanism have different original species banks",
            ));
        }
        let domain = flow.gas().domain();
        if chemistry.minimum_temperature_k < domain.minimum_temperature
            || chemistry.maximum_temperature_k > domain.maximum_temperature
        {
            return Err(ReactiveError::Invalid(
                "chemistry interval exceeds the flow caloric path",
            ));
        }
        for i in 0..flow.cell_count() {
            let p = flow
                .gas()
                .primitive(flow.cell(i).expect("admitted cell index"), &mut cancelled)?;
            if !(chemistry.minimum_temperature_k..=chemistry.maximum_temperature_k)
                .contains(&p.temperature)
            {
                return Err(ReactiveError::Invalid(
                    "initial flow temperature is outside chemistry interval",
                ));
            }
        }
        let trial = flow.duplicate(&mut cancelled)?;
        let amounts = allocated(flow.gas().species_count(), Amount::new(0.0))?;
        let reactor = IsochoricReactor::new(mechanism, chemistry, &mut cancelled)
            .map_err(|source| chemistry_error(0, source))?;
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: flow,
            trial,
            reactor,
            amounts,
        })
    }
    /// Read only the accepted flow, including its original sources and mesh.
    pub const fn flow(&self) -> &MixtureEuler1d {
        &self.accepted
    }
    /// Original source-bound finite-rate mechanism.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }
    /// Declared chemistry tolerances and work bounds.
    pub fn chemistry_controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls()
    }
    /// Release the complete accepted conservative field.
    pub fn into_flow(self) -> MixtureEuler1d {
        self.accepted
    }

    /// Advance chemistry/transport/chemistry for one common accepted duration.
    ///
    /// Heat release can increase acoustic speed. If transport reduces the trial
    /// duration, the complete composition restarts from the original state using
    /// that shorter duration; an already integrated chemistry half step is never
    /// silently paired with a different transport duration.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(ReactiveError::Invalid(
                "reacting duration, CFL or retry controls",
            ));
        }
        let before = self.accepted.totals();
        let mut duration = maximum_duration.value();
        // Start within the original acoustic scale before paying for chemistry.
        for i in 0..self.accepted.cell_count() {
            let p = self.accepted.gas().primitive(
                self.accepted.cell(i).expect("admitted index"),
                &mut cancelled,
            )?;
            let area = 0.5 * self.accepted.mesh.face_areas()[i]
                + 0.5 * self.accepted.mesh.face_areas()[i + 1];
            let rate =
                (p.velocity.abs() + p.sound_speed) * (area / self.accepted.mesh.volumes()[i]);
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
        let mut last = ReactiveError::Invalid("no representable reacting time increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "reacting time increment cannot advance the clock",
                ));
            }
            self.trial.copy_physical_from(&self.accepted)?;
            let mut report = ChemistryReport::default();
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut report,
                    &mut cancelled,
                )?;
                let transport =
                    self.trial
                        .step(Time::new(duration), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != duration {
                    return Ok((transport, false));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut report,
                    &mut cancelled,
                )?;
                Ok::<_, ReactiveError>((transport, true))
            })();
            match attempt {
                Ok((transport, false)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveError::Invalid("transport shortened the coupled duration");
                }
                Ok((transport, true)) => {
                    let after = self.trial.totals();
                    let n = self.trial.gas.species_count();
                    let mut residual = [0.0; MAX_COMPONENTS];
                    for k in 0..n {
                        residual[k] = after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                            - report.species_change_kg[k];
                    }
                    residual[n] = after.momentum
                        - before.momentum
                        - transport.boundary_change.momentum
                        - transport.area_momentum_impulse;
                    residual[n + 1] = after.energy
                        - before.energy
                        - transport.boundary_change.energy
                        - report.referenced_energy_change_j;
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
                            * self.trial.gas.energy_references()[k]
                    })));
                    if !residual[..n + 2].iter().all(|x| x.is_finite())
                        || !mass_residual_kg.is_finite()
                        || !absolute_energy_residual_j.is_finite()
                        || !report.entropy_change_j_k.is_finite()
                    {
                        return Err(ReactiveError::Invalid(
                            "non-finite coupled extensive balance",
                        ));
                    }
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    return Ok(ReactiveStepReport {
                        elapsed: Time::new(duration),
                        coupling_retries: retry,
                        transport,
                        chemistry: report,
                        before,
                        after,
                        conservation_residual: Totals::from_values(&residual[..n + 2], n),
                        mass_residual_kg,
                        absolute_energy_residual_j,
                    });
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error) => {
                    last = error;
                    duration *= 0.5;
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn react_cells(
    flow: &mut MixtureEuler1d,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    let n = flow.gas.species_count();
    let width = n + 2;
    let controls = reactor.controls();
    let domain = flow.gas.domain();
    let mut next = [0.0; MAX_COMPONENTS];
    for i in 0..flow.cell_count() {
        poll(cancelled)?;
        let q = &flow.states[i * width..(i + 1) * width];
        let primitive = flow.gas.decode(q, None, cancelled)?;
        let volume = flow.mesh.volumes()[i];
        let thermal = q[n + 1] - 0.5 * q[n] * primitive.velocity;
        let reference = sum(q[..n]
            .iter()
            .zip(flow.gas.energy_references())
            .map(|(rho, e)| rho * e));
        let internal = (thermal + reference) * volume;
        let tolerance = controls.energy_tolerance_j
            + volume
                * (primitive.density * domain.absolute_energy_tolerance
                    + domain.relative_energy_tolerance * thermal.abs());
        let initial_roundtrip = (internal / volume - reference - thermal) * volume;
        if !internal.is_finite() || !tolerance.is_finite() || initial_roundtrip.abs() > tolerance {
            return Err(ReactiveError::EnergyConsistency {
                cell: i,
                residual_j: initial_roundtrip,
                tolerance_j: tolerance,
            });
        }
        for k in 0..n {
            amounts[k] =
                Amount::new((q[k] * volume) / flow.gas.source().models()[k].molar_mass().value());
        }
        let endpoint = reactor
            .advance_adiabatic_energy(
                amounts,
                Energy::new(internal),
                Volume::new(volume),
                Temperature::new(primitive.temperature),
                Time::new(duration),
                &mut *cancelled,
            )
            .map_err(|source| chemistry_error(i, source))?;
        for k in 0..n {
            next[k] = (endpoint.amounts()[k].value()
                * flow.gas.source().models()[k].molar_mass().value())
                / volume;
        }
        next[n] = q[n];
        // A linear conserved-variable gauge: never add a separate heat-release
        // source to the original chemical total energy.
        let reference_change =
            sum((0..n).map(|k| (next[k] - q[k]) * flow.gas.energy_references()[k]));
        next[n + 1] = q[n + 1] - reference_change;
        let resolved = flow.gas.decode(
            &next[..width],
            Some(endpoint.temperature().value()),
            cancelled,
        )?;
        let next_thermal = next[n + 1] - 0.5 * next[n] * resolved.velocity;
        let balance = endpoint.balance();
        let mut source_thermal = [0.0; MAX_SPECIES];
        for k in 0..n {
            source_thermal[k] = (next[k] / flow.gas.source().models()[k].molar_mass().value())
                * flow.gas.source().models()[k]
                    .molar_internal_energy_change(
                        Temperature::new(domain.minimum_temperature),
                        endpoint.temperature(),
                    )
                    .map_err(FlowError::from)?
                    .value();
        }
        // Independently reconstruct sensible energy at the reactor's endpoint
        // temperature. Checking only the gauge update would be an algebraic
        // identity and could hide a thermodynamically inconsistent endpoint.
        let reconstruction_residual =
            (next_thermal - sum(source_thermal[..n].iter().copied())) * volume;
        let physical_residual = ((next_thermal - thermal) + reference_change) * volume;
        let reconstruction_tolerance = tolerance + balance.energy_residual_j.abs();
        if !reconstruction_residual.is_finite()
            || reconstruction_residual.abs() > reconstruction_tolerance
        {
            return Err(ReactiveError::EnergyConsistency {
                cell: i,
                residual_j: reconstruction_residual,
                tolerance_j: reconstruction_tolerance,
            });
        }
        // A molecular-mass mismatch may also change kinetic energy at fixed
        // momentum. It must fit the same actual energy budget, not be hidden.
        if !physical_residual.is_finite() || physical_residual.abs() > tolerance {
            return Err(ReactiveError::EnergyConsistency {
                cell: i,
                residual_j: physical_residual,
                tolerance_j: tolerance,
            });
        }
        report.cell_solves += 1;
        report.accepted_steps += balance.accepted_steps;
        report.rejected_steps += balance.rejected_steps;
        report.nonlinear_iterations += balance.nonlinear_iterations;
        report.maximum_energy_residual_j = report
            .maximum_energy_residual_j
            .max(balance.energy_residual_j.abs());
        report.maximum_reconstruction_residual_j = report
            .maximum_reconstruction_residual_j
            .max(reconstruction_residual.abs());
        report.entropy_change_j_k += balance.entropy_change_j_k;
        for k in 0..n {
            report.species_change_kg[k] += (next[k] - q[k]) * volume;
        }
        report.referenced_energy_change_j += (next[n + 1] - q[n + 1]) * volume;
        flow.states[i * width..(i + 1) * width].copy_from_slice(&next[..width]);
    }
    flow.validate_states(cancelled)?;
    poll(cancelled)
}

#[cfg(test)]
mod tests;
