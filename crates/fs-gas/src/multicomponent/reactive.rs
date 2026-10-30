//! Atomic coupling of conservative species transport and finite-rate chemistry.
//!
//! Chemistry advances each physical cell at fixed volume, momentum and total
//! chemical energy. A symmetric chemistry/transport/chemistry composition uses
//! the actual accepted transport duration. Optional molecular half steps enclose
//! that composition, and their final state must retain the chemistry interval.
//! The original chemistry owner retains its selected time method; symmetric
//! composition alone is not a second-order or entropy-stability certificate.

use super::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, DissipationOperator,
    DissipationReport, TransportLaw,
};
use super::{
    Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, StepControl,
    StepReport as TransportStepReport, Totals, allocated, positive, sum,
};
use crate::chemistry::{IsochoricKineticsError, IsochoricReactor};
pub use crate::chemistry::{IsochoricKineticsMethod, IsochoricKineticsSolve};
use fs_qty::{Amount, Time};
pub use fs_thermochem::kinetics::IdealGasKineticsV1;

mod cell;
pub use cell::CellChemistry;

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

/// Molecular transport on both sides of the chemistry/flow composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitDissipationReport {
    /// Complete molecular half step before chemistry and advection.
    pub first: DissipationReport,
    /// Complete molecular half step after chemistry and advection.
    pub second: DissipationReport,
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
    /// Molecular species, heat and momentum transport when explicitly enabled.
    pub dissipation: Option<SplitDissipationReport>,
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
    /// Molecular transport refused before the complete step could publish.
    Dissipation(DissipationError),
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
        write!(f, "reacting multicomponent gas: {self:?}")
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
impl From<DissipationError> for ReactiveError {
    fn from(value: DissipationError) -> Self {
        match value {
            DissipationError::Cancelled => Self::Cancelled,
            other => Self::Dissipation(other),
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            cancelled,
        )
    }

    /// Advance molecular diffusion, viscosity and heat conduction together with
    /// conservative gas transport and finite-rate chemistry.
    ///
    /// The symmetric composition is D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2).
    /// Each molecular half step may take bounded parabolic substeps but must
    /// complete its exact duration. A shorter acoustic step restarts all five
    /// stages. Only the accepted advection step advances the physical clock.
    pub fn step_with_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
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
        mut operator: Option<&mut DissipationOperator<L>>,
        dissipation_controls: DissipationControl,
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
                let first = match operator.as_deref_mut() {
                    Some(operator) => Some(operator.advance_for(
                        &mut self.trial,
                        Time::new(0.5 * duration),
                        dissipation_controls,
                        &mut cancelled,
                    )?),
                    None => None,
                };
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
                    return Ok((transport, false, None));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut report,
                    &mut cancelled,
                )?;
                let second = match operator.as_deref_mut() {
                    Some(operator) => Some(operator.advance_for(
                        &mut self.trial,
                        Time::new(0.5 * duration),
                        dissipation_controls,
                        &mut cancelled,
                    )?),
                    None => None,
                };
                if second.is_some() {
                    // The molecular source bank can admit a wider temperature
                    // path than the reactor. In particular, viscous heating
                    // after the last chemistry solve can leave that path.
                    // Reject this whole attempt before publishing such a state.
                    let chemistry = self.reactor.controls();
                    for i in 0..self.trial.cell_count() {
                        poll(&mut cancelled)?;
                        let primitive = self.trial.gas().primitive(
                            self.trial.cell(i).expect("admitted cell index"),
                            &mut cancelled,
                        )?;
                        if !(chemistry.minimum_temperature_k..=chemistry.maximum_temperature_k)
                            .contains(&primitive.temperature)
                        {
                            return Err(ReactiveError::Invalid(
                                "final molecular state is outside chemistry interval",
                            ));
                        }
                    }
                }
                let dissipation = first
                    .zip(second)
                    .map(|(first, second)| SplitDissipationReport { first, second });
                Ok::<_, ReactiveError>((transport, true, dissipation))
            })();
            match attempt {
                Ok((transport, false, _)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveError::Invalid("transport shortened the coupled duration");
                }
                Ok((transport, true, dissipation)) => {
                    let after = self.trial.totals();
                    let n = self.trial.gas.species_count();
                    let mut residual = [0.0; MAX_COMPONENTS];
                    let mut molecular_boundary = [0.0; MAX_COMPONENTS];
                    if let Some(d) = dissipation {
                        for k in 0..n + 2 {
                            molecular_boundary[k] =
                                d.first.boundary_change.get(k) + d.second.boundary_change.get(k);
                        }
                    }
                    for k in 0..n {
                        residual[k] = after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                            - molecular_boundary[k]
                            - report.species_change_kg[k];
                    }
                    residual[n] = after.momentum
                        - before.momentum
                        - transport.boundary_change.momentum
                        - molecular_boundary[n]
                        - transport.area_momentum_impulse;
                    residual[n + 1] = after.energy
                        - before.energy
                        - transport.boundary_change.energy
                        - molecular_boundary[n + 1]
                        - report.referenced_energy_change_j;
                    let mass_residual_kg = sum((0..n).map(|k| {
                        after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                            - molecular_boundary[k]
                    }));
                    let absolute_energy_residual_j = sum(core::iter::once(
                        after.energy
                            - before.energy
                            - transport.boundary_change.energy
                            - molecular_boundary[n + 1],
                    )
                    .chain((0..n).map(|k| {
                        (after.species_masses[k]
                            - before.species_masses[k]
                            - transport.boundary_change.species_masses[k]
                            - molecular_boundary[k])
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
                        dissipation,
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
    let width = flow.gas.species_count() + 2;
    for i in 0..flow.cell_count() {
        react_cell(
            &flow.gas,
            &mut flow.states[i * width..(i + 1) * width],
            flow.mesh.volumes()[i],
            i,
            reactor,
            amounts,
            duration,
            report,
            cancelled,
        )?;
    }
    flow.validate_states(cancelled)?;
    poll(cancelled)
}

/// One physical control volume at fixed momentum and chemical total energy.
/// Branch and compact-volume owners reuse the identical original-source update.
#[allow(clippy::too_many_arguments)]
pub(super) fn react_cell(
    gas: &super::MixtureGas,
    q: &mut [f64],
    volume: f64,
    i: usize,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    cell::react_components::<1>(
        gas, q, volume, i, reactor, amounts, duration, report, cancelled,
    )
}

#[cfg(test)]
mod tests;
