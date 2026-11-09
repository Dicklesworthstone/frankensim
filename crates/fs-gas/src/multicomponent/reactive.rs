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
    Cell, Error as FlowError, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, StepControl,
    StepReport as TransportStepReport, Totals, allocated, positive, sum,
};
use crate::chemistry::{IsochoricKineticsError, IsochoricReactor};
use crate::body_force::{BodyForce, ForceChange, ForceControl};
pub use crate::chemistry::{IsochoricKineticsMethod, IsochoricKineticsSolve};
use fs_qty::{Amount, Time};
pub use fs_thermochem::kinetics::IdealGasKineticsV1;

mod cell;
pub use cell::CellChemistry;
mod entropy;
pub use entropy::{
    ReactiveEntropyReport, ReactiveEntropyStepReport, ReactiveExteriorEntropyReport,
    ReactiveExteriorEntropyStepReport, SplitExteriorMolecularEntropyReport,
    StabilizedReactiveExteriorEntropyStepReport, StabilizedReactiveExteriorEntropyIntervalReport,
    SplitMolecularEntropyReport,
};
mod reservoir;
mod stabilized;
pub use stabilized::{
    SplitStabilizedDissipationReport, StabilizedReactiveIntervalControl,
    StabilizedReactiveEntropyIntervalReport, StabilizedReactiveEntropyStepReport,
    StabilizedReactiveIntervalReport,
    StabilizedReactiveStepReport,
};
mod inertial_walls;
pub use inertial_walls::ReactiveInertialWallsStep;
mod piston;
pub use piston::{
    ReactiveSharedPistonGasReport, ReactiveSharedPistonStep, ReactiveSharedPistonThermalStep,
};
pub use piston::{
    PistonThermalEnds, PistonThermalReport, ReactivePistonStep, SplitPistonThermalReport,
    ThermalReactiveStepReport,
};

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

impl ReactiveStepControl {
    /// Admit the original acoustic and whole-composition retry controls.
    pub fn validate(self) -> Result<(), ReactiveError> {
        if self.maximum_coupling_retries > 32
            || !positive(self.transport.cfl)
            || self.transport.cfl > 0.45
            || self.transport.maximum_halvings > 32
        {
            return Err(ReactiveError::Invalid("reacting duration, CFL or retry controls"));
        }
        Ok(())
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
    /// Physical total-energy change minus boundary transport and body work [J].
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
    /// The original resolved pressure/piston owner refused a coupled trial.
    Piston(crate::dynamics::piston::PistonError),
    /// Original end-wall heat storage or thermal exchange refused.
    Wall(super::wall_exchange::WallExchangeError),
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
    /// A complete-duration update exhausted its coupled-step budget.
    IntervalWorkLimit {
        elapsed_time_s: f64,
        completed_steps: usize,
        attempts: usize,
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
impl From<crate::dynamics::piston::PistonError> for ReactiveError {
    fn from(value: crate::dynamics::piston::PistonError) -> Self {
        match value {
            crate::dynamics::piston::PistonError::Cancelled => Self::Cancelled,
            other => Self::Piston(other),
        }
    }
}
impl From<super::wall_exchange::WallExchangeError> for ReactiveError {
    fn from(value: super::wall_exchange::WallExchangeError) -> Self {
        match value {
            super::wall_exchange::WallExchangeError::Cancelled => Self::Cancelled,
            other => Self::Wall(other),
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
    // Only complete-duration evolution lazily reserves this physical rollback.
    interval_initial: Vec<f64>,
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
            interval_initial: Vec::new(),
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
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            [None; 2],
            None,
            cancelled,
        ).map(|(report, _)| report)
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
            stabilized::MolecularMethod::Explicit(dissipation_controls),
            [None; 2],
            None,
            cancelled,
        ).map(|(report, _)| report)
    }

    /// Advance finite-rate chemistry, body acceleration and conservative flow.
    ///
    /// The existing exact mechanical kicks enclose the original transport step
    /// inside the chemistry half steps: C(h/2), F(h/2), A(h), F(h/2), C(h/2).
    /// Acceleration uses the actual mixture density at each force stage; force
    /// density retains its authored units. Mechanical work changes kinetic and
    /// total energy together, while chemistry retains the original species and
    /// formation-energy references. The accepted impulse/work are reported in
    /// `transport.source_change` and enter the complete balance exactly once.
    ///
    /// Heat release and acceleration share the original acoustic retry owner.
    /// A shortened transport step restarts the entire composition, including
    /// chemistry and both kicks, from the accepted field. No partial chemistry,
    /// force, clock or report is published. A zero source takes the unforced
    /// numerical path. Forcing is held over the accepted interval; this method
    /// does not infer a time-dependent law, moving walls or entropy authority.
    pub fn step_with_body_force(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        self.step_inner_forced::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            [None; 2],
            None,
            None,
            Some((force, force_controls)),
            cancelled,
        ).map(|(report, _)| report)
    }

    /// Compose the original molecular half steps with reacting forced flow.
    ///
    /// D(h/2) and D(h/2) surround the complete chemistry/force/transport map in
    /// [`Self::step_with_body_force`]. Molecular stress, species diffusion and
    /// heat retain the supplied source-bound transport law. Both source kicks
    /// and both molecular halves use the same finally accepted duration; final
    /// caloric and Mach admission precede atomic publication. The original
    /// chemistry and molecular time methods retain their own accuracy limits.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_body_force<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        self.step_inner_forced(
            maximum_duration,
            controls,
            Some(operator),
            stabilized::MolecularMethod::Explicit(dissipation_controls),
            [None; 2],
            None,
            None,
            Some((force, force_controls)),
            cancelled,
        ).map(|(report, _)| report)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_impl<L: TransportLaw>(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        operator: Option<&mut DissipationOperator<L>>,
        molecular: stabilized::MolecularMethod<'_>,
        reservoirs: [Option<Cell<'_>>; 2],
        interval_budget: Option<&mut stabilized::CouplingBudget>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(ReactiveStepReport, Option<SplitStabilizedDissipationReport>), ReactiveError> {
        self.step_inner(
            maximum_duration,
            controls,
            operator,
            molecular,
            reservoirs,
            interval_budget,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_inner<L: TransportLaw>(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        operator: Option<&mut DissipationOperator<L>>,
        molecular: stabilized::MolecularMethod<'_>,
        reservoirs: [Option<Cell<'_>>; 2],
        interval_budget: Option<&mut stabilized::CouplingBudget>,
        entropy: Option<&mut entropy::CoupledEntropy>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(ReactiveStepReport, Option<SplitStabilizedDissipationReport>), ReactiveError> {
        self.step_inner_forced(
            maximum_duration, controls, operator, molecular, reservoirs,
            interval_budget, entropy, None, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_inner_forced<L: TransportLaw>(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut operator: Option<&mut DissipationOperator<L>>,
        mut molecular: stabilized::MolecularMethod<'_>,
        reservoirs: [Option<Cell<'_>>; 2],
        mut interval_budget: Option<&mut stabilized::CouplingBudget>,
        mut entropy: Option<&mut entropy::CoupledEntropy>,
        mut forcing: Option<(BodyForce<'_>, ForceControl)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(ReactiveStepReport, Option<SplitStabilizedDissipationReport>), ReactiveError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(maximum_duration.value()) {
            return Err(ReactiveError::Invalid(
                "reacting duration, CFL or retry controls",
            ));
        }
        if let Some((force, force_controls)) = forcing {
            force_controls.admit::<FlowError>()?;
            if !force.admit::<FlowError>(self.accepted.cell_count(), &mut cancelled)? {
                forcing = None;
            }
        }
        if forcing.is_some()
            && (entropy.is_some() || reservoirs.iter().any(Option::is_some))
        {
            return Err(ReactiveError::Invalid(
                "reacting body forcing does not compose prescribed reservoirs or entropy accounting",
            ));
        }
        if entropy.as_ref().is_some_and(|account| account.has_exterior_states()) {
            // Only the dedicated entropy entry interprets these borrowed cells
            // as moving exterior traces. Legacy reservoir admission is unchanged.
            if operator.is_some() {
                // Molecular ports retain the existing stationary bath model.
                // The same source-compatible states feed both D halves and A.
                reservoir::admit(&self.accepted, reservoirs, &mut cancelled)?;
            }
        } else {
            reservoir::admit(&self.accepted, reservoirs, &mut cancelled)?;
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
            if let Some(budget) = interval_budget.as_deref_mut() {
                budget.reserve(self.accepted.time)?;
            }
            if !positive(duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "reacting time increment cannot advance the clock",
                ));
            }
            self.trial.copy_physical_with_motion_from(&self.accepted, &mut cancelled)?;
            if let Some(account) = entropy.as_deref_mut() {
                account.begin_attempt();
            }
            let mut report = ChemistryReport::default();
            let attempt = (|| {
                let first = match operator.as_deref_mut() {
                    Some(operator) => Some(molecular.advance(
                        operator,
                        &mut self.trial,
                        Time::new(0.5 * duration),
                        reservoirs,
                        entropy.as_deref_mut(),
                        &mut cancelled,
                    )?),
                    None => None,
                };
                react_cells_inner(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut report,
                    entropy.as_deref_mut(),
                    &mut cancelled,
                )?;
                let (transport, transport_entropy) = if let Some((force, force_controls)) = forcing {
                    (
                        self.trial.step_with_body_force(
                            Time::new(duration),
                            controls.transport,
                            force,
                            force_controls,
                            &mut cancelled,
                        )?,
                        None,
                    )
                } else if let Some(account) = entropy.as_deref_mut() {
                    account.advance_transport(
                        &mut self.trial,
                        Time::new(duration),
                        controls.transport,
                        reservoirs,
                        &mut cancelled,
                    )?
                } else {
                    (
                        self.trial.step_with_reservoirs(
                            Time::new(duration),
                            controls.transport,
                            reservoirs,
                            &mut cancelled,
                        )?,
                        None,
                    )
                };
                if transport.elapsed.value() != duration {
                    return Ok((transport, false, None, None));
                }
                react_cells_inner(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut report,
                    entropy.as_deref_mut(),
                    &mut cancelled,
                )?;
                let second = match operator.as_deref_mut() {
                    Some(operator) => Some(molecular.advance(
                        operator,
                        &mut self.trial,
                        Time::new(0.5 * duration),
                        reservoirs,
                        entropy.as_deref_mut(),
                        &mut cancelled,
                    )?),
                    None => None,
                };
                if (second.is_some() || forcing.is_some()) && entropy.is_none() {
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
                        if let Some((_, force_controls)) = forcing {
                            force_controls.admit_velocity::<FlowError>(
                                primitive.velocity, primitive.sound_speed,
                            )?;
                        }
                    }
                }
                if let Some(account) = entropy.as_deref_mut() {
                    let temperature_bounds = second.as_ref().map(|_| {
                        let chemistry = self.reactor.controls();
                        [
                            chemistry.minimum_temperature_k,
                            chemistry.maximum_temperature_k,
                        ]
                    });
                    if account.has_exterior_states() {
                        account.check_exterior_candidate(
                            &self.accepted,
                            &self.trial,
                            temperature_bounds,
                            [
                                first.as_ref().and_then(|half| half.exterior_entropy),
                                second.as_ref().and_then(|half| half.exterior_entropy),
                            ],
                            &mut cancelled,
                        )?;
                    } else {
                        account.check_candidate(
                            &self.accepted,
                            &self.trial,
                            temperature_bounds,
                            transport_entropy.ok_or(ReactiveError::Invalid(
                                "missing reacting advection entropy receipt",
                            ))?,
                            [
                                first.as_ref().and_then(|half| half.entropy),
                                second.as_ref().and_then(|half| half.entropy),
                            ],
                            &mut cancelled,
                        )?;
                    }
                }
                let (dissipation, stabilized) = stabilized::split_halves(first, second);
                Ok::<_, ReactiveError>((transport, true, dissipation, stabilized))
            })();
            match attempt {
                Ok((transport, false, _, _)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveError::Invalid("transport shortened the coupled duration");
                }
                Ok((transport, true, dissipation, stabilized)) => {
                    let after = self.trial.totals();
                    let complete = coupled_report(&self.trial.gas, before, after,
                        transport, report, dissipation, retry)?;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    return Ok((complete, stabilized));
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error @ ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { .. })) => {
                    return Err(error);
                }
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

// Shared accounting for fixed and moving-volume reactive compositions. Actual
// endpoint species changes retain the original physical chemical-energy datum.
fn coupled_report(
    gas: &super::MixtureGas,
    before: Totals,
    after: Totals,
    transport: TransportStepReport,
    chemistry: ChemistryReport,
    dissipation: Option<SplitDissipationReport>,
    coupling_retries: usize,
) -> Result<ReactiveStepReport, ReactiveError> {
    let n = gas.species_count();
    let mut molecular_boundary = [0.0; MAX_COMPONENTS];
    if let Some(d) = dissipation {
        for k in 0..n + 2 {
            molecular_boundary[k] =
                d.first.boundary_change.get(k) + d.second.boundary_change.get(k);
        }
    }
    let (conservation_residual, mass_residual_kg, absolute_energy_residual_j) = coupled_balance(
        gas, before, after, transport.boundary_change,
        Totals::from_values(&molecular_boundary[..n + 2], n),
        transport.area_momentum_impulse, transport.source_change, &chemistry,
    )?;
    Ok(ReactiveStepReport {
        elapsed: transport.elapsed,
        coupling_retries,
        transport,
        chemistry,
        dissipation,
        before,
        after,
        conservation_residual,
        mass_residual_kg,
        absolute_energy_residual_j,
    })
}


// One endpoint balance owner for accepted steps and complete stationary windows.
#[allow(clippy::too_many_arguments)]
fn coupled_balance(
    gas: &super::MixtureGas,
    before: Totals,
    after: Totals,
    transport_boundary: Totals,
    molecular_boundary: Totals,
    area_momentum_impulse: f64,
    source_change: ForceChange,
    chemistry: &ChemistryReport,
) -> Result<(Totals, f64, f64), ReactiveError> {
    let n = gas.species_count();
    let mut residual = [0.0; MAX_COMPONENTS];
    for k in 0..n {
        residual[k] = after.species_masses[k]
            - before.species_masses[k]
            - transport_boundary.species_masses[k]
            - molecular_boundary.get(k)
            - chemistry.species_change_kg[k];
    }
    residual[n] = after.momentum
        - before.momentum
        - transport_boundary.momentum
        - molecular_boundary.get(n)
        - area_momentum_impulse
        - source_change.impulse;
    residual[n + 1] = after.energy
        - before.energy
        - transport_boundary.energy
        - molecular_boundary.get(n + 1)
        - chemistry.referenced_energy_change_j
        - source_change.work;
    let mass_residual_kg = sum((0..n).map(|k| {
        after.species_masses[k]
            - before.species_masses[k]
            - transport_boundary.species_masses[k]
            - molecular_boundary.get(k)
    }));
    let absolute_energy_residual_j = sum(core::iter::once(
        after.energy - before.energy - transport_boundary.energy - molecular_boundary.get(n + 1)
            - source_change.work,
    )
    .chain((0..n).map(|k| {
        (after.species_masses[k]
            - before.species_masses[k]
            - transport_boundary.species_masses[k]
            - molecular_boundary.get(k))
            * gas.energy_references()[k]
    })));
    if !residual[..n + 2].iter().all(|x| x.is_finite())
        || !mass_residual_kg.is_finite()
        || !absolute_energy_residual_j.is_finite()
        || !chemistry.entropy_change_j_k.is_finite()
    {
        return Err(ReactiveError::Invalid(
            "non-finite coupled extensive balance",
        ));
    }
    Ok((Totals::from_values(&residual[..n + 2], n), mass_residual_kg, absolute_energy_residual_j))
}

fn react_cells(
    flow: &mut MixtureEuler1d,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    react_cells_inner(flow, reactor, amounts, duration, report, None, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn react_cells_inner(
    flow: &mut MixtureEuler1d,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    mut entropy: Option<&mut entropy::CoupledEntropy>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    let width = flow.gas.species_count() + 2;
    let mut previous = [0.0; MAX_COMPONENTS];
    for i in 0..flow.cell_count() {
        if entropy.is_some() {
            previous[..width].copy_from_slice(&flow.states[i * width..(i + 1) * width]);
        }
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
        if let Some(account) = entropy.as_deref_mut() {
            account.observe_chemistry_row(
                &previous[..width],
                &flow.states[i * width..(i + 1) * width],
                flow.mesh.volumes()[i],
                cancelled,
            )?;
        }
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
