use super::stabilized::StepWork;
use super::*;
use crate::compressible::mixture::molecular::stabilized::IntervalFailure;

/// One physical interval shared by diffusion, shock transport and diffusion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepControl {
    /// Original Euler reconstruction, CFL and retry policy.
    pub advection: StepControl,
    /// Source-owned transport substep and positivity budgets.
    pub molecular: MolecularControl,
    /// Maximum complete-step restarts when advection shortens its interval.
    pub maximum_coupling_retries: usize,
}
impl Default for MolecularStepControl {
    fn default() -> Self {
        Self {
            advection: StepControl::default(),
            molecular: MolecularControl::default(),
            maximum_coupling_retries: 16,
        }
    }
}

/// One completely accepted Fourier/Newtonian/Euler interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport {
    /// Same physical duration represented by the complete split.
    pub elapsed: Time,
    /// Number of rejected whole-step trials.
    pub coupling_retries: usize,
    /// First molecular half step.
    pub first: MolecularReport,
    /// Original advection report between molecular half steps.
    pub advection: StepReport,
    /// Second molecular half step.
    pub second: MolecularReport,
    /// Complete original field.
    pub before: Totals,
    /// Complete accepted field.
    pub after: Totals,
    /// Combined inward molecular and advective boundary transfer.
    pub boundary_change: Totals,
    /// Applied body impulse/work over the same complete physical interval.
    pub source_change: ForceChange2d,
    /// `after - before - boundary_change - source_change`.
    pub conservation_residual: Totals,
}

/// Complete coupled physical interval and actual RKL2 work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedMolecularStepReport {
    /// Original accepted molecular/advection/molecular balances.
    pub transport: MolecularStepReport,
    /// All molecular evaluations across both halves and discarded whole trials.
    pub operator_evaluations: usize,
    /// Largest stage count in the final accepted pair of molecular halves.
    pub maximum_stages_used: usize,
    /// Largest scaled parabolic number in the final accepted molecular halves.
    pub maximum_scaled_parabolic_number: f64,
}

/// Atomic material-driven two-dimensional compressible viscous/thermal flow.
///
/// A symmetric molecular/advection/molecular composition uses the exact duration
/// accepted by the original Euler owner. If its CFL changes after heating, the
/// complete trial restarts from the original field. Both molecular stages use
/// SSP RK2 by default or opt-in RKL2; MUSCL Euler permits second-order splitting.
/// Limiting and first-order retries reduce order. This is not an entropy,
/// turbulence, shock-thickness or experimental-validation certificate.
pub struct MolecularEuler2d<L> {
    accepted: Euler2d,
    trial: Euler2d,
    molecular: MolecularTransport2d<L>,
    boundary_impulses: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
    last_step: Option<MolecularStepReport>,
}
impl<L: TransportLaw> MolecularEuler2d<L> {
    /// Apply one authored wall to both molecular half steps. Its shear, work
    /// and heat join the original Euler wall pressure in the accepted account.
    /// A changed condition invalidates the previous complete-step report;
    /// identical reconfiguration preserves the accepted interval and loads.
    pub fn with_wall(
        mut self,
        side: BoundarySide,
        wall: MolecularWall2d,
    ) -> Result<Self, MolecularError> {
        let changed = self.molecular.wall(side) != wall;
        self.molecular = self.molecular.with_wall(side, wall)?;
        if changed {
            self.last_step = None;
        }
        Ok(self)
    }

    /// Authored molecular wall on one original physical side.
    pub fn wall(&self, side: BoundarySide) -> MolecularWall2d {
        self.molecular.wall(side)
    }

    /// Preserve the original accepted Euler field and bind its identical source law.
    pub fn new(
        flow: Euler2d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let molecular = MolecularTransport2d::new(&flow, law, &mut cancelled)?;
        let mut cells = allocated(flow.states.len(), State::default())?;
        cells.copy_from_slice(&flow.states);
        let mut trial = Euler2d::new(
            flow.gas.clone(),
            flow.grid,
            cells,
            flow.boundaries,
            &mut cancelled,
        )?;
        if let Some(balance) = &flow.hydrostatic {
            trial.bind_hydrostatic_reference(
                balance.axis,
                balance.reference,
                balance.quadrature,
                &mut cancelled,
            )?;
        }
        let faces = flow.boundary_impulses.len();
        Ok(Self {
            accepted: flow,
            trial,
            molecular,
            boundary_impulses: allocated(faces, Totals::default())?,
            boundary_trial: allocated(faces, Totals::default())?,
            last_step: None,
        })
    }
    /// Complete accepted physical field and clock. Its Euler step metadata
    /// describes the middle advective stage; use this owner's report for totals.
    pub const fn flow(&self) -> &Euler2d {
        &self.accepted
    }
    /// Exact immutable material/collision source bound to this field.
    pub fn law(&self) -> &L {
        self.molecular.law()
    }
    /// Last completely accepted coupled interval.
    pub const fn last_step(&self) -> Option<&MolecularStepReport> {
        self.last_step.as_ref()
    }
    /// Combined outward transfer, including the stationary-wall viscous impulse.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.accepted.grid.boundary_range(side)]
    }

    /// Advance at most the requested duration, publishing only a complete split.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, MolecularError> {
        self.step_impl(maximum_duration, controls, None, None, cancelled)
    }

    /// Drive the original viscous/thermal flow with acceleration or volumetric
    /// force. Molecular halves surround the same conservative forced Euler
    /// interval; a shortened or refused central step restarts the whole split
    /// from the accepted state. Source work and wall transfer remain distinct.
    /// Allocate reusable scratch through `flow().forcing_workspace(maximum_bytes)`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        force: BodyForce2d<'_>,
        force_controls: ForceControl2d,
        workspace: &mut ForcingWorkspace2d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step(maximum_duration, controls, cancelled);
        }
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            None,
            cancelled,
        )
    }

    /// Compose bounded RKL2 halves around the original Euler step.
    ///
    /// One evaluation budget covers both halves and every discarded complete
    /// trial. Hyperbolic CFL, source gates and complete-step conservation retain
    /// their original owners. A failure leaves field, clock, loads and report
    /// unchanged; `last_step()` retains the accepted physical balances.
    pub fn step_stabilized(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport, MolecularError> {
        let mut work = StepWork::new(stabilization);
        let transport =
            self.step_impl(maximum_duration, controls, None, Some(&mut work), cancelled)?;
        Ok(StabilizedMolecularStepReport {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    /// Apply the same stabilization around conservative forced Euler evolution.
    /// The force workspace, Mach gate and source-work accounting are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_stabilized_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        stabilization: MolecularStabilization,
        force: BodyForce2d<'_>,
        force_controls: ForceControl2d,
        workspace: &mut ForcingWorkspace2d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedMolecularStepReport, MolecularError> {
        force_controls.admit::<Error>()?;
        let active = force.admit::<Error>(self.accepted.states.len(), &mut cancelled)?;
        workspace.admit(self.accepted.states.len())?;
        if !active {
            return self.step_stabilized(maximum_duration, controls, stabilization, cancelled);
        }
        let mut work = StepWork::new(stabilization);
        let transport = self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            Some(&mut work),
            cancelled,
        )?;
        Ok(StabilizedMolecularStepReport {
            transport,
            operator_evaluations: work.evaluations,
            maximum_stages_used: work.maximum_stages_used,
            maximum_scaled_parabolic_number: work.maximum_scaled_parabolic_number,
        })
    }

    fn molecular_half(
        &mut self,
        duration: Time,
        controls: MolecularControl,
        stabilization: Option<&mut StepWork>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, MolecularError> {
        match stabilization {
            None => self
                .molecular
                .advance(&mut self.trial, duration, controls, cancelled),
            Some(work) => {
                let report = self.molecular.advance_stabilized_with_budget(
                    &mut self.trial,
                    duration,
                    controls,
                    work.controls,
                    &mut work.evaluations,
                    cancelled,
                )?;
                work.record(
                    report.maximum_stages_used,
                    report.maximum_scaled_parabolic_number,
                );
                Ok(report.transport)
            }
        }
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut forcing: Option<(BodyForce2d<'_>, ForceControl2d, &mut ForcingWorkspace2d)>,
        mut stabilization: Option<&mut StepWork>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport, MolecularError> {
        poll(&mut cancelled)?;
        controls.molecular.validate()?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
        {
            return Err(MolecularError::Invalid(
                "coupled molecular duration or controls",
            ));
        }
        if let Some(work) = stabilization.as_deref() {
            work.controls.controller(controls.molecular)?;
            self.molecular
                .prepare_stabilized(work.controls.maximum_extra_workspace_bytes)?;
        }
        let before = checked_totals(
            &self.accepted.states,
            self.accepted.grid.cell_volume(),
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
        for q in &self.accepted.states {
            let p = q.primitive(&self.accepted.gas, &mut cancelled)?;
            if let Some((_, force_controls, _)) = &forcing {
                force_controls.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            rate = rate.max(
                (p.velocity[0].abs() + p.sound_speed) / self.accepted.grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / self.accepted.grid.spacing[1],
            );
        }
        rate = rate.max(super::super::pressure_outlet_rate(
            &self.accepted.gas,
            &self.accepted.grid,
            self.accepted.boundaries,
            &self.accepted.states,
            self.accepted.hydrostatic.as_ref(),
            &mut cancelled,
        )?);
        let mut dt = maximum_duration.value().min(controls.advection.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if let Some(work) = stabilization.as_deref_mut() {
                work.reset_extrema();
            }
            if !positive(dt * 0.5) {
                return Err(MolecularError::Invalid("unrepresentable coupled duration"));
            }
            self.trial.states.copy_from_slice(&self.accepted.states);
            self.trial.time = self.accepted.time;
            self.trial.last_step = self.accepted.last_step;
            self.trial
                .boundary_impulses
                .copy_from_slice(&self.accepted.boundary_impulses);
            let mut shortened = None;
            let attempt = (|| {
                let first = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    &mut cancelled,
                )?;
                self.boundary_trial
                    .copy_from_slice(&self.molecular.boundary_impulses);
                let advection = match forcing.as_mut() {
                    Some((force, force_controls, workspace)) => self.trial.step_impl(
                        Time::new(dt),
                        controls.advection,
                        Some((*force, *force_controls, &mut **workspace)),
                        &mut cancelled,
                    )?,
                    None => self
                        .trial
                        .step(Time::new(dt), controls.advection, &mut cancelled)?,
                };
                if advection.elapsed.value() != dt {
                    shortened = Some(advection.elapsed.value());
                    return Err(MolecularError::Invalid(
                        "Euler CFL shortened the split interval",
                    ));
                }
                for (total, impulse) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.trial.boundary_impulses)
                {
                    *total = add_totals(*total, *impulse);
                }
                let second = self.molecular_half(
                    Time::new(0.5 * dt),
                    controls.molecular,
                    stabilization.as_deref_mut(),
                    &mut cancelled,
                )?;
                for (total, impulse) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.molecular.boundary_impulses)
                {
                    *total = add_totals(*total, *impulse);
                }
                if let Some((_, force_controls, _)) = &forcing {
                    let mut endpoint_rate = 0.0_f64;
                    for state in &self.trial.states {
                        let p = state.primitive(&self.trial.gas, &mut cancelled)?;
                        force_controls.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
                        let local: f64 = (0..2)
                            .map(|axis| {
                                (p.velocity[axis].abs() + p.sound_speed)
                                    / self.trial.grid.spacing[axis]
                            })
                            .sum();
                        endpoint_rate = endpoint_rate.max(local);
                    }
                    cartesian::check_cfl::<Error>(endpoint_rate * dt, controls.advection.cfl)?;
                }
                let after = checked_totals(
                    &self.trial.states,
                    self.accepted.grid.cell_volume(),
                    &mut cancelled,
                )?;
                let outward = sum_arrays(
                    self.boundary_trial.iter().map(|q| q.array()),
                    &mut cancelled,
                )?;
                let boundary_change = Totals::from_array(outward.map(|v| -v));
                let source_change = advection.source_change;
                let source_totals = Totals {
                    mass: 0.0,
                    momentum: source_change.impulse,
                    energy: source_change.work,
                };
                let conservation_residual =
                    balance(before, after, add_totals(boundary_change, source_totals))?;
                poll(&mut cancelled)?;
                Ok(MolecularStepReport {
                    elapsed: Time::new(dt),
                    coupling_retries: retry,
                    first,
                    advection,
                    second,
                    before,
                    after,
                    boundary_change,
                    source_change,
                    conservation_residual,
                })
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    self.boundary_impulses.copy_from_slice(&self.boundary_trial);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(error) if error.is_terminal() => return Err(error),
                Err(error) => {
                    if retry == controls.maximum_coupling_retries {
                        return Err(MolecularError::CouplingRetries {
                            attempts: retry + 1,
                            last: Box::new(error),
                        });
                    }
                    dt = shortened.unwrap_or(dt * 0.5);
                }
            }
        }
        unreachable!("bounded coupled retries return on the final attempt")
    }
}
