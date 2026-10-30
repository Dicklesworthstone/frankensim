use super::*;

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
    /// `after - before - boundary_change` for all conservative components.
    pub conservation_residual: Totals,
}

/// Atomic material-driven two-dimensional compressible viscous/thermal flow.
///
/// A symmetric molecular/advection/molecular composition uses the exact duration
/// accepted by the original Euler owner. If its CFL changes after heating, the
/// complete trial restarts from the original field. Both molecular stages use
/// SSP RK2; MUSCL Euler permits second-order splitting for smooth admitted data.
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
        let trial = Euler2d::new(
            flow.gas.clone(),
            flow.grid,
            cells,
            flow.boundaries,
            &mut cancelled,
        )?;
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
        let before = checked_totals(
            &self.accepted.states,
            self.accepted.grid.cell_volume(),
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
        for q in &self.accepted.states {
            let p = q.primitive(&self.accepted.gas, &mut cancelled)?;
            rate = rate.max(
                (p.velocity[0].abs() + p.sound_speed) / self.accepted.grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / self.accepted.grid.spacing[1],
            );
        }
        let mut dt = maximum_duration.value().min(controls.advection.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
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
                let first = self.molecular.advance(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    &mut cancelled,
                )?;
                self.boundary_trial
                    .copy_from_slice(&self.molecular.boundary_impulses);
                let advection =
                    self.trial
                        .step(Time::new(dt), controls.advection, &mut cancelled)?;
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
                let second = self.molecular.advance(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    &mut cancelled,
                )?;
                for (total, impulse) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.molecular.boundary_impulses)
                {
                    *total = add_totals(*total, *impulse);
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
                let conservation_residual = balance(before, after, boundary_change)?;
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
                Err(MolecularError::Cancelled) => return Err(MolecularError::Cancelled),
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
