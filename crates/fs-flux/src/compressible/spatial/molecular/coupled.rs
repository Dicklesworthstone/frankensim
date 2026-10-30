use super::*;

/// One completely accepted three-dimensional molecular/Euler/molecular step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularStepReport3d {
    /// Physical duration shared by all stages [s].
    pub elapsed: Time,
    /// Rejected complete trials before acceptance.
    pub coupling_retries: usize,
    /// First molecular half step.
    pub first: MolecularReport3d,
    /// Original Euler transport report between molecular half steps.
    pub advection: StepReport,
    /// Second molecular half step.
    pub second: MolecularReport3d,
    /// Original complete field totals.
    pub before: Totals,
    /// Accepted complete field totals.
    pub after: Totals,
    /// Combined inward molecular and advective boundary transfer.
    pub boundary_change: Totals,
    /// `after - before - boundary_change` for all five conservative components.
    pub conservation_residual: Totals,
}

/// Atomic source-driven three-dimensional compressible viscous/thermal flow.
///
/// A symmetric molecular/advection/molecular composition retains one actual
/// accepted duration. A shorter Euler step or a refused later stage restarts
/// the complete trial from the accepted field. Both molecular stages use SSP
/// RK2; smooth MUSCL Euler permits second-order splitting, while limiting or
/// first-order fallback reduces order. No entropy, turbulence, physical shock
/// thickness or experimental-accuracy certificate is implied.
pub struct MolecularEuler3d<L> {
    accepted: Euler3d,
    trial: Euler3d,
    molecular: MolecularTransport3d<L>,
    boundary_impulses: Vec<BoundaryImpulse>,
    boundary_trial: Vec<BoundaryImpulse>,
    last_step: Option<MolecularStepReport3d>,
}

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Bind the identical original source while preserving the complete field.
    pub fn new(
        flow: Euler3d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let molecular = MolecularTransport3d::new(&flow, law, &mut cancelled)?;
        let mut cells = allocated(flow.states.len(), State::default())?;
        cells.copy_from_slice(&flow.states);
        let trial = Euler3d::new(
            flow.gas.clone(),
            flow.grid,
            cells,
            flow.boundaries,
            &mut cancelled,
        )?;
        let faces = flow.boundary_impulses.len();
        let result = Self {
            accepted: flow,
            trial,
            molecular,
            boundary_impulses: allocated(faces, Totals::default())?,
            boundary_trial: allocated(faces, Totals::default())?,
            last_step: None,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Accepted field and physical clock. Its Euler report describes the middle
    /// advective stage; this owner's report and loads cover the complete split.
    pub const fn flow(&self) -> &Euler3d {
        &self.accepted
    }

    /// Exact immutable source law, including retained material/collision data.
    pub fn law(&self) -> &L {
        self.molecular.law()
    }

    /// Last completely accepted physical interval.
    pub const fn last_step(&self) -> Option<&MolecularStepReport3d> {
        self.last_step.as_ref()
    }

    /// Combined outward transfer, including both molecular support impulses.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.accepted.grid.boundary_range(side)]
    }

    /// Advance at most the requested duration, publishing only a complete split.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularStepReport3d, MolecularError> {
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
            let local: f64 = (0..3)
                .map(|a| (p.velocity[a].abs() + p.sound_speed) / self.accepted.grid.spacing[a])
                .sum();
            rate = rate.max(local);
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
                for (i, (total, impulse)) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.trial.boundary_impulses)
                    .enumerate()
                {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    *total = add_totals(*total, *impulse);
                }
                let second = self.molecular.advance(
                    &mut self.trial,
                    Time::new(0.5 * dt),
                    controls.molecular,
                    &mut cancelled,
                )?;
                for (i, (total, impulse)) in self
                    .boundary_trial
                    .iter_mut()
                    .zip(&self.molecular.boundary_impulses)
                    .enumerate()
                {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
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
                Ok(MolecularStepReport3d {
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
