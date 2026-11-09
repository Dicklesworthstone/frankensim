//! Native molecular stress and heat on the actual moving Cartesian gas volume.

use super::*;

/// A complete molecular/ALE/molecular interval on one accepted moving volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingMolecularStepReport3d {
    /// Original conservative transport and source-admissibility accounts.
    pub transport: MolecularStepReport3d,
    /// Molecular traction work into the gas [J], in [`BoundarySide::ALL`] order.
    /// Includes all three laboratory wall-velocity components. Pressure work
    /// belongs to `transport.advection` and is not included here.
    pub molecular_wall_work_j: [f64; 6],
    /// Nonmechanical energy into the gas [J], in [`BoundarySide::ALL`] order.
    /// Fourier heat and the Navier heat partition use the same SSP quadrature
    /// as the face flux. Heat sent to a wall is negative. These values are
    /// accumulated directly, without subtracting work from total energy.
    pub wall_heat_j: [f64; 6],
}

/// A fixed-duration mechanics proposal must refuse an internally shortened
/// Euler interval so that both molecular halves can be retried consistently.
pub(super) struct MovingTrialFailure {
    pub(super) error: MolecularError,
    pub(super) shortened: Option<f64>,
}

impl From<MolecularError> for MovingTrialFailure {
    fn from(error: MolecularError) -> Self {
        Self {
            error,
            shortened: None,
        }
    }
}

impl From<Error> for MovingTrialFailure {
    fn from(error: Error) -> Self {
        MolecularError::from(error).into()
    }
}

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Advance source-owned Newtonian stress, Fourier heat and ALE advection
    /// on a translating or single-axis deforming box.
    ///
    /// The first SSP molecular half uses the original cell volumes and face
    /// areas; the second uses the actual accepted ALE geometry. Each moving
    /// wall's normal velocity comes from `motion`, while its authored model
    /// retains the tangential slip and thermal laws. The same laboratory
    /// traction supplies all momentum components and mechanical work.
    ///
    /// A shorter Euler interval, source refusal or endpoint CFL violation
    /// restarts the complete trial. Cells, geometry, time, reports and loads
    /// publish together after the final cancellation check. Existing cell and
    /// face storage is reused. Zero motion retains the stationary kernels and
    /// their arithmetic, including stationary thermal and tangential wall work.
    ///
    /// This entry admits periodic directions and reflecting molecular walls.
    /// Exterior profiles, retained hydrostatics, open ends and SGS closures
    /// require their own moving-domain laws and are refused here.
    pub fn step_moving(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        motion: MeshMotion,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingMolecularStepReport3d, MolecularError> {
        poll(&mut cancelled)?;
        self.admit_moving(maximum_duration.value(), controls, motion)?;
        let before = checked_totals(
            &self.accepted.states,
            self.accepted.grid.cell_volume(),
            &mut cancelled,
        )?;
        let rate = relative_rate(&self.accepted, motion, &mut cancelled)?;
        let mut dt = maximum_duration.value().min(controls.advection.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            match self.trial_moving(dt, controls, motion, before, retry, &mut cancelled) {
                Ok(report) => {
                    poll(&mut cancelled)?;
                    self.publish_moving(report);
                    return Ok(report);
                }
                Err(failure) if failure.error.is_terminal() => return Err(failure.error),
                Err(failure) => {
                    if retry == controls.maximum_coupling_retries {
                        return Err(MolecularError::CouplingRetries {
                            attempts: retry + 1,
                            last: Box::new(failure.error),
                        });
                    }
                    dt = failure.shortened.unwrap_or(dt * 0.5);
                }
            }
        }
        unreachable!("bounded moving split retries return on the final attempt")
    }

    pub(super) fn admit_moving(
        &self,
        duration: f64,
        controls: MolecularStepControl,
        motion: MeshMotion,
    ) -> Result<(), MolecularError> {
        controls.molecular.validate()?;
        if !positive(duration * 0.5)
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
            || controls.maximum_coupling_retries > 32
        {
            return Err(MolecularError::Invalid(
                "moving molecular duration or controls",
            ));
        }
        motion.admit(&self.accepted)?;
        if self.molecular.reservoir_exchange
            || !self.molecular.uses_compatible_molecular_stress()
            || (0..3).any(|axis| match self.accepted.boundaries.axis(axis) {
                AxisBoundaries::Periodic => false,
                AxisBoundaries::Ends { lower, upper } => {
                    lower != BoundaryCondition::Reflecting || upper != BoundaryCondition::Reflecting
                }
            })
        {
            return Err(MolecularError::Invalid(
                "moving molecular transport requires closed compatible walls or periodicity",
            ));
        }
        #[cfg(feature = "experimental-sgs")]
        if self.molecular.sgs.is_some() {
            return Err(MolecularError::Invalid(
                "moving molecular SGS is not admitted",
            ));
        }
        Ok(())
    }

    /// Evaluate one held-duration moving gas proposal without publication.
    /// The finite-inertia owner also uses this exact transaction and its six
    /// complete boundary-face accounts, including transverse support loads.
    pub(super) fn trial_moving(
        &mut self,
        duration: f64,
        controls: MolecularStepControl,
        motion: MeshMotion,
        before: Totals,
        retry: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingMolecularStepReport3d, MovingTrialFailure> {
        poll(cancelled)?;
        self.admit_moving(duration, controls, motion)?;
        if !self.accepted.time.is_finite()
            || !((self.accepted.time + duration).is_finite())
            || self.accepted.time + duration <= self.accepted.time
        {
            return Err(MolecularError::Invalid("unrepresentable moving molecular clock").into());
        }
        // Admit the entire affine geometry before spending work on diffusion.
        motion.grid_after(&self.accepted.grid, duration)?;
        self.trial.grid = self.accepted.grid;
        self.trial.states.copy_from_slice(&self.accepted.states);
        self.trial.time = self.accepted.time;
        self.trial.last_step = self.accepted.last_step;
        self.trial
            .boundary_impulses
            .copy_from_slice(&self.accepted.boundary_impulses);
        let (first, mut transfer) = self.molecular.advance_at_geometry(
            &mut self.trial,
            Time::new(0.5 * duration),
            controls.molecular,
            timed::Clock::split(self.accepted.time, duration, false),
            None,
            Some(motion),
            false,
            &mut *cancelled,
        )?;
        self.boundary_trial
            .copy_from_slice(&self.molecular.work.boundary_trial);
        let advection = self.trial.step_moving(
            Time::new(duration),
            controls.advection,
            motion,
            &mut *cancelled,
        )?;
        if advection.elapsed.value() != duration {
            return Err(MovingTrialFailure {
                error: MolecularError::Invalid("Euler CFL shortened the moving split interval"),
                shortened: Some(advection.elapsed.value()),
            });
        }
        accumulate_boundaries(
            &mut self.boundary_trial,
            &self.trial.boundary_impulses,
            cancelled,
        )?;
        let (second, second_transfer) = self.molecular.advance_at_geometry(
            &mut self.trial,
            Time::new(0.5 * duration),
            controls.molecular,
            timed::Clock::split(self.accepted.time, duration, true),
            None,
            Some(motion),
            false,
            &mut *cancelled,
        )?;
        transfer.add_scaled(second_transfer, 1.0)?;
        accumulate_boundaries(
            &mut self.boundary_trial,
            &self.molecular.work.boundary_trial,
            cancelled,
        )?;
        // Source heating can increase the endpoint acoustic speed after ALE.
        check_cfl::<Error>(
            relative_rate(&self.trial, motion, cancelled)? * duration,
            controls.advection.cfl,
        )?;
        let after = checked_totals(&self.trial.states, self.trial.grid.cell_volume(), cancelled)?;
        let outward = sum_arrays(self.boundary_trial.iter().map(|q| q.array()), cancelled)?;
        let boundary_change = Totals::from_array(outward.map(|q| -q));
        let conservation_residual = balance(before, after, boundary_change)?;
        poll(cancelled)?;
        Ok(MovingMolecularStepReport3d {
            transport: MolecularStepReport3d {
                elapsed: Time::new(duration),
                coupling_retries: retry,
                first,
                advection,
                second,
                before,
                after,
                boundary_change,
                source_change: ForceChange3d::default(),
                conservation_residual,
            },
            molecular_wall_work_j: transfer.work,
            wall_heat_j: transfer.heat,
        })
    }

    /// No fallible operation may follow the outer owner's final checkpoint.
    pub(super) fn publish_moving(&mut self, report: MovingMolecularStepReport3d) {
        std::mem::swap(&mut self.accepted, &mut self.trial);
        self.molecular.grid = self.accepted.grid;
        self.molecular
            .boundary_impulses
            .copy_from_slice(&self.molecular.work.boundary_trial);
        self.boundary_impulses.copy_from_slice(&self.boundary_trial);
        self.last_step = Some(report.transport);
    }
}

fn accumulate_boundaries(
    total: &mut [BoundaryImpulse],
    increment: &[BoundaryImpulse],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for (i, (sum, value)) in total.iter_mut().zip(increment).enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        *sum = add_totals(*sum, *value);
    }
    Ok(())
}

fn relative_rate(
    flow: &Euler3d,
    motion: MeshMotion,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut rate = 0.0_f64;
    for (cell, state) in flow.states.iter().enumerate() {
        let primitive = state.primitive(&flow.gas, &mut *cancelled)?;
        let coordinate = flow.grid.coordinates(cell);
        let local: f64 = (0..3)
            .map(|axis| {
                let lower = motion.velocity(&flow.grid, axis, coordinate[axis]);
                let upper = motion.velocity(&flow.grid, axis, coordinate[axis] + 1);
                ((primitive.velocity[axis] - lower)
                    .abs()
                    .max((primitive.velocity[axis] - upper).abs())
                    + primitive.sound_speed)
                    / flow.grid.spacing[axis]
            })
            .sum();
        if !local.is_finite() {
            return Err(Error::Invalid("moving molecular relative CFL overflow"));
        }
        rate = rate.max(local);
    }
    Ok(rate)
}
