//! Source-owned molecular transport on the actual old and new ALE volumes.

use super::*;
use crate::compressible::spatial::MeshMotion;

/// Complete evolving-species molecular/ALE interval and direct wall transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingMolecularStepReport3d {
    /// Original species, three momenta, energy and complete boundary balances.
    pub transport: MolecularStepReport3d,
    /// Signed molecular stress work into the gas on XMin, XMax, YMin, YMax,
    /// ZMin, ZMax [J]. Pressure work is already in `transport.advection`.
    pub molecular_wall_work_j: [f64; 6],
    /// Signed Fourier heat and the Navier wall heat partition into the gas [J],
    /// in the same side order. Evaluated before adding mechanical wall work.
    pub wall_heat_j: [f64; 6],
}

/// A complete fixed-duration trial must be restarted if native ALE shortens it.
pub(super) struct MovingTrialFailure {
    pub(super) error: DissipationError,
    pub(super) shortened: Option<f64>,
}
impl From<DissipationError> for MovingTrialFailure {
    fn from(error: DissipationError) -> Self {
        Self {
            error,
            shortened: None,
        }
    }
}

impl From<Error> for MovingTrialFailure {
    fn from(error: Error) -> Self {
        Self::from(DissipationError::from(error))
    }
}

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Borrow instantaneous source-owned mixture stress on the moving volume.
    ///
    /// `motion` supplies the current normal plane velocities. Tangential slip,
    /// temperature and mixture composition use the accepted original owners.
    /// After finite-piston stepping, pass endpoint piston velocities rather
    /// than the interval's midpoint speeds. No historical velocity is inferred
    /// or retained, and [`Self::stress_field`] retains stationary semantics.
    ///
    /// The original Maxwell--Stefan/Fourier/full-stress kernel is evaluated in
    /// prepared buffers under `controls`. Cells, geometry, physical time and
    /// all accepted interval loads remain unchanged. The same closed/periodic
    /// domain as [`Self::step_moving`] is required; SGS, exterior profiles and
    /// finite tangential wall inertia are refused. Failure returns no view;
    /// reusable source caches and scratch may change.
    pub fn stress_field_moving(
        &mut self,
        controls: DissipationControl,
        motion: MeshMotion,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::compressible::spatial::molecular::StressField3d<'_>, DissipationError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        self.admit_moving_domain(motion)?;
        let operator = &mut self.molecular;
        spatial_with_motion(
            &mut operator.source,
            &operator.grid,
            operator.boundaries,
            operator.pressure_outlets,
            reservoir::View::from(&self.accepted.reservoirs),
            &operator.walls,
            &self.accepted.states,
            &mut operator.work.primitive,
            &mut operator.work.gradient,
            &mut operator.work.viscosity,
            &mut operator.work.flux,
            #[cfg(feature = "experimental-sgs")]
            None,
            controls,
            Some(motion),
            None,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(
            crate::compressible::spatial::molecular::StressField3d::mixture(
                &operator.grid,
                operator.boundaries,
                self.accepted.time,
                &operator.work.primitive,
                &operator.work.flux,
                &self.accepted.states,
                self.accepted.gas.species_count(),
            ),
        )
    }

    /// Advance the original evolving mixture with molecular transport on a
    /// translating or single-axis deforming Cartesian box.
    ///
    /// An SSP molecular half on the old geometry, the original conservative
    /// ALE interval, and a molecular half on the new geometry share one trial.
    /// Maxwell--Stefan species transport, enthalpy diffusion, Fourier heat and
    /// all nine velocity derivatives use the original source and prepared
    /// buffers. Wall species flow is zero. Full moving-wall velocity supplies
    /// normal stress and signed mechanical work, including on free-slip walls.
    ///
    /// Closed reflecting and periodic domains are supported, with prescribed
    /// tangential no-slip/Navier motion and held normal mesh velocities. Open
    /// reservoirs, retained SGS histories, finite tangential wall inertia and
    /// simultaneous deformation of multiple axes need their own composition.
    /// A stationary mesh uses the same numerical operators and reports actual
    /// prescribed belt work and wall heat. It does not waive this admission.
    ///
    /// Native caloric, diffusive and relative acoustic limits remain active.
    /// Any shortened phase restarts both halves. Rejection or cancellation
    /// preserves cells, grid, clock, source binding, loads and prior reports.
    pub fn step_moving(
        &mut self,
        maximum_duration: Time,
        controls: MolecularStepControl,
        motion: MeshMotion,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingMolecularStepReport3d, DissipationError> {
        poll(&mut cancelled)?;
        self.admit_moving(maximum_duration.value(), controls, motion)?;
        let n = self.accepted.gas.species_count();
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            n,
            &mut cancelled,
        )?;
        let rate = relative_rate(&self.accepted, motion, &mut cancelled)?;
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.advection.cfl / rate);
        let mut last = DissipationError::Invalid("no moving molecular mixture trial");
        let mut last_dt = dt;
        for retry in 0..=controls.maximum_coupling_retries {
            match self.trial_moving(dt, controls, motion, before, retry, &mut cancelled) {
                Ok(report) => {
                    self.publish_moving(report);
                    return Ok(report);
                }
                Err(MovingTrialFailure {
                    error: DissipationError::Cancelled,
                    ..
                }) => {
                    return Err(DissipationError::Cancelled);
                }
                Err(failure) => {
                    last = failure.error;
                    last_dt = dt;
                    dt = failure.shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        Err(DissipationError::StepFailure {
            elapsed_time_s: 0.0,
            step_s: last_dt,
            last: Box::new(last),
        })
    }

    pub(super) fn admit_moving(
        &self,
        maximum_duration: f64,
        controls: MolecularStepControl,
        motion: MeshMotion,
    ) -> Result<(), DissipationError> {
        controls.molecular.validate()?;
        if !positive(maximum_duration)
            || controls.maximum_coupling_retries > 32
            || !positive(controls.advection.cfl)
            || controls.advection.cfl > 0.45
            || controls.advection.maximum_halvings > 32
        {
            return Err(DissipationError::Invalid(
                "moving molecular mixture duration or controls",
            ));
        }
        self.admit_moving_domain(motion)
    }

    fn admit_moving_domain(&self, motion: MeshMotion) -> Result<(), DissipationError> {
        admit_workspace(self.workspace_bytes(), self.maximum_workspace_bytes)?;
        self.molecular.admit_moving_field(&self.accepted, motion)
    }

    /// Fixed-duration, unpublished complete split. Finite-body coupling may
    /// probe this repeatedly; neither an accepted source binding nor an
    /// accepted boundary impulse is replaced until `publish_moving`.
    pub(super) fn trial_moving(
        &mut self,
        dt: f64,
        controls: MolecularStepControl,
        motion: MeshMotion,
        before: Totals,
        retry: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingMolecularStepReport3d, MovingTrialFailure> {
        poll(cancelled)?;
        if !positive(0.5 * dt)
            || !(self.accepted.time + dt).is_finite()
            || self.accepted.time + dt <= self.accepted.time
        {
            return Err(DissipationError::Invalid(
                "unrepresentable moving molecular mixture increment",
            )
            .into());
        }
        // Admit the endpoint before spending either half on a collapsing box.
        // The original ALE transaction still owns its exact stage geometry.
        if motion != MeshMotion::default() {
            motion
                .grid_after(&self.accepted.grid, dt)
                .map_err(|error| {
                    DissipationError::Invalid(match error {
                        fs_gas::dynamics::Error::Invalid(reason) => reason,
                        _ => "moving molecular geometry",
                    })
                })?;
        }
        copy_physical(&self.accepted, &mut self.trial, cancelled)?;
        let (first, first_transfer) = self.molecular.advance_moving_unpublished(
            &mut self.trial,
            Time::new(0.5 * dt),
            controls.molecular,
            self.maximum_workspace_bytes,
            motion,
            cancelled,
        )?;
        copy_values(
            self.molecular.trial_impulses(),
            &mut self.boundary_trial,
            cancelled,
        )?;
        let advection = self
            .trial
            .step_moving(Time::new(dt), controls.advection, motion, &mut *cancelled)
            .map_err(DissipationError::from)?;
        if advection.elapsed.value() != dt {
            return Err(MovingTrialFailure {
                error: DissipationError::Invalid("ALE shortened moving molecular mixture interval"),
                shortened: Some(advection.elapsed.value()),
            });
        }
        add_impulses(
            &self.trial.boundary_impulses,
            &mut self.boundary_trial,
            cancelled,
        )?;
        let (second, second_transfer) = self.molecular.advance_moving_unpublished(
            &mut self.trial,
            Time::new(0.5 * dt),
            controls.molecular,
            self.maximum_workspace_bytes,
            motion,
            cancelled,
        )?;
        add_impulses(
            self.molecular.trial_impulses(),
            &mut self.boundary_trial,
            cancelled,
        )?;
        // The final molecular half can heat or accelerate the gas after the
        // last wave stage. Admit its actual new-volume acoustic rate too.
        let final_rate = relative_rate(&self.trial, motion, cancelled)?;
        stepping::check_cfl::<Error>(final_rate * dt, controls.advection.cfl)
            .map_err(DissipationError::from)?;
        let n = self.trial.gas.species_count();
        let width = n + 4;
        let after = checked_totals(&self.trial.states, &self.trial.grid, n, cancelled)?;
        let outward = sum_rows(&self.boundary_trial, width, 1.0, cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|x| -x), n);
        let source_change = advection.source_change;
        let mut residual = [0.0; WIDTH_LIMIT];
        for (d, value) in residual[..width].iter_mut().enumerate() {
            *value = after.component(d) - before.component(d) - boundary_change.component(d);
        }
        for axis in 0..3 {
            residual[n + axis] -= source_change.impulse[axis];
        }
        residual[n + 3] -= source_change.work;
        if !residual[..width].iter().all(|x| x.is_finite()) {
            return Err(
                DissipationError::Invalid("moving molecular mixture balance overflow").into(),
            );
        }
        let conservation_residual = Totals::from_row(&residual[..width], n);
        let absolute_energy_residual_j = conservation_residual
            .absolute_energy(&self.trial.gas)
            .map_err(DissipationError::from)?;
        let mut transfer = first_transfer;
        transfer.add_scaled(second_transfer, 1.0, &mut WallTransfer::default())?;
        poll(cancelled)?;
        Ok(MovingMolecularStepReport3d {
            transport: MolecularStepReport3d {
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
                absolute_energy_residual_j,
            },
            molecular_wall_work_j: transfer.mechanical,
            wall_heat_j: transfer.heat,
        })
    }

    pub(super) fn publish_moving(&mut self, report: MovingMolecularStepReport3d) {
        std::mem::swap(&mut self.accepted, &mut self.trial);
        std::mem::swap(&mut self.boundary_impulses, &mut self.boundary_trial);
        self.molecular
            .publish_moving(self.accepted.grid, report.transport.second);
        self.last_step = Some(report.transport);
    }
}

fn relative_rate(
    flow: &MixtureEuler3d,
    motion: MeshMotion,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let n = flow.gas.species_count();
    let mut rate: f64 = 0.0;
    for (index, row) in flow.states.chunks_exact(n + 4).enumerate() {
        let p = cell(row, n).primitive(&flow.gas, &mut *cancelled)?;
        let local = if motion == MeshMotion::default() {
            forcing::acoustic_rate(p, &flow.grid)?
        } else {
            let coordinates = flow.grid.coordinates(index);
            (0..3)
                .map(|axis| {
                    let lo = motion.velocity(&flow.grid, axis, coordinates[axis]);
                    let hi = motion.velocity(&flow.grid, axis, coordinates[axis] + 1);
                    ((p.velocity[axis] - lo)
                        .abs()
                        .max((p.velocity[axis] - hi).abs())
                        + p.sound_speed)
                        / flow.grid.spacing[axis]
                })
                .sum::<f64>()
        };
        if !local.is_finite() {
            return Err(DissipationError::Invalid(
                "moving molecular relative CFL overflow",
            ));
        }
        rate = rate.max(local);
    }
    Ok(rate)
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Shared source/topology admission for molecular and reacting moving owners.
    pub(in crate::compressible::spatial::mixture) fn admit_moving_field(
        &self,
        flow: &MixtureEuler3d,
        motion: MeshMotion,
    ) -> Result<(), DissipationError> {
        flow.admit_motion(motion)?;
        if self.grid != flow.grid || self.gas() != flow.gas() || self.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "moving molecular source binding differs",
            ));
        }
        if self.inertia.is_some() {
            return Err(DissipationError::Invalid(
                "moving molecular volumes do not compose finite tangential wall inertia",
            ));
        }
        #[cfg(feature = "experimental-sgs")]
        if self.sgs.is_some() {
            return Err(DissipationError::Invalid(
                "moving molecular volumes do not compose retained subgrid models",
            ));
        }
        for axis in 0..3 {
            if let AxisBoundaries::Ends { lower, upper } = flow.boundaries.axis(axis)
                && [lower, upper]
                    .into_iter()
                    .any(|bc| bc != BoundaryCondition::Reflecting)
            {
                return Err(DissipationError::Invalid(
                    "moving molecular volumes require closed walls or periodic boundaries",
                ));
            }
        }
        Ok(())
    }

    /// The composing owner has admitted the complete trial and final checkpoint.
    pub(in crate::compressible::spatial::mixture) fn publish_moving(
        &mut self,
        grid: CartesianGrid,
        report: MolecularReport,
    ) {
        self.grid = grid;
        self.publish(report, None);
    }

    /// One native SSP interval at the split's current geometry. Only private
    /// trial cells change here; the accepted molecular binding never moves.
    pub(in crate::compressible::spatial::mixture) fn advance_moving_unpublished(
        &mut self,
        flow: &mut MixtureEuler3d,
        duration: Time,
        controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(MolecularReport, WallTransfer), DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        let before =
            self.prepare_interval_on_grid(flow, duration, maximum_workspace_bytes, cancelled)?;
        let grid = flow.grid;
        let width = self.gas().species_count() + 4;
        let mut first_rates = WallTransfer::default();
        let mut transfer = WallTransfer::default();
        let mut transfer_correction = WallTransfer::default();
        let interval = interval::advance(
            duration.value(),
            controls,
            cancelled,
            |action, _, cancelled| -> Result<Stage, DissipationError> {
                match action {
                    Action::First => self.moving_spatial(
                        &grid,
                        &flow.reservoirs,
                        motion,
                        false,
                        controls,
                        &mut first_rates,
                        cancelled,
                    ),
                    Action::Trial { dt, first, cached } => {
                        let first = if cached {
                            first
                        } else {
                            self.moving_spatial(
                                &grid,
                                &flow.reservoirs,
                                motion,
                                false,
                                controls,
                                &mut first_rates,
                                cancelled,
                            )?
                        };
                        parabolic(first.rate * dt, controls.parabolic_cfl)?;
                        update(
                            self.source.gas(),
                            &grid,
                            &self.work.current,
                            &self.work.current,
                            &self.work.flux,
                            &mut self.work.stage,
                            dt,
                            false,
                            None,
                            cancelled,
                        )?;
                        capture(
                            &grid,
                            width,
                            &self.work.flux,
                            0.5 * dt,
                            false,
                            &mut self.work.boundary_step,
                            cancelled,
                        )?;
                        let mut second_rates = WallTransfer::default();
                        let second = self.moving_spatial(
                            &grid,
                            &flow.reservoirs,
                            motion,
                            true,
                            controls,
                            &mut second_rates,
                            cancelled,
                        )?;
                        parabolic(second.rate * dt, controls.parabolic_cfl)?;
                        update(
                            self.source.gas(),
                            &grid,
                            &self.work.current,
                            &self.work.stage,
                            &self.work.flux,
                            &mut self.work.candidate,
                            dt,
                            true,
                            None,
                            cancelled,
                        )?;
                        capture(
                            &grid,
                            width,
                            &self.work.flux,
                            0.5 * dt,
                            true,
                            &mut self.work.boundary_step,
                            cancelled,
                        )?;
                        // Calculate both ledgers privately before an accepted
                        // substep can change interval accumulation or cells.
                        let mut next_transfer = transfer;
                        let mut next_correction = transfer_correction;
                        next_transfer.add_scaled(first_rates, 0.5 * dt, &mut next_correction)?;
                        next_transfer.add_scaled(second_rates, 0.5 * dt, &mut next_correction)?;
                        for (index, ((total, correction), change)) in self
                            .work
                            .boundary_trial
                            .iter_mut()
                            .zip(&mut self.work.boundary_correction)
                            .zip(&self.work.boundary_step)
                            .enumerate()
                        {
                            if index % 4096 == 0 {
                                poll(cancelled)?;
                            }
                            let adjusted = change - *correction;
                            let next = *total + adjusted;
                            *correction = (next - *total) - adjusted;
                            *total = next;
                        }
                        transfer = next_transfer;
                        transfer_correction = next_correction;
                        std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                        Ok(first.merge(second))
                    }
                }
            },
        )?;
        let report = self.complete_interval(flow, duration, before, interval, None, cancelled)?;
        Ok((report, transfer))
    }

    #[allow(clippy::too_many_arguments)]
    fn moving_spatial(
        &mut self,
        grid: &CartesianGrid,
        reservoirs: &Reservoirs,
        motion: MeshMotion,
        second_stage: bool,
        controls: DissipationControl,
        rates: &mut WallTransfer,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Stage, DissipationError> {
        spatial_with_motion(
            &mut self.source,
            grid,
            self.boundaries,
            self.pressure_outlets,
            reservoirs,
            &self.walls,
            if second_stage {
                &self.work.stage
            } else {
                &self.work.current
            },
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.viscosity,
            &mut self.work.flux,
            #[cfg(feature = "experimental-sgs")]
            None,
            controls,
            Some(motion),
            Some(rates),
            cancelled,
        )
    }
}
