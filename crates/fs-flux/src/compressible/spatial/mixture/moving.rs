//! Original-species momentum and energy on translating and piston-driven grids.
//!
//! The shared Cartesian topology, composition reconstruction and face storage
//! carry one source-owned ALE flux. SSP RK2 evolves extensive species masses,
//! all momenta and referenced energy with actual stage volumes and face areas.
//! The supported single-axis deformation has an affine volume, so this same
//! quadrature satisfies the discrete geometric conservation law.

use super::super::super::positive;
use super::*;

impl MixtureEuler3d {
    /// Advance this field on the same moving box as the frozen gas owner.
    ///
    /// Original species calorics, both transverse momenta and chemical energy
    /// references enter each moving-face Riemann solve. Reflecting planes use
    /// the source's nonlinear wall pressure with exactly zero species transfer
    /// and pressure work at the actual plane velocity. Interior neighbors use
    /// the same stored flux with opposite signs; periodic ends share its bits.
    ///
    /// Translation is allowed on all axes and at most one axis may expand or
    /// contract. Existing accepted/stage/face buffers are reused. Rejection or
    /// cancellation preserves cells, grid, clock, report and boundary loads.
    /// Zero motion follows [`Self::step`] including its cancellation sequence.
    ///
    /// This inviscid entry admits periodic translation, moving slip walls and
    /// unbound zero-gradient openings. It does not compose molecular transport,
    /// chemistry, prescribed exterior states/profiles, pressure outlets or
    /// contact-energy projection on a moving domain. Those require the same
    /// geometry and stage-time conventions before they can be combined.
    pub fn step_moving(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        motion: MeshMotion,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        if motion == MeshMotion::default() {
            return self.step(maximum_duration, controls, cancelled);
        }
        self.step_moving_impl(maximum_duration, controls, motion, None, cancelled)
    }

    /// Compose held mechanical source kicks with conservative moving volumes.
    ///
    /// The original exact kinetic-work kicks surround native ALE transport.
    /// Acceleration uses the actual species mass at each half; force density
    /// uses that half's actual old or new volume. Physical Mach admission and
    /// mesh-relative acoustic limits apply to source, RK and face states.
    /// Species and their energy references are unchanged by either kick, and
    /// the source impulse/work remain separate from boundary exchange.
    ///
    /// Every retry restarts both kicks and ALE from accepted storage. The
    /// supplied workspace is prepared once and retains the first kicked field.
    /// Zero force follows the unforced moving path; zero motion follows the
    /// original forced stationary path. Publication remains atomic.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        motion: MeshMotion,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        if motion == MeshMotion::default() {
            return self.step_with_forcing(
                maximum_duration,
                controls,
                force,
                force_controls,
                workspace,
                cancelled,
            );
        }
        stepping::poll::<Error>(&mut cancelled)?;
        force_controls.admit::<Error>()?;
        workspace.admit(self)?;
        if !force.admit::<Error>(self.cell_count(), &mut cancelled)? {
            return self.step_moving(maximum_duration, controls, motion, cancelled);
        }
        self.step_moving_impl(
            maximum_duration,
            controls,
            motion,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }

    fn step_moving_impl(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        motion: MeshMotion,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        self.admit_motion(motion)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("volumetric duration, CFL or retry controls"));
        }
        let n = self.gas.species_count();
        let width = n + 4;
        let grid = self.grid;
        let mut rate: f64 = 0.0;
        for (k, (q, p)) in self
            .states
            .chunks_exact(width)
            .zip(&mut self.work.primitive)
            .enumerate()
        {
            *p = cell(q, n).primitive(&self.gas, &mut cancelled)?;
            if let Some((_, control, _)) = forcing.as_ref() {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            let c = grid.coordinates(k);
            let local: f64 = (0..3)
                .map(|axis| {
                    let lo = motion.velocity(&grid, axis, c[axis]);
                    let hi = motion.velocity(&grid, axis, c[axis] + 1);
                    ((p.velocity[axis] - lo)
                        .abs()
                        .max((p.velocity[axis] - hi).abs())
                        + p.sound_speed)
                        / grid.spacing[axis]
                })
                .sum();
            if !local.is_finite() {
                return Err(Error::Invalid(
                    "moving mixture initial relative CFL overflow",
                ));
            }
            rate = rate.max(local);
        }
        let before = checked_totals(&self.states, &grid, n, &mut cancelled)?;
        let mut next_grid = grid;
        let mut source_change = ForceChange3d::default();
        let accepted = stepping::advance::<Error, _>(
            maximum_duration.value(),
            self.time,
            rate,
            controls,
            &mut cancelled,
            |dt, reconstruction, cancelled| {
                let next = motion.grid_after(&grid, dt).map_err(|error| match error {
                    fs_gas::dynamics::Error::Invalid(reason) => Error::Invalid(reason),
                    _ => Error::Invalid("moving mixture geometry"),
                })?;
                let stats = if let Some((force, control, workspace)) = forcing.as_mut() {
                    copy_values(&self.states, &mut workspace.kicked, cancelled)?;
                    let mut source = forcing::Accumulator::default();
                    forcing::kick_field(
                        &self.gas,
                        &grid,
                        &mut workspace.kicked,
                        *force,
                        *control,
                        0.5 * dt,
                        &mut source,
                        false,
                        cancelled,
                    )?;
                    let mut stats = attempt_with_forcing(
                        &self.gas,
                        &grid,
                        &next,
                        self.boundaries,
                        &self.reservoirs,
                        motion,
                        &workspace.kicked,
                        &mut self.work,
                        dt,
                        controls.cfl,
                        reconstruction,
                        Some(*control),
                        cancelled,
                    )?;
                    forcing::kick_field(
                        &self.gas,
                        &next,
                        &mut self.work.candidate,
                        *force,
                        *control,
                        0.5 * dt,
                        &mut source,
                        false,
                        cancelled,
                    )?;
                    let final_rate = forcing::admit_moving_field(
                        &self.gas,
                        &next,
                        &self.work.candidate,
                        *control,
                        motion,
                        cancelled,
                    )?;
                    stepping::check_cfl::<Error>(final_rate * dt, controls.cfl)?;
                    stats.rate = stats.rate.max(final_rate);
                    source_change = source.finish()?;
                    stats
                } else {
                    attempt(
                        &self.gas,
                        &grid,
                        &next,
                        self.boundaries,
                        &self.reservoirs,
                        motion,
                        &self.states,
                        &mut self.work,
                        dt,
                        controls.cfl,
                        reconstruction,
                        cancelled,
                    )?
                };
                next_grid = next;
                Ok(stats)
            },
        )?;
        let after = checked_totals(&self.work.candidate, &next_grid, n, &mut cancelled)?;
        let mut report = reported_step(
            &self.gas,
            before,
            after,
            &self.work.boundary_trial,
            accepted.dt,
            accepted.stats,
            accepted.halvings,
            accepted.reconstruction != controls.reconstruction,
            &mut cancelled,
        )?;
        if forcing.is_some() {
            report.source_change = source_change;
            for axis in 0..3 {
                report.conservation_residual.momentum[axis] -= source_change.impulse[axis];
            }
            report.conservation_residual.energy -= source_change.work;
            report.absolute_energy_residual_j =
                report.conservation_residual.absolute_energy(&self.gas)?;
            if !report
                .conservation_residual
                .momentum
                .iter()
                .all(|v| v.is_finite())
                || !report.conservation_residual.energy.is_finite()
            {
                return Err(Error::Invalid("moving mixture forced balance overflow"));
            }
        }
        stepping::poll::<Error>(&mut cancelled)?;
        std::mem::swap(&mut self.states, &mut self.work.candidate);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        self.grid = next_grid;
        self.time += accepted.dt;
        self.last_step = Some(report);
        Ok(report)
    }

    pub(super) fn admit_motion(&self, motion: MeshMotion) -> Result<(), Error> {
        motion.admit_boundaries::<Error>(self.boundaries)?;
        if self
            .prepared_reservoir_profile_sides()
            .into_iter()
            .any(|x| x)
            || BoundarySide::ALL
                .into_iter()
                .any(|side| self.reservoir(side).is_some())
        {
            return Err(Error::Invalid(
                "moving mixture does not admit prescribed exterior states or profiles",
            ));
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn reported_step(
    gas: &MixtureGas,
    before: Totals,
    after: Totals,
    boundary_trial: &[f64],
    dt: f64,
    stats: stepping::StageStats,
    halvings: usize,
    first_order_retry: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StepReport, Error> {
    let n = gas.species_count();
    let width = n + 4;
    let outward = sum_rows(boundary_trial, width, 1.0, cancelled)?;
    let boundary_change = Totals::from_row(&outward.map(|v| -v), n);
    let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
        if d < width {
            after.component(d) - before.component(d) - boundary_change.component(d)
        } else {
            0.0
        }
    });
    if !residual.iter().all(|v| v.is_finite()) {
        return Err(Error::Invalid("moving mixture extensive balance overflow"));
    }
    let conservation_residual = Totals::from_row(&residual[..width], n);
    Ok(StepReport {
        elapsed: Time::new(dt),
        maximum_cfl: stats.rate * dt,
        halvings,
        hll_fallback_faces: stats.fallback,
        limited_cells: stats.limited,
        first_order_retry,
        before,
        after,
        boundary_change,
        source_change: ForceChange3d::default(),
        absolute_energy_residual_j: conservation_residual.absolute_energy(gas)?,
        conservation_residual,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &MixtureGas,
    original_grid: &CartesianGrid,
    next_grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    motion: MeshMotion,
    states: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    attempt_with_forcing(
        gas,
        original_grid,
        next_grid,
        boundaries,
        reservoirs,
        motion,
        states,
        work,
        dt,
        cfl,
        reconstruction,
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn attempt_with_forcing(
    gas: &MixtureGas,
    original_grid: &CartesianGrid,
    next_grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    motion: MeshMotion,
    states: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    let width = gas.species_count() + 4;
    stepping::rk2(cancelled, |final_stage, cancelled| {
        let grid = if final_stage {
            next_grid
        } else {
            original_grid
        };
        let input = if final_stage { &work.stage[..] } else { states };
        let mut stats = spatial_with_motion(
            gas,
            grid,
            boundaries,
            reservoirs,
            false,
            FiniteStage::Accepted,
            input,
            &mut work.primitive,
            &mut work.slopes,
            &mut work.flux,
            &mut work.metadata,
            reconstruction,
            force_control,
            None,
            Some(motion),
            cancelled,
        )?;
        stepping::check_cfl::<Error>(stats.rate * dt, cfl)?;
        let updated_rate = if final_stage {
            update(
                gas,
                original_grid,
                grid,
                next_grid,
                states,
                &work.stage,
                &work.flux,
                &mut work.candidate,
                dt,
                true,
                force_control,
                motion,
                cancelled,
            )?
        } else {
            update(
                gas,
                original_grid,
                grid,
                next_grid,
                states,
                states,
                &work.flux,
                &mut work.stage,
                dt,
                false,
                force_control,
                motion,
                cancelled,
            )?
        };
        if force_control.is_some() {
            stepping::check_cfl::<Error>(updated_rate * dt, cfl)?;
            stats.rate = stats.rate.max(updated_rate);
        }
        cartesian::capture::<Error>(
            grid,
            width,
            0.5 * dt,
            final_stage,
            |axis, face, d| work.flux[axis][face * width + d],
            |boundary, d, delta, add| {
                let out = &mut work.boundary_trial[boundary * width + d];
                *out = (if add { *out } else { 0.0 }) + delta;
            },
            cancelled,
        )?;
        Ok(stats)
    })
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &MixtureGas,
    original_grid: &CartesianGrid,
    stage_grid: &CartesianGrid,
    next_grid: &CartesianGrid,
    original: &[f64],
    stage: &[f64],
    flux: &[Vec<f64>; 3],
    output: &mut [f64],
    dt: f64,
    final_stage: bool,
    force_control: Option<ForceControl3d>,
    motion: MeshMotion,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let n = gas.species_count();
    let width = n + 4;
    let v0 = original_grid.cell_volume();
    let vs = stage_grid.cell_volume();
    let v1 = next_grid.cell_volume();
    let mut row = [0.0; WIDTH_LIMIT];
    let mut rate: f64 = 0.0;
    cartesian::divergence::<WIDTH_LIMIT, Error, _>(
        stage_grid,
        width,
        dt,
        |axis, face, d| flux[axis][face * width + d],
        cancelled,
        |k, delta, cancelled| {
            let start = k * width;
            for d in 0..width {
                let forward = vs * stage[start + d] + vs * delta[d];
                row[d] = if final_stage {
                    (0.5 * v0 * original[start + d] + 0.5 * forward) / v1
                } else {
                    forward / v1
                };
            }
            let p = cell(&row[..width], n).primitive(gas, &mut *cancelled)?;
            if let Some(control) = force_control {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
                rate = rate.max(forcing::moving_acoustic_rate(p, next_grid, k, motion)?);
            }
            output[start..start + width].copy_from_slice(&row[..width]);
            Ok(())
        },
    )?;
    Ok(rate)
}

#[cfg(test)]
mod tests;
