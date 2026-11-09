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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        if motion == MeshMotion::default() {
            return self.step(maximum_duration, controls, cancelled);
        }
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
                let stats = attempt(
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
                )?;
                next_grid = next;
                Ok(stats)
            },
        )?;
        let after = checked_totals(&self.work.candidate, &next_grid, n, &mut cancelled)?;
        let report = reported_step(
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
    let width = gas.species_count() + 4;
    stepping::rk2(cancelled, |final_stage, cancelled| {
        let grid = if final_stage {
            next_grid
        } else {
            original_grid
        };
        let input = if final_stage { &work.stage[..] } else { states };
        let stats = spatial_with_motion(
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
            None,
            None,
            Some(motion),
            cancelled,
        )?;
        stepping::check_cfl::<Error>(stats.rate * dt, cfl)?;
        if final_stage {
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
                cancelled,
            )?;
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
                cancelled,
            )?;
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
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let n = gas.species_count();
    let width = n + 4;
    let v0 = original_grid.cell_volume();
    let vs = stage_grid.cell_volume();
    let v1 = next_grid.cell_volume();
    let mut row = [0.0; WIDTH_LIMIT];
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
            cell(&row[..width], n).primitive(gas, &mut *cancelled)?;
            output[start..start + width].copy_from_slice(&row[..width]);
            Ok(())
        },
    )
}

#[cfg(test)]
mod tests;
