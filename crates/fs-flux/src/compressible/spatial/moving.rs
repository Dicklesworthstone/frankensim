//! Conservative translation and piston motion of the original volumetric field.
//!
//! One axis may change width while the box translates in all three directions.
//! Cell volumes are then affine in time. SSP RK2 advances extensive V*U using
//! the actual stage face areas and the same stage volumes; its geometric
//! conservation law is exact up to floating-point arithmetic. General
//! simultaneous multi-axis deformation needs a different geometric quadrature.

use super::super::cartesian::advance;
use super::*;

/// Held laboratory velocities of the six Cartesian bounding planes [m/s].
///
/// Each component is the velocity along the positive coordinate axis, even on
/// a minimum side. Interior planes interpolate their two endpoints. At most
/// one axis may deform, and periodic axes must translate without changing
/// width. Reflecting boundaries move with their planes and transfer pressure
/// impulse and work; transmissive boundaries retain zero-gradient gas traces.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MeshMotion {
    /// Velocities of the minimum x/y/z planes [m/s].
    pub lower_velocity: [f64; 3],
    /// Velocities of the maximum x/y/z planes [m/s].
    pub upper_velocity: [f64; 3],
}

impl MeshMotion {
    pub(super) fn admit(self, flow: &Euler3d) -> Result<(), Error> {
        self.admit_boundaries::<Error>(flow.boundaries)?;
        if flow.hydrostatic.is_some()
            || flow
                .prepared_reservoir_profile_sides()
                .into_iter()
                .any(|x| x)
            || BoundarySide::ALL
                .into_iter()
                .any(|side| flow.reservoir(side).is_some() || flow.pressure_outlet(side).is_some())
        {
            return Err(Error::Invalid(
                "moving box does not admit retained gravity or prescribed exterior profiles",
            ));
        }
        Ok(())
    }

    /// Shared kinematic and topological admission; each field owner separately
    /// admits its retained source and exterior data.
    pub(super) fn admit_boundaries<E: super::super::cartesian::Failure>(
        self,
        boundaries: Boundaries3d,
    ) -> Result<(), E> {
        let rates =
            std::array::from_fn::<_, 3, _>(|a| self.upper_velocity[a] - self.lower_velocity[a]);
        if !self
            .lower_velocity
            .into_iter()
            .chain(self.upper_velocity)
            .chain(rates)
            .all(f64::is_finite)
            || rates.into_iter().filter(|x| *x != 0.0).count() > 1
        {
            return Err(E::invalid(
                "moving box needs finite translation and at most one deforming axis",
            ));
        }
        for (axis, rate) in rates.into_iter().enumerate() {
            match boundaries.axis(axis) {
                AxisBoundaries::Periodic if rate != 0.0 => {
                    return Err(E::invalid(
                        "periodic moving planes must have identical velocity",
                    ));
                }
                AxisBoundaries::Ends { lower, upper } => {
                    if [lower, upper].into_iter().any(|condition| {
                        !matches!(
                            condition,
                            BoundaryCondition::Reflecting | BoundaryCondition::Transmissive
                        )
                    }) {
                        return Err(E::invalid(
                            "moving box requires wall or zero-gradient boundary traces",
                        ));
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    pub(super) fn velocity(self, grid: &CartesianGrid, axis: usize, plane: usize) -> f64 {
        if plane == 0 {
            self.lower_velocity[axis]
        } else if plane == grid.shape[axis] {
            self.upper_velocity[axis]
        } else {
            self.lower_velocity[axis]
                + (self.upper_velocity[axis] - self.lower_velocity[axis])
                    * (plane as f64 / grid.shape[axis] as f64)
        }
    }

    pub(super) fn face_velocity(self, grid: &CartesianGrid, axis: usize, index: usize) -> f64 {
        let mut shape = grid.shape;
        shape[axis] += 1;
        let plane = match axis {
            0 => index % shape[0],
            1 => (index / shape[0]) % shape[1],
            _ => index / (shape[0] * shape[1]),
        };
        self.velocity(grid, axis, plane)
    }

    pub(super) fn grid_after(self, grid: &CartesianGrid, dt: f64) -> Result<CartesianGrid, Error> {
        let spacing = std::array::from_fn::<_, 3, _>(|a| {
            grid.spacing[a]
                + dt * ((self.upper_velocity[a] - self.lower_velocity[a]) / grid.shape[a] as f64)
        });
        let mut next = CartesianGrid::new(
            grid.shape[0],
            grid.shape[1],
            grid.shape[2],
            spacing[0],
            spacing[1],
            spacing[2],
        )?;
        next.origin = std::array::from_fn(|a| grid.origin[a] + dt * self.lower_velocity[a]);
        if (0..3).any(|a| {
            !next.origin[a].is_finite()
                || (self.lower_velocity[a] != 0.0 && next.origin[a] == grid.origin[a])
                || !(next.origin[a] + next.shape[a] as f64 * next.spacing[a]).is_finite()
                || next.origin[a] + next.spacing[a] <= next.origin[a]
                || (self.upper_velocity[a] != self.lower_velocity[a]
                    && next.spacing[a] == grid.spacing[a])
        }) {
            return Err(Error::Invalid(
                "unrepresentable moving box coordinates or deformation",
            ));
        }
        // Every width is affine. Positive old and new widths exclude collapse
        // everywhere in the interval, not merely at the quadrature endpoints.
        Ok(next)
    }
}

impl Euler3d {
    /// Advance this same gas field on a translating or single-axis deforming box.
    ///
    /// The held plane speeds define actual geometry throughout the accepted
    /// duration. Each original face carries one laboratory ALE flux; adjacent
    /// cells use it with opposite signs. Extensive mass, all momenta and energy
    /// are advanced before division by the actual new volume. Moving reflecting
    /// walls exchange zero mass and the original nonlinear pressure work.
    ///
    /// Both SSP stages use relative face-wave CFL admission and the source's
    /// caloric domain. Failure or cancellation preserves cells, geometry, time,
    /// loads and the previous report. A zero motion uses [`Self::step`] exactly.
    /// The original field/face scratch is reused without step-time allocation.
    ///
    /// This entry admits periodic translation, moving slip walls and unbound
    /// zero-gradient openings. Prescribed exterior states/profiles, pressure
    /// outlets, retained hydrostatics, simultaneous deformation of multiple
    /// axes, and molecular/FSI composition require their own moving-domain laws.
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
        poll(&mut cancelled)?;
        motion.admit(self)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("volumetric duration, CFL or retry controls"));
        }
        let grid = self.grid;
        let mut rate: f64 = 0.0;
        for (cell, (p, q)) in self.work.primitive.iter_mut().zip(&self.states).enumerate() {
            *p = q.primitive(&self.gas, &mut cancelled)?;
            let c = grid.coordinates(cell);
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
                return Err(Error::Invalid("moving box initial relative CFL overflow"));
            }
            rate = rate.max(local);
        }
        let before = checked_totals(&self.states, grid.cell_volume(), &mut cancelled)?;
        let mut next_grid = grid;
        let accepted = advance::<Error, _>(
            maximum_duration.value(),
            self.time,
            rate,
            controls,
            &mut cancelled,
            |dt, reconstruction, cancelled| {
                let next = motion.grid_after(&grid, dt)?;
                let stats = attempt(
                    &self.gas,
                    &grid,
                    &next,
                    self.boundaries,
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
        let after = checked_totals(
            &self.work.candidate,
            next_grid.cell_volume(),
            &mut cancelled,
        )?;
        let report = reported_step(
            before,
            after,
            &self.work.boundary_trial,
            accepted.dt,
            accepted.stats,
            accepted.halvings,
            accepted.reconstruction != controls.reconstruction,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        std::mem::swap(&mut self.states, &mut self.work.candidate);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        self.grid = next_grid;
        self.time += accepted.dt;
        self.last_step = Some(report);
        Ok(report)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn reported_step(
    before: Totals,
    after: Totals,
    boundary_trial: &[BoundaryImpulse],
    dt: f64,
    stats: StageStats,
    halvings: usize,
    first_order_retry: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StepReport, Error> {
    let outward = sum_arrays(boundary_trial.iter().map(|x| x.array()), cancelled)?;
    let boundary_change = Totals::from_array(outward.map(|x| -x));
    let a = after.array();
    let b = before.array();
    let c = boundary_change.array();
    let residual = std::array::from_fn(|k| a[k] - b[k] - c[k]);
    if !a
        .into_iter()
        .chain(b)
        .chain(c)
        .chain(residual)
        .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "unrepresentable moving-volume extensive balance",
        ));
    }
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
        conservation_residual: Totals::from_array(residual),
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &FrozenGas,
    original_grid: &CartesianGrid,
    next_grid: &CartesianGrid,
    boundaries: Boundaries3d,
    motion: MeshMotion,
    states: &[State],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    super::super::cartesian::rk2(cancelled, |final_stage, cancelled| {
        let grid = if final_stage {
            next_grid
        } else {
            original_grid
        };
        let input = if final_stage { &work.stage[..] } else { states };
        let stats = spatial(
            gas,
            grid,
            boundaries,
            motion,
            input,
            &mut work.primitive,
            &mut work.slopes,
            &mut work.flux,
            reconstruction,
            cancelled,
        )?;
        check_cfl::<Error>(stats.rate * dt, cfl)?;
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
        evolve::capture_boundaries(
            grid,
            &work.flux,
            0.5 * dt,
            &mut work.boundary_trial,
            final_stage,
            cancelled,
        )?;
        Ok(stats)
    })
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    motion: MeshMotion,
    states: &[State],
    primitive: &mut [Primitive],
    slopes: &mut [[[f64; 5]; 3]],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    for (p, q) in primitive.iter_mut().zip(states) {
        *p = q.primitive(gas, &mut *cancelled)?;
    }
    let trace = |cell: usize, sign: f64, slope: &[f64]| {
        let p = primitive[cell].array();
        std::array::from_fn::<_, 5, _>(|d| p[d] + sign * slope[d])
    };
    let limited = cartesian::slopes_with_wall_velocity::<Error>(
        grid,
        boundaries,
        5,
        [1, 2, 3],
        reconstruction,
        std::array::from_fn(|a| [motion.lower_velocity[a], motion.upper_velocity[a]]),
        |cell, d| primitive[cell].array()[d],
        |_, _, _, _| Ok(0),
        slopes.as_flattened_mut().as_flattened_mut(),
        |cell, _, sign, slope| {
            let p = trace(cell, sign, slope);
            let temperature = p[4] / (gas.gas_constant() * p[0]);
            p.into_iter().all(f64::is_finite)
                && positive(p[0])
                && positive(p[4])
                && temperature >= gas.domain().minimum_temperature
                && temperature <= gas.domain().maximum_temperature
        },
        cancelled,
    )?;
    let mut stats = StageStats {
        limited,
        ..StageStats::default()
    };
    cartesian::faces::<Error, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            fluxes[axis][face.index] = fluxes[axis][first];
            return Ok(());
        }
        let [a, b] = tangents(axis);
        let reconstruct = |cell: usize, sign: f64| {
            let p = trace(cell, sign, &slopes[cell][axis]);
            SpatialPrimitive {
                density: p[0],
                normal_velocity: p[axis + 1],
                transverse_velocity: [p[a + 1], p[b + 1]],
                pressure: p[4],
            }
        };
        let speed = motion.face_velocity(grid, axis, face.index);
        let flux = match face.kind {
            cartesian::FaceKind::Interior { left, right } => gas
                .moving_hllc_spatial_flux_from_primitive(
                    reconstruct(left, 0.5),
                    reconstruct(right, -0.5),
                    speed,
                    &mut *cancelled,
                )?,
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition,
            } => {
                let p = reconstruct(cell, 0.5 * sign);
                match condition {
                    BoundaryCondition::Reflecting => gas.moving_wall_spatial_flux_from_primitive(
                        p,
                        sign,
                        speed,
                        &mut *cancelled,
                    )?,
                    BoundaryCondition::Transmissive => {
                        gas.moving_hllc_spatial_flux_from_primitive(p, p, speed, &mut *cancelled)?
                    }
                    _ => unreachable!("moving-boundary conditions admitted before the trial"),
                }
            }
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        stats.fallback += usize::from(flux.used_hll_fallback);
        fluxes[axis][face.index] = flux;
        Ok(())
    })?;
    stats.rate = cartesian::rate::<Error>(
        grid,
        |axis, face| fluxes[axis][face].signal_speed,
        cancelled,
    )?;
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &FrozenGas,
    original_grid: &CartesianGrid,
    stage_grid: &CartesianGrid,
    next_grid: &CartesianGrid,
    original: &[State],
    stage: &[State],
    fluxes: &[Vec<SpatialFaceFlux>; 3],
    output: &mut [State],
    dt: f64,
    final_stage: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let v0 = original_grid.cell_volume();
    let vs = stage_grid.cell_volume();
    let v1 = next_grid.cell_volume();
    cartesian::divergence::<5, Error, _>(
        stage_grid,
        5,
        dt,
        |axis, face, d| evolve::flux_array(fluxes[axis][face], axis)[d],
        cancelled,
        |cell, delta, cancelled| {
            let q0 = original[cell].array();
            let qs = stage[cell].array();
            let values = std::array::from_fn(|d| {
                let forward = vs * qs[d] + vs * delta[d];
                if final_stage {
                    (0.5 * v0 * q0[d] + 0.5 * forward) / v1
                } else {
                    forward / v1
                }
            });
            let candidate = State::from_array(values);
            candidate.primitive(gas, &mut *cancelled)?;
            output[cell] = candidate;
            Ok(())
        },
    )
}

#[cfg(test)]
mod tests;
