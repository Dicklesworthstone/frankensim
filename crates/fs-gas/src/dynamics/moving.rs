//! Spatial piston dynamics on the original source-owned Euler field.
//!
//! The affine, constant-area geometry and extensive RK2 state use the same
//! moving-face velocities. No finite-reservoir or reduced pressure law replaces
//! the original resolved wave and caloric owners.

use super::evolve::{Stage, Workspace, operator_with_face_motion, totals};
use super::{
    Boundaries, Boundary, Conserved, Error, Euler1d, FrozenGas, Mesh1d, Reconstruction,
    StepControl, StepReport, Totals, allocated, poll, positive,
};
use fs_qty::Time;

pub(crate) struct MovingWorkspace {
    pub(crate) mesh: Mesh1d,
    pub(crate) velocity: Vec<f64>,
}

impl MovingWorkspace {
    pub(crate) fn new(mesh: &Mesh1d) -> Result<Self, Error> {
        let area = mesh.uniform_cross_section().ok_or(Error::Invalid(
            "moving interval requires constant cross section",
        ))?;
        // Uniform construction intentionally retains the declared width: on
        // fine meshes, subtracting rounded face coordinates loses more bits
        // than the original width*area volume. Keep that same geometry owner.
        let declared_width = mesh.uniform_cell_width().map(|width| width.value());
        for (i, &volume) in mesh.volumes().iter().enumerate() {
            let geometric = area * declared_width.unwrap_or_else(|| mesh.width(i));
            if !positive(geometric)
                || (volume - geometric).abs() > 64.0 * f64::EPSILON * volume.max(geometric)
            {
                return Err(Error::Invalid(
                    "moving interval volumes must match constant-area geometry",
                ));
            }
        }
        Ok(Self {
            mesh: mesh.try_copy()?,
            velocity: allocated(mesh.cell_count() + 1, 0.0)?,
        })
    }

    /// Shared geometry-only preparation for frozen and evolving mixtures.
    pub(crate) fn prepare_velocities(
        &mut self,
        mesh: &Mesh1d,
        wall_velocity: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), Error> {
        let n = mesh.cell_count();
        let positions = mesh.positions();
        let length = positions[n] - positions[0];
        let difference = wall_velocity[1] - wall_velocity[0];
        if !difference.is_finite() {
            return Err(Error::Invalid("moving-wall velocity difference overflow"));
        }
        for (i, speed) in self.velocity.iter_mut().enumerate() {
            poll(cancelled)?;
            *speed = if i == 0 {
                wall_velocity[0]
            } else if i == n {
                wall_velocity[1]
            } else {
                wall_velocity[0] + ((positions[i] - positions[0]) / length) * difference
            };
            if !speed.is_finite() {
                return Err(Error::Invalid("affine face velocity overflow"));
            }
        }
        Ok((length, difference))
    }
}

impl Euler1d {
    /// Advance the resolved gas between two prescribed impermeable moving walls.
    ///
    /// Both speeds [m/s] use the increasing axial coordinate and are held fixed
    /// over the accepted step. The interior faces move by affine interpolation
    /// between these endpoints, preserving mesh connectivity and nonuniform cell
    /// fractions. The caller may change the speeds between accepted steps.
    /// A returned duration can be smaller than `maximum_duration` under the
    /// original relative-wave CFL, source-domain and bounded-retry controls.
    ///
    /// This evolves extensive mass, laboratory momentum and total energy on the
    /// changing volumes. MUSCL reconstruction and HLLC/HLL retain the original
    /// source calorics; every Riemann problem is solved in its face's frame and
    /// transformed back to the laboratory ALE flux. Each end uses the original
    /// one-sided moving-wall shock/rarefaction pressure. No mass crosses either
    /// wall, while pressure impulse and work enter the same two-stage boundary
    /// exchanges in [`StepReport`]. In particular `boundary_change.energy` is
    /// the accepted wall work, and `wall_heat` remains zero.
    ///
    /// The RK2 state is advanced as `Q = volume * U` together with the exact
    /// linear volume trajectory `dV/dt = area * (w_right-w_left)`. The report's
    /// CFL is relative to moving faces and checked on both stage geometries.
    /// State, geometry and clock publish together after all endpoint checks.
    /// Geometry scratch is allocated fallibly on the first nonzero-motion call
    /// and reused thereafter. Zero speeds use the existing stationary step.
    ///
    /// This bounded lane requires constant area and two `Reflecting` ends.
    /// It refuses periodic/open/reservoir ends and bound gravitational reference
    /// fields. Piston inertia, wall heat, evolving mixture composition, radial
    /// deformation and multidimensional moving geometry remain separate owners.
    /// Boundary cells keep the existing unreconstructed trace, so smooth global
    /// accuracy near the walls is not asserted to be second order.
    ///
    /// # Errors
    /// Refuses unsupported geometry/boundaries, nonfinite speeds, invalid time
    /// controls, crossed/unrepresentable faces, source limits and exhausted
    /// retries without changing accepted cells, mesh or physical time.
    pub fn step_with_moving_walls(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.hydrostatic.is_some()
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(Error::Invalid(
                "moving walls require a constant-area closed interval without a gravitational reference",
            ));
        }
        if wall_velocity.iter().any(|speed| !speed.is_finite())
            || !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("moving-wall speeds or step controls"));
        }
        if wall_velocity == [0.0; 2] {
            return self.step(maximum_duration, controls, cancelled);
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh)?);
        }
        let moving = self.moving.as_mut().expect("moving workspace initialized");
        let (length, difference) =
            moving.prepare_velocities(&self.mesh, wall_velocity, &mut cancelled)?;
        let area = self.mesh.face_areas()[0];
        let mut signal_rate = 0.0_f64;
        for (i, &state) in self.states.iter().enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            let relative = (p.velocity - moving.velocity[i])
                .abs()
                .max((p.velocity - moving.velocity[i + 1]).abs());
            let rate = (relative + p.sound_speed) * (area / self.mesh.volumes()[i]);
            if !positive(rate) {
                return Err(Error::Invalid("moving interval acoustic CFL overflow"));
            }
            signal_rate = signal_rate.max(rate);
        }
        let mut dt = maximum_duration
            .value()
            .min((0.9 * controls.cfl) / signal_rate);
        if difference < 0.0 {
            // Retain a positive trial volume before either stage is attempted.
            // Both actual stage face-speed CFLs still determine acceptance.
            dt = dt.min((0.9 * controls.cfl) * (length / -difference));
        }
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        let mut last = Error::Invalid("no admissible moving time increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(Error::Invalid("moving step cannot advance the clock"));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            match attempt(
                &self.gas,
                &self.states,
                &self.mesh,
                self.boundaries,
                dt,
                controls.cfl,
                reconstruction,
                &mut self.work,
                moving,
                &mut cancelled,
            ) {
                Ok((a, b)) => {
                    let after =
                        totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                    let factor = 0.5 * dt * area;
                    let left = Totals::from_array(core::array::from_fn(|k| {
                        factor * (a.left.array()[k] + b.left.array()[k])
                    }));
                    let right = Totals::from_array(core::array::from_fn(|k| {
                        -factor * (a.right.array()[k] + b.right.array()[k])
                    }));
                    let boundary = Totals::from_array(core::array::from_fn(|k| {
                        left.array()[k] + right.array()[k]
                    }));
                    let residual = Totals::from_array(core::array::from_fn(|k| {
                        after.array()[k] - before.array()[k] - boundary.array()[k]
                    }));
                    if !after
                        .array()
                        .iter()
                        .chain(&left.array())
                        .chain(&right.array())
                        .chain(&boundary.array())
                        .chain(&residual.array())
                        .all(|value| value.is_finite())
                    {
                        return Err(Error::Invalid("moving extensive balance overflow"));
                    }
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    return Ok(StepReport {
                        elapsed: Time::new(dt),
                        maximum_cfl: a.signal_rate.max(b.signal_rate) * dt,
                        halvings,
                        hll_fallback_faces: a.fallbacks + b.fallbacks,
                        limited_cells: a.limited + b.limited,
                        first_order_retry: first_order_retry
                            && controls.reconstruction == Reconstruction::Muscl,
                        before,
                        after,
                        boundary_change: boundary,
                        left_boundary_change: left,
                        right_boundary_change: right,
                        area_momentum_impulse: 0.0,
                        gravity_momentum_impulse: 0.0,
                        gravity_work: 0.0,
                        wall_momentum_impulse: 0.0,
                        wall_heat: 0.0,
                        conservation_residual: residual,
                    });
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => {
                    if !matches!(error, Error::CflExceeded { .. }) {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(Error::StepRetries {
            attempts: controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &FrozenGas,
    states: &[Conserved],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    work: &mut Workspace,
    moving: &mut MovingWorkspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage), Error> {
    moving
        .mesh
        .prepare_motion(mesh, &moving.velocity, dt, cancelled)?;
    let a = operator_with_face_motion(
        gas,
        states,
        mesh,
        boundaries,
        reconstruction,
        None,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        [None; 2],
        Some(&moving.velocity),
        cancelled,
    )?;
    check_cfl(a.signal_rate, dt, cfl)?;
    let area = mesh.face_areas()[0];
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        let before_volume = mesh.volumes()[i];
        let after_volume = moving.mesh.volumes()[i];
        work.stage[i] = Conserved::from_array(core::array::from_fn(|k| {
            (before_volume * q.array()[k]
                + dt * area * (work.faces[i].array()[k] - work.faces[i + 1].array()[k]))
                / after_volume
        }));
        if !gas.admissible_state(work.stage[i]) {
            return Err(Error::Invalid("moving predictor outside source domain"));
        }
    }
    let b = operator_with_face_motion(
        gas,
        &work.stage,
        &moving.mesh,
        boundaries,
        reconstruction,
        None,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        [None; 2],
        Some(&moving.velocity),
        cancelled,
    )?;
    check_cfl(b.signal_rate, dt, cfl)?;
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        let before_volume = mesh.volumes()[i];
        let after_volume = moving.mesh.volumes()[i];
        let candidate = Conserved::from_array(core::array::from_fn(|k| {
            (0.5 * before_volume * q.array()[k]
                + 0.5
                    * (after_volume * work.stage[i].array()[k]
                        + dt * area * (work.faces[i].array()[k] - work.faces[i + 1].array()[k])))
                / after_volume
        }));
        gas.decode(candidate, Some(work.primitive[i].temperature), cancelled)?;
        work.candidate[i] = candidate;
    }
    Ok((a, b))
}

fn check_cfl(rate: f64, dt: f64, cfl: f64) -> Result<(), Error> {
    if rate * dt > cfl * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded {
            observed: rate * dt,
            maximum: cfl,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
