//! Resolved moving-wall pressure on the original evolving species field.

use super::evolve::{Stage, Workspace, operator_with_face_motion, totals};
use super::{
    Boundaries, Boundary, Error, MAX_COMPONENTS, Mesh1d, MixtureEuler1d, MixtureGas,
    Reconstruction, StepControl, StepReport, Totals, allocated, poll, positive,
};
use crate::body_force::{Accumulator, BodyForce, ForceChange, ForceControl};
use crate::dynamics::MovingWorkspace;
use fs_qty::Time;

impl MixtureEuler1d {
    /// Advance all resolved species between two prescribed moving end walls.
    ///
    /// Axial wall speeds [m/s] are held over the accepted step; interior faces
    /// follow the same affine motion and constant-area geometry as the frozen
    /// Euler owner. Each canonical species, laboratory momentum and referenced
    /// total energy evolve extensively as `Q=V*U`, with a shared RK2 volume
    /// trajectory and relative-face CFL at both stages. Zero speeds use the
    /// original stationary step and preserve its arithmetic.
    ///
    /// Current partial densities and source calorics supply every face pressure
    /// and wave. The existing compatible mass-fraction reconstruction and
    /// HLLC/HLL owner remain in charge of species transport. At moving walls
    /// every species flux is exactly zero, and the same-stage pressure impulse
    /// and `p_star*w` work appear in the returned boundary exchanges. Species
    /// energy references are transported with their own species fluxes; no
    /// clipping, reference replacement or composition normalization is added.
    ///
    /// Mesh, complete species field and clock publish together. Geometry scratch
    /// is fallibly allocated on first nonzero motion and reused thereafter.
    /// Cancellation from the shared geometry owner remains cancellation, never
    /// a retryable source refusal.
    ///
    /// This entry requires constant area and two `Reflecting` ends. It does not
    /// advance chemistry, molecular/wall transport, finite reservoirs or piston
    /// inertia. Their separately bound stationary integrators are unchanged.
    /// General moving unequal-caloric material contacts retain the existing
    /// conservative scheme's pressure-oscillation limitation; the moving frame
    /// does not establish exact pressure equilibrium or an entropy theorem.
    ///
    /// # Errors
    /// Refuses unsupported geometry or boundaries, invalid speeds/controls,
    /// unrepresentable motion, source-domain failures and exhausted retries
    /// without publishing partial physical state.
    pub fn step_with_moving_walls(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        self.advance_moving(maximum_duration, controls, wall_velocity, None, cancelled)
    }

    /// Compose distributed body work and moving-wall pressure on every species.
    ///
    /// Exact source half-steps surround the original extensive ALE SSP-RK2 map.
    /// Each source uses its actual mixture density and cell volume: accepted
    /// geometry first, admitted moved geometry second. All species retain their
    /// original transport and formation-energy references. The returned
    /// `source_change` contains actual impulse and conjugate work, independently
    /// of the original boundary pressure impulse and `p_star * wall_velocity`
    /// work. No source work is relabeled as heat or a balance correction.
    ///
    /// The original EOS, authored Mach bound and full-interval CFL relative to
    /// moving faces admit both source endpoints. Rejected trials restart both
    /// sources, transport and geometry from the accepted field. State, complete
    /// mesh, clock and accounting publish once, after every check. The existing
    /// flat body-source buffer is reused without cloning fields per retry.
    /// Zero force follows the original moving method; zero wall speeds follow
    /// the stationary forced method with its original arithmetic.
    ///
    /// Sources are held in original cell order for the accepted interval.
    /// Constant area and reflecting ends remain required. Chemistry, piston
    /// inertia, radial heat transport and changing source histories retain
    /// separate owners. This does not change the underlying scheme's admitted
    /// pressure-equilibrium or entropy qualifications.
    ///
    /// # Errors
    /// Invalid source/controls, geometry, source-domain or CFL refusal, exhausted
    /// retries and cancellation leave all accepted physical state unchanged.
    pub fn step_with_moving_walls_and_body_force(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        force: BodyForce<'_>,
        force_controls: ForceControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        force_controls.admit::<Error>()?;
        let source = force
            .admit::<Error>(self.cell_count(), &mut cancelled)?
            .then_some((force, force_controls));
        self.advance_moving(maximum_duration, controls, wall_velocity, source, cancelled)
    }

    fn advance_moving(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        source: Option<(BodyForce<'_>, ForceControl)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(Error::Invalid(
                "moving mixture walls require a constant-area closed interval",
            ));
        }
        if wall_velocity.iter().any(|speed| !speed.is_finite())
            || !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("moving mixture speeds or step controls"));
        }
        if wall_velocity == [0.0; 2] {
            return match source {
                Some((force, force_controls)) => self.step_with_body_force(
                    maximum_duration,
                    controls,
                    force,
                    force_controls,
                    cancelled,
                ),
                None => self.step(maximum_duration, controls, cancelled),
            };
        }
        if source.is_some() && self.body_force_stage.is_none() {
            self.body_force_stage = Some(allocated(self.states.len(), 0.0)?);
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh).map_err(geometry_error)?);
        }
        let moving = self
            .moving
            .as_mut()
            .expect("moving mixture workspace initialized");
        let (length, difference) = moving
            .prepare_velocities(&self.mesh, wall_velocity, &mut cancelled)
            .map_err(geometry_error)?;
        let n = self.gas.species_count();
        let m = n + 2;
        let area = self.mesh.face_areas()[0];
        let mut signal_rate = 0.0_f64;
        for (i, state) in self.states.chunks_exact(m).enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            let rate = moving.relative_signal_rate(&self.mesh, i, p.velocity, p.sound_speed);
            if !positive(rate) {
                return Err(Error::Invalid("moving mixture acoustic CFL overflow"));
            }
            signal_rate = signal_rate.max(rate);
        }
        let mut dt = maximum_duration
            .value()
            .min((0.9 * controls.cfl) / signal_rate);
        if difference < 0.0 {
            dt = dt.min((0.9 * controls.cfl) * (length / -difference));
        }
        let before = totals(&self.states, self.mesh.volumes(), n, &mut cancelled)?;
        let mut last = Error::Invalid("no admissible moving mixture increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(Error::Invalid(
                    "moving mixture step cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            let trial = if let Some((force, force_controls)) = source {
                attempt_with_body_force(
                    &self.gas,
                    &self.states,
                    &self.mesh,
                    self.boundaries,
                    dt,
                    controls.cfl,
                    reconstruction,
                    &mut self.work,
                    moving,
                    self.body_force_stage
                        .as_mut()
                        .expect("source workspace initialized"),
                    force,
                    force_controls,
                    &mut cancelled,
                )
            } else {
                attempt(
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
                )
                .map(|(a, b)| (a, b, ForceChange::default(), 0.0))
            };
            match trial {
                Ok((a, b, source_change, source_signal_rate)) => {
                    let after = totals(
                        &self.work.candidate,
                        moving.mesh.volumes(),
                        n,
                        &mut cancelled,
                    )?;
                    let mut report = reported_step_with_source(
                        before,
                        after,
                        [a, b],
                        dt,
                        area,
                        n,
                        source_change,
                        source_signal_rate,
                    )?;
                    report.halvings = halvings;
                    report.first_order_retry =
                        first_order_retry && controls.reconstruction == Reconstruction::Muscl;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    return Ok(report);
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

pub(super) fn reported_step(
    before: Totals,
    after: Totals,
    [a, b]: [Stage; 2],
    dt: f64,
    area: f64,
    n: usize,
) -> Result<StepReport, Error> {
    reported_step_with_source(
        before,
        after,
        [a, b],
        dt,
        area,
        n,
        ForceChange::default(),
        0.0,
    )
}

#[allow(clippy::too_many_arguments)]
fn reported_step_with_source(
    before: Totals,
    after: Totals,
    [a, b]: [Stage; 2],
    dt: f64,
    area: f64,
    n: usize,
    source_change: ForceChange,
    source_signal_rate: f64,
) -> Result<StepReport, Error> {
    let m = n + 2;
    let factor = 0.5 * dt * area;
    let mut left = [0.0; MAX_COMPONENTS];
    let mut right = [0.0; MAX_COMPONENTS];
    let mut boundary = [0.0; MAX_COMPONENTS];
    let mut residual = [0.0; MAX_COMPONENTS];
    for k in 0..m {
        left[k] = factor * (a.left[k] + b.left[k]);
        right[k] = -factor * (a.right[k] + b.right[k]);
        boundary[k] = left[k] + right[k];
        residual[k] = after.get(k) - before.get(k) - boundary[k];
        let source = if k == n {
            source_change.impulse
        } else if k == n + 1 {
            source_change.work
        } else {
            0.0
        };
        if source != 0.0 {
            residual[k] -= source;
        }
    }
    if !left[..m]
        .iter()
        .chain(&right[..m])
        .chain(&boundary[..m])
        .chain(&residual[..m])
        .all(|value| value.is_finite())
    {
        return Err(Error::Invalid("moving mixture extensive balance overflow"));
    }
    Ok(StepReport {
        elapsed: Time::new(dt),
        maximum_cfl: a.signal_rate.max(b.signal_rate).max(source_signal_rate) * dt,
        halvings: 0,
        hll_fallback_faces: a.fallbacks + b.fallbacks,
        limited_cells: a.limited + b.limited,
        first_order_retry: false,
        before,
        after,
        boundary_change: Totals::from_values(&boundary[..m], n),
        left_boundary_change: Totals::from_values(&left[..m], n),
        right_boundary_change: Totals::from_values(&right[..m], n),
        area_momentum_impulse: 0.0,
        source_change,
        conservation_residual: Totals::from_values(&residual[..m], n),
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &MixtureGas,
    states: &[f64],
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
        .prepare_motion(mesh, &moving.velocity, dt, cancelled)
        .map_err(geometry_error)?;
    let n = gas.species_count();
    let m = n + 2;
    let a = operator_with_face_motion(
        gas,
        states,
        mesh,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        &mut work.face_data,
        [None; 2],
        Some(&moving.velocity),
        cancelled,
    )?;
    check_cfl(a.signal_rate * dt, cfl)?;
    let area = mesh.face_areas()[0];
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        let before_volume = mesh.volumes()[i];
        let after_volume = moving.mesh.volumes()[i];
        let stage = &mut work.stage[i * m..(i + 1) * m];
        for k in 0..m {
            stage[k] = (before_volume * q[k]
                + dt * area * (work.faces[i * m + k] - work.faces[(i + 1) * m + k]))
                / after_volume;
        }
        if !gas.admissible(stage) {
            return Err(Error::Invalid(
                "moving mixture predictor outside source domain",
            ));
        }
    }
    let b = operator_with_face_motion(
        gas,
        &work.stage,
        &moving.mesh,
        boundaries,
        reconstruction,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.faces,
        &mut work.face_data,
        [None; 2],
        Some(&moving.velocity),
        cancelled,
    )?;
    check_cfl(b.signal_rate * dt, cfl)?;
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        let before_volume = mesh.volumes()[i];
        let after_volume = moving.mesh.volumes()[i];
        let candidate = &mut work.candidate[i * m..(i + 1) * m];
        for k in 0..m {
            candidate[k] = (0.5 * before_volume * q[k]
                + 0.5
                    * (after_volume * work.stage[i * m + k]
                        + dt * area * (work.faces[i * m + k] - work.faces[(i + 1) * m + k])))
                / after_volume;
        }
        gas.decode(candidate, Some(work.primitive[i].temperature), cancelled)?;
    }
    Ok((a, b))
}

#[allow(clippy::too_many_arguments)]
fn attempt_with_body_force(
    gas: &MixtureGas,
    states: &[f64],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    work: &mut Workspace,
    moving: &mut MovingWorkspace,
    forced: &mut [f64],
    force: BodyForce<'_>,
    controls: ForceControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage, ForceChange, f64), Error> {
    let m = gas.species_count() + 2;
    let mut exchange = Accumulator::default();
    let mut signal_rate = 0.0_f64;
    for (i, (q, next)) in states
        .chunks_exact(m)
        .zip(forced.chunks_exact_mut(m))
        .enumerate()
    {
        next.copy_from_slice(q);
        super::evolve::forcing::update(
            gas,
            next,
            &mut work.primitive[i],
            i,
            force,
            0.5 * dt,
            controls,
            mesh.volumes()[i],
            &mut exchange,
            cancelled,
        )?;
        let p = work.primitive[i];
        let rate = moving.relative_signal_rate(mesh, i, p.velocity, p.sound_speed);
        check_cfl(rate * dt, cfl)?;
        signal_rate = signal_rate.max(rate);
    }
    let (a, b) = attempt(
        gas,
        forced,
        mesh,
        boundaries,
        dt,
        cfl,
        reconstruction,
        work,
        moving,
        cancelled,
    )?;
    for (i, q) in work.candidate.chunks_exact_mut(m).enumerate() {
        super::evolve::forcing::update(
            gas,
            q,
            &mut work.primitive[i],
            i,
            force,
            0.5 * dt,
            controls,
            moving.mesh.volumes()[i],
            &mut exchange,
            cancelled,
        )?;
        let p = work.primitive[i];
        let rate = moving.relative_signal_rate(&moving.mesh, i, p.velocity, p.sound_speed);
        check_cfl(rate * dt, cfl)?;
        signal_rate = signal_rate.max(rate);
    }
    Ok((a, b, exchange.finish(), signal_rate))
}

fn check_cfl(observed: f64, maximum: f64) -> Result<(), Error> {
    if !observed.is_finite() || observed > maximum * (1.0 + 64.0 * f64::EPSILON) {
        Err(Error::CflExceeded { observed, maximum })
    } else {
        Ok(())
    }
}

pub(super) fn geometry_error(error: crate::dynamics::Error) -> Error {
    match error {
        crate::dynamics::Error::Cancelled => Error::Cancelled,
        crate::dynamics::Error::Allocation => Error::Allocation,
        other => Error::Geometry(other),
    }
}

#[cfg(test)]
mod tests;
