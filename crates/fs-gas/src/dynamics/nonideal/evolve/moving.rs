//! Prescribed impermeable moving boundaries on the original real-gas field.

use super::*;
use crate::dynamics::MovingWorkspace;

mod flux;
mod kinematics;
pub(super) use flux::{interface, wall_flux};

impl RealGasEuler1d {
    /// Advance the resolved original real gas between two prescribed end walls.
    ///
    /// Both velocities [m/s] use the increasing axial coordinate and are held
    /// over this accepted interval. Interior faces move by affine interpolation
    /// through the original nonuniform mesh. The extensive SSP-RK2 state and
    /// native geometric-conservation law use those SAME face velocities.
    ///
    /// Reflecting ends retain the original source-checked mirrored Riemann
    /// traction; PhysicalReflecting ends retain their original shock/isentrope
    /// controls and physical pressure. Mass transfer is exactly zero and energy
    /// transfer is the actual pressure times wall velocity. Boundary impulse
    /// and work use the same stage fluxes in StepReport.
    ///
    /// Both stages and the second Euler endpoint resolve through the original
    /// NASA/Peng--Robinson branch. CFL is relative to moving faces and is
    /// checked on both geometries. Retry starts from accepted cells and mesh;
    /// cells, cached primitives, geometry and time publish together.
    ///
    /// This entry requires constant area and two impermeable ends. It supports
    /// first-order or original-source MUSCL traces, retaining unreconstructed
    /// end cells. Zero velocities delegate to the stationary method. Geometry
    /// scratch is allocated once and reused. It does not infer piston inertia,
    /// wall heat, radial deformation, a phase flash or a global accuracy theorem.
    pub fn step_with_moving_walls(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        self.advance_moving(maximum_duration, controls, wall_velocity, None, cancelled)
    }

    /// Compose held mechanical-force half steps with moving real-gas transport.
    ///
    /// The first kick uses the accepted density and physical volume. The
    /// second uses the actual transported density and moved volume. Both use
    /// original calorics, laboratory kinetic work, the caller's Mach limit and
    /// full-interval relative-face CFL. Pressure work stays in boundary_change;
    /// external body impulse/work stays in source_change, counted once.
    ///
    /// Force values follow the original cell order for this whole interval.
    /// Refusal or cancellation publishes neither source nor motion. Zero force
    /// uses the source-free path and zero motion uses the existing forced step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_moving_walls_and_body_force(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        force: BodyForce<'_>,
        force_controls: ForceControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        poll(&mut cancelled)?;
        force_controls.admit::<RealGasFlowError>()?;
        let source = force
            .admit::<RealGasFlowError>(self.states.len(), &mut cancelled)?
            .then_some((force, force_controls));
        self.advance_moving(maximum_duration, controls, wall_velocity, source, cancelled)
    }

    #[allow(clippy::too_many_lines)]
    fn advance_moving(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        wall_velocity: [f64; 2],
        source: Option<(BodyForce<'_>, ForceControl)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        poll(&mut cancelled)?;
        let impermeable = |boundary| {
            matches!(
                boundary,
                RealGasBoundary::Reflecting | RealGasBoundary::PhysicalReflecting(_)
            )
        };
        if !matches!(self.boundaries, RealGasBoundaries::Ends { left, right }
            if impermeable(left) && impermeable(right))
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(RealGasFlowError::Invalid(
                "moving real gas requires constant area and two impermeable ends",
            ));
        }
        if wall_velocity.iter().any(|speed| !speed.is_finite())
            || !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(RealGasFlowError::Invalid(
                "finite moving-wall velocities, positive duration, CFL <= 0.45 and bounded retries required",
            ));
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
        if controls.reconstruction == Reconstruction::Muscl && self.work.traces.is_none() {
            self.work.traces = Some(allocated(
                self.states.len(),
                [(Conserved::default(), Primitive::default()); 2],
            )?);
        }
        if source.is_some() && self.force_work.is_none() {
            self.force_work = Some(forcing::Workspace::new(self.states.len())?);
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh)?);
        }
        let moving = self.moving.as_mut().expect("admitted moving workspace");
        let (length, difference) =
            moving.prepare_velocities(&self.mesh, wall_velocity, &mut cancelled)?;
        let mut signal_rate = 0.0_f64;
        for (i, p) in self.primitive.iter().enumerate() {
            poll(&mut cancelled)?;
            let local = moving.relative_signal_rate(&self.mesh, i, p.velocity, p.sound_speed);
            if !positive(local) {
                return Err(RealGasFlowError::Invalid(
                    "unrepresentable real-gas relative wave rate",
                ));
            }
            signal_rate = signal_rate.max(local);
        }
        let mut dt = maximum_duration
            .value()
            .min(0.9 * controls.cfl / signal_rate);
        if difference < 0.0 {
            dt = dt.min((0.9 * controls.cfl) * (length / -difference));
        }
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        if before.array().iter().any(|value| !value.is_finite()) {
            return Err(RealGasFlowError::Invalid(
                "nonfinite moving real-gas inventory",
            ));
        }
        let area = self.mesh.face_areas()[0];
        let mut last = RealGasFlowError::Invalid("no admissible moving real-gas interval");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || !(self.time + dt).is_finite() || self.time + dt <= self.time {
                return Err(RealGasFlowError::Invalid(
                    "moving real-gas interval cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            let trial = (|| {
                let mut exchange = Accumulator::default();
                let mut source_rate = 0.0_f64;
                let (states, primitive) = if let Some((force, force_controls)) = source {
                    let forced = self.force_work.as_mut().expect("admitted source workspace");
                    for (i, &q) in self.states.iter().enumerate() {
                        (forced.states[i], forced.primitive[i]) = forcing::update(
                            &self.gas,
                            q,
                            self.primitive[i],
                            force,
                            force_controls,
                            i,
                            0.5 * dt,
                            self.mesh.volumes()[i],
                            &mut exchange,
                            &mut cancelled,
                        )?;
                        source_cfl(
                            forced.primitive[i],
                            &self.mesh,
                            moving,
                            i,
                            dt,
                            controls.cfl,
                            &mut source_rate,
                        )?;
                    }
                    (&forced.states[..], &forced.primitive[..])
                } else {
                    (&self.states[..], &self.primitive[..])
                };
                let stages = attempt(
                    &self.gas,
                    states,
                    primitive,
                    &self.mesh,
                    self.boundaries,
                    self.boundary_primitive,
                    reconstruction,
                    dt,
                    controls.cfl,
                    &mut self.work,
                    moving,
                    &mut cancelled,
                )?;
                if let Some((force, force_controls)) = source {
                    for (i, q) in self.work.candidate.iter_mut().enumerate() {
                        (*q, self.work.candidate_primitive[i]) = forcing::update(
                            &self.gas,
                            *q,
                            self.work.candidate_primitive[i],
                            force,
                            force_controls,
                            i,
                            0.5 * dt,
                            moving.mesh.volumes()[i],
                            &mut exchange,
                            &mut cancelled,
                        )?;
                        source_cfl(
                            self.work.candidate_primitive[i],
                            &moving.mesh,
                            moving,
                            i,
                            dt,
                            controls.cfl,
                            &mut source_rate,
                        )?;
                    }
                }
                Ok((stages, exchange.finish(), source_rate))
            })();
            match trial {
                Ok(([a, b], source_change, source_rate)) => {
                    let after =
                        totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                    let mut report = crate::dynamics::moving::reported_step_with_source(
                        before,
                        after,
                        [a, b],
                        dt,
                        area,
                        source_change,
                        source_rate,
                    )?;
                    report.halvings = halvings;
                    report.first_order_retry =
                        first_order_retry && controls.reconstruction == Reconstruction::Muscl;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.primitive, &mut self.work.candidate_primitive);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    return Ok(report);
                }
                Err(RealGasFlowError::Cancelled) => return Err(RealGasFlowError::Cancelled),
                Err(error) => {
                    if !matches!(
                        error,
                        RealGasFlowError::Numerical(Error::CflExceeded { .. })
                    ) {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(RealGasFlowError::StepRetries {
            attempts: controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &RealGas,
    states: &[Conserved],
    primitive: &[Primitive],
    mesh: &Mesh1d,
    boundaries: RealGasBoundaries,
    boundary_primitive: [Option<Primitive>; 2],
    reconstruction: Reconstruction,
    dt: f64,
    cfl: f64,
    work: &mut Workspace,
    moving: &mut MovingWorkspace,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Stage; 2], RealGasFlowError> {
    moving
        .mesh
        .prepare_motion(mesh, &moving.velocity, dt, cancelled)?;
    let a = super::operator(
        gas,
        states,
        primitive,
        mesh,
        boundaries,
        boundary_primitive,
        Some(&moving.velocity),
        &mut work.faces,
        work.traces
            .as_deref_mut()
            .filter(|_| reconstruction == Reconstruction::Muscl),
        cancelled,
    )?;
    cfl_gate(a.signal_rate * dt, cfl)?;
    let dt_area = dt * mesh.face_areas()[0];
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        let next = Conserved::from_array(kinematics::predictor(
            q.array(),
            mesh.volumes()[i],
            moving.mesh.volumes()[i],
            dt_area,
            [work.faces[i].array(), work.faces[i + 1].array()],
        ));
        work.stage_primitive[i] = if next == q {
            primitive[i]
        } else {
            gas.decode(next, cancelled)?
        };
        work.stage[i] = next;
    }
    let b = super::operator(
        gas,
        &work.stage,
        &work.stage_primitive,
        &moving.mesh,
        boundaries,
        boundary_primitive,
        Some(&moving.velocity),
        &mut work.faces,
        work.traces
            .as_deref_mut()
            .filter(|_| reconstruction == Reconstruction::Muscl),
        cancelled,
    )?;
    cfl_gate(b.signal_rate * dt, cfl)?;
    for (i, &q) in states.iter().enumerate() {
        poll(cancelled)?;
        // The second Euler endpoint belongs to the twice-displaced trial
        // volume. It is source-admitted before the extensive SSP average,
        // matching the original real-gas stage discipline.
        let volume = moving.mesh.volumes()[i];
        let twice = volume + dt_area * (moving.velocity[i + 1] - moving.velocity[i]);
        if !positive(twice) {
            return Err(RealGasFlowError::Invalid(
                "nonpositive second moving Euler trial volume",
            ));
        }
        let fluxes = [work.faces[i].array(), work.faces[i + 1].array()];
        let second = Conserved::from_array(kinematics::predictor(
            work.stage[i].array(),
            volume,
            twice,
            dt_area,
            fluxes,
        ));
        let second_primitive = if second == work.stage[i] {
            work.stage_primitive[i]
        } else {
            gas.decode(second, cancelled)?
        };
        let candidate = Conserved::from_array(kinematics::corrector(
            q.array(),
            work.stage[i].array(),
            mesh.volumes()[i],
            volume,
            dt_area,
            fluxes,
        ));
        work.candidate_primitive[i] = if candidate == q {
            primitive[i]
        } else if candidate == second {
            second_primitive
        } else {
            gas.decode(candidate, cancelled)?
        };
        work.candidate[i] = candidate;
    }
    Ok([a, b])
}

#[allow(clippy::too_many_arguments)]
fn source_cfl(
    primitive: Primitive,
    mesh: &Mesh1d,
    moving: &MovingWorkspace,
    cell: usize,
    dt: f64,
    maximum: f64,
    rate: &mut f64,
) -> Result<(), RealGasFlowError> {
    let local = moving.relative_signal_rate(mesh, cell, primitive.velocity, primitive.sound_speed);
    if !positive(local) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable moving source CFL",
        ));
    }
    cfl_gate(local * dt, maximum)?;
    *rate = rate.max(local);
    Ok(())
}
