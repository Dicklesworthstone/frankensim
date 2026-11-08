//! Finite piston inertia coupled to original NASA/Peng--Robinson pressure work.

use super::*;
use crate::dynamics::MovingWorkspace;
use crate::dynamics::moving::reported_step;
use crate::dynamics::piston::{
    MechanicalExchange, Piston, PistonControl, PistonError, PistonLoad, PistonStep,
    PistonSuspension, SuspendedPistonStep, mechanical_report, solve, suspension, validate,
};

mod thermal;

impl RealGasEuler1d {
    /// Advance the original real-gas wave field and a finite-mass right piston.
    ///
    /// The shared midpoint mechanics solves the piston impulse equation using
    /// the pressure impulse from each original moving-volume RK2 trial. That
    /// SAME wall traction supplies the gas pressure work; no ideal-gas pressure,
    /// frozen gamma, or independent compression-energy estimate is substituted.
    /// Both numerical reflecting and original-source physical shock/isentrope
    /// walls retain their declared boundary law and its admission controls.
    ///
    /// The left wall is stationary and the right wall moves at the solved
    /// midpoint speed. Source, geometry, clock, piston velocity and cached
    /// thermodynamic states publish together only after independent impulse and
    /// gas-plus-piston energy checks. All nonlinear proposals share the same
    /// accepted initial field and duration. Retries restart the whole coupled
    /// interval; numerical workspace is reused instead of cloning the field.
    ///
    /// Requires constant area and two impermeable ends. The original homogeneous
    /// real-gas branch restrictions remain; this does not compose phase change,
    /// heat, chemistry, distributed forcing or additional moving masses.
    ///
    /// # Errors
    /// Invalid controls, unsupported boundaries, source or CFL refusal, exhausted
    /// nonlinear/retry budgets and cancellation leave accepted state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        self.step_with_piston_suspension(
            piston,
            maximum_duration,
            gas_controls,
            load,
            None,
            controls,
            None,
            cancelled,
        )
        .map(|step| step.piston)
    }

    /// Couple a stationary spring/damper to the same real-gas piston solve.
    ///
    /// The existing discrete spring law uses actual midpoint displacement and
    /// the damper uses piston midpoint velocity. Its storage and nonnegative
    /// dissipation are separately retained by [`SuspendedPistonStep`]; damping
    /// is not injected again as gas heat. Pressure work and every candidate gas
    /// state still resolve through the original real-gas owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_suspended_piston(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        spring: PistonSuspension,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep, PistonError> {
        self.step_with_piston_suspension(
            piston,
            maximum_duration,
            gas_controls,
            load,
            Some(spring),
            controls,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_with_piston_suspension(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        spring: Option<PistonSuspension>,
        controls: PistonControl,
        mut walls: Option<(&mut RealGasWallOperator<'_>, WallControl)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep, PistonError> {
        poll(&mut cancelled)?;
        validate(*piston, maximum_duration, gas_controls, load, controls)?;
        if let Some((walls, wall_controls)) = &walls {
            walls.validate_binding(self)?;
            wall_controls.validate().map_err(RealGasFlowError::from)?;
        }
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
            return Err(PistonError::Invalid(
                "real-gas piston requires constant area and two impermeable ends",
            ));
        }
        let end = self.states.len();
        let spring = spring
            .map(|law| law.bind(self.mesh.positions()[end]))
            .transpose()?;
        let area = self.mesh.face_areas()[0];
        if !(load.external_force_n - area * load.exterior_pressure_pa).is_finite() {
            return Err(PistonError::Invalid("external piston force overflow"));
        }
        if gas_controls.reconstruction == Reconstruction::Muscl && self.work.traces.is_none() {
            self.work.traces = Some(allocated(
                end,
                [(Conserved::default(), Primitive::default()); 2],
            )?);
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh)?);
        }
        if walls.is_some() && self.force_work.is_none() {
            self.force_work = Some(forcing::Workspace::new(end)?);
        }
        let moving = self
            .moving
            .as_mut()
            .expect("admitted real-gas piston workspace");
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s;
        let mut rate = 0.0_f64;
        for (i, p) in self.primitive.iter().enumerate() {
            poll(&mut cancelled)?;
            let relative = p
                .velocity
                .abs()
                .max((p.velocity - lo).abs())
                .max((p.velocity - hi).abs());
            let local = (relative + p.sound_speed) * area / self.mesh.volumes()[i];
            if !positive(local) {
                return Err(PistonError::Invalid(
                    "real-gas piston relative-wave CFL overflow",
                ));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        if lo < 0.0 {
            let length = self.mesh.positions()[end] - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / -lo));
        }
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        if before.array().iter().any(|value| !value.is_finite()) {
            return Err(PistonError::Invalid("nonfinite real-gas piston inventory"));
        }
        let mut last = PistonError::Invalid("no admissible coupled real-gas piston increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || !(self.time + dt).is_finite() || self.time + dt <= self.time {
                return Err(PistonError::Invalid(
                    "real-gas piston step cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let mut wave_failure = false;
            let result = (|| {
                let first = if let Some((walls, wall_controls)) = &mut walls {
                    let source = self
                        .force_work
                        .as_mut()
                        .expect("prepared gas source buffers");
                    for (i, &state) in self.states.iter().enumerate() {
                        poll(&mut cancelled)?;
                        source.states[i] = state;
                        source.primitive[i] = self.primitive[i];
                    }
                    Some(walls.prepare_piston_first(
                        &mut source.states,
                        &mut source.primitive,
                        Time::new(0.5 * dt),
                        *wall_controls,
                        &mut cancelled,
                    )?)
                } else {
                    None
                };
                let (states, primitives, transport_before) = if first.is_some() {
                    let source = self
                        .force_work
                        .as_ref()
                        .expect("prepared gas source buffers");
                    let total = totals(&source.states, self.mesh.volumes(), &mut cancelled)?;
                    (&source.states[..], &source.primitive[..], total)
                } else {
                    (&self.states[..], &self.primitive[..], before)
                };
                let report = solve(piston.velocity_m_s(), controls, |speed| {
                    gas_trials += 1;
                    moving.prepare_velocities(&self.mesh, [0.0, speed], &mut cancelled)?;
                    let stages = moving::attempt(
                        &self.gas,
                        states,
                        primitives,
                        &self.mesh,
                        self.boundaries,
                        self.boundary_primitive,
                        reconstruction,
                        dt,
                        gas_controls.cfl,
                        &mut self.work,
                        moving,
                        &mut cancelled,
                    )
                    .map_err(|error| {
                        wave_failure = !matches!(
                            error,
                            RealGasFlowError::Cancelled
                                | RealGasFlowError::Numerical(Error::CflExceeded { .. })
                        );
                        error
                    })?;
                    let after =
                        totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                    let gas = reported_step(transport_before, after, stages, dt, area)?;
                    let exchange = MechanicalExchange::from(&gas);
                    let displacement = moving.mesh.positions()[end] - self.mesh.positions()[end];
                    let trial_load =
                        spring.map_or(Ok(load), |spring| spring.load(load, speed, displacement))?;
                    mechanical_report(
                        gas,
                        exchange,
                        *piston,
                        trial_load,
                        speed,
                        displacement,
                        area,
                    )
                })?;
                let mut step = suspension::finish(spring, report, load, dt, controls)?;
                let thermal = if let (Some((walls, wall_controls)), Some(first)) =
                    (&mut walls, first)
                {
                    let second = walls.prepare_piston_second(
                        &mut self.work.candidate,
                        &mut self.work.candidate_primitive,
                        &moving.mesh,
                        Time::new(0.5 * dt),
                        *wall_controls,
                        &mut cancelled,
                    )?;
                    for (i, p) in self.work.candidate_primitive.iter().enumerate() {
                        poll(&mut cancelled)?;
                        let rate =
                            moving.relative_signal_rate(&moving.mesh, i, p.velocity, p.sound_speed);
                        if !positive(rate) {
                            return Err(PistonError::Invalid(
                                "thermal piston endpoint CFL overflow",
                            ));
                        }
                        cfl_gate(rate * dt, gas_controls.cfl)?;
                        step.piston.gas.maximum_cfl = step.piston.gas.maximum_cfl.max(rate * dt);
                    }
                    let after =
                        totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                    Some(thermal::finish(
                        &mut step,
                        before,
                        after,
                        [first, second],
                        *piston,
                        load,
                        area,
                        controls,
                    )?)
                } else {
                    None
                };
                Ok((step, thermal))
            })();
            match result {
                Ok((mut step, thermal)) => {
                    let report = &mut step.piston;
                    report.gas.halvings = halvings;
                    report.gas.first_order_retry =
                        first_order_retry && gas_controls.reconstruction == Reconstruction::Muscl;
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.primitive, &mut self.work.candidate_primitive);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    piston.accept(report);
                    if let (Some((walls, _)), Some(thermal)) = (&mut walls, thermal) {
                        walls.commit_piston(thermal);
                    }
                    return Ok(step);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::Gas(Error::Allocation)) => {
                    return Err(PistonError::Gas(Error::Allocation));
                }
                Err(error) => {
                    // A refused heat half or endpoint CFL needs a shorter
                    // complete interval, not a less accurate pressure field.
                    if wave_failure {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(PistonError::StepRetries {
            attempts: gas_controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}
