//! Two original real-gas wave fields transferring work through one finite piston.

use super::*;
use crate::multicomponent::SharedPistonGasReport;

impl RealGasEuler1d {
    /// Advance two resolved real-gas chambers and their shared inertial piston.
    ///
    /// `self` lies left of `right` in one axial coordinate system. Their shared
    /// face position, constant cross section and accepted clock must match
    /// exactly. Each sealed chamber keeps its own original NASA/Peng--Robinson
    /// source, frozen composition, energy reference, grid and wall pressure law.
    /// Neither chamber pressure is prescribed or equated to the other side.
    /// The piston has negligible thickness; `external_force_n` excludes both
    /// resolved gas pressures and is positive to the right.
    ///
    /// Every midpoint-speed proposal advances BOTH native moving-volume RK2
    /// fields at one duration, with end speeds `[0,w]` and `[w,0]`. Positive
    /// speed expands the left chamber and compresses the right. The same two
    /// interface impulses supply `M*(v1-v0)` and both gas pressure works.
    /// Actual gas-energy increments, piston kinetic change and applied-force
    /// displacement work must close independently. The two fixed-end reactions
    /// close total momentum and do no work; their sum is the reported fixed-wall
    /// impulse. No exterior-pressure work is added for an already resolved gas.
    ///
    /// Acoustic and contraction limits cover both fields and the complete speed
    /// bracket. Refused trials restart both fields at a shorter common duration.
    /// Native state, thermodynamic and geometry buffers are reused; no field or
    /// source table is cloned per nonlinear proposal. Both fields, caches,
    /// meshes, clocks and the piston publish together after the final poll.
    /// `gas_trials` counts paired proposals, each with at most two native field
    /// evaluations; an early left-field refusal skips the right-field trial.
    ///
    /// This entry is adiabatic, impermeable, constant-area real-gas Euler flow.
    /// Thermal walls, body forcing, chemistry, phase flash and leakage are not
    /// composed here. Original source and near-wall accuracy limits remain;
    /// the scalar mechanical solve is not a general FSI stability theorem.
    ///
    /// # Errors
    /// Invalid bindings or controls, original source/geometry/CFL refusal,
    /// exhausted iteration/retry budgets and cancellation leave all accepted
    /// gas, cache, geometry, clock and piston state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston(
        &mut self,
        right: &mut RealGasEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonGasReport<StepReport>>, PistonError> {
        self.step_with_shared_piston_suspension(
            right,
            piston,
            maximum_duration,
            gas_controls,
            external_force_n,
            None,
            controls,
            cancelled,
        )
        .map(|step| step.piston)
    }

    /// Couple both real-gas pressures to one stationary spring/damper suspension.
    ///
    /// The original midpoint spring and damper act within every joint pressure
    /// trial. Actual spring storage and independently integrated nonnegative
    /// damping join both gas energies and piston kinetic energy in the accepted
    /// balance. Suspension impulse has an opposite stationary-anchor reaction.
    /// Neither gas receives the damping heat. All bindings, source limits and
    /// atomic refusal semantics of [`Self::step_with_shared_piston`] apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_suspended_piston(
        &mut self,
        right: &mut RealGasEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        spring: PistonSuspension,
        controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<SharedPistonGasReport<StepReport>>, PistonError> {
        self.step_with_shared_piston_suspension(
            right,
            piston,
            maximum_duration,
            gas_controls,
            external_force_n,
            Some(spring),
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_with_shared_piston_suspension(
        &mut self,
        right: &mut RealGasEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        spring: Option<PistonSuspension>,
        controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<SharedPistonGasReport<StepReport>>, PistonError> {
        poll(&mut cancelled)?;
        let load = PistonLoad {
            exterior_pressure_pa: 0.0,
            external_force_n,
        };
        validate(*piston, maximum_duration, gas_controls, load, controls)?;
        let area = self
            .mesh
            .uniform_cross_section()
            .ok_or(PistonError::Invalid(
                "shared real-gas piston requires constant-area chambers",
            ))?;
        let interface = self.mesh.positions()[self.states.len()];
        if !closed(self.boundaries)
            || !closed(right.boundaries)
            || right.mesh.uniform_cross_section() != Some(area)
            || right.mesh.positions()[0] != interface
            || self.time != right.time
        {
            return Err(PistonError::Invalid(
                "shared real-gas piston requires impermeable chambers with equal areas, interface coordinates and clocks",
            ));
        }
        let spring = spring.map(|law| law.bind(interface)).transpose()?;
        let original_piston = *piston;
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s;
        let (left_before, left_rate) = prepare(self, gas_controls, [lo, hi], &mut cancelled)?;
        let (right_before, right_rate) = prepare(right, gas_controls, [lo, hi], &mut cancelled)?;
        let mut dt = maximum_duration
            .value()
            .min(0.9 * gas_controls.cfl / left_rate.max(right_rate));
        if lo < 0.0 {
            dt = dt.min(0.9 * gas_controls.cfl * ((interface - self.mesh.positions()[0]) / -lo));
        }
        if hi > 0.0 {
            let length = right.mesh.positions()[right.states.len()] - interface;
            dt = dt.min(0.9 * gas_controls.cfl * (length / hi));
        }
        let mut last = PistonError::Invalid("no admissible shared real-gas piston increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || !(self.time + dt).is_finite() || self.time + dt <= self.time {
                return Err(PistonError::Invalid(
                    "shared real-gas piston cannot advance both clocks",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let mut wave_failure = false;
            let result = solve(original_piston.velocity_m_s(), controls, |speed| {
                gas_trials += 1;
                let left = trial(
                    self,
                    left_before,
                    dt,
                    gas_controls.cfl,
                    reconstruction,
                    [0.0, speed],
                    &mut wave_failure,
                    &mut cancelled,
                )?;
                let right_report = trial(
                    right,
                    right_before,
                    dt,
                    gas_controls.cfl,
                    reconstruction,
                    [speed, 0.0],
                    &mut wave_failure,
                    &mut cancelled,
                )?;
                let left_mesh = &self.moving.as_ref().expect("prepared left motion").mesh;
                let right_mesh = &right.moving.as_ref().expect("prepared right motion").mesh;
                let next_interface = left_mesh.positions()[self.states.len()];
                if next_interface != right_mesh.positions()[0] {
                    return Err(PistonError::Invalid(
                        "shared real-gas candidate endpoints do not represent one piston",
                    ));
                }
                let exchange = MechanicalExchange {
                    duration_s: dt,
                    momentum_change: sum([
                        left.after.momentum,
                        -left.before.momentum,
                        right_report.after.momentum,
                        -right_report.before.momentum,
                    ]),
                    energy_change: sum([
                        left.after.energy,
                        -left.before.energy,
                        right_report.after.energy,
                        -right_report.before.energy,
                    ]),
                    left_impulse: sum([
                        left.left_boundary_change.momentum,
                        right_report.right_boundary_change.momentum,
                    ]),
                    right_impulse: sum([
                        left.right_boundary_change.momentum,
                        right_report.left_boundary_change.momentum,
                    ]),
                    right_work: sum([
                        left.right_boundary_change.energy,
                        right_report.left_boundary_change.energy,
                    ]),
                };
                let displacement = next_interface - interface;
                let trial_load =
                    spring.map_or(Ok(load), |law| law.load(load, speed, displacement))?;
                mechanical_report(
                    SharedPistonGasReport {
                        left,
                        right: right_report,
                    },
                    exchange,
                    original_piston,
                    trial_load,
                    speed,
                    displacement,
                    area,
                )
            })
            .and_then(|report| suspension::finish(spring, report, load, dt, controls));
            match result {
                Ok(mut step) => {
                    let report = &mut step.piston;
                    for side in [&mut report.gas.left, &mut report.gas.right] {
                        side.halvings = halvings;
                        side.first_order_retry = first_order_retry
                            && gas_controls.reconstruction == Reconstruction::Muscl;
                    }
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    for flow in [self, right] {
                        let moving = flow.moving.as_mut().expect("prepared accepted motion");
                        core::mem::swap(&mut flow.states, &mut flow.work.candidate);
                        core::mem::swap(&mut flow.primitive, &mut flow.work.candidate_primitive);
                        core::mem::swap(&mut flow.mesh, &mut moving.mesh);
                        flow.time += dt;
                    }
                    piston.accept(report);
                    return Ok(step);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::Gas(Error::Allocation)) => {
                    return Err(PistonError::Gas(Error::Allocation));
                }
                Err(error) => {
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

fn closed(boundaries: RealGasBoundaries) -> bool {
    let wall = |end| {
        matches!(
            end,
            RealGasBoundary::Reflecting | RealGasBoundary::PhysicalReflecting(_)
        )
    };
    matches!(boundaries, RealGasBoundaries::Ends { left, right } if wall(left) && wall(right))
}

fn prepare(
    flow: &mut RealGasEuler1d,
    controls: StepControl,
    bracket: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Totals, f64), PistonError> {
    if controls.reconstruction == Reconstruction::Muscl && flow.work.traces.is_none() {
        flow.work.traces = Some(allocated(
            flow.states.len(),
            [(Conserved::default(), Primitive::default()); 2],
        )?);
    }
    if flow.moving.is_none() {
        flow.moving = Some(MovingWorkspace::new(&flow.mesh)?);
    }
    let area = flow.mesh.face_areas()[0];
    let mut rate = 0.0_f64;
    for (i, p) in flow.primitive.iter().enumerate() {
        poll(cancelled)?;
        let relative = p
            .velocity
            .abs()
            .max((p.velocity - bracket[0]).abs())
            .max((p.velocity - bracket[1]).abs());
        let local = (relative + p.sound_speed) * area / flow.mesh.volumes()[i];
        if !positive(local) {
            return Err(PistonError::Invalid(
                "shared real-gas relative-wave CFL overflow",
            ));
        }
        rate = rate.max(local);
    }
    let before = totals(&flow.states, flow.mesh.volumes(), cancelled)?;
    if before.array().iter().any(|value| !value.is_finite()) {
        return Err(PistonError::Invalid("nonfinite shared real-gas inventory"));
    }
    Ok((before, rate))
}

#[allow(clippy::too_many_arguments)]
fn trial(
    flow: &mut RealGasEuler1d,
    before: Totals,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    wall_velocity: [f64; 2],
    wave_failure: &mut bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StepReport, PistonError> {
    let moving = flow
        .moving
        .as_mut()
        .expect("prepared shared real-gas motion");
    moving.prepare_velocities(&flow.mesh, wall_velocity, cancelled)?;
    let stages = moving::attempt(
        &flow.gas,
        &flow.states,
        &flow.primitive,
        &flow.mesh,
        flow.boundaries,
        flow.boundary_primitive,
        reconstruction,
        dt,
        cfl,
        &mut flow.work,
        moving,
        cancelled,
    )
    .map_err(|error| {
        *wave_failure = !matches!(
            error,
            RealGasFlowError::Cancelled | RealGasFlowError::Numerical(Error::CflExceeded { .. })
        );
        error
    })?;
    let after = totals(&flow.work.candidate, moving.mesh.volumes(), cancelled)?;
    Ok(reported_step(
        before,
        after,
        stages,
        dt,
        flow.mesh.face_areas()[0],
    )?)
}

fn sum<const N: usize>(values: [f64; N]) -> f64 {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + correction
}
