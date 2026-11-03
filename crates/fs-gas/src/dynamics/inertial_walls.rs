//! Resolved gas exchanging momentum and work with two finite end walls.
//!
//! Both midpoint wall speeds are solved together on the original ALE RK2
//! field. Gas momentum remains resolved: no reduced wall mass, uniform cavity
//! pressure or stationary support replaces it. Every trial reuses the original
//! field/source workspace; accepted gas, mesh, clock and both walls publish
//! together. This is a bounded engineering coupling, not an FSI stability proof.

use super::evolve::totals;
use super::moving::{attempt, reported_step};
use super::piston::{Piston, PistonControl, PistonError, PistonLoad, validate};
use super::{
    Boundaries, Boundary, Error, Euler1d, MovingWorkspace, Reconstruction, StepControl, StepReport,
    poll, positive,
};
use fs_qty::Time;

/// Loads on one wall, held over the same accepted gas interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct InertialWallLoad {
    /// Absolute exterior pressure [Pa]. Its force is +A*p on the LEFT wall
    /// and -A*p on the RIGHT wall. Zero explicitly means vacuum.
    pub exterior_pressure_pa: f64,
    /// Other force [N], positive along the increasing axis for EITHER wall.
    /// Excludes the resolved gas pressure and the exterior pressure above.
    pub external_force_n: f64,
}

/// Bounded simultaneous solve and independent physical acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialWallsControl {
    /// [left, right] finite midpoint-speed intervals [m/s], each containing
    /// that wall's accepted velocity. They bound trials, not clamped outputs.
    /// Both intervals enter the relative-wave and mesh-contraction CFL.
    pub midpoint_velocity_bracket_m_s: [[f64; 2]; 2],
    /// Absolute tolerance for EACH wall and combined gas/walls impulse [N s].
    pub momentum_tolerance_n_s: f64,
    /// Absolute tolerance for EACH wall work and combined total energy [J].
    pub energy_tolerance_j: f64,
    /// Damped Newton iterations per complete attempt, in 1..=128. Each uses
    /// two one-sided Jacobian proposals and at most eight backtracking trials.
    pub maximum_iterations: usize,
}

pub(crate) fn validate_controls(
    walls: [Piston; 2],
    maximum_duration: Time,
    gas_controls: StepControl,
    loads: [InertialWallLoad; 2],
    controls: InertialWallsControl,
) -> Result<(), PistonError> {
    for i in 0..2 {
        validate(
            walls[i],
            maximum_duration,
            gas_controls,
            PistonLoad {
                exterior_pressure_pa: loads[i].exterior_pressure_pa,
                external_force_n: loads[i].external_force_n,
            },
            PistonControl {
                midpoint_velocity_bracket_m_s: controls.midpoint_velocity_bracket_m_s[i],
                momentum_tolerance_n_s: controls.momentum_tolerance_n_s,
                energy_tolerance_j: controls.energy_tolerance_j,
                maximum_iterations: controls.maximum_iterations,
            },
        )?;
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s[i];
        if !positive(hi - lo) {
            return Err(PistonError::Invalid(
                "wall speed interval width is unrepresentable",
            ));
        }
    }
    Ok(())
}

/// Original accepted field and actual extensive exchanges; all arrays are
/// [left, right]. There is no unmodeled fixed support supplying momentum.
/// `R` retains the original frozen or evolving-species report, including all
/// species balances and the source owner's energy-reference convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialWallsStep<R = StepReport> {
    /// Original extensive RK2 field and per-boundary exchange report.
    pub gas: R,
    /// Accepted input wall velocities along the common axial coordinate [m/s].
    pub velocity_before_m_s: [f64; 2],
    /// Jointly accepted endpoint wall velocities [m/s].
    pub velocity_after_m_s: [f64; 2],
    /// Solved constant mesh speeds, equal to the wall midpoint velocities [m/s].
    pub midpoint_velocity_m_s: [f64; 2],
    /// Actual represented endpoint-face displacement, not assumed dt*w [m].
    pub displacement_m: [f64; 2],
    /// Negatives of the original gas-side boundary impulses [N s].
    pub gas_impulse_n_s: [f64; 2],
    /// Independently evaluated changes in each wall's kinetic storage [J].
    pub kinetic_energy_change_j: [f64; 2],
    /// Other applied force times the actual represented wall displacement [J].
    pub external_force_work_j: [f64; 2],
    /// Work AGAINST exterior pressure: [-p_L*A*dx_L, p_R*A*dx_R] [J].
    pub exterior_pressure_work_j: [f64; 2],
    /// Each wall momentum change minus its gas and exterior/load impulses [N s].
    pub momentum_residual_n_s: [f64; 2],
    /// Actual gas plus both wall momentum change minus external impulse [N s].
    pub combined_momentum_residual_n_s: f64,
    /// Each wall kinetic change plus gas-side work minus exterior/load work [J].
    pub mechanical_energy_residual_j: [f64; 2],
    /// Actual gas plus both wall energy change minus exterior/load work [J].
    pub total_energy_residual_j: f64,
    /// Completed Newton iterations in the accepted coupled attempt.
    pub iterations: usize,
    /// Complete gas proposals, including Jacobian, line-search and retry work.
    pub gas_trials: usize,
}

impl<R> InertialWallsStep<R> {
    fn residual_norm(&self) -> f64 {
        self.momentum_residual_n_s[0]
            .abs()
            .max(self.momentum_residual_n_s[1].abs())
    }
    fn defects(&self) -> (f64, f64) {
        (
            self.residual_norm()
                .max(self.combined_momentum_residual_n_s.abs()),
            self.mechanical_energy_residual_j[0]
                .abs()
                .max(self.mechanical_energy_residual_j[1].abs())
                .max(self.total_energy_residual_j.abs()),
        )
    }
    pub(crate) fn admitted(&self, controls: InertialWallsControl) -> bool {
        let (momentum, energy) = self.defects();
        momentum <= controls.momentum_tolerance_n_s && energy <= controls.energy_tolerance_j
    }
}

impl Euler1d {
    /// Advance a sealed constant-area field and BOTH finite end walls.
    ///
    /// Both wall velocities use the increasing axial coordinate. The left
    /// wall moves the first mesh face and the right wall the last. At each
    /// simultaneous midpoint-speed trial, the original moving-volume operator
    /// supplies integrated pressure impulses AND their conjugate gas work.
    /// Wall kinetic changes, each work balance, combined momentum and total
    /// energy are checked independently; no state or energy correction is made.
    ///
    /// A row-scaled two-variable Newton solve uses bounded finite differences
    /// and backtracking. Every proposal starts from the SAME accepted field
    /// at one duration. Failure halves the entire interval under the original
    /// gas controls. Scratch is reused; no field or source table is cloned.
    ///
    /// This closed, adiabatic, frozen-composition Euler lane has no gravity,
    /// heat transfer, contact stops, leakage, transverse wall motion or housing
    /// mass beyond the two explicit walls. Each wall must remain associated
    /// with its own end of this field. Cancellation/refusal publishes nothing.
    ///
    /// # Errors
    /// Refuses unsupported geometry, invalid controls or loads, unrepresentable
    /// states, exhausted physical/root retries and cancellation before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<InertialWallsStep, PistonError> {
        poll(&mut cancelled)?;
        validate_controls(*walls, maximum_duration, gas_controls, loads, controls)?;
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.hydrostatic.is_some()
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(PistonError::Invalid(
                "two inertial walls require a constant-area closed interval without gravity",
            ));
        }
        let area = self.mesh.face_areas()[0];
        let force = [
            loads[0].external_force_n + area * loads[0].exterior_pressure_pa,
            loads[1].external_force_n - area * loads[1].exterior_pressure_pa,
        ];
        if force.iter().any(|value| !value.is_finite()) {
            return Err(PistonError::Invalid("external wall force overflow"));
        }
        let bounds = controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let mut rate = 0.0_f64;
        for (i, &state) in self.states.iter().enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            let relative = (p.velocity - low).abs().max((p.velocity - high).abs());
            let local = (relative + p.sound_speed) * (area / self.mesh.volumes()[i]);
            if !positive(local) {
                return Err(PistonError::Invalid("two-wall relative-wave CFL overflow"));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(PistonError::Invalid("two-wall contraction bound overflow"));
        }
        if contraction > 0.0 {
            let length = self.mesh.positions()[self.states.len()] - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / contraction));
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh)?);
        }
        let moving = self.moving.as_mut().expect("moving workspace initialized");
        let original = *walls;
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        let mut last = PistonError::Invalid("no admissible two-wall increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(PistonError::Invalid(
                    "two-wall step cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve(
                original.map(|wall| wall.velocity_m_s()),
                controls,
                |speed| {
                    gas_trials += 1;
                    moving.prepare_velocities(&self.mesh, speed, &mut cancelled)?;
                    let (a, b) = attempt(
                        &self.gas,
                        &self.states,
                        &self.mesh,
                        self.boundaries,
                        dt,
                        gas_controls.cfl,
                        reconstruction,
                        &mut self.work,
                        moving,
                        &mut cancelled,
                    )?;
                    let after =
                        totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                    let gas = reported_step(before, after, [a, b], dt, area)?;
                    let exchange = MechanicalExchange::from(&gas);
                    let displacement = [0, self.states.len()]
                        .map(|end| moving.mesh.positions()[end] - self.mesh.positions()[end]);
                    mechanical_report(gas, exchange, original, loads, speed, displacement, area)
                },
            );
            match result {
                Ok(mut report) => {
                    // Validate both endpoints before the final cancellation gate;
                    // nothing fallible may follow publication of any owner.
                    let next = [
                        Piston::new(original[0].mass_kg(), report.velocity_after_m_s[0])?,
                        Piston::new(original[1].mass_kg(), report.velocity_after_m_s[1])?,
                    ];
                    report.gas.halvings = halvings;
                    report.gas.first_order_retry =
                        first_order_retry && gas_controls.reconstruction == Reconstruction::Muscl;
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    *walls = next;
                    return Ok(report);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::Gas(Error::Allocation)) => {
                    return Err(PistonError::Gas(Error::Allocation));
                }
                Err(error) => {
                    if matches!(&error, PistonError::Gas(cause) if !matches!(cause, Error::CflExceeded { .. }))
                    {
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

// Only these original extensive exchanges enter the shared wall mechanics.
// The complete species/source report remains on the accepted result.
pub(crate) struct MechanicalExchange {
    duration_s: f64,
    momentum_change: f64,
    energy_change: f64,
    boundary_impulse: [f64; 2],
    boundary_work: [f64; 2],
}

impl From<&StepReport> for MechanicalExchange {
    fn from(report: &StepReport) -> Self {
        Self {
            duration_s: report.elapsed.value(),
            momentum_change: report.after.momentum - report.before.momentum,
            energy_change: report.after.energy - report.before.energy,
            boundary_impulse: [
                report.left_boundary_change.momentum,
                report.right_boundary_change.momentum,
            ],
            boundary_work: [
                report.left_boundary_change.energy,
                report.right_boundary_change.energy,
            ],
        }
    }
}

impl From<&crate::multicomponent::StepReport> for MechanicalExchange {
    fn from(report: &crate::multicomponent::StepReport) -> Self {
        Self {
            duration_s: report.elapsed.value(),
            momentum_change: report.after.momentum - report.before.momentum,
            energy_change: report.after.energy - report.before.energy,
            boundary_impulse: [
                report.left_boundary_change.momentum,
                report.right_boundary_change.momentum,
            ],
            boundary_work: [
                report.left_boundary_change.energy,
                report.right_boundary_change.energy,
            ],
        }
    }
}

pub(crate) fn mechanical_report<R>(
    gas: R,
    exchange: MechanicalExchange,
    walls: [Piston; 2],
    loads: [InertialWallLoad; 2],
    speed: [f64; 2],
    displacement: [f64; 2],
    area: f64,
) -> Result<InertialWallsStep<R>, PistonError> {
    let before = walls.map(|wall| wall.velocity_m_s());
    let after = core::array::from_fn(|i| speed[i] + (speed[i] - before[i]));
    let momentum: [f64; 2] = core::array::from_fn(|i| walls[i].mass_kg() * (after[i] - before[i]));
    let kinetic = core::array::from_fn(|i| {
        (0.5 * walls[i].mass_kg() * (after[i] - before[i])) * (after[i] + before[i])
    });
    let exterior = core::array::from_fn(|i| {
        [-1.0, 1.0][i] * area * loads[i].exterior_pressure_pa * displacement[i]
    });
    let external = core::array::from_fn(|i| loads[i].external_force_n * displacement[i]);
    let external_impulse: [f64; 2] = core::array::from_fn(|i| {
        exchange.duration_s
            * (loads[i].external_force_n - [-1.0, 1.0][i] * area * loads[i].exterior_pressure_pa)
    });
    let residual =
        core::array::from_fn(|i| momentum[i] + exchange.boundary_impulse[i] - external_impulse[i]);
    let work = core::array::from_fn(|i| {
        kinetic[i] + exchange.boundary_work[i] + exterior[i] - external[i]
    });
    let combined = exchange.momentum_change
        + (momentum[0] - external_impulse[0])
        + (momentum[1] - external_impulse[1]);
    let total = exchange.energy_change
        + (kinetic[0] + exterior[0] - external[0])
        + (kinetic[1] + exterior[1] - external[1]);
    if [
        after,
        momentum,
        kinetic,
        exterior,
        external,
        external_impulse,
        residual,
        work,
        [combined, total],
    ]
    .into_iter()
    .flatten()
    .any(|value| !value.is_finite())
    {
        return Err(PistonError::Invalid("two-wall momentum or energy overflow"));
    }
    Ok(InertialWallsStep {
        gas,
        velocity_before_m_s: before,
        velocity_after_m_s: after,
        midpoint_velocity_m_s: speed,
        displacement_m: displacement,
        gas_impulse_n_s: exchange.boundary_impulse.map(|impulse| -impulse),
        kinetic_energy_change_j: kinetic,
        external_force_work_j: external,
        exterior_pressure_work_j: exterior,
        momentum_residual_n_s: residual,
        combined_momentum_residual_n_s: combined,
        mechanical_energy_residual_j: work,
        total_energy_residual_j: total,
        iterations: 0,
        gas_trials: 0,
    })
}

pub(crate) fn solve<R>(
    initial: [f64; 2],
    controls: InertialWallsControl,
    mut trial: impl FnMut([f64; 2]) -> Result<InertialWallsStep<R>, PistonError>,
) -> Result<InertialWallsStep<R>, PistonError> {
    let mut speed = initial;
    let mut current = trial(speed)?;
    for iteration in 0..=controls.maximum_iterations {
        // Only this CURRENT proposal may be returned: it must still own the
        // native candidate buffers, including after derivative/line searches.
        if current.admitted(controls) {
            current.iterations = iteration;
            return Ok(current);
        }
        let norm = current.residual_norm();
        if norm == 0.0 {
            let (momentum, energy) = current.defects();
            return Err(PistonError::Balance {
                momentum_residual_n_s: momentum,
                energy_residual_j: energy,
            });
        }
        if iteration == controls.maximum_iterations {
            break;
        }
        let mut jacobian = [[0.0; 2]; 2];
        for column in 0..2 {
            let [lo, hi] = controls.midpoint_velocity_bracket_m_s[column];
            let h = f64::EPSILON.sqrt() * (hi - lo).max(speed[column].abs()).max(1.0);
            let mut probe = speed;
            probe[column] = if speed[column] <= f64::midpoint(lo, hi) {
                (speed[column] + h).min(hi)
            } else {
                (speed[column] - h).max(lo)
            };
            let delta = probe[column] - speed[column];
            if delta == 0.0 {
                return Err(PistonError::NotConverged {
                    iterations: iteration,
                    impulse_residual_n_s: norm,
                });
            }
            let value = trial(probe)?;
            for row in 0..2 {
                jacobian[row][column] =
                    (value.momentum_residual_n_s[row] - current.momentum_residual_n_s[row]) / delta;
            }
        }
        let correction = newton_correction(jacobian, current.momentum_residual_n_s)?;
        let mut factor = 1.0_f64;
        for i in 0..2 {
            let [lo, hi] = controls.midpoint_velocity_bracket_m_s[i];
            if correction[i] > 0.0 {
                factor = factor.min((hi - speed[i]) / correction[i]);
            }
            if correction[i] < 0.0 {
                factor = factor.min((lo - speed[i]) / correction[i]);
            }
        }
        let mut accepted = None;
        let mut last = PistonError::NotConverged {
            iterations: iteration,
            impulse_residual_n_s: norm,
        };
        for _ in 0..8 {
            let proposed = core::array::from_fn(|i| speed[i] + factor * correction[i]);
            if proposed == speed || factor <= 0.0 {
                break;
            }
            if (0..2).all(|i| {
                proposed[i] >= controls.midpoint_velocity_bracket_m_s[i][0]
                    && proposed[i] <= controls.midpoint_velocity_bracket_m_s[i][1]
            }) {
                match trial(proposed) {
                    Ok(value) if value.admitted(controls) || value.residual_norm() < norm => {
                        accepted = Some((proposed, value));
                        break;
                    }
                    Ok(_) => {}
                    Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                    Err(PistonError::Gas(Error::Allocation)) => {
                        return Err(PistonError::Gas(Error::Allocation));
                    }
                    Err(PistonError::MixtureGas(crate::multicomponent::Error::Allocation)) => {
                        return Err(PistonError::MixtureGas(
                            crate::multicomponent::Error::Allocation,
                        ));
                    }
                    Err(error) => last = error,
                }
            }
            factor *= 0.5;
        }
        (speed, current) = accepted.ok_or(last)?;
    }
    Err(PistonError::NotConverged {
        iterations: controls.maximum_iterations,
        impulse_residual_n_s: current.residual_norm(),
    })
}

fn newton_correction(
    mut matrix: [[f64; 2]; 2],
    mut residual: [f64; 2],
) -> Result<[f64; 2], PistonError> {
    for i in 0..2 {
        let scale = matrix[i][0].abs().max(matrix[i][1].abs());
        if !positive(scale) || matrix[i].iter().any(|value| !value.is_finite()) {
            return Err(PistonError::Invalid(
                "unrepresentable two-wall impulse Jacobian",
            ));
        }
        matrix[i] = matrix[i].map(|value| value / scale);
        residual[i] /= scale;
    }
    let [[a, b], [c, d]] = matrix;
    let determinant = a * d - b * c;
    if determinant.abs() <= 64.0 * f64::EPSILON {
        return Err(PistonError::Invalid("singular two-wall impulse Jacobian"));
    }
    let step = [
        (b * residual[1] - d * residual[0]) / determinant,
        (c * residual[0] - a * residual[1]) / determinant,
    ];
    if step.iter().any(|value| !value.is_finite()) {
        return Err(PistonError::Invalid("two-wall Newton increment overflow"));
    }
    Ok(step)
}

#[cfg(test)]
mod tests;
