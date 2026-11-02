//! A resolved Euler field driving its inertial right end wall.
//!
//! Each scalar trial reuses the original moving-volume RK2 operator at one
//! common duration. If `w = (v0+v1)/2`, the discrete identity
//! `K1-K0 = w*M*(v1-v0)` pairs its mechanical update with the SAME integrated
//! pressure impulse and work that changed the gas. The left wall is fixed.
//! This is a bounded engineering coupling, not an FSI stability theorem.

use super::evolve::totals;
use super::moving::{attempt, reported_step};
use super::{
    Boundaries, Boundary, Error as GasError, Euler1d, MovingWorkspace, Reconstruction, StepControl,
    StepReport, poll, positive,
};
use fs_qty::Time;

/// Inertial right wall. Its position and area belong to the gas mesh alone.
/// Positive axial velocity expands the gas interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Piston {
    mass_kg: f64,
    velocity_m_s: f64,
}

impl Piston {
    /// Bind an explicit positive inertial mass and finite initial wall velocity.
    pub fn new(mass_kg: f64, velocity_m_s: f64) -> Result<Self, PistonError> {
        if !positive(mass_kg)
            || !velocity_m_s.is_finite()
            || !(0.5 * mass_kg * velocity_m_s * velocity_m_s).is_finite()
        {
            return Err(PistonError::Invalid(
                "positive mass and finite wall kinetic energy required",
            ));
        }
        Ok(Self {
            mass_kg,
            velocity_m_s,
        })
    }

    /// Accepted wall mass [kg], distinct from the conserved gas mass.
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Accepted axial wall velocity [m/s].
    #[must_use]
    pub const fn velocity_m_s(&self) -> f64 {
        self.velocity_m_s
    }
}

/// Loads held constant throughout one accepted coupled step.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PistonLoad {
    /// Absolute pressure outside the right wall [Pa]; zero means vacuum.
    pub exterior_pressure_pa: f64,
    /// Additional force [N], positive towards gas expansion. Gas and exterior
    /// pressure forces are supplied by the coupling and must not be added here.
    pub external_force_n: f64,
}

/// Explicit scalar-solve limits and independently admitted physical balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonControl {
    /// Finite bracket for the midpoint wall velocity [m/s], containing the
    /// currently accepted velocity. All queried wall states must be supported
    /// by the original source. This is not a clamp on the final velocity.
    pub midpoint_velocity_bracket_m_s: [f64; 2],
    /// Maximum absolute piston and combined gas/piston impulse residual [N s].
    pub momentum_tolerance_n_s: f64,
    /// Maximum absolute mechanical and combined gas/piston energy residual [J].
    pub energy_tolerance_j: f64,
    /// Interior safeguarded secant/bisection trials in 1..=128. Up to three
    /// additional proposals check the current speed and bracket endpoints.
    pub maximum_iterations: usize,
}

/// A refusal never advances accepted gas, mesh, clock or wall velocity.
#[derive(Clone, Debug, PartialEq)]
pub enum PistonError {
    /// Invalid physical data, controls or unrepresentable balance.
    Invalid(&'static str),
    /// Original gas/source/geometry refusal with its concrete cause.
    Gas(GasError),
    /// Caller cancellation, propagated before any retry or publication.
    Cancelled,
    /// No impulse-residual sign change in the authored midpoint-speed bracket.
    OutsideBracket {
        /// Residual at the lower midpoint speed [N s].
        lower_impulse_n_s: f64,
        /// Residual at the upper midpoint speed [N s].
        upper_impulse_n_s: f64,
    },
    /// Bounded nonlinear work or floating-point progress was exhausted.
    NotConverged {
        /// Number of interior proposals already evaluated.
        iterations: usize,
        /// Last evaluated physical piston impulse residual [N s].
        impulse_residual_n_s: f64,
    },
    /// Even a zero impulse residual failed independent accepted-state balances.
    Balance {
        /// Maximum absolute mechanical or combined impulse residual [N s].
        momentum_residual_n_s: f64,
        /// Maximum absolute mechanical or combined energy residual [J].
        energy_residual_j: f64,
    },
    /// Every complete coupled attempt was refused. No trial was published.
    StepRetries {
        /// Number of attempts including the initial one.
        attempts: usize,
        /// Last concrete refusal.
        last: Box<PistonError>,
    },
}

impl core::fmt::Display for PistonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved gas/piston coupling: {self:?}")
    }
}
impl std::error::Error for PistonError {}
impl From<GasError> for PistonError {
    fn from(error: GasError) -> Self {
        if error == GasError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}

/// Actual common-duration gas and wall exchange, measured before publication.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonStep {
    /// The original extensive RK2 field and boundary-exchange report.
    pub gas: StepReport,
    /// Initial and accepted final wall velocity [m/s].
    pub velocity_m_s: [f64; 2],
    /// Constant mesh speed for this step, equal to the wall's midpoint speed.
    pub midpoint_velocity_m_s: f64,
    /// Actual represented displacement of the mesh's right endpoint [m].
    pub displacement_m: f64,
    /// Accepted two-stage pressure impulse ON the piston [N s].
    pub gas_impulse_n_s: f64,
    /// Opposite of the left-face impulse on the gas, delivered to its fixed
    /// support [N s]. This reaction does no work but carries momentum.
    pub fixed_wall_impulse_n_s: f64,
    /// Independently measured change of wall kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Additional force times accepted displacement [J].
    pub external_force_work_j: f64,
    /// Work against exterior pressure, `p_ext*A*displacement` [J].
    pub exterior_pressure_work_j: f64,
    /// `M*(v1-v0) - I_gas - dt*(F-A*p_ext)` [N s].
    pub momentum_residual_n_s: f64,
    /// Gas plus piston momentum change, including fixed-wall reaction and
    /// external pressure/force impulse [N s].
    pub combined_momentum_residual_n_s: f64,
    /// Kinetic change plus actual gas-side right-wall work, exterior work,
    /// and negative additional-force work [J].
    pub mechanical_energy_residual_j: f64,
    /// Actual gas-energy change plus wall kinetic change and exterior work,
    /// minus additional-force work [J]. No artificial energy correction.
    pub total_energy_residual_j: f64,
    /// Interior scalar proposals in the accepted coupled attempt.
    pub iterations: usize,
    /// All original spatial gas proposals, including rejected attempts.
    pub gas_trials: usize,
}

impl PistonStep {
    fn admitted(&self, controls: PistonControl) -> bool {
        self.momentum_residual_n_s
            .abs()
            .max(self.combined_momentum_residual_n_s.abs())
            <= controls.momentum_tolerance_n_s
            && self
                .mechanical_energy_residual_j
                .abs()
                .max(self.total_energy_residual_j.abs())
                <= controls.energy_tolerance_j
    }
}

impl Euler1d {
    /// Advance resolved pressure waves and one freely translating right piston.
    ///
    /// Both ends must be reflecting, the left wall is fixed, and cross section
    /// is constant. The original moving-volume RK2 operator supplies pressure
    /// impulse `I` and work for each scalar midpoint-speed trial. Mechanics
    /// solves `M*(v1-v0)=I+dt*(F-A*p_ext)` with `w=(v0+v1)/2`; geometry moves
    /// by `dt*w`. All trials in one nonlinear solve use the SAME duration and
    /// start from the same accepted field. CFL/source failure halves and resets
    /// the entire coupled attempt under `gas_controls.maximum_halvings`.
    ///
    /// Each accepted step independently closes piston impulse, combined momentum
    /// including the fixed-left-wall reaction, mechanical work and actual gas
    /// plus wall kinetic energy. Only the final evaluated trial's scratch state
    /// is published, together with geometry, time and piston velocity. Workspace
    /// is allocated on first use and reused; no full field is cloned per trial.
    ///
    /// The velocity bracket also bounds the initial relative-wave CFL. A wide
    /// bracket can therefore reduce the accepted duration. The caller retains
    /// the association between this field and its piston, and may change loads
    /// or controls between steps. No spring, friction, collision, valve, heat,
    /// gravity, multiple piston, or general FSI stability claim is implied.
    ///
    /// # Errors
    /// Invalid data, unsupported geometry/source states, failed balance or root
    /// admission, exhausted retries and cancellation leave physical state intact.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s;
        if !positive(maximum_duration.value())
            || !positive(gas_controls.cfl)
            || gas_controls.cfl > 0.45
            || gas_controls.maximum_halvings > 32
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || piston.velocity_m_s < lo
            || piston.velocity_m_s > hi
            || !positive(controls.momentum_tolerance_n_s)
            || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
            || !load.exterior_pressure_pa.is_finite()
            || load.exterior_pressure_pa < 0.0
            || !load.external_force_n.is_finite()
        {
            return Err(PistonError::Invalid(
                "invalid piston load, speed bracket, tolerances or step budget",
            ));
        }
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.hydrostatic.is_some()
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(PistonError::Invalid(
                "piston requires a constant-area closed interval without gravity",
            ));
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh)?);
        }
        let moving = self.moving.as_mut().expect("moving workspace initialized");
        let area = self.mesh.face_areas()[0];
        let force = load.external_force_n - area * load.exterior_pressure_pa;
        if !force.is_finite() {
            return Err(PistonError::Invalid("external piston force overflow"));
        }
        let mut rate = 0.0_f64;
        for (i, &q) in self.states.iter().enumerate() {
            let p = self
                .gas
                .decode(q, Some(self.work.primitive[i].temperature), &mut cancelled)?;
            self.work.primitive[i] = p;
            let relative = p
                .velocity
                .abs()
                .max((p.velocity - lo).abs())
                .max((p.velocity - hi).abs());
            let local = (relative + p.sound_speed) * area / self.mesh.volumes()[i];
            if !positive(local) {
                return Err(PistonError::Invalid("piston relative-wave CFL overflow"));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        if lo < 0.0 {
            let length = self.mesh.positions()[self.states.len()] - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / -lo));
        }
        let before = totals(&self.states, self.mesh.volumes(), &mut cancelled)?;
        let mut last = PistonError::Invalid("no admissible coupled piston increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(PistonError::Invalid("piston step cannot advance the clock"));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve(piston.velocity_m_s, controls, |speed| {
                gas_trials += 1;
                moving.prepare_velocities(&self.mesh, [0.0, speed], &mut cancelled)?;
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
                let after = totals(&self.work.candidate, moving.mesh.volumes(), &mut cancelled)?;
                let gas = reported_step(before, after, [a, b], dt, area)?;
                let end = self.states.len();
                let displacement = moving.mesh.positions()[end] - self.mesh.positions()[end];
                mechanical_report(gas, *piston, load, speed, displacement, area)
            });
            match result {
                Ok(mut report) => {
                    report.gas.halvings = halvings;
                    report.gas.first_order_retry =
                        first_order_retry && gas_controls.reconstruction == Reconstruction::Muscl;
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    piston.velocity_m_s = report.velocity_m_s[1];
                    return Ok(report);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::Gas(GasError::Allocation)) => {
                    return Err(PistonError::Gas(GasError::Allocation));
                }
                Err(error) => {
                    if matches!(&error, PistonError::Gas(cause) if !matches!(cause, GasError::CflExceeded { .. }))
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

fn mechanical_report(
    gas: StepReport,
    piston: Piston,
    load: PistonLoad,
    speed: f64,
    displacement: f64,
    area: f64,
) -> Result<PistonStep, PistonError> {
    let v0 = piston.velocity_m_s;
    let v1 = speed + (speed - v0);
    let momentum = piston.mass_kg * (v1 - v0);
    let kinetic = (0.5 * piston.mass_kg * (v1 - v0)) * (v1 + v0);
    let impulse = -gas.right_boundary_change.momentum;
    let fixed = -gas.left_boundary_change.momentum;
    let external_impulse =
        gas.elapsed.value() * (load.external_force_n - area * load.exterior_pressure_pa);
    // Use represented endpoint displacement, not an assumed unrounded dt*w.
    // A geometry/work discrepancy must pass the same independent energy gate.
    let external_work = load.external_force_n * displacement;
    let exterior_work = (area * load.exterior_pressure_pa) * displacement;
    let residual = momentum - impulse - external_impulse;
    let combined = (gas.after.momentum - gas.before.momentum) + momentum + fixed - external_impulse;
    let mechanical = (kinetic + gas.right_boundary_change.energy) + exterior_work - external_work;
    let total = ((gas.after.energy - gas.before.energy) + kinetic) + exterior_work - external_work;
    if ![
        v1,
        momentum,
        kinetic,
        impulse,
        fixed,
        external_impulse,
        external_work,
        exterior_work,
        residual,
        combined,
        mechanical,
        total,
    ]
    .iter()
    .all(|x| x.is_finite())
        || !(0.5 * piston.mass_kg * v1 * v1).is_finite()
    {
        return Err(PistonError::Invalid(
            "piston candidate kinetic energy, impulse or work overflow",
        ));
    }
    Ok(PistonStep {
        gas,
        velocity_m_s: [v0, v1],
        midpoint_velocity_m_s: speed,
        displacement_m: displacement,
        gas_impulse_n_s: impulse,
        fixed_wall_impulse_n_s: fixed,
        kinetic_energy_change_j: kinetic,
        external_force_work_j: external_work,
        exterior_pressure_work_j: exterior_work,
        momentum_residual_n_s: residual,
        combined_momentum_residual_n_s: combined,
        mechanical_energy_residual_j: mechanical,
        total_energy_residual_j: total,
        iterations: 0,
        gas_trials: 0,
    })
}

fn solve(
    initial: f64,
    controls: PistonControl,
    mut trial: impl FnMut(f64) -> Result<PistonStep, PistonError>,
) -> Result<PistonStep, PistonError> {
    let mut evaluate = |speed| {
        let value = trial(speed)?;
        if value.momentum_residual_n_s == 0.0 && !value.admitted(controls) {
            return Err(PistonError::Balance {
                momentum_residual_n_s: value.combined_momentum_residual_n_s.abs(),
                energy_residual_j: value
                    .mechanical_energy_residual_j
                    .abs()
                    .max(value.total_energy_residual_j.abs()),
            });
        }
        Ok(value)
    };
    let current = evaluate(initial)?;
    if current.admitted(controls) {
        return Ok(current);
    }
    let [mut lo, mut hi] = controls.midpoint_velocity_bracket_m_s;
    let low = evaluate(lo)?;
    if low.admitted(controls) {
        return Ok(low);
    }
    let high = evaluate(hi)?;
    if high.admitted(controls) {
        return Ok(high);
    }
    let mut flo = low.momentum_residual_n_s;
    let mut fhi = high.momentum_residual_n_s;
    if flo.is_sign_positive() == fhi.is_sign_positive() {
        return Err(PistonError::OutsideBracket {
            lower_impulse_n_s: flo,
            upper_impulse_n_s: fhi,
        });
    }
    if initial > lo && initial < hi {
        if current.momentum_residual_n_s.is_sign_positive() == flo.is_sign_positive() {
            lo = initial;
            flo = current.momentum_residual_n_s;
        } else {
            hi = initial;
            fhi = current.momentum_residual_n_s;
        }
    }
    let mut residual = high.momentum_residual_n_s;
    for iterations in 1..=controls.maximum_iterations {
        let scale = flo.abs().max(fhi.abs());
        let fraction = (flo.abs() / scale) / (flo.abs() / scale + fhi.abs() / scale);
        // Periodic bisection retains bracket contraction without forcing a
        // full field proposal far from a nearly converged secant every time.
        // Scaled weights avoid overflow in the secant's residual denominator.
        let mut speed = if iterations % 4 != 0 {
            (1.0 - fraction) * lo + fraction * hi
        } else {
            f64::midpoint(lo, hi)
        };
        if speed <= lo || speed >= hi {
            speed = f64::midpoint(lo, hi);
        }
        if speed <= lo || speed >= hi {
            return Err(PistonError::NotConverged {
                iterations: iterations - 1,
                impulse_residual_n_s: residual,
            });
        }
        let mut value = evaluate(speed)?;
        residual = value.momentum_residual_n_s;
        if value.admitted(controls) {
            value.iterations = iterations;
            return Ok(value);
        }
        if residual.is_sign_positive() == flo.is_sign_positive() {
            lo = speed;
            flo = residual;
        } else {
            hi = speed;
            fhi = residual;
        }
    }
    Err(PistonError::NotConverged {
        iterations: controls.maximum_iterations,
        impulse_residual_n_s: residual,
    })
}

#[cfg(test)]
mod tests;
