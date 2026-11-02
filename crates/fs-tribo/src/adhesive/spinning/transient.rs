//! Curved uniform-disk relaxation with physical-power quadrature and finite arrest.
//!
//! The original homogeneous disk law drives `q' = W*w(q)` in the complete
//! radius-scaled coordinates. Positive-stage Bogacki--Shampine 3(2) updates
//! integrate actual force, spin moment, travel and native sliding power. The
//! embedded difference controls numerical error; no numerical energy defect is
//! renamed frictional heat. An independent original-matrix work gate remains.
//!
//! With `W=L*L^T`, `r=|L^-1*q|` and kinetic force capacity C, the real disk obeys
//! `2*C*|q|/sqrt(13) <= H(q) <= sqrt(13)*C*|q|/3`. Thus positive computable
//! bounds alpha,beta give `-beta <= r' <= -alpha`: the remainder stops within
//! `[r/beta,r/alpha]`. Once this entire tail is before the requested horizon
//! and its time/travel bounds fit the declared accuracy, its total impulse is
//! exactly `-W^-1*q`, and its physical heat is `q dot W^-1*q/2`. Its unresolved
//! travel is returned as a bound, never fabricated. Bodies are not clipped.
//!
//! The stop interval bounds the analytic tail from the numerically resolved
//! prefix; embedded prefix errors are estimates, not validated-flow enclosures.
//! Area, strengths, mobility and geometry remain fixed, with no external drive.

mod metric;

use super::contact::{AdhesiveWrenchConfig, AdhesiveWrenchContact, AdhesiveWrenchError};
use super::relaxation::{
    AdhesiveRelaxationError, AdhesiveRelaxationInput, AdhesiveRelaxationReport,
};
use fs_math::det;
use metric::Metric;

/// Physical-time accuracy and bounded numerical work, independent of law budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTransientControls {
    /// Largest moving physical substep [s].
    pub maximum_step_s: f64,
    /// Smallest unresolved moving substep [s]; refusal replaces hidden clipping.
    pub minimum_step_s: f64,
    /// Accepted plus rejected moving attempts, including local work refusals.
    pub maximum_attempts: usize,
    /// Total original support evaluations, including the exact-ray probe.
    pub maximum_evaluations: usize,
    /// Accumulated estimated error in `[u_x,u_y,a*omega]` [m/s].
    pub scaled_rate_tolerance_m_s: f64,
    /// Accumulated estimated error in `[J_x,J_y,J_spin/a]` [N s].
    pub scaled_impulse_tolerance_ns: f64,
    /// Estimated quadrature error plus unresolved tail in `[du_x,du_y,a*dtheta]` [m].
    pub scaled_travel_tolerance_m: f64,
    /// Accumulated estimated error of original physical-power quadrature [J].
    pub heat_tolerance_j: f64,
    /// Largest admitted remaining stopping-time interval width [s].
    pub stopping_time_tolerance_s: f64,
}

/// Integrated original traction, resolved travel and explicitly bounded tail.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveTransientReport {
    /// Endpoint target; zero only after an exact ray event or finite extinction.
    pub velocity_after_m_s: [f64; 2],
    /// Endpoint relative spin target [rad/s].
    pub spin_rate_after_rad_s: f64,
    /// Integrated actual acting force [N s].
    pub force_impulse_ns: [f64; 2],
    /// Integrated actual acting free moment [N m s].
    pub spin_angular_impulse_nms: f64,
    /// Resolved centre travel; any uncomputed tail lies in the reported bounds.
    pub slip_increment_m: [f64; 2],
    /// Resolved relative rotation; any tail uncertainty is reported separately.
    pub spin_increment_rad: f64,
    /// Integrated original sliding/spin heat, with no numerical-loss addition.
    pub dissipated_energy_j: f64,
    /// Original physical-power quadrature before any analytic extinction tail.
    pub prefix_dissipated_energy_j: f64,
    /// Exact remaining relative kinetic energy consumed by the extinction tail.
    pub stopping_tail_heat_j: f64,
    /// Absorbed work from the unchanged full matrix and integrated impulse [J].
    pub contact_work_j: f64,
    /// Full physical heat minus original-matrix absorbed work [J].
    pub work_residual_j: f64,
    /// Independent prefix heat/work residual, before analytical tail closure [J].
    pub prefix_work_residual_j: f64,
    /// Sum of embedded local rate-error estimates in physical scaled units.
    pub estimated_scaled_rate_error_m_s: f64,
    /// Sum of embedded local integrated-wrench error estimates [N s].
    pub estimated_scaled_impulse_error_ns: f64,
    /// Sum of embedded local travel-error estimates [m].
    pub estimated_scaled_travel_error_m: f64,
    /// Sum of embedded local physical-power quadrature error estimates [J].
    pub estimated_heat_error_j: f64,
    /// Componentwise absolute uncomputed centre travel bounds [m].
    pub tail_slip_bound_m: [f64; 2],
    /// Absolute uncomputed relative-spin bound [rad].
    pub tail_spin_bound_rad: f64,
    /// Stop interval within this horizon; equal ends only for an exact ray.
    pub stopping_time_bracket_s: Option<[f64; 2]>,
    /// Physical time explicitly integrated before a tail, or the entire horizon.
    pub integrated_until_s: f64,
    /// Original-matrix endpoint minus returned endpoint target [m/s].
    pub scaled_velocity_residual_m_s: [f64; 3],
    /// Original exact-ray path used instead of numerical substeps.
    pub exact_ray: bool,
    /// Accepted physical moving substeps, excluding analytic tail/rest.
    pub accepted_steps: usize,
    /// Rejected local numerical candidates; no physical work is accumulated.
    pub rejected_steps: usize,
    /// Original support evaluations including rejected candidates.
    pub evaluations: usize,
}

/// Numerical and physical refusals preserve the complete original input.
#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveTransientError {
    /// Original exact-ray/law admission; curved motion alone uses this integrator.
    Exact(AdhesiveRelaxationError),
    /// Original matrix, material or effort-control admission.
    Contact(AdhesiveWrenchError),
    /// Invalid controls or unrepresentable finite mechanics.
    Invalid(&'static str),
    /// Original instantaneous force or moment exceeded the declared budget.
    EffortBudget,
    /// Unknown tail directions cannot be covered by the declared effort budget.
    TailEffortBudget,
    /// Accepted plus rejected attempts exhausted their declared bound.
    AttemptBudget,
    /// Original support evaluations exhausted their declared bound.
    EvaluationBudget,
    /// Further physical-time refinement is not representable or admitted.
    StepUnderflow,
    /// Original independent prefix or full work/endpoint gates did not close.
    UnresolvedWork,
    /// Cooperative cancellation observed before publication.
    Cancelled,
}
impl From<AdhesiveWrenchError> for AdhesiveTransientError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for AdhesiveTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive transient: {self:?}")
    }
}
impl std::error::Error for AdhesiveTransientError {}
type Error = AdhesiveTransientError;
type Vector = [f64; 3];

impl AdhesiveWrenchContact {
    /// Integrate curved passive slip/spin, retaining the existing exact-ray fast path.
    /// A remaining tail is accepted only when its complete stop interval precedes
    /// the horizon and its unresolved travel fits the declared accuracy.
    #[allow(clippy::too_many_lines)]
    pub fn advance_transient(
        &self,
        input: &AdhesiveRelaxationInput,
        config: AdhesiveWrenchConfig,
        controls: AdhesiveTransientControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdhesiveTransientReport, Error> {
        poll(&mut cancelled)?;
        controls.validate()?;
        match self.relax_velocity(input, config, &mut cancelled) {
            Ok(exact) => {
                let residual = norm(exact.scaled_velocity_residual_m_s);
                if !residual.is_finite() || residual > controls.scaled_rate_tolerance_m_s {
                    return Err(Error::UnresolvedWork);
                }
                poll(&mut cancelled)?;
                return Ok(from_exact(exact, input.dt_s));
            }
            Err(AdhesiveRelaxationError::CurvedTrajectory { .. }) => {}
            Err(AdhesiveRelaxationError::Cancelled) => return Err(Error::Cancelled),
            Err(error) => return Err(Error::Exact(error)),
        }
        let radius = self.parameters().contact_radius_m;
        let matrix = input.scaled_mobility_per_kg;
        let initial = [
            input.velocity_m_s[0],
            input.velocity_m_s[1],
            radius * input.spin_rate_rad_s,
        ];
        let metric = Metric::new(matrix, self.kinetic_patch().sliding_capacity_n(), config)?;
        let horizon = input.dt_s.min(metric.tail(initial)?.0[1]);
        let mut work = Work {
            law: self,
            radius,
            config,
            controls,
            evaluations: 1,
            cancelled: &mut cancelled,
        };
        let mut impulse = [0.0; 3];
        let mut travel = [0.0; 3];
        let mut heat = 0.0;
        let mut errors = [0.0; 4];
        let mut elapsed = 0.0;
        let mut next_step = controls.maximum_step_s;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut cached = None;
        for _ in 0..controls.maximum_attempts {
            poll(work.cancelled)?;
            let q = add(initial, multiply(matrix, impulse));
            let remaining = input.dt_s - elapsed;
            let (stop, tail_travel) = metric.tail(q)?;
            if remaining == 0.0 {
                let candidate = finish(
                    input,
                    config,
                    controls,
                    radius,
                    impulse,
                    travel,
                    heat,
                    0.0,
                    errors,
                    [0.0; 3],
                    None,
                    elapsed,
                    accepted,
                    rejected,
                    work.evaluations,
                )?;
                poll(work.cancelled)?;
                return Ok(candidate);
            }
            if stop[1] <= remaining
                && stop[1] - stop[0] <= controls.stopping_time_tolerance_s
                && errors[2] + norm(tail_travel) <= controls.scaled_travel_tolerance_m
            {
                if self.kinetic_patch().sliding_capacity_n() > config.maximum_force_n
                    || self.kinetic_patch().spinning_capacity_nm() > config.maximum_spin_moment_nm
                {
                    return Err(Error::TailEffortBudget);
                }
                let prefix_work = absorbed(initial, matrix, impulse);
                if !prefix_work.is_finite()
                    || (heat - prefix_work).abs() > config.energy_tolerance_j
                {
                    return Err(Error::UnresolvedWork);
                }
                let tail_impulse = metric.stop_impulse(q)?;
                let tail_heat = finite(-0.5 * dot(q, tail_impulse))?;
                if tail_heat < 0.0 {
                    return Err(Error::UnresolvedWork);
                }
                let candidate = finish(
                    input,
                    config,
                    controls,
                    radius,
                    add(impulse, tail_impulse),
                    travel,
                    heat,
                    tail_heat,
                    errors,
                    tail_travel,
                    Some([elapsed + stop[0], elapsed + stop[1]]),
                    elapsed,
                    accepted,
                    rejected,
                    work.evaluations,
                )?;
                poll(work.cancelled)?;
                return Ok(AdhesiveTransientReport {
                    prefix_work_residual_j: heat - prefix_work,
                    ..candidate
                });
            }
            // Every positive embedded stage stays within h*beta of its initial
            // mass-metric state, strictly before any possible zero crossing.
            let safe_step = 0.25 * metric.radius(q)?.0 / metric.maximum_speed;
            let h = next_step.min(safe_step).min(remaining);
            if !h.is_finite()
                || h <= 0.0
                || elapsed + h == elapsed
                || (h < controls.minimum_step_s && h != remaining)
            {
                return Err(Error::StepUnderflow);
            }
            let first = match cached {
                Some(value) => value,
                None => work.sample(q)?,
            };
            let step = step(&mut work, matrix, q, first, h)?;
            let local_errors = [
                metric.rate_error(step.impulse_error)?,
                norm(step.impulse_error),
                norm(step.travel_error),
                step.heat_error.abs(),
            ];
            let limits = [
                controls.scaled_rate_tolerance_m_s,
                controls.scaled_impulse_tolerance_ns,
                controls.scaled_travel_tolerance_m,
                controls.heat_tolerance_j,
            ];
            let fraction = h / horizon;
            let ratio = local_errors
                .into_iter()
                .zip(limits)
                .map(|(e, t)| ratio(e, t * fraction))
                .fold(0.0, f64::max)
                .max(ratio(
                    (step.heat - absorbed(q, matrix, step.impulse)).abs(),
                    config.energy_tolerance_j * fraction,
                ));
            if !ratio.is_finite() || ratio > 1.0 || step.heat < 0.0 {
                rejected += 1;
                next_step = 0.5 * h;
                cached = Some(first);
                continue;
            }
            impulse = add(impulse, step.impulse);
            travel = add(travel, step.travel);
            heat = finite(heat + step.heat)?;
            for i in 0..4 {
                errors[i] = finite(errors[i] + local_errors[i])?;
            }
            elapsed = if h == remaining {
                input.dt_s
            } else {
                elapsed + h
            };
            accepted += 1;
            if elapsed == input.dt_s {
                let candidate = finish(
                    input,
                    config,
                    controls,
                    radius,
                    impulse,
                    travel,
                    heat,
                    0.0,
                    errors,
                    [0.0; 3],
                    None,
                    elapsed,
                    accepted,
                    rejected,
                    work.evaluations,
                )?;
                poll(work.cancelled)?;
                return Ok(candidate);
            }
            // Recompute at the actual accumulated-impulse endpoint. Reusing a
            // different rounded stage endpoint would obscure native work closure.
            cached = None;
            next_step = if ratio < 0.125 {
                (2.0 * h).min(controls.maximum_step_s)
            } else {
                h
            };
        }
        Err(Error::AttemptBudget)
    }
}

#[derive(Clone, Copy)]
struct Sample {
    wrench: Vector,
    velocity: Vector,
    power: f64,
}
struct Work<'a, C> {
    law: &'a AdhesiveWrenchContact,
    radius: f64,
    config: AdhesiveWrenchConfig,
    controls: AdhesiveTransientControls,
    evaluations: usize,
    cancelled: &'a mut C,
}
impl<C: FnMut() -> bool> Work<'_, C> {
    fn sample(&mut self, q: Vector) -> Result<Sample, Error> {
        poll(self.cancelled)?;
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(Error::EvaluationBudget);
        }
        self.evaluations += 1;
        let response = self
            .law
            .kinetic_patch()
            .evaluate([q[0], q[1]], finite(q[2] / self.radius)?)
            .map_err(AdhesiveWrenchError::Law)?;
        if det::hypot(response.force_n[0], response.force_n[1]) > self.config.maximum_force_n
            || response.spin_moment_nm.abs() > self.config.maximum_spin_moment_nm
        {
            return Err(Error::EffortBudget);
        }
        if response.at_rest {
            return Err(Error::StepUnderflow);
        }
        Ok(Sample {
            wrench: [
                response.force_n[0],
                response.force_n[1],
                response.spin_moment_nm / self.radius,
            ],
            velocity: q,
            power: response.dissipated_power_w,
        })
    }
}
struct Step {
    impulse: Vector,
    impulse_error: Vector,
    travel: Vector,
    travel_error: Vector,
    heat: f64,
    heat_error: f64,
}
fn step<C: FnMut() -> bool>(
    work: &mut Work<'_, C>,
    matrix: [[f64; 3]; 3],
    q: Vector,
    a: Sample,
    h: f64,
) -> Result<Step, Error> {
    let b = work.sample(add(q, multiply(matrix, a.wrench.map(|v| 0.5 * h * v))))?;
    let c = work.sample(add(q, multiply(matrix, b.wrench.map(|v| 0.75 * h * v))))?;
    let combine =
        |values: [f64; 3]| h * (2.0 * values[0] / 9.0 + values[1] / 3.0 + 4.0 * values[2] / 9.0);
    let impulse = core::array::from_fn(|i| combine([a.wrench[i], b.wrench[i], c.wrench[i]]));
    let d = work.sample(add(q, multiply(matrix, impulse)))?;
    // Third minus second order weights: [-5/72, 1/12, 1/9, -1/8].
    let error = |values: [f64; 4]| {
        h * (-5.0 * values[0] / 72.0 + values[1] / 12.0 + values[2] / 9.0 - values[3] / 8.0)
    };
    Ok(Step {
        impulse,
        impulse_error: core::array::from_fn(|i| {
            error([a.wrench[i], b.wrench[i], c.wrench[i], d.wrench[i]])
        }),
        travel: core::array::from_fn(|i| combine([a.velocity[i], b.velocity[i], c.velocity[i]])),
        travel_error: core::array::from_fn(|i| {
            error([a.velocity[i], b.velocity[i], c.velocity[i], d.velocity[i]])
        }),
        heat: finite(combine([a.power, b.power, c.power]))?,
        heat_error: finite(error([a.power, b.power, c.power, d.power]))?,
    })
}

#[allow(clippy::too_many_arguments)]
fn finish(
    input: &AdhesiveRelaxationInput,
    config: AdhesiveWrenchConfig,
    controls: AdhesiveTransientControls,
    radius: f64,
    impulse: Vector,
    travel: Vector,
    prefix_heat: f64,
    tail_heat: f64,
    errors: [f64; 4],
    tail_travel: Vector,
    stop: Option<[f64; 2]>,
    elapsed: f64,
    accepted: usize,
    rejected: usize,
    evaluations: usize,
) -> Result<AdhesiveTransientReport, Error> {
    let initial = [
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        radius * input.spin_rate_rad_s,
    ];
    let actual = add(initial, multiply(input.scaled_mobility_per_kg, impulse));
    let endpoint = if stop.is_some() { [0.0; 3] } else { actual };
    let heat = finite(prefix_heat + tail_heat)?;
    let contact_work = finite(absorbed(initial, input.scaled_mobility_per_kg, impulse))?;
    let residual = finite(heat - contact_work)?;
    if residual.abs() > config.energy_tolerance_j
        || norm(add(actual, endpoint.map(|v| -v))) > controls.scaled_rate_tolerance_m_s
        || errors[0] > controls.scaled_rate_tolerance_m_s
        || errors[1] > controls.scaled_impulse_tolerance_ns
        || errors[2] + norm(tail_travel) > controls.scaled_travel_tolerance_m
        || errors[3] > controls.heat_tolerance_j
    {
        return Err(Error::UnresolvedWork);
    }
    Ok(AdhesiveTransientReport {
        velocity_after_m_s: [endpoint[0], endpoint[1]],
        spin_rate_after_rad_s: finite(endpoint[2] / radius)?,
        force_impulse_ns: [finite(impulse[0])?, finite(impulse[1])?],
        spin_angular_impulse_nms: finite(radius * impulse[2])?,
        slip_increment_m: [finite(travel[0])?, finite(travel[1])?],
        spin_increment_rad: finite(travel[2] / radius)?,
        dissipated_energy_j: heat,
        prefix_dissipated_energy_j: prefix_heat,
        stopping_tail_heat_j: tail_heat,
        contact_work_j: contact_work,
        work_residual_j: residual,
        prefix_work_residual_j: residual,
        estimated_scaled_rate_error_m_s: errors[0],
        estimated_scaled_impulse_error_ns: errors[1],
        estimated_scaled_travel_error_m: errors[2],
        estimated_heat_error_j: errors[3],
        tail_slip_bound_m: [tail_travel[0], tail_travel[1]],
        tail_spin_bound_rad: finite(tail_travel[2] / radius)?,
        stopping_time_bracket_s: stop,
        integrated_until_s: elapsed,
        scaled_velocity_residual_m_s: add(actual, endpoint.map(|v| -v)),
        exact_ray: false,
        accepted_steps: accepted,
        rejected_steps: rejected,
        evaluations,
    })
}

fn from_exact(r: AdhesiveRelaxationReport, dt: f64) -> AdhesiveTransientReport {
    AdhesiveTransientReport {
        velocity_after_m_s: r.velocity_after_m_s,
        spin_rate_after_rad_s: r.spin_rate_after_rad_s,
        force_impulse_ns: r.force_impulse_ns,
        spin_angular_impulse_nms: r.spin_angular_impulse_nms,
        slip_increment_m: r.slip_increment_m,
        spin_increment_rad: r.spin_increment_rad,
        dissipated_energy_j: r.dissipated_energy_j,
        prefix_dissipated_energy_j: r.dissipated_energy_j,
        stopping_tail_heat_j: 0.0,
        contact_work_j: r.contact_work_j,
        work_residual_j: r.work_residual_j,
        prefix_work_residual_j: r.work_residual_j,
        estimated_scaled_rate_error_m_s: 0.0,
        estimated_scaled_impulse_error_ns: 0.0,
        estimated_scaled_travel_error_m: 0.0,
        estimated_heat_error_j: 0.0,
        tail_slip_bound_m: [0.0; 2],
        tail_spin_bound_rad: 0.0,
        stopping_time_bracket_s: r.kinetic_phase_end_s.map(|t| [t, t]),
        integrated_until_s: r.kinetic_phase_end_s.unwrap_or(dt),
        scaled_velocity_residual_m_s: r.scaled_velocity_residual_m_s,
        exact_ray: true,
        accepted_steps: 0,
        rejected_steps: 0,
        evaluations: 1,
    }
}
impl AdhesiveTransientControls {
    fn validate(self) -> Result<(), Error> {
        if !self.maximum_step_s.is_finite()
            || !self.minimum_step_s.is_finite()
            || self.minimum_step_s <= 0.0
            || self.maximum_step_s < self.minimum_step_s
            || self.maximum_attempts == 0
            || self.maximum_evaluations == 0
            || [
                self.scaled_rate_tolerance_m_s,
                self.scaled_impulse_tolerance_ns,
                self.scaled_travel_tolerance_m,
                self.heat_tolerance_j,
                self.stopping_time_tolerance_s,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(Error::Invalid(
                "finite positive transient tolerances and bounded work required",
            ));
        }
        Ok(())
    }
}
fn ratio(error: f64, tolerance: f64) -> f64 {
    if error == 0.0 { 0.0 } else { error / tolerance }
}
fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: Vector) -> f64 {
    det::hypot(det::hypot(a[0], a[1]), a[2])
}
fn add(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] + b[i])
}
fn multiply(matrix: [[f64; 3]; 3], vector: Vector) -> Vector {
    matrix.map(|row| dot(row, vector))
}
fn absorbed(q: Vector, matrix: [[f64; 3]; 3], impulse: Vector) -> f64 {
    -dot(impulse, add(q, multiply(matrix, impulse).map(|v| 0.5 * v)))
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("unrepresentable adhesive transient"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
