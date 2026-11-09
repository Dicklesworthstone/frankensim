//! Driven physical time under the original local-speed Hertz pressure law.
//!
//! Positive RK4 step doubling integrates `q'=a+W*w(q)`, actual pressure-field
//! power and travel. It follows the initial-value branch of weakening laws;
//! no nonlinear algebraic endpoint root or effective patch coefficient replaces
//! that evolution. Spatial pressure quadrature remains an independent axis.
//!
//! Rest is retained explicitly. The original joint static Hertz surface holds
//! or releases it; the first moving stage uses the original zero-speed kinetic
//! right limit, not the patch evaluator's zero-at-rest selection. A separating
//! support plane proves positive departure over that first interval. Later
//! intervals exclude zero by their original-matrix speed bound. Nonconstant-law
//! arrest that cannot be resolved this way refuses explicitly: passive stopping
//! tails, speed capture and fabricated terminal heat are never used with drive.
//! Genuine Coulomb delegates the existing driven owner, including its exact
//! rays and explicitly first-order curved events.
//! Original affine invariant translation/spin additionally has analytic event,
//! travel and heat integration, including static dwell and overdriven reversal.

use core::f64::consts::{PI, SQRT_2};
use fs_math::det;

use super::transient::{HertzRateTransientError, coefficient_bounds};
use super::{HertzFrictionPatch, HertzFrictionQuadrature};
use crate::adhesive::spinning::transient::metric::Metric;
use crate::adhesive::spinning::transient::{Sample, TransientError, doubled_rk4_with_drive};
use crate::sliding_patch::contact::driven::{
    HertzDrivenControls, HertzDrivenError, HertzDrivenInput, HertzDrivenMode, HertzDrivenReport,
};
use crate::sliding_patch::contact::transient::{HertzTransientControls, HertzTransientError};
use crate::sliding_patch::contact::{
    Factor, HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput,
    HertzWrenchMode, HertzWrenchParameters, HertzWrenchPhase, factor_circular_response,
};
use crate::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};
use crate::{FrictionLaw, TriboError};

type V = [f64; 3];
type M = [[f64; 3]; 3];
const K: f64 = 3.0 * PI / 16.0;
const PAD: f64 = 1024.0 * f64::EPSILON;

mod invariant;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzRateDrivenControls {
    /// Adaptive physical-time errors/budgets. The stopping-time tolerance does
    /// not authorize an unresolved event; admitted affine events are analytic.
    pub integration: HertzTransientControls,
    /// Admission only for an explicitly retained resting constraint [m/s].
    pub rest_scaled_rate_tolerance_m_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HertzRateDrivenReport {
    pub trajectory: HertzDrivenReport,
    /// Analytic invariant affine solution, separately from adaptive estimates.
    pub exact_invariant: bool,
    /// First resolved arrest interval. Analytic events give the same computed
    /// time at both ends; this is not a validated floating-point enclosure.
    pub stopping_time_bracket_s: Option<[f64; 2]>,
    pub local_speed_bound_m_s: f64,
    pub kinetic_coefficient_bounds: [f64; 2],
    /// Includes original static support as well as the moving pressure bound.
    pub maximum_scaled_wrench_norm_n: f64,
    /// False for analytic invariant motion and delegated Coulomb: zero fields
    /// do not assert measured adaptive error. Coulomb curved continuation still
    /// retains its refinement contract; `exact_invariant` identifies this path.
    pub adaptive_error_estimates_available: bool,
    pub estimated_scaled_rate_error_m_s: f64,
    pub estimated_scaled_impulse_error_ns: f64,
    pub estimated_scaled_travel_error_m: f64,
    pub estimated_heat_error_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HertzRateDrivenError {
    Law(TriboError),
    Driven(HertzDrivenError),
    Transient(HertzTransientError),
    /// No admitted nonzero continuation or proved driven arrest remains.
    /// Caller state, including contact mode, has not been published.
    UnresolvedDrivenEvent {
        elapsed_s: f64,
        scaled_rate_m_s: [f64; 3],
    },
}
impl From<TriboError> for HertzRateDrivenError {
    fn from(e: TriboError) -> Self {
        Self::Law(e)
    }
}
impl From<HertzDrivenError> for HertzRateDrivenError {
    fn from(e: HertzDrivenError) -> Self {
        Self::Driven(e)
    }
}
impl From<HertzWrenchError> for HertzRateDrivenError {
    fn from(e: HertzWrenchError) -> Self {
        Self::Driven(e.into())
    }
}
impl From<TransientError> for HertzRateDrivenError {
    fn from(e: TransientError) -> Self {
        Self::Transient(e.with_exact())
    }
}
impl From<HertzRateTransientError> for HertzRateDrivenError {
    fn from(e: HertzRateTransientError) -> Self {
        match e {
            HertzRateTransientError::Law(e) => Self::Law(e),
            HertzRateTransientError::Transient(e) => Self::Transient(e),
            HertzRateTransientError::NonCoercive { .. } => {
                invalid("unrepresentable coefficient bound")
            }
        }
    }
}
impl core::fmt::Display for HertzRateDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven rate-dependent Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for HertzRateDrivenError {}
type Error = HertzRateDrivenError;

impl HertzFrictionPatch {
    /// Follow the original moving initial-value branch with held drive, source,
    /// normal load, geometry and full SPD mobility. No candidate is published
    /// before the original impulse/work, numerical error and cancellation gates.
    #[allow(clippy::too_many_lines)]
    pub fn advance_driven(
        &self,
        input: &HertzDrivenInput,
        config: HertzWrenchConfig,
        quadrature: HertzFrictionQuadrature,
        controls: HertzRateDrivenControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzRateDrivenReport, Error> {
        poll(&mut cancelled)?;
        config.validate()?;
        quadrature.validate()?;
        controls.integration.validate()?;
        if !controls.rest_scaled_rate_tolerance_m_s.is_finite()
            || controls.rest_scaled_rate_tolerance_m_s < 0.0
        {
            return Err(invalid("finite retained-rest tolerance required"));
        }
        let radius = self.contact_radius_m();
        let normal = self.normal_force_n();
        let zero = self.law.zero_speed_kinetic_coefficient()?;
        let static_mu = self.law.static_coefficient()?;
        let limit = HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n: normal,
            contact_radius_m: radius,
            static_mu,
            kinetic_mu: zero,
        })?;
        let upper_norm = det::sqrt(1.0 + K * K);
        if matches!(self.law, FrictionLaw::Coulomb { .. }) || normal == 0.0 {
            let c = controls.integration;
            let trajectory = limit.advance_driven(
                input,
                config,
                HertzDrivenControls {
                    maximum_step_s: c.maximum_step_s,
                    minimum_step_s: c.minimum_step_s,
                    maximum_attempts: c.maximum_attempts,
                    maximum_evaluations: c.maximum_evaluations,
                    scaled_rate_tolerance_m_s: c.scaled_rate_tolerance_m_s,
                    rest_scaled_rate_tolerance_m_s: controls.rest_scaled_rate_tolerance_m_s,
                },
                &mut cancelled,
            )?;
            let (speed_bound, bounds) = if matches!(self.law, FrictionLaw::Coulomb { .. }) {
                (f64::INFINITY, [zero; 2])
            } else {
                // An unloaded curve has zero traction but retains its original
                // nonconstant coefficient. The free path is affine in time;
                // ||u||+|a*omega| is convex, so its endpoint maximum bounds it.
                let start = [
                    input.velocity_m_s[0],
                    input.velocity_m_s[1],
                    radius * input.spin_rate_rad_s,
                ];
                let acceleration = [
                    input.driving_acceleration_m_s2[0],
                    input.driving_acceleration_m_s2[1],
                    radius * input.driving_spin_acceleration_rad_s2,
                ];
                let end = vector(add(start, scale(acceleration, input.dt_s)))?;
                let speed = checked(
                    (det::hypot(start[0], start[1]) + start[2].abs())
                        .max(det::hypot(end[0], end[1]) + end[2].abs())
                        * (1.0 + PAD),
                )?;
                (speed, coefficient_bounds(&self.law, speed)?)
            };
            return Ok(report(
                trajectory,
                speed_bound,
                bounds,
                checked(normal * static_mu * upper_norm)?,
                false,
                [0.0; 4],
            ));
        }
        if !input.dt_s.is_finite() || input.dt_s <= 0.0 {
            return Err(invalid("positive finite driven duration required"));
        }
        let original = vector([
            input.velocity_m_s[0],
            input.velocity_m_s[1],
            radius * input.spin_rate_rad_s,
        ])?;
        let drive = vector([
            input.driving_acceleration_m_s2[0],
            input.driving_acceleration_m_s2[1],
            radius * input.driving_spin_acceleration_rad_s2,
        ])?;
        if !input.spin_rate_rad_s.is_finite()
            || !input.driving_spin_acceleration_rad_s2.is_finite()
            || (input.spin_rate_rad_s != 0.0 && original[2] == 0.0)
            || (input.driving_spin_acceleration_rad_s2 != 0.0 && drive[2] == 0.0)
        {
            return Err(invalid("unrepresentable scaled spin port"));
        }
        let initial = match input.mode {
            HertzDrivenMode::Resting => {
                if norm(original) > controls.rest_scaled_rate_tolerance_m_s {
                    return Err(invalid("retained rest exceeds original rate gate"));
                }
                [0.0; 3]
            }
            HertzDrivenMode::Sliding => {
                if original == [0.0; 3] {
                    return Err(invalid("zero rate requires explicit resting mode"));
                }
                original
            }
        };
        let matrix = input.scaled_mobility_per_kg;
        let metric = Metric::new(matrix, 1.0, [0.0, 1.0], config)?;
        let factor = factor_circular_response(matrix, config)?;
        let demand = solve(factor, drive)?;
        let r0 = metric.radius(original)?.1;
        let drive_radius = metric.radius(drive)?.1;
        // Exact flow obeys r<=r0+D*t. Numerical positive RK stages need a
        // separate work allowance: r^2<=r0^2+2*D*T*(1.25*R)+2*epsilon.
        // The positive root below admits every endpoint; every stage is checked
        // against 1.25*R. No final residual cancellation hides an excursion.
        let forcing = checked(1.25 * input.dt_s * drive_radius)?;
        let mass_bound = checked(
            (forcing
                + det::hypot(
                    det::hypot(r0, forcing),
                    SQRT_2 * det::sqrt(config.energy_tolerance_j),
                ))
                * (1.0 + PAD),
        )?;
        let speed_bound = checked(1.25 * SQRT_2 * metric.maximum_speed * mass_bound * (1.0 + PAD))?;
        let bounds = coefficient_bounds(&self.law, speed_bound)?;
        let moving_wrench_bound = checked(normal * bounds[1] * upper_norm)?;
        let maximum_wrench =
            checked(moving_wrench_bound.max(normal * static_mu * upper_norm) * (1.0 + PAD))?;
        let slope = match self.law {
            FrictionLaw::VelocityDependent {
                slope_per_speed, ..
            } => slope_per_speed.max(0.0),
            FrictionLaw::Stribeck {
                viscous_per_speed, ..
            } => viscous_per_speed,
            _ => 0.0,
        };
        let affine_bounds = [
            checked(normal * zero * upper_norm * (1.0 + PAD))?,
            checked(normal * slope * (1.0 + PAD))?,
        ];
        if let Some(exact) = invariant::advance(
            self,
            input,
            original,
            initial,
            drive,
            demand,
            factor,
            config,
            quadrature,
            controls,
            &limit,
            &mut cancelled,
        )? {
            return Ok(exact);
        }
        let mut evaluations = 0;
        let mut projection_iterations = 0;
        let mut release = None;
        let mut first = None;
        let mut attempts = 0;
        if input.mode == HertzDrivenMode::Resting {
            attempts += 1;
            let mut probe_config = config;
            probe_config.maximum_evaluations = probe_config
                .maximum_evaluations
                .min(controls.integration.maximum_evaluations);
            let support = limit.advance_phase(
                &HertzWrenchInput {
                    free_displacement_m: [drive[0], drive[1]],
                    free_spin_angle_rad: checked(drive[2] / radius)?,
                    scaled_compliance_m_per_n: matrix,
                },
                HertzWrenchPhase::StaticAdmissible,
                probe_config,
                &mut cancelled,
            )?;
            evaluations += support.evaluations;
            projection_iterations += support.iterations;
            let w = vector([
                support.force_n[0],
                support.force_n[1],
                support.spin_moment_nm / radius,
            ])?;
            if support.mode == HertzWrenchMode::Sticking || drive == [0.0; 3] {
                let trajectory = finish(
                    input,
                    original,
                    [0.0; 3],
                    drive,
                    demand,
                    factor,
                    radius,
                    scale(w, input.dt_s),
                    [0.0; 3],
                    0.0,
                    [0.0; 4],
                    HertzDrivenMode::Resting,
                    config,
                    controls,
                    1,
                    0,
                    evaluations,
                    projection_iterations,
                )?;
                poll(&mut cancelled)?;
                return Ok(report(
                    trajectory,
                    speed_bound,
                    bounds,
                    maximum_wrench,
                    true,
                    [0.0; 4],
                ));
            }
            let direction = vector([
                support.slip_displacement_m[0],
                support.slip_displacement_m[1],
                radius * support.spin_angle_rad,
            ])?;
            if direction == [0.0; 3] {
                return Err(event(0.0, initial));
            }
            // One original right-limit support and one unit-pressure support
            // are explicit counted probes, not zero-speed regularization.
            if controls
                .integration
                .maximum_evaluations
                .saturating_sub(evaluations)
                < 2
            {
                return Err(TransientError::EvaluationBudget.into());
            }
            evaluations += 2;
            let right = self
                .zero_speed_patch
                .evaluate([direction[0], direction[1]], direction[2] / radius)?;
            if det::hypot(right.force_n[0], right.force_n[1]) > config.maximum_force_n
                || right.spin_moment_nm.abs() > config.maximum_spin_moment_nm
            {
                return Err(TransientError::EffortBudget.into());
            }
            let unit = HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: normal,
                contact_radius_m: radius,
                friction_coefficient: 1.0,
            })?
            .evaluate([direction[0], direction[1]], direction[2] / radius)?;
            let separator = solve(factor, direction)?;
            let projected_drive = dot(separator, drive);
            let allowed = checked(0.5 * (zero + projected_drive / unit.dissipated_power_w))?;
            if allowed <= zero || unit.dissipated_power_w <= 0.0 {
                return Err(event(0.0, initial));
            }
            let low_speed = if slope > 0.0 {
                speed_bound.min(checked(0.5 * (allowed - zero) / slope)?)
            } else {
                speed_bound
            };
            if low_speed <= 0.0 || coefficient_bounds(&self.law, low_speed)?[1] >= allowed {
                return Err(event(0.0, initial));
            }
            release = Some((low_speed, separator));
            first = Some(Sample {
                wrench: [
                    right.force_n[0],
                    right.force_n[1],
                    right.spin_moment_nm / radius,
                ],
                velocity: [0.0; 3],
                power: 0.0,
            });
            poll(&mut cancelled)?;
        }
        let mut work = Work {
            law: self,
            quadrature,
            config,
            controls,
            metric: &metric,
            radius,
            speed_bound,
            mass_bound,
            moving_wrench_bound,
            affine_bounds,
            release,
            evaluations,
            cancelled: &mut cancelled,
        };
        let mut impulse = [0.0; 3];
        let mut travel = [0.0; 3];
        let mut heat = 0.0;
        let mut errors = [0.0; 4];
        let mut elapsed = 0.0;
        let mut next = controls.integration.maximum_step_s;
        let mut accepted = 0;
        let mut rejected = 0;
        while attempts < controls.integration.maximum_attempts {
            attempts += 1;
            poll(work.cancelled)?;
            let q = vector(add(
                add(initial, scale(drive, elapsed)),
                multiply(matrix, impulse),
            ))?;
            if metric.radius(q)?.1 > mass_bound {
                return Err(invalid(
                    "driven endpoint left its admitted mass-radius domain",
                ));
            }
            let remaining = input.dt_s - elapsed;
            let global_beta = checked(drive_radius + metric.maximum_speed * moving_wrench_bound)?;
            let safe = if let Some((low_speed, _)) = work.release {
                if global_beta == 0.0 {
                    return Err(event(elapsed, q));
                }
                0.5 * low_speed / (SQRT_2 * metric.maximum_speed * global_beta)
            } else {
                let (lower, upper) = metric.radius(q)?;
                let local = checked(
                    metric.maximum_speed * affine_bounds[0]
                        + 1.25
                            * metric.maximum_speed
                            * metric.maximum_speed
                            * affine_bounds[1]
                            * upper,
                )?;
                let beta = checked(
                    drive_radius + (metric.maximum_speed * moving_wrench_bound).min(local),
                )?;
                if beta == 0.0 {
                    f64::INFINITY
                } else {
                    0.25 * lower / beta
                }
            };
            let h = next.min(safe).min(remaining);
            if !h.is_finite()
                || h <= 0.0
                || elapsed + h == elapsed
                || (h < controls.integration.minimum_step_s && h != remaining)
            {
                return Err(event(elapsed, q));
            }
            let sample = match first {
                Some(s) => s,
                None => work.sample(q)?,
            };
            let step = match doubled_rk4_with_drive(
                &mut |q| work.sample(q),
                matrix,
                q,
                drive,
                sample,
                h,
            ) {
                Ok(step) => step,
                Err(TransientError::StepUnderflow) => return Err(event(elapsed, q)),
                Err(error) => return Err(error.into()),
            };
            let local_errors = [
                metric.rate_error(step.impulse_error)?,
                norm(step.impulse_error),
                norm(step.travel_error),
                step.heat_error.abs(),
            ];
            let limits = error_limits(controls);
            let fraction = h / input.dt_s;
            let change = vector(add(scale(drive, h), multiply(matrix, step.impulse)))?;
            let kinetic = checked(dot(solve(factor, change)?, add(q, scale(change, 0.5))))?;
            let balance = checked(kinetic + step.heat - dot(demand, step.travel))?;
            let ratio = local_errors
                .into_iter()
                .zip(limits)
                .map(|(e, l)| e / (l * fraction))
                .fold(0.0_f64, f64::max)
                .max(balance.abs() / (config.energy_tolerance_j * fraction));
            if !ratio.is_finite() || ratio > 1.0 || step.heat < 0.0 {
                rejected += 1;
                next = 0.5 * h;
                first = Some(sample);
                continue;
            }
            let new_impulse = vector(add(impulse, step.impulse))?;
            let new_elapsed = if h == remaining {
                input.dt_s
            } else {
                checked(elapsed + h)?
            };
            let after = vector(add(
                add(initial, scale(drive, new_elapsed)),
                multiply(matrix, new_impulse),
            ))?;
            if after == [0.0; 3] {
                return Err(event(elapsed, q));
            }
            if metric.radius(after)?.1 > mass_bound {
                return Err(invalid(
                    "driven endpoint left its admitted mass-radius domain",
                ));
            }
            impulse = new_impulse;
            travel = vector(add(travel, step.travel))?;
            heat = checked(heat + step.heat)?;
            for i in 0..4 {
                errors[i] = checked(errors[i] + local_errors[i])?;
            }
            elapsed = new_elapsed;
            accepted += 1;
            first = None;
            work.release = None;
            if elapsed == input.dt_s {
                let trajectory = finish(
                    input,
                    original,
                    after,
                    drive,
                    demand,
                    factor,
                    radius,
                    impulse,
                    travel,
                    heat,
                    errors,
                    HertzDrivenMode::Sliding,
                    config,
                    controls,
                    accepted,
                    rejected,
                    work.evaluations,
                    projection_iterations,
                )?;
                poll(work.cancelled)?;
                return Ok(report(
                    trajectory,
                    speed_bound,
                    bounds,
                    maximum_wrench,
                    true,
                    errors,
                ));
            }
            next = if ratio < 1.0 / 32.0 {
                (2.0 * h).min(controls.integration.maximum_step_s)
            } else {
                h
            };
        }
        Err(TransientError::AttemptBudget.into())
    }
}

struct Work<'a, C> {
    law: &'a HertzFrictionPatch,
    quadrature: HertzFrictionQuadrature,
    config: HertzWrenchConfig,
    controls: HertzRateDrivenControls,
    metric: &'a Metric,
    radius: f64,
    speed_bound: f64,
    mass_bound: f64,
    moving_wrench_bound: f64,
    affine_bounds: [f64; 2],
    release: Option<(f64, V)>,
    evaluations: usize,
    cancelled: &'a mut C,
}
impl<C: FnMut() -> bool> Work<'_, C> {
    fn sample(&mut self, q: V) -> Result<Sample, TransientError> {
        if (self.cancelled)() {
            return Err(TransientError::Cancelled);
        }
        if self.evaluations >= self.controls.integration.maximum_evaluations {
            return Err(TransientError::EvaluationBudget);
        }
        if q == [0.0; 3] {
            return Err(TransientError::StepUnderflow);
        }
        let speed = det::hypot(q[0], q[1]) + q[2].abs();
        if !speed.is_finite()
            || speed > self.speed_bound
            || self.metric.radius(q)?.1 > 1.25 * self.mass_bound
        {
            return Err(TransientError::Invalid(
                "driven RK stage left its admitted rate domain",
            ));
        }
        if let Some((limit, separator)) = self.release {
            if speed > limit || dot(separator, q) <= 0.0 {
                return Err(TransientError::StepUnderflow);
            }
        }
        self.evaluations += 1;
        let response = self
            .law
            .evaluate([q[0], q[1]], q[2] / self.radius, self.quadrature)
            .map_err(HertzWrenchError::Law)?;
        let wrench = [
            response.force_n[0],
            response.force_n[1],
            response.spin_moment_nm / self.radius,
        ];
        if det::hypot(wrench[0], wrench[1]) > self.config.maximum_force_n
            || response.spin_moment_nm.abs() > self.config.maximum_spin_moment_nm
        {
            return Err(TransientError::EffortBudget);
        }
        if norm(wrench) > self.moving_wrench_bound
            || norm(wrench) > self.affine_bounds[0] + self.affine_bounds[1] * norm(q)
        {
            return Err(TransientError::Invalid(
                "original driven traction exceeded its admitted bound",
            ));
        }
        if response.at_rest || response.dissipated_power_w < 0.0 {
            return Err(TransientError::StepUnderflow);
        }
        if (self.cancelled)() {
            return Err(TransientError::Cancelled);
        }
        Ok(Sample {
            wrench,
            velocity: q,
            power: response.dissipated_power_w,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn finish(
    input: &HertzDrivenInput,
    original: V,
    endpoint: V,
    drive: V,
    demand: V,
    factor: Factor,
    radius: f64,
    impulse: V,
    travel: V,
    heat: f64,
    errors: [f64; 4],
    mode: HertzDrivenMode,
    config: HertzWrenchConfig,
    controls: HertzRateDrivenControls,
    accepted: usize,
    rejected: usize,
    evaluations: usize,
    iterations: usize,
) -> Result<HertzDrivenReport, Error> {
    let actual = vector(add(
        add(original, scale(drive, input.dt_s)),
        multiply(input.scaled_mobility_per_kg, impulse),
    ))?;
    let residual = vector(add(actual, scale(endpoint, -1.0)))?;
    let kinetic = checked(dot(
        solve(factor, add(actual, scale(original, -1.0)))?,
        scale(add(actual, original), 0.5),
    ))?;
    let drive_work = checked(dot(demand, travel))?;
    let balance = checked(kinetic + heat - drive_work)?;
    if heat < 0.0
        || norm(residual) > controls.integration.scaled_rate_tolerance_m_s
        || balance.abs() > config.energy_tolerance_j
        || errors
            .into_iter()
            .zip(error_limits(controls))
            .any(|(e, l)| e > l)
    {
        return Err(HertzDrivenError::UnresolvedWork {
            scaled_velocity_residual_m_s: residual,
            energy_residual_j: balance,
        }
        .into());
    }
    Ok(HertzDrivenReport {
        mode,
        velocity_after_m_s: [endpoint[0], endpoint[1]],
        spin_rate_after_rad_s: checked(endpoint[2] / radius)?,
        driving_scaled_wrench_n: demand,
        force_impulse_ns: [impulse[0], impulse[1]],
        spin_angular_impulse_nms: checked(radius * impulse[2])?,
        slip_increment_m: [travel[0], travel[1]],
        spin_increment_rad: checked(travel[2] / radius)?,
        dissipated_energy_j: heat,
        drive_work_j: drive_work,
        kinetic_energy_change_j: kinetic,
        energy_residual_j: balance,
        scaled_velocity_residual_m_s: residual,
        first_arrest_time_s: None,
        arrest_time_is_exact: false,
        accepted_steps: accepted,
        rejected_steps: rejected,
        evaluations,
        projection_iterations: iterations,
    })
}
fn report(
    trajectory: HertzDrivenReport,
    speed: f64,
    bounds: [f64; 2],
    wrench: f64,
    available: bool,
    errors: [f64; 4],
) -> HertzRateDrivenReport {
    let stopping_time_bracket_s = if trajectory.arrest_time_is_exact {
        trajectory.first_arrest_time_s.map(|t| [t, t])
    } else {
        None
    };
    HertzRateDrivenReport {
        trajectory,
        exact_invariant: false,
        stopping_time_bracket_s,
        local_speed_bound_m_s: speed,
        kinetic_coefficient_bounds: bounds,
        maximum_scaled_wrench_norm_n: wrench,
        adaptive_error_estimates_available: available,
        estimated_scaled_rate_error_m_s: errors[0],
        estimated_scaled_impulse_error_ns: errors[1],
        estimated_scaled_travel_error_m: errors[2],
        estimated_heat_error_j: errors[3],
    }
}
fn error_limits(c: HertzRateDrivenControls) -> [f64; 4] {
    [
        c.integration.scaled_rate_tolerance_m_s,
        c.integration.scaled_impulse_tolerance_ns,
        c.integration.scaled_travel_tolerance_m,
        c.integration.heat_tolerance_j,
    ]
}
fn event(t: f64, q: V) -> Error {
    Error::UnresolvedDrivenEvent {
        elapsed_s: t,
        scaled_rate_m_s: q,
    }
}
fn invalid(reason: &'static str) -> Error {
    Error::Transient(HertzTransientError::Invalid(reason))
}
fn checked(v: f64) -> Result<f64, Error> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("unrepresentable driven rate quantity"))
    }
}
fn vector(v: V) -> Result<V, Error> {
    for x in v {
        checked(x)?;
    }
    Ok(v)
}
fn solve(f: Factor, v: V) -> Result<V, Error> {
    vector(
        f.solve(v)
            .ok_or_else(|| invalid("unrepresentable original mobility solve"))?,
    )
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(TransientError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn norm(v: V) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn add(a: V, b: V) -> V {
    core::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: f64) -> V {
    a.map(|v| v * s)
}
fn multiply(m: M, v: V) -> V {
    m.map(|r| dot(r, v))
}
