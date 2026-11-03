//! Held external loads on the physical-exhaust trajectory.
//!
//! During emission, `v(s)=v0+d*phi(s)+g*b*s`, where `s=t/b`,
//! `phi=-ln(1-x*s)/x`, `x=dm/M0`, and `d=F*b/M0-w*x`.
//! This normalization never forms the potentially overflowing `F/q`.
//! Three scalar moments integrate particle momentum, energy and gravity work;
//! a positive variance decomposition preserves opposing force/gravity limits.

use super::{ExhaustForcing, ExhaustStep};
use crate::Vec3;
use crate::rocket::{RocketError, RocketForcing, RocketState, mass_integrals, poll};
use fs_math::det;

/// External loads held in the same inertial frame as the vehicle and exhaust.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExhaustEnvironment {
    /// Uniform gravitational acceleration [m/s2]. It acts on retained mass.
    pub gravity_m_s2: Vec3,
    /// Other applied force [N], excluding exhaust recoil and gravity.
    /// Remains active through any coast after burnout.
    pub external_force_n: Vec3,
}

impl ExhaustEnvironment {
    /// Free translation, using the original optimized exhaust calculation.
    pub const ZERO: Self = Self {
        gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::ZERO,
    };
}

// Return mean velocity increment, velocity variance, and <s*phi>.
// Center the variance on s before squaring. For x -> 0, phi -> s and its
// covariance becomes rank one; subtracting raw second moments would erase
// the remaining O(x^2) variance when force acceleration cancels gravity.
pub(super) fn velocity_moments(
    x: f64,
    a: f64,
    c: f64,
    final_mass_ratio: f64,
    d: Vec3,
    gravity_increment: Vec3,
) -> Result<(Vec3, f64, f64), RocketError> {
    let sqrt_twelve = det::sqrt(12.0);
    let (mean, linear, nonlinear, time_log) = if x < 0.125 {
        let (mut mean_tail, mut covariance_tail, mut time_log) = (0.0_f64, 0.0_f64, 0.0_f64);
        for n in (0..=16).rev() {
            let k = f64::from(n);
            mean_tail = mean_tail.mul_add(x, 1.0 / ((k + 2.0) * (k + 3.0)));
            covariance_tail = covariance_tail.mul_add(x, 0.5 / ((k + 3.0) * (k + 4.0)));
            time_log = time_log.mul_add(x, 1.0 / ((k + 1.0) * (k + 3.0)));
        }
        // K=Var(phi)-12*Cov(phi,s)^2. All coefficients of K/x^2
        // retained here are positive; no floor or artificial variance is used.
        let coefficients = [
            1.0 / 720.0,
            1.0 / 360.0,
            2.0 / 525.0,
            113.0 / 25200.0,
            41.0 / 8400.0,
            8.0 / 1575.0,
            2993.0 / 582120.0,
            29761.0 / 5821200.0,
            56323.0 / 11211200.0,
            741203.0 / 151351200.0,
            1901.0 / 400400.0,
            39659.0 / 8648640.0,
            4871359.0 / 1102701600.0,
            630677.0 / 148440600.0,
            121943009.0 / 29875045200.0,
            702472649.0 / 179250271200.0,
            53940071.0 / 14340021696.0,
        ];
        let mut remainder = 0.0_f64;
        for coefficient in coefficients.into_iter().rev() {
            remainder = remainder.mul_add(x, coefficient);
        }
        let combined = d.add(gravity_increment);
        (
            combined.scale(0.5).add(d.scale(x * mean_tail)),
            combined
                .scale(sqrt_twelve.recip())
                .add(d.scale(sqrt_twelve * x * covariance_tail)),
            d.scale(x * det::sqrt(remainder)),
            time_log,
        )
    } else {
        let time_log = ((0.5 * (1.0 + x)).mul_add(c, -0.25)) / x;
        let second = (-(final_mass_ratio * a)).mul_add(a, 2.0 * c) / x;
        let covariance = time_log - 0.5 * c;
        let remainder = (-12.0 * covariance).mul_add(covariance, (-c).mul_add(c, second));
        if !remainder.is_finite() || remainder < 0.0 {
            return Err(RocketError::Unrepresentable(
                "forced exhaust velocity variance",
            ));
        }
        (
            d.scale(c).add(gravity_increment.scale(0.5)),
            d.scale(sqrt_twelve * covariance)
                .add(gravity_increment.scale(sqrt_twelve.recip())),
            d.scale(det::sqrt(remainder)),
            time_log,
        )
    };
    let variance = linear.norm_squared() + nonlinear.norm_squared();
    if !mean.is_finite() || !variance.is_finite() || !time_log.is_finite() {
        return Err(RocketError::Unrepresentable(
            "forced exhaust velocity moments",
        ));
    }
    Ok((mean, variance, time_log))
}

impl RocketState {
    /// Integrate physical exhaust with uniform gravity and a held external force.
    ///
    /// Uses the original variable-mass translation and exact burnout split.
    /// Exported particle momentum and energy follow the changing vehicle speed.
    /// Gravity impulse integrates retained mass; its work integrates the
    /// correlated mass and velocity, rather than multiplying their averages.
    /// The other force contributes `F dot displacement`, including coast.
    ///
    /// Residuals use the actual returned endpoint and retain rounding defects.
    /// The constant-load trajectory is analytic; varying thermochemistry or
    /// loads require caller-controlled time advancement. All work is O(1),
    /// allocation-free and cancellation-safe. No rotational claim is made.
    pub fn advance_with_exhaust_forced(
        &self,
        forcing: ExhaustForcing,
        environment: ExhaustEnvironment,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExhaustStep, RocketError> {
        poll(&mut cancelled)?;
        if environment == ExhaustEnvironment::ZERO {
            return self.advance_with_exhaust(forcing, duration_s, cancelled);
        }
        if !environment.gravity_m_s2.is_finite() || !environment.external_force_n.is_finite() {
            return Err(RocketError::NonFiniteInput("exhaust environment"));
        }
        if !forcing.mass_flow_kg_s.is_finite()
            || !forcing.relative_momentum_rate_n.is_finite()
            || !forcing.relative_energy_rate_w.is_finite()
        {
            return Err(RocketError::NonFiniteInput("physical exhaust moments"));
        }
        if forcing.mass_flow_kg_s < 0.0
            || (forcing.mass_flow_kg_s == 0.0
                && (forcing.relative_momentum_rate_n != Vec3::ZERO
                    || forcing.relative_energy_rate_w != 0.0))
        {
            return Err(RocketError::InvalidInput(
                "physical exhaust mass and moments",
            ));
        }
        let q = forcing.mass_flow_kg_s;
        let w = if q == 0.0 {
            Vec3::ZERO
        } else {
            Vec3::new(
                forcing.relative_momentum_rate_n.x / q,
                forcing.relative_momentum_rate_n.y / q,
                forcing.relative_momentum_rate_n.z / q,
            )
        };
        if !w.is_finite() {
            return Err(RocketError::Unrepresentable("relative exhaust velocity"));
        }
        let mechanics = self.advance(
            RocketForcing {
                mass_flow_kg_s: q,
                effective_exhaust_velocity_m_s: w.scale(-1.0),
                gravity_m_s2: environment.gravity_m_s2,
                external_force_n: environment.external_force_n,
            },
            duration_s,
            &mut cancelled,
        )?;
        if duration_s == 0.0 {
            return Ok(ExhaustStep::empty(mechanics));
        }
        let dm = mechanics.burned_propellant_mass_kg;
        let b = mechanics.powered_duration_s;
        let coast = mechanics.coast_duration_s;
        let m0 = self.dry_mass_kg + self.propellant_mass_kg;
        let m1 = self.dry_mass_kg + mechanics.state_after.propellant_mass_kg;
        let g = environment.gravity_m_s2;
        let force = environment.external_force_n;
        let v0 = self.velocity_m_s;
        let x = dm / m0;
        let (log_ratio, a, c) = if dm > 0.0 {
            mass_integrals(x, m0, m1)
        } else {
            (0.0, 1.0, 0.5)
        };
        let gravity_increment = g.scale(b);
        let d = force.scale(b / m0).sub(w.scale(x));
        let (mean_increment, velocity_variance, time_log) =
            velocity_moments(x, a, c, m1 / m0, d, gravity_increment)?;
        let mean_squared_increment = mean_increment.norm_squared() + velocity_variance;
        let burn_increment = w
            .scale(-log_ratio)
            .add(force.scale((b / m0) * a))
            .add(gravity_increment);
        let coast_increment = force.scale(coast / m1).add(g.scale(coast));
        let coast_displacement = burn_increment
            .scale(coast)
            .add(coast_increment.scale(0.5 * coast));
        // Integrate displacement analytically. Subtracting huge absolute
        // endpoint positions would erase small but physically real work.
        let displacement_relative = mean_increment.scale(b).add(coast_displacement);
        let gravity_impulse = gravity_increment
            .scale(m0 * (1.0 - 0.5 * x))
            .add(g.scale(coast).scale(m1));
        let external_impulse = force.scale(duration_s).add(gravity_impulse);
        let weighted_velocity = d
            .scale((-x).mul_add(time_log, c))
            .add(gravity_increment.scale(0.5 - x / 3.0));
        let gravity_work_relative = gravity_increment.dot(weighted_velocity.scale(m0))
            + g.dot(coast_displacement.scale(m1));
        let external_work_relative = force.dot(displacement_relative) + gravity_work_relative;
        let external_work = v0.dot(external_impulse) + external_work_relative;
        let relative_momentum = mechanics.thrust_impulse_n_s.scale(-1.0);
        let relative_energy = forcing.relative_energy_rate_w * b;
        let export_momentum_relative = mean_increment.scale(dm).add(relative_momentum);
        let export_kinetic_relative =
            relative_momentum.dot(mean_increment) + (0.5 * dm) * mean_squared_increment;
        let initial_speed_squared = if dm == 0.0 { 0.0 } else { v0.norm_squared() };
        let export_kinetic = (0.5 * dm) * initial_speed_squared
            + v0.dot(export_momentum_relative)
            + export_kinetic_relative;
        let actual_increment = mechanics.state_after.velocity_m_s.sub(v0);
        let actual_mass_loss = self.propellant_mass_kg - mechanics.state_after.propellant_mass_kg;
        let mass_residual = dm - actual_mass_loss;
        let momentum_residual_relative = actual_increment
            .scale(m1)
            .add(export_momentum_relative)
            .sub(external_impulse);
        let kinetic_change = v0.dot(actual_increment.scale(m1))
            + (0.5 * m1) * actual_increment.norm_squared()
            - (0.5 * actual_mass_loss) * initial_speed_squared;
        let energy_residual = v0.dot(momentum_residual_relative)
            + (0.5 * m1) * actual_increment.norm_squared()
            + export_kinetic_relative
            - external_work_relative
            + (0.5 * mass_residual) * initial_speed_squared;
        let mean_velocity = if dm == 0.0 {
            Vec3::ZERO
        } else {
            v0.add(mean_increment)
        };
        // A change of inertial frame changes only the mean, not the variance.
        let mean_speed_squared = if dm == 0.0 {
            0.0
        } else {
            mean_velocity.norm_squared() + velocity_variance
        };
        let step = ExhaustStep {
            mechanics,
            relative_momentum_n_s: relative_momentum,
            relative_energy_j: relative_energy,
            exhaust_momentum_world_kg_m_s: v0.scale(dm).add(export_momentum_relative),
            exhaust_energy_world_j: relative_energy + export_kinetic,
            relative_kinetic_transport_energy_j: export_kinetic,
            mean_velocity_world_m_s: mean_velocity,
            mean_speed_squared_m2_s2: mean_speed_squared,
            kinetic_energy_change_j: kinetic_change,
            external_impulse_n_s: external_impulse,
            external_work_j: external_work,
            mass_residual_kg: mass_residual,
            momentum_residual_kg_m_s: momentum_residual_relative.add(v0.scale(mass_residual)),
            kinetic_energy_residual_j: energy_residual,
        };
        if !step.relative_momentum_n_s.is_finite()
            || !step.exhaust_momentum_world_kg_m_s.is_finite()
            || !step.mean_velocity_world_m_s.is_finite()
            || !step.external_impulse_n_s.is_finite()
            || !step.momentum_residual_kg_m_s.is_finite()
            || ![
                step.relative_energy_j,
                step.exhaust_energy_world_j,
                step.relative_kinetic_transport_energy_j,
                step.mean_speed_squared_m2_s2,
                step.kinetic_energy_change_j,
                step.external_work_j,
                step.mass_residual_kg,
                step.kinetic_energy_residual_j,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(RocketError::Unrepresentable(
                "forced inertial exhaust moments",
            ));
        }
        poll(&mut cancelled)?;
        Ok(step)
    }
}

#[cfg(test)]
#[path = "forced/tests.rs"]
mod tests;
