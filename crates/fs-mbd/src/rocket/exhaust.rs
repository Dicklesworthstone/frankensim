//! Physical exhaust in a common inertial frame, on the native rocket trajectory.
//!
//! Let `q` be positive expelled mass rate, `w=P_relative/q` the outward mean
//! velocity and `L(t)=ln(M0/M(t))`. With no external force the native trajectory
//! is `v(t)=v0-w*L(t)`. Its exported physical fluxes are
//! `q*v+P_relative` and `E_relative+v.dot(P_relative)+q*|v|²/2`.
//! Integrating these along the CHANGING velocity gives, with `dv=-w*L(end)`,
//! `P_out=dm*v0-M1*dv` and
//! `E_out=E_relative*burn_time+dm*|v0|²/2-M1*v0.dot(dv)-M1*|dv|²/2`.
//! Thus mass loss carries inertial momentum even when the relative thrust is zero.
//!
//! These are physical particle fluxes. A continuum pressure-thrust correction
//! must not be folded into `P_relative`; it is a separate surface interaction.
//! Relative energy includes all transported random/internal energy in exactly
//! one source datum and can be signed. The thermodynamic owner must debit it;
//! this mechanics owner does not create a second propellant caloric model.
//! Molecular beams with different velocities can be summed before this step:
//! their full second moment remains in the supplied relative energy.
//!
//! MIT 16.522, Lecture 2, eqs. (1)--(5), gives the variable-mass momentum
//! integral: https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/.
//! The energy expression above follows the Galilean particle-energy transform.
//! The step is exact for held relative fluxes, including burnout and coast;
//! changing gas states require caller-controlled time refinement. The forced
//! extension includes held gravity and external mechanical work. Rotation and
//! a shifting centre of mass are excluded.

use super::{RocketError, RocketForcing, RocketState, RocketStep, mass_integrals, poll};
use crate::Vec3;
use fs_math::det;

#[path = "exhaust/forced.rs"]
mod forced;
pub use forced::ExhaustEnvironment;

#[path = "exhaust/exchange.rs"]
mod exchange;
pub use exchange::{
    ParticleExchangeForcing, ParticleExchangeStep, ParticleStreamIntegral, ParticleStreamMoments,
};

/// Physical comoving exhaust moments, held during the powered interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExhaustForcing {
    /// Nonnegative actual outward mass rate [kg/s].
    pub mass_flow_kg_s: f64,
    /// Outward relative momentum rate [N], opposite the vehicle recoil.
    /// This is a particle moment, with no additional continuum pressure thrust.
    pub relative_momentum_rate_n: Vec3,
    /// Complete comoving exported energy [W], including its original datum.
    pub relative_energy_rate_w: f64,
}

/// Native trajectory and physical exported moments from the same interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExhaustStep {
    /// Original variable-mass trajectory; its thrust points opposite the exhaust.
    pub mechanics: RocketStep,
    /// Integrated comoving outward momentum [N s].
    pub relative_momentum_n_s: Vec3,
    /// Complete comoving exported energy [J], before a frame transformation.
    pub relative_energy_j: f64,
    /// Actual exported momentum in the initial state's inertial frame [kg m/s].
    pub exhaust_momentum_world_kg_m_s: Vec3,
    /// Complete exported energy in that same inertial frame [J].
    pub exhaust_energy_world_j: f64,
    /// Inertial kinetic correction to comoving energy [J], retained separately
    /// so a large arbitrary formation datum cannot erase recoil accounting.
    pub relative_kinetic_transport_energy_j: f64,
    /// Mean vehicle velocity during emission [m/s]; zero when no mass emits.
    /// A composing owner can use this to account for accepted-rate rounding.
    pub mean_velocity_world_m_s: Vec3,
    /// Mean squared vehicle speed during emission [m2/s2], including its variance.
    pub mean_speed_squared_m2_s2: f64,
    /// Vehicle kinetic change, evaluated with its actual endpoint velocity [J].
    pub kinetic_energy_change_j: f64,
    /// Impulse from the declared external force and gravity [N s], including
    /// coast. Gravity acts on the changing retained vehicle mass.
    pub external_impulse_n_s: Vec3,
    /// Work of those external loads on the moving vehicle [J], including coast.
    /// No gravitational potential or propellant caloric datum is introduced.
    pub external_work_j: f64,
    /// Integrated emitted mass minus actual native propellant decrement [kg].
    /// Its inertial momentum/energy effect is included in the residuals below.
    pub mass_residual_kg: f64,
    /// Actual retained momentum change plus exported momentum minus external
    /// impulse [kg m/s], in a stable chart about the initial velocity. This
    /// exposes lost small recoil.
    pub momentum_residual_kg_m_s: Vec3,
    /// Actual kinetic change plus the inertial export correction minus external
    /// work [J], without subtracting large common kinetic or formation baselines.
    pub kinetic_energy_residual_j: f64,
}

impl ExhaustStep {
    fn empty(mechanics: RocketStep) -> Self {
        Self {
            mechanics,
            relative_momentum_n_s: Vec3::ZERO,
            relative_energy_j: 0.0,
            exhaust_momentum_world_kg_m_s: Vec3::ZERO,
            exhaust_energy_world_j: 0.0,
            relative_kinetic_transport_energy_j: 0.0,
            mean_velocity_world_m_s: Vec3::ZERO,
            mean_speed_squared_m2_s2: 0.0,
            kinetic_energy_change_j: 0.0,
            external_impulse_n_s: Vec3::ZERO,
            external_work_j: 0.0,
            mass_residual_kg: 0.0,
            momentum_residual_kg_m_s: Vec3::ZERO,
            kinetic_energy_residual_j: 0.0,
        }
    }
}

// Average L(t)^2 over a linearly decreasing mass interval. The positive series
// is 2*sum_{n>=2} H_(n-1)*x^n/[n(n+1)]; it avoids cancellation for tiny burns.
fn mean_log_squared(x: f64, log_ratio: f64, mean_log: f64, mass_ratio: f64) -> f64 {
    if x < 0.125 {
        let mut coefficients = [0.0_f64; 19];
        let mut harmonic = 1.0;
        for (n, coefficient) in coefficients.iter_mut().enumerate().skip(2) {
            let k = n as f64;
            *coefficient = 2.0 * harmonic / (k * (k + 1.0));
            harmonic += k.recip();
        }
        let mut value = 0.0_f64;
        for coefficient in coefficients[2..].iter().rev() {
            value = value.mul_add(x, *coefficient);
        }
        (value * x) * x
    } else {
        (-(mass_ratio * log_ratio)).mul_add(log_ratio, 2.0 * mean_log)
    }
}

impl RocketState {
    /// Advance a free translating vehicle using physical exhaust moments.
    ///
    /// The original native rocket owner supplies mass depletion, position,
    /// velocity and exact burnout. Fluxes use the changing vehicle velocity,
    /// not a frozen initial or endpoint inertial exhaust speed. A zero mass
    /// rate requires zero particle momentum and energy; an exhausted vehicle
    /// emits nothing and coasts. Cancellation returns no changed state.
    ///
    /// Finite physical residuals are returned for the composing owner's absolute
    /// momentum/energy budgets. In particular a boost so large that an actual
    /// velocity increment rounds away is exposed, not certified as conservative.
    pub fn advance_with_exhaust(
        &self,
        forcing: ExhaustForcing,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExhaustStep, RocketError> {
        poll(&mut cancelled)?;
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
                gravity_m_s2: Vec3::ZERO,
                external_force_n: Vec3::ZERO,
            },
            duration_s,
            &mut cancelled,
        )?;
        let mut step = ExhaustStep::empty(mechanics);
        let dm = mechanics.burned_propellant_mass_kg;
        if dm == 0.0 {
            poll(&mut cancelled)?;
            return Ok(step);
        }
        let m0 = self.dry_mass_kg + self.propellant_mass_kg;
        let m1 = self.dry_mass_kg + mechanics.state_after.propellant_mass_kg;
        let x = dm / m0;
        let (log_ratio, _, position_integral) = mass_integrals(x, m0, m1);
        let ideal_dv = w.scale(-log_ratio);
        let actual_dv = mechanics.state_after.velocity_m_s.sub(self.velocity_m_s);
        let v0 = self.velocity_m_s;
        let relative_export = ideal_dv.scale(-m1);
        let mean_log = x * position_integral;
        let mean_log_squared = mean_log_squared(x, log_ratio, mean_log, m1 / dm);
        let variance_log = mean_log_squared - mean_log * mean_log;
        if !variance_log.is_finite() || variance_log < 0.0 {
            return Err(RocketError::Unrepresentable("exhaust velocity variance"));
        }
        step.mean_velocity_world_m_s = v0.sub(w.scale(mean_log));
        let velocity_standard_deviation = w.scale(det::sqrt(variance_log));
        step.mean_speed_squared_m2_s2 = step.mean_velocity_world_m_s.norm_squared()
            + velocity_standard_deviation.norm_squared();
        step.relative_momentum_n_s = mechanics.thrust_impulse_n_s.scale(-1.0);
        step.relative_energy_j = forcing.relative_energy_rate_w * mechanics.powered_duration_s;
        step.exhaust_momentum_world_kg_m_s = v0.scale(dm).add(relative_export);
        let initial_advected_kinetic = (0.5 * dm) * v0.norm_squared();
        let ideal_recoil_kinetic = (0.5 * m1) * ideal_dv.norm_squared();
        step.relative_kinetic_transport_energy_j =
            initial_advected_kinetic + v0.dot(relative_export) - ideal_recoil_kinetic;
        step.exhaust_energy_world_j =
            step.relative_energy_j + step.relative_kinetic_transport_energy_j;
        let actual_mass_loss = self.propellant_mass_kg - mechanics.state_after.propellant_mass_kg;
        step.mass_residual_kg = dm - actual_mass_loss;
        step.kinetic_energy_change_j = v0.dot(actual_dv.scale(m1))
            + (0.5 * m1) * actual_dv.norm_squared()
            - (0.5 * actual_mass_loss) * v0.norm_squared();
        // Cancel the common velocity analytically, but use the ACTUAL accepted
        // velocity increment. Endpoint rounding therefore remains observable.
        let dv_defect = actual_dv.sub(ideal_dv);
        let recoil_momentum_defect = dv_defect.scale(m1);
        step.momentum_residual_kg_m_s = recoil_momentum_defect.add(v0.scale(step.mass_residual_kg));
        step.kinetic_energy_residual_j = v0.dot(recoil_momentum_defect)
            + (0.5 * m1) * dv_defect.dot(actual_dv.add(ideal_dv))
            + (0.5 * step.mass_residual_kg) * v0.norm_squared();
        if !step.mean_velocity_world_m_s.is_finite()
            || !step.relative_momentum_n_s.is_finite()
            || !step.exhaust_momentum_world_kg_m_s.is_finite()
            || !step.momentum_residual_kg_m_s.is_finite()
            || ![
                step.relative_energy_j,
                step.exhaust_energy_world_j,
                step.relative_kinetic_transport_energy_j,
                step.mean_speed_squared_m2_s2,
                step.kinetic_energy_change_j,
                step.mass_residual_kg,
                step.kinetic_energy_residual_j,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(RocketError::Unrepresentable("inertial exhaust moments"));
        }
        poll(&mut cancelled)?;
        Ok(step)
    }
}

#[cfg(test)]
#[path = "exhaust/tests.rs"]
mod tests;
