//! Continuous inelastic capture of an inertial particle stream.
//!
//! Incoming velocity `u` is held in the world frame, independently of the
//! changing body velocity. With mass rate `q >= 0` and no external loads,
//! `M' = q`, `(M v)' = q u`, and `v(t)-u = M0*(v0-u)/M(t)`.
//! Hence capture produces ram drag even for stationary incoming particles.
//! This is the mass-gain law of MIT 16.07 L14, eq. (11):
//! <https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/98d19cb2b9983c14c6907052499e7a98_MIT16_07F09_Lec14.pdf>.
//!
//! The original native mass-integral kernel also gives the exact trajectory.
//! For `dm=q*h`, `x=dm/M1`, and its kernel `C`, the velocity increment is
//! `(u-v0)*x` and the displacement is `h*(v0+(u-v0)*x*C)`.
//! Inelastic heat is integrated independently as
//! `Q = integral(q*|u-v|^2/2 dt) = M0*x*|u-v0|^2/2`.
//! It is not manufactured from a residual or added to the bulk particle energy.
//!
//! Particles join at the instantaneous point-body position. The orbital influx
//! therefore uses that moving position, not either endpoint. Spin, finite-size
//! capture, outgoing particles and an inferred capture rate are outside this
//! operator. Held external force and gravity use the same variable-mass
//! trajectory through `advance_with_inertial_capture_forced`. Source inventory
//! and caloric admission remain with the composing material owner.

use super::ExhaustEnvironment;
use crate::Vec3;
use crate::rocket::{RocketError, RocketForcing, RocketState, mass_integrals, poll};

#[path = "capture/forced.rs"]
mod forced;

/// A physical incoming stream, held in the body's declared inertial frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleCaptureForcing {
    /// Nonnegative captured mass rate [kg/s], supplied by the source owner.
    pub mass_flow_kg_s: f64,
    /// Incoming mean velocity [m/s], before capture, held in world coordinates.
    pub velocity_world_m_s: Vec3,
    /// Incoming internal/random energy [W], excluding mean bulk kinetic energy.
    /// A formation datum may make this signed. It must be zero at zero mass rate.
    pub internal_energy_rate_w: f64,
}

/// Native endpoint and physical momentum, angular momentum and energy influx.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleCaptureStep {
    /// Complete new native state; the captured mass joins its existing inventory.
    pub state_after: RocketState,
    /// Complete held interval [s].
    pub duration_s: f64,
    /// Integrated incoming mass [kg], before native inventory rounding.
    pub captured_mass_kg: f64,
    /// Incoming physical momentum in world coordinates [kg m/s].
    pub incoming_momentum_world_kg_m_s: Vec3,
    /// Incoming orbital angular momentum about the world origin [kg m2/s].
    /// The capture location follows the continuous native trajectory.
    pub incoming_angular_momentum_world_kg_m2_s: Vec3,
    /// Applied force plus gravity impulse, integrated on retained mass [N s].
    pub external_momentum_world_kg_m_s: Vec3,
    /// Applied force plus gravity orbital impulse about the world origin [kg m2/s].
    pub external_angular_momentum_world_kg_m2_s: Vec3,
    /// Work of the held nongravitational force along the trajectory [J].
    pub external_work_j: f64,
    /// Work of uniform gravity on the continuously changing retained mass [J].
    pub gravity_work_j: f64,
    /// Incoming complete energy in the world frame [J].
    pub incoming_energy_world_j: f64,
    /// Incoming internal/random energy in its original source datum [J].
    pub incoming_internal_energy_j: f64,
    /// Nonnegative kinetic energy thermalized by inelastic capture [J].
    pub thermalization_energy_j: f64,
    /// Required material credit: incoming internal energy plus capture heat [J].
    /// The composing material owner applies this exactly once.
    pub retained_internal_energy_j: f64,
    /// Time mean of the body's velocity during the complete interval [m/s].
    pub mean_velocity_world_m_s: Vec3,
    /// Actual retained bulk kinetic energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Actual native inventory increase minus integrated incoming mass [kg].
    pub mass_residual_kg: f64,
    /// Actual momentum change minus incoming and external impulses [N s].
    pub momentum_residual_kg_m_s: Vec3,
    /// Actual orbital change minus incoming and external angular impulses [kg m2/s].
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Actual kinetic change plus heat minus incoming kinetic energy and load work.
    /// Internal/formation energies are excluded from this mechanical check [J].
    pub kinetic_energy_residual_j: f64,
}

impl RocketState {
    /// Capture particles whose incoming velocity is held in the inertial frame.
    ///
    /// This exact force-free operator distinguishes ambient capture from a held
    /// body-relative stream: its relative velocity changes as the body slows or
    /// accelerates. Position integrates that changing velocity, and inelastic
    /// heat is a separate, Galilean-invariant material credit. A zero mass rate
    /// gives free flight; zero duration is the admitted identity.
    ///
    /// All inputs remain borrowed. Cancellation or invalid/unrepresentable
    /// arithmetic publishes no state. Finite residuals expose actual rounding,
    /// including a velocity increment lost under a large inertial boost; they
    /// are not a conservation certificate or a repaired endpoint.
    pub fn advance_with_inertial_capture(
        &self,
        forcing: ParticleCaptureForcing,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ParticleCaptureStep, RocketError> {
        self.advance_with_inertial_capture_forced(
            forcing,
            ExhaustEnvironment::ZERO,
            duration_s,
            cancelled,
        )
    }

    /// Capture an inertial stream with held external force and uniform gravity.
    ///
    /// Integrates `M'=q`, `(M*v)'=q*u+F+M*g` on one native trajectory.
    /// Incoming particles, external angular impulse, work and inelastic heat
    /// follow that same trajectory. Gravity acts on retained mass; incoming
    /// velocity and source inventory remain the source owner's responsibility.
    ///
    /// All work is allocation-free with bounded scalar moment evaluations.
    /// A zero rate gives constant-mass forced motion; zero duration is an
    /// admitted identity. Cancellation and rejection publish no state.
    pub fn advance_with_inertial_capture_forced(
        &self,
        forcing: ParticleCaptureForcing,
        environment: ExhaustEnvironment,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ParticleCaptureStep, RocketError> {
        poll(&mut cancelled)?;
        if !environment.gravity_m_s2.is_finite() || !environment.external_force_n.is_finite() {
            return Err(RocketError::NonFiniteInput("inertial capture environment"));
        }
        if !forcing.mass_flow_kg_s.is_finite()
            || !forcing.velocity_world_m_s.is_finite()
            || !forcing.internal_energy_rate_w.is_finite()
            || !duration_s.is_finite()
        {
            return Err(RocketError::NonFiniteInput("inertial particle capture"));
        }
        if forcing.mass_flow_kg_s < 0.0
            || duration_s < 0.0
            || (forcing.mass_flow_kg_s == 0.0 && forcing.internal_energy_rate_w != 0.0)
        {
            return Err(RocketError::InvalidInput("inertial particle capture"));
        }
        // The existing native owner admits state, masses and inertial axes.
        self.advance(
            RocketForcing {
                mass_flow_kg_s: 0.0,
                effective_exhaust_velocity_m_s: Vec3::ZERO,
                gravity_m_s2: Vec3::ZERO,
                external_force_n: Vec3::ZERO,
            },
            0.0,
            &mut cancelled,
        )?;
        let captured = forcing.mass_flow_kg_s * duration_s;
        let internal = forcing.internal_energy_rate_w * duration_s;
        let mut next = *self;
        next.propellant_mass_kg = forcing
            .mass_flow_kg_s
            .mul_add(duration_s, self.propellant_mass_kg);
        let m0 = self.dry_mass_kg + self.propellant_mass_kg;
        let m1 = self.dry_mass_kg + next.propellant_mass_kg;
        if !captured.is_finite() || !m1.is_finite() || !internal.is_finite() {
            return Err(RocketError::Unrepresentable("captured particle inventory"));
        }
        if duration_s > 0.0
            && ((forcing.mass_flow_kg_s > 0.0
                && (captured == 0.0 || next.propellant_mass_kg == self.propellant_mass_kg))
                || (forcing.internal_energy_rate_w != 0.0 && internal == 0.0))
        {
            return Err(RocketError::Unrepresentable("captured particle increment"));
        }
        let x = captured / m1;
        if !x.is_finite() || x > 1.0 || (captured > 0.0 && x == 0.0) {
            return Err(RocketError::Unrepresentable("captured mass fraction"));
        }
        let (_, a, c) = mass_integrals(x, m1, m0);
        let v0 = self.velocity_m_s;
        let slip = if captured > 0.0 {
            forcing.velocity_world_m_s.sub(v0)
        } else {
            Vec3::ZERO
        };
        let mut mean_velocity = if x < 0.5 {
            let mean_increment = slip.scale(x * c);
            next.velocity_m_s = v0.add(slip.scale(x));
            next.position_m = self
                .position_m
                .add(v0.scale(duration_s))
                .add(mean_increment.scale(duration_s));
            v0.add(mean_increment)
        } else {
            // Keep the original body's residual velocity and displacement
            // when incoming mass dominates and x itself rounds to one.
            let k = m0 / m1;
            next.velocity_m_s = forcing.velocity_world_m_s.sub(slip.scale(k));
            let correction = slip.scale(k * a);
            next.position_m = self
                .position_m
                .add(forcing.velocity_world_m_s.scale(duration_s))
                .sub(correction.scale(duration_s));
            forcing.velocity_world_m_s.sub(correction)
        };
        poll(&mut cancelled)?;

        let incoming_momentum = forcing.velocity_world_m_s.scale(captured);
        let incoming_kinetic = if captured > 0.0 {
            (0.5 * captured) * forcing.velocity_world_m_s.norm_squared()
        } else {
            0.0
        };
        let slip_squared = slip.norm_squared();
        let mut thermalization = (0.5 * (m0 * x)) * slip_squared;
        // Integral of q*(position-position0) cross u. The large-gain form
        // retains the finite M0*log(M1/M0) contribution when C rounds to one.
        let angular_mass = if x < 0.125 {
            captured * (1.0 - c)
        } else {
            m0 * (a - 1.0)
        };
        let mut centered_angular = if captured > 0.0 {
            v0.cross(forcing.velocity_world_m_s)
                .scale(angular_mass * duration_s)
        } else {
            Vec3::ZERO
        };
        // Use the initial-velocity chart for small capture and the incoming
        // chart when capture dominates. In the latter, the original momentum
        // remains directly visible even when M0/M1 is below machine epsilon.
        // Keep the separately stored native masses in actual endpoint balances.
        let actual_mass = next.propellant_mass_kg - self.propellant_mass_kg;
        let mass_residual = actual_mass - captured;
        let boost = if x < 0.5 {
            v0
        } else {
            forcing.velocity_world_m_s
        };
        let loads = if environment == ExhaustEnvironment::ZERO || duration_s == 0.0 {
            forced::CaptureLoads::ZERO
        } else {
            let loads = forced::CaptureInterval {
                duration_s,
                captured_mass_kg: captured,
                initial_mass_kg: m0,
                final_mass_kg: m1,
                x,
                a,
                c,
            }
            .loads(self, forcing, environment, boost)?;
            next.velocity_m_s = next.velocity_m_s.add(loads.velocity_increment);
            next.position_m = next
                .position_m
                .add(loads.mean_velocity_increment.scale(duration_s));
            mean_velocity = mean_velocity.add(loads.mean_velocity_increment);
            centered_angular = centered_angular.add(loads.incoming_centered_angular_increment);
            thermalization = loads.thermalization_energy_j;
            loads
        };
        let incoming_angular = self
            .position_m
            .cross(incoming_momentum)
            .add(centered_angular);
        let before_relative = v0.sub(boost);
        let after_relative = next.velocity_m_s.sub(boost);
        let incoming_relative = if captured > 0.0 {
            forcing.velocity_world_m_s.sub(boost)
        } else {
            Vec3::ZERO
        };
        let momentum_before_relative = before_relative
            .scale(self.dry_mass_kg)
            .add(before_relative.scale(self.propellant_mass_kg));
        let momentum_after_relative = after_relative
            .scale(self.dry_mass_kg)
            .add(after_relative.scale(next.propellant_mass_kg));
        let momentum_change_relative = momentum_after_relative.sub(momentum_before_relative);
        let centered_momentum_residual = momentum_change_relative
            .sub(incoming_relative.scale(captured))
            .sub(loads.external_impulse);
        let momentum_residual = centered_momentum_residual.add(boost.scale(mass_residual));
        let kinetic_before_relative = (0.5 * self.dry_mass_kg) * before_relative.norm_squared()
            + (0.5 * self.propellant_mass_kg) * before_relative.norm_squared();
        let kinetic_after_relative = (0.5 * self.dry_mass_kg) * after_relative.norm_squared()
            + (0.5 * next.propellant_mass_kg) * after_relative.norm_squared();
        let kinetic_change_relative = kinetic_after_relative - kinetic_before_relative;
        let boost_squared = if actual_mass == 0.0 && mass_residual == 0.0 {
            0.0
        } else {
            boost.norm_squared()
        };
        let final_momentum = next
            .velocity_m_s
            .scale(self.dry_mass_kg)
            .add(next.velocity_m_s.scale(next.propellant_mass_kg));
        let angular_residual = next
            .position_m
            .sub(self.position_m)
            .cross(final_momentum)
            .sub(centered_angular)
            .sub(loads.external_centered_angular)
            .add(self.position_m.cross(momentum_residual));
        let result = ParticleCaptureStep {
            state_after: next,
            duration_s,
            captured_mass_kg: captured,
            incoming_momentum_world_kg_m_s: incoming_momentum,
            incoming_angular_momentum_world_kg_m2_s: incoming_angular,
            external_momentum_world_kg_m_s: loads.external_impulse,
            external_angular_momentum_world_kg_m2_s: self
                .position_m
                .cross(loads.external_impulse)
                .add(loads.external_centered_angular),
            external_work_j: loads.external_work_j,
            gravity_work_j: loads.gravity_work_j,
            incoming_energy_world_j: internal + incoming_kinetic,
            incoming_internal_energy_j: internal,
            thermalization_energy_j: thermalization,
            retained_internal_energy_j: internal + thermalization,
            mean_velocity_world_m_s: mean_velocity,
            kinetic_energy_change_j: boost.dot(momentum_change_relative)
                + kinetic_change_relative
                + (0.5 * actual_mass) * boost_squared,
            mass_residual_kg: mass_residual,
            momentum_residual_kg_m_s: momentum_residual,
            angular_momentum_residual_kg_m2_s: angular_residual,
            kinetic_energy_residual_j: boost.dot(centered_momentum_residual)
                + kinetic_change_relative
                + thermalization
                - (0.5 * captured) * incoming_relative.norm_squared()
                + (0.5 * mass_residual) * boost_squared
                - loads.work_relative_j,
        };
        if ![
            result.state_after.position_m,
            result.state_after.velocity_m_s,
            result.incoming_momentum_world_kg_m_s,
            result.incoming_angular_momentum_world_kg_m2_s,
            result.external_momentum_world_kg_m_s,
            result.external_angular_momentum_world_kg_m2_s,
            result.mean_velocity_world_m_s,
            result.momentum_residual_kg_m_s,
            result.angular_momentum_residual_kg_m2_s,
        ]
        .into_iter()
        .all(Vec3::is_finite)
            || ![
                result.incoming_energy_world_j,
                result.thermalization_energy_j,
                result.retained_internal_energy_j,
                result.kinetic_energy_change_j,
                result.mass_residual_kg,
                result.kinetic_energy_residual_j,
                result.external_work_j,
                result.gravity_work_j,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(RocketError::Unrepresentable("inertial capture moments"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}
