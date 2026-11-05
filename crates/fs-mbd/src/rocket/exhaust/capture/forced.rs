//! Exact held-load moments on the native inelastic-capture trajectory.
//!
//! With s=t/h, x=dm/M1, k=M0/M1 and psi=s/(k+x*s), velocity is
//! v=v0+B*psi+G*s/2, B=(u-v0)*x+F*h/M1+k*G/2, G=g*h.
//! C=<psi>, K=<(1-s)*psi>=(C-1/2)/x, and
//! J=<(k+x*s)*integral_0^s psi>=(C+k*K)/2-1/6.
//! These moments integrate the moving capture point and the mass-weighted
//! gravity torque. Heat uses the squared mean slip plus an orthogonal,
//! positive velocity variance; it is never inferred from energy imbalance.

use super::{ExhaustEnvironment, ParticleCaptureForcing};
use crate::Vec3;
use crate::rocket::{RocketError, RocketState};
use fs_math::det;

pub(super) struct CaptureInterval {
    pub duration_s: f64,
    pub captured_mass_kg: f64,
    pub initial_mass_kg: f64,
    pub final_mass_kg: f64,
    pub x: f64,
    pub a: f64,
    pub c: f64,
}

pub(super) struct CaptureLoads {
    pub velocity_increment: Vec3,
    pub mean_velocity_increment: Vec3,
    pub incoming_centered_angular_increment: Vec3,
    pub external_impulse: Vec3,
    pub external_centered_angular: Vec3,
    pub external_work_j: f64,
    pub gravity_work_j: f64,
    pub work_relative_j: f64,
    pub thermalization_energy_j: f64,
}

impl CaptureLoads {
    pub const ZERO: Self = Self {
        velocity_increment: Vec3::ZERO,
        mean_velocity_increment: Vec3::ZERO,
        incoming_centered_angular_increment: Vec3::ZERO,
        external_impulse: Vec3::ZERO,
        external_centered_angular: Vec3::ZERO,
        external_work_j: 0.0,
        gravity_work_j: 0.0,
        work_relative_j: 0.0,
        thermalization_energy_j: 0.0,
    };
}

impl CaptureInterval {
    pub fn loads(
        &self,
        initial: &RocketState,
        forcing: ParticleCaptureForcing,
        environment: ExhaustEnvironment,
        boost: Vec3,
    ) -> Result<CaptureLoads, RocketError> {
        let Self {
            duration_s: h,
            captured_mass_kg: dm,
            initial_mass_kg: m0,
            final_mass_kg: m1,
            x,
            a,
            c,
        } = *self;
        // Keep k independently when x rounds to one: the original body still
        // carries momentum and thermalizes a finite amount of kinetic energy.
        let k = m0 / m1;
        let force_impulse = environment.external_force_n.scale(h);
        let f = Vec3::new(
            force_impulse.x / m1,
            force_impulse.y / m1,
            force_impulse.z / m1,
        );
        let g = environment.gravity_m_s2.scale(h);
        let gravity_impulse = g.scale(m0 + 0.5 * dm);
        let slip = if dm > 0.0 {
            forcing.velocity_world_m_s.sub(initial.velocity_m_s)
        } else {
            Vec3::ZERO
        };
        let (position_kernel, position_tail, covariance_tail) = if x < 0.125 {
            let (mut position, mut tail, mut covariance) = (0.0_f64, 0.0_f64, 0.0_f64);
            for n in (0..=16).rev() {
                let n = f64::from(n);
                position = position.mul_add(x, 1.0 / ((n + 2.0) * (n + 3.0)));
                tail = tail.mul_add(x, 1.0 / ((n + 3.0) * (n + 4.0)));
                covariance =
                    covariance.mul_add(x, -(n + 1.0) / (2.0 * (n + 3.0) * (n + 4.0) * (n + 5.0)));
            }
            (position, tail, covariance)
        } else {
            let position = (c - 0.5) / x;
            (position, (position - 1.0 / 6.0) / x, 0.0)
        };
        // Center the load corrections on f+g. This retains the O(x) response
        // when a constant applied force nearly balances gravity at small q.
        let combined = f.add(g);
        let velocity_increment = combined.sub(g.scale(0.5 * x));
        let mean_velocity_increment = combined
            .scale(c)
            .sub(g.scale(0.5 * x * (position_kernel + c)));
        let mean_position_load = combined
            .scale(position_kernel)
            .sub(g.scale(0.5 * x * (position_tail + position_kernel)))
            .scale(h);
        let mean_position_free = if x < 0.5 {
            initial
                .velocity_m_s
                .scale(0.5)
                .add(slip.scale(x * position_kernel))
        } else {
            forcing
                .velocity_world_m_s
                .scale(0.5)
                .sub(slip.scale(k * (a - c)))
        }
        .scale(h);
        let incoming_centered_angular_increment = if dm > 0.0 {
            mean_position_load
                .cross(forcing.velocity_world_m_s)
                .scale(dm)
        } else {
            Vec3::ZERO
        };
        let force_centered_angular = mean_position_free
            .add(mean_position_load)
            .cross(force_impulse);
        let j = 0.5 * (c + k * position_kernel) - 1.0 / 6.0;
        let gravity_position_moment = if x < 0.5 {
            initial
                .velocity_m_s
                .scale(0.5 * m0 + dm / 3.0)
                .add(slip.scale(dm * j))
        } else {
            // In the incoming chart this retains the M0*log(M1/M0) term,
            // instead of subtracting two O(dm) position moments.
            let retained = 0.5 * (1.0 + k) * a - 0.25 - 0.5 * k * c;
            forcing
                .velocity_world_m_s
                .scale(0.5 * m0 + dm / 3.0)
                .sub(slip.scale(m0 * retained))
        }
        .add(force_impulse.scale(j));
        // The remaining position contribution is parallel to gravity.
        let gravity_centered_angular = gravity_position_moment.scale(h).cross(g);

        let relative_initial = initial.velocity_m_s.sub(boost);
        let relative_incoming = if dm > 0.0 {
            forcing.velocity_world_m_s.sub(boost)
        } else {
            Vec3::ZERO
        };
        let mean_relative_free = if x < 0.5 {
            slip.scale(x * c)
        } else {
            slip.scale(-k * a)
        };
        let mean_relative_velocity = mean_relative_free.add(mean_velocity_increment);
        let force_work_relative = force_impulse.dot(mean_relative_velocity);
        // M(t)*v(t) is a quadratic polynomial even though v(t) is rational.
        // Integrate that polynomial directly for gravity work.
        let initial_load_impulse = force_impulse.add(g.scale(m0));
        let gravity_work_relative = g.dot(
            relative_initial
                .scale(m0)
                .add(relative_incoming.scale(0.5 * dm))
                .add(initial_load_impulse.scale(0.5))
                .add(g.scale(dm / 6.0)),
        );
        // Integrate the total load directly as well. Adding large opposite
        // force/gravity works would erase the mechanical work near balance.
        let time_velocity_kernel = c - position_kernel;
        let time_relative_free = if x < 0.5 {
            slip.scale(x * time_velocity_kernel)
        } else {
            slip.scale(-k * c)
        };
        let time_relative_velocity = time_relative_free
            .add(combined.scale(time_velocity_kernel))
            .sub(g.scale(0.5 * x * (time_velocity_kernel + position_kernel - position_tail)));
        let work_relative_j =
            initial_load_impulse.dot(mean_relative_velocity) + dm * g.dot(time_relative_velocity);

        let thermalization_energy_j = if dm == 0.0 {
            0.0
        } else {
            let b = slip.scale(x).add(f).add(g.scale(0.5 * k));
            let mean_relative = slip.scale(-k * a).add(mean_velocity_increment);
            let sqrt_twelve = det::sqrt(12.0);
            let sqrt_half_mass = det::sqrt(dm) * core::f64::consts::FRAC_1_SQRT_2;
            let (linear, nonlinear) = if x < 0.125 {
                // R=Var(psi)-12*Cov(psi,s)^2. The first two coefficients
                // vanish identically; evaluate R/x^2 without cancellation.
                let coefficients = [
                    1.0 / 180.0,
                    1.0 / 180.0,
                    13.0 / 3150.0,
                    17.0 / 6300.0,
                    1.0 / 630.0,
                    1.0 / 1260.0,
                    181.0 / 727650.0,
                    -169.0 / 1455300.0,
                    -2441.0 / 6879600.0,
                    -7663.0 / 15135120.0,
                    -905.0 / 1513512.0,
                    -7013.0 / 10810800.0,
                    -370327.0 / 551350800.0,
                    -1304309.0 / 1929727800.0,
                    -5385503.0 / 8066262204.0,
                    -5837081.0 / 8962513560.0,
                    -56443523.0 / 89625135600.0,
                ];
                let mut remainder = 0.0_f64;
                for coefficient in coefficients.into_iter().rev() {
                    remainder = remainder.mul_add(x, coefficient);
                }
                let linear = combined
                    .add(slip.sub(g.scale(0.5)).scale(x))
                    .scale(sqrt_twelve.recip())
                    .add(b.scale(sqrt_twelve * x * x * covariance_tail));
                (
                    linear.scale(sqrt_half_mass),
                    b.scale(sqrt_half_mass * x * det::sqrt(remainder)),
                )
            } else {
                let covariance_over_k = (1.0 + k).mul_add(a, -2.0) / (2.0 * x * x);
                // Factor out k before evaluating R. Multiplying it back as
                // dm*k=M0*x retains heat even if M0/M1 underflows to zero.
                let remainder_over_k = (-12.0 * k * covariance_over_k)
                    .mul_add(covariance_over_k, (-k * a).mul_add(a, 1.0) / (x * x));
                if !remainder_over_k.is_finite() || remainder_over_k < 0.0 {
                    return Err(RocketError::Unrepresentable(
                        "forced inertial capture variance",
                    ));
                }
                (
                    b.scale(sqrt_twelve * k * covariance_over_k)
                        .add(g.scale(0.5 / sqrt_twelve))
                        .scale(sqrt_half_mass),
                    b.scale(det::sqrt(m0) * det::sqrt(0.5 * x * remainder_over_k)),
                )
            };
            mean_relative.scale(sqrt_half_mass).norm_squared()
                + linear.norm_squared()
                + nonlinear.norm_squared()
        };
        Ok(CaptureLoads {
            velocity_increment,
            mean_velocity_increment,
            incoming_centered_angular_increment,
            external_impulse: force_impulse.add(gravity_impulse),
            external_centered_angular: force_centered_angular.add(gravity_centered_angular),
            external_work_j: force_work_relative + boost.dot(force_impulse),
            gravity_work_j: gravity_work_relative + boost.dot(gravity_impulse),
            work_relative_j,
            thermalization_energy_j,
        })
    }
}
