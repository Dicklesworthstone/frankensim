//! Curved Hertz sliding/spinning with the original pressure-resolved traction.
//!
//! In scaled coordinates `q=[u_x,u_y,a*omega]`, the complete original mobility
//! drives `q'=W*w(q)`. The shared circular-traction integrator evaluates this
//! Hertz law at every positive classical RK4 stage. Two half-steps supply the
//! physical candidate; their difference from one full step estimates error.
//! No Richardson-extrapolated physical state or negative heat weight is used.
//! The independent power quadrature never absorbs the numerical work residual.
//!
//! Let `C=mu*N`, `s=a*omega`, and `k=3*pi/16`, the Hertz pressure's mean radius
//! divided by the outer radius. The dissipation integral is even and convex in
//! either u or s; angular symmetry and Jensen's inequality give
//! `H(q)>=C*max(|u|,k*|s|)>=C*k/sqrt(1+k*k)*|q|`. The triangle inequality gives
//! `H(q)<=C*(|u|+k*|s|)<=C*sqrt(1+k*k)*|q|`. These pressure-specific bounds,
//! and the original positive-definite mobility, bound finite extinction.
//!
//! A terminal tail is used only when its whole stopping interval is before the
//! requested horizon and its time/travel uncertainty fits the supplied limits.
//! Its impulse is the actual eventual-rest impulse `-W^-1*q`; its heat is the
//! remaining relative kinetic energy. Uncomputed travel stays an explicit bound.
//! Prefix errors are numerical estimates, not validated trajectory enclosures.
//! Geometry, load, coefficients and mobility remain fixed, without external drive.

use super::relaxation::{CircularRelaxationInput, HertzRelaxationError, HertzRelaxationInput};
use super::{HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, SPIN_RADIUS};
use crate::adhesive::spinning::transient::{
    CircularTransientError, CircularTransientLaw, CircularTransientQuadrature,
    CircularTransientTraction, TransientError, advance_circular_transient,
};
use fs_math::det;

/// Common physical-time accuracy and bounded numerical work controls.
pub use crate::adhesive::spinning::transient::AdhesiveTransientControls as HertzTransientControls;
/// Integrated actual Hertz traction, resolved travel and bounded stopping tail.
pub use crate::adhesive::spinning::transient::AdhesiveTransientReport as HertzTransientReport;
/// Hertz admission and the common numerical, physical and cancellation refusals.
pub type HertzTransientError = CircularTransientError<HertzRelaxationError>;

impl HertzWrenchContact {
    /// Resolve passive curved slip/spin under the complete native mobility.
    /// Exact invariant directions retain the one-evaluation analytic fast path.
    /// Every other accepted stage uses the original Hertz pressure integral;
    /// static friction never captures a moving state before physical extinction.
    pub fn advance_transient(
        &self,
        input: &HertzRelaxationInput,
        config: HertzWrenchConfig,
        controls: HertzTransientControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzTransientReport, HertzTransientError> {
        if cancelled() {
            return Err(HertzTransientError::Cancelled);
        }
        controls.validate().map_err(TransientError::with_exact)?;
        match self.relax_velocity(input, config, &mut cancelled) {
            Ok(exact) => {
                let residual = det::hypot(
                    det::hypot(
                        exact.scaled_velocity_residual_m_s[0],
                        exact.scaled_velocity_residual_m_s[1],
                    ),
                    exact.scaled_velocity_residual_m_s[2],
                );
                if !residual.is_finite() || residual > controls.scaled_rate_tolerance_m_s {
                    return Err(HertzTransientError::UnresolvedWork);
                }
                if cancelled() {
                    return Err(HertzTransientError::Cancelled);
                }
                return Ok(HertzTransientReport {
                    velocity_after_m_s: exact.velocity_after_m_s,
                    spin_rate_after_rad_s: exact.spin_rate_after_rad_s,
                    force_impulse_ns: exact.force_impulse_ns,
                    spin_angular_impulse_nms: exact.spin_angular_impulse_nms,
                    slip_increment_m: exact.slip_increment_m,
                    spin_increment_rad: exact.spin_increment_rad,
                    dissipated_energy_j: exact.dissipated_energy_j,
                    prefix_dissipated_energy_j: exact.dissipated_energy_j,
                    stopping_tail_heat_j: 0.0,
                    contact_work_j: exact.contact_work_j,
                    work_residual_j: exact.work_residual_j,
                    prefix_work_residual_j: exact.work_residual_j,
                    estimated_scaled_rate_error_m_s: 0.0,
                    estimated_scaled_impulse_error_ns: 0.0,
                    estimated_scaled_travel_error_m: 0.0,
                    estimated_heat_error_j: 0.0,
                    tail_slip_bound_m: [0.0; 2],
                    tail_spin_bound_rad: 0.0,
                    stopping_time_bracket_s: exact.kinetic_phase_end_s.map(|t| [t, t]),
                    integrated_until_s: exact.kinetic_phase_end_s.unwrap_or(input.dt_s),
                    scaled_velocity_residual_m_s: exact.scaled_velocity_residual_m_s,
                    exact_ray: true,
                    accepted_steps: 0,
                    rejected_steps: 0,
                    evaluations: 1,
                });
            }
            Err(HertzRelaxationError::CurvedTrajectory { .. }) => {}
            Err(HertzRelaxationError::Cancelled) => return Err(HertzTransientError::Cancelled),
            Err(error) => return Err(HertzTransientError::Exact(error)),
        }
        let upper = det::sqrt(1.0 + SPIN_RADIUS * SPIN_RADIUS);
        advance_circular_transient(
            &CircularRelaxationInput {
                velocity_m_s: input.velocity_m_s,
                spin_rate_rad_s: input.spin_rate_rad_s,
                scaled_mobility_per_kg: input.scaled_mobility_per_kg,
                dt_s: input.dt_s,
            },
            config,
            controls,
            CircularTransientLaw {
                radius: self.parameters().contact_radius_m,
                sliding_capacity_n: self.kinetic_patch.sliding_capacity_n(),
                spinning_capacity_nm: self.kinetic_patch.spinning_capacity_nm(),
                power_norm_bounds: [
                    (SPIN_RADIUS / upper) * (1.0 - 32.0 * f64::EPSILON),
                    upper * (1.0 + 32.0 * f64::EPSILON),
                ],
                quadrature: CircularTransientQuadrature::Rk4Doubling,
                evaluate: |velocity, spin| {
                    let native = self
                        .kinetic_patch
                        .evaluate(velocity, spin)
                        .map_err(HertzWrenchError::Law)?;
                    Ok(CircularTransientTraction {
                        force_n: native.force_n,
                        spin_moment_nm: native.spin_moment_nm,
                        dissipated_power_w: native.dissipated_power_w,
                        at_rest: native.at_rest,
                    })
                },
            },
            cancelled,
        )
        .map_err(TransientError::with_exact)
    }
}
