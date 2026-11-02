//! Exact first-order-hold forcing for the shared scalar oscillator.
//!
//! A force record supplies values at interval endpoints. This owner integrates
//! their linear interpolant, retaining the phase of the continuous oscillator
//! instead of introducing an additional time-discretization frequency shift.
//! The interpolation remains a model of the supplied record, not a claim about
//! unresolved forcing between measurements.

use super::{DampedTransition, ModalAcousticState, ModalAcousticTimeError};

/// Cached transition for `q'' + 2*zeta*omega*q' + omega^2*q = F(t)` with
/// `F(t) = F_start + (F_end-F_start)*t/duration`.
///
/// Coordinates and forces use the same mass-normalized units as
/// [`ModalAcousticState`]. This is also the unit-mass acceleration-response
/// equation. Preparation is independent of state and forcing; advancing a
/// prepared interval allocates nothing and calls no transcendental functions.
///
/// The entire-series/doubling transition is shared with the modal ZOH owner.
/// Its additional force moment avoids dividing by stiffness, a damped natural
/// frequency, or the separation of two poles. Undamped, critical, overdamped
/// and explicitly free coordinates therefore share one continuous formula.
///
/// This API returns motion only: it does not publish an acoustic observation,
/// force work, viscous heat, or a certified floating-point enclosure. Exact
/// refers to the linear differential equation under this forcing interpolant;
/// finite arithmetic still contributes roundoff, including phase error for
/// very long oscillatory intervals. Existing ZOH energy accounts are unchanged.
#[derive(Clone, Copy, Debug)]
pub struct LinearForceTransition {
    duration_s: f64,
    angular_frequency_rad_s: f64,
    transition: DampedTransition,
    ramp_displacement: f64,
}

impl LinearForceTransition {
    /// Prepare one positive-duration interval with nonnegative finite damping.
    /// A zero frequency explicitly denotes a free coordinate and requires zero
    /// damping ratio; no drag coefficient is inferred from that singular ratio.
    /// There is no sample-rate/Nyquist admission here: callers must resolve their
    /// physical forcing and any requested intermediate observations themselves.
    ///
    /// # Errors
    /// Invalid parameters or unrepresentable dimensionless transition values.
    pub fn try_new(
        angular_frequency_rad_s: f64,
        damping_ratio: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        if !angular_frequency_rad_s.is_finite()
            || angular_frequency_rad_s < 0.0
            || !damping_ratio.is_finite()
            || damping_ratio < 0.0
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || (angular_frequency_rad_s == 0.0 && damping_ratio != 0.0)
        {
            return Err(invalid(
                "linear-force transition requires finite nonnegative frequency/damping, positive duration, and zero damping for a free coordinate",
            ));
        }
        let mut angle = angular_frequency_rad_s * duration_s;
        let mut drag = 2.0 * (damping_ratio * angle);
        if !angle.is_finite() || !drag.is_finite() {
            return Err(invalid(
                "linear-force dimensionless frequency/damping must be representable",
            ));
        }
        let mut doublings = 0;
        while angle.max(drag) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            doublings += 1;
        }
        let mut transition = DampedTransition::small(angle, drag);
        let mut ramp_displacement = small_ramp_displacement(angle, drag);
        for _ in 0..doublings {
            // R(h)=integral_0^h G(t) dt, where G is the held-force
            // displacement response. Composition of two linear-force halves:
            // R(2h)=(1+E11)*R(h)+B(h)*G(h)+h*G(h).
            // Keep R/h^3 scaled; it multiplies h^2*(F_end-F_start).
            ramp_displacement = ((1.0 + transition.position_retention) * ramp_displacement
                + (1.0 + transition.b) * transition.g)
                * 0.125;
            transition.double();
        }
        if ![
            transition.b,
            transition.g,
            transition.position_retention,
            transition.velocity_retention,
            transition.position_decrement,
            transition.velocity_decrement,
            ramp_displacement,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid(
                "linear-force transition coefficients must be representable",
            ));
        }
        Ok(Self {
            duration_s,
            angular_frequency_rad_s,
            transition,
            ramp_displacement,
        })
    }

    /// Advance from a supplied boundary state under the two endpoint forces.
    /// No state is retained or mutated, so rejected trials and checkpoint replay
    /// cannot advance accepted motion accidentally.
    ///
    /// Small increments are formed about the current displacement/velocity.
    /// Long decays instead retain the independently computed homogeneous tails;
    /// subtracting a nearly equal whole state would erase those tails.
    ///
    /// # Errors
    /// Nonfinite state/forces or unrepresentable scaled forcing or response.
    pub fn advance(
        &self,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        let q0 = state.displacement_m_sqrt_kg;
        let v0 = state.velocity_m_sqrt_kg_per_s;
        if ![q0, v0, force_start, force_end]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "linear-force state and endpoint forces must be finite",
            ));
        }
        let h = self.duration_s;
        let flow = self.transition;
        let force_change = force_end - force_start;
        let scaled_change = if force_change.is_finite() {
            h * force_change
        } else {
            // Opposite finite endpoints can overflow their difference even
            // though a short interval has a representable velocity impulse.
            h * force_end - h * force_start
        };
        let scaled_force = h * force_start;
        // Form h*omega^2*q without first squaring a very soft/stiff omega.
        let scaled_restoring = (self.angular_frequency_rad_s * q0) * flow.angle;
        let residual = scaled_force - scaled_restoring;
        let driven_velocity = flow.b * residual + flow.g * scaled_change;
        let q = if flow.position_decrement.abs() < 0.5 {
            q0 + h * (flow.b * v0 + flow.g * residual + self.ramp_displacement * scaled_change)
        } else {
            flow.position_retention * q0
                + h * (flow.b * v0 + flow.g * scaled_force + self.ramp_displacement * scaled_change)
        };
        let v = if flow.velocity_decrement.abs() < 0.5 {
            v0 + (driven_velocity - flow.velocity_decrement * v0)
        } else {
            flow.velocity_retention * v0 + driven_velocity
        };
        if ![scaled_change, scaled_force, scaled_restoring, q, v]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "linear-force scaled inputs and response must be representable",
            ));
        }
        Ok(ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        })
    }
}

fn small_ramp_displacement(angle: f64, drag: f64) -> f64 {
    let stiffness = angle * angle;
    let mut previous = 0.0;
    let mut power = 1.0;
    let mut reciprocal_factorial = 1.0 / 6.0;
    let mut ramp = reciprocal_factorial;
    // The existing impulse-response series has p_j/(j+1)!. Two time
    // integrations give p_j/(j+3)! for R/h^3 on the resolved interval.
    for j in 1..=18 {
        let next = -drag * power - stiffness * previous;
        previous = power;
        power = next;
        reciprocal_factorial /= f64::from(j + 3);
        ramp += power * reciprocal_factorial;
    }
    ramp
}

fn invalid(what: &'static str) -> ModalAcousticTimeError {
    ModalAcousticTimeError::InvalidInput { what }
}
