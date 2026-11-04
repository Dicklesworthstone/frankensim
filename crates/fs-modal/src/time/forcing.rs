//! Exact first-order-hold forcing for the shared scalar oscillator.
//!
//! A force record supplies values at interval endpoints. This owner integrates
//! their linear interpolant, retaining the phase of the continuous oscillator
//! instead of introducing an additional time-discretization frequency shift.
//! The interpolation remains a model of the supplied record, not a claim about
//! unresolved forcing between measurements.

use super::{DampedTransition, ModalAcousticState, ModalAcousticTimeError};
use super::heat::LinearViscousIntegral;

/// Cached transition for `q'' + c*q' + omega^2*q = F(t)` with
/// `F(t) = F_start + (F_end-F_start)*t/duration`.
/// Supply `c = 2*zeta*omega` through [`Self::try_new`], or specify the physical
/// damping rate directly through [`Self::try_new_with_damping_rate`].
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
/// [`Self::advance_with_energy`] also integrates force work and nonnegative
/// physical viscous heat along that same continuous trajectory. No acoustic
/// observation or certified floating-point enclosure is inferred. Exact
/// refers to the linear differential equation under this forcing interpolant;
/// finite arithmetic still contributes roundoff, including phase error for
/// very long oscillatory intervals.
#[derive(Clone, Copy, Debug)]
pub struct LinearForceTransition {
    duration_s: f64,
    angular_frequency_rad_s: f64,
    transition: DampedTransition,
    ramp_displacement: f64,
    ramp_integral: f64,
    viscous_integral: LinearViscousIntegral,
}

/// Endpoint and independent physical work integrals for one linear force ramp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearForceFrame {
    /// Motion at the interval's endpoint in the original modal coordinates.
    pub state: ModalAcousticState,
    /// Integral of `F(t)*v(t)` [J], including either sign of supplied work.
    pub input_work_j: f64,
    /// Integral of `c*v(t)^2` [J], evaluated as a sum of squares.
    pub viscous_dissipation_j: f64,
    /// Scale-aware floating-point allowance for the viscous integral [J].
    pub dissipation_roundoff_tolerance_j: f64,
}

/// Cotangents of a scalar objective through one fixed oscillator interval.
///
/// These are derivatives, not a physical modal state. Each has the objective's
/// units divided by the corresponding input's units; callers supply the scalar
/// objective and any additional parameterization chain rules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearForcePullback {
    /// Derivative with respect to initial displacement, dJ/dq_start.
    pub initial_displacement: f64,
    /// Derivative with respect to initial velocity, dJ/dv_start.
    pub initial_velocity: f64,
    /// Derivative with respect to the first endpoint force, dJ/dF_start.
    pub force_start: f64,
    /// Derivative with respect to the second endpoint force, dJ/dF_end.
    pub force_end: f64,
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
        let angle = angular_frequency_rad_s * duration_s;
        let drag = 2.0 * (damping_ratio * angle);
        Self::prepare(angular_frequency_rad_s, duration_s, angle, drag)
    }

    /// Prepare `q'' + c*q' + omega^2*q = F(t)` with an explicit damping rate
    /// `c` [1/s]. For a dimensional coordinate of mass `m` with viscous force
    /// `-b*v`, supply `c = b/m` and mass-normalize the state and forcing.
    ///
    /// Unlike a damping ratio, this remains defined at zero stiffness. A free
    /// coordinate with positive `c` decays under viscous drag and dissipates
    /// physical heat through [`Self::advance_with_energy`]. Nonzero `omega`
    /// uses the same transition; no small artificial stiffness or division by
    /// the natural frequency is introduced.
    ///
    /// # Errors
    /// Invalid parameters or unrepresentable dimensionless transition values.
    pub fn try_new_with_damping_rate(
        angular_frequency_rad_s: f64,
        damping_rate_per_s: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        if ![angular_frequency_rad_s, damping_rate_per_s]
            .into_iter()
            .all(|value| value.is_finite() && value >= 0.0)
            || !duration_s.is_finite()
            || duration_s <= 0.0
        {
            return Err(invalid(
                "linear-force transition requires finite nonnegative frequency/damping rate and positive duration",
            ));
        }
        Self::prepare(
            angular_frequency_rad_s,
            duration_s,
            angular_frequency_rad_s * duration_s,
            damping_rate_per_s * duration_s,
        )
    }

    fn prepare(
        angular_frequency_rad_s: f64,
        duration_s: f64,
        mut angle: f64,
        mut drag: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        if !angle.is_finite() || !drag.is_finite() {
            return Err(invalid(
                "linear-force dimensionless frequency/damping must be representable",
            ));
        }
        let original_angle = angle;
        let original_drag = drag;
        let mut doublings = 0;
        while angle.max(drag) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            doublings += 1;
        }
        let mut transition = DampedTransition::small(angle, drag);
        let (mut ramp_displacement, mut ramp_integral) = small_ramp_moments(angle, drag);
        for _ in 0..doublings {
            // S(h)=integral_0^h R(t) dt. Its fourth scaled moment is
            // needed by integral t*v(t), hence the physical ramp work.
            let scaled_ramp = transition.angle * ramp_displacement;
            ramp_integral = (2.0 * ramp_integral + 2.0 * ramp_displacement
                - scaled_ramp * scaled_ramp
                + transition.g * transition.g) * 0.0625;
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
            ramp_integral,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid(
                "linear-force transition coefficients must be representable",
            ));
        }
        let viscous_integral = LinearViscousIntegral::new(original_angle, original_drag);
        if !viscous_integral.is_finite() {
            return Err(invalid("linear-force viscous integral must be representable"));
        }
        Ok(Self {
            duration_s,
            angular_frequency_rad_s,
            transition,
            ramp_displacement,
            ramp_integral,
            viscous_integral,
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
        self.motion(state, force_start, force_end).map(|(next, _, _)| next)
    }

    fn motion(
        &self,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
    ) -> Result<(ModalAcousticState, f64, f64), ModalAcousticTimeError> {
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
        Ok((ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        }, residual, scaled_change))
    }

    /// Advance with independently integrated supplied work and physical heat.
    ///
    /// The viscous integral is a prepared positive square-root Gramian, so
    /// high-Q heat smaller than stored-energy subtraction roundoff survives.
    /// Exact force moments integrate the ramp's work without subtracting two
    /// large positions. Preparation does all quadrature, composition and roots;
    /// this call allocates nothing and uses only fixed scalar arithmetic.
    ///
    /// # Errors
    /// Every refusal from [`Self::advance`], or nonfinite work/heat arithmetic.
    pub fn advance_with_energy(
        &self,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
    ) -> Result<LinearForceFrame, ModalAcousticTimeError> {
        let (next, residual, change) = self.motion(state, force_start, force_end)?;
        let v = state.velocity_m_sqrt_kg_per_s;
        let flow = self.transition;
        // Integrals of v and (t/h)*v, each divided by h. All operands are
        // velocities; multiplying by h*F avoids overflow in the force slope.
        let mean_velocity = flow.b * v + flow.g * residual + self.ramp_displacement * change;
        let first_moment = (flow.b - flow.g) * v
            + (flow.g - self.ramp_displacement) * residual
            + (self.ramp_displacement - self.ramp_integral) * change;
        let work = (self.duration_s * force_start) * mean_velocity + change * first_moment;
        let (heat, allowance) = self.viscous_integral.evaluate(v, residual, change, flow.angle);
        if ![work, heat, allowance].into_iter().all(f64::is_finite) {
            return Err(invalid("linear-force work and viscous heat must be representable"));
        }
        Ok(LinearForceFrame {
            state: next,
            input_work_j: work,
            viscous_dissipation_j: heat,
            dissipation_roundoff_tolerance_j: allowance,
        })
    }

    /// Pull endpoint cotangents `(dJ/dq_end, dJ/dv_end)` through this linear map.
    ///
    /// Frequency, damping ratio and duration are fixed. The transpose reuses
    /// the prepared flow and force moments, with no matrix, allocation,
    /// transcendental call or differentiation through preparation iterations.
    /// Small state corrections use the same independently stored decrements as
    /// the forward map; long decays retain the same homogeneous tails.
    ///
    /// For a force sample shared by two adjacent intervals, sum this interval's
    /// `force_end` and the next interval's `force_start`. A trajectory reverse
    /// sweep must propagate both state cotangents. Derivatives of physical
    /// coefficients, sample times and nonsmooth peak selection are not supplied.
    ///
    /// # Errors
    /// Nonfinite endpoint cotangents or unrepresentable scaled products/results.
    /// No input or prepared transition is mutated and no partial gradient returns.
    pub fn pullback(
        &self,
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
    ) -> Result<LinearForcePullback, ModalAcousticTimeError> {
        let (q_bar, v_bar) = (endpoint_displacement_cotangent, endpoint_velocity_cotangent);
        if ![q_bar, v_bar].into_iter().all(f64::is_finite) {
            return Err(invalid("linear-force endpoint cotangents must be finite"));
        }
        let flow = self.transition;
        let h = self.duration_s;
        let restoring_gain = self.angular_frequency_rad_s * (flow.angle * flow.b);
        let restoring = if restoring_gain.is_finite() && restoring_gain != 0.0 {
            restoring_gain * v_bar
        } else {
            // Avoid squaring omega; a finite seed can make a very large/small
            // homogeneous gain representable after the cotangent is included.
            self.angular_frequency_rad_s * ((flow.angle * flow.b) * v_bar)
        };
        let initial_displacement = if flow.position_decrement.abs() < 0.5 {
            q_bar - (flow.position_decrement * q_bar + restoring)
        } else {
            flow.position_retention * q_bar - restoring
        };
        let scaled_q_bar = h * q_bar;
        let initial_velocity = if flow.velocity_decrement.abs() < 0.5 {
            v_bar + (flow.b * scaled_q_bar - flow.velocity_decrement * v_bar)
        } else {
            flow.b * scaled_q_bar + flow.velocity_retention * v_bar
        };
        let force_start =
            h * ((flow.g - self.ramp_displacement) * scaled_q_bar + (flow.b - flow.g) * v_bar);
        let force_end = h * (self.ramp_displacement * scaled_q_bar + flow.g * v_bar);
        if ![
            restoring,
            scaled_q_bar,
            initial_displacement,
            initial_velocity,
            force_start,
            force_end,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid(
                "linear-force scaled cotangents and pullback must be representable",
            ));
        }
        Ok(LinearForcePullback {
            initial_displacement,
            initial_velocity,
            force_start,
            force_end,
        })
    }
}

fn small_ramp_moments(angle: f64, drag: f64) -> (f64, f64) {
    let stiffness = angle * angle;
    let mut previous = 0.0;
    let mut power = 1.0;
    let mut reciprocal_factorial = 1.0 / 6.0;
    let mut ramp = reciprocal_factorial;
    let mut integral = reciprocal_factorial * 0.25;
    // The existing impulse-response series has p_j/(j+1)!. Two time
    // integrations give p_j/(j+3)! for R/h^3 on the resolved interval.
    for j in 1..=18 {
        let next = -drag * power - stiffness * previous;
        previous = power;
        power = next;
        reciprocal_factorial /= f64::from(j + 3);
        ramp += power * reciprocal_factorial;
        integral += power * reciprocal_factorial / f64::from(j + 4);
    }
    (ramp, integral)
}

fn invalid(what: &'static str) -> ModalAcousticTimeError {
    ModalAcousticTimeError::InvalidInput { what }
}
