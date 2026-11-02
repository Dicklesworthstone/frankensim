//! Exact sinusoidal forcing of the shared mass-normalized oscillator.
//!
//! The force is a continuous physical tone, not held samples of that tone.
//! An entire convolution and block composition retain secular resonant growth
//! without dividing by detuning or subtracting nearly equal oscillator poles.

use super::{DampedTransition, ModalAcousticState, ModalAcousticTimeError};
use fs_math::det;

/// Cached motion for `q'' + 2*zeta*omega*q' + omega^2*q = F(t)`, where
/// `F(t) = cosine_force*cos(Omega*t) + sine_force*sin(Omega*t)` and time is
/// measured from the supplied initial state. Forces and state use the same
/// mass-normalized units as [`ModalAcousticState`]. Under `exp(-i Omega t)`,
/// the complex force amplitude is `cosine_force + i*sine_force`.
///
/// Preparation uses the existing scalar homogeneous transition and an entire
/// force convolution. Undamped resonance, both sides of critical damping and
/// zero-frequency free coordinates therefore need no singular special case.
/// Once prepared, all methods allocate nothing and call no transcendentals.
/// Multiple harmonics superpose through [`Self::forced_increment`]; advance
/// the homogeneous state only once before summing their zero-state responses.
///
/// This API supplies motion, not harmonic-balance convergence, force work,
/// viscous heat, or a broadband acoustic observation. Exact refers to this
/// linear differential equation and continuous forcing law. Finite arithmetic,
/// unresolved intermediate observations and nonlinear generated harmonics
/// retain their own accuracy requirements. A sampled force record should use
/// its actual interpolation law instead of being silently treated as a tone.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicForceTransition {
    duration_s: f64,
    angular_frequency_rad_s: f64,
    transition: DampedTransition,
    // [q_cos/h^2, q_sin/h^2, v_cos/h, v_sin/h]. When remainder is true,
    // the cosine column integrates cos(Omega*t)-1; the ZOH owner supplies 1.
    force: [f64; 4],
    cosine_remainder: bool,
    phase_cosine: f64,
    phase_sine: f64,
}

impl HarmonicForceTransition {
    /// Prepare a positive-duration interval with finite nonnegative natural
    /// and forcing frequencies and damping. A free coordinate (`omega = 0`)
    /// requires zero damping ratio, matching the shared scalar owner.
    /// No Nyquist cutoff is imposed: this computes an endpoint of a continuous
    /// prescribed load, and the caller owns observation/time resolution.
    ///
    /// # Errors
    /// Invalid parameters or unrepresentable dimensionless coefficients.
    pub fn try_new(
        angular_frequency_rad_s: f64,
        damping_ratio: f64,
        forcing_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        if ![
            angular_frequency_rad_s,
            damping_ratio,
            forcing_angular_frequency_rad_s,
        ]
        .into_iter()
        .all(|value| value.is_finite() && value >= 0.0)
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || (angular_frequency_rad_s == 0.0 && damping_ratio != 0.0)
        {
            return Err(invalid(
                "harmonic transition requires finite nonnegative frequencies/damping, positive duration, and zero damping for a free coordinate",
            ));
        }
        let mut angle = angular_frequency_rad_s * duration_s;
        let mut drag = 2.0 * (damping_ratio * angle);
        let full_phase = forcing_angular_frequency_rad_s * duration_s;
        if ![angle, drag, full_phase].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "harmonic dimensionless frequencies/damping must be representable",
            ));
        }
        let mut phase = full_phase;
        let mut doublings = 0;
        while angle.max(drag).max(phase) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            phase *= 0.5;
            doublings += 1;
        }
        // A separately evaluated cos-1 convolution retains the tiny motion of
        // an almost-static drive about equilibrium. At larger phase use the
        // complete gains, avoiding cancellation against the held-force gain.
        let cosine_remainder = full_phase <= 0.25;
        let mut transition = DampedTransition::small(angle, drag);
        let mut force = small_force(angle, drag, phase, cosine_remainder);
        for _ in 0..doublings {
            let (sine, cosine) = det::sin_cos(phase);
            let [qc, qs, vc, vs] = force;
            let restoring_gain = -(transition.angle * transition.b) * transition.angle;
            let (constant_cos, constant_sin) = if cosine_remainder {
                let half_sine = det::sin(0.5 * phase);
                (-2.0 * half_sine * half_sine, sine)
            } else {
                (0.0, 0.0)
            };
            // U(2h) = E(h) U(h) + U(h) R(Omega*h). The position
            // gains are scaled by h^2 and velocity gains by h. For a cos-1
            // column, add U_ZOH (R-I) before rescaling to the doubled h.
            force = [
                0.25 * (transition.position_retention * qc + transition.b * vc + qc * cosine
                    - qs * sine
                    + transition.g * constant_cos),
                0.25 * (transition.position_retention * qs
                    + transition.b * vs
                    + qc * sine
                    + qs * cosine
                    + transition.g * constant_sin),
                0.5 * (restoring_gain * qc + transition.velocity_retention * vc + vc * cosine
                    - vs * sine
                    + transition.b * constant_cos),
                0.5 * (restoring_gain * qs
                    + transition.velocity_retention * vs
                    + vc * sine
                    + vs * cosine
                    + transition.b * constant_sin),
            ];
            transition.double();
            phase *= 2.0;
        }
        if !force.into_iter().all(f64::is_finite)
            || ![
                transition.b,
                transition.g,
                transition.position_retention,
                transition.velocity_retention,
                transition.position_decrement,
                transition.velocity_decrement,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "harmonic transition coefficients must be representable",
            ));
        }
        let (phase_sine, phase_cosine) = det::sin_cos(full_phase);
        Ok(Self {
            duration_s,
            angular_frequency_rad_s,
            transition,
            force,
            cosine_remainder,
            phase_cosine,
            phase_sine,
        })
    }

    /// Advance a boundary state under one continuous force tone. No accepted
    /// state is retained or mutated; a failed trial cannot advance motion.
    ///
    /// # Errors
    /// Nonfinite inputs or unrepresentable scaled force/state arithmetic.
    pub fn advance(
        &self,
        state: ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        let q0 = state.displacement_m_sqrt_kg;
        let v0 = state.velocity_m_sqrt_kg_per_s;
        if ![q0, v0, cosine_force, sine_force]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "harmonic state and force coefficients must be finite",
            ));
        }
        let h = self.duration_s;
        let flow = self.transition;
        let fc = h * cosine_force;
        let fs = h * sine_force;
        let restoring = (self.angular_frequency_rad_s * q0) * flow.angle;
        let [qc, qs, vc, vs] = self.force;
        let driven_q = qc * fc + qs * fs;
        let driven_v = vc * fc + vs * fs;
        let residual = if self.cosine_remainder {
            fc - restoring
        } else {
            -restoring
        };
        let q = if flow.position_decrement.abs() < 0.5 {
            q0 + h * (flow.b * v0 + flow.g * residual + driven_q)
        } else {
            let constant = if self.cosine_remainder {
                flow.g * fc
            } else {
                0.0
            };
            flow.position_retention * q0 + h * (flow.b * v0 + constant + driven_q)
        };
        let velocity_increment = flow.b * residual + driven_v;
        let v = if flow.velocity_decrement.abs() < 0.5 {
            v0 + (velocity_increment - flow.velocity_decrement * v0)
        } else {
            flow.velocity_retention * v0 + velocity_increment
        };
        if ![fc, fs, restoring, q, v].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "harmonic scaled inputs and response must be representable",
            ));
        }
        Ok(ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        })
    }

    /// Zero-state response contributed by this harmonic over the interval.
    /// Sum these contributions for several frequencies after advancing the
    /// homogeneous state once with zero forces. Never sum several copies of
    /// the same initial state. This operation is linear in both force entries.
    ///
    /// # Errors
    /// The same finite-force and representability refusals as [`Self::advance`].
    pub fn forced_increment(
        &self,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        self.advance(ModalAcousticState::default(), cosine_force, sine_force)
    }

    /// Shift interval-local force coefficients to the next interval without
    /// trigonometric calls: `(Fc cos + Fs sin, Fs cos - Fc sin)`. This carries
    /// the continuous tone's phase; it does not reset it at each sample.
    /// Repeated floating-point rotations accumulate roundoff. A caller may
    /// re-evaluate phase from its own absolute clock when that error matters;
    /// no amplitude normalization or artificial loss is applied here.
    ///
    /// # Errors
    /// Nonfinite force coefficients or an unrepresentable rotated force.
    pub fn advance_force_phase(
        &self,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<(f64, f64), ModalAcousticTimeError> {
        let cosine = self.phase_cosine * cosine_force + self.phase_sine * sine_force;
        let sine = self.phase_cosine * sine_force - self.phase_sine * cosine_force;
        if ![cosine_force, sine_force, cosine, sine]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("harmonic force phase must remain finite"));
        }
        Ok((cosine, sine))
    }
}

fn invalid(what: &'static str) -> ModalAcousticTimeError {
    ModalAcousticTimeError::InvalidInput { what }
}

fn small_force(angle: f64, drag: f64, phase: f64, cosine_remainder: bool) -> [f64; 4] {
    let mut gains = [0.0; 4];
    for column in 0..2 {
        // Exp(action) on [q, h*v, h^2*Fc, h^2*Fs]. All quantities are
        // dimensionless; term_n = G term_(n-1)/n with zero initial motion.
        let mut term = [
            0.0,
            0.0,
            if column == 0 { 1.0 } else { 0.0 },
            if column == 1 { 1.0 } else { 0.0 },
        ];
        for order in 1..=24 {
            let inverse = 1.0 / f64::from(order);
            let source = if cosine_remainder && column == 0 && order == 1 {
                0.0
            } else {
                term[2]
            };
            term = [
                term[1] * inverse,
                (-angle * angle * term[0] - drag * term[1] + source) * inverse,
                phase * term[3] * inverse,
                -phase * term[2] * inverse,
            ];
            gains[column] += term[0];
            gains[2 + column] += term[1];
        }
    }
    gains
}
