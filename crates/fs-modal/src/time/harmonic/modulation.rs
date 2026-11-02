//! Phase-integrated reconstruction of a prescribed changing-frequency force.
use super::{ModalAcousticTimeError, det, invalid};

/// One immutable force-modulation proposal on a caller-owned physical clock.
///
/// The desired force has a linearly interpolated signed gain and instantaneous
/// angular frequency over the interval. An existing affine tone propagator
/// integrates its reconstructed cosine/sine envelopes about a fixed carrier.
/// This preparation supplies those endpoints and integrates the desired phase
/// exactly; it does not replace physical mode frequencies or alter material.
///
/// The envelope reconstruction is second order in duration, not an exact chirp
/// solver. [`Self::coefficient_error_bound`] bounds that interpolation defect in
/// real arithmetic. Roundoff and accumulated caller-owned phase error are not
/// enclosed, so this estimate is not an interval certificate. The caller owns
/// source/observation bandwidth admission and the acceptable physical force
/// error. Multiplying the bound by `hypot(cosine, sine)` gives an absolute force
/// error estimate in the original coefficient units.
///
/// Construction uses deterministic trigonometry and allocates nothing. Both
/// evaluation methods are pure, allocate nothing and call no transcendentals.
/// Publish phase only after the associated physical and observer step accepts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ToneModulationStep {
    gain: [f64; 2],
    relative_phase: [f64; 2],
    force_phase: [f64; 2],
    error_bound: f64,
}

impl ToneModulationStep {
    /// Prepare one positive-duration interval. Frequencies are finite and
    /// nonnegative in rad/s; gains are finite signed dimensionless multipliers.
    ///
    /// # Errors
    /// Invalid inputs or unrepresentable phase, interpolation-bound arithmetic,
    /// or deterministic trigonometric values. No Nyquist limit is guessed here;
    /// source users must apply their original physical clock's admission rule.
    pub fn try_new(
        duration_s: f64,
        carrier_omega_rad_s: f64,
        instantaneous_omega_rad_s: [f64; 2],
        gain: [f64; 2],
    ) -> Result<Self, ModalAcousticTimeError> {
        let [w0, w1] = instantaneous_omega_rad_s;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || ![carrier_omega_rad_s, w0, w1]
                .into_iter()
                .all(|w| w.is_finite() && w >= 0.0)
            || !gain.into_iter().all(f64::is_finite)
        {
            return Err(invalid(
                "tone modulation requires positive finite duration, nonnegative finite frequencies and finite gains",
            ));
        }
        // Subtract frequencies before scaling, retaining small detuning about a
        // large carrier. Nonnegative frequency differences cannot overflow.
        let d0 = (w0 - carrier_omega_rad_s) * duration_s;
        let d1 = (w1 - carrier_omega_rad_s) * duration_s;
        let frequency_change = (w1 - w0) * duration_s;
        let absolute_angle = average(w0, w1) * duration_s;
        if ![d0, d1, frequency_change, absolute_angle]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("tone modulation phase must be representable"));
        }
        let d = d0.abs().max(d1.abs());
        let g = gain[0].abs().max(gain[1].abs());
        // For z=g*exp(-i theta), theta'=omega-carrier:
        // |z''| <= 2*|g'|*max|theta'| + max|g|*(|theta''|+max|theta'|^2).
        // Linear interpolation has a sup-norm error <= h^2*sup|z''|/8.
        // Evaluate dimensionless terms without first squaring a large d or
        // overflowing the difference between opposite large signed gains.
        let error_bound = if w0 == carrier_omega_rad_s && w1 == carrier_omega_rad_s {
            0.0 // The existing affine gain law is recovered exactly.
        } else {
            let change = gain[1] - gain[0];
            let mixed = if change.is_finite() {
                scaled_product(change.abs(), d, 0.25)
            } else {
                (gain[1] * 0.25 - gain[0] * 0.25).abs() * d
            };
            let curvature = scaled_product(g, frequency_change.abs(), 0.125);
            let rotation = if (1.0..8.0).contains(&d) {
                scaled_product(g, d * d, 0.125)
            } else {
                scaled_product(g, d, 0.125) * d
            };
            mixed + curvature + rotation
        };
        if !error_bound.is_finite() {
            return Err(invalid(
                "tone modulation coefficient error bound must be representable",
            ));
        }
        let relative_phase = rotation(average(d0, d1))?;
        let force_phase = rotation(absolute_angle)?;
        Ok(Self {
            gain,
            relative_phase,
            force_phase,
            error_bound,
        })
    }

    /// Absolute dimensionless complex-envelope interpolation bound. Multiply
    /// by the unscaled original force-coefficient norm for force units.
    #[must_use]
    pub const fn coefficient_error_bound(&self) -> f64 {
        self.error_bound
    }

    /// Return `[cosine_start, sine_start, cosine_end, sine_end]` for the
    /// existing affine tone propagator, all in its SAME interval-start carrier
    /// basis. Inputs retain the desired source's unscaled phase at this start.
    /// A zero end gain does not destroy the phase needed by a later re-attack.
    ///
    /// # Errors
    /// Nonfinite input coefficients or unrepresentable rotated/scaled values.
    pub fn endpoints(&self, cosine: f64, sine: f64) -> Result<[f64; 4], ModalAcousticTimeError> {
        let (end_cosine, end_sine) = rotate(cosine, sine, self.relative_phase)?;
        let endpoints = [
            cosine * self.gain[0],
            sine * self.gain[0],
            end_cosine * self.gain[1],
            end_sine * self.gain[1],
        ];
        if !endpoints.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "modulated physical force endpoints must be representable",
            ));
        }
        Ok(endpoints)
    }

    /// Advance the unscaled force by the exact integral of the desired
    /// linearly varying frequency. Call only after accepting the physical step;
    /// pass the original coefficients, not gain-scaled envelope endpoints.
    /// Repeated rotations accumulate ordinary floating-point phase error.
    ///
    /// # Errors
    /// Nonfinite input coefficients or an unrepresentable rotated force.
    pub fn advance_force_phase(
        &self,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), ModalAcousticTimeError> {
        rotate(cosine, sine, self.force_phase)
    }
}

fn average(a: f64, b: f64) -> f64 {
    if a == b { a } else { a * 0.5 + b * 0.5 }
}

// Divide the larger factor first: this avoids avoidable overflow without
// underflowing a tiny factor before multiplication by a large one.
fn scaled_product(a: f64, b: f64, scale: f64) -> f64 {
    if a >= b {
        (a * scale) * b
    } else {
        a * (b * scale)
    }
}

fn rotation(angle: f64) -> Result<[f64; 2], ModalAcousticTimeError> {
    let (sine, cosine) = if angle == 0.0 {
        (0.0, 1.0)
    } else {
        det::sin_cos(angle)
    };
    if ![sine, cosine].into_iter().all(f64::is_finite) {
        return Err(invalid("tone modulation rotation must remain finite"));
    }
    Ok([sine, cosine])
}

fn rotate(
    cosine: f64,
    sine: f64,
    [phase_sine, phase_cosine]: [f64; 2],
) -> Result<(f64, f64), ModalAcousticTimeError> {
    if ![cosine, sine].into_iter().all(f64::is_finite) {
        return Err(invalid("tone modulation coefficients must remain finite"));
    }
    if phase_sine == 0.0 && phase_cosine == 1.0 {
        return Ok((cosine, sine));
    }
    let next_cosine = phase_cosine * cosine + phase_sine * sine;
    let next_sine = phase_cosine * sine - phase_sine * cosine;
    if ![next_cosine, next_sine].into_iter().all(f64::is_finite) {
        return Err(invalid("modulated force phase must be representable"));
    }
    Ok((next_cosine, next_sine))
}
