//! Prescribed source clock derivatives through the original tone convolution.

use super::{
    HarmonicForceTransition, ModalAcousticTimeError, PreparedHarmonicForceSensitivity, invalid,
};

/// Transpose of the SAME continuous source-phase rotation used by
/// [`HarmonicForceTransition::advance_force_phase`]. Combine these force-input
/// cotangents with the motion pullback at each step of a reverse trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicForcePhasePullback {
    /// dJ/dFc at the interval start, from the two advanced source coefficients.
    pub cosine_force: f64,
    /// dJ/dFs at the interval start, from the two advanced source coefficients.
    pub sine_force: f64,
    /// dJ/dOmega from phase advancement alone [objective seconds/radian].
    /// Add the motion pullback's forcing-frequency contribution separately.
    pub forcing_angular_frequency: f64,
}

impl HarmonicForceTransition {
    /// Pull back the advanced source quadratures AND their physical clock.
    ///
    /// For R = [[cos,sin],[-sin,cos]], the initial-force cotangent is R^T times
    /// the supplied endpoint-force cotangent. The forcing-frequency contribution
    /// is h*(next_sine*bar_next_cosine-next_cosine*bar_next_sine), since h is
    /// fixed. This source-state chain rule is required when many intervals share
    /// one authored frequency: each motion partial otherwise keeps its local
    /// input phase fixed. The final source cotangent is zero if the objective
    /// depends only on physical motion, then accumulates through the reverse sweep.
    ///
    /// Runtime evaluation allocates nothing and calls no transcendental
    /// functions. Neither source coefficients nor physical motion are mutated.
    ///
    /// # Errors
    /// Nonfinite inputs/cotangents or unrepresentable physical cotangents.
    pub fn pullback_force_phase(
        &self,
        cosine_force: f64,
        sine_force: f64,
        endpoint_cosine_cotangent: f64,
        endpoint_sine_cotangent: f64,
    ) -> Result<HarmonicForcePhasePullback, ModalAcousticTimeError> {
        if ![endpoint_cosine_cotangent, endpoint_sine_cotangent]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("harmonic source-phase cotangents must be finite"));
        }
        let (cosine, sine) = self.advance_force_phase(cosine_force, sine_force)?;
        let cosine_force = self.phase_cosine * endpoint_cosine_cotangent
            - self.phase_sine * endpoint_sine_cotangent;
        let sine_force = self.phase_sine * endpoint_cosine_cotangent
            + self.phase_cosine * endpoint_sine_cotangent;
        let forcing_angular_frequency = self.duration_s
            * sine.mul_add(endpoint_cosine_cotangent, -cosine * endpoint_sine_cotangent);
        if ![cosine_force, sine_force, forcing_angular_frequency]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "harmonic source-phase pullback must be representable",
            ));
        }
        Ok(HarmonicForcePhasePullback {
            cosine_force,
            sine_force,
            forcing_angular_frequency,
        })
    }
}

impl PreparedHarmonicForceSensitivity {
    pub(in crate::time::harmonic) fn clock_partials(
        &self,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<[f64; 4], ModalAcousticTimeError> {
        let h = self.motion.duration_s;
        let fc = h * cosine_force;
        let fs = h * sine_force;
        let [qc, qs, vc, vs] = self.force_envelope;
        // delta_gains integrates (t/h)*(Dc*cos+Ds*sin). The frequency source is
        // t*(Fs*cos-Fc*sin): one extra h, with the original rotated quadratures.
        let qw = h * (h * qc.mul_add(fs, -qs * fc));
        let vw = h * vc.mul_add(fs, -vs * fc);
        // The phase source has no t factor. Reuse the constant-tone owner,
        // including its separately evaluated near-DC cosine remainder.
        let phase = self.motion.forced_increment(sine_force, -cosine_force)?;
        let result = [
            qw,
            vw,
            phase.displacement_m_sqrt_kg,
            phase.velocity_m_sqrt_kg_per_s,
        ];
        if !result.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "harmonic physical clock partials must be representable",
            ));
        }
        Ok(result)
    }
}
