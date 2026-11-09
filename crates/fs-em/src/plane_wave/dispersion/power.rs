//! Normal Poynting-flux fractions and their original coherent frequency slope.
//!
//! Both fields use the same incident power port: E/E_inc and H/(Y_inc E_inc).
//! Consequently T=Re(E conj(H)), including dispersive exit AND incident
//! admittances. A transmitted electric-amplitude slope alone omits that work.

use fs_math::{c64::C64, dd::Dd, det};

/// A native directed power fraction and its analytic angular-frequency slope.
/// These are monochromatic normal energy-flux ratios, not field intensities,
/// material heat, finite-band pulse energy or an optical phase measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerFractionDispersion {
    fraction: f64,
    omega: f64,
    scaled_derivative: f64,
}

impl PowerFractionDispersion {
    /// Original native power fraction, without clipping or renormalization.
    pub const fn fraction(self) -> f64 {
        self.fraction
    }
    /// omega d(fraction)/domega; it can have either sign.
    pub const fn scaled_derivative(self) -> f64 {
        self.scaled_derivative
    }
    /// d(fraction)/domega, per rad/s.
    pub fn derivative_s(self) -> f64 {
        self.scaled_derivative / self.omega
    }

    /// The primal power comes from its native positive-flux calculation. The
    /// tangent retains BOTH coherent fields and their original port slopes.
    pub(in crate::plane_wave) fn from_fields(
        fraction: f64,
        omega: f64,
        electric: C64,
        magnetic: C64,
        electric_slope: C64,
        magnetic_slope: C64,
    ) -> Option<Self> {
        if !fraction.is_finite()
            || fraction < 0.0
            || !omega.is_finite()
            || omega <= 0.0
            || [electric, magnetic, electric_slope, magnetic_slope]
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return None;
        }
        let slope = Dd::from_f64(electric_slope.re) * Dd::from_f64(magnetic.re)
            + Dd::from_f64(electric_slope.im) * Dd::from_f64(magnetic.im)
            + Dd::from_f64(electric.re) * Dd::from_f64(magnetic_slope.re)
            + Dd::from_f64(electric.im) * Dd::from_f64(magnetic_slope.im);
        let scaled_derivative = slope.to_f64();
        if !scaled_derivative.is_finite() || !(scaled_derivative / omega).is_finite() {
            return None;
        }
        Some(Self {
            fraction,
            omega,
            scaled_derivative,
        })
    }

    /// A phase-free search observation A=sqrt(P), D A=D P/(2 sqrt(P)).
    /// Its squared amplitude is the ACTUAL native power fraction. Zero power
    /// has no admitted square-root tangent here; no phase or slope is invented.
    pub(in crate::plane_wave) fn square_root_observation(self) -> Option<(C64, C64)> {
        if self.fraction == 0.0 {
            return None;
        }
        let amplitude = det::sqrt(self.fraction);
        let slope = (self.scaled_derivative / amplitude) * 0.5;
        if !slope.is_finite() || !(slope / self.omega).is_finite() {
            return None;
        }
        Some((C64::from_re(amplitude), C64::from_re(slope)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 3e-15 * b.abs().max(1.0), "{a} != {b}");
    }

    #[test]
    fn dispersive_admittance_and_coherent_phase_enter_the_original_flux_slope() {
        let e = C64::new(0.5, -0.25);
        let de = C64::new(-0.125, 0.375);
        let y = C64::new(2.0, -3.0);
        let dy = C64::new(0.25, 4.0);
        let h = y * e;
        let dh = dy * e + y * de;
        let p = y.re * e.norm_sq();
        let expected = dy.re * e.norm_sq() + 2.0 * y.re * (de * e.conj()).re;
        let jet = PowerFractionDispersion::from_fields(p, 4.0, e, h, de, dh).unwrap();
        close(jet.fraction(), p);
        close(jet.scaled_derivative(), expected);
        close(jet.derivative_s(), expected / 4.0);
        let (amplitude, slope) = jet.square_root_observation().unwrap();
        close(amplitude.norm_sq(), p);
        close(2.0 * (amplitude.conj() * slope).re, expected);
        // The field-amplitude objective misses the exit-admittance direction.
        assert!((jet.scaled_derivative() - 2.0 * (de * e.conj()).re).abs() > 0.1);
    }

    #[test]
    fn common_phase_motion_cancels_without_losing_small_real_flux_change() {
        let e = C64::new(0.5, 0.25);
        let h = C64::new(0.125, -0.5);
        let phase = C64::new(0.0, 1_099_511_627_776.0);
        let de = phase * e + C64::new(0.125, 0.0);
        let dh = phase * h;
        // Negative signed reference flux is not a nonnegative power fraction;
        // select the positive-flux conjugate port for this fixture instead.
        let h = -h;
        let dh = -dh;
        let expected = 0.125 * h.re;
        let power = (e * h.conj()).re;
        let result = PowerFractionDispersion::from_fields(power, 2.0, e, h, de, dh).unwrap();
        close(result.scaled_derivative(), expected);
        let multiplier = 1_048_576.0;
        let scaled = PowerFractionDispersion::from_fields(
            power,
            2.0,
            e.scale(multiplier),
            h.scale(1.0 / multiplier),
            de.scale(multiplier),
            dh.scale(1.0 / multiplier),
        )
        .unwrap();
        assert_eq!(scaled, result);
    }

    #[test]
    fn zero_power_has_no_sqrt_search_jet_and_invalid_data_refuse() {
        let zero = PowerFractionDispersion::from_fields(
            0.0,
            1.0,
            C64::ONE,
            C64::new(0.0, 1.0),
            C64::ZERO,
            C64::ZERO,
        )
        .unwrap();
        assert_eq!(zero.scaled_derivative(), 0.0);
        assert!(zero.square_root_observation().is_none());
        for (power, omega) in [(-1.0, 1.0), (1.0, 0.0), (f64::NAN, 1.0)] {
            assert!(
                PowerFractionDispersion::from_fields(
                    power,
                    omega,
                    C64::ONE,
                    C64::ONE,
                    C64::ZERO,
                    C64::ZERO
                )
                .is_none()
            );
        }
    }
}
