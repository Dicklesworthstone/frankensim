//! Analytic rectangular-window moments of real, peak-phasor harmonics.
//!
//! For `x(t) = Re[a exp(-i omega t)]`, a detector measures moments over
//! its actual observation interval. Different frequencies retain their beat
//! terms; even a single frequency retains its sum-frequency term when the
//! interval is shorter than a carrier period. Positive-time conventions use
//! conjugated amplitudes. This is integration of authored harmonic fields,
//! not a source-coherence model or an error certificate.

use fs_math::{c64::C64, dd::Dd, det};

/// A rejected observation window or unrepresentable moment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowError {
    /// The start/end must be finite and the duration finite and positive.
    InvalidWindow,
    /// Angular frequencies must be finite (signed frequencies and DC work).
    InvalidFrequency,
    /// Both components of each amplitude must be finite.
    NonFiniteAmplitude,
    /// A phase product, amplitude rotation or result exceeds floating range.
    Unrepresentable,
}

impl core::fmt::Display for WindowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidWindow => "harmonic window needs a finite start/end and positive duration",
            Self::InvalidFrequency => "harmonic angular frequency must be finite",
            Self::NonFiniteAmplitude => "harmonic peak amplitude must be finite",
            Self::Unrepresentable => "harmonic window moment is not representable",
        })
    }
}

impl std::error::Error for WindowError {}

/// A uniform observation on `[start_s, start_s + duration_s]` in seconds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicWindow {
    start_s: f64,
    duration_s: f64,
}

impl HarmonicWindow {
    /// Admit an interval without rounding a short duration into its start.
    pub fn new(start_s: f64, duration_s: f64) -> Result<Self, WindowError> {
        if !start_s.is_finite()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !(start_s + duration_s).is_finite()
        {
            return Err(WindowError::InvalidWindow);
        }
        Ok(Self {
            start_s,
            duration_s,
        })
    }

    /// Beginning of the observation, in seconds on the source clock.
    #[must_use]
    pub const fn start_s(self) -> f64 {
        self.start_s
    }

    /// Duration in seconds, retained independently of the start timestamp.
    #[must_use]
    pub const fn duration_s(self) -> f64 {
        self.duration_s
    }

    /// Average of `exp(-i omega t)`, including the exact DC limit.
    pub fn mean_phase(self, omega: f64) -> Result<C64, WindowError> {
        let half_width = self.half_width(omega)?;
        let scale = sinc(half_width).to_f64();
        Ok(self.rotation(omega, half_width)?.scale(scale))
    }

    /// Cache the time integration for a frequency pair. No allocation; reuse
    /// it for every field component or spatial point at those frequencies.
    pub fn product_kernel(
        self,
        omega_a: f64,
        omega_b: f64,
    ) -> Result<HarmonicProduct, WindowError> {
        let a = self.half_width(omega_a)?;
        let b = self.half_width(omega_b)?;
        let difference = checked(a - b)?;
        let sum = checked(a + b)?;
        let half = Dd::from_f64(0.5);
        // About the midpoint, cos*sin is odd and integrates to zero.
        // C = (sinc(a-b)+sinc(a+b))/2, S = (sinc(a-b)-sinc(a+b))/2.
        // Subtract sinc-1 directly for short shutters: subtracting two rounded
        // numbers near one would erase the positive sin² contribution.
        let (cos_cos, sin_sin) = if difference.hi.abs().max(sum.hi.abs()) <= 0.5 {
            let d = sinc_minus_one(difference);
            let s = sinc_minus_one(sum);
            (Dd::ONE + half * (d + s), half * (d - s))
        } else {
            let d = sinc(difference);
            let s = sinc(sum);
            let (small, large) = if a.hi.abs() <= b.hi.abs() {
                (a, b)
            } else {
                (b, a)
            };
            let ratio = small / large;
            let ss = if ratio.hi.abs() < 1e-4 {
                // Stable divided difference when one frequency tends to DC.
                let ca = phase(small).re;
                let cb = phase(large).re;
                ratio * (sinc(large) * Dd::from_f64(ca) - sinc(small) * Dd::from_f64(cb))
                    / (Dd::ONE - ratio * ratio)
            } else {
                half * (d - s)
            };
            (half * (d + s), ss)
        };
        Ok(HarmonicProduct {
            rotation_a: self.rotation(omega_a, a)?,
            rotation_b: self.rotation(omega_b, b)?,
            cos_cos: checked(cos_cos)?,
            sin_sin: checked(sin_sin)?,
        })
    }

    /// Average of the actual real product, retaining sum and difference
    /// frequencies. Amplitudes are peak values in the negative-time convention.
    pub fn mean_real_product(
        self,
        a: C64,
        omega_a: f64,
        b: C64,
        omega_b: f64,
    ) -> Result<f64, WindowError> {
        self.product_kernel(omega_a, omega_b)?
            .mean_real_product(a, b)
    }

    fn half_width(self, omega: f64) -> Result<Dd, WindowError> {
        if !omega.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        // Halve the larger operand, so neither a tiny duration nor a huge
        // frequency loses a representable product through scaling order.
        let (a, b) = if omega.abs() >= self.duration_s {
            (0.5 * omega, self.duration_s)
        } else {
            (omega, 0.5 * self.duration_s)
        };
        checked(Dd::from_f64(a) * Dd::from_f64(b))
    }

    fn rotation(self, omega: f64, half_width: Dd) -> Result<C64, WindowError> {
        let start = checked(Dd::from_f64(omega) * Dd::from_f64(self.start_s))?;
        // Do not form start + duration/2: a short window on a large clock
        // still has a resolved duration. TwoProd also retains product phase
        // bits that a single rounded omega*time multiplication would discard.
        Ok(phase(start) * phase(half_width))
    }
}

/// Prepared constant-work real-product moment for one frequency pair/window.
///
/// Rotations and centered Gram weights are retained, so repeated contractions
/// need no trigonometry or allocation. Scalar products use double-double
/// arithmetic; this is not an interval or a relative guarantee near nulls.
#[derive(Debug, Clone, Copy)]
pub struct HarmonicProduct {
    rotation_a: C64,
    rotation_b: C64,
    cos_cos: Dd,
    sin_sin: Dd,
}

impl HarmonicProduct {
    /// Average `Re[a exp(-i omega_a t)] * Re[b exp(-i omega_b t)]`.
    pub fn mean_real_product(self, a: C64, b: C64) -> Result<f64, WindowError> {
        self.mean_real_product_dd(a, b).map(Dd::to_f64)
    }

    /// Retain the double-double product until a whole field's interference
    /// sum is complete. Rounding each large pair to f64 before cancellation
    /// can discard a small beat signal even with an accurate outer sum.
    /// Cached trigonometry remains f64; no relative error bound near an
    /// arbitrary null is implied by this extended-precision accumulator.
    pub fn mean_real_product_dd(self, a: C64, b: C64) -> Result<Dd, WindowError> {
        let (ar, ai) = rotate(a, self.rotation_a)?;
        let (br, bi) = rotate(b, self.rotation_b)?;
        let cosine = checked(checked(self.cos_cos * ar)? * br)?;
        let sine = checked(checked(self.sin_sin * ai)? * bi)?;
        checked(cosine + sine)
    }
}

fn checked(x: Dd) -> Result<Dd, WindowError> {
    if x.hi.is_finite() && x.lo.is_finite() {
        Ok(x)
    } else {
        Err(WindowError::Unrepresentable)
    }
}

fn rotate(a: C64, r: C64) -> Result<(Dd, Dd), WindowError> {
    if !a.re.is_finite() || !a.im.is_finite() {
        return Err(WindowError::NonFiniteAmplitude);
    }
    let ar = Dd::from_f64(a.re);
    let ai = Dd::from_f64(a.im);
    let rr = Dd::from_f64(r.re);
    let ri = Dd::from_f64(r.im);
    Ok((checked(ar * rr - ai * ri)?, checked(ar * ri + ai * rr)?))
}

fn phase(x: Dd) -> C64 {
    let (sh, ch) = det::sin_cos(x.hi);
    let (sl, cl) = det::sin_cos(x.lo);
    C64::new(ch.mul_add(cl, -sh * sl), -sh.mul_add(cl, ch * sl))
}

fn sinc(x: Dd) -> Dd {
    if x.hi.abs() <= 0.5 {
        Dd::ONE + sinc_minus_one(x)
    } else {
        Dd::from_f64(-phase(x).im) / x
    }
}

fn sinc_minus_one(x: Dd) -> Dd {
    let minus_square = -(x * x);
    let mut term = minus_square / Dd::from_f64(6.0);
    let mut result = term;
    // At |x| <= 1/2 the omitted term is below 2e-36. In particular the
    // leading x²/6 survives even when rounded sinc(x) equals exactly one.
    for n in 2..=12 {
        term = term * minus_square / Dd::from_f64(f64::from((2 * n) * (2 * n + 1)));
        result = result + term;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, TAU};

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    // G1 independent integration of REAL signals, with platform trigonometry
    // confined to this test oracle. It never calls a phasor-product formula.
    fn simpson(window: HarmonicWindow, a: C64, wa: f64, b: C64, wb: f64) -> f64 {
        let n = 20_000_usize;
        let dt = window.duration_s() / n as f64;
        let mut total = 0.0;
        for i in 0..=n {
            let t = window.start_s() + i as f64 * dt;
            let x = a.re * (wa * t).cos() + a.im * (wa * t).sin();
            let y = b.re * (wb * t).cos() + b.im * (wb * t).sin();
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            total += weight * x * y;
        }
        total / (3 * n) as f64
    }

    #[test]
    fn g0_dc_and_complete_cycle_peak_normalization() {
        let window = HarmonicWindow::new(0.0, TAU).unwrap();
        assert_eq!(window.mean_phase(0.0).unwrap(), C64::ONE);
        close(
            window
                .mean_real_product(C64::new(2.0, 7.0), 0.0, C64::new(3.0, -4.0), 0.0)
                .unwrap(),
            6.0,
            0.0,
        );
        let a = C64::new(2.0, -3.0);
        let b = C64::new(4.0, 5.0);
        close(
            window.mean_real_product(a, 1.0, b, 1.0).unwrap(),
            -3.5,
            2e-15,
        );
        close(
            window.mean_real_product(a, 1.0, b, 2.0).unwrap(),
            0.0,
            2e-15,
        );
    }

    #[test]
    fn g1_short_shutter_retains_sum_frequency_and_sine_energy() {
        let window = HarmonicWindow::new(0.0, PI / 4.0).unwrap();
        let cosine = C64::ONE;
        let sine = C64::new(0.0, 1.0);
        close(
            window.mean_real_product(cosine, 1.0, cosine, 1.0).unwrap(),
            0.5 + 1.0 / PI,
            3e-16,
        );
        close(
            window.mean_real_product(cosine, 1.0, sine, 1.0).unwrap(),
            1.0 / PI,
            3e-16,
        );
        // Standard 0.5*(1-sinc(2*omega*h)) rounds to zero here. The
        // centered Gram weights must retain the physical omega²*T²/12.
        for duration in [1e-4, 1e-8, 1e-12, 1e-50] {
            let w = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
            let expected = duration * duration / 12.0;
            let actual = w.mean_real_product(sine, 1.0, sine, 1.0).unwrap();
            close(actual, expected, 1e-9 * expected);
            assert!(actual > 0.0);
        }
    }

    #[test]
    fn g1_arbitrary_phases_beats_and_signed_frequencies_match_real_quadrature() {
        let a = C64::new(-0.7, 1.3);
        let b = C64::new(2.1, -0.4);
        for (start, duration, wa, wb) in [
            (-0.13, 0.027, 13.0, 17.0),
            (0.25, 0.3, 200.0, 201.0),
            (-0.2, 0.8, -19.0, 23.0),
            (0.15, 0.5, 0.0, 31.0),
        ] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            let got = window.mean_real_product(a, wa, b, wb).unwrap();
            close(got, simpson(window, a, wa, b, wb), 3e-11);
            close(window.mean_real_product(b, wb, a, wa).unwrap(), got, 2e-15);
        }
    }

    #[test]
    fn g1_small_frequency_divided_difference_retains_cross_energy() {
        let window = HarmonicWindow::new(-1.0, 2.0).unwrap();
        let a = C64::new(0.0, 1.0);
        let small = 1e-16;
        let large = 3.0_f64;
        // sin(epsilon*t) = epsilon*t + O(epsilon³); the integral of
        // t*sin(b*t) on [-1,1], divided by 2, is (sin(b)-b*cos(b))/b².
        let expected = small * (large.sin() - large * large.cos()) / (large * large);
        let got = window.mean_real_product(a, small, a, large).unwrap();
        close(got, expected, expected.abs() * 1e-14);
    }

    #[test]
    fn g3_short_duration_survives_a_large_clock_origin() {
        let start = 9_007_199_254_740_992.0;
        let window = HarmonicWindow::new(start, 1.0).unwrap();
        assert_eq!(start + 0.5, start); // A rounded midpoint loses the window.
        let c = start.cos();
        let s = start.sin();
        let expected_re = (c * 0.5_f64.cos() - s * 0.5_f64.sin()) * (0.5_f64.sin() / 0.5);
        let expected_im = -(s * 0.5_f64.cos() + c * 0.5_f64.sin()) * (0.5_f64.sin() / 0.5);
        let mean = window.mean_phase(1.0).unwrap();
        close(mean.re, expected_re, 3e-16);
        close(mean.im, expected_im, 3e-16);
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_are_errors() {
        for (start, duration) in [
            (f64::NAN, 1.0),
            (0.0, 0.0),
            (0.0, -1.0),
            (0.0, f64::INFINITY),
            (f64::MAX, f64::MAX),
        ] {
            assert_eq!(
                HarmonicWindow::new(start, duration),
                Err(WindowError::InvalidWindow)
            );
        }
        let window = HarmonicWindow::new(0.0, 1.0).unwrap();
        assert_eq!(
            window.mean_phase(f64::NAN),
            Err(WindowError::InvalidFrequency)
        );
        assert_eq!(
            window.mean_real_product(C64::new(f64::NAN, 0.0), 1.0, C64::ONE, 1.0),
            Err(WindowError::NonFiniteAmplitude)
        );
        assert_eq!(
            window.mean_real_product(C64::from_re(f64::MAX), 1.0, C64::from_re(f64::MAX), 1.0),
            Err(WindowError::Unrepresentable)
        );
        assert_eq!(
            HarmonicWindow::new(1e200, 1.0).unwrap().mean_phase(1e200),
            Err(WindowError::Unrepresentable)
        );
    }
}
