//! Analytic rectangular-window moments of real, peak-phasor harmonics.
//!
//! For `x(t) = Re[a exp(-i omega t)]`, a detector measures moments over
//! its actual observation interval. Different frequencies retain their beat
//! terms; even a single frequency retains its sum-frequency term when the
//! interval is shorter than a carrier period. Positive-time conventions use
//! conjugated amplitudes. [`HarmonicProduct`] integrates deterministic fields;
//! [`ProperHarmonicProduct`] integrates explicitly authored, jointly proper
//! complex covariance factors. Neither path infers source statistics or gives
//! an error certificate.

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

    /// Cache a real second moment for one explicitly proper complex factor.
    ///
    /// For `E[z] = E[z*z] = 0`, `E[|z|²] = 1`, this integrates
    /// `E[Re(a*z*exp(-i*omega_a*t)) * Re(b*z*exp(-i*omega_b*t))]`.
    /// The result includes the peak-amplitude factor 1/2. Only the difference
    /// frequency survives: properness removes the sum term even for a short
    /// window. Correlated frequencies can still describe a nonstationary field.
    ///
    /// Multiple factors require JOINT moments `E[z_r*conj(z_s)] = delta_rs`
    /// and `E[z_r*z_s] = 0`; their second moments add. No Gaussian distribution,
    /// stationarity, physical source covariance or PSD normalization is inferred.
    /// Deterministic mean fields use [`Self::product_kernel`] separately.
    ///
    /// The real-part covariance identity follows Neeser and Massey (1993),
    /// equations (2)–(5) and Definition 1:
    /// <https://www.isiweb.ee.ethz.ch/papers/arch/fneese-mass-inspec-1993-1.pdf>.
    #[allow(clippy::float_cmp)] // Equal physical frequencies have exactly zero beat.
    pub fn proper_product_kernel(
        self,
        omega_a: f64,
        omega_b: f64,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        if !omega_a.is_finite() || !omega_b.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let half = Dd::from_f64(0.5);
        if omega_a == omega_b {
            return Ok(ProperHarmonicProduct {
                real: half,
                imaginary: Dd::ZERO,
            });
        }
        // Subtract before applying the physical clock. Two independently
        // rounded carrier phases can lose a well-resolved, slow beat.
        let difference = checked(Dd::from_f64(omega_a) - Dd::from_f64(omega_b))?;
        let width = if difference.hi.abs() >= self.duration_s {
            checked((difference * half) * Dd::from_f64(self.duration_s))?
        } else {
            checked(difference * Dd::from_f64(0.5 * self.duration_s))?
        };
        let start = checked(difference * Dd::from_f64(self.start_s))?;
        // Keep start and half-width as separate rotations, just as in the
        // deterministic path: even their DD sum can lose a tiny duration on
        // an enormous clock. Preserve cos-1 as well as sinc-1 near cancellation.
        let (sr, si) = proper_phase(start);
        let (wr, wi) = proper_phase(width);
        let scale = half * sinc(width);
        Ok(ProperHarmonicProduct {
            real: checked(checked(sr * wr - si * wi)? * scale)?,
            imaginary: checked(checked(sr * wi + si * wr)? * scale)?,
        })
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

/// Prepared real second-moment contraction of a proper complex factor.
///
/// Stores `(1/2) * mean(exp(-i*(omega_a-omega_b)*t))` in double-double
/// arithmetic. Cache once per frequency pair and reuse across field components
/// and spatial points with no trigonometry or allocation in the contraction.
/// This is the exposure average of the instantaneous second moment, not the
/// variance of a window-averaged field and not detector uncertainty.
#[derive(Debug, Clone, Copy)]
pub struct ProperHarmonicProduct {
    real: Dd,
    imaginary: Dd,
}

impl ProperHarmonicProduct {
    /// Average proper-ensemble real product, including peak normalization 1/2.
    pub fn mean_real_product(self, a: C64, b: C64) -> Result<f64, WindowError> {
        self.mean_real_product_dd(a, b).map(Dd::to_f64)
    }

    /// Retain low bits through a complete covariance-factor contraction.
    /// Cross products may be signed; only a fully contracted self-moment is
    /// nonnegative. Floating trigonometry supplies no relative guarantee at an
    /// arbitrary null. Equal-frequency fields should be combined first.
    pub fn mean_real_product_dd(self, a: C64, b: C64) -> Result<Dd, WindowError> {
        if [a.re, a.im, b.re, b.im].iter().any(|x| !x.is_finite()) {
            return Err(WindowError::NonFiniteAmplitude);
        }
        let ar = Dd::from_f64(a.re);
        let ai = Dd::from_f64(a.im);
        let br = Dd::from_f64(b.re);
        let bi = Dd::from_f64(b.im);
        let real = checked(ar * br + ai * bi)?;
        let imaginary = checked(ai * br - ar * bi)?;
        checked(checked(real * self.real)? - checked(imaginary * self.imaginary)?)
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

// Exp(-i*x), retaining the small real decrement near ±1 instead of rounding
// cos(x) to a unit value before a destructive-interference sum. The half-angle
// form also retains this decrement near full turns. Trig leaves remain the
// existing deterministic f64 implementation, not a new transcendental owner.
fn proper_phase(x: Dd) -> (Dd, Dd) {
    fn component(x: f64) -> (Dd, Dd) {
        let (s, c) = det::sin_cos(x);
        let cosine = if c > 0.5 {
            let sine_half = Dd::from_f64(det::sin(0.5 * x));
            Dd::ONE - Dd::from_f64(2.0) * sine_half * sine_half
        } else if c < -0.5 {
            let cosine_half = Dd::from_f64(det::cos(0.5 * x));
            -Dd::ONE + Dd::from_f64(2.0) * cosine_half * cosine_half
        } else {
            Dd::from_f64(c)
        };
        (cosine, Dd::from_f64(-s))
    }
    let (hr, hi) = component(x.hi);
    let (lr, li) = component(x.lo);
    (hr * lr - hi * li, hr * li + hi * lr)
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

    #[test]
    fn g0_proper_single_line_is_half_power_for_every_shutter() {
        let a = C64::new(2.0, -3.0);
        for (start, duration, omega) in [
            (0.0, 1e-30, 0.0),
            (-0.3, 0.02, 17.0),
            (1e200, 1.0, f64::MAX),
        ] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            let kernel = window.proper_product_kernel(omega, omega).unwrap();
            close(kernel.mean_real_product(a, a).unwrap(), 6.5, 0.0);
            close(
                kernel
                    .mean_real_product(C64::new(0.0, 2.0), C64::new(0.0, 2.0))
                    .unwrap(),
                2.0,
                0.0,
            );
        }
    }

    #[test]
    fn g1_proper_products_match_four_phase_real_field_quadrature() {
        let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
        let a = C64::new(-0.7, 1.3);
        let b = C64::new(2.1, -0.4);
        for (start, duration, wa, wb) in [
            (-0.13, 0.027, 13.0, 17.0),
            (0.25, 0.3, 200.0, 201.0),
            (-0.2, 0.8, -19.0, 23.0),
            (0.15, 0.5, 0.0, 31.0),
        ] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            // The finite ensemble has E[z]=E[z²]=0 and E[|z|²]=1 exactly.
            // Its oracle integrates REAL signals, without Hermitian products.
            let expected = phases
                .iter()
                .map(|&z| simpson(window, a * z, wa, b * z, wb))
                .sum::<f64>()
                * 0.25;
            let value = window
                .proper_product_kernel(wa, wb)
                .unwrap()
                .mean_real_product(a, b)
                .unwrap();
            close(value, expected, 3e-11);
            close(
                window
                    .proper_product_kernel(wb, wa)
                    .unwrap()
                    .mean_real_product(b, a)
                    .unwrap(),
                value,
                3e-15,
            );
        }
    }

    #[test]
    fn g1_proper_destructive_beat_retains_sinc_and_midpoint_decrements() {
        let delta = 2.0_f64.powi(-20);
        let duration = 2.0_f64.powi(-20);
        let frequencies = [1.0, 1.0 + delta];
        let amplitudes = [C64::ONE, -C64::ONE];
        for (start, divisor) in [(0.0, 6.0), (-duration * 0.5, 24.0)] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            let mut sum = Dd::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    sum = sum
                        + window
                            .proper_product_kernel(frequencies[i], frequencies[j])
                            .unwrap()
                            .mean_real_product_dd(amplitudes[i], amplitudes[j])
                            .unwrap();
                }
            }
            // 1-sinc(delta*T) on [0,T], or 1-sinc(delta*T/2) on [-T/2,T/2].
            // The next series term is negligible at this exactly dyadic scale.
            let expected = (delta * duration).powi(2) / divisor;
            assert!(sum.to_f64() > 0.0);
            close(sum.to_f64(), expected, expected * 3e-14);
        }
    }

    #[test]
    fn g3_proper_kernel_retains_beat_clock_and_ignores_common_carrier() {
        let a = C64::new(0.7, -1.2);
        let b = C64::new(-0.2, 0.8);
        let window = HarmonicWindow::new(9_007_199_254_740_992.0, 1.0).unwrap();
        let carrier = 2.0_f64.powi(40);
        let carrier_value = window
            .proper_product_kernel(carrier + 1.0, carrier)
            .unwrap()
            .mean_real_product_dd(a, b)
            .unwrap();
        assert_eq!(
            carrier_value,
            window
                .proper_product_kernel(1.0, 0.0)
                .unwrap()
                .mean_real_product_dd(a, b)
                .unwrap()
        );
        let original = HarmonicWindow::new(-0.13, 0.027).unwrap();
        let shift = 0.4_f64;
        let rotate_clock = |z, omega: f64| {
            let phase = omega * shift;
            z * C64::new(phase.cos(), -phase.sin())
        };
        let expected = original
            .proper_product_kernel(13.0, 17.0)
            .unwrap()
            .mean_real_product(a, b)
            .unwrap();
        let shifted = HarmonicWindow::new(original.start_s() - shift, original.duration_s())
            .unwrap()
            .proper_product_kernel(13.0, 17.0)
            .unwrap()
            .mean_real_product(rotate_clock(a, 13.0), rotate_clock(b, 17.0))
            .unwrap();
        close(shifted, expected, 2e-14);
    }

    #[test]
    fn g0_proper_kernel_refuses_invalid_or_unrepresentable_moments() {
        let window = HarmonicWindow::new(0.0, 1.0).unwrap();
        assert!(matches!(
            window.proper_product_kernel(f64::NAN, 0.0),
            Err(WindowError::InvalidFrequency)
        ));
        let kernel = window.proper_product_kernel(1.0, 1.0).unwrap();
        assert_eq!(
            kernel.mean_real_product(C64::new(f64::NAN, 0.0), C64::ONE),
            Err(WindowError::NonFiniteAmplitude)
        );
        assert_eq!(
            kernel.mean_real_product(C64::from_re(f64::MAX), C64::from_re(f64::MAX)),
            Err(WindowError::Unrepresentable)
        );
        assert!(matches!(
            HarmonicWindow::new(1e200, 1.0)
                .unwrap()
                .proper_product_kernel(1e200, 0.0),
            Err(WindowError::Unrepresentable)
        ));
    }
}
