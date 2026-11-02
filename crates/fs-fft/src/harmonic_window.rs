//! Analytic normalized-window moments of real, peak-phasor harmonics.
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

mod spectrum;
pub use spectrum::{
    AngularHarmonicProductInput, HarmonicProductInput, HarmonicProductLaw,
    HarmonicProductSpectrumLimits, HarmonicProductTone, MAX_PRODUCT_SPECTRUM_TERMS,
    PreparedHarmonicProductSpectrum, ProductSpectrumError,
};
pub use spectrum::{
    AngularHarmonicSpectrumTerm, HarmonicCausalIntegral, HarmonicSpectrumTermLimits,
};

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

/// Nonnegative, unit-integral temporal sensitivity of a detector.
///
/// These profiles weight the instantaneous quantity being observed. They do
/// not taper each field amplitude before multiplication: that would square the
/// sensitivity and describe a different detector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HarmonicWindowProfile {
    /// Constant sensitivity 1/T over the physical interval.
    Uniform,
    /// Linear opening and closing, zero at either endpoint and 2/T at the
    /// midpoint. For u = (t-start)/T in [0,1], w(t) = 2*(1-|2*u-1|)/T.
    Triangular,
}

/// A normalized observation on `[start_s, start_s + duration_s]` in seconds.
///
/// Every returned mean is integral(w(t) * quantity(t) dt), where w has unit
/// integral. Multiplication by duration gives a sensitivity-weighted exposure.
/// Only the uniform profile gives the unweighted interval integral or makes
/// a mean storage rate equal to the endpoint storage difference divided by T.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicWindow {
    start_s: f64,
    duration_s: f64,
    profile: HarmonicWindowProfile,
}

impl HarmonicWindow {
    /// Admit a uniform interval without rounding a short duration into its start.
    pub fn new(start_s: f64, duration_s: f64) -> Result<Self, WindowError> {
        Self::with_profile(start_s, duration_s, HarmonicWindowProfile::Uniform)
    }

    /// Admit an explicitly shaped, normalized detector sensitivity.
    ///
    /// Profiles integrate products of the physical fields with one temporal
    /// weight. They neither rescale peak amplitudes nor infer source statistics.
    pub fn with_profile(
        start_s: f64,
        duration_s: f64,
        profile: HarmonicWindowProfile,
    ) -> Result<Self, WindowError> {
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
            profile,
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

    /// Original, explicitly selected detector sensitivity.
    #[must_use]
    pub const fn profile(self) -> HarmonicWindowProfile {
        self.profile
    }

    /// Weighted average of `exp(-i omega t)`, including the exact DC limit.
    pub fn mean_phase(self, omega: f64) -> Result<C64, WindowError> {
        if !omega.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        self.mean_phase_extended(Dd::from_f64(omega))
    }

    /// Mean phase for a frequency authored in cycles/second. Convert to rad/s
    /// in double-double precision before applying the physical source clock.
    /// This avoids an independently rounded carrier shifting a resolved beat.
    pub fn mean_phase_hz(self, frequency_hz: f64) -> Result<C64, WindowError> {
        self.mean_phase_extended(angular_frequency(frequency_hz)?)
    }

    fn mean_phase_extended(self, omega: Dd) -> Result<C64, WindowError> {
        let half_width = self.half_width_extended(omega)?;
        let scale = self.centered_transform(half_width);
        let [real, imaginary] = self.rotation(omega, half_width)?;
        Ok(C64::new(
            checked(real * scale)?.to_f64(),
            checked(imaginary * scale)?.to_f64(),
        ))
    }

    /// Cache the time integration for a frequency pair. No allocation; reuse
    /// it for every field component or spatial point at those frequencies.
    pub fn product_kernel(
        self,
        omega_a: f64,
        omega_b: f64,
    ) -> Result<HarmonicProduct, WindowError> {
        if !omega_a.is_finite() || !omega_b.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        self.product_kernel_extended(Dd::from_f64(omega_a), Dd::from_f64(omega_b))
    }

    /// Deterministic real-product kernel for frequencies authored in Hz.
    ///
    /// Retains the Hz-to-radian conversion through both carrier rotations and
    /// the centered sum/difference Gram weights. Separately rounded carriers
    /// can otherwise turn a slow beat into a different physical frequency.
    pub fn product_kernel_hz(
        self,
        frequency_a_hz: f64,
        frequency_b_hz: f64,
    ) -> Result<HarmonicProduct, WindowError> {
        self.product_kernel_extended(
            angular_frequency(frequency_a_hz)?,
            angular_frequency(frequency_b_hz)?,
        )
    }

    fn product_kernel_extended(
        self,
        omega_a: Dd,
        omega_b: Dd,
    ) -> Result<HarmonicProduct, WindowError> {
        let a = self.half_width_extended(omega_a)?;
        let b = self.half_width_extended(omega_b)?;
        let (cos_cos, sin_sin) = match self.profile {
            HarmonicWindowProfile::Uniform => rectangular_gram(a, b)?,
            HarmonicWindowProfile::Triangular => {
                // A triangular clock displacement is the sum of two uniform
                // displacements with half the support. If their Gram entries
                // are C,S, independence and symmetry give C²+S² and 2CS.
                // This avoids subtracting two almost equal sinc² transforms
                // for a short shutter or vastly disparate frequencies.
                let half = Dd::from_f64(0.5);
                let (c, s) = rectangular_gram(a * half, b * half)?;
                (checked(c * c + s * s)?, checked(Dd::from_f64(2.0) * c * s)?)
            }
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
    pub fn proper_product_kernel(
        self,
        omega_a: f64,
        omega_b: f64,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        if !omega_a.is_finite() || !omega_b.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        // Subtract before applying the physical clock. Two independently
        // rounded carrier phases can lose a well-resolved, slow beat.
        self.proper_difference_kernel(checked(Dd::from_f64(omega_a) - Dd::from_f64(omega_b))?)
    }

    /// Proper second-moment kernel for frequencies authored in cycles/second.
    ///
    /// The difference is formed in double-double BEFORE conversion to rad/s.
    /// Separately rounding `TAU*f_a` and `TAU*f_b` can corrupt a resolved beat
    /// between close carriers; rounding the difference before clock evaluation
    /// can also discard the phase of a much smaller component. This method
    /// retains both operations through the same difference-only kernel.
    pub fn proper_product_kernel_hz(
        self,
        frequency_a_hz: f64,
        frequency_b_hz: f64,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        if !frequency_a_hz.is_finite() || !frequency_b_hz.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let difference = checked(Dd::from_f64(frequency_a_hz) - Dd::from_f64(frequency_b_hz))?;
        self.proper_difference_kernel(checked(difference * Dd::from_f64(core::f64::consts::TAU))?)
    }

    #[allow(clippy::float_cmp)] // Exact zero is the normalized DD identity.
    fn proper_difference_kernel(
        self,
        difference: Dd,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        let half = Dd::from_f64(0.5);
        if difference.hi == 0.0 && difference.lo == 0.0 {
            return Ok(ProperHarmonicProduct {
                real: half,
                imaginary: Dd::ZERO,
            });
        }
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
        let scale = half * self.centered_transform(width);
        Ok(ProperHarmonicProduct {
            real: checked(checked(sr * wr - si * wi)? * scale)?,
            imaginary: checked(checked(sr * wi + si * wr)? * scale)?,
        })
    }

    fn centered_transform(self, half_width: Dd) -> Dd {
        match self.profile {
            HarmonicWindowProfile::Uniform => sinc(half_width),
            HarmonicWindowProfile::Triangular => {
                let half_sinc = sinc(half_width * Dd::from_f64(0.5));
                half_sinc * half_sinc
            }
        }
    }

    fn half_width_extended(self, omega: Dd) -> Result<Dd, WindowError> {
        // Halve the larger operand to preserve representable tiny products.
        if omega.hi.abs() >= self.duration_s {
            checked((omega * Dd::from_f64(0.5)) * Dd::from_f64(self.duration_s))
        } else {
            checked(omega * Dd::from_f64(0.5 * self.duration_s))
        }
    }

    fn rotation(self, omega: Dd, half_width: Dd) -> Result<[Dd; 2], WindowError> {
        let start = checked(omega * Dd::from_f64(self.start_s))?;
        // Preserve a short duration on a large clock and small cosine
        // decrements near cancellation through the full quadratic contraction.
        let (sr, si) = proper_phase(start);
        let (wr, wi) = proper_phase(half_width);
        Ok([checked(sr * wr - si * wi)?, checked(sr * wi + si * wr)?])
    }
}

// Centered cosine/cosine and sine/sine means under a uniform clock.
fn rectangular_gram(a: Dd, b: Dd) -> Result<(Dd, Dd), WindowError> {
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
    Ok((checked(cos_cos)?, checked(sin_sin)?))
}

/// Prepared constant-work real-product moment for one frequency pair/window.
///
/// Rotations and centered Gram weights are retained, so repeated contractions
/// need no trigonometry or allocation. Scalar products use double-double
/// arithmetic; this is not an interval or a relative guarantee near nulls.
#[derive(Debug, Clone, Copy)]
pub struct HarmonicProduct {
    rotation_a: [Dd; 2],
    rotation_b: [Dd; 2],
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
        self.mean_real_product_parts_dd(
            [Dd::from_f64(a.re), Dd::from_f64(a.im)],
            [Dd::from_f64(b.re), Dd::from_f64(b.im)],
        )
    }

    /// Contract retained real/imaginary double-double field components.
    ///
    /// Equal-frequency field sums can carry a small low component that only
    /// becomes visible after a later cancellation. Keeping it here also avoids
    /// expanding one product into four separate C64 contractions.
    pub fn mean_real_product_parts_dd(self, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, WindowError> {
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
        self.mean_real_product_parts_dd(
            [Dd::from_f64(a.re), Dd::from_f64(a.im)],
            [Dd::from_f64(b.re), Dd::from_f64(b.im)],
        )
    }

    /// Contract retained real/imaginary double-double field components.
    pub fn mean_real_product_parts_dd(self, a: [Dd; 2], b: [Dd; 2]) -> Result<Dd, WindowError> {
        for value in a.into_iter().chain(b) {
            check_amplitude(value)?;
        }
        let [ar, ai] = a;
        let [br, bi] = b;
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

fn angular_frequency(frequency_hz: f64) -> Result<Dd, WindowError> {
    if !frequency_hz.is_finite() {
        return Err(WindowError::InvalidFrequency);
    }
    checked(Dd::from_f64(frequency_hz) * Dd::from_f64(core::f64::consts::TAU))
}

fn check_amplitude(value: Dd) -> Result<(), WindowError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(())
    } else {
        Err(WindowError::NonFiniteAmplitude)
    }
}

fn rotate(a: [Dd; 2], r: [Dd; 2]) -> Result<(Dd, Dd), WindowError> {
    check_amplitude(a[0])?;
    check_amplitude(a[1])?;
    let [ar, ai] = a;
    let [rr, ri] = r;
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
            let sensitivity = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => {
                    let u = i as f64 / n as f64;
                    2.0 * (1.0 - (2.0 * u - 1.0).abs())
                }
            };
            total += weight * sensitivity * x * y;
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

    #[test]
    fn g3_proper_hz_conversion_preserves_authored_frequency_differences() {
        let carrier = 2.0_f64.powi(52);
        // Converting the two carriers first yields an 8 rad/s difference,
        // despite an exactly represented physical difference of one hertz.
        assert_eq!(TAU * (carrier + 1.0) - TAU * carrier, 8.0);
        let window = HarmonicWindow::new(0.17, 0.63).unwrap();
        let a = C64::new(0.7, -1.2);
        let b = C64::new(-0.2, 0.8);
        let large = window
            .proper_product_kernel_hz(carrier + 1.0, carrier)
            .unwrap()
            .mean_real_product_dd(a, b)
            .unwrap();
        let small = window
            .proper_product_kernel_hz(1.0, 0.0)
            .unwrap()
            .mean_real_product_dd(a, b)
            .unwrap();
        assert_eq!(large, small);
        close(
            large.to_f64(),
            window
                .proper_product_kernel(TAU, 0.0)
                .unwrap()
                .mean_real_product(a, b)
                .unwrap(),
            0.0,
        );
        assert!(matches!(
            window.proper_product_kernel_hz(f64::INFINITY, 0.0),
            Err(WindowError::InvalidFrequency)
        ));
    }

    #[test]
    fn g1_hz_kernels_match_real_signals_and_keep_a_resolved_high_carrier_beat() {
        let window = HarmonicWindow::new(-0.23, 0.41).unwrap();
        let a = C64::new(0.7, -0.4);
        let b = C64::new(-0.3, 1.2);
        for (fa, fb) in [(3.25, 4.125), (-1.5, 0.0), (0.0, 0.0)] {
            close(
                window
                    .product_kernel_hz(fa, fb)
                    .unwrap()
                    .mean_real_product(a, b)
                    .unwrap(),
                simpson(window, a, TAU * fa, b, TAU * fb),
                2e-12,
            );
        }
        // The sum cos(a)cos(b)+sin(a)sin(b) removes the common carrier
        // exactly. This low-frequency oracle never rotates either carrier.
        let fa = 2.0_f64.powi(48);
        let fb = fa + 0.125;
        let shutter = HarmonicWindow::new(-0.5, 1.0).unwrap();
        let kernel = shutter.product_kernel_hz(fa, fb).unwrap();
        let got = kernel.mean_real_product_dd(C64::ONE, C64::ONE).unwrap()
            + kernel
                .mean_real_product_dd(C64::new(0.0, 1.0), C64::new(0.0, 1.0))
                .unwrap();
        let half_beat = 0.5 * TAU * (fa - fb);
        let expected = half_beat.sin() / half_beat;
        close(got.to_f64(), expected, 3e-14);
        let rounded = shutter.product_kernel(TAU * fa, TAU * fb).unwrap();
        let wrong = rounded.mean_real_product(C64::ONE, C64::ONE).unwrap()
            + rounded
                .mean_real_product(C64::new(0.0, 1.0), C64::new(0.0, 1.0))
                .unwrap();
        assert!((wrong - expected).abs() > 1e-3);
    }

    #[test]
    fn g1_hz_mean_phase_retains_a_fractional_carrier_on_its_clock() {
        let base = 2.0_f64.powi(48);
        let fraction = 0.125;
        let start = 2.0_f64.powi(-20);
        let duration = 2.0_f64.powi(-30);
        let window = HarmonicWindow::new(start, duration).unwrap();
        // Products by powers of two are exact; evaluate the large and small
        // physical frequency contributions independently using platform trig.
        let rotate = |x: f64| C64::new(x.cos(), -x.sin());
        let base_start = TAU * (base * start);
        let base_half = TAU * (base * (0.5 * duration));
        let small_start = TAU * (fraction * start);
        let small_half = TAU * (fraction * (0.5 * duration));
        let expected_rotation =
            rotate(base_start) * rotate(small_start) * rotate(base_half) * rotate(small_half);
        let half = Dd::from_f64(base_half) + Dd::from_f64(small_half);
        let expected_sine = (rotate(base_half) * rotate(small_half)).im;
        let expected = expected_rotation.scale(-expected_sine / half.to_f64());
        let actual = window.mean_phase_hz(base + fraction).unwrap();
        close(actual.re, expected.re, 2e-20);
        close(actual.im, expected.im, 2e-20);
    }

    #[test]
    fn g0_extended_amplitudes_keep_grouped_low_components() {
        let big = Dd::from_f64(2.0_f64.powi(54));
        let weak = Dd::from_f64(0.25);
        let a = [big + weak, Dd::ZERO];
        let b = [Dd::ONE, Dd::ZERO];
        let window = HarmonicWindow::new(0.0, 1.0).unwrap();
        let deterministic = window.product_kernel(0.0, 0.0).unwrap();
        let low = deterministic.mean_real_product_parts_dd(a, b).unwrap()
            - deterministic
                .mean_real_product_parts_dd([big, Dd::ZERO], b)
                .unwrap();
        close(low.to_f64(), 0.25, 0.0);
        let proper = window.proper_product_kernel(7.0, 7.0).unwrap();
        let low = proper.mean_real_product_parts_dd(a, b).unwrap()
            - proper
                .mean_real_product_parts_dd([big, Dd::ZERO], b)
                .unwrap();
        close(low.to_f64(), 0.125, 0.0);
        let invalid = [
            Dd {
                hi: 1.0,
                lo: f64::NAN,
            },
            Dd::ZERO,
        ];
        assert_eq!(
            deterministic.mean_real_product_parts_dd(invalid, b),
            Err(WindowError::NonFiniteAmplitude)
        );
        assert_eq!(
            proper.mean_real_product_parts_dd(invalid, b),
            Err(WindowError::NonFiniteAmplitude)
        );
        for frequency in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                window.product_kernel_hz(frequency, 1.0),
                Err(WindowError::InvalidFrequency)
            ));
            assert_eq!(
                window.mean_phase_hz(frequency),
                Err(WindowError::InvalidFrequency)
            );
        }
    }

    #[test]
    fn g1_deterministic_rotation_retains_a_sub_ulp_cosine_decrement() {
        let omega = 2.0_f64.powi(-40);
        let window = HarmonicWindow::new(1.0, 2.0_f64.powi(-40)).unwrap();
        let one = window
            .product_kernel(0.0, 0.0)
            .unwrap()
            .mean_real_product_dd(C64::ONE, C64::ONE)
            .unwrap();
        let cosine = window
            .product_kernel(omega, 0.0)
            .unwrap()
            .mean_real_product_dd(C64::ONE, C64::ONE)
            .unwrap();
        let decrement = (one - cosine).to_f64();
        let expected = 0.5 * omega * omega;
        assert!(decrement > 0.0);
        close(decrement, expected, expected * 3e-11);
    }
    #[test]
    fn g0_shaped_windows_are_normalized_and_uniform_is_unchanged() {
        let uniform = HarmonicWindow::new(-0.3, 0.7).unwrap();
        let explicit =
            HarmonicWindow::with_profile(-0.3, 0.7, HarmonicWindowProfile::Uniform).unwrap();
        assert_eq!(uniform, explicit);
        assert_eq!(uniform.profile(), HarmonicWindowProfile::Uniform);
        let triangular =
            HarmonicWindow::with_profile(-0.3, 0.7, HarmonicWindowProfile::Triangular).unwrap();
        assert_eq!(triangular.mean_phase(0.0).unwrap(), C64::ONE);
        assert_eq!(
            triangular
                .mean_real_product(C64::new(2.0, 17.0), 0.0, C64::new(3.0, -9.0), 0.0,)
                .unwrap(),
            6.0,
        );
        assert_eq!(triangular.start_s(), -0.3);
        assert_eq!(triangular.duration_s(), 0.7);
        for (start, duration) in [(0.0, 0.0), (f64::NAN, 1.0), (0.0, -1.0)] {
            assert_eq!(
                HarmonicWindow::with_profile(start, duration, HarmonicWindowProfile::Triangular,),
                Err(WindowError::InvalidWindow),
            );
        }
    }

    #[test]
    fn g1_triangular_products_match_weighted_real_signal_quadrature() {
        let a = C64::new(-0.7, 1.3);
        let b = C64::new(2.1, -0.4);
        let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
        for (start, duration, wa, wb) in [
            (-0.13, 0.027, 13.0, 17.0),
            (0.25, 0.3, 200.0, 201.0),
            (-0.2, 0.8, -19.0, 23.0),
            (0.15, 0.5, 0.0, 31.0),
        ] {
            let window =
                HarmonicWindow::with_profile(start, duration, HarmonicWindowProfile::Triangular)
                    .unwrap();
            let expected = simpson(window, a, wa, b, wb);
            let got = window.mean_real_product(a, wa, b, wb).unwrap();
            close(got, expected, 3e-11);
            close(window.mean_real_product(b, wb, a, wa).unwrap(), got, 3e-15);
            let proper_expected = phases
                .iter()
                .map(|&z| simpson(window, a * z, wa, b * z, wb))
                .sum::<f64>()
                * 0.25;
            let proper = window.proper_product_kernel(wa, wb).unwrap();
            close(
                proper.mean_real_product(a, b).unwrap(),
                proper_expected,
                3e-11,
            );
            // A unit DC partner extracts the weighted linear field mean.
            let mean = window.mean_phase(wa).unwrap();
            close(
                window.mean_real_product(a, wa, C64::ONE, 0.0).unwrap(),
                (a * mean).re,
                3e-15,
            );
        }
    }

    #[test]
    fn g1_triangular_shutter_preserves_weak_sine_and_disparate_frequency_energy() {
        let sine = C64::new(0.0, 1.0);
        for duration in [1e-4, 1e-8, 1e-12, 1e-50] {
            let window = HarmonicWindow::with_profile(
                -0.5 * duration,
                duration,
                HarmonicWindowProfile::Triangular,
            )
            .unwrap();
            // The normalized triangular clock has variance T²/24.
            let expected = duration * duration / 24.0;
            let actual = window.mean_real_product(sine, 1.0, sine, 1.0).unwrap();
            assert!(actual > 0.0);
            close(actual, expected, 1e-9 * expected);
        }
        let window =
            HarmonicWindow::with_profile(-1.0, 2.0, HarmonicWindowProfile::Triangular).unwrap();
        let small = 1e-40;
        let large = 3.0_f64;
        let x = 0.5 * large;
        // -d/db sinc(b/2)² = sinc(x)*(sin(x)-x*cos(x))/x².
        let expected = small * (x.sin() / x) * (x.sin() - x * x.cos()) / (x * x);
        close(
            window.mean_real_product(sine, small, sine, large).unwrap(),
            expected,
            2e-14 * expected.abs(),
        );
    }

    #[test]
    fn g1_triangular_proper_null_retains_sub_ulp_shutter_and_clock_decrements() {
        let delta = 2.0_f64.powi(-20);
        let duration = 2.0_f64.powi(-20);
        let omega = [1.0, 1.0 + delta];
        let amplitudes = [C64::ONE, -C64::ONE];
        for (start, coefficient) in [(-0.5 * duration, 1.0 / 48.0), (0.0, 7.0 / 48.0)] {
            let window =
                HarmonicWindow::with_profile(start, duration, HarmonicWindowProfile::Triangular)
                    .unwrap();
            let mut moment = Dd::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    moment = moment
                        + window
                            .proper_product_kernel(omega[i], omega[j])
                            .unwrap()
                            .mean_real_product_dd(amplitudes[i], amplitudes[j])
                            .unwrap();
                }
            }
            let expected = coefficient * (delta * duration).powi(2);
            assert!(moment.to_f64() > 0.0);
            close(moment.to_f64(), expected, 3e-14 * expected);
        }
    }

    #[test]
    fn g3_triangular_transform_clock_shift_and_high_carrier_hz_are_consistent() {
        let window =
            HarmonicWindow::with_profile(-0.5, 1.0, HarmonicWindowProfile::Triangular).unwrap();
        for omega in [0.0, 0.3, 4.0 * PI, 13.0] {
            let x = 0.25 * omega;
            let expected = if x == 0.0 { 1.0 } else { (x.sin() / x).powi(2) };
            let mean = window.mean_phase(omega).unwrap();
            close(mean.re, expected, 4e-16);
            close(mean.im, 0.0, 2e-16);
        }
        let carrier = 2.0_f64.powi(48);
        let difference = 0.125;
        let kernel = window
            .product_kernel_hz(carrier, carrier + difference)
            .unwrap();
        let value = kernel.mean_real_product_dd(C64::ONE, C64::ONE).unwrap()
            + kernel
                .mean_real_product_dd(C64::new(0.0, 1.0), C64::new(0.0, 1.0))
                .unwrap();
        let quarter_beat = 0.25 * TAU * difference;
        let expected = (quarter_beat.sin() / quarter_beat).powi(2);
        close(value.to_f64(), expected, 3e-14);
        assert_eq!(
            window
                .proper_product_kernel_hz(carrier + difference, carrier)
                .unwrap()
                .mean_real_product_dd(C64::ONE, C64::ONE)
                .unwrap(),
            window
                .proper_product_kernel_hz(difference, 0.0)
                .unwrap()
                .mean_real_product_dd(C64::ONE, C64::ONE)
                .unwrap(),
        );

        let shifted =
            HarmonicWindow::with_profile(-0.5 - 0.37, 1.0, HarmonicWindowProfile::Triangular)
                .unwrap();
        let a = C64::new(0.7, -1.2);
        let b = C64::new(-0.2, 0.8);
        let rotate = |z, omega: f64| z * C64::new((omega * 0.37).cos(), -(omega * 0.37).sin());
        close(
            shifted
                .mean_real_product(rotate(a, 3.0), 3.0, rotate(b, 7.0), 7.0)
                .unwrap(),
            window.mean_real_product(a, 3.0, b, 7.0).unwrap(),
            3e-15,
        );
    }

    #[test]
    fn g1_triangular_weights_the_observable_once_and_is_not_endpoint_storage() {
        let window =
            HarmonicWindow::with_profile(0.0, 1.0, HarmonicWindowProfile::Triangular).unwrap();
        let uniform = HarmonicWindow::new(0.0, 1.0).unwrap();
        // A constant physical field has unit mean square for BOTH detectors;
        // squaring the triangular sensitivity would incorrectly give 4/3.
        assert_eq!(
            window
                .mean_real_product(C64::ONE, 0.0, C64::ONE, 0.0)
                .unwrap(),
            1.0,
        );
        // U(t) = cos(t), dU/dt = -sin(t). The uniform integral is the
        // endpoint change; the shaped detector instead averages the local rate.
        let rate = C64::new(0.0, -1.0);
        let endpoint = 1.0_f64.cos() - 1.0;
        close(
            uniform.mean_real_product(rate, 1.0, C64::ONE, 0.0).unwrap(),
            endpoint,
            3e-16,
        );
        let weighted = window.mean_real_product(rate, 1.0, C64::ONE, 0.0).unwrap();
        close(weighted, simpson(window, rate, 1.0, C64::ONE, 0.0), 3e-13);
        assert!((weighted - endpoint).abs() > 1e-3);
    }
}
