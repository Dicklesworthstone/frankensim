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

mod band;
pub use band::{
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES, PreparedUniformHarmonicBand,
    UniformHarmonicBandRef, UniformHarmonicBandWorkspace,
};
mod spectrum;
pub use spectrum::{
    AngularHarmonicProductInput, HarmonicProductInput, HarmonicProductLaw,
    HarmonicProductSpectrumLimits, HarmonicProductTone, MAX_PRODUCT_SPECTRUM_TERMS,
    PreparedHarmonicProductSpectrum, ProductSpectrumError,
};
pub use spectrum::{
    AngularHarmonicSpectrumTerm, HarmonicCausalIntegral, HarmonicSpectrumTermLimits,
};

pub use spectrum::{
    HarmonicDecayIntervalLimits, HarmonicDecayMoments, HarmonicPhaseInterval,
    MAX_HARMONIC_DECAY_TONES, PreparedHarmonicDecayInterval, phase_angular_interval_from_clock_dd,
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
    /// Field or scratch dimensions do not match the prepared coherent band.
    InvalidBandShape,
    /// A phase product, amplitude rotation or result exceeds floating range.
    Unrepresentable,
}

impl core::fmt::Display for WindowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidWindow => "harmonic window needs a finite start/end and positive duration",
            Self::InvalidFrequency => "harmonic angular frequency must be finite",
            Self::NonFiniteAmplitude => "harmonic peak amplitude must be finite",
            Self::InvalidBandShape => "harmonic band field/scratch dimensions do not match",
            Self::Unrepresentable => "harmonic window moment is not representable",
        })
    }
}

impl std::error::Error for WindowError {}

/// Instantaneous `exp(-i 2 pi f t)` on an explicitly retained source clock.
///
/// Signed Hz and DC are admitted. Reduce the retained cycle count before
/// conversion to radians, so the approximation to 2*pi cannot accumulate a
/// spurious phase on a large absolute source clock. The returned real/imaginary
/// parts retain small cancellation residuals until the final field reduction.
pub fn phase_hz_from_clock_dd(frequency_hz: f64, time_s: Dd) -> Result<[Dd; 2], WindowError> {
    if !frequency_hz.is_finite() {
        return Err(WindowError::InvalidFrequency);
    }
    if !time_s.hi.is_finite() || !time_s.lo.is_finite() {
        return Err(WindowError::InvalidWindow);
    }
    let time = checked(Dd::from_pair(time_s.hi, time_s.lo))?;
    phase_cycles(Dd::from_f64(frequency_hz), time, 1.0)
}

/// Instantaneous `exp(-i omega t)` for an originally angular-frequency source.
///
/// Signed rad/s and DC are admitted. The source clock and complex phase remain
/// in DD; no finite shutter or rounded conversion through Hz is introduced.
pub fn phase_rad_s_from_clock_dd(
    angular_frequency_rad_s: f64,
    time_s: Dd,
) -> Result<[Dd; 2], WindowError> {
    phase_rad_s_dd_from_clock_dd(Dd::from_f64(angular_frequency_rad_s), time_s)
}

/// Native angular phase retaining a DD frequency, including an exact retained
/// difference between two authored source frequencies. Proper field sums can
/// remove a common carrier before evaluating relative phases; actual real
/// fields must also retain that carrier's final rotation.
///
/// This uses the original deterministic DD radian reduction. It is not a
/// relative-error certificate for arbitrarily cancelled coherent fields.
pub fn phase_rad_s_dd_from_clock_dd(
    angular_frequency_rad_s: Dd,
    time_s: Dd,
) -> Result<[Dd; 2], WindowError> {
    let omega = admit_frequency_dd(angular_frequency_rad_s)?;
    if !time_s.hi.is_finite() || !time_s.lo.is_finite() {
        return Err(WindowError::InvalidWindow);
    }
    let time = checked(Dd::from_pair(time_s.hi, time_s.lo))?;
    let (real, imaginary) = proper_phase(checked(omega * time)?);
    Ok([checked(real)?, checked(imaginary)?])
}

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

    /// Mean phase for a frequency authored in cycles/second. Reduce exact
    /// cycle products before conversion to radians, retaining the original
    /// source clock without accumulating the error in the 2*pi constant.
    pub fn mean_phase_hz(self, frequency_hz: f64) -> Result<C64, WindowError> {
        if !frequency_hz.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let frequency = Dd::from_f64(frequency_hz);
        let scale = self.centered_transform_hz(frequency)?;
        let [real, imaginary] = self.rotation_hz(frequency)?;
        Ok(C64::new(
            checked(real * scale)?.to_f64(),
            checked(imaginary * scale)?.to_f64(),
        ))
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

    /// Real-product kernel for angular frequencies retained in double-double
    /// precision. A common carrier and a much smaller resolved offset must
    /// not be rounded to one binary64 value before the sum/difference Gram
    /// weights are formed. Inputs are normalized without discarding either
    /// finite component; the original clock, shutter and contraction apply.
    pub fn product_kernel_dd(
        self,
        omega_a: Dd,
        omega_b: Dd,
    ) -> Result<HarmonicProduct, WindowError> {
        let omega_a = admit_frequency_dd(omega_a)?;
        let omega_b = admit_frequency_dd(omega_b)?;
        if omega_a.lo == 0.0 && omega_b.lo == 0.0 {
            return self.product_kernel_extended(omega_a, omega_b);
        }
        let a = self.half_width_extended(omega_a)?;
        let b = self.half_width_extended(omega_b)?;
        // Subtract BEFORE multiplying by the interval. The small beat can
        // otherwise disappear into the multiplication error of the common
        // carrier, even though both authored DD frequencies retain it.
        let difference = self.half_width_extended(checked(omega_a - omega_b)?)?;
        let sum = self.half_width_extended(checked(omega_a + omega_b)?)?;
        let (cos_cos, sin_sin) = match self.profile {
            HarmonicWindowProfile::Uniform => rectangular_gram_parts(a, b, difference, sum)?,
            HarmonicWindowProfile::Triangular => {
                let half = Dd::from_f64(0.5);
                let (c, s) =
                    rectangular_gram_parts(a * half, b * half, difference * half, sum * half)?;
                (checked(c * c + s * s)?, checked(Dd::from_f64(2.0) * c * s)?)
            }
        };
        Ok(HarmonicProduct {
            rotation_a: self.rotation_parts(omega_a)?,
            rotation_b: self.rotation_parts(omega_b)?,
            cos_cos,
            sin_sin,
        })
    }

    /// Deterministic product for two offsets on one common angular carrier.
    ///
    /// Each physical frequency is common_offset + omega_a/b. Add in
    /// double-double precision before the clock rotation and centered Gram
    /// contraction, so a small envelope offset is not rounded out of a large
    /// carrier. Both sum and difference terms are retained. Signed offsets,
    /// carriers and DC are admitted; all three inputs must be finite.
    ///
    /// This is the same analytic owner as product_kernel, not a carrier
    /// average or an approximation that removes the sum-frequency term.
    pub fn product_kernel_shifted(
        self,
        omega_a: f64,
        omega_b: f64,
        common_offset: f64,
    ) -> Result<HarmonicProduct, WindowError> {
        if !omega_a.is_finite() || !omega_b.is_finite() || !common_offset.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let carrier = Dd::from_f64(common_offset);
        self.product_kernel_dd(
            checked(carrier + Dd::from_f64(omega_a))?,
            checked(carrier + Dd::from_f64(omega_b))?,
        )
    }

    /// Deterministic real-product kernel for frequencies authored in Hz.
    ///
    /// Forms sum/difference frequencies in Hz and reduces clock cycles before
    /// converting to radians. Separately converted carriers can otherwise
    /// corrupt a resolved beat or drift in phase on a large source epoch.
    pub fn product_kernel_hz(
        self,
        frequency_a_hz: f64,
        frequency_b_hz: f64,
    ) -> Result<HarmonicProduct, WindowError> {
        if !frequency_a_hz.is_finite() || !frequency_b_hz.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let a = Dd::from_f64(frequency_a_hz);
        let b = Dd::from_f64(frequency_b_hz);
        let (cos_cos, sin_sin) = match self.profile {
            HarmonicWindowProfile::Uniform => self.rectangular_gram_hz(a, b, 2.0)?,
            HarmonicWindowProfile::Triangular => {
                let (c, s) = self.rectangular_gram_hz(a, b, 4.0)?;
                (checked(c * c + s * s)?, checked(Dd::from_f64(2.0) * c * s)?)
            }
        };
        Ok(HarmonicProduct {
            rotation_a: self.rotation_hz(a)?,
            rotation_b: self.rotation_hz(b)?,
            cos_cos,
            sin_sin,
        })
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

    /// Proper second-moment kernel for double-double angular frequencies.
    /// Subtract the retained frequencies before applying the physical clock,
    /// preserving beats smaller than a binary64 step of their common carrier.
    /// The result includes the same peak-amplitude factor 1/2 as
    /// [`Self::proper_product_kernel`].
    pub fn proper_product_kernel_dd(
        self,
        omega_a: Dd,
        omega_b: Dd,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        self.proper_difference_kernel(checked(
            admit_frequency_dd(omega_a)? - admit_frequency_dd(omega_b)?,
        )?)
    }

    /// Proper second-moment kernel for frequencies authored in cycles/second.
    ///
    /// The difference is formed in double-double BEFORE conversion to rad/s.
    /// Separately rounding `TAU*f_a` and `TAU*f_b` can corrupt a resolved beat
    /// between close carriers; rounding the difference before clock evaluation
    /// can also discard the phase of a much smaller component. This method
    /// retains the difference and reduces its cycles before radian conversion.
    pub fn proper_product_kernel_hz(
        self,
        frequency_a_hz: f64,
        frequency_b_hz: f64,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        if !frequency_a_hz.is_finite() || !frequency_b_hz.is_finite() {
            return Err(WindowError::InvalidFrequency);
        }
        let difference = checked(Dd::from_f64(frequency_a_hz) - Dd::from_f64(frequency_b_hz))?;
        let scale = checked(Dd::from_f64(0.5) * self.centered_transform_hz(difference)?)?;
        let [real, imaginary] = self.rotation_hz(difference)?;
        Ok(ProperHarmonicProduct {
            real: checked(real * scale)?,
            imaginary: checked(imaginary * scale)?,
        })
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

    fn centered_transform_hz(self, frequency: Dd) -> Result<Dd, WindowError> {
        match self.profile {
            HarmonicWindowProfile::Uniform => self.sinc_hz(frequency, 2.0),
            HarmonicWindowProfile::Triangular => {
                let half_sinc = self.sinc_hz(frequency, 4.0)?;
                checked(half_sinc * half_sinc)
            }
        }
    }

    // The divisor is 2 for a rectangular interval and 4 for each half-width
    // interval whose convolution gives the triangular profile. It is applied
    // after the cycle product, preserving subnormal authored frequencies/times.
    fn width_hz(self, frequency: Dd, divisor: f64) -> Result<Dd, WindowError> {
        checked(self.half_width_extended(frequency)? * (TAU_DD * Dd::from_f64(2.0 / divisor)))
    }

    fn sinc_hz(self, frequency: Dd, divisor: f64) -> Result<Dd, WindowError> {
        let width = self.width_hz(frequency, divisor)?;
        if width.hi.abs() <= 0.5 {
            Ok(Dd::ONE + sinc_minus_one(width))
        } else {
            let phase = phase_cycles(frequency, Dd::from_f64(self.duration_s), divisor)?;
            checked(-phase[1] / width)
        }
    }

    fn rectangular_gram_hz(
        self,
        frequency_a: Dd,
        frequency_b: Dd,
        divisor: f64,
    ) -> Result<(Dd, Dd), WindowError> {
        let difference = checked(frequency_a - frequency_b)?;
        let sum = checked(frequency_a + frequency_b)?;
        let a = self.width_hz(frequency_a, divisor)?;
        let b = self.width_hz(frequency_b, divisor)?;
        let d = self.width_hz(difference, divisor)?;
        let s = self.width_hz(sum, divisor)?;
        if d.hi.abs().max(s.hi.abs()) <= 0.5 {
            // Retain the sinc-1 difference that resolves a short shutter's
            // small positive sine energy. Only the large-phase path changes.
            return rectangular_gram_parts(a, b, d, s);
        }
        let d = self.sinc_hz(difference, divisor)?;
        let s = self.sinc_hz(sum, divisor)?;
        let (small, large, fs, fl) = if a.hi.abs() <= b.hi.abs() {
            (a, b, frequency_a, frequency_b)
        } else {
            (b, a, frequency_b, frequency_a)
        };
        let ratio = small / large;
        let half = Dd::from_f64(0.5);
        let sin_sin = if ratio.hi.abs() < 1e-4 {
            let time = Dd::from_f64(self.duration_s);
            let cs = phase_cycles(fs, time, divisor)?[0];
            let cl = phase_cycles(fl, time, divisor)?[0];
            ratio * (self.sinc_hz(fl, divisor)? * cs - self.sinc_hz(fs, divisor)? * cl)
                / (Dd::ONE - ratio * ratio)
        } else {
            half * (d - s)
        };
        Ok((checked(half * (d + s))?, checked(sin_sin)?))
    }

    fn rotation_hz(self, frequency: Dd) -> Result<[Dd; 2], WindowError> {
        let a = phase_cycles(frequency, Dd::from_f64(self.start_s), 1.0)?;
        let b = phase_cycles(frequency, Dd::from_f64(self.duration_s), 2.0)?;
        Ok([
            checked(a[0] * b[0] - a[1] * b[1])?,
            checked(a[0] * b[1] + a[1] * b[0])?,
        ])
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

    fn rotation_parts(self, omega: Dd) -> Result<[Dd; 2], WindowError> {
        let high = Dd::from_f64(omega.hi);
        let a = self.rotation(high, self.half_width_extended(high)?)?;
        if omega.lo == 0.0 {
            return Ok(a);
        }
        // The high-frequency product already needs both DD components.
        // Rotate the retained offset separately before combining phases;
        // adding its small angle to that product can erase the offset.
        let low = Dd::from_f64(omega.lo);
        let b = self.rotation(low, self.half_width_extended(low)?)?;
        Ok([
            checked(a[0] * b[0] - a[1] * b[1])?,
            checked(a[0] * b[1] + a[1] * b[0])?,
        ])
    }
}

// Centered cosine/cosine and sine/sine means under a uniform clock.
fn rectangular_gram(a: Dd, b: Dd) -> Result<(Dd, Dd), WindowError> {
    let difference = checked(a - b)?;
    let sum = checked(a + b)?;
    rectangular_gram_parts(a, b, difference, sum)
}

fn rectangular_gram_parts(a: Dd, b: Dd, difference: Dd, sum: Dd) -> Result<(Dd, Dd), WindowError> {
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
            let ca = det::sin_cos_dd(small).1;
            let cb = det::sin_cos_dd(large).1;
            ratio * (sinc(large) * ca - sinc(small) * cb) / (Dd::ONE - ratio * ratio)
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
    /// Centered cosine/cosine and sine/sine Gram weights, before the two
    /// carrier rotations. Symmetry makes centered mixed cosine/sine means zero.
    ///
    /// Separable spatial/temporal observers may compose these weights and
    /// apply one joint center phase to retained DD field amplitudes. The
    /// original source clock and each interval midpoint must still be applied;
    /// these weights alone do not describe the original uncentered product.
    pub fn centered_gram_dd(self) -> [Dd; 2] {
        [self.cos_cos, self.sin_sin]
    }

    /// Average `Re[a exp(-i omega_a t)] * Re[b exp(-i omega_b t)]`.
    pub fn mean_real_product(self, a: C64, b: C64) -> Result<f64, WindowError> {
        self.mean_real_product_dd(a, b).map(Dd::to_f64)
    }

    /// Retain the double-double product until a whole field's interference
    /// sum is complete. Rounding each large pair to f64 before cancellation
    /// can discard a small beat signal even with an accurate outer sum.
    /// Cached trigonometry remains DD; no relative error bound near an
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
        let real = ar * br + ai * bi;
        let imaginary = ai * br - ar * bi;
        if let (Ok(real), Ok(imaginary)) = (checked(real), checked(imaginary)) {
            return checked(checked(real * self.real)? - checked(imaginary * self.imaginary)?);
        }
        // The kernel already contains the physical peak factor 1/2. Move it
        // into an overflowing field product, preserving the ordinary path's
        // DD contraction and decorrelation when both components fit.
        let real = half_product_sum(ar, br, ai, bi, false)?;
        let imaginary = half_product_sum(ai, br, ar, bi, true)?;
        let two = Dd::from_f64(2.0);
        checked(checked(real * (two * self.real))? - checked(imaginary * (two * self.imaginary))?)
    }
}

fn admit_frequency_dd(frequency: Dd) -> Result<Dd, WindowError> {
    if !frequency.hi.is_finite() || !frequency.lo.is_finite() {
        return Err(WindowError::InvalidFrequency);
    }
    checked(Dd::from_pair(frequency.hi, frequency.lo))
}

fn checked(x: Dd) -> Result<Dd, WindowError> {
    if x.hi.is_finite() && x.lo.is_finite() {
        Ok(x)
    } else {
        Err(WindowError::Unrepresentable)
    }
}

// Physical peak normalization of (a*b +/- c*d). Keep a representable sum
// before halving: two half-subnormal products can have a representable sum.
// Only an overflowing intermediate needs the half moved into each product;
// halve its larger operand so a tiny operand beside a large one survives.
fn half_product_sum(a: Dd, b: Dd, c: Dd, d: Dd, subtract: bool) -> Result<Dd, WindowError> {
    let half = Dd::from_f64(0.5);
    let sum = if subtract {
        a * b - c * d
    } else {
        a * b + c * d
    };
    if let Ok(sum) = checked(sum) {
        return checked(half * sum);
    }
    let product = |left: Dd, right: Dd| {
        checked(if left.abs().lt(right.abs()) {
            left * (half * right)
        } else {
            (half * left) * right
        })
    };
    let first = product(a, b)?;
    let second = product(c, d)?;
    checked(if subtract {
        first - second
    } else {
        first + second
    })
}

// Mathematical 2*pi to DD precision, using fs-math's existing pi-tail
// convention. Lifting only f64::TAU leaves its ~2.45e-16 radian error in every
// cycle; multiplying by an absolute source epoch then creates a spurious phase.
pub(crate) const TAU_DD: Dd = Dd {
    hi: core::f64::consts::TAU,
    lo: 2.449_293_598_294_706_4e-16,
};

/// Convert authored hertz to angular frequency with the mathematical 2*pi
/// constant retained in DD. A source already authored in rad/s must instead
/// retain that original value; converting through f64 Hz changes its clock.
pub fn angular_frequency(frequency_hz: f64) -> Result<Dd, WindowError> {
    if !frequency_hz.is_finite() {
        return Err(WindowError::InvalidFrequency);
    }
    checked(Dd::from_f64(frequency_hz) * TAU_DD)
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

// Exp(-i*2*pi*f*t/divisor), where divisor is an exact power of two (1, 2, 4).
// Each scalar product is an exact two-product before range reduction. Reduce
// its high/low components independently: adding a subperiod clock offset to a
// huge carrier product first would erase it, even in DD. Radian-authored inputs
// deliberately do not use this cycle identity or change their physical clock.
#[allow(clippy::float_cmp)] // Exact zero skips a mathematically absent product.
fn phase_cycles(frequency: Dd, time: Dd, divisor: f64) -> Result<[Dd; 2], WindowError> {
    let remainder = |x: f64| x - (x / divisor).trunc() * divisor;
    let mut cycles = Dd::ZERO;
    for f in [frequency.hi, frequency.lo] {
        for t in [time.hi, time.lo] {
            if f == 0.0 || t == 0.0 {
                continue;
            }
            let product = checked(Dd::from_f64(f) * Dd::from_f64(t))?;
            cycles = cycles + Dd::from_pair(remainder(product.hi), remainder(product.lo));
        }
    }
    let cycles = Dd::from_pair(remainder(cycles.hi), remainder(cycles.lo));
    let phase = checked(cycles * (TAU_DD * Dd::from_f64(1.0 / divisor)))?;
    let (real, imaginary) = proper_phase(phase);
    Ok([checked(real)?, checked(imaginary)?])
}

// Exp(-i*x), retaining the small real decrement near ±1 instead of rounding
// independently evaluated sinc/carrier terms to f64 before interference.
// The shared fs-math owner also preserves low sine components: DD products
// alone cannot repair a Gram matrix whose inputs were already rounded.
fn proper_phase(x: Dd) -> (Dd, Dd) {
    let (sine, cosine) = det::sin_cos_dd(x);
    (cosine, -sine)
}

fn sinc(x: Dd) -> Dd {
    if x.hi.abs() <= 0.5 {
        Dd::ONE + sinc_minus_one(x)
    } else {
        det::sin_cos_dd(x).0 / x
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

    #[test]
    fn g1_native_angular_instant_phase_retains_frequency_and_clock_parts() {
        let compose =
            |a: [Dd; 2], b: [Dd; 2]| [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
        for (omega, epoch, offset) in [
            (512.0, 2.0_f64.powi(40), 2.0_f64.powi(-18)),
            (2.0_f64.powi(60), 1.0, 2.0_f64.powi(-61)),
        ] {
            let actual = phase_rad_s_from_clock_dd(omega, Dd::from_pair(epoch, offset)).unwrap();
            let expected = compose(
                phase_rad_s_from_clock_dd(omega, Dd::from_f64(epoch)).unwrap(),
                phase_rad_s_from_clock_dd(omega, Dd::from_f64(offset)).unwrap(),
            );
            for c in 0..2 {
                close((actual[c] - expected[c]).to_f64(), 0.0, 2e-15);
            }
            let rounded = phase_rad_s_from_clock_dd(omega, Dd::from_f64(epoch)).unwrap();
            assert!(
                (actual[0] - rounded[0]).abs().to_f64() + (actual[1] - rounded[1]).abs().to_f64()
                    > 1e-4
            );
        }
        let carrier = 2.0_f64.powi(60);
        let difference = Dd::from_f64(carrier) - Dd::ONE;
        let actual = phase_rad_s_dd_from_clock_dd(difference, Dd::ONE).unwrap();
        let expected = compose(
            phase_rad_s_from_clock_dd(carrier, Dd::ONE).unwrap(),
            [Dd::from_f64(1.0_f64.cos()), Dd::from_f64(1.0_f64.sin())],
        );
        for c in 0..2 {
            close((actual[c] - expected[c]).to_f64(), 0.0, 2e-15);
        }
    }

    #[test]
    fn g0_native_angular_instant_phase_sign_small_residual_and_admission() {
        let positive = phase_rad_s_from_clock_dd(1.0, Dd::from_f64(0.25)).unwrap();
        close(positive[0].to_f64(), 0.25_f64.cos(), 2e-16);
        close(positive[1].to_f64(), -0.25_f64.sin(), 2e-16);
        let negative = phase_rad_s_from_clock_dd(-1.0, Dd::from_f64(0.25)).unwrap();
        assert_eq!(positive[0], negative[0]);
        assert_eq!(positive[1], -negative[1]);
        let tiny = phase_rad_s_from_clock_dd(1.0, Dd::from_f64(1e-12)).unwrap();
        close((Dd::ONE - tiny[0]).to_f64() / 5e-25, 1.0, 1e-7);
        close(tiny[1].to_f64() / -1e-12, 1.0, 2e-16);
        assert_eq!(
            phase_rad_s_from_clock_dd(0.0, Dd::from_pair(1e100, 1.0)).unwrap()[0],
            Dd::ONE
        );
        assert!(matches!(
            phase_rad_s_from_clock_dd(f64::NAN, Dd::ZERO),
            Err(WindowError::InvalidFrequency)
        ));
        assert!(matches!(
            phase_rad_s_from_clock_dd(1.0, Dd::from_f64(f64::INFINITY)),
            Err(WindowError::InvalidWindow)
        ));
        assert!(matches!(
            phase_rad_s_from_clock_dd(f64::MAX, Dd::from_f64(2.0)),
            Err(WindowError::Unrepresentable)
        ));
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
    fn g1_shifted_kernel_retains_offsets_below_carrier_spacing() {
        let carrier = 1e20;
        let offsets = [-1.0, 1.0];
        assert_eq!(carrier + offsets[0], carrier + offsets[1]);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.13, 0.7, profile).unwrap();
            let kernel = window
                .product_kernel_shifted(offsets[0], offsets[1], carrier)
                .unwrap();
            // cos(A)cos(B) + sin(A)sin(B) = cos(A-B), pointwise.
            // Independent difference-frequency reference has no huge carrier.
            let resolved = kernel.mean_real_product_dd(C64::ONE, C64::ONE).unwrap()
                + kernel
                    .mean_real_product_dd(C64::new(0.0, 1.0), C64::new(0.0, 1.0))
                    .unwrap();
            let reference = window.mean_phase(offsets[0] - offsets[1]).unwrap().re;
            // DD retains the beat, with finite absolute phase precision at 1e20.
            let tolerance = 4e-14 + 8.0 * carrier * f64::EPSILON * f64::EPSILON;
            close(resolved.to_f64(), reference, tolerance);
            assert!((resolved.to_f64() - 1.0).abs() > 0.1);
        }
    }

    #[test]
    fn g0_shifted_kernel_zero_offset_and_refusals() {
        let window = HarmonicWindow::new(-0.3, 0.8).unwrap();
        let a = C64::new(0.7, -0.4);
        let b = C64::new(-0.2, 0.9);
        assert_eq!(
            window
                .product_kernel_shifted(-2.0, 3.0, 0.0)
                .unwrap()
                .mean_real_product_dd(a, b)
                .unwrap(),
            window
                .product_kernel(-2.0, 3.0)
                .unwrap()
                .mean_real_product_dd(a, b)
                .unwrap(),
        );
        for (wa, wb, carrier) in [
            (f64::NAN, 1.0, 0.0),
            (1.0, f64::INFINITY, 0.0),
            (1.0, 2.0, f64::NEG_INFINITY),
        ] {
            assert_eq!(
                window.product_kernel_shifted(wa, wb, carrier).unwrap_err(),
                WindowError::InvalidFrequency,
            );
        }
        assert_eq!(
            window
                .product_kernel_shifted(f64::MAX, 1.0, f64::MAX)
                .unwrap_err(),
            WindowError::Unrepresentable,
        );
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
    fn g1_nearby_coherent_tones_retain_positive_sub_ulp_beat_energy() {
        // Integrate the REAL difference signal through the stable identity
        // cos(a)-cos(b) = 2*sin((a+b)/2)*sin((b-a)/2). This oracle never
        // subtracts squared carrier powers or calls the moment kernels.
        for (start, duration, separation) in [
            (0.0, 0.0001, 2.0_f64.powi(-20)),
            (0.0, 0.0011, 2.0_f64.powi(-30)),
            (0.2345, 0.0011, 2.0_f64.powi(-30)),
        ] {
            let frequencies = [1000.0, 1000.0 + separation];
            for profile in [
                HarmonicWindowProfile::Uniform,
                HarmonicWindowProfile::Triangular,
            ] {
                let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
                let mut actual = Dd::ZERO;
                for i in 0..2 {
                    for j in 0..2 {
                        let amplitude = |index| if index == 0 { C64::ONE } else { -C64::ONE };
                        actual = actual
                            + window
                                .product_kernel_hz(frequencies[i], frequencies[j])
                                .unwrap()
                                .mean_real_product_dd(amplitude(i), amplitude(j))
                                .unwrap();
                    }
                }
                let mut expected = 0.0;
                let intervals = 20_000;
                for i in 0..=intervals {
                    let u = f64::from(i) / f64::from(intervals);
                    let time = start + duration * u;
                    let difference = 2.0
                        * (TAU * (1000.0 + 0.5 * separation) * time).sin()
                        * (PI * separation * time).sin();
                    let sensitivity = match profile {
                        HarmonicWindowProfile::Uniform => 1.0,
                        HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
                    };
                    let weight = if i == 0 || i == intervals {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    expected += weight * sensitivity * difference * difference;
                }
                expected /= 3.0 * f64::from(intervals);
                assert!(
                    actual.to_f64() > 0.0,
                    "negative coherent energy: {actual:?}"
                );
                close(actual.to_f64(), expected, expected * 2e-8);
            }
        }
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
    fn g1_dd_product_preserves_a_beat_below_the_carrier_spacing() {
        let duration = 2.0_f64.powi(60);
        let delta = 1.0 / duration;
        assert_eq!(1.0 + delta, 1.0);
        let window = HarmonicWindow::new(0.0, duration).unwrap();
        let actual = window
            .product_kernel_dd(Dd::from_pair(1.0, delta), Dd::ONE)
            .unwrap()
            .mean_real_product(C64::ONE, C64::ONE)
            .unwrap();
        // cos((1+delta)t)*cos(t) has a resolved beat delta*T=1.
        // Its sum-frequency mean has magnitude <= 1/(2*(2+delta)*T),
        // below 3e-19 here; no high-frequency quadrature is needed.
        close(actual, 0.5 * 1.0_f64.sin(), 3e-15);
        let rounded = window
            .product_kernel(1.0 + delta, 1.0)
            .unwrap()
            .mean_real_product(C64::ONE, C64::ONE)
            .unwrap();
        assert!((rounded - actual).abs() > 0.07);
    }

    #[test]
    fn g1_dd_real_beats_survive_the_carrier_multiplication_error() {
        for power in [60, 80, 90, 100] {
            let duration = 0.7 * 2.0_f64.powi(power);
            let delta = 2.0_f64.powi(-power);
            for profile in [
                HarmonicWindowProfile::Uniform,
                HarmonicWindowProfile::Triangular,
            ] {
                let window = HarmonicWindow::with_profile(0.0, duration, profile).unwrap();
                let actual = window
                    .product_kernel_dd(Dd::from_pair(0.3, delta), Dd::from_f64(0.3))
                    .unwrap()
                    .mean_real_product(C64::ONE, C64::ONE)
                    .unwrap();
                assert_eq!(
                    actual,
                    window
                        .product_kernel_shifted(delta, 0.0, 0.3)
                        .unwrap()
                        .mean_real_product(C64::ONE, C64::ONE)
                        .unwrap(),
                );
                // delta*T=0.7 exactly at the authored binary64 precision.
                // The real cross product is half the mean cosine of that
                // beat plus a fast sum term bounded below 2e-18 here.
                let expected = match profile {
                    HarmonicWindowProfile::Uniform => 0.5 * 0.7_f64.sin() / 0.7,
                    HarmonicWindowProfile::Triangular => {
                        0.5 * 0.35_f64.cos() * (0.175_f64.sin() / 0.175).powi(2)
                    }
                };
                close(actual, expected, 3e-15);
            }
        }
    }

    #[test]
    fn g3_dd_proper_beats_cancel_the_retained_common_carrier() {
        let delta = 2.0_f64.powi(-60);
        let window = HarmonicWindow::new(2.0_f64.powi(90), 2.0_f64.powi(59)).unwrap();
        let a = C64::new(0.7, -1.2);
        let b = C64::new(-0.2, 0.8);
        let expected = window
            .proper_product_kernel(delta, 0.0)
            .unwrap()
            .mean_real_product_dd(a, b)
            .unwrap();
        for carrier in [-1.0, 0.0, 1.0, 2.0_f64.powi(30)] {
            let actual = window
                .proper_product_kernel_dd(Dd::from_pair(carrier, delta), Dd::from_f64(carrier))
                .unwrap()
                .mean_real_product_dd(a, b)
                .unwrap();
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn g0_dd_product_admission_and_binary64_compatibility() {
        let a = C64::new(0.7, -0.3);
        let b = C64::new(-0.4, 1.1);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(-0.17, 0.31, profile).unwrap();
            for (wa, wb) in [(0.0, 0.0), (-1.7, 2.3), (1.0, 1.0), (1e12, 1e12 + 0.25)] {
                assert_eq!(
                    window
                        .product_kernel_dd(Dd::from_f64(wa), Dd::from_f64(wb))
                        .unwrap()
                        .mean_real_product_dd(a, b)
                        .unwrap(),
                    window
                        .product_kernel(wa, wb)
                        .unwrap()
                        .mean_real_product_dd(a, b)
                        .unwrap(),
                );
                assert_eq!(
                    window
                        .proper_product_kernel_dd(Dd::from_f64(wa), Dd::from_f64(wb))
                        .unwrap()
                        .mean_real_product_dd(a, b)
                        .unwrap(),
                    window
                        .proper_product_kernel(wa, wb)
                        .unwrap()
                        .mean_real_product_dd(a, b)
                        .unwrap(),
                );
            }
            assert_eq!(
                window
                    .product_kernel_dd(Dd { hi: 1.0, lo: 2.0 }, Dd::ONE)
                    .unwrap()
                    .mean_real_product_dd(a, b)
                    .unwrap(),
                window
                    .product_kernel(3.0, 1.0)
                    .unwrap()
                    .mean_real_product_dd(a, b)
                    .unwrap(),
            );
            for invalid in [
                Dd {
                    hi: f64::NAN,
                    lo: 0.0,
                },
                Dd {
                    hi: 0.0,
                    lo: f64::INFINITY,
                },
            ] {
                assert!(matches!(
                    window.product_kernel_dd(invalid, Dd::ONE),
                    Err(WindowError::InvalidFrequency)
                ));
                assert!(matches!(
                    window.proper_product_kernel_dd(Dd::ONE, invalid),
                    Err(WindowError::InvalidFrequency)
                ));
            }
        }
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
        // This rad/s input is the rounded f64 approximation to 2*pi. The Hz
        // owner retains the mathematical constant's low part, so the two
        // independently authored frequencies agree numerically, not bitwise.
        close(
            large.to_f64(),
            window
                .proper_product_kernel(TAU, 0.0)
                .unwrap()
                .mean_real_product(a, b)
                .unwrap(),
            2e-15,
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
        // The base completes an EXACT integer number of cycles at both the
        // start and half-width. Its phase is therefore one; only the fractional
        // frequency rotates. Applying platform trig to TAU*large_integer would
        // instead validate the growing phase error of a rounded 2*pi constant.
        let rotate = |x: f64| C64::new(x.cos(), -x.sin());
        assert_eq!(base * start, 2.0_f64.powi(28));
        assert_eq!(base * (0.5 * duration), 2.0_f64.powi(17));
        let small_start = TAU * (fraction * start);
        let small_half = TAU * (fraction * (0.5 * duration));
        let expected_rotation = rotate(small_start) * rotate(small_half);
        let half = TAU * ((base + fraction) * (0.5 * duration));
        let expected = expected_rotation.scale(small_half.sin() / half);
        let actual = window.mean_phase_hz(base + fraction).unwrap();
        close(actual.re, expected.re, 2e-20);
        close(actual.im, expected.im, 2e-20);
    }

    #[test]
    fn g3_instantaneous_hz_phase_retains_integer_cycles_and_sub_ulp_quarter_turns() {
        for epoch in [2.0_f64.powi(40), -2.0_f64.powi(40)] {
            let time = Dd::from_f64(epoch);
            let phase = phase_hz_from_clock_dd(4096.0, time).unwrap();
            close(phase[0].to_f64(), 1.0, 2e-15);
            close(phase[1].to_f64(), 0.0, 2e-15);
            // A quarter period lies below one f64 ULP of the absolute epoch.
            let quarter = time + Dd::from_f64(1.0 / 16384.0);
            for (frequency, imaginary) in [(4096.0, -1.0), (-4096.0, 1.0)] {
                let phase = phase_hz_from_clock_dd(frequency, quarter).unwrap();
                close(phase[0].to_f64(), 0.0, 2e-15);
                close(phase[1].to_f64(), imaginary, 2e-15);
            }
        }
    }

    #[test]
    fn g3_hz_phase_reduces_cycles_before_the_irrational_conversion() {
        for power in [40, 80, 120, 600] {
            for sign in [-1.0, 1.0] {
                let epoch = Dd::from_f64(sign * 2.0_f64.powi(power));
                let full_cycles = phase_hz_from_clock_dd(4096.0, epoch).unwrap();
                assert_eq!(full_cycles, [Dd::ONE, Dd::ZERO]);
                // The quarter turn must not be added to an enormous radian
                // phase before range reduction. It remains physically real
                // even when neither f64 nor the DD carrier sum can retain it.
                let quarter = Dd::from_pair(epoch.hi, 1.0 / 16384.0);
                for (frequency, imaginary) in [(4096.0, -1.0), (-4096.0, 1.0)] {
                    let phase = phase_hz_from_clock_dd(frequency, quarter).unwrap();
                    close(phase[0].to_f64(), 0.0, 2e-30);
                    close(phase[1].to_f64(), imaginary, 2e-30);
                }
            }
        }
    }

    #[test]
    fn g1_hz_optical_carrier_beats_remain_on_the_original_clock() {
        let carrier = 1e15;
        let duration = 0.001;
        let half_beat = PI * duration;
        for power in [20, 30, 40, 50, 80] {
            let start = 0.7 * 2.0_f64.powi(power);
            for profile in [
                HarmonicWindowProfile::Uniform,
                HarmonicWindowProfile::Triangular,
            ] {
                let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
                let kernel = window.product_kernel_hz(carrier + 1.0, carrier).unwrap();
                // A quadrature pair eliminates the sum frequency identically:
                // cos(a)cos(b)+sin(a)sin(b)=cos(a-b). The remaining 1 Hz beat
                // can be integrated without any high carrier in the oracle.
                let actual = kernel.mean_real_product_dd(C64::ONE, C64::ONE).unwrap()
                    + kernel
                        .mean_real_product_dd(C64::new(0.0, 1.0), C64::new(0.0, 1.0))
                        .unwrap();
                let envelope = match profile {
                    HarmonicWindowProfile::Uniform => half_beat.sin() / half_beat,
                    HarmonicWindowProfile::Triangular => {
                        ((0.5 * half_beat).sin() / (0.5 * half_beat)).powi(2)
                    }
                };
                let expected = (TAU * (start.fract() + 0.5 * duration)).cos() * envelope;
                close(actual.to_f64(), expected, 3e-15);
                let proper = window
                    .proper_product_kernel_hz(carrier + 1.0, carrier)
                    .unwrap()
                    .mean_real_product(C64::ONE, C64::ONE)
                    .unwrap();
                close(2.0 * proper, expected, 3e-15);
            }
        }
    }

    #[test]
    fn g3_hz_complete_cycles_and_subperiod_exposures_are_epoch_invariant() {
        let carrier = 1e15;
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            // Both profiles integrate an integer number of these cycles to
            // zero, regardless of the absolute source-clock epoch.
            let full = HarmonicWindow::with_profile(0.125, 2.0_f64.powi(30), profile).unwrap();
            assert_eq!(full.mean_phase_hz(carrier).unwrap(), C64::ZERO);
            for power in [30, 60, 100] {
                let late =
                    HarmonicWindow::with_profile(2.0_f64.powi(power), 2.0_f64.powi(-60), profile)
                        .unwrap();
                let early = HarmonicWindow::with_profile(0.0, late.duration_s(), profile).unwrap();
                // The integral-duration epoch adds only complete cycles to
                // both integer-Hz sources. Retain even subperiod interference.
                assert_eq!(
                    late.mean_phase_hz(carrier).unwrap(),
                    early.mean_phase_hz(carrier).unwrap()
                );
                let a = C64::new(0.7, -0.4);
                let b = C64::new(-0.2, 0.9);
                let observe = |w: HarmonicWindow| {
                    w.product_kernel_hz(carrier, carrier + 1.0)
                        .unwrap()
                        .mean_real_product_dd(a, b)
                        .unwrap()
                };
                assert_eq!(observe(late), observe(early));
            }
        }
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
