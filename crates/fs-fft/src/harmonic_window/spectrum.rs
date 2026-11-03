//! Exact spectral bookkeeping for quadratic real-field observables.
//!
//! Sum and difference frequencies remain explicit so a beat can drive a
//! downstream time-dependent model instead of being replaced by its DC power.
//! This is a finite line spectrum, not a sampled FFT, PSD or bandwidth model.

use super::*;

mod causal;
pub use causal::HarmonicCausalIntegral;
mod increment;
mod ripple;

/// Hard bound on generated tones BEFORE frequency coalescing.
pub const MAX_PRODUCT_SPECTRUM_TERMS: usize = 131_072;

/// Second-moment law of one bilinear contribution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicProductLaw {
    /// Product of deterministic real waveforms, or the expected product of
    /// waveforms sharing one zero-mean REAL coefficient with unit variance.
    Real,
    /// Both waveforms share one complex coefficient z with E[z]=E[z²]=0 and
    /// E[|z|²]=1. Only the difference-frequency contribution survives.
    ProperComplex,
}

/// One term in an explicitly authored bilinear field contraction.
///
/// Operands use peak amplitudes and exp(-i omega t). Signed frequencies work.
/// The caller must combine coherent fields at each frequency BEFORE forming
/// these products to retain weak residuals beside large cancelling sources.
/// Statistical factors contribute separately only under their declared joint
/// moments; unknown phases do not imply independence.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicProductInput {
    pub left_frequency_hz: f64,
    pub right_frequency_hz: f64,
    pub left_amplitude: [Dd; 2],
    pub right_amplitude: [Dd; 2],
    pub law: HarmonicProductLaw,
}

/// One contraction whose source frequencies are authored in radians/second.
///
/// This retains the original angular clock without first rounding omega/TAU
/// to one f64 Hz value. Amplitudes and the negative-time convention have the
/// same meaning as in [`HarmonicProductInput`].
#[derive(Clone, Copy, Debug)]
pub struct AngularHarmonicProductInput {
    pub left_angular_frequency_rad_s: f64,
    pub right_angular_frequency_rad_s: f64,
    pub left_amplitude: [Dd; 2],
    pub right_amplitude: [Dd; 2],
    pub law: HarmonicProductLaw,
}

/// One already-contracted term of a REAL observable using exp(-i omega t).
///
/// Spatial quadratures can contract original field products before emitting
/// this term. The peak amplitude already includes every physical weight and
/// any one-half product normalization. Zero angular frequency contributes
/// only its real part to DC; negative frequencies are conjugated. The input
/// DD pairs are normalized without first rounding their sums to f64.
///
/// This declaration does not infer source statistics, spatial quadrature
/// accuracy, units, dissipation, or a random realization from an expected
/// waveform. Those belong to the original physical contraction.
#[derive(Clone, Copy, Debug)]
pub struct AngularHarmonicSpectrumTerm {
    pub angular_frequency_rad_s: Dd,
    pub peak_amplitude: [Dd; 2],
}

/// Admission for already-contracted terms, before frequency coalescing.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSpectrumTermLimits {
    pub maximum_terms: usize,
    /// Logical vector payload at its admitted capacity. Caller inputs,
    /// iterator state, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

struct NativeHarmonicProductInput {
    left_frequency: f64,
    right_frequency: f64,
    left_amplitude: [Dd; 2],
    right_amplitude: [Dd; 2],
    law: HarmonicProductLaw,
}

impl From<HarmonicProductInput> for NativeHarmonicProductInput {
    fn from(value: HarmonicProductInput) -> Self {
        Self {
            left_frequency: value.left_frequency_hz,
            right_frequency: value.right_frequency_hz,
            left_amplitude: value.left_amplitude,
            right_amplitude: value.right_amplitude,
            law: value.law,
        }
    }
}

impl From<AngularHarmonicProductInput> for NativeHarmonicProductInput {
    fn from(value: AngularHarmonicProductInput) -> Self {
        Self {
            left_frequency: value.left_angular_frequency_rad_s,
            right_frequency: value.right_angular_frequency_rad_s,
            left_amplitude: value.left_amplitude,
            right_amplitude: value.right_amplitude,
            law: value.law,
        }
    }
}

/// Admission bounds including the temporary unsimplified line spectrum.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicProductSpectrumLimits {
    pub maximum_products: usize,
    pub maximum_generated_terms: usize,
    /// Logical retained vector payload at its admitted capacity. Inputs,
    /// iterator state, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// A positive-frequency peak phasor of the resulting REAL observable.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicProductTone {
    frequency_hz: Dd,
    angular_frequency_rad_s: Dd,
    amplitude: [Dd; 2],
    order: usize,
}

impl HarmonicProductTone {
    /// Preserve low frequency bits: a sum of two authored f64 carriers may
    /// not itself be exactly representable as one f64.
    pub const fn frequency_hz_dd(self) -> Dd {
        self.frequency_hz
    }
    /// Original angular frequency [rad/s] used by temporal kernels. Angular
    /// inputs retain their exact DD sum/difference without a round trip via Hz.
    pub const fn angular_frequency_dd(self) -> Dd {
        self.angular_frequency_rad_s
    }
    /// Real/imaginary peak amplitude in the product of the operand units.
    pub const fn peak_amplitude_dd(self) -> [Dd; 2] {
        self.amplitude
    }
}

/// A refused spectral preparation or observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProductSpectrumError {
    Window(WindowError),
    Budget,
    Allocation,
    Cancelled,
}

impl core::fmt::Display for ProductSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Window(error) => write!(f, "quadratic spectrum: {error}"),
            Self::Budget => f.write_str("quadratic spectrum work/storage budget exceeded"),
            Self::Allocation => f.write_str("quadratic spectrum allocation failed"),
            Self::Cancelled => f.write_str("quadratic spectrum cancelled"),
        }
    }
}

impl std::error::Error for ProductSpectrumError {}

impl From<WindowError> for ProductSpectrumError {
    fn from(value: WindowError) -> Self {
        Self::Window(value)
    }
}

/// Prepared DC and distinct positive-frequency terms of a quadratic quantity.
///
/// Real products retain both sums and differences; proper products retain
/// only differences. Negative output frequencies are conjugated exactly.
/// Equal DD frequencies in the authored unit coalesce in input order, with no
/// tolerance-based merging or amplitude cutoff. Preparation uses O(K log K) work and O(K)
/// storage for generated K terms; scalar/window queries are allocation-free
/// O(distinct output tones). No sampled time grid or aliasing is introduced.
///
/// This is an expected waveform when proper/real random factors were supplied,
/// not a random realization or a fourth-moment noise spectrum of the product.
/// Double-double algebra improves cancellation, without a relative-error
/// certificate near arbitrary nulls or increased source/propagator accuracy.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicProductSpectrum {
    dc: Dd,
    tones: Vec<HarmonicProductTone>,
}

impl PreparedHarmonicProductSpectrum {
    pub fn new(
        products: impl IntoIterator<Item = HarmonicProductInput>,
        limits: HarmonicProductSpectrumLimits,
    ) -> Result<Self, ProductSpectrumError> {
        Self::new_cancellable(products, limits, || false)
    }

    /// Iterator work is bounded even if it supplies more items than expected.
    /// Cancellation/refusal discards preparation without publishing a spectrum.
    pub fn new_cancellable(
        products: impl IntoIterator<Item = HarmonicProductInput>,
        limits: HarmonicProductSpectrumLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        Self::prepare(
            products.into_iter().map(NativeHarmonicProductInput::from),
            false,
            limits,
            cancelled,
        )
    }

    /// Prepare source-authored angular frequencies without converting each
    /// carrier to f64 Hz. Coalescing compares exact DD angular sums/differences.
    pub fn from_angular_frequencies(
        products: impl IntoIterator<Item = AngularHarmonicProductInput>,
        limits: HarmonicProductSpectrumLimits,
    ) -> Result<Self, ProductSpectrumError> {
        Self::from_angular_frequencies_cancellable(products, limits, || false)
    }

    pub fn from_angular_frequencies_cancellable(
        products: impl IntoIterator<Item = AngularHarmonicProductInput>,
        limits: HarmonicProductSpectrumLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        Self::prepare(
            products.into_iter().map(NativeHarmonicProductInput::from),
            true,
            limits,
            cancelled,
        )
    }

    /// Prepare already-contracted real-observable terms without expanding
    /// physical quadrature samples into separate retained tones.
    pub fn from_angular_terms(
        terms: impl IntoIterator<Item = AngularHarmonicSpectrumTerm>,
        limits: HarmonicSpectrumTermLimits,
    ) -> Result<Self, ProductSpectrumError> {
        Self::from_angular_terms_cancellable(terms, limits, || false)
    }

    pub fn from_angular_terms_cancellable(
        terms: impl IntoIterator<Item = AngularHarmonicSpectrumTerm>,
        limits: HarmonicSpectrumTermLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        let capacity = limits.maximum_terms.min(MAX_PRODUCT_SPECTRUM_TERMS);
        let bytes = capacity
            .checked_mul(core::mem::size_of::<HarmonicProductTone>())
            .ok_or(ProductSpectrumError::Budget)?;
        if bytes > limits.maximum_workspace_bytes {
            return Err(ProductSpectrumError::Budget);
        }
        let mut tones = Vec::new();
        tones
            .try_reserve_exact(capacity)
            .map_err(|_| ProductSpectrumError::Allocation)?;
        for (index, input) in terms.into_iter().enumerate() {
            spectrum_poll(&mut cancelled)?;
            if index >= capacity {
                return Err(ProductSpectrumError::Budget);
            }
            let frequency = input.angular_frequency_rad_s;
            if !frequency.hi.is_finite() || !frequency.lo.is_finite() {
                return Err(WindowError::InvalidFrequency.into());
            }
            let frequency = checked(Dd::from_pair(frequency.hi, frequency.lo))?;
            let mut amplitude = input.peak_amplitude;
            for value in &mut amplitude {
                check_amplitude(*value)?;
                *value = checked(Dd::from_pair(value.hi, value.lo))?;
            }
            tones.push(tone(frequency, amplitude, index, true)?);
        }
        Self::finish(tones, true, cancelled)
    }

    fn prepare(
        products: impl Iterator<Item = NativeHarmonicProductInput>,
        angular: bool,
        limits: HarmonicProductSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        let capacity = limits
            .maximum_generated_terms
            .min(MAX_PRODUCT_SPECTRUM_TERMS);
        let bytes = capacity
            .checked_mul(core::mem::size_of::<HarmonicProductTone>())
            .ok_or(ProductSpectrumError::Budget)?;
        if bytes > limits.maximum_workspace_bytes {
            return Err(ProductSpectrumError::Budget);
        }
        let mut tones = Vec::new();
        tones
            .try_reserve_exact(capacity)
            .map_err(|_| ProductSpectrumError::Allocation)?;
        for (index, input) in products.enumerate() {
            spectrum_poll(&mut cancelled)?;
            let required = match input.law {
                HarmonicProductLaw::Real => 2,
                HarmonicProductLaw::ProperComplex => 1,
            };
            if index >= limits.maximum_products || required > capacity - tones.len() {
                return Err(ProductSpectrumError::Budget);
            }
            if !input.left_frequency.is_finite() || !input.right_frequency.is_finite() {
                return Err(WindowError::InvalidFrequency.into());
            }
            let a = input.left_amplitude;
            let b = input.right_amplitude;
            for value in a.into_iter().chain(b) {
                check_amplitude(value)?;
            }
            let fa = Dd::from_f64(input.left_frequency);
            let fb = Dd::from_f64(input.right_frequency);
            let half = Dd::from_f64(0.5);
            let difference = [
                checked(half * checked(a[0] * b[0] + a[1] * b[1])?)?,
                checked(half * checked(a[1] * b[0] - a[0] * b[1])?)?,
            ];
            let order = tones.len();
            tones.push(tone(checked(fa - fb)?, difference, order, angular)?);
            if input.law == HarmonicProductLaw::Real {
                let sum = [
                    checked(half * checked(a[0] * b[0] - a[1] * b[1])?)?,
                    checked(half * checked(a[0] * b[1] + a[1] * b[0])?)?,
                ];
                let order = tones.len();
                tones.push(tone(checked(fa + fb)?, sum, order, angular)?);
            }
        }
        Self::finish(tones, angular, cancelled)
    }

    // One coalescing path serves field products and physical quadrature
    // contractions, preserving exact authored frequencies and input order.
    fn finish(
        mut tones: Vec<HarmonicProductTone>,
        angular: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        // The fixed generated-term ceiling bounds this synchronous sort.
        tones.sort_unstable_by(|a, b| {
            frequency_order(native_frequency(*a, angular), native_frequency(*b, angular))
                .then_with(|| a.order.cmp(&b.order))
        });
        spectrum_poll(&mut cancelled)?;
        let mut dc = Dd::ZERO;
        let mut written = 0usize;
        for index in 0..tones.len() {
            if index % 64 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let next = tones[index];
            if is_zero(native_frequency(next, angular)) {
                dc = checked(dc + next.amplitude[0])?;
            } else if written > 0
                && frequency_order(
                    native_frequency(tones[written - 1], angular),
                    native_frequency(next, angular),
                )
                .is_eq()
            {
                let previous = &mut tones[written - 1];
                previous.amplitude[0] = checked(previous.amplitude[0] + next.amplitude[0])?;
                previous.amplitude[1] = checked(previous.amplitude[1] + next.amplitude[1])?;
            } else {
                tones[written] = next;
                written += 1;
            }
        }
        tones.truncate(written);
        spectrum_poll(&mut cancelled)?;
        Ok(Self { dc, tones })
    }

    /// Constant part, in the product of the operand units.
    pub const fn dc_dd(&self) -> Dd {
        self.dc
    }

    pub fn tones(&self) -> &[HarmonicProductTone] {
        &self.tones
    }

    pub fn evaluate(&self, time_s: f64) -> Result<f64, ProductSpectrumError> {
        self.evaluate_cancellable(time_s, || false)
    }

    /// Evaluate the expected REAL quadratic observable at one source-clock time.
    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, ProductSpectrumError> {
        spectrum_result(self.evaluate_from_clock_dd_cancellable(Dd::from_f64(time_s), cancelled)?)
    }

    pub fn mean(&self, window: HarmonicWindow) -> Result<f64, ProductSpectrumError> {
        self.mean_cancellable(window, || false)
    }

    /// Integrate the complete output spectrum with the selected normalized
    /// detector profile, retaining low carrier-frequency and amplitude parts.
    pub fn mean_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        let mut sum = self.dc;
        for (index, tone) in self.tones.iter().enumerate() {
            if index % 64 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let omega = tone.angular_frequency_rad_s;
            let kernel = window.product_kernel_extended(omega, Dd::ZERO)?;
            sum = checked(
                sum + kernel.mean_real_product_parts_dd(tone.amplitude, [Dd::ONE, Dd::ZERO])?,
            )?;
        }
        spectrum_poll(&mut cancelled)?;
        spectrum_result(sum)
    }
}

fn tone(
    mut frequency: Dd,
    mut amplitude: [Dd; 2],
    order: usize,
    angular: bool,
) -> Result<HarmonicProductTone, ProductSpectrumError> {
    if frequency.hi < 0.0 || (frequency.hi == 0.0 && frequency.lo < 0.0) {
        frequency = -frequency;
        amplitude[1] = -amplitude[1];
    }
    if is_zero(frequency) {
        frequency = Dd::ZERO;
    }
    if frequency.lo == 0.0 {
        frequency.lo = 0.0;
    }
    let tau = Dd::from_f64(core::f64::consts::TAU);
    let (frequency_hz, angular_frequency_rad_s) = if angular {
        (checked(frequency / tau)?, frequency)
    } else {
        (frequency, checked(frequency * tau)?)
    };
    Ok(HarmonicProductTone {
        frequency_hz,
        angular_frequency_rad_s,
        amplitude,
        order,
    })
}

// Each preparation has one authored unit. Equality is established before
// its display conversion, preserving distinct low-frequency carrier beats.
fn native_frequency(tone: HarmonicProductTone, angular: bool) -> Dd {
    if angular {
        tone.angular_frequency_rad_s
    } else {
        tone.frequency_hz
    }
}

fn frequency_order(a: Dd, b: Dd) -> core::cmp::Ordering {
    a.hi.total_cmp(&b.hi).then_with(|| a.lo.total_cmp(&b.lo))
}

fn is_zero(value: Dd) -> bool {
    value.hi == 0.0 && value.lo == 0.0
}

fn spectrum_result(value: Dd) -> Result<f64, ProductSpectrumError> {
    let rounded = checked(value)?.to_f64();
    if rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(WindowError::Unrepresentable.into())
    }
}

fn spectrum_poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ProductSpectrumError> {
    if cancelled() {
        Err(ProductSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
