//! Time-dependent spectra of pressure/velocity products from actual fields.
//!
//! Prepared sum/difference tones expose acoustic beats for downstream
//! temporal models. The quantity is neither PSD/Hz nor an inferred heat load.

use super::*;
use fs_fft::harmonic_window::{
    HarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
    MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum, ProductSpectrumError,
};

/// Physical field coordinate in the original observation frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AcousticFieldComponent {
    Pressure,
    Velocity { axis: usize },
}

/// Expected real waveform of one pressure/velocity product.
///
/// Pressure/pressure has units Pa², velocity/velocity m²/s², and either
/// pressure/velocity ordering W/m². A pressure/velocity product is signed
/// acoustic energy flux, including temporary reactive storage exchange; it
/// is not a local heating rate or a plane-wave pressure-squared estimate.
///
/// A deterministic field returns its actual product waveform. For an ensemble,
/// this is the expected product under the authored joint second moments,
/// including any nonstationary beats. It is not a fluctuation realization or
/// a spectrum of fourth-moment fluctuations of acoustic intensity.
#[derive(Clone, Debug)]
pub struct AcousticProductSpectrum {
    left: AcousticFieldComponent,
    right: AcousticFieldComponent,
    spectrum: PreparedHarmonicProductSpectrum,
}

impl AcousticProductSpectrum {
    pub const fn components(&self) -> [AcousticFieldComponent; 2] {
        [self.left, self.right]
    }
    /// DC plus positive-frequency peak tones. Frequencies/amplitudes retain DD
    /// parts; sum frequencies must not be rounded back to an unrelated carrier.
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    pub fn evaluate(&self, time_s: f64) -> Result<f64, AeroacError> {
        self.spectrum.evaluate(time_s).map_err(spectrum_error)
    }
    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        self.spectrum
            .evaluate_cancellable(time_s, cancelled)
            .map_err(spectrum_error)
    }
    pub fn observe_window(&self, window: HarmonicWindow) -> Result<f64, AeroacError> {
        self.spectrum.mean(window).map_err(spectrum_error)
    }
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        self.spectrum
            .mean_cancellable(window, cancelled)
            .map_err(spectrum_error)
    }
}

impl<const D: usize> PreparedHarmonicField<D> {
    /// Prepare one physical quadratic observable from already grouped coherent
    /// fields. No propagation or field approximation is added.
    pub fn product_spectrum(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        self.product_spectrum_cancellable(left, right, limits, || false)
    }
    pub fn product_spectrum_cancellable(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        build(
            core::iter::once((HarmonicProductLaw::Real, self)),
            left,
            right,
            limits,
            cancelled,
        )
    }
}

impl<const D: usize> PreparedProperHarmonicField<D> {
    pub fn product_spectrum(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        self.product_spectrum_cancellable(left, right, limits, || false)
    }
    pub fn product_spectrum_cancellable(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        build(
            self.modes
                .iter()
                .map(|field| (HarmonicProductLaw::ProperComplex, field)),
            left,
            right,
            limits,
            cancelled,
        )
    }
}

impl<const D: usize> PreparedHarmonicEnsemble<D> {
    pub fn product_spectrum(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        self.product_spectrum_cancellable(left, right, limits, || false)
    }
    pub fn product_spectrum_cancellable(
        &self,
        left: AcousticFieldComponent,
        right: AcousticFieldComponent,
        limits: HarmonicProductSpectrumLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticProductSpectrum, AeroacError> {
        build(
            core::iter::once((HarmonicProductLaw::Real, &self.mean)).chain(
                self.factors.iter().map(|(law, field)| {
                    (
                        match law {
                            HarmonicFactorStatistics::Real => HarmonicProductLaw::Real,
                            HarmonicFactorStatistics::ProperComplex => {
                                HarmonicProductLaw::ProperComplex
                            }
                        },
                        field,
                    )
                }),
            ),
            left,
            right,
            limits,
            cancelled,
        )
    }
}

fn build<'a, const D: usize>(
    banks: impl Iterator<Item = (HarmonicProductLaw, &'a PreparedHarmonicField<D>)> + Clone,
    left: AcousticFieldComponent,
    right: AcousticFieldComponent,
    limits: HarmonicProductSpectrumLimits,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AcousticProductSpectrum, AeroacError> {
    poll(&mut cancelled)?;
    for component in [left, right] {
        if let AcousticFieldComponent::Velocity { axis } = component {
            if axis >= D {
                return Err(invalid(
                    "quadratic acoustic axis is outside the physical frame",
                ));
            }
        }
    }
    let mut pairs = 0usize;
    let mut terms = 0usize;
    for (law, field) in banks.clone() {
        poll(&mut cancelled)?;
        let count = field
            .lines
            .len()
            .checked_mul(field.lines.len())
            .ok_or_else(|| invalid("quadratic acoustic pair count overflows"))?;
        pairs = pairs
            .checked_add(count)
            .ok_or_else(|| invalid("quadratic acoustic work overflows"))?;
        terms = count
            .checked_mul(if law == HarmonicProductLaw::Real {
                2
            } else {
                1
            })
            .and_then(|n| terms.checked_add(n))
            .ok_or_else(|| invalid("quadratic acoustic spectrum count overflows"))?;
    }
    if pairs > limits.maximum_products
        || terms > limits.maximum_generated_terms
        || terms > MAX_PRODUCT_SPECTRUM_TERMS
    {
        return Err(invalid(
            "quadratic acoustic spectrum exceeds its work budget",
        ));
    }
    // A lazy contraction iterator avoids materializing the O(M²) product
    // inputs. Source fields were coalesced before any of these products, with
    // their low DD components retained through the shared spectrum engine.
    let products = banks.flat_map(move |(law, field)| {
        field.lines.iter().flat_map(move |a| {
            field.lines.iter().map(move |b| HarmonicProductInput {
                left_frequency_hz: a.frequency_hz,
                right_frequency_hz: b.frequency_hz,
                left_amplitude: parts(a, left),
                right_amplitude: parts(b, right),
                law,
            })
        })
    });
    let spectrum = PreparedHarmonicProductSpectrum::new_cancellable(
        products,
        HarmonicProductSpectrumLimits {
            maximum_products: pairs,
            maximum_generated_terms: terms,
            maximum_workspace_bytes: limits.maximum_workspace_bytes,
        },
        &mut cancelled,
    )
    .map_err(spectrum_error)?;
    poll(&mut cancelled)?;
    Ok(AcousticProductSpectrum {
        left,
        right,
        spectrum,
    })
}

fn parts<const D: usize>(line: &PreparedLine<D>, component: AcousticFieldComponent) -> [Dd; 2] {
    match component {
        AcousticFieldComponent::Pressure => line.pressure.parts(),
        AcousticFieldComponent::Velocity { axis } => line.velocity[axis].parts(),
    }
}

fn spectrum_error(error: ProductSpectrumError) -> AeroacError {
    match error {
        ProductSpectrumError::Cancelled => AeroacError::Cancelled {
            what: "quadratic acoustic spectrum",
        },
        _ => invalid("quadratic acoustic spectrum is outside its work/storage or numerical domain"),
    }
}
