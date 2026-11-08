//! Physical harmonic levels and distortion from a joint, known-frequency fit.
//!
//! [`HarmonicAnalysisPlan`] reuses `fs_fft::harmonic::HarmonicFitPlan`: DC and
//! every declared sine/cosine pair are fitted together, including on records
//! containing a noninteger number of periods. Frequencies are supplied in Hz,
//! never snapped to FFT bins. The phasor convention is `exp(-i omega t)`, with
//! phase referenced to the first sample.
//!
//! Distortion uses the fitted periodic model: `THD = sqrt(sum_{h>1}|P_h|²)/|P_1|`.
//! The common peak-to-RMS factor cancels. These are period-averaged harmonic
//! contributions, not a decomposition of the raw finite-record mean square;
//! finite records can have phase-dependent cross terms. Only explicitly
//! retained orders enter the sum. DC and the unexplained sampled residual are
//! reported separately, and no noise, alias, frequency-identification or
//! instrument-realism verdict is inferred. For inharmonic or decaying partials,
//! use the tone or damped-tone owners in `fs_fft::harmonic` instead.
//!
//! A caller-supplied fundamental amplitude floor prevents normalization by an
//! unresolved reference. It is an admission threshold in the input's units,
//! not an estimated uncertainty or an automatic noise floor. No harmonic is
//! clipped to that threshold, and zero levels are represented without a dB
//! floor. These physical descriptors are not perceptual or listening verdicts.

use fs_fft::C64;
pub use fs_fft::harmonic::{HarmonicFitBudget, HarmonicFitError};
use fs_fft::harmonic::{HarmonicFitPlan, HarmonicFitSummary};
use fs_math::det;

/// Whether the fitted fundamental can serve as the declared amplitude reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FundamentalStatus {
    /// Order one was not included in the fit; its amplitude is unknown.
    NotRetained,
    /// Both fitted fundamental quadratures are exactly zero.
    Zero,
    /// The nonzero peak is at or below the caller's resolution floor.
    BelowMinimum,
    /// The fitted peak exceeds the declared floor. This is not a noise test.
    Resolved,
}

/// A logarithmic amplitude ratio, with exact zero cases kept explicit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HarmonicLevelRatio {
    /// `20 log10(numerator / denominator)` for positive amplitudes.
    Decibels(f64),
    /// Zero numerator and positive reference: mathematically minus infinity.
    NegativeInfinity,
    /// Positive numerator and zero reference: mathematically plus infinity.
    PositiveInfinity,
    /// Both numerator and reference are zero: the ratio is undefined.
    Undefined,
}

/// One explicitly retained harmonic in the record's original units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicComponent {
    /// Positive integer multiple of the supplied fundamental.
    pub order: usize,
    /// Physical frequency in Hz, without a spectral-bin offset.
    pub frequency_hz: f64,
    /// Peak `exp(-i omega t)` phasor, referenced to the first sample.
    pub peak_phasor: C64,
    /// Level relative to order one, absent when that reference is unresolved.
    /// A retained exact-zero component has `NegativeInfinity`, not absence.
    pub relative_level: Option<HarmonicLevelRatio>,
}

/// Harmonic RMS ratios relative to the fundamental RMS, excluding DC/residual.
///
/// A ratio of 0.01 is 1%, not 0.01%. `None` means that no orders in that
/// category were retained; `Some(0.0)` means they were fitted and all were zero.
/// Even and odd ratios are contributions to THD: their squares sum to the
/// total's square when both categories are present.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicDistortion {
    /// Root-sum-square ratio over every retained order above one.
    pub total_ratio: Option<f64>,
    /// Root-sum-square ratio over retained even orders.
    pub even_ratio: Option<f64>,
    /// Root-sum-square ratio over retained odd orders above one.
    pub odd_ratio: Option<f64>,
    /// Largest retained order; no claim is made about higher frequencies.
    pub highest_order: usize,
    /// All orders 1 through `highest_order` were retained. Otherwise the ratios
    /// cover only the explicitly listed subset, not the complete finite band.
    pub complete_band: bool,
}

/// Prepared, immutable joint fit shared across records with one physical clock.
///
/// Storage and conditioning admission, preparation, and fit work belong to the
/// existing FFT harmonic owner. Analysis adds O(K) scalar descriptor work and
/// no heap allocation. Preparation and fixed-block analysis are synchronous;
/// this front adds no streaming or cancellation guarantee to that owner.
#[derive(Debug, Clone)]
pub struct HarmonicAnalysisPlan {
    fit: HarmonicFitPlan,
    minimum_fundamental_peak: f64,
}

impl HarmonicAnalysisPlan {
    /// Prepare positive, strictly increasing orders below Nyquist, with a
    /// finite nonnegative reference floor in the input samples' amplitude units.
    /// The supplied frequency is not estimated; a wrong frequency or omitted
    /// signal component can affect both fitted coefficients and the residual.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        fundamental_hz: f64,
        harmonics: &[usize],
        minimum_fundamental_peak: f64,
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if !(minimum_fundamental_peak.is_finite() && minimum_fundamental_peak >= 0.0) {
            return Err(HarmonicFitError::InvalidParameter {
                what: "fundamental amplitude floor must be finite and nonnegative",
            });
        }
        Ok(Self {
            fit: HarmonicFitPlan::new(sample_count, timestep_s, fundamental_hz, harmonics, budget)?,
            minimum_fundamental_peak,
        })
    }

    /// Shared numerical plan, including clock, orders, workspace size and
    /// reciprocal-condition diagnostic. No second harmonic estimator is used.
    #[must_use]
    pub const fn fit_plan(&self) -> &HarmonicFitPlan {
        &self.fit
    }

    /// Caller-declared minimum resolved fundamental peak in sample units.
    #[must_use]
    pub const fn minimum_fundamental_peak(&self) -> f64 {
        self.minimum_fundamental_peak
    }

    /// Fit without allocation using N scratch scalars and K output phasors.
    ///
    /// The returned read-only view keeps the phasor buffer and its plan paired.
    /// Drop the view before reusing that buffer for another record. Refusals
    /// preserve the output phasors; scratch may change, as in the shared fitter.
    pub fn analyze_into<'a>(
        &'a self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &'a mut [C64],
    ) -> Result<HarmonicAnalysis<'a>, HarmonicFitError> {
        let summary = self.fit.fit_into(samples, scratch, phasors)?;
        Ok(HarmonicAnalysis {
            plan: self,
            phasors,
            summary,
        })
    }
}

/// Fitted physical spectrum paired with its declared frequencies and clock.
#[derive(Debug)]
pub struct HarmonicAnalysis<'a> {
    plan: &'a HarmonicAnalysisPlan,
    phasors: &'a [C64],
    summary: HarmonicFitSummary,
}

impl HarmonicAnalysis<'_> {
    /// Retained positive orders in phasor order. Unlisted orders are unknown.
    #[must_use]
    pub fn harmonics(&self) -> &[usize] {
        self.plan.fit.harmonics()
    }

    /// Peak phasors in [`Self::harmonics`] order, in the input's original units.
    #[must_use]
    pub const fn peak_phasors(&self) -> &[C64] {
        self.phasors
    }

    /// Jointly fitted constant, not the arithmetic mean of a partial period.
    #[must_use]
    pub const fn dc(&self) -> f64 {
        self.summary.dc
    }

    /// RMS of the unexplained sampled record; not a noise or alias certificate.
    #[must_use]
    pub const fn residual_rms(&self) -> f64 {
        self.summary.residual_rms
    }

    /// Whether order one can serve as the caller-declared reference.
    #[must_use]
    pub fn fundamental_status(&self) -> FundamentalStatus {
        if self.harmonics()[0] != 1 {
            return FundamentalStatus::NotRetained;
        }
        let fundamental = ScaledSquares::of(self.phasors[0]);
        if fundamental.scale == 0.0 {
            FundamentalStatus::Zero
        } else if fundamental.at_most(self.plan.minimum_fundamental_peak) {
            FundamentalStatus::BelowMinimum
        } else {
            FundamentalStatus::Resolved
        }
    }

    /// Describe a retained harmonic. `None` means the order was not fitted;
    /// an unresolved fundamental only withholds its relative level.
    #[must_use]
    pub fn component(&self, order: usize) -> Option<HarmonicComponent> {
        let index = self.harmonics().binary_search(&order).ok()?;
        let peak_phasor = self.phasors[index];
        let relative_level =
            (self.fundamental_status() == FundamentalStatus::Resolved).then(|| {
                level_ratio(
                    ScaledSquares::of(peak_phasor).log_norm(),
                    ScaledSquares::of(self.phasors[0]).log_norm(),
                )
            });
        Some(HarmonicComponent {
            order,
            frequency_hz: self.plan.fit.fundamental_hz() * order as f64,
            peak_phasor,
            relative_level,
        })
    }

    /// Period-averaged THD and even/odd contributions over the retained orders.
    ///
    /// `Ok(None)` means the fundamental reference is unresolved. Missing
    /// categories remain absent. Ratios outside the nonzero finite f64 range
    /// refuse instead of turning a nonzero contribution into zero or infinity;
    /// individual logarithmic levels remain available through [`Self::component`].
    pub fn distortion(&self) -> Result<Option<HarmonicDistortion>, HarmonicFitError> {
        if self.fundamental_status() != FundamentalStatus::Resolved {
            return Ok(None);
        }
        let fundamental = ScaledSquares::of(self.phasors[0]);
        let mut total = ScaledSquares::default();
        let mut even = ScaledSquares::default();
        let mut odd = ScaledSquares::default();
        let mut have_even = false;
        let mut have_odd = false;
        for (&order, &phasor) in self.harmonics().iter().zip(self.phasors).skip(1) {
            total.add_phasor(phasor);
            if order % 2 == 0 {
                have_even = true;
                even.add_phasor(phasor);
            } else {
                have_odd = true;
                odd.add_phasor(phasor);
            }
        }
        let ratio = |present: bool, norm: ScaledSquares| {
            present.then(|| norm.ratio(fundamental)).transpose()
        };
        let highest_order = self.harmonics()[self.harmonics().len() - 1];
        Ok(Some(HarmonicDistortion {
            total_ratio: ratio(have_even || have_odd, total)?,
            even_ratio: ratio(have_even, even)?,
            odd_ratio: ratio(have_odd, odd)?,
            highest_order,
            complete_band: highest_order == self.harmonics().len(),
        }))
    }

    /// `L_h - (L_{h-1} + L_{h+1}) / 2` in dB for a positive even order h.
    ///
    /// The reference amplitude is the geometric mean of the adjacent odd
    /// amplitudes, so the result does not need a resolved fundamental. `None`
    /// means the order is not positive/even or one of the three orders was
    /// not retained. Exact zeros yield explicit infinite/undefined cases.
    /// This is the even-deficit measurement needed by the panpipe X-Pp1
    /// protocol; it does not supply that protocol's physical or listening proof.
    #[must_use]
    pub fn even_harmonic_deficit(&self, order: usize) -> Option<HarmonicLevelRatio> {
        if order == 0 || order % 2 != 0 {
            return None;
        }
        let log_at = |h| {
            let index = self.harmonics().binary_search(&h).ok()?;
            Some(ScaledSquares::of(self.phasors[index]).log_norm())
        };
        let center = log_at(order)?;
        let lower = log_at(order - 1)?;
        let upper = log_at(order.checked_add(1)?)?;
        let reference = lower.zip(upper).map(|(a, b)| 0.5 * (a + b));
        Some(level_ratio(center, reference))
    }
}

fn level_ratio(numerator: Option<f64>, denominator: Option<f64>) -> HarmonicLevelRatio {
    match (numerator, denominator) {
        (Some(a), Some(b)) => {
            HarmonicLevelRatio::Decibels(20.0 / core::f64::consts::LN_10 * (a - b))
        }
        (None, Some(_)) => HarmonicLevelRatio::NegativeInfinity,
        (Some(_), None) => HarmonicLevelRatio::PositiveInfinity,
        (None, None) => HarmonicLevelRatio::Undefined,
    }
}

// Keep separate scales for the reference and each harmonic category. Squaring
// a physical pressure directly would lose weak partials or overflow loud ones;
// scaling everything by the largest category would also lose a weak reference.
#[derive(Debug, Default, Clone, Copy)]
struct ScaledSquares {
    scale: f64,
    sum: f64,
}

impl ScaledSquares {
    fn of(phasor: C64) -> Self {
        let mut norm = Self::default();
        norm.add_phasor(phasor);
        norm
    }

    fn add_phasor(&mut self, phasor: C64) {
        for value in [phasor.re.abs(), phasor.im.abs()] {
            if value > self.scale {
                let ratio = self.scale / value;
                self.sum = 1.0 + self.sum * ratio * ratio;
                self.scale = value;
            } else if value > 0.0 {
                let ratio = value / self.scale;
                self.sum += ratio * ratio;
            }
        }
    }

    fn log_norm(self) -> Option<f64> {
        (self.scale > 0.0).then(|| det::ln(self.scale) + 0.5 * det::ln(self.sum))
    }

    fn at_most(self, threshold: f64) -> bool {
        self.scale <= threshold && self.scale / threshold * det::sqrt(self.sum) <= 1.0
    }

    fn ratio(self, reference: Self) -> Result<f64, HarmonicFitError> {
        if self.scale == 0.0 {
            return Ok(0.0);
        }
        let mut value = (self.scale / reference.scale) * det::sqrt(self.sum / reference.sum);
        if !(value.is_finite() && value > 0.0) {
            // A scale quotient can overflow even when the norm quotient fits.
            // Log reconstruction is the cold extreme-range fallback only.
            value = det::exp(self.log_norm().unwrap() - reference.log_norm().unwrap());
        }
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(HarmonicFitError::NonFinite {
                what: "nonzero harmonic distortion ratio is outside the finite f64 range",
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_zero_levels_have_no_invented_floor() {
        assert_eq!(
            level_ratio(None, Some(0.0)),
            HarmonicLevelRatio::NegativeInfinity
        );
        assert_eq!(
            level_ratio(Some(0.0), None),
            HarmonicLevelRatio::PositiveInfinity
        );
        assert_eq!(level_ratio(None, None), HarmonicLevelRatio::Undefined);
    }

    #[test]
    fn extreme_norm_ratios_do_not_square_away_the_reference() {
        let weak = ScaledSquares::of(C64::new(1e-200, 1e-200));
        let loud = ScaledSquares::of(C64::new(0.3, -0.4));
        let ratio = loud.ratio(weak).unwrap();
        assert!((ratio / (0.5 / (1e-200 * 2.0_f64.sqrt())) - 1.0).abs() < 3e-15);
        let subnormal = ScaledSquares::of(C64::new(f64::from_bits(1), 0.0));
        assert!(loud.ratio(subnormal).is_err());
        assert!(
            subnormal
                .ratio(ScaledSquares::of(C64::new(1e100, 0.0)))
                .is_err()
        );
        assert_eq!(ScaledSquares::default().ratio(weak), Ok(0.0));
    }
}
