//! Joint, known-frequency tone and harmonic estimation for sampled records.
//!
//! The model is `x[j] = dc + sum_h (a_h cos(w_h j dt) + b_h sin(w_h j dt))`.
//! All columns, including DC, are fitted together with the shared Householder
//! QR owner. Unlike picking a nearest FFT bin or independently projecting
//! each tone, this does not assume an integer number of periods or orthogonal
//! sampled tones. [`ToneFitPlan`] retains arbitrary physical frequencies,
//! including inharmonic partials and beating sources. [`HarmonicFitPlan`]
//! supplies integer multiples of one fundamental to the same numerical owner.
//! Neither plan estimates frequency or removes unmodeled tones. For a supplied
//! local pitch seed, [`refine_harmonic_frequency`] jointly refines the common
//! fundamental and amplitudes within explicit bounds and evaluation budgets.
//!
//! The returned peak phasor is `a_h + i b_h` for `exp(-i w_h t)`, referenced
//! to the first sample. It is twice the conjugate of the positive-frequency
//! normalized DFT coefficient for a bin-centered real sinusoid.

use crate::C64;
use fs_la::factor::{Qr, qr};
use fs_math::det;

#[path = "harmonic/frequency.rs"]
mod frequency;
pub use frequency::{
    HarmonicFrequencyBudget, HarmonicFrequencyError, HarmonicFrequencyFit,
    refine_harmonic_frequency,
};

/// Explicit storage and conditioning limits for a prepared harmonic fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicFitBudget {
    /// Conservative peak bytes for numeric arrays during preparation or one
    /// fit, including the QR input/copy, scratch, and harmonic coefficients.
    /// Allocator metadata and caller-owned input records are excluded.
    pub max_workspace_bytes: usize,
    /// Smallest admitted `1 / (||R||_1 ||R^-1||_1)` of the joint design's QR.
    /// This is a numerical conditioning diagnostic, not an uncertainty bound.
    pub min_reciprocal_condition: f64,
}

/// A refused fit never produces partial amplitude or phase results.
#[derive(Debug, Clone, PartialEq)]
pub enum HarmonicFitError {
    /// Sampling, harmonics, tolerances, or buffer lengths are invalid.
    InvalidParameter {
        /// The violated admission condition.
        what: &'static str,
    },
    /// The prepared design exceeds the caller's storage budget.
    BudgetExceeded {
        /// Conservative peak numeric-array storage.
        required_bytes: usize,
        /// Caller-supplied maximum.
        max_bytes: usize,
    },
    /// The sampled DC and harmonic columns cannot be reliably separated.
    IllConditioned {
        /// Measured reciprocal 1-norm condition of the triangular QR factor.
        reciprocal_condition: f64,
        /// Caller-supplied minimum.
        minimum: f64,
    },
    /// Input samples or a computed result are not representable finite values.
    NonFinite {
        /// The offending quantity.
        what: &'static str,
    },
}

impl core::fmt::Display for HarmonicFitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidParameter { what } => write!(f, "harmonic fit: {what}"),
            Self::BudgetExceeded {
                required_bytes,
                max_bytes,
            } => write!(
                f,
                "harmonic fit requires {required_bytes} bytes, budget is {max_bytes}"
            ),
            Self::IllConditioned {
                reciprocal_condition,
                minimum,
            } => write!(
                f,
                "harmonic fit reciprocal condition {reciprocal_condition:e} is below {minimum:e}"
            ),
            Self::NonFinite { what } => write!(f, "non-finite harmonic fit: {what}"),
        }
    }
}

impl std::error::Error for HarmonicFitError {}

fn invalid(what: &'static str) -> HarmonicFitError {
    HarmonicFitError::InvalidParameter { what }
}

/// Joint fit result in the input record's physical units.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicFit {
    /// Fitted constant term, not the arithmetic mean of a partial-period record.
    pub dc: f64,
    /// Peak `exp(-i omega t)` phasors in the plan's tone or harmonic order.
    pub phasors: Vec<C64>,
    /// RMS of the unexplained sampled record, with denominator sample count.
    /// This is not a noise, aliasing, or physical model-error certificate.
    pub residual_rms: f64,
}

/// Allocation-free fit diagnostics; phasors are written into caller storage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicFitSummary {
    /// Jointly fitted constant term in input units.
    pub dc: f64,
    /// RMS of the unexplained sampled record in input units.
    pub residual_rms: f64,
}

/// Immutable QR plan shared across force components, sensors, or successive
/// records with the same sampling clock and declared harmonic frequencies.
///
/// For N samples and K harmonics, preparation takes O(N K^2 + K^3) work and
/// O(N K) storage; each fit takes O(N K + K^2) work. The K^3 term measures
/// triangular conditioning. This small-known-spectrum path complements FFTs;
/// it is not a replacement for broadband O(N log N) spectral analysis.
#[derive(Debug, Clone)]
pub struct HarmonicFitPlan {
    fundamental_hz: f64,
    harmonics: Vec<usize>,
    design: FitDesign,
}

/// A joint DC/sinusoid fit at explicit, possibly inharmonic frequencies.
///
/// Frequencies come from the caller's physical model or independent frequency
/// estimate. Close tones are fitted simultaneously, preserving their relative
/// phase and finite-record cross terms; unresolved tones are refused at the
/// declared conditioning floor. There is no rounding to harmonics or FFT bins.
/// Preparation and fitting have the same costs as [`HarmonicFitPlan`].
#[derive(Debug, Clone)]
pub struct ToneFitPlan {
    frequencies_hz: Vec<f64>,
    design: FitDesign,
}

// One factorization, scaling, and solve owner for both frequency declarations.
#[derive(Debug, Clone)]
struct FitDesign {
    samples: usize,
    timestep_s: f64,
    tones: usize,
    factor: Qr,
    reciprocal_condition: f64,
    workspace_bytes: usize,
}

impl HarmonicFitPlan {
    /// Prepare the joint DC/cosine/sine design. Harmonics must be nonempty,
    /// positive, strictly increasing, and strictly below Nyquist. Record length
    /// need not be a power of two or contain a whole number of periods.
    ///
    /// The fundamental is supplied, never inferred or rounded to an FFT bin.
    /// Short/ill-conditioned records refuse at the caller's conditioning floor.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        fundamental_hz: f64,
        harmonics: &[usize],
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if !(timestep_s.is_finite()
            && timestep_s > 0.0
            && fundamental_hz.is_finite()
            && fundamental_hz > 0.0)
        {
            return Err(invalid(
                "timestep and fundamental must be positive and finite",
            ));
        }
        if harmonics.is_empty()
            || harmonics[0] == 0
            || harmonics.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "harmonics must be positive, nonempty and strictly increasing",
            ));
        }
        let cycles_per_sample = fundamental_hz * timestep_s;
        for &order in harmonics {
            if order as u128 > (1_u128 << 53) || !(fundamental_hz * order as f64).is_finite() {
                return Err(invalid(
                    "harmonic frequencies must be finite, positive and below Nyquist",
                ));
            }
        }
        let design = FitDesign::new(
            sample_count,
            timestep_s,
            harmonics.len(),
            |k| cycles_per_sample * harmonics[k] as f64,
            budget,
        )?;
        Ok(Self {
            fundamental_hz,
            harmonics: harmonics.to_vec(),
            design,
        })
    }

    /// Required record and scratch length.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.design.samples
    }

    /// Uniform physical sampling interval in seconds.
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.design.timestep_s
    }

    /// Caller-supplied fundamental frequency in hertz.
    #[must_use]
    pub const fn fundamental_hz(&self) -> f64 {
        self.fundamental_hz
    }

    /// Retained positive harmonic orders, matching the output phasor ordering.
    #[must_use]
    pub fn harmonics(&self) -> &[usize] {
        &self.harmonics
    }

    /// Measured reciprocal 1-norm condition of the triangular QR factor.
    #[must_use]
    pub const fn reciprocal_condition(&self) -> f64 {
        self.design.reciprocal_condition
    }

    /// Conservative numeric-array storage admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.design.workspace_bytes
    }

    /// Fit one record, allocating its scratch and returned phasors.
    pub fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        self.design.fit(samples)
    }

    /// Fit without allocation, reusing a length-N scratch slice and K phasors.
    /// On error, phasors remain unchanged; scratch may have been modified.
    pub fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        self.design.fit_into(samples, scratch, phasors)
    }
}

impl ToneFitPlan {
    /// Prepare a joint fit at positive, strictly increasing finite frequencies
    /// below Nyquist. The clock is uniform, but record length and tone spacing
    /// need not coincide with any FFT grid. The complete design, including DC,
    /// must pass the declared conditioning and storage limits before use.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        frequencies_hz: &[f64],
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if frequencies_hz.is_empty()
            || frequencies_hz.iter().any(|f| !f.is_finite() || *f <= 0.0)
            || frequencies_hz.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "tone frequencies must be finite, positive and strictly increasing",
            ));
        }
        let design = FitDesign::new(
            sample_count,
            timestep_s,
            frequencies_hz.len(),
            |k| frequencies_hz[k] * timestep_s,
            budget,
        )?;
        Ok(Self {
            frequencies_hz: frequencies_hz.to_vec(),
            design,
        })
    }

    /// Physical tone frequencies in the same order as the returned phasors.
    #[must_use]
    pub fn frequencies_hz(&self) -> &[f64] {
        &self.frequencies_hz
    }

    /// Required record and scratch length.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.design.samples
    }

    /// Uniform physical sampling interval in seconds.
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.design.timestep_s
    }

    /// Measured reciprocal 1-norm condition of the joint triangular factor.
    #[must_use]
    pub const fn reciprocal_condition(&self) -> f64 {
        self.design.reciprocal_condition
    }

    /// Conservative numeric-array storage admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.design.workspace_bytes
    }

    /// Fit one record, allocating scratch and peak phasors in input units.
    pub fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        self.design.fit(samples)
    }

    /// Fit with caller-owned scratch and output, without allocation. On error,
    /// phasors remain unchanged; scratch may have been modified.
    pub fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        self.design.fit_into(samples, scratch, phasors)
    }
}

impl FitDesign {
    fn new(
        sample_count: usize,
        timestep_s: f64,
        tones: usize,
        cycles_per_sample: impl Fn(usize) -> f64,
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if !timestep_s.is_finite() || timestep_s <= 0.0 {
            return Err(invalid("timestep must be positive and finite"));
        }
        if !(budget.min_reciprocal_condition.is_finite()
            && budget.min_reciprocal_condition > 0.0
            && budget.min_reciprocal_condition <= 1.0)
        {
            return Err(invalid("reciprocal-condition floor must be in (0, 1]"));
        }
        let columns = tones
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("sinusoid design size overflows"))?;
        if sample_count < columns {
            return Err(invalid("need at least 1 + 2 * tone_count samples"));
        }
        for k in 0..tones {
            let cycles = cycles_per_sample(k);
            if !(cycles > 0.0 && cycles < 0.5) {
                return Err(invalid(
                    "tone frequencies must be positive and strictly below Nyquist",
                ));
            }
        }
        // QR keeps a copy of the N*p design while constructing its factor.
        // The conservative sum also covers tau, condition scratch, one fit's
        // sample scratch/output, and retained tone frequencies or harmonic
        // indices (usize <= u64). Admission precedes every numeric allocation.
        let entries = sample_count
            .checked_mul(columns)
            .ok_or_else(|| invalid("harmonic design size overflows"))?;
        let workspace_bytes = fit_workspace_bytes(sample_count, tones)?;
        if workspace_bytes > budget.max_workspace_bytes {
            return Err(HarmonicFitError::BudgetExceeded {
                required_bytes: workspace_bytes,
                max_bytes: budget.max_workspace_bytes,
            });
        }
        let mut design = vec![0.0; entries];
        for (j, row) in design.chunks_exact_mut(columns).enumerate() {
            row[0] = 1.0;
            for k in 0..tones {
                let theta = core::f64::consts::TAU * (cycles_per_sample(k) * j as f64);
                let (s, c) = det::sin_cos(theta);
                row[2 * k + 1] = c;
                row[2 * k + 2] = s;
            }
        }
        let factor = qr(&design, sample_count, columns);
        let reciprocal_condition = triangular_reciprocal_condition(&factor, columns);
        if reciprocal_condition < budget.min_reciprocal_condition {
            return Err(HarmonicFitError::IllConditioned {
                reciprocal_condition,
                minimum: budget.min_reciprocal_condition,
            });
        }
        Ok(Self {
            samples: sample_count,
            timestep_s,
            tones,
            factor,
            reciprocal_condition,
            workspace_bytes,
        })
    }

    fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        if samples.len() != self.samples {
            return Err(invalid(
                "record length differs from the prepared sampling grid",
            ));
        }
        let mut scratch = vec![0.0; self.samples];
        let mut phasors = vec![C64::default(); self.tones];
        let summary = self.fit_into(samples, &mut scratch, &mut phasors)?;
        Ok(HarmonicFit {
            dc: summary.dc,
            phasors,
            residual_rms: summary.residual_rms,
        })
    }

    /// Fit without allocation, reusing a length-N scratch slice and K phasors.
    /// All inputs are admitted before output mutation. On any later error the
    /// scratch is unspecified; output phasors remain unchanged.
    fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        if samples.len() != self.samples
            || scratch.len() != self.samples
            || phasors.len() != self.tones
        {
            return Err(invalid(
                "record/scratch/phasor lengths differ from the prepared design",
            ));
        }
        if samples.iter().any(|x| !x.is_finite()) {
            return Err(HarmonicFitError::NonFinite {
                what: "input sample",
            });
        }
        // Subtract a representable central value before scaling. This retains
        // small harmonic variation on a large preload and avoids overflow for
        // records spanning both signs near f64::MAX. DC is still jointly fitted.
        let low = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let high = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let center = low.midpoint(high);
        let scale = samples
            .iter()
            .map(|&x| (x - center).abs())
            .fold(0.0, f64::max);
        if scale == 0.0 {
            phasors.fill(C64::default());
            return Ok(HarmonicFitSummary {
                dc: center,
                residual_rms: 0.0,
            });
        }
        for (out, &x) in scratch.iter_mut().zip(samples) {
            *out = (x - center) / scale;
        }
        self.factor.apply_qt(scratch);
        let columns = 2 * self.tones + 1;
        let residual_rms =
            det::sqrt(scratch[columns..].iter().map(|v| v * v).sum::<f64>() / self.samples as f64)
                * scale;
        for i in (0..columns).rev() {
            let mut x = scratch[i];
            for (j, &value) in scratch.iter().enumerate().take(columns).skip(i + 1) {
                x = (-self.factor.r(i, j)).mul_add(value, x);
            }
            scratch[i] = x / self.factor.r(i, i);
        }
        let dc = scratch[0].mul_add(scale, center);
        if !dc.is_finite()
            || !residual_rms.is_finite()
            || scratch[1..columns]
                .iter()
                .any(|&x| !(x * scale).is_finite())
        {
            return Err(HarmonicFitError::NonFinite {
                what: "fitted coefficient or residual",
            });
        }
        for (k, phasor) in phasors.iter_mut().enumerate() {
            *phasor = C64::new(scratch[2 * k + 1] * scale, scratch[2 * k + 2] * scale);
        }
        Ok(HarmonicFitSummary { dc, residual_rms })
    }
}

// Shared admission calculation for fixed-frequency plans and local refinement.
// Refinement adds its retained sample/derivative buffers before any allocation.
fn fit_workspace_bytes(samples: usize, tones: usize) -> Result<usize, HarmonicFitError> {
    let columns = tones.checked_mul(2).and_then(|n| n.checked_add(1));
    columns
        .and_then(|p| samples.checked_mul(p).and_then(|n| n.checked_mul(2)))
        .and_then(|n| n.checked_add(samples))
        .and_then(|n| columns?.checked_mul(3).and_then(|p| n.checked_add(p)))
        .and_then(|n| tones.checked_mul(3).and_then(|h| n.checked_add(h)))
        .and_then(|n| n.checked_mul(core::mem::size_of::<f64>()))
        .ok_or_else(|| invalid("harmonic workspace size overflows"))
}

// The complete triangular inverse norm costs O(p^3), with O(p) scratch.
// This avoids treating a diagonal-ratio heuristic as a condition estimate.
fn triangular_reciprocal_condition(factor: &Qr, columns: usize) -> f64 {
    let mut r_norm = 0.0_f64;
    let mut inverse_norm = 0.0_f64;
    let mut inverse_column = vec![0.0; columns];
    for j in 0..columns {
        let column_sum: f64 = (0..=j).map(|i| factor.r(i, j).abs()).sum();
        if !column_sum.is_finite() || factor.r(j, j) == 0.0 {
            return 0.0;
        }
        r_norm = r_norm.max(column_sum);
        inverse_column.fill(0.0);
        for i in (0..=j).rev() {
            let mut value = if i == j { 1.0 } else { 0.0 };
            for (k, &entry) in inverse_column.iter().enumerate().take(j + 1).skip(i + 1) {
                value = (-factor.r(i, k)).mul_add(entry, value);
            }
            inverse_column[i] = value / factor.r(i, i);
        }
        let sum: f64 = inverse_column[..=j].iter().map(|x| x.abs()).sum();
        if !sum.is_finite() {
            return 0.0;
        }
        inverse_norm = inverse_norm.max(sum);
    }
    (1.0 / r_norm) / inverse_norm
}
