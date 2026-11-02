//! Deterministic means with explicit real or proper-complex fluctuations.
//!
//! Real source factors retain nonzero pseudocovariance and sum-frequency
//! interference. They can describe phase-locked amplitude noise, or a real
//! covariance square root of the stacked real/imaginary source components.
//! No distribution, random phase, source validation or stationarity is inferred.

use super::*;

/// Joint second-moment model of one zero-mean source factor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicFactorStatistics {
    /// Real coefficient xi, E[xi]=0 and E[xi²]=1. Its waveform keeps its
    /// authored phase, so both sum- and difference-frequency terms remain.
    Real,
    /// Complex coefficient z, E[z]=E[z²]=0 and E[|z|²]=1.
    /// Properness removes this factor's sum-frequency contribution.
    ProperComplex,
}

/// One coherent field per explicitly authored random coefficient.
///
/// All coefficients in a bank are jointly uncorrelated: real-real covariance
/// is identity, proper-proper covariance is identity and pseudocovariance zero,
/// and real-proper mixed moments vanish. Those assumptions, with zero means,
/// make factors additive in second moments without Gaussian assumptions.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicFieldFactor<'a, const D: usize> {
    /// Actual propagated peak fields at the same receiver/frame/source clock.
    pub lines: &'a [HarmonicFieldLine<D>],
    pub statistics: HarmonicFactorStatistics,
}

/// Preparation/application bounds including the deterministic mean.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicEnsembleBudget {
    /// Supplied fluctuation factors, including empty factors.
    pub max_factors: usize,
    /// Sum of line counts in the mean and every factor.
    pub max_total_lines: usize,
    /// Sum of SQUARED line counts per mean/factor bank, before grouping.
    pub max_line_pairs: usize,
}

/// Mean and uncentered second moments over the ensemble and physical window.
///
/// The velocity tensor includes deterministic motion; it is not the centered
/// covariance of fluctuations or uncertainty in a detector's estimate. Signed
/// intensity is the actual mean p*u, never a pressure-only plane-wave proxy.
#[derive(Clone, Debug)]
pub struct WindowedAcousticEnsemble<const D: usize> {
    pub window: HarmonicWindow,
    pub mean_pressure_pa: f64,
    pub mean_particle_velocity_m_s: [f64; D],
    pub mean_square_pressure_pa2: f64,
    pub velocity_second_moment_m2_s2: [[f64; D]; D],
    pub mean_intensity_w_m2: [f64; D],
}

impl<const D: usize> WindowedAcousticEnsemble<D> {
    #[must_use]
    pub fn pressure_rms_pa(&self) -> f64 {
        det::sqrt(self.mean_square_pressure_pa2)
    }
}

/// Prepared deterministic mean plus explicit real/proper source factors.
///
/// Each factor retains coherent phases across sources AND frequencies. A real
/// factor can represent phase-locked amplitude variability; independent real
/// factors for the stacked real/imaginary coefficients represent general
/// positive-semidefinite nonproper complex covariance by construction. Proper
/// factors remain useful for circular complex fluctuations without doubling
/// their rank. Inputs already carry statistical square-root weights.
///
/// Preparation is O(sum N_r²) and stores O(sum N_r) fields for fixed dimension.
/// Equal frequencies group in source order with retained DD components.
/// Each observation allocates nothing and uses O(sum M_r²) temporal kernels,
/// reused across the full p/u moment tensor. Resource refusal or cancellation
/// publishes nothing; repeated observations do not mutate the prepared state.
///
/// A single real factor with unit variance has the same SECOND moments as its
/// deterministic waveform, but zero ensemble mean. A single proper factor
/// generally has different second moments over a short window. Unknown source
/// phase is never permission to substitute either statistical assumption.
///
/// This owner adds no medium, source calibration, spatial propagation, PSD/Hz
/// scaling, sampling bandwidth or error enclosure. All input fields must
/// already refer to one receiver, Cartesian frame and physical source clock.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicEnsemble<const D: usize> {
    mean: PreparedHarmonicField<D>,
    factors: Vec<(HarmonicFactorStatistics, PreparedHarmonicField<D>)>,
}

impl<const D: usize> PreparedHarmonicEnsemble<D> {
    pub fn new(
        deterministic_lines: &[HarmonicFieldLine<D>],
        factors: &[HarmonicFieldFactor<'_, D>],
        budget: HarmonicEnsembleBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(deterministic_lines, factors, budget, || false)
    }

    pub fn new_cancellable(
        deterministic_lines: &[HarmonicFieldLine<D>],
        factors: &[HarmonicFieldFactor<'_, D>],
        budget: HarmonicEnsembleBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if factors.len() > budget.max_factors {
            return Err(invalid("harmonic ensemble exceeds its factor budget"));
        }
        let mut total_lines = deterministic_lines.len();
        let mut pairs = total_lines
            .checked_mul(total_lines)
            .ok_or_else(|| invalid("harmonic ensemble pair count overflows"))?;
        for factor in factors {
            poll(&mut cancelled)?;
            total_lines = total_lines
                .checked_add(factor.lines.len())
                .ok_or_else(|| invalid("harmonic ensemble line count overflows"))?;
            pairs = factor
                .lines
                .len()
                .checked_mul(factor.lines.len())
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(|| invalid("harmonic ensemble pair count overflows"))?;
        }
        if total_lines > budget.max_total_lines || pairs > budget.max_line_pairs {
            return Err(invalid("harmonic ensemble exceeds its line/work budget"));
        }
        let field_budget = HarmonicObservationBudget {
            max_lines: budget.max_total_lines,
            max_line_pairs: budget.max_line_pairs,
        };
        let mean = PreparedHarmonicField::new_cancellable(
            deterministic_lines,
            field_budget,
            &mut cancelled,
        )?;
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(factors.len())
            .map_err(|_| invalid("harmonic ensemble allocation exceeds available memory"))?;
        for factor in factors {
            prepared.push((
                factor.statistics,
                PreparedHarmonicField::new_cancellable(factor.lines, field_budget, &mut cancelled)?,
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            mean,
            factors: prepared,
        })
    }

    #[must_use]
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }

    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Includes full velocity second moments and signed acoustic energy flux.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        let mut sums = FieldMoments::new();
        self.mean.accumulate(
            window,
            FactorKind::Deterministic,
            true,
            true,
            &mut sums,
            &mut cancelled,
        )?;
        for (statistics, field) in &self.factors {
            let kind = match statistics {
                HarmonicFactorStatistics::Real => FactorKind::Deterministic,
                HarmonicFactorStatistics::ProperComplex => FactorKind::Proper,
            };
            field.accumulate(window, kind, false, true, &mut sums, &mut cancelled)?;
        }
        poll(&mut cancelled)?;
        sums.finish(window)
    }
}
