//! Finite-duration observations of harmonic acoustic fields and source ensembles.
//!
//! All lines refer to one receiver, coordinate frame and physical time origin,
//! using peak phasors `Re(A exp(-i omega t))`. The caller establishes that clock
//! and geometry; neither unknown source phases nor stochastic coherence are
//! inferred here. Propagation delays must already be in the observed phasors.
//!
//! A short observation contains both difference-frequency beats and
//! sum-frequency terms. Adding `|p|²/2` across frequencies loses them. We
//! integrate the actual real pressure and velocity over an explicit rectangular
//! window through the shared `fs_fft::harmonic_window` kernel. No sampling rate,
//! aliasing, or carrier-period averaging is hidden in this operation.
//! [`PreparedProperHarmonicField`] separately accepts explicitly authored
//! zero-mean proper complex coherence factors. Its ensemble average removes
//! sum-frequency terms by the declared statistics, while cross-frequency beats
//! within each factor remain. It never treats unknown phase as independence.
//! [`PreparedHarmonicEnsemble`] also admits a deterministic mean and explicitly
//! real factors, retaining phase-locked amplitude fluctuations and their
//! nonzero pseudocovariance through the same field-moment implementation.

use crate::AeroacError;
use fs_fft::harmonic_window::{
    HarmonicProduct, HarmonicWindow, ProperHarmonicProduct, WindowError, angular_frequency,
    phase_hz_from_clock_dd,
};
use fs_math::{c64::C64, dd::Dd, det};

mod band;
mod quadratic;
mod sampled;
pub use fs_fft::harmonic_window::HarmonicProductSpectrumLimits;
pub use quadratic::{AcousticFieldComponent, AcousticProductSpectrum};
pub use sampled::{HarmonicFieldSample, HarmonicSampleBudget, HarmonicSampleClock};

mod ensemble;
pub use ensemble::{
    HarmonicEnsembleBudget, HarmonicFactorStatistics, HarmonicFieldFactor,
    PreparedHarmonicEnsemble, WindowedAcousticEnsemble,
};

/// Spatial correlations and physical readouts on one authored source clock.
pub mod array;

/// A spectral line of an already propagated acoustic field.
///
/// Pressure and velocity share one physical clock and the same peak-amplitude,
/// negative-time harmonic convention. Different frequencies must retain their
/// actual phases. [`PreparedHarmonicField`] interprets these as deterministic
/// amplitudes; [`PreparedProperHarmonicField`] interprets each factor's lines
/// as coefficients of an explicitly declared proper complex random variable.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicFieldLine<const D: usize> {
    /// Positive frequency [Hz]. Exactly equal frequencies are coalesced.
    pub frequency_hz: f64,
    /// Observed complex peak acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Observed complex peak particle velocity [m/s] in one Cartesian frame.
    pub particle_velocity_m_s: [C64; D],
}

/// An already propagated line retaining a source-authored angular frequency.
///
/// This has the same physical peak-phasor convention as [`HarmonicFieldLine`],
/// but no round trip through f64 hertz is applied to its original rad/s clock.
#[derive(Clone, Copy, Debug)]
pub struct AngularHarmonicFieldLine<const D: usize> {
    /// Positive original angular frequency [rad/s]. Exact equals coalesce.
    pub angular_frequency_rad_s: f64,
    /// Observed complex peak acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Observed complex peak particle velocity [m/s] in one Cartesian frame.
    pub particle_velocity_m_s: [C64; D],
}

/// Admission bounds for preparation and each finite-window observation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicObservationBudget {
    /// Maximum supplied line count, checked before allocation.
    pub max_lines: usize,
    /// Maximum square of the supplied line count, with overflow refused.
    /// This bounds coalescing and the fallback pairwise interference work.
    /// Short uniform bands instead use at most 64 tones and 17 temporal modes
    /// per field component, with fixed scratch and bounded projection work.
    pub max_line_pairs: usize,
}

/// Real-field averages on one explicit finite observation interval.
///
/// The intensity is `<p(t) u(t)>`, including interference. It is not a
/// plane-wave pressure-squared approximation, and Cartesian components may be
/// negative. A finite window does not define a separate reactive intensity.
/// These are floating-point evaluations, not certified error enclosures or
/// calibrated/source-validated SPL, and they preserve every input model limit.
#[derive(Clone, Debug)]
pub struct WindowedAcousticField<const D: usize> {
    /// Rectangular observer-time interval [s].
    pub window: HarmonicWindow,
    /// Mean real pressure [Pa]; finite windows need not have zero mean.
    pub mean_pressure_pa: f64,
    /// Mean square of the real pressure [Pa²], without subtracting its mean.
    pub mean_square_pressure_pa2: f64,
    /// Mean real particle velocity [m/s].
    pub mean_particle_velocity_m_s: [f64; D],
    /// Mean physical energy-flux vector `<p(t) u(t)>` [W/m²].
    pub mean_intensity_w_m2: [f64; D],
}

impl<const D: usize> WindowedAcousticField<D> {
    /// Pressure RMS on this interval [Pa], not a sound-power estimate.
    #[must_use]
    pub fn pressure_rms_pa(&self) -> f64 {
        det::sqrt(self.mean_square_pressure_pa2)
    }
}

/// Prepared deterministic harmonic field with equal frequencies combined.
///
/// Preparation is bounded O(N²) and allocates at most N spectral lines. Equal
/// frequencies are summed in input order using double-double accumulators
/// before any quadratic observation, retaining small cancellation residuals.
/// Frequencies that are merely close remain distinct so their beats survive.
/// Each window is allocation-free. Short uniform bands project pressure and
/// velocity through one positive temporal factor (at most 64 analytic or 32
/// real-waveform tones and 17 modes). Other windows retain O(M²) pair kernels
/// for M distinct frequencies, reused across pressure and velocity components.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicField<const D: usize> {
    lines: Vec<PreparedLine<D>>,
}

impl<const D: usize> PreparedHarmonicField<D> {
    /// Prepare lines at one receiver from known coherent source clocks.
    ///
    /// Empty input represents silence. The caller must not mix receivers,
    /// coordinate frames, amplitude conventions or unrelated random phases.
    /// Existing 2D/3D observers supply correctly normalized lines through their
    /// evaluate_harmonic_line methods; their source/propagation limits remain.
    pub fn new(
        lines: &[HarmonicFieldLine<D>],
        budget: HarmonicObservationBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(lines, budget, || false)
    }

    /// Admit and coalesce fields with bounded cancellation polls.
    /// No partially prepared field is published after cancellation or refusal.
    pub fn new_cancellable(
        lines: &[HarmonicFieldLine<D>],
        budget: HarmonicObservationBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            lines.iter().map(|line| {
                (
                    PreparedFrequency::Hertz(line.frequency_hz),
                    line.pressure_pa,
                    line.particle_velocity_m_s,
                )
            }),
            budget,
            cancelled,
        )
    }

    /// Preserve source-authored rad/s through samples, window moments and
    /// quadratic spectra. Frequencies are never reconstructed from rounded Hz.
    pub fn new_angular(
        lines: &[AngularHarmonicFieldLine<D>],
        budget: HarmonicObservationBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_angular_cancellable(lines, budget, || false)
    }

    /// Admit native angular lines with bounded cancellation and no partial result.
    pub fn new_angular_cancellable(
        lines: &[AngularHarmonicFieldLine<D>],
        budget: HarmonicObservationBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            lines.iter().map(|line| {
                (
                    PreparedFrequency::Angular(line.angular_frequency_rad_s),
                    line.pressure_pa,
                    line.particle_velocity_m_s,
                )
            }),
            budget,
            cancelled,
        )
    }

    fn prepare(
        lines: impl ExactSizeIterator<Item = (PreparedFrequency, C64, [C64; D])>,
        budget: HarmonicObservationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if lines.len() > budget.max_lines
            || lines
                .len()
                .checked_mul(lines.len())
                .is_none_or(|pairs| pairs > budget.max_line_pairs)
        {
            return Err(invalid("harmonic observation exceeds the line/work budget"));
        }
        let mut combined: Vec<PreparedLine<D>> = Vec::new();
        combined
            .try_reserve_exact(lines.len())
            .map_err(|_| invalid("harmonic field allocation exceeds available memory"))?;
        for (frequency, pressure_pa, particle_velocity_m_s) in lines {
            poll(&mut cancelled)?;
            let omega = frequency
                .angular()
                .map_err(|_| invalid("harmonic frequencies must be positive finite"))?;
            if !omega.hi.is_finite() || !omega.lo.is_finite() || omega.hi <= 0.0 {
                return Err(invalid("harmonic frequencies must be positive finite"));
            }
            if !finite(pressure_pa) || particle_velocity_m_s.iter().any(|value| !finite(*value)) {
                return Err(nonfinite("harmonic pressure or particle velocity"));
            }
            let mut group = None;
            for (index, previous) in combined.iter().enumerate() {
                poll(&mut cancelled)?;
                if previous.frequency == frequency {
                    group = Some(index);
                    break;
                }
            }
            let index = group.unwrap_or(combined.len());
            if index == combined.len() {
                combined.push(PreparedLine {
                    frequency,
                    pressure: ComplexSum::default(),
                    velocity: [ComplexSum::default(); D],
                });
            }
            combined[index].pressure.add(pressure_pa);
            for (sum, value) in combined[index]
                .velocity
                .iter_mut()
                .zip(particle_velocity_m_s)
            {
                sum.add(value);
            }
        }
        for line in &combined {
            poll(&mut cancelled)?;
            if !line.pressure.is_finite() || line.velocity.iter().any(|value| !value.is_finite()) {
                return Err(nonfinite("coalesced harmonic field"));
            }
        }
        poll(&mut cancelled)?;
        Ok(Self { lines: combined })
    }

    /// Number of distinct physical frequencies; close frequencies stay separate.
    #[must_use]
    pub fn frequency_count(&self) -> usize {
        self.lines.len()
    }

    /// Integrate actual real pressure and energy flux on the source clock.
    ///
    /// Both sum- and difference-frequency contributions survive. Short uniform
    /// bands use positive temporal field factors before forming any moments.
    /// Grouped fields, original frequencies and any Hz conversion remain in
    /// extended precision. This is not a relative certificate at a null.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticField<D>, AeroacError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Allocation-free observation with bounded mode/pair cancellation polls.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticField<D>, AeroacError> {
        let mut sums = FieldMoments::new();
        self.accumulate(
            window,
            FactorKind::Deterministic,
            true,
            false,
            &mut sums,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let value = sums.finish(window)?;
        Ok(WindowedAcousticField {
            window,
            mean_pressure_pa: value.mean_pressure_pa,
            mean_square_pressure_pa2: value.mean_square_pressure_pa2,
            mean_particle_velocity_m_s: value.mean_particle_velocity_m_s,
            mean_intensity_w_m2: value.mean_intensity_w_m2,
        })
    }

    /// Also retain the full uncentered real velocity second-moment tensor.
    pub fn observe_window_moments(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        self.observe_window_moments_cancellable(window, || false)
    }

    /// Full deterministic moments with atomic cancellation and no allocation.
    pub fn observe_window_moments_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        let mut sums = FieldMoments::new();
        self.accumulate(
            window,
            FactorKind::Deterministic,
            true,
            true,
            &mut sums,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        sums.finish(window)
    }

    fn accumulate(
        &self,
        window: HarmonicWindow,
        kind: FactorKind,
        include_mean: bool,
        full_velocity: bool,
        sums: &mut FieldMoments<D>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(cancelled)?;
        if band::accumulate(
            self,
            window,
            kind,
            include_mean,
            full_velocity,
            sums,
            cancelled,
        )? {
            return Ok(());
        }
        if include_mean {
            for line in &self.lines {
                poll(cancelled)?;
                // A unit DC field integrates each retained component without
                // rounding a tiny mean-phase decrement to C64 first.
                let kernel = line.frequency.mean_kernel(window).map_err(window_error)?;
                let unit = [Dd::ONE, Dd::ZERO];
                sums.pressure_mean = sums.pressure_mean
                    + kernel
                        .mean_real_product_parts_dd(line.pressure.parts(), unit)
                        .map_err(window_error)?;
                for (sum, value) in sums.velocity_mean.iter_mut().zip(line.velocity) {
                    *sum = *sum
                        + kernel
                            .mean_real_product_parts_dd(value.parts(), unit)
                            .map_err(window_error)?;
                }
            }
        }
        for (i, a) in self.lines.iter().enumerate() {
            for (j, b) in self.lines[..=i].iter().enumerate() {
                poll(cancelled)?;
                let kernel = match kind {
                    FactorKind::Deterministic => MomentKernel::Deterministic(
                        a.frequency
                            .product_kernel(b.frequency, window)
                            .map_err(window_error)?,
                    ),
                    FactorKind::Proper => MomentKernel::Proper(
                        a.frequency
                            .proper_kernel(b.frequency, window)
                            .map_err(window_error)?,
                    ),
                };
                let pp = kernel.product(a.pressure, b.pressure)?;
                sums.pressure_square = sums.pressure_square + pp;
                if i != j {
                    sums.pressure_square = sums.pressure_square + pp;
                }
                for (axis, sum) in sums.intensity.iter_mut().enumerate() {
                    *sum = *sum + kernel.product(a.pressure, b.velocity[axis])?;
                    if i != j {
                        *sum = *sum + kernel.product(a.velocity[axis], b.pressure)?;
                    }
                }
                if full_velocity {
                    for (row, values) in sums.velocity.iter_mut().enumerate() {
                        for (column, sum) in values.iter_mut().enumerate().skip(row) {
                            *sum = *sum + kernel.product(a.velocity[row], b.velocity[column])?;
                            if i != j {
                                *sum =
                                    *sum + kernel.product(a.velocity[column], b.velocity[row])?;
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)
    }
}

/// Bounds on an explicitly supplied bank of proper complex coherence factors.
#[derive(Clone, Copy, Debug)]
pub struct ProperHarmonicObservationBudget {
    /// Maximum factor count, including factors with no lines.
    pub max_modes: usize,
    /// Maximum sum of supplied line counts across all factors.
    pub max_total_lines: usize,
    /// Maximum sum of squared supplied line counts, checked for overflow.
    /// Bounds preparation and fallback pairwise observation before allocation.
    /// Short uniform bands have the fixed modal bound of HarmonicObservationBudget.
    pub max_line_pairs: usize,
}

/// Ensemble and rectangular-window second moments of a zero-mean proper field.
///
/// Both ensemble-mean pressure and ensemble-mean velocity are zero at every
/// time. The entries below average the actual real-field products, not products
/// of window-mean fields. They are neither PSD/Hz nor statistical confidence
/// intervals, and retain every original source/propagation model limitation.
#[derive(Clone, Debug)]
pub struct WindowedAcousticStatistics<const D: usize> {
    /// Rectangular receiver-time interval in seconds on the common source clock.
    pub window: HarmonicWindow,
    /// `<E[p(t)²]>` [Pa²], including authored cross-frequency interference.
    pub mean_square_pressure_pa2: f64,
    /// `<E[u_i(t) u_j(t)]>` [m²/s²], a real symmetric velocity covariance.
    pub velocity_covariance_m2_s2: [[f64; D]; D],
    /// `<E[p(t) u(t)]>` [W/m²]; Cartesian components may be negative.
    pub mean_intensity_w_m2: [f64; D],
}

impl<const D: usize> WindowedAcousticStatistics<D> {
    /// Root mean square over the ensemble and interval [Pa]. This is not the
    /// expectation of each realization's RMS, nor a sound-power estimate.
    #[must_use]
    pub fn pressure_rms_pa(&self) -> f64 {
        det::sqrt(self.mean_square_pressure_pa2)
    }
}

/// Prepared, explicitly authored proper complex source-coherence factors.
///
/// Factor r's field is `sum_j A[r,j] exp(-i omega[r,j] t)`. The physical field
/// is the real part of the sum of those factors multiplied by random variables
/// xi_r with `E[xi_r]=0`, `E[xi_r conj(xi_s)]=delta_rs` and `E[xi_r xi_s]=0`.
/// Thus the joint frequency/source covariance is positive semidefinite by
/// construction; no dense covariance, invented phase, or positivity repair is
/// needed. The second moments, not a Gaussian distribution, are the assumption.
///
/// Each factor coherently combines sources AND frequencies. Distinct factors
/// have no cross terms. A one-factor bank is not a deterministic phased tone:
/// properness removes its sum-frequency terms even during a subperiod window.
/// Cross-frequency coherence may produce time-dependent ensemble power; no
/// stationarity is claimed. Factors use peak units: a one-sided PSD bin of width
/// df needs coefficients whose covariance is `2 df S`, not the unscaled PSD.
///
/// Use the existing 2D/3D observers to propagate each source factor at each
/// frequency before preparation, retaining the same receiver, material, frame
/// and source clock. Equal frequencies coalesce within each factor before
/// squaring; merely close frequencies remain distinct. A deterministic mean
/// field can be observed with [`PreparedHarmonicField`]: its mean-square and
/// intensity add to these zero-mean fluctuation contributions without cross
/// terms. Nonzero pseudocovariance is supported separately by explicitly real
/// factors in [`PreparedHarmonicEnsemble`]; it is never inferred from these modes.
#[derive(Clone, Debug)]
pub struct PreparedProperHarmonicField<const D: usize> {
    modes: Vec<PreparedHarmonicField<D>>,
}

impl<const D: usize> PreparedProperHarmonicField<D> {
    /// Admit factors and total work. Empty factors count against max_modes.
    pub fn new(
        modes: &[&[HarmonicFieldLine<D>]],
        budget: ProperHarmonicObservationBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(modes, budget, || false)
    }

    /// Cancellable preparation with no partially prepared factor bank.
    pub fn new_cancellable(
        modes: &[&[HarmonicFieldLine<D>]],
        budget: ProperHarmonicObservationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if modes.len() > budget.max_modes {
            return Err(invalid("proper harmonic field exceeds the mode budget"));
        }
        let mut lines = 0_usize;
        let mut pairs = 0_usize;
        for mode in modes {
            poll(&mut cancelled)?;
            lines = lines
                .checked_add(mode.len())
                .ok_or_else(|| invalid("proper harmonic line count overflows"))?;
            pairs = mode
                .len()
                .checked_mul(mode.len())
                .and_then(|count| pairs.checked_add(count))
                .ok_or_else(|| invalid("proper harmonic pair count overflows"))?;
        }
        if lines > budget.max_total_lines || pairs > budget.max_line_pairs {
            return Err(invalid(
                "proper harmonic field exceeds the line/work budget",
            ));
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(modes.len())
            .map_err(|_| invalid("proper harmonic mode allocation exceeds available memory"))?;
        for mode in modes {
            prepared.push(PreparedHarmonicField::new_cancellable(
                mode,
                HarmonicObservationBudget {
                    max_lines: budget.max_total_lines,
                    max_line_pairs: budget.max_line_pairs,
                },
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(Self { modes: prepared })
    }

    /// Number of authored uncorrelated factors, including empty ones.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }

    /// Allocate nothing; reuse one temporal factor or each pair kernel across
    /// all field moments. Proper factors retain the real-peak factor of one half.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticStatistics<D>, AeroacError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Bounded cancellation polls without partial publication or owner mutation.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticStatistics<D>, AeroacError> {
        let mut sums = FieldMoments::new();
        poll(&mut cancelled)?;
        for mode in &self.modes {
            mode.accumulate(
                window,
                FactorKind::Proper,
                false,
                true,
                &mut sums,
                &mut cancelled,
            )?;
        }
        poll(&mut cancelled)?;
        let value = sums.finish(window)?;
        Ok(WindowedAcousticStatistics {
            window,
            mean_square_pressure_pa2: value.mean_square_pressure_pa2,
            velocity_covariance_m2_s2: value.velocity_second_moment_m2_s2,
            mean_intensity_w_m2: value.mean_intensity_w_m2,
        })
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct ComplexSum {
    re: Dd,
    im: Dd,
}

impl ComplexSum {
    fn add(&mut self, value: C64) {
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
    }

    fn parts(self) -> [Dd; 2] {
        [self.re, self.im]
    }

    fn is_finite(self) -> bool {
        [self.re.hi, self.re.lo, self.im.hi, self.im.lo]
            .iter()
            .all(|x| x.is_finite())
    }
}

#[derive(Clone, Debug)]
struct PreparedLine<const D: usize> {
    frequency: PreparedFrequency,
    pressure: ComplexSum,
    velocity: [ComplexSum; D],
}

// Preserve the unit in which the source clock was authored. The Hertz branch
// delegates to the existing Hz kernels, retaining their conversion semantics.
#[derive(Clone, Copy, Debug, PartialEq)]
enum PreparedFrequency {
    Hertz(f64),
    Angular(f64),
}

impl PreparedFrequency {
    fn angular(self) -> Result<Dd, WindowError> {
        match self {
            Self::Hertz(hz) => angular_frequency(hz),
            Self::Angular(omega) if omega.is_finite() => Ok(Dd::from_f64(omega)),
            Self::Angular(_) => Err(WindowError::InvalidFrequency),
        }
    }

    fn native(self) -> f64 {
        match self {
            Self::Hertz(value) | Self::Angular(value) => value,
        }
    }

    fn is_angular(self) -> bool {
        matches!(self, Self::Angular(_))
    }

    fn phase(self, time: Dd) -> Result<[Dd; 2], WindowError> {
        if let Self::Hertz(hz) = self {
            return phase_hz_from_clock_dd(hz, time);
        }
        if !time.hi.is_finite() || !time.lo.is_finite() {
            return Err(WindowError::InvalidWindow);
        }
        let angle = self.angular()? * Dd::from_pair(time.hi, time.lo);
        if !angle.hi.is_finite() || !angle.lo.is_finite() {
            return Err(WindowError::Unrepresentable);
        }
        let (sin, cos) = det::sin_cos_dd(angle);
        Ok([cos, -sin])
    }

    fn below_nyquist(self, sample_rate_hz: u32) -> Result<bool, WindowError> {
        let limit_hz = 0.5 * f64::from(sample_rate_hz);
        match self {
            Self::Hertz(hz) => Ok(hz < limit_hz),
            Self::Angular(_) => Ok(self.angular()?.lt(angular_frequency(limit_hz)?)),
        }
    }

    fn mean_kernel(self, window: HarmonicWindow) -> Result<HarmonicProduct, WindowError> {
        match self {
            Self::Hertz(hz) => window.product_kernel_hz(hz, 0.0),
            Self::Angular(omega) => window.product_kernel(omega, 0.0),
        }
    }

    fn product_kernel(
        self,
        other: Self,
        window: HarmonicWindow,
    ) -> Result<HarmonicProduct, WindowError> {
        match (self, other) {
            (Self::Hertz(a), Self::Hertz(b)) => window.product_kernel_hz(a, b),
            _ => window.product_kernel_dd(self.angular()?, other.angular()?),
        }
    }

    fn proper_kernel(
        self,
        other: Self,
        window: HarmonicWindow,
    ) -> Result<ProperHarmonicProduct, WindowError> {
        match (self, other) {
            (Self::Hertz(a), Self::Hertz(b)) => window.proper_product_kernel_hz(a, b),
            _ => window.proper_product_kernel_dd(self.angular()?, other.angular()?),
        }
    }
}

#[derive(Clone, Copy)]
enum FactorKind {
    Deterministic,
    Proper,
}

#[derive(Clone, Copy)]
enum MomentKernel {
    Deterministic(HarmonicProduct),
    Proper(ProperHarmonicProduct),
}

impl MomentKernel {
    fn product(self, a: ComplexSum, b: ComplexSum) -> Result<Dd, AeroacError> {
        match self {
            Self::Deterministic(kernel) => kernel.mean_real_product_parts_dd(a.parts(), b.parts()),
            Self::Proper(kernel) => kernel.mean_real_product_parts_dd(a.parts(), b.parts()),
        }
        .map_err(window_error)
    }
}

struct FieldMoments<const D: usize> {
    pressure_mean: Dd,
    velocity_mean: [Dd; D],
    pressure_square: Dd,
    velocity: [[Dd; D]; D],
    intensity: [Dd; D],
}

impl<const D: usize> FieldMoments<D> {
    fn new() -> Self {
        Self {
            pressure_mean: Dd::ZERO,
            velocity_mean: [Dd::ZERO; D],
            pressure_square: Dd::ZERO,
            velocity: [[Dd::ZERO; D]; D],
            intensity: [Dd::ZERO; D],
        }
    }

    fn finish(self, window: HarmonicWindow) -> Result<WindowedAcousticEnsemble<D>, AeroacError> {
        let value = WindowedAcousticEnsemble {
            window,
            mean_pressure_pa: self.pressure_mean.to_f64(),
            mean_particle_velocity_m_s: self.velocity_mean.map(Dd::to_f64),
            mean_square_pressure_pa2: self.pressure_square.to_f64(),
            velocity_second_moment_m2_s2: core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    self.velocity[row.min(column)][row.max(column)].to_f64()
                })
            }),
            mean_intensity_w_m2: self.intensity.map(Dd::to_f64),
        };
        if [value.mean_pressure_pa, value.mean_square_pressure_pa2]
            .iter()
            .chain(&value.mean_particle_velocity_m_s)
            .chain(&value.mean_intensity_w_m2)
            .chain(value.velocity_second_moment_m2_s2.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(nonfinite("finite-window acoustic field moments"));
        }
        if value.mean_square_pressure_pa2 < 0.0
            || (0..D).any(|axis| value.velocity_second_moment_m2_s2[axis][axis] < 0.0)
        {
            return Err(invalid(
                "finite-window acoustic self moment lost numerical positivity",
            ));
        }
        Ok(value)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "harmonic acoustic observation",
        })
    } else {
        Ok(())
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

fn nonfinite(what: &'static str) -> AeroacError {
    AeroacError::NonFinite { what }
}

fn window_error(_: WindowError) -> AeroacError {
    invalid("harmonic window integral is not representable")
}
