//! Finite-duration observations of harmonic acoustic fields and proper ensembles.
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

use crate::AeroacError;
use fs_fft::harmonic_window::{HarmonicWindow, WindowError};
use fs_math::{c64::C64, dd::Dd, det};

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

/// Admission bounds for preparation and each finite-window observation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicObservationBudget {
    /// Maximum supplied line count, checked before allocation.
    pub max_lines: usize,
    /// Maximum square of the supplied line count, with overflow refused.
    /// This bounds coalescing and all pairwise interference work.
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
/// Each window is allocation-free O(M²) for M distinct frequencies, caching
/// each frequency-pair integral across pressure and velocity components.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicField<const D: usize> {
    lines: Vec<HarmonicFieldLine<D>>,
}

impl<const D: usize> PreparedHarmonicField<D> {
    /// Prepare lines at one receiver from known coherent source clocks.
    ///
    /// Empty input represents silence. The caller must not mix receivers,
    /// coordinate frames, amplitude conventions or unrelated random phases.
    /// Existing 2D/3D observers supply correctly normalized lines through their
    /// `evaluate_harmonic_line` methods; the original source authority remains.
    #[allow(clippy::float_cmp)] // Only exact frequency equality permits coherent coalescing.
    pub fn new(
        lines: &[HarmonicFieldLine<D>],
        budget: HarmonicObservationBudget,
    ) -> Result<Self, AeroacError> {
        if lines.len() > budget.max_lines
            || lines
                .len()
                .checked_mul(lines.len())
                .is_none_or(|pairs| pairs > budget.max_line_pairs)
        {
            return Err(invalid("harmonic observation exceeds the line/work budget"));
        }
        for line in lines {
            let omega = core::f64::consts::TAU * line.frequency_hz;
            if !omega.is_finite() || omega <= 0.0 {
                return Err(invalid("harmonic frequencies must be positive finite"));
            }
            if !finite(line.pressure_pa) || line.particle_velocity_m_s.iter().any(|v| !finite(*v)) {
                return Err(nonfinite("harmonic pressure or particle velocity"));
            }
        }
        let mut combined = Vec::new();
        combined
            .try_reserve_exact(lines.len())
            .map_err(|_| invalid("harmonic field allocation exceeds available memory"))?;
        for (i, line) in lines.iter().enumerate() {
            if lines[..i]
                .iter()
                .any(|previous| previous.frequency_hz == line.frequency_hz)
            {
                continue;
            }
            let mut pressure = ComplexSum::default();
            let mut velocity = [ComplexSum::default(); D];
            for same in lines[i..]
                .iter()
                .filter(|same| same.frequency_hz == line.frequency_hz)
            {
                pressure.add(same.pressure_pa);
                for (sum, value) in velocity.iter_mut().zip(same.particle_velocity_m_s) {
                    sum.add(value);
                }
            }
            let pressure_pa = pressure.value();
            let particle_velocity_m_s = velocity.map(ComplexSum::value);
            if !finite(pressure_pa) || particle_velocity_m_s.iter().any(|v| !finite(*v)) {
                return Err(nonfinite("coalesced harmonic field"));
            }
            combined.push(HarmonicFieldLine {
                frequency_hz: line.frequency_hz,
                pressure_pa,
                particle_velocity_m_s,
            });
        }
        Ok(Self { lines: combined })
    }

    /// Number of distinct physical frequencies; close frequencies stay separate.
    #[must_use]
    pub fn frequency_count(&self) -> usize {
        self.lines.len()
    }

    /// Integrate pressure and energy flux over a finite observer-time window.
    ///
    /// Both sum- and difference-frequency contributions are retained, including
    /// for a window shorter than one period. Severe cancellation between close
    /// distinct frequencies is still floating-point limited; a negative computed
    /// pressure square refuses rather than being silently clamped to zero.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticField<D>, AeroacError> {
        let mut pressure_mean = Dd::from_f64(0.0);
        let mut pressure_square = Dd::from_f64(0.0);
        let mut velocity_mean = [Dd::from_f64(0.0); D];
        let mut intensity = [Dd::from_f64(0.0); D];
        for (i, a) in self.lines.iter().enumerate() {
            let omega_a = core::f64::consts::TAU * a.frequency_hz;
            let phase = window.mean_phase(omega_a).map_err(window_error)?;
            pressure_mean = pressure_mean + Dd::from_f64((a.pressure_pa * phase).re);
            for (sum, value) in velocity_mean.iter_mut().zip(a.particle_velocity_m_s) {
                *sum = *sum + Dd::from_f64((value * phase).re);
            }
            for (j, b) in self.lines[..=i].iter().enumerate() {
                let kernel = window
                    .product_kernel(omega_a, core::f64::consts::TAU * b.frequency_hz)
                    .map_err(window_error)?;
                let pp = kernel
                    .mean_real_product_dd(a.pressure_pa, b.pressure_pa)
                    .map_err(window_error)?;
                pressure_square = pressure_square + pp;
                if i != j {
                    pressure_square = pressure_square + pp;
                }
                for (axis, sum) in intensity.iter_mut().enumerate() {
                    let pu = kernel
                        .mean_real_product_dd(a.pressure_pa, b.particle_velocity_m_s[axis])
                        .map_err(window_error)?;
                    *sum = *sum + pu;
                    if i != j {
                        // Swap amplitudes AND their frequencies: this symmetric
                        // real product uses u_a as its first argument.
                        let up = kernel
                            .mean_real_product_dd(a.particle_velocity_m_s[axis], b.pressure_pa)
                            .map_err(window_error)?;
                        *sum = *sum + up;
                    }
                }
            }
        }
        let mean_pressure_pa = pressure_mean.to_f64();
        let mean_square_pressure_pa2 = pressure_square.to_f64();
        let mean_particle_velocity_m_s = velocity_mean.map(Dd::to_f64);
        let mean_intensity_w_m2 = intensity.map(Dd::to_f64);
        if [mean_pressure_pa, mean_square_pressure_pa2]
            .iter()
            .chain(&mean_particle_velocity_m_s)
            .chain(&mean_intensity_w_m2)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("finite-window acoustic field"));
        }
        if mean_square_pressure_pa2 < 0.0 {
            return Err(invalid(
                "finite-window pressure square lost numerical positivity",
            ));
        }
        Ok(WindowedAcousticField {
            window,
            mean_pressure_pa,
            mean_square_pressure_pa2,
            mean_particle_velocity_m_s,
            mean_intensity_w_m2,
        })
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
    /// Bounds preparation and pairwise observation before any allocation.
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
/// terms. Nonzero pseudocovariance requires a different model and is not implied.
#[derive(Clone, Debug)]
pub struct PreparedProperHarmonicField<const D: usize> {
    modes: Vec<PreparedHarmonicField<D>>,
}

impl<const D: usize> PreparedProperHarmonicField<D> {
    /// Admit supplied factors and their total preparation/observation work.
    ///
    /// Preparation has the same bounded synchronous behavior as the deterministic
    /// owner. Empty factors are valid but count against `max_modes`; no factors
    /// represent silence. Shape/range failures return no partially prepared bank.
    pub fn new(
        modes: &[&[HarmonicFieldLine<D>]],
        budget: ProperHarmonicObservationBudget,
    ) -> Result<Self, AeroacError> {
        if modes.len() > budget.max_modes {
            return Err(invalid("proper harmonic field exceeds the mode budget"));
        }
        let mut lines = 0_usize;
        let mut pairs = 0_usize;
        for mode in modes {
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
            prepared.push(PreparedHarmonicField::new(
                mode,
                HarmonicObservationBudget {
                    max_lines: budget.max_total_lines,
                    max_line_pairs: budget.max_line_pairs,
                },
            )?);
        }
        Ok(Self { modes: prepared })
    }

    /// Number of authored uncorrelated factors, including empty ones.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }

    /// Observe the declared ensemble over a finite receiver-time interval.
    ///
    /// Work is O(sum_r M_r²) for M_r distinct frequencies per factor and fixed
    /// spatial dimension, without allocation. One temporal kernel per pair is
    /// reused across pressure, velocity covariance and intensity. Double-double
    /// contractions survive until the entire bank is accumulated. Numerical
    /// negativity refuses; relative accuracy at an arbitrary null is not assured.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticStatistics<D>, AeroacError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Observe with bounded cancellation polls and no partial publication.
    ///
    /// The callback is polled before work, for each factor and frequency pair,
    /// and before publication. Preparation is separately bounded and synchronous.
    /// The immutable prepared bank remains reusable after a cancelled observation.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut is_cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticStatistics<D>, AeroacError> {
        let mut pressure_square = Dd::ZERO;
        let mut velocity = [[Dd::ZERO; D]; D];
        let mut intensity = [Dd::ZERO; D];
        let mut poll = || {
            if is_cancelled() {
                Err(AeroacError::Cancelled {
                    what: "proper harmonic observation",
                })
            } else {
                Ok(())
            }
        };
        poll()?;
        for mode in &self.modes {
            poll()?;
            for (i, a) in mode.lines.iter().enumerate() {
                for (j, b) in mode.lines[..=i].iter().enumerate() {
                    poll()?;
                    // The shared Hz path retains the frequency difference and
                    // radian conversion together in extended precision; two
                    // separately rounded carriers can lose the actual beat.
                    let kernel = window
                        .proper_product_kernel_hz(a.frequency_hz, b.frequency_hz)
                        .map_err(window_error)?;
                    let pp = kernel
                        .mean_real_product_dd(a.pressure_pa, b.pressure_pa)
                        .map_err(window_error)?;
                    pressure_square = pressure_square + pp;
                    if i != j {
                        pressure_square = pressure_square + pp;
                    }
                    for (axis, sum) in intensity.iter_mut().enumerate() {
                        *sum = *sum
                            + kernel
                                .mean_real_product_dd(a.pressure_pa, b.particle_velocity_m_s[axis])
                                .map_err(window_error)?;
                        if i != j {
                            *sum = *sum
                                + kernel
                                    .mean_real_product_dd(
                                        a.particle_velocity_m_s[axis],
                                        b.pressure_pa,
                                    )
                                    .map_err(window_error)?;
                        }
                    }
                    for (row, values) in velocity.iter_mut().enumerate() {
                        for (column, sum) in values.iter_mut().enumerate().skip(row) {
                            *sum = *sum
                                + kernel
                                    .mean_real_product_dd(
                                        a.particle_velocity_m_s[row],
                                        b.particle_velocity_m_s[column],
                                    )
                                    .map_err(window_error)?;
                            if i != j {
                                *sum = *sum
                                    + kernel
                                        .mean_real_product_dd(
                                            a.particle_velocity_m_s[column],
                                            b.particle_velocity_m_s[row],
                                        )
                                        .map_err(window_error)?;
                            }
                        }
                    }
                }
            }
        }
        poll()?;
        let mean_square_pressure_pa2 = pressure_square.to_f64();
        let mean_intensity_w_m2 = intensity.map(Dd::to_f64);
        let velocity_covariance_m2_s2: [[f64; D]; D] = core::array::from_fn(|row| {
            core::array::from_fn(|column| velocity[row.min(column)][row.max(column)].to_f64())
        });
        if !mean_square_pressure_pa2.is_finite()
            || mean_intensity_w_m2.iter().any(|value| !value.is_finite())
            || velocity_covariance_m2_s2
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(nonfinite("finite-window proper acoustic statistics"));
        }
        if mean_square_pressure_pa2 < 0.0
            || (0..D).any(|axis| velocity_covariance_m2_s2[axis][axis] < 0.0)
        {
            return Err(invalid(
                "proper acoustic variance lost numerical positivity",
            ));
        }
        Ok(WindowedAcousticStatistics {
            window,
            mean_square_pressure_pa2,
            velocity_covariance_m2_s2,
            mean_intensity_w_m2,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct ComplexSum {
    re: Dd,
    im: Dd,
}

impl ComplexSum {
    fn add(&mut self, value: C64) {
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
    }

    fn value(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
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
