//! Finite-duration observations of deterministic harmonic acoustic fields.
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

use crate::AeroacError;
use fs_fft::harmonic_window::{HarmonicWindow, WindowError};
use fs_math::{c64::C64, dd::Dd, det};

/// A deterministic spectral line of an already propagated acoustic field.
///
/// Pressure and velocity share one physical clock and the same peak-amplitude,
/// negative-time harmonic convention. Different frequencies must retain their
/// actual phases; independent random source modes are not deterministic tones.
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
