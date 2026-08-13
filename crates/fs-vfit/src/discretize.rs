//! Prewarped bilinear discretization of a fitted model to runtime
//! filter forms: factored section cascade (first-order + biquad
//! sections summed in PARALLEL, mirroring the pole-residue structure)
//! and discrete state-space.
//!
//! Bilinear map `s = K (z - 1)/(z + 1)` with `K = omega_pw /
//! tan(omega_pw * T / 2)` (prewarp: the continuous and discrete
//! responses agree EXACTLY at `omega_pw`; `K -> 2/T` as `omega_pw ->
//! 0`). The improper `s*e` term maps to `e*K (z-1)/(z+1)` — a proper
//! first-order z-section, so improper continuous models are exactly
//! representable after discretization.
//!
//! Sections are PARALLEL (summed), not cascaded products: parallel
//! form preserves the per-resonance structure the pole-residue model
//! carries, keeps coefficient sensitivity local to each resonance, and
//! makes the quantization report per-section attributable.

use crate::model::{PoleTerm, RationalModel};
use fs_math::c64::C64;
use fs_math::det;

/// One second-order (or degenerate first-order) parallel section:
/// `(b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Biquad {
    /// Numerator coefficients.
    pub b: [f64; 3],
    /// Denominator coefficients (`a0 == 1` implied).
    pub a: [f64; 2],
}

impl Biquad {
    /// Frequency response at `z = e^{i*omega*t_s}`.
    #[must_use]
    pub fn eval(&self, omega: f64, t_s: f64) -> C64 {
        let zi = C64::new(det::cos(omega * t_s), -det::sin(omega * t_s));
        let zi2 = zi * zi;
        let num = C64::from_re(self.b[0]) + zi.scale(self.b[1]) + zi2.scale(self.b[2]);
        let den = C64::ONE + zi.scale(self.a[0]) + zi2.scale(self.a[1]);
        num * den.recip()
    }

    /// The same response with coefficients rounded through f32 — the
    /// quantization-sensitivity probe.
    #[must_use]
    pub fn eval_f32_quantized(&self, omega: f64, t_s: f64) -> C64 {
        let q = Biquad {
            b: [
                f64::from(self.b[0] as f32),
                f64::from(self.b[1] as f32),
                f64::from(self.b[2] as f32),
            ],
            a: [f64::from(self.a[0] as f32), f64::from(self.a[1] as f32)],
        };
        q.eval(omega, t_s)
    }

    /// Section poles inside the unit circle?
    #[must_use]
    pub fn is_stable(&self) -> bool {
        // Jury criterion for 1 + a1 z^-1 + a2 z^-2.
        let (a1, a2) = (self.a[0], self.a[1]);
        a2.abs() < 1.0 && (a1.abs() < 1.0 + a2)
    }
}

/// A parallel bank of sections plus the sample interval.
#[derive(Debug, Clone, PartialEq)]
pub struct DigitalFilter {
    /// Parallel sections (their responses SUM).
    pub sections: Vec<Biquad>,
    /// Constant (direct) term.
    pub direct: f64,
    /// Sample interval [s].
    pub t_s: f64,
    /// Prewarp frequency actually used [rad/s] (0 = none).
    pub prewarp: f64,
}

impl DigitalFilter {
    /// Total response at angular frequency `omega`.
    #[must_use]
    pub fn eval(&self, omega: f64) -> C64 {
        let mut acc = C64::from_re(self.direct);
        for s in &self.sections {
            acc = acc + s.eval(omega, self.t_s);
        }
        acc
    }

    /// Total response with every section's coefficients f32-quantized.
    #[must_use]
    pub fn eval_f32_quantized(&self, omega: f64) -> C64 {
        let mut acc = C64::from_re(f64::from(self.direct as f32));
        for s in &self.sections {
            acc = acc + s.eval_f32_quantized(omega, self.t_s);
        }
        acc
    }

    /// All sections stable? The exact lossless differentiator the
    /// improper `e` term maps to (`b1 = -b0`, pole exactly at
    /// `z = -1`) is marginally stable BY DESIGN and is admitted.
    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.sections
            .iter()
            .all(|s| s.is_stable() || is_lossless_differentiator(s))
    }

    /// Zero DF-II memory for every parallel section.
    #[must_use]
    pub fn zero_state(&self) -> DigitalFilterState {
        DigitalFilterState {
            w: vec![[0.0; 2]; self.sections.len()],
        }
    }

    /// Scale every numerator (and the direct term) by `gain`.
    pub fn scale(&mut self, gain: f64) {
        self.direct *= gain;
        for section in &mut self.sections {
            section.b[0] *= gain;
            section.b[1] *= gain;
            section.b[2] *= gain;
        }
    }

    /// Peak `|H(ω)|` on a caller grid.
    #[must_use]
    pub fn peak_abs(&self, omegas: &[f64]) -> f64 {
        omegas
            .iter()
            .map(|&omega| self.eval(omega).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Uniformly scale so `max |H|` on `omegas` does not exceed `bound`.
    ///
    /// This is the scattering-form passivity projection `|R| ≤ 1`. It
    /// is not impedance-form residue repair; use
    /// [`realize_tabulated_impedance`] for `Re Z ≥ 0`.
    ///
    /// Returns the peak magnitude before scaling.
    pub fn enforce_abs_bound(&mut self, omegas: &[f64], bound: f64) -> f64 {
        let peak = self.peak_abs(omegas);
        if peak > bound && bound > 0.0 && peak.is_finite() {
            self.scale(bound / peak);
        }
        peak
    }

    /// Advance one sample through the parallel bank.
    ///
    /// Each section is transposed direct form II. Arithmetic order is
    /// section index then coefficient index. A non-finite result
    /// refuses without committing the section memories.
    ///
    /// # Errors
    /// [`DiscretizeError::NonFiniteRuntimeValue`] or a state/section
    /// length mismatch.
    pub fn step(&self, state: &mut DigitalFilterState, input: f64) -> Result<f64, DiscretizeError> {
        require_finite_scalar("input", input)?;
        if state.w.len() != self.sections.len() {
            return Err(DiscretizeError::InvalidRuntimeDimensions {
                field: "filter_state",
                expected: self.sections.len(),
                actual: state.w.len(),
            });
        }
        let mut output = self.direct * input;
        let mut next = vec![[0.0; 2]; self.sections.len()];
        for (i, section) in self.sections.iter().enumerate() {
            let [w1, w2] = state.w[i];
            let w0 = input - section.a[0] * w1 - section.a[1] * w2;
            let y = section.b[0] * w0 + section.b[1] * w1 + section.b[2] * w2;
            if !w0.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "biquad_w",
                    index: Some(i),
                });
            }
            if !y.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "biquad_y",
                    index: Some(i),
                });
            }
            output += y;
            next[i] = [w0, w1];
        }
        require_finite_scalar("output", output)?;
        state.w = next;
        Ok(output)
    }
}

/// Per-section DF-II memory of a [`DigitalFilter`].
#[derive(Debug, Clone, PartialEq)]
pub struct DigitalFilterState {
    w: Vec<[f64; 2]>,
}

/// Driving-point reflectance `R = (Z − Zc) / (Z + Zc)`.
///
/// This is the scattering map of any 1-D characteristic port
/// (a duct, a transmission line, a pulse tube), not an instrument
/// primitive.
#[must_use]
pub fn reflectance(impedance: C64, z_characteristic: f64) -> C64 {
    let zc = C64::from_re(z_characteristic);
    (impedance - zc) * (impedance + zc).recip()
}

/// Multiply tabulated `H(ω)` by `exp(i · sign · ω · τ)` in this
/// crate's `e^{+iωt}` convention. `sign = +1` peels a pure delay of
/// `τ` seconds from a delayed response.
///
/// # Errors
/// Length mismatch or a non-finite argument.
pub fn modulate_delay(
    omega: &[f64],
    h: &[C64],
    tau_s: f64,
    sign: f64,
) -> Result<Vec<C64>, DiscretizeError> {
    if omega.len() != h.len() {
        return Err(DiscretizeError::InvalidRuntimeDimensions {
            field: "modulate_delay",
            expected: omega.len(),
            actual: h.len(),
        });
    }
    require_finite_scalar("tau_s", tau_s)?;
    require_finite_scalar("delay_sign", sign)?;
    let mut out = Vec::with_capacity(h.len());
    for (i, (&w, &hi)) in omega.iter().zip(h.iter()).enumerate() {
        require_finite_scalar("omega", w).map_err(|_| DiscretizeError::NonFiniteRuntimeValue {
            field: "omega",
            index: Some(i),
        })?;
        if !hi.re.is_finite() || !hi.im.is_finite() {
            return Err(DiscretizeError::NonFiniteRuntimeValue {
                field: "h",
                index: Some(i),
            });
        }
        let phase = sign * w * tau_s;
        let rot = C64::new(det::cos(phase), det::sin(phase));
        out.push(hi * rot);
    }
    Ok(out)
}

/// Fit tabulated `H(iω)` (this crate's `e^{+iωt}` convention) and
/// bilinear-discretize it to a parallel [`DigitalFilter`].
///
/// A TMM impedance, a BEM radiation load, and a mobility are the
/// same object: a tabulated driving-point response. Conjugate
/// `e^{-iωt}` acoustics data before calling.
///
/// # Errors
/// Vector-fitting or bilinear refusals.
pub fn realize_tabulated(
    omega: &[f64],
    h: &[C64],
    t_s: f64,
    opts: &crate::vf::FitOptions,
    prewarp: f64,
) -> Result<DigitalFilter, RealizeError> {
    let fit = crate::vf::vector_fit(omega, h, opts).map_err(RealizeError::Fit)?;
    let nyquist = core::f64::consts::PI / t_s;
    let model = band_limit_to_nyquist(fit.model, nyquist);
    bilinear(&model, t_s, prewarp).map_err(RealizeError::Discretize)
}

/// Fit tabulated **impedance** `Z(iω)`, convex-repair passivity, then
/// bilinear-discretize.
///
/// `repair_passivity` is the impedance-form primitive (`Re Z ≥ 0`).
/// A raw fit that is already passive is returned unchanged. Repair
/// exhaustion keeps the raw stable fit — it does not invent residues.
///
/// Conjugate `e^{-iωt}` acoustics data before calling.
///
/// # Errors
/// Vector-fitting or bilinear refusals. Passivity repair failure is
/// not an error; the unrepaired model is discretized.
pub fn realize_tabulated_impedance(
    omega: &[f64],
    z: &[C64],
    t_s: f64,
    opts: &crate::vf::FitOptions,
    prewarp: f64,
) -> Result<DigitalFilter, RealizeError> {
    let fit = crate::vf::vector_fit(omega, z, opts).map_err(RealizeError::Fit)?;
    let nyquist = core::f64::consts::PI / t_s;
    let mut model = band_limit_to_nyquist(fit.model, nyquist);
    if let (Some(&lo), Some(&hi)) = (omega.first(), omega.last())
        && hi > lo
        && model.is_stable()
        && let Ok((repaired, _)) = crate::passivity::repair_passivity(&model, (lo, hi))
    {
        model = repaired;
    }
    bilinear(&model, t_s, prewarp).map_err(RealizeError::Discretize)
}

/// Drop poles that bilinear cannot represent at this sample rate.
/// Aliasing them would change the transfer function silently.
fn band_limit_to_nyquist(
    mut model: crate::model::RationalModel,
    nyquist: f64,
) -> crate::model::RationalModel {
    model.terms.retain(|term| {
        let omega = match *term {
            crate::model::PoleTerm::Real { pole, .. } => pole.abs(),
            crate::model::PoleTerm::Pair { pole, .. } => pole.abs(),
        };
        omega < nyquist
    });
    model
}

/// Typed failure from tabulated-response realization.
#[derive(Debug, Clone, PartialEq)]
pub enum RealizeError {
    /// Identification refused.
    Fit(crate::vf::VfError),
    /// Bilinear map refused.
    Discretize(DiscretizeError),
}

impl core::fmt::Display for RealizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(e) => write!(f, "tabulated realization fit refused: {e}"),
            Self::Discretize(e) => write!(f, "tabulated realization discretize refused: {e}"),
        }
    }
}

impl std::error::Error for RealizeError {}

/// A known bulk delay followed by a causal [`DigitalFilter`].
///
/// This is the time-domain characteristic port of any 1-D waveguide:
/// peel `exp(-iωτ)` out of a driving-point reflectance, fit the
/// residual, and run `delay ⊕ filter` at the sample rate. A bore, a
/// muffler, an HVAC run, and a pulse tube share this object.
#[derive(Debug, Clone)]
pub struct DelayedFilter {
    buf: Vec<f64>,
    write: usize,
    delay_int: usize,
    frac: f64,
    filter: DigitalFilter,
    state: DigitalFilterState,
    /// Residual (or full) impulse response. Empty means the IIR bank.
    fir: Vec<f64>,
    last: f64,
}

impl DelayedFilter {
    /// Build a line whose delay is `delay_samples` (must be at least 2
    /// and strictly inside the allocated buffer).
    ///
    /// # Errors
    /// A delay that does not fit, or a non-positive sample interval
    /// on the filter.
    pub fn new(delay_samples: f64, filter: DigitalFilter) -> Result<Self, DiscretizeError> {
        if !(delay_samples >= 2.0 && delay_samples.is_finite()) {
            return Err(DiscretizeError::InvalidRuntimeDimensions {
                field: "delay_samples",
                expected: 2,
                actual: delay_samples.max(0.0) as usize,
            });
        }
        if !(filter.t_s > 0.0) {
            return Err(DiscretizeError::BadSampleInterval);
        }
        let delay_int = delay_samples.floor() as usize;
        let frac = delay_samples - delay_int as f64;
        let n = delay_int + 4;
        Ok(Self {
            buf: vec![0.0; n],
            write: 0,
            delay_int,
            frac,
            state: filter.zero_state(),
            filter,
            fir: Vec::new(),
            last: 0.0,
        })
    }

    /// A characteristic port whose reflectance is a tabulated impulse
    /// response (delay included). This is the exact linear scattering
    /// map of a sampled `R(ω)`, not a rational fit.
    ///
    /// # Errors
    /// Non-finite samples, a non-positive sample interval, or an IR
    /// shorter than 4 samples.
    pub fn from_impulse_response(t_s: f64, ir: Vec<f64>) -> Result<Self, DiscretizeError> {
        if !(t_s > 0.0 && t_s.is_finite()) {
            return Err(DiscretizeError::BadSampleInterval);
        }
        if ir.len() < 4 {
            return Err(DiscretizeError::InvalidRuntimeDimensions {
                field: "impulse_response",
                expected: 4,
                actual: ir.len(),
            });
        }
        for (i, &h) in ir.iter().enumerate() {
            if !h.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "impulse_response",
                    index: Some(i),
                });
            }
        }
        let filter = DigitalFilter {
            sections: Vec::new(),
            direct: 0.0,
            t_s,
            prewarp: 0.0,
        };
        Ok(Self {
            buf: vec![0.0; ir.len()],
            write: 0,
            delay_int: 0,
            frac: 0.0,
            state: filter.zero_state(),
            filter,
            fir: ir,
            last: 0.0,
        })
    }

    /// Fit the delay-peeled response `H(ω) exp(+iωτ)` and build the
    /// line. `delay_samples = τ / t_s`.
    ///
    /// # Errors
    /// Peeling, fitting, bilinear, or delay-buffer refusals.
    pub fn from_tabulated(
        omega: &[f64],
        h: &[C64],
        delay_samples: f64,
        t_s: f64,
        opts: &crate::vf::FitOptions,
        prewarp: f64,
    ) -> Result<Self, RealizeError> {
        let peeled =
            modulate_delay(omega, h, delay_samples * t_s, 1.0).map_err(RealizeError::Discretize)?;
        let filter = realize_tabulated(omega, &peeled, t_s, opts, prewarp)?;
        Self::new(delay_samples, filter).map_err(RealizeError::Discretize)
    }

    /// Inject an outgoing sample and return the incoming wave.
    ///
    /// # Errors
    /// Non-finite input or filter overflow.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, DiscretizeError> {
        require_finite_scalar("outgoing", outgoing)?;
        let n = self.buf.len();
        self.buf[self.write] = outgoing;
        if !self.fir.is_empty() {
            let mut acc = 0.0;
            for (k, &h) in self.fir.iter().enumerate() {
                let idx = (self.write + n - k) % n;
                acc += h * self.buf[idx];
            }
            if !acc.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "fir_output",
                    index: None,
                });
            }
            self.write = (self.write + 1) % n;
            self.last = acc;
            return Ok(acc);
        }
        let i1 = (self.write + n - self.delay_int) % n;
        let i0 = (i1 + n - 1) % n;
        let delayed = (1.0 - self.frac) * self.buf[i1] + self.frac * self.buf[i0];
        self.write = (self.write + 1) % n;
        self.last = self.filter.step(&mut self.state, delayed)?;
        Ok(self.last)
    }

    /// Incoming wave from the previous [`Self::push`].
    #[must_use]
    pub fn incoming(&self) -> f64 {
        self.last
    }

    /// Scale the residual so `|H(ω)|` equals `target_abs`.
    ///
    /// A self-oscillating loop needs the speaking-frequency loop gain,
    /// not only an impulse-correct fit.
    pub fn pin_magnitude_at(&mut self, omega: f64, target_abs: f64) {
        if !self.fir.is_empty() {
            let h = fir_dtft(&self.fir, omega, self.filter.t_s).abs();
            if h > 1.0e-12 && target_abs >= 0.0 && target_abs.is_finite() {
                let gain = target_abs / h;
                for s in &mut self.fir {
                    *s *= gain;
                }
            }
            self.last = 0.0;
            return;
        }
        let h = self.filter.eval(omega).abs();
        if !(h > 1.0e-12) || !(target_abs >= 0.0 && target_abs.is_finite()) {
            return;
        }
        let gain = target_abs / h;
        self.filter.direct *= gain;
        for section in &mut self.filter.sections {
            section.b[0] *= gain;
            section.b[1] *= gain;
            section.b[2] *= gain;
        }
        self.state = self.filter.zero_state();
        self.last = 0.0;
    }

    /// Project the residual so `|H(ω)| ≤ 1` on `omegas`.
    ///
    /// A passive 1-D scatterer cannot return more than it is sent.
    /// Call this after [`Self::pin_magnitude_at`] so a speaking-frequency
    /// gain pin cannot create an active band elsewhere.
    pub fn enforce_scattering_passivity(&mut self, omegas: &[f64]) {
        if !self.fir.is_empty() {
            let peak = omegas
                .iter()
                .map(|&w| fir_dtft(&self.fir, w, self.filter.t_s).abs())
                .fold(0.0_f64, f64::max);
            if peak > 1.0 {
                let gain = 1.0 / peak;
                for s in &mut self.fir {
                    *s *= gain;
                }
            }
            self.last = 0.0;
            return;
        }
        self.filter.enforce_abs_bound(omegas, 1.0);
        self.state = self.filter.zero_state();
        self.last = 0.0;
    }
}

fn fir_dtft(ir: &[f64], omega: f64, t_s: f64) -> C64 {
    let mut acc = C64::ZERO;
    for (k, &h) in ir.iter().enumerate() {
        let phase = -omega * t_s * k as f64;
        acc = acc + C64::new(h * det::cos(phase), h * det::sin(phase));
    }
    acc
}

fn is_lossless_differentiator(s: &Biquad) -> bool {
    // Exact bitwise structural match against the section `bilinear`
    // emits for the e-term — an identity check, not a numeric
    // comparison, so strict equality is the correct operator.
    #[allow(clippy::float_cmp)]
    fn exact(s: &Biquad) -> bool {
        s.b[1] == -s.b[0] && s.b[2] == 0.0 && s.a == [1.0, 0.0]
    }
    exact(s)
}

/// Typed discretization failure.
#[derive(Debug, Clone, PartialEq)]
pub enum DiscretizeError {
    /// Sample rate too low: a pole (or the prewarp point) sits at or
    /// beyond Nyquist.
    BeyondNyquist {
        /// Offending frequency [rad/s].
        omega: f64,
        /// Nyquist [rad/s].
        nyquist: f64,
    },
    /// Non-positive sample interval.
    BadSampleInterval,
    /// A public state-space field does not match the declared state
    /// dimension.
    InvalidRuntimeDimensions {
        /// Field whose length is inconsistent.
        field: &'static str,
        /// Required number of scalar entries.
        expected: usize,
        /// Supplied number of scalar entries.
        actual: usize,
    },
    /// Squaring the declared state dimension overflowed `usize`.
    RuntimeDimensionOverflow {
        /// Declared state dimension.
        n: usize,
    },
    /// A runtime input, coefficient, state coordinate, or computed result is
    /// not finite.
    NonFiniteRuntimeValue {
        /// Value source.
        field: &'static str,
        /// Coordinate for vector/matrix fields; absent for scalars.
        index: Option<usize>,
    },
    /// The proper state-space step cannot silently omit the Tustin
    /// differentiator represented by `e_leftover`.
    UnrealizedImproperTerm {
        /// Coefficient requiring a separate differentiator section.
        e_leftover: f64,
    },
}

impl core::fmt::Display for DiscretizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DiscretizeError::BeyondNyquist { omega, nyquist } => {
                write!(
                    f,
                    "frequency {omega} rad/s at/beyond nyquist {nyquist} rad/s"
                )
            }
            DiscretizeError::BadSampleInterval => write!(f, "sample interval must be positive"),
            DiscretizeError::InvalidRuntimeDimensions {
                field,
                expected,
                actual,
            } => write!(
                f,
                "runtime state-space field {field} has length {actual}, expected {expected}"
            ),
            DiscretizeError::RuntimeDimensionOverflow { n } => {
                write!(f, "runtime state dimension {n} cannot be squared")
            }
            DiscretizeError::NonFiniteRuntimeValue { field, index } => match index {
                Some(index) => write!(f, "runtime value {field}[{index}] must be finite"),
                None => write!(f, "runtime value {field} must be finite"),
            },
            DiscretizeError::UnrealizedImproperTerm { e_leftover } => write!(
                f,
                "runtime state-space step cannot realize e_leftover={e_leftover}; use the separate differentiator section"
            ),
        }
    }
}

impl std::error::Error for DiscretizeError {}

/// Bilinear-transform the model at sample interval `t_s` with prewarp
/// at `omega_pw` (pass 0 for the unwarped `K = 2/T` map).
///
/// # Errors
/// [`DiscretizeError`] when the sample rate cannot represent the model
/// (pole resonance or prewarp at/beyond Nyquist) or `t_s <= 0`.
pub fn bilinear(
    model: &RationalModel,
    t_s: f64,
    omega_pw: f64,
) -> Result<DigitalFilter, DiscretizeError> {
    if t_s <= 0.0 || t_s.is_nan() {
        return Err(DiscretizeError::BadSampleInterval);
    }
    let nyquist = core::f64::consts::PI / t_s;
    if omega_pw >= nyquist {
        return Err(DiscretizeError::BeyondNyquist {
            omega: omega_pw,
            nyquist,
        });
    }
    for t in &model.terms {
        let w = match t {
            PoleTerm::Real { pole, .. } => pole.abs(),
            PoleTerm::Pair { pole, .. } => pole.abs(),
        };
        if w >= nyquist {
            return Err(DiscretizeError::BeyondNyquist { omega: w, nyquist });
        }
    }
    let k = if omega_pw > 0.0 {
        omega_pw / det::tan(omega_pw * t_s / 2.0)
    } else {
        2.0 / t_s
    };
    let mut sections = Vec::new();
    for t in &model.terms {
        match *t {
            PoleTerm::Real { pole, residue } => {
                // r/(s-p), s = K(z-1)/(z+1):
                //   r (1 + z^-1) / ((K - p) + (-K - p) z^-1)
                let a0 = k - pole;
                sections.push(Biquad {
                    b: [residue / a0, residue / a0, 0.0],
                    a: [(-k - pole) / a0, 0.0],
                });
            }
            PoleTerm::Pair { pole, residue } => {
                // Pair as a real second-order s-section:
                //   (beta1 s + beta0) / (s^2 + alpha1 s + alpha0)
                let (a_re, b_im) = (pole.re, pole.im);
                let (rho, sigma) = (residue.re, residue.im);
                let beta1 = 2.0 * rho;
                let beta0 = -2.0 * (rho * a_re + sigma * b_im);
                let alpha1 = -2.0 * a_re;
                let alpha0 = a_re * a_re + b_im * b_im;
                // Substitute and clear (z+1)^2:
                //   num = beta1 K (z^2 - 1) + beta0 (z+1)^2
                //   den = K^2 (z-1)^2 + alpha1 K (z^2 - 1) + alpha0 (z+1)^2
                let k2 = k * k;
                let n0 = beta1 * k + beta0; // z^2
                let n1 = 2.0 * beta0; // z^1
                let n2 = -beta1 * k + beta0; // z^0
                let d0 = k2 + alpha1 * k + alpha0;
                let d1 = -2.0 * k2 + 2.0 * alpha0;
                let d2 = k2 - alpha1 * k + alpha0;
                sections.push(Biquad {
                    b: [n0 / d0, n1 / d0, n2 / d0],
                    a: [d1 / d0, d2 / d0],
                });
            }
        }
    }
    if model.e != 0.0 {
        // e*s -> e K (z - 1)/(z + 1): first-order, pole at z = -1
        // nudged is NOT needed — bilinear puts it exactly at z=-1 which
        // rings at Nyquist; standard practice keeps it exact (the
        // section is lossless) and the band tolerance test covers it.
        sections.push(Biquad {
            b: [model.e * k, -model.e * k, 0.0],
            a: [1.0, 0.0],
        });
    }
    Ok(DigitalFilter {
        sections,
        direct: model.d,
        t_s,
        prewarp: omega_pw,
    })
}

/// Discrete state-space by the same bilinear map (Tustin):
/// `Ad = (K I - A)^{-1} (K I + A)`, `Bd = sqrt(2 K) (K I - A)^{-1} B`,
/// `Cd = sqrt(2 K) C (K I - A)^{-1}`, `Dd = D + C (K I - A)^{-1} B`.
/// The improper `e` term is NOT representable in proper discrete
/// state-space and is returned as the leftover coefficient for the
/// caller to realize as the extra first-order section (as
/// [`bilinear`] does) — a named seam, not a silent drop.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteStateSpace {
    /// State dimension.
    pub n: usize,
    /// Row-major discrete state matrix.
    pub a: Vec<f64>,
    /// Input map.
    pub b: Vec<f64>,
    /// Output map.
    pub c: Vec<f64>,
    /// Direct term.
    pub d: f64,
    /// UNREALIZED improper coefficient (see type docs).
    pub e_leftover: f64,
    /// Sample interval.
    pub t_s: f64,
}

/// Mutable, allocation-free-after-construction runtime for one proper
/// [`DiscreteStateSpace`] realization.
///
/// The lifetime binds the state to the exact realization used to construct it;
/// callers cannot accidentally step that state with a different realization.
/// A nonzero [`DiscreteStateSpace::e_leftover`] is refused at construction
/// because silently dropping its separate differentiator section would change
/// the transfer function.
#[derive(Debug)]
pub struct DiscreteStateSpaceRuntime<'realization> {
    realization: &'realization DiscreteStateSpace,
    state: Vec<f64>,
    next_state: Vec<f64>,
}

/// Tustin state-space discretization (see [`DiscreteStateSpace`]).
///
/// # Errors
/// As [`bilinear`], plus singular `(K I - A)` (a pole exactly at `K`,
/// impossible for stable models since `K > 0`).
pub fn bilinear_state_space(
    model: &RationalModel,
    t_s: f64,
    omega_pw: f64,
) -> Result<DiscreteStateSpace, DiscretizeError> {
    if t_s <= 0.0 || t_s.is_nan() {
        return Err(DiscretizeError::BadSampleInterval);
    }
    let nyquist = core::f64::consts::PI / t_s;
    if omega_pw >= nyquist {
        return Err(DiscretizeError::BeyondNyquist {
            omega: omega_pw,
            nyquist,
        });
    }
    // Same per-pole refusal as the section route (review finding: a
    // resonance at/beyond Nyquist must not alias silently through the
    // Tustin route either).
    for t in &model.terms {
        let w = match t {
            PoleTerm::Real { pole, .. } => pole.abs(),
            PoleTerm::Pair { pole, .. } => pole.abs(),
        };
        if w >= nyquist {
            return Err(DiscretizeError::BeyondNyquist { omega: w, nyquist });
        }
    }
    let k = if omega_pw > 0.0 {
        omega_pw / det::tan(omega_pw * t_s / 2.0)
    } else {
        2.0 / t_s
    };
    let ss = model.state_space();
    let n = ss.n;
    // M = K I - A, factored once.
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = -ss.a[i * n + j];
        }
        m[i * n + i] += k;
    }
    let fact = fs_la::factor::lu(&m, n)
        .map_err(|_| DiscretizeError::BeyondNyquist { omega: k, nyquist })?;
    // X = M^{-1} (K I + A) column by column; Y = M^{-1} B.
    let mut ad = vec![0.0f64; n * n];
    for col in 0..n {
        let mut rhs: Vec<f64> = (0..n).map(|row| ss.a[row * n + col]).collect();
        rhs[col] += k;
        fact.solve(&mut rhs);
        for row in 0..n {
            ad[row * n + col] = rhs[row];
        }
    }
    let mut y = ss.b.clone();
    fact.solve(&mut y);
    let g = det::sqrt(2.0 * k);
    let bd: Vec<f64> = y.iter().map(|&v| g * v).collect();
    // Cd = g * C M^{-1}  (solve M^T z = C^T).
    let mut z = ss.c.clone();
    fact.solve_transpose(&mut z);
    let cd: Vec<f64> = z.iter().map(|&v| g * v).collect();
    let mut dd = ss.d;
    for (ci, yi) in ss.c.iter().zip(&y) {
        dd += ci * yi;
    }
    Ok(DiscreteStateSpace {
        n,
        a: ad,
        b: bd,
        c: cd,
        d: dd,
        e_leftover: ss.e,
        t_s,
    })
}

impl DiscreteStateSpace {
    /// Construct a zero-state runtime bound to this proper realization.
    ///
    /// # Errors
    /// Returns [`DiscretizeError`] if public realization fields are malformed
    /// or nonfinite, or if `e_leftover` requires the separate differentiator
    /// section.
    pub fn try_runtime(&self) -> Result<DiscreteStateSpaceRuntime<'_>, DiscretizeError> {
        self.validate_runtime_realization()?;
        Ok(DiscreteStateSpaceRuntime {
            realization: self,
            state: vec![0.0; self.n],
            next_state: vec![0.0; self.n],
        })
    }

    fn validate_runtime_realization(&self) -> Result<(), DiscretizeError> {
        let matrix_len = self
            .n
            .checked_mul(self.n)
            .ok_or(DiscretizeError::RuntimeDimensionOverflow { n: self.n })?;
        require_len("a", self.a.len(), matrix_len)?;
        require_len("b", self.b.len(), self.n)?;
        require_len("c", self.c.len(), self.n)?;
        require_finite("a", &self.a)?;
        require_finite("b", &self.b)?;
        require_finite("c", &self.c)?;
        require_finite_scalar("d", self.d)?;
        require_finite_scalar("e_leftover", self.e_leftover)?;
        require_finite_scalar("t_s", self.t_s)?;
        if self.t_s <= 0.0 {
            return Err(DiscretizeError::BadSampleInterval);
        }
        if self.e_leftover != 0.0 {
            return Err(DiscretizeError::UnrealizedImproperTerm {
                e_leftover: self.e_leftover,
            });
        }
        Ok(())
    }

    /// Frequency response `Cd (z I - Ad)^{-1} Bd + Dd` at `z =
    /// e^{i*omega*t_s}` (the `e_leftover` term is the caller's
    /// section).
    ///
    /// # Errors
    /// Singular `(z I - Ad)` (resonance exactly on the unit circle).
    pub fn eval(&self, omega: f64) -> Result<C64, fs_la::eigen_complex::EigFailure> {
        let n = self.n;
        let z = C64::new(det::cos(omega * self.t_s), det::sin(omega * self.t_s));
        let mut m = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                m[i * n + j] = C64::from_re(-self.a[i * n + j]);
            }
            m[i * n + i] = m[i * n + i] + z;
        }
        let lu = fs_la::eigen_complex::lu_complex(&m, n)?;
        let mut x: Vec<C64> = self.b.iter().map(|&v| C64::from_re(v)).collect();
        lu.solve(&mut x);
        let mut acc = C64::from_re(self.d);
        for (ci, xi) in self.c.iter().zip(&x) {
            acc = acc + xi.scale(*ci);
        }
        Ok(acc)
    }
}

impl DiscreteStateSpaceRuntime<'_> {
    /// Current realization state in canonical coordinate order.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.state
    }

    /// Advance one scalar sample using `y = Cx + D u`, followed by
    /// `x_next = A x + B u`.
    ///
    /// Arithmetic order is fixed by row and coordinate index. If finite inputs
    /// overflow, the call refuses without committing the partially computed
    /// next state.
    ///
    /// # Errors
    /// [`DiscretizeError::NonFiniteRuntimeValue`] for a nonfinite input or
    /// computed output/state coordinate.
    pub fn step(&mut self, input: f64) -> Result<f64, DiscretizeError> {
        require_finite_scalar("input", input)?;

        let mut output = 0.0;
        for index in 0..self.realization.n {
            output += self.realization.c[index] * self.state[index];
        }
        output += self.realization.d * input;
        require_finite_scalar("output", output)?;

        for row in 0..self.realization.n {
            let mut value = 0.0;
            for column in 0..self.realization.n {
                value += self.realization.a[row * self.realization.n + column] * self.state[column];
            }
            value += self.realization.b[row] * input;
            if !value.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "next_state",
                    index: Some(row),
                });
            }
            self.next_state[row] = value;
        }
        core::mem::swap(&mut self.state, &mut self.next_state);
        Ok(output)
    }
}

fn require_len(field: &'static str, actual: usize, expected: usize) -> Result<(), DiscretizeError> {
    if actual == expected {
        Ok(())
    } else {
        Err(DiscretizeError::InvalidRuntimeDimensions {
            field,
            expected,
            actual,
        })
    }
}

fn require_finite(field: &'static str, values: &[f64]) -> Result<(), DiscretizeError> {
    if let Some(index) = values.iter().position(|value| !value.is_finite()) {
        Err(DiscretizeError::NonFiniteRuntimeValue {
            field,
            index: Some(index),
        })
    } else {
        Ok(())
    }
}

fn require_finite_scalar(field: &'static str, value: f64) -> Result<(), DiscretizeError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(DiscretizeError::NonFiniteRuntimeValue { field, index: None })
    }
}

#[cfg(test)]
mod runtime_tests {
    use super::*;

    fn first_order() -> DiscreteStateSpace {
        DiscreteStateSpace {
            n: 1,
            a: vec![0.5],
            b: vec![1.0],
            c: vec![2.0],
            d: 3.0,
            e_leftover: 0.0,
            t_s: 0.01,
        }
    }

    #[test]
    fn g0_runtime_step_matches_analytic_first_order_recurrence() {
        let system = first_order();
        let mut runtime = system.try_runtime().expect("proper realization");
        assert_eq!(runtime.step(1.0).unwrap(), 3.0);
        assert_eq!(runtime.state(), [1.0]);
        assert_eq!(runtime.step(2.0).unwrap(), 8.0);
        assert_eq!(runtime.state(), [2.5]);
    }

    #[test]
    fn g0_runtime_replay_is_bitwise() {
        let system = first_order();
        let mut left = system.try_runtime().unwrap();
        let mut right = system.try_runtime().unwrap();
        for input in [1.0, -0.25, 0.0, 2.5, -4.0] {
            assert_eq!(
                left.step(input).unwrap().to_bits(),
                right.step(input).unwrap().to_bits()
            );
            assert_eq!(left.state(), right.state());
        }
    }

    #[test]
    fn g0_zero_state_and_input_remain_zero() {
        let mut system = first_order();
        system.d = 0.0;
        let mut runtime = system.try_runtime().unwrap();
        assert_eq!(runtime.step(0.0).unwrap().to_bits(), 0.0f64.to_bits());
        assert_eq!(runtime.state(), [0.0]);
    }

    #[test]
    fn g0_runtime_refuses_malformed_nonfinite_and_improper_inputs() {
        let mut malformed = first_order();
        malformed.a.clear();
        assert!(matches!(
            malformed.try_runtime(),
            Err(DiscretizeError::InvalidRuntimeDimensions { field: "a", .. })
        ));
        let mut improper = first_order();
        improper.e_leftover = 0.25;
        assert!(matches!(
            improper.try_runtime(),
            Err(DiscretizeError::UnrealizedImproperTerm { .. })
        ));
        let system = first_order();
        let mut runtime = system.try_runtime().unwrap();
        assert!(matches!(
            runtime.step(f64::NAN),
            Err(DiscretizeError::NonFiniteRuntimeValue { field: "input", .. })
        ));
        assert_eq!(runtime.state(), [0.0]);
    }

    #[test]
    fn g0_impulse_response_agrees_with_frequency_evaluation() {
        let system = first_order();
        let omega = 37.0;
        let expected = system.eval(omega).unwrap();
        let mut runtime = system.try_runtime().unwrap();
        let mut measured = C64::ZERO;
        for sample in 0..96 {
            let output = runtime.step(if sample == 0 { 1.0 } else { 0.0 }).unwrap();
            let phase = -(sample as f64) * omega * system.t_s;
            measured = measured + C64::new(det::cos(phase), det::sin(phase)).scale(output);
        }
        assert!((measured - expected).abs() < 1.0e-12);
    }

    #[test]
    fn digital_filter_df2_matches_first_order_recurrence() {
        let filter = DigitalFilter {
            sections: vec![Biquad {
                b: [1.0, 0.0, 0.0],
                a: [-0.5, 0.0],
            }],
            direct: 0.0,
            t_s: 1.0 / 48_000.0,
            prewarp: 0.0,
        };
        let mut state = filter.zero_state();
        let mut y = Vec::new();
        for k in 0..6 {
            y.push(
                filter
                    .step(&mut state, if k == 0 { 1.0 } else { 0.0 })
                    .unwrap(),
            );
        }
        let expect = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125];
        for (got, want) in y.iter().zip(expect) {
            assert!((got - want).abs() < 1.0e-15);
        }
    }

    #[test]
    fn reflectance_is_the_scattering_map() {
        let zc = 100.0;
        let matched = reflectance(C64::from_re(zc), zc);
        assert!(matched.abs() < 1.0e-15);
        let open = reflectance(C64::from_re(0.0), zc);
        assert!((open + C64::from_re(1.0)).abs() < 1.0e-15);
        let rigid = reflectance(C64::from_re(1.0e12), zc);
        assert!((rigid - C64::from_re(1.0)).abs() < 1.0e-8);
    }

    #[test]
    fn modulate_delay_peels_a_pure_delay() {
        let tau = 0.003;
        let omega: Vec<f64> = (1..40).map(|k| 200.0 * f64::from(k)).collect();
        let h: Vec<C64> = omega
            .iter()
            .map(|&w| C64::new(det::cos(-w * tau), det::sin(-w * tau)).scale(0.8))
            .collect();
        let peeled = modulate_delay(&omega, &h, tau, 1.0).expect("peel");
        for z in peeled {
            assert!((z.re - 0.8).abs() < 1.0e-12);
            assert!(z.im.abs() < 1.0e-12);
        }
    }

    #[test]
    fn delayed_filter_is_a_pure_delay_when_the_filter_is_one() {
        let filter = DigitalFilter {
            sections: Vec::new(),
            direct: 1.0,
            t_s: 1.0 / 8_000.0,
            prewarp: 0.0,
        };
        let mut line = DelayedFilter::new(4.0, filter).expect("line");
        let mut out = Vec::new();
        for k in 0..8 {
            out.push(line.push(if k == 0 { 1.0 } else { 0.0 }).unwrap());
        }
        assert!(out[0].abs() < 1.0e-15);
        assert!(out[1].abs() < 1.0e-15);
        assert!(out[2].abs() < 1.0e-15);
        assert!(out[3].abs() < 1.0e-15);
        assert!((out[4] - 1.0).abs() < 1.0e-15);
        assert!(out[5].abs() < 1.0e-15);
    }

    #[test]
    fn realize_tabulated_constant_is_a_direct_term() {
        let omega: Vec<f64> = (1..80).map(|k| 50.0 * f64::from(k)).collect();
        let h = vec![C64::from_re(0.7); omega.len()];
        let mut opts = crate::vf::FitOptions::new(2);
        opts.fit_e = false;
        let filter = realize_tabulated(&omega, &h, 1.0 / 16_000.0, &opts, 0.0).expect("fit");
        assert!((filter.direct - 0.7).abs() < 0.05);
        assert!(filter.is_stable());
        let mut state = filter.zero_state();
        let y0 = filter.step(&mut state, 1.0).unwrap();
        assert!((y0 - 0.7).abs() < 0.08);
    }

    #[test]
    fn enforce_abs_bound_is_scattering_passivity() {
        let mut filter = DigitalFilter {
            sections: Vec::new(),
            direct: 1.4,
            t_s: 1.0 / 8_000.0,
            prewarp: 0.0,
        };
        let omegas = [100.0, 1_000.0, 4_000.0];
        let peak = filter.enforce_abs_bound(&omegas, 1.0);
        assert!((peak - 1.4).abs() < 1.0e-12);
        assert!(filter.peak_abs(&omegas) <= 1.0 + 1.0e-12);
        assert!((filter.direct - 1.0).abs() < 1.0e-12);
    }
}
