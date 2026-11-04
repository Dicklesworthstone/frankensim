//! Seeded turbulence noise: the shared core of the opt-in glottal aspiration source of the voice, written so
//! that any other flow-driven noise source (a flute's jet noise, say) can use the same generator and
//! filter and supply its own amplitude law.
//!
//! THE CORE (no physics): [`ColouredGaussian`], a unit-variance Gaussian passed through a first-order low-pass
//! whose cutoff the caller may change every step. The Gaussian numbers are the repository's own: `fs_rand::Stream`
//! (counter-based Philox keyed by `(seed, kernel, tile)`, normals by Box-Muller on `fs_math::det`), under this
//! module's own kernel id [`GLOTTAL_NOISE_KERNEL`], so a declared seed gives the same bits on every platform and
//! the stream never touches another consumer's. The filter ORDER is a modelling choice of this code: the source
//! below says "a low-pass filter" and names no order.
//!
//! THE VOICE DRIVER ([`AspirationCard`], [`GlottalAspiration`]): a pressure source in series with the glottal
//! channel (Sondhi and Schroeter 1987 put the aspiration source there), gated by a threshold Reynolds number,
//! with the amplitude law of Elie and Laprie, "Acoustic impact of the glottal chink on the production of
//! fricatives: A numerical study" (preprint dated 2016-12-28, Role of glottal chink in fricatives), Eq. (4):
//!
//! ```text
//! P_n(t) = xi * w_c(t) * (Re^2 - Re_c^2) * U^3 / a^(3/2)      for Re > Re_c, and 0 otherwise
//! Re     = rho u d / mu,  u = U / a,  d = 2 sqrt(a / pi)
//! f_c    = 0.15 sqrt(pi) U / a^(3/2)        (the cutoff of the low-pass that colours w_c; Elie and Laprie, p. 6)
//! ```
//!
//! with `U` the glottal volume flow, `a` the glottal area, `w_c` unit-variance coloured Gaussian noise.
//! WHAT IS AND IS NOT FROM THE SOURCE:
//! - `Re_c = 1700` is the value Elie and Laprie use, citing Sondhi and Schroeter 1987, and they call it
//!   "arbitrarily chosen": a published value, not a measured one. It is the card's optional third field.
//! - `xi` is "an arbitrarily adjustable real constant used to control the noise level" and was tuned there
//!   to natural-speech fricative/vowel ratios (0.61e-7 .. 4.3e-7 by place of articulation, fricatives, not
//!   aspiration). No source gives it for glottal aspiration, so it has NO DEFAULT: a card without it is
//!   refused by name. Its unit here is Pa s^3 m^-6 (`P_n` in Pa, `U^3 / a^(3/2)` in m^6 s^-3).
//! - READING 1, the Reynolds number: the preprint prints `Re = 2 rho U / (mu pi sqrt(a))`. With `d = 2 sqrt(a / pi)`
//!   (the diameter behind the cutoff above) and `u = U / a`, the definition `Re = rho u d / mu` gives
//!   `2 rho U / (mu sqrt(pi a))`, which differs from the printed form by a factor `sqrt(pi)`. This code uses the
//!   definition.
//! - READING 2, the gate: printed, `max{0, xi w_c (Re^2 - Re_c^2) ...}` would pass noise BELOW `Re_c` (a negative
//!   `Re^2 - Re_c^2` times a negative `w_c`) and half-wave rectify it above. The text says the source "is
//!   activated when the Reynolds number is above a certain critical value", so this code gates on `Re > Re_c`
//!   and keeps the noise signed (zero-mean).
//! - READING 3, the colour: `w_c` is unit variance, whatever the cutoff, so that `xi` times the flow terms is
//!   the rms pressure; the cutoff changes with the flow every step, and the normalisation is that of the
//!   stationary one-pole filter at the present cutoff (quasi-stationary).
//! - READING 4, the sample: the law is quasi-steady in the island's own previous-step flow and area (one
//!   mechanics step of lag, no implicit loop).
//!
//! With no `aspiration_noise` key the island never touches any of this and a performance is bit-identical to
//! one written before the key existed.
//!
//! THE VOICE JET NOISE ([`JetNoiseCard`], [`GlottalJetNoise`]): the second opt-in glottal source, beside the one above and
//! not a replacement for it. It is Stevens' form (a series pressure source proportional to the pressure drop across the
//! constriction) with Verge's measured fraction and band from a flue jet, applied to the glottal jet as an extrapolation;
//! see [`JetNoiseCard`] for the law, what each part rests on, and the label a page must carry. It draws through the
//! flue row's band-pass [`ButterworthBandGaussian`] under its own kernel id [`GLOTTAL_JET_NOISE_KERNEL`] and, like the
//! flue's `jet_noise` row, has five required fields and no default.
//!
//! THE JET NOISE FILTER ([`ButterworthBandGaussian`]): the flue's, the voice's and the brass lip's `jet_noise` rows draw
//! through one 2-pole Butterworth band-pass, Verge's prescribed shape, normalised so that the rms inside the card's band
//! is the card's gain. A band that reaches the step's Nyquist frequency keeps only its part below it, at the whole band's
//! density, and a band wholly above it is silent ([`BandCoverage`]); see that type for the source and what it does not
//! claim. [`BandedGaussian`], the band these rows drew through before, is kept unchanged only as the tests' planted
//! reference.

use fs_math::det;
use fs_rand::{Stream, StreamKey};

/// `Re_c` of Elie and Laprie, citing Sondhi and Schroeter (1987), who call it arbitrary [-].
pub const RE_CRIT_SONDHI_SCHROETER: f64 = 1700.0;

/// Kernel id of this module's noise stream in `fs_rand`'s `(seed, kernel, tile)` key: "glot". Another driver of the
/// same core takes its own kernel id, so no two consumers share a stream.
pub const GLOTTAL_NOISE_KERNEL: u32 = 0x676c_6f74;

/// Gaussian noise through a first-order low-pass of changing cutoff, normalised to unit variance.
#[derive(Clone, Debug)]
pub struct ColouredGaussian {
    white: Stream,
    y: f64,
}

impl ColouredGaussian {
    /// A stream at `seed` under `kernel` (tile 0), filter at rest.
    #[must_use]
    pub fn new(seed: u64, kernel: u32) -> Self {
        Self {
            white: StreamKey { seed, kernel, tile: 0 }.stream(),
            y: 0.0,
        }
    }

    /// One filtered sample at `cutoff_hz` and step `dt` seconds. The pole is matched,
    /// `alpha = 1 - exp(-2 pi f_c dt)`; for a white input of unit variance the stationary output variance is
    /// `alpha / (2 - alpha)`, so the output is scaled by its inverse square root. A cutoff at or above the
    /// matched-pole limit (`alpha = 1`) is the white sample itself; a non-positive or non-finite cutoff counts as
    /// closed and is [`ColouredGaussian::skip`].
    pub fn step(&mut self, cutoff_hz: f64, dt: f64) -> f64 {
        if !(cutoff_hz > 0.0 && cutoff_hz.is_finite() && dt > 0.0 && dt.is_finite()) {
            return self.skip();
        }
        let x = self.white.next_normal();
        let alpha = (1.0 - det::exp(-core::f64::consts::TAU * cutoff_hz * dt)).clamp(0.0, 1.0);
        if alpha <= 0.0 {
            return self.skip_with(x);
        }
        self.y += alpha * (x - self.y);
        self.y * det::sqrt((2.0 - alpha) / alpha)
    }

    /// A step with the source closed: the generator advances by one white sample (so the noise at a given
    /// step does not depend on how long the gate was shut), the filter is cleared and the output is zero.
    pub fn skip(&mut self) -> f64 {
        let x = self.white.next_normal();
        self.skip_with(x)
    }

    fn skip_with(&mut self, _white: f64) -> f64 {
        self.y = 0.0;
        0.0
    }
}

/// Kernel id of the flue jet-noise stream in `fs_rand`'s `(seed, kernel, tile)` key: "jetn". The flue driver
/// ([`ButterworthBandGaussian`] under this id) never touches the glottal stream or the flue island's onset-seed stream.
pub const FLUE_JET_NOISE_KERNEL: u32 = 0x6a65_746e;

/// The matched-pole coefficient of a first-order low-pass of cutoff `cutoff_hz` at step `dt`.
fn one_pole_alpha(cutoff_hz: f64, dt: f64) -> f64 {
    (1.0 - det::exp(-core::f64::consts::TAU * cutoff_hz * dt)).clamp(0.0, 1.0)
}

/// Stationary variance of `y_hi - y_lo`, two one-pole low-passes `y_k += alpha_k (x - y_k)` driven by the SAME unit
/// white sample: `v_k = alpha_k / (2 - alpha_k)` each and the covariance `alpha_lo alpha_hi / (1 - (1 - alpha_lo)(1 - alpha_hi))`
/// (the fixed point of `c = (1 - alpha_lo)(1 - alpha_hi) c + alpha_lo alpha_hi`).
fn band_variance(alpha_lo: f64, alpha_hi: f64) -> f64 {
    let denominator = alpha_lo + alpha_hi - alpha_lo * alpha_hi;
    if denominator <= 0.0 {
        return 0.0;
    }
    alpha_hi / (2.0 - alpha_hi) + alpha_lo / (2.0 - alpha_lo) - 2.0 * alpha_lo * alpha_hi / denominator
}

/// Gaussian noise through a band, normalised to unit variance: the difference of two first-order low-passes (cutoffs
/// `f_lo < f_hi`, both changeable every step) driven by ONE white sample per step. The band edges are a modelling
/// choice of this code (a soft first-order edge on each side), not of any source that names a band.
#[derive(Clone, Debug)]
pub struct BandedGaussian {
    white: Stream,
    y_lo: f64,
    y_hi: f64,
}

impl BandedGaussian {
    /// A stream at `seed` under `kernel` (tile 0), both filters at rest.
    #[must_use]
    pub fn new(seed: u64, kernel: u32) -> Self {
        Self {
            white: StreamKey { seed, kernel, tile: 0 }.stream(),
            y_lo: 0.0,
            y_hi: 0.0,
        }
    }

    /// One band sample for edges `f_lo_hz < f_hi_hz` and step `dt` seconds, scaled by the closed-form stationary
    /// standard deviation at the present edges (quasi-stationary, like [`ColouredGaussian::step`]). Edges that are not
    /// ordered, positive and finite count as closed and are [`BandedGaussian::skip`]; edges so close or so high that
    /// the two poles coincide give exactly zero (the difference vanishes identically).
    pub fn step(&mut self, f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> f64 {
        if !(f_lo_hz > 0.0 && f_hi_hz > f_lo_hz && f_hi_hz.is_finite() && dt > 0.0 && dt.is_finite()) {
            return self.skip();
        }
        let x = self.white.next_normal();
        let (a_lo, a_hi) = (one_pole_alpha(f_lo_hz, dt), one_pole_alpha(f_hi_hz, dt));
        self.y_lo += a_lo * (x - self.y_lo);
        self.y_hi += a_hi * (x - self.y_hi);
        let variance = band_variance(a_lo, a_hi);
        if !(variance > 0.0 && variance.is_finite()) {
            return 0.0;
        }
        (self.y_hi - self.y_lo) / det::sqrt(variance)
    }

    /// A step with the source closed: exactly one white sample is drawn (the stream index does not depend on how long
    /// the gate was shut), both filters are cleared and the output is zero.
    pub fn skip(&mut self) -> f64 {
        let _ = self.white.next_normal();
        self.y_lo = 0.0;
        self.y_hi = 0.0;
        0.0
    }
}

/// The 16-point Gauss-Legendre rule on `[-1, 1]`, positive half `(node, weight)`; the rule is symmetric.
const GAUSS_LEGENDRE_16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_6),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_6),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_65),
    (0.755_404_408_355_003_1, 0.124_628_971_255_533_95),
    (0.865_631_202_387_831_6, 0.095_158_511_682_493),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_61),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_055),
];

/// The variance INSIDE the -3 dB edges of the band-pass `b s / (s^2 + b s + c)` taken through the bilinear map
/// `s = (1 - z^-1) / (1 + z^-1)` (prewarped: `u = tan(pi f dt)`, `b = u_hi - u_lo`, `c = u_lo u_hi`), for a unit white
/// input: `(1/pi) int_{w_lo}^{w_hi} |H(e^{jw})|^2 dw`. With `u = tan(w / 2)` and the low-pass prototype frequency
/// `x = (u^2 - c) / (b u)` (`x = -1` at `u_lo`, `+1` at `u_hi`, `|H|^2 = 1 / (1 + x^2)`, `du / dx = b u / r`,
/// `r = sqrt(b^2 x^2 + 4 c)`) it is `(2/pi) int_{-1}^{1} b u / (r (1 + x^2) (1 + u^2)) dx`, a smooth integrand, here by
/// 16-point Gauss-Legendre (about 1e-12 relative for the flue band; a test checks it against brute force). Without
/// warping (`u << 1`) it is `b / 2`: an analog 2-pole Butterworth band-pass puts exactly half its power inside its edges.
fn butterworth_in_band_variance(b: f64, c: f64) -> f64 {
    let mut sum = 0.0;
    for &(node, weight) in &GAUSS_LEGENDRE_16 {
        for x in [-node, node] {
            let r = det::sqrt(b * b * x * x + 4.0 * c);
            // The positive root of u^2 - b x u - c = 0, in the form without cancellation on either side of x = 0.
            let u = if x >= 0.0 { 0.5 * (b * x + r) } else { 2.0 * c / (r - b * x) };
            sum += weight * (b * u / r) / ((1.0 + x * x) * (1.0 + u * u));
        }
    }
    sum * core::f64::consts::FRAC_2_PI
}

/// One 2-pole Butterworth band-pass section and the variance of its output inside its edges for a unit white input.
/// `y[n] = b0 (x[n] - x[n-2]) - a1 y[n-1] - a2 y[n-2]`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct BandBiquad {
    b0: f64,
    a1: f64,
    a2: f64,
    in_band_variance: f64,
}

impl BandBiquad {
    /// The section with -3 dB at exactly `f_lo_hz` and `f_hi_hz` at step `dt` (prewarped bilinear map of
    /// `b s / (s^2 + b s + c)`: `a0 = 1 + b + c`, `b0 = b / a0`, `a1 = 2 (c - 1) / a0`, `a2 = (1 - b + c) / a0`). None for edges
    /// that are not ordered and positive, or an upper edge at or above the Nyquist frequency `1 / (2 dt)`.
    fn design(f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> Option<Self> {
        if !(f_lo_hz > 0.0 && f_hi_hz > f_lo_hz && dt > 0.0 && dt.is_finite() && f_hi_hz * dt < 0.5) {
            return None;
        }
        let pi = core::f64::consts::PI;
        let (u_lo, u_hi) = (det::tan(pi * f_lo_hz * dt), det::tan(pi * f_hi_hz * dt));
        let (b, c) = (u_hi - u_lo, u_lo * u_hi);
        if !(b > 0.0 && c > 0.0 && b.is_finite() && c.is_finite()) {
            return None;
        }
        let in_band_variance = butterworth_in_band_variance(b, c);
        if !(in_band_variance > 0.0 && in_band_variance.is_finite()) {
            return None;
        }
        let a0 = 1.0 + b + c;
        Some(Self { b0: b / a0, a1: 2.0 * (c - 1.0) / a0, a2: (1.0 - b + c) / a0, in_band_variance })
    }
}

/// Gaussian noise through ONE 2-pole Butterworth band-pass (a biquad), normalised so that its rms INSIDE the -3 dB edges
/// is one: the filter of the flue's, the voice's and the brass lip's `jet_noise` rows. It is Verge's own recipe for the flue
/// source, white noise through "an order two Butterworth pass-band digital filter" (M. P. Verge, PhD thesis, TU Eindhoven
/// 1995, sec. 6.5, p. 196). The recipe names no edges, and Verge never ran it: "In the results presented later, we,
/// however, did not use a turbulence noise source" (p. 198). One biquad has 6 dB/oct skirts, the slope of the pedestal
/// his measured source shows outside its peak (Fig. 5.15a, p. 160, read approximately: 8-10 x 1e-5 per analyser line over
/// St_h 0.055-0.10 and 2.5-5 x 1e-5 on either side; nothing is printed above St_h 0.135). The edges are the card's
/// Strouhal band (p. 196 prefers edges fixed in Hz; that is a separate step, not taken here), and the gain is defined
/// over that band, so the normalisation is in-band: the analog section puts
/// exactly half its power inside its edges, the bilinear one a little more (0.51-0.56 at the flue's settings), and the
/// total rms is about 1.4 times the in-band rms. [`BandedGaussian`] (two first-order low-passes, normalised on the
/// whole band) put only 0.26 of its power inside the same edges, 5.8 dB below the card's gain; it is kept unchanged,
/// used by no row, as the tests' planted reference. The coefficients and the normalisation follow the edges every step
/// (quasi-stationary, like [`BandedGaussian::step`]); the filter is direct form I, so a coefficient change acts on
/// stored inputs and outputs, not on a hidden state.
///
/// AT NYQUIST (the voice and lip rows). The edges follow the jet, and a fast jet through a
/// narrow gap puts them above the step's Nyquist frequency `f_N = 1 / (2 dt)`. A band with `f_lo < f_N <= f_hi` keeps
/// only `[f_lo, f_N)`, at the whole band's mean in-band density (unit power over `f_hi - f_lo`, so power
/// `(f_N - f_lo) / (f_hi - f_lo)` there). It is drawn through a first-order high-pass at `f_lo` (prewarped: -3 dB at
/// exactly `f_lo`), the band-pass's own lower skirt (6 dB/oct). Nothing is renormalised into the narrower band:
/// Verge gives no spectrum above his band, and moving its power below Nyquist would be a number made here. A band with
/// `f_lo >= f_N` is silent. A band wholly below Nyquist is the band-pass above, bit for bit what it was before this
/// rule. [`ButterworthBandGaussian::step_with_coverage`] says which of the three a step was.
#[derive(Clone, Debug)]
pub struct ButterworthBandGaussian {
    white: Stream,
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
}

impl ButterworthBandGaussian {
    /// A stream at `seed` under `kernel` (tile 0), the section at rest.
    #[must_use]
    pub fn new(seed: u64, kernel: u32) -> Self {
        Self { white: StreamKey { seed, kernel, tile: 0 }.stream(), x1: 0.0, x2: 0.0, y1: 0.0, y2: 0.0 }
    }

    /// One sample for -3 dB edges `f_lo_hz < f_hi_hz` and step `dt` seconds, of unit rms inside the edges. An upper edge at
    /// or above the Nyquist frequency keeps the band's part below it (the type's AT NYQUIST); edges that are not ordered,
    /// positive and finite, or a lower edge at or above the Nyquist frequency, count as closed and are
    /// [`ButterworthBandGaussian::skip`].
    pub fn step(&mut self, f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> f64 {
        self.step_with_coverage(f_lo_hz, f_hi_hz, dt).0
    }

    /// [`Self::step`], and how much of the band this step could represent.
    pub fn step_with_coverage(&mut self, f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> (f64, BandCoverage) {
        if let Some(q) = BandBiquad::design(f_lo_hz, f_hi_hz, dt) {
            let x = self.white.next_normal();
            let y = q.b0 * (x - self.x2) - q.a1 * self.y1 - q.a2 * self.y2;
            (self.x2, self.x1, self.y2, self.y1) = (self.x1, x, self.y1, y);
            return (y / det::sqrt(q.in_band_variance), BandCoverage::Whole);
        }
        let Some(h) = TruncatedHighPass::design(f_lo_hz, f_hi_hz, dt) else {
            return (self.skip(), BandCoverage::Silent);
        };
        let x = self.white.next_normal();
        let y = h.b0 * (x - self.x1) - h.a1 * self.y1;
        (self.x2, self.x1, self.y2, self.y1) = (self.x1, x, self.y1, y);
        (y * det::sqrt(h.fraction / h.in_band_variance), BandCoverage::Truncated)
    }

    /// A step with the source closed: exactly one white sample is drawn (the stream index does not depend on how long
    /// the gate was shut), the section is cleared and the output is zero.
    pub fn skip(&mut self) -> f64 {
        let _ = self.white.next_normal();
        (self.x1, self.x2, self.y1, self.y2) = (0.0, 0.0, 0.0, 0.0);
        0.0
    }
}

/// How much of a jet-noise band one step can represent at its clock (Nyquist `f_N = 1 / (2 dt)`); see
/// [`ButterworthBandGaussian`]'s AT NYQUIST.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BandCoverage {
    /// The whole band lies below `f_N`: the 2-pole band-pass.
    Whole,
    /// `f_lo < f_N <= f_hi`: only `[f_lo, f_N)` is emitted, at the whole band's mean in-band density.
    Truncated,
    /// The lower edge is at or above `f_N` (or the edges are unusable): nothing in the band can be represented.
    Silent,
}

/// A jet-noise source's open-gate steps counted by [`BandCoverage`]. A step with the gate shut is not counted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JetNoiseCoverage {
    /// Open-gate steps whose whole band lies below Nyquist.
    pub whole: u64,
    /// Open-gate steps whose band straddles Nyquist (only the part below it emitted).
    pub truncated: u64,
    /// Open-gate steps whose band lies wholly at or above Nyquist (silent although the gate was open).
    pub silent: u64,
}

impl JetNoiseCoverage {
    fn count(&mut self, coverage: BandCoverage) {
        match coverage {
            BandCoverage::Whole => self.whole += 1,
            BandCoverage::Truncated => self.truncated += 1,
            BandCoverage::Silent => self.silent += 1,
        }
    }

    /// Every open-gate step.
    #[must_use]
    pub const fn open_gate(&self) -> u64 {
        self.whole + self.truncated + self.silent
    }
}

/// The variance on `[f_lo, f_N)` of the prewarped first-order high-pass `s / (s + a)`, `a = u_lo = tan theta_lo`, for a
/// unit white input: `(2/pi) int_{theta_lo}^{pi/2} sin^2 / (sin^2 + a^2 cos^2) d theta` (`theta = pi f dt`,
/// `|H|^2 = u^2 / (u^2 + a^2)` with `u = tan theta`). In `u` it is `(2/pi) int_a^inf u^2 / ((u^2 + a^2)(1 + u^2)) du
/// = (2/pi) (pi/2 - atan a - a pi / 4) / (1 - a^2)`; near `a = 1`, where that is 0/0, the integrand in `theta` is smooth and
/// 16-point Gauss-Legendre is used instead (a test checks both against brute force).
fn high_pass_variance_to_nyquist(theta_lo: f64, u_lo: f64) -> f64 {
    let (pi, a) = (core::f64::consts::PI, u_lo);
    if (1.0 - a * a).abs() > 0.05 {
        return core::f64::consts::FRAC_2_PI * (0.5 * pi - det::atan(a) - 0.25 * pi * a) / (1.0 - a * a);
    }
    let half_pi = core::f64::consts::FRAC_PI_2;
    let (mid, half) = (0.5 * (theta_lo + half_pi), 0.5 * (half_pi - theta_lo));
    let mut sum = 0.0;
    for &(node, weight) in &GAUSS_LEGENDRE_16 {
        for x in [-node, node] {
            let theta = mid + half * x;
            let (s, c) = (det::sin(theta), det::cos(theta));
            sum += weight * (s * s) / (s * s + u_lo * u_lo * c * c);
        }
    }
    sum * half * core::f64::consts::FRAC_2_PI
}

/// The truncated band, `f_lo < f_N <= f_hi`: one first-order high-pass with -3 dB at exactly `f_lo` (prewarped bilinear
/// map of `s / (s + u_lo)`: `b0 = 1 / (1 + u_lo)`, `a1 = (u_lo - 1) / (1 + u_lo)`;
/// `y[n] = b0 (x[n] - x[n-1]) - a1 y[n-1]`), its variance on `[f_lo, f_N)` and the power that part of the band carries at
/// the whole band's density, `fraction = (f_N - f_lo) / (f_hi - f_lo)`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct TruncatedHighPass {
    b0: f64,
    a1: f64,
    in_band_variance: f64,
    fraction: f64,
}

impl TruncatedHighPass {
    /// None unless the edges are ordered, positive and finite, the lower edge below Nyquist and the upper edge at or
    /// above it.
    fn design(f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> Option<Self> {
        if !(f_lo_hz > 0.0
            && f_hi_hz > f_lo_hz
            && f_hi_hz.is_finite()
            && dt > 0.0
            && dt.is_finite()
            && f_lo_hz * dt < 0.5
            && f_hi_hz * dt >= 0.5)
        {
            return None;
        }
        let theta_lo = core::f64::consts::PI * f_lo_hz * dt;
        let u_lo = det::tan(theta_lo);
        let in_band_variance = high_pass_variance_to_nyquist(theta_lo, u_lo);
        let fraction = (0.5 / dt - f_lo_hz) / (f_hi_hz - f_lo_hz);
        if !(u_lo > 0.0 && u_lo.is_finite() && in_band_variance > 0.0 && in_band_variance.is_finite() && fraction > 0.0) {
            return None;
        }
        let a0 = 1.0 + u_lo;
        Some(Self { b0: 1.0 / a0, a1: (u_lo - 1.0) / a0, in_band_variance, fraction })
    }
}

/// Test estimator: the mean square of `xs` inside `[f_lo_hz, f_hi_hz]`, from Hann-windowed 4096-sample periodograms
/// (Goertzel on each bin whose centre lies in the band, both signs of frequency, averaged over the whole segments).
#[cfg(test)]
pub(crate) fn in_band_mean_square(xs: &[f64], f_lo_hz: f64, f_hi_hz: f64, dt: f64) -> f64 {
    const N: usize = 4096;
    let tau = core::f64::consts::TAU;
    let window: Vec<f64> = (0..N).map(|i| 0.5 - 0.5 * (tau * i as f64 / N as f64).cos()).collect();
    let w2: f64 = window.iter().map(|w| w * w).sum();
    let df = 1.0 / (N as f64 * dt);
    let (k_lo, k_hi) = ((f_lo_hz / df).ceil() as usize, (f_hi_hz / df).floor() as usize);
    let (mut total, mut segments) = (0.0, 0usize);
    for segment in xs.chunks_exact(N) {
        for k in k_lo..=k_hi {
            let coeff = 2.0 * (tau * k as f64 / N as f64).cos();
            let (mut s1, mut s2) = (0.0f64, 0.0f64);
            for (x, w) in segment.iter().zip(&window) {
                (s2, s1) = (s1, x * w + coeff * s1 - s2);
            }
            total += 2.0 * (s1 * s1 + s2 * s2 - coeff * s1 * s2) / (N as f64 * w2);
        }
        segments += 1;
    }
    total / segments as f64
}

/// The voice's opt-in aspiration card: `aspiration_noise <seed> <xi> [<re_crit>]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AspirationCard {
    /// Declared seed of the noise stream.
    pub seed: u64,
    /// The gain `xi` [Pa s^3 m^-6]. No default: no source gives it for glottal aspiration.
    pub xi: f64,
    /// The threshold Reynolds number `Re_c` [-]; [`RE_CRIT_SONDHI_SCHROETER`] when the card omits it.
    pub re_crit: f64,
}

impl AspirationCard {
    /// Refuse a card whose gain or threshold cannot be used.
    ///
    /// # Errors
    /// A named reason: `xi` must be finite and positive, `re_crit` finite and non-negative.
    pub fn validate(&self) -> Result<(), &'static str> {
        if !(self.xi.is_finite() && self.xi > 0.0) {
            return Err("aspiration_noise xi must be finite and positive");
        }
        if !(self.re_crit.is_finite() && self.re_crit >= 0.0) {
            return Err("aspiration_noise re_crit must be finite and non-negative");
        }
        Ok(())
    }
}

/// The Reynolds number of the glottal jet, `rho u d / mu` with `u = U / a` and the equivalent-circle diameter
/// `d = 2 sqrt(a / pi)`: `2 rho U / (mu sqrt(pi a))` (module doc, READING 1).
#[must_use]
pub fn reynolds_number(flow_m3_s: f64, area_m2: f64, rho: f64, mu: f64) -> f64 {
    if !(flow_m3_s > 0.0 && area_m2 > 0.0 && rho > 0.0 && mu > 0.0) {
        return 0.0;
    }
    2.0 * rho * flow_m3_s / (mu * det::sqrt(core::f64::consts::PI * area_m2))
}

/// The flow at which the jet reaches `re_crit`: the inverse of [`reynolds_number`] in the flow [m^3/s].
#[must_use]
pub fn onset_flow_m3_s(area_m2: f64, rho: f64, mu: f64, re_crit: f64) -> f64 {
    re_crit * mu * det::sqrt(core::f64::consts::PI * area_m2) / (2.0 * rho)
}

/// The cutoff of the low-pass that colours `w_c`, `0.15 sqrt(pi) U / a^(3/2)` [Hz] (Elie and Laprie, p. 6).
#[must_use]
pub fn colour_cutoff_hz(flow_m3_s: f64, area_m2: f64) -> f64 {
    if !(flow_m3_s > 0.0 && area_m2 > 0.0) {
        return 0.0;
    }
    0.15 * det::sqrt(core::f64::consts::PI) * flow_m3_s / (area_m2 * det::sqrt(area_m2))
}

/// The gated amplitude law for one sample of unit-variance coloured noise `w`: zero unless `re > re_crit`.
#[must_use]
pub fn gated_pressure_pa(xi: f64, w: f64, re: f64, re_crit: f64, flow_m3_s: f64, area_m2: f64) -> f64 {
    if !(re > re_crit && flow_m3_s > 0.0 && area_m2 > 0.0) {
        return 0.0;
    }
    xi * w * (re * re - re_crit * re_crit) * flow_m3_s * flow_m3_s * flow_m3_s / (area_m2 * det::sqrt(area_m2))
}

/// The aspiration source of one glottal island: the card, its noise stream and the last step's flow and area.
#[derive(Clone, Debug)]
pub struct GlottalAspiration {
    card: AspirationCard,
    noise: ColouredGaussian,
    last_flow_m3_s: f64,
    last_area_m2: f64,
}

impl GlottalAspiration {
    /// The source for `card`.
    ///
    /// # Errors
    /// The reason [`AspirationCard::validate`] gives.
    pub fn new(card: AspirationCard) -> Result<Self, &'static str> {
        card.validate()?;
        Ok(Self {
            noise: ColouredGaussian::new(card.seed, GLOTTAL_NOISE_KERNEL),
            card,
            last_flow_m3_s: 0.0,
            last_area_m2: 0.0,
        })
    }

    /// The card.
    #[must_use]
    pub const fn card(&self) -> &AspirationCard {
        &self.card
    }

    /// The series pressure [Pa] of this step from the last step's flow and area. The stream advances by exactly one
    /// white sample per call, gate open or shut.
    pub fn pressure_pa(&mut self, rho: f64, mu: f64, dt: f64) -> f64 {
        let (u, a) = (self.last_flow_m3_s, self.last_area_m2);
        let re = reynolds_number(u, a, rho, mu);
        if !(re > self.card.re_crit) {
            self.noise.skip();
            return 0.0;
        }
        let w = self.noise.step(colour_cutoff_hz(u, a), dt);
        gated_pressure_pa(self.card.xi, w, re, self.card.re_crit, u, a)
    }

    /// Remember this step's flow [m^3/s] and glottal area [m^2] for the next call of [`Self::pressure_pa`].
    pub fn record(&mut self, flow_m3_s: f64, area_m2: f64) {
        self.last_flow_m3_s = flow_m3_s;
        self.last_area_m2 = area_m2;
    }
}

/// Kernel id of the voice's jet-noise stream in `fs_rand`'s `(seed, kernel, tile)` key: "gjet". It is neither the
/// aspiration stream ("glot") nor the flue's jet-noise stream ("jetn"), so a card may carry both voice rows.
pub const GLOTTAL_JET_NOISE_KERNEL: u32 = 0x676a_6574;

/// The voice's opt-in jet-noise card: `jet_noise <seed> <gain> <st_lo> <st_hi> <re_on>`, every field required, no default.
///
/// THE LAW (the form is Stevens', the constant is Verge's, and putting the two on the glottis is an EXTRAPOLATION).
/// A pressure source in series with the glottal channel, proportional to the pressure drop across the constriction, of rms
/// `gain * rho U_j^2 / 2` inside the band `St_h = f h / U_j` in `[st_lo, st_hi]`, shut while `Re_h = U_j h / nu` is below
/// `re_on`. `U_j = Q / (L h)` is the speed of the glottal jet (`Q` the glottal flow, `L` the fold length, `h` the gap), so
/// `Re_h = Q / (L nu)`: the gap cancels and the gate is a flow threshold. The form, a sound-pressure source proportional to
/// the pressure drop across a constriction with weak dependence on its area, is Stevens (1971, JASA 50:1180); Krane (2005,
/// JASA 118:410) derives the same `p' ~ rho U^2` for plane-wave radiation in a tube. The CONSTANT, band and gate are
/// Verge's measurement (PhD thesis, TU Eindhoven 1995, sec. 5.5, Fig. 5.15a): a flue jet at a labium, slit 1 mm,
/// `gain` 1.3e-3 (0.92e-3 for the other analyser line width, +-3 dB), `St_h` 0.045-0.115, not determinable below
/// `Re_h` about 1700. No source measures any of it on a glottis, and the page that shows a render says so (MB's label:
/// "constriction noise from a flue-jet measurement (Verge 1995: flue jet at a labium, slit 1 mm), applied to the glottal
/// jet as an extrapolation; Stevens 1971 gives the form (a pressure source proportional to the pressure drop), not this
/// constant"). The numbers are the card's; nothing here supplies, tunes or reads one from a recording.
///
/// MEASURED SHORTFALL (SunnySnow, 2026-11-03, rendered on the 14 straight-tone site cards with this constant, through the
/// band this row used before, [`BandedGaussian`], whose rms inside the band was 2.8 to 5.8 dB under the gain at the
/// tests' settings; not re-measured on the present band): the noise the
/// law adds is about 53-64 dB under the harmonic power on the F4 notes (grader hnr, ours side, before and after), against a
/// total inter-harmonic power 28-33 dB under it in the recordings the grader compares with; no number from those recordings
/// enters the code. So at Verge's constant the law is about 20-30 dB too weak on the glottis: it leaves the grader's
/// spectral flux low (F4-a 0.37 to 0.51 against 2.8) and the hnr far from the recordings' (F4-a 72 to 58 dB against 28). It
/// is kept as the opt-in place for a constriction constant that is measured on a constriction like this one; only the
/// constant, band and gate would change, with their source.
///
/// The jet speed cannot exceed the Bernoulli speed of the whole driving pressure, `sqrt(2 p_drive / rho)`; the caller
/// passes that bound and `U_j` is clamped to it (a numerically tiny gap with a finite flow would otherwise give an
/// unbounded speed). With no `jet_noise` key the island never touches any of this and a performance is bit-identical to
/// one written before the key existed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JetNoiseCard {
    /// Declared seed of the noise stream.
    pub seed: u64,
    /// rms of the series pressure as a fraction of `rho U_j^2 / 2` [-]. No default.
    pub gain: f64,
    /// Lower edge of the band in `St_h = f h / U_j` [-].
    pub strouhal_lo: f64,
    /// Upper edge of the band in `St_h` [-].
    pub strouhal_hi: f64,
    /// `Re_h = U_j h / nu` below which the source is shut [-].
    pub reynolds_onset: f64,
}

impl JetNoiseCard {
    /// Refuse a card whose fields cannot be used.
    ///
    /// # Errors
    /// A named reason: the gain must be finite and positive, the band `0 < st_lo < st_hi` finite, the onset finite and non-negative.
    pub fn validate(&self) -> Result<(), &'static str> {
        if !(self.gain.is_finite() && self.gain > 0.0) {
            return Err("jet_noise gain must be finite and positive");
        }
        if !(self.strouhal_lo.is_finite() && self.strouhal_hi.is_finite() && self.strouhal_lo > 0.0 && self.strouhal_hi > self.strouhal_lo) {
            return Err("jet_noise needs a Strouhal band 0 < ST_LO < ST_HI");
        }
        if !(self.reynolds_onset.is_finite() && self.reynolds_onset >= 0.0) {
            return Err("jet_noise reynolds onset must be finite and non-negative");
        }
        Ok(())
    }
}

/// The glottal jet-noise source of one island: the card, its band-limited stream, the last step's flow and gap, and the
/// open-gate steps counted by how much of the band each could represent.
#[derive(Clone, Debug)]
pub struct GlottalJetNoise {
    card: JetNoiseCard,
    noise: ButterworthBandGaussian,
    last_flow_m3_s: f64,
    last_gap_m: f64,
    coverage: JetNoiseCoverage,
}

impl GlottalJetNoise {
    /// The source for `card`.
    ///
    /// # Errors
    /// The reason [`JetNoiseCard::validate`] gives.
    pub fn new(card: JetNoiseCard) -> Result<Self, &'static str> {
        card.validate()?;
        Ok(Self {
            noise: ButterworthBandGaussian::new(card.seed, GLOTTAL_JET_NOISE_KERNEL),
            card,
            last_flow_m3_s: 0.0,
            last_gap_m: 0.0,
            coverage: JetNoiseCoverage::default(),
        })
    }

    /// The open-gate steps so far, by how much of the band each could represent.
    #[must_use]
    pub const fn coverage(&self) -> JetNoiseCoverage {
        self.coverage
    }

    /// The card.
    #[must_use]
    pub const fn card(&self) -> &JetNoiseCard {
        &self.card
    }

    /// The series pressure [Pa] of this step from the last step's flow and gap (the law is quasi-steady in the island's
    /// own previous step, one mechanics step of lag and no implicit loop, as READING 4 of the aspiration source).
    /// `fold_length_m` is `L`, `max_speed_m_s` the Bernoulli bound `sqrt(2 p_drive / rho)` on the jet speed, `rho` and
    /// `mu` the gas density and dynamic viscosity. The stream advances by exactly one white sample per call, gate open
    /// or shut, so its index never depends on how long the gate was shut.
    pub fn pressure_pa(&mut self, fold_length_m: f64, max_speed_m_s: f64, rho: f64, mu: f64, dt: f64) -> f64 {
        let (q, h) = (self.last_flow_m3_s, self.last_gap_m);
        if !(q > 0.0 && h > 0.0 && fold_length_m > 0.0 && max_speed_m_s > 0.0 && rho > 0.0 && mu > 0.0) {
            let _ = self.noise.skip();
            return 0.0;
        }
        let speed = (q / (fold_length_m * h)).min(max_speed_m_s);
        let reynolds = speed * h * rho / mu;
        if !(reynolds >= self.card.reynolds_onset) {
            let _ = self.noise.skip();
            return 0.0;
        }
        let per_hz = speed / h;
        let (w, coverage) = self.noise.step_with_coverage(self.card.strouhal_lo * per_hz, self.card.strouhal_hi * per_hz, dt);
        self.coverage.count(coverage);
        self.card.gain * 0.5 * rho * speed * speed * w
    }

    /// Remember this step's flow [m^3/s] and gap [m] for the next call of [`Self::pressure_pa`].
    pub fn record(&mut self, flow_m3_s: f64, gap_m: f64) {
        self.last_flow_m3_s = flow_m3_s;
        self.last_gap_m = gap_m;
    }
}

/// Kernel id of the brass lip's jet-noise stream in `fs_rand`'s `(seed, kernel, tile)` key: "ljet". It is neither the voice's
/// ("gjet", "glot") nor the flue's ("jetn") stream.
pub const LIP_JET_NOISE_KERNEL: u32 = 0x6c6a_6574;

/// The brass lip's opt-in jet-noise source: the SAME card and the SAME law as the voice row ([`JetNoiseCard`]; Stevens' form,
/// Verge's constant, band and gate, an extrapolation, with the label and the measured shortfall written there), on the lip
/// channel instead of the glottis. `width` is the lip slit width `w`, `h` the lip gap, `Q` the jet flow the lip meters:
/// `U_j = Q / (w h)` is exactly the Bernoulli speed `sqrt(2 (p_blow - p_cup) / rho)` of the lip loop, and `Re_h = U_j h / nu =
/// Q / (w nu)`. The gate is therefore a flow threshold `Q >= Re_on w nu`, and it sees FORWARD flow only: a lip whose
/// mouth pressure is at or below the cup pressure (reverse or zero flow) has no jet and no source.
///
/// Measured on SilentBluff's V1 Bb4 hold (Cullen embouchure 19, 11.42 kPa, `--mouthpiece-wav` tap and the engine's own lip equation
/// run offline): `Re_h >= 1700` on 56 % of the steady hold, in every cycle, peaking above 5400; the source is 11.7 Pa rms while the
/// gate is open, 61 dB under the cup's AC rms (through the band this row used before, [`BandedGaussian`], whose whole-band
/// rms was the gain). Kept additive on purpose: it repeats the arithmetic of [`GlottalJetNoise`] so that
/// the landed voice path is untouched.
#[derive(Clone, Debug)]
pub struct LipJetNoise {
    card: JetNoiseCard,
    noise: ButterworthBandGaussian,
    last_flow_m3_s: f64,
    last_gap_m: f64,
    coverage: JetNoiseCoverage,
}

impl LipJetNoise {
    /// The source for `card`.
    ///
    /// # Errors
    /// The reason [`JetNoiseCard::validate`] gives.
    pub fn new(card: JetNoiseCard) -> Result<Self, &'static str> {
        card.validate()?;
        Ok(Self {
            noise: ButterworthBandGaussian::new(card.seed, LIP_JET_NOISE_KERNEL),
            card,
            last_flow_m3_s: 0.0,
            last_gap_m: 0.0,
            coverage: JetNoiseCoverage::default(),
        })
    }

    /// The open-gate samples so far, by how much of the band each could represent.
    #[must_use]
    pub const fn coverage(&self) -> JetNoiseCoverage {
        self.coverage
    }

    /// The card.
    #[must_use]
    pub const fn card(&self) -> &JetNoiseCard {
        &self.card
    }

    /// The series pressure [Pa] of this sample from the last sample's flow and gap (one sample of lag, no implicit loop).
    /// `width_m` is the lip slit width, `max_speed_m_s` the Bernoulli bound `sqrt(2 |p_blow - p_cup| / rho)` of the present
    /// pressure drop on the jet speed, `rho` and `mu` the gas density and dynamic viscosity. The stream advances by exactly
    /// one white sample per call, gate open or shut.
    pub fn pressure_pa(&mut self, width_m: f64, max_speed_m_s: f64, rho: f64, mu: f64, dt: f64) -> f64 {
        let (q, h) = (self.last_flow_m3_s, self.last_gap_m);
        if !(q > 0.0 && h > 0.0 && width_m > 0.0 && max_speed_m_s > 0.0 && rho > 0.0 && mu > 0.0) {
            let _ = self.noise.skip();
            return 0.0;
        }
        let speed = (q / (width_m * h)).min(max_speed_m_s);
        let reynolds = speed * h * rho / mu;
        if !(reynolds >= self.card.reynolds_onset) {
            let _ = self.noise.skip();
            return 0.0;
        }
        let per_hz = speed / h;
        let (w, coverage) = self.noise.step_with_coverage(self.card.strouhal_lo * per_hz, self.card.strouhal_hi * per_hz, dt);
        self.coverage.count(coverage);
        self.card.gain * 0.5 * rho * speed * speed * w
    }

    /// Remember this sample's jet flow [m^3/s] (signed: negative is reverse flow) and lip gap [m].
    pub fn record(&mut self, flow_m3_s: f64, gap_m: f64) {
        self.last_flow_m3_s = flow_m3_s;
        self.last_gap_m = gap_m;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_seed_gives_one_stream_other_seeds_and_other_kernels_other_streams() {
        let run = |seed, kernel| {
            let mut c = ColouredGaussian::new(seed, kernel);
            (0..64).map(|_| c.step(900.0, 1.0 / 96_000.0).to_bits()).collect::<Vec<u64>>()
        };
        assert_eq!(run(7, GLOTTAL_NOISE_KERNEL), run(7, GLOTTAL_NOISE_KERNEL));
        assert_ne!(run(7, GLOTTAL_NOISE_KERNEL), run(8, GLOTTAL_NOISE_KERNEL));
        assert_ne!(run(7, GLOTTAL_NOISE_KERNEL), run(7, GLOTTAL_NOISE_KERNEL + 1), "a driver's kernel id keeps its stream apart");
    }

    #[test]
    fn the_repository_stream_has_zero_mean_and_unit_variance_here() {
        let mut g = StreamKey { seed: 12_345, kernel: GLOTTAL_NOISE_KERNEL, tile: 0 }.stream();
        let n = 1usize << 17;
        let (mut s, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let z = g.next_normal();
            s += z;
            s2 += z * z;
        }
        let mean = s / n as f64;
        let var = s2 / n as f64 - mean * mean;
        assert!(mean.abs() < 0.02, "mean {mean}");
        assert!((var - 1.0).abs() < 0.03, "variance {var}");
    }

    #[test]
    fn coloured_noise_has_unit_variance_and_the_one_pole_correlation() {
        let (fc, dt) = (1000.0, 1.0 / 96_000.0);
        let mut c = ColouredGaussian::new(99, GLOTTAL_NOISE_KERNEL);
        // Settle the filter, then measure.
        for _ in 0..4096 {
            c.step(fc, dt);
        }
        let n = 1usize << 18;
        let xs: Vec<f64> = (0..n).map(|_| c.step(fc, dt)).collect();
        let var = xs.iter().map(|x| x * x).sum::<f64>() / n as f64;
        let r1 = xs.windows(2).map(|w| w[0] * w[1]).sum::<f64>() / (n - 1) as f64 / var;
        let rho = det::exp(-core::f64::consts::TAU * fc * dt);
        assert!((var - 1.0).abs() < 0.06, "variance {var}");
        assert!((r1 - rho).abs() < 0.02, "lag-1 correlation {r1} against {rho}");
    }

    #[test]
    fn the_gate_is_exactly_shut_at_and_below_the_onset_flow_and_open_above_it() {
        let (area, rho, mu) = (5.0e-6, 1.14, 1.9e-5);
        let onset = onset_flow_m3_s(area, rho, mu, RE_CRIT_SONDHI_SCHROETER);
        assert!((reynolds_number(onset, area, rho, mu) - 1700.0).abs() < 1.0e-9 * 1700.0);
        for u in [0.0, 0.5 * onset, onset * (1.0 - 1.0e-9)] {
            let re = reynolds_number(u, area, rho, mu);
            assert_eq!(gated_pressure_pa(1.0e-3, 1.0, re, 1700.0, u, area), 0.0, "flow {u}");
        }
        // At the onset itself the law is continuous: rounding may put Re one ulp above Re_c, and then (Re^2 - Re_c^2) is
        // about 1e-12 of Re^2, so the pressure is below 1e-12 Pa, not exactly zero.
        let at = gated_pressure_pa(1.0e-3, 1.0, reynolds_number(onset, area, rho, mu), 1700.0, onset, area);
        assert!(at.abs() < 1.0e-12, "pressure at the onset flow {at}");
        let u = onset * 1.001;
        let re = reynolds_number(u, area, rho, mu);
        assert!(gated_pressure_pa(1.0e-3, 1.0, re, 1700.0, u, area) > 0.0);
        assert!(gated_pressure_pa(1.0e-3, -1.0, re, 1700.0, u, area) < 0.0, "the noise stays signed");
    }

    #[test]
    fn the_law_is_linear_in_the_gain_and_in_the_noise_and_has_the_printed_flow_terms() {
        let (u, a, rho, mu) = (3.0e-4, 6.0e-6, 1.14, 1.9e-5);
        let re = reynolds_number(u, a, rho, mu);
        let p = gated_pressure_pa(2.0e-5, 0.7, re, 1700.0, u, a);
        let by_hand = 2.0e-5 * 0.7 * (re * re - 1700.0 * 1700.0) * u.powi(3) / a.powf(1.5);
        assert!((p - by_hand).abs() < 1.0e-12 * by_hand.abs());
        assert!((gated_pressure_pa(4.0e-5, 0.7, re, 1700.0, u, a) - 2.0 * p).abs() < 1.0e-12 * p.abs());
        assert!((gated_pressure_pa(2.0e-5, 1.4, re, 1700.0, u, a) - 2.0 * p).abs() < 1.0e-12 * p.abs());
        // 0.15 sqrt(pi) U / a^(3/2)
        let fc = colour_cutoff_hz(u, a);
        assert!((fc - 0.15 * core::f64::consts::PI.sqrt() * u / a.powf(1.5)).abs() < 1.0e-9 * fc);
    }

    #[test]
    fn the_card_without_a_usable_gain_or_threshold_is_refused() {
        let ok = AspirationCard { seed: 1, xi: 1.0e-5, re_crit: 1700.0 };
        assert!(GlottalAspiration::new(ok).is_ok());
        for bad in [
            AspirationCard { xi: 0.0, ..ok },
            AspirationCard { xi: -1.0e-5, ..ok },
            AspirationCard { xi: f64::NAN, ..ok },
            AspirationCard { xi: f64::INFINITY, ..ok },
            AspirationCard { re_crit: -1.0, ..ok },
            AspirationCard { re_crit: f64::NAN, ..ok },
        ] {
            assert!(GlottalAspiration::new(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn the_stream_advances_one_sample_per_step_whether_the_gate_is_open_or_shut() {
        let card = AspirationCard { seed: 5, xi: 1.0e-5, re_crit: 1700.0 };
        let (area, rho, mu, dt) = (5.0e-6, 1.14, 1.9e-5, 1.0 / 96_000.0);
        let u_open = 2.0 * onset_flow_m3_s(area, rho, mu, 1700.0);
        // A: gate open from the start. B: gate shut for 10 steps, then open. After 10 steps both consumed 10 white samples, so
        // the first filtered samples differ only by the filter state (cleared in B) and the streams agree afterwards in the white draw.
        let mut a = GlottalAspiration::new(card).unwrap();
        let mut b = GlottalAspiration::new(card).unwrap();
        a.record(u_open, area);
        for _ in 0..10 {
            a.pressure_pa(rho, mu, dt);
        }
        b.record(0.0, area);
        for _ in 0..10 {
            assert_eq!(b.pressure_pa(rho, mu, dt), 0.0);
        }
        assert_eq!(a.noise.white.index(), b.noise.white.index(), "one white sample per step either way");
        assert!(a.noise.white.index() > 0, "the stream advanced");
    }

    #[test]
    fn banded_noise_has_unit_variance_no_dc_and_one_normal_per_step_gated_or_not() {
        let (f_lo, f_hi, dt) = (1000.0, 2500.0, 1.0 / 48_000.0);
        let mut g = BandedGaussian::new(41, FLUE_JET_NOISE_KERNEL);
        for _ in 0..4096 {
            g.step(f_lo, f_hi, dt);
        }
        let n = 1usize << 18;
        let xs: Vec<f64> = (0..n).map(|_| g.step(f_lo, f_hi, dt)).collect();
        let var = xs.iter().map(|x| x * x).sum::<f64>() / n as f64;
        assert!((var - 1.0).abs() < 0.06, "variance {var}");
        // A band-pass has no DC: the variance of 256-sample block means is far below the 1/256 of white noise.
        let means: Vec<f64> = xs.chunks(256).map(|c| c.iter().sum::<f64>() / 256.0).collect();
        let var_mean = means.iter().map(|m| m * m).sum::<f64>() / means.len() as f64;
        assert!(var_mean < 0.25 / 256.0, "block-mean variance {var_mean} against white's {}", 1.0 / 256.0);
        // One white sample per step, band open or shut, and the shut band clears the filters.
        let mut open = BandedGaussian::new(7, FLUE_JET_NOISE_KERNEL);
        let mut shut = BandedGaussian::new(7, FLUE_JET_NOISE_KERNEL);
        for _ in 0..6 {
            open.step(f_lo, f_hi, dt);
        }
        for _ in 0..6 {
            assert_eq!(shut.skip(), 0.0);
        }
        assert_eq!(open.white.index(), shut.white.index());
        assert!(open.white.index() > 0);
        assert_eq!((shut.y_lo, shut.y_hi), (0.0, 0.0));
        // A closed band (unordered or non-positive edges) is a skip as well: the same number of raw draws per step as an
        // open step (`next_normal` consumes more than one raw draw, so the per-step count is measured, not assumed).
        let before = open.white.index();
        let per_step = before / 6;
        assert_eq!(per_step * 6, before, "six open steps must advance the stream evenly");
        assert_eq!(open.step(2500.0, 1000.0, dt), 0.0);
        assert_eq!(open.step(0.0, 1000.0, dt), 0.0);
        assert_eq!(open.white.index(), before + 2 * per_step);
    }

    #[test]
    fn the_banded_variance_formula_matches_a_brute_force_run_and_the_seed_and_kernel_matter() {
        // Two one-pole low-passes driven by the same white stream, the difference measured directly.
        let (f_lo, f_hi, dt) = (700.0, 1900.0, 1.0 / 44_100.0);
        let (a_lo, a_hi) = (one_pole_alpha(f_lo, dt), one_pole_alpha(f_hi, dt));
        let mut white = StreamKey { seed: 3, kernel: FLUE_JET_NOISE_KERNEL, tile: 0 }.stream();
        let (mut y_lo, mut y_hi) = (0.0f64, 0.0f64);
        for _ in 0..8192 {
            let x = white.next_normal();
            y_lo += a_lo * (x - y_lo);
            y_hi += a_hi * (x - y_hi);
        }
        let n = 1usize << 19;
        let mut s2 = 0.0;
        for _ in 0..n {
            let x = white.next_normal();
            y_lo += a_lo * (x - y_lo);
            y_hi += a_hi * (x - y_hi);
            s2 += (y_hi - y_lo) * (y_hi - y_lo);
        }
        let measured = s2 / n as f64;
        let formula = band_variance(a_lo, a_hi);
        assert!((measured - formula).abs() < 0.03 * formula, "measured {measured} against {formula}");
        let run = |seed, kernel| {
            let mut g = BandedGaussian::new(seed, kernel);
            (0..64).map(|_| g.step(f_lo, f_hi, dt).to_bits()).collect::<Vec<u64>>()
        };
        assert_eq!(run(7, FLUE_JET_NOISE_KERNEL), run(7, FLUE_JET_NOISE_KERNEL));
        assert_ne!(run(7, FLUE_JET_NOISE_KERNEL), run(8, FLUE_JET_NOISE_KERNEL));
        assert_ne!(run(7, FLUE_JET_NOISE_KERNEL), run(7, GLOTTAL_NOISE_KERNEL), "the flue's kernel keeps its stream apart from the glottal one");
    }

    // |H(e^{jw})|^2 of a section, evaluated from its coefficients (independent of the closed forms in `design`).
    fn biquad_power(q: &BandBiquad, f_hz: f64, dt: f64) -> f64 {
        let w = core::f64::consts::TAU * f_hz * dt;
        let (c1, s1, c2, s2) = (w.cos(), w.sin(), (2.0 * w).cos(), (2.0 * w).sin());
        let (nr, ni) = (q.b0 * (1.0 - c2), q.b0 * s2);
        let (dr, di) = (1.0 + q.a1 * c1 + q.a2 * c2, -(q.a1 * s1 + q.a2 * s2));
        (nr * nr + ni * ni) / (dr * dr + di * di)
    }

    // (1/pi) int |H|^2 dw over [f_lo, f_hi], composite Simpson on `intervals` (even) intervals.
    fn biquad_variance_brute(q: &BandBiquad, f_lo: f64, f_hi: f64, dt: f64, intervals: usize) -> f64 {
        let h = (f_hi - f_lo) / intervals as f64;
        let mut s = biquad_power(q, f_lo, dt) + biquad_power(q, f_hi, dt);
        for i in 1..intervals {
            s += if i % 2 == 1 { 4.0 } else { 2.0 } * biquad_power(q, f_lo + i as f64 * h, dt);
        }
        s * h / 3.0 * 2.0 * dt
    }

    #[test]
    fn butterworth_band_is_minus_3_db_at_both_edges_unity_at_its_centre_and_6_db_per_octave_outside() {
        let (f_lo, f_hi, dt) = (1000.0, 2500.0, 1.0 / 48_000.0);
        let q = BandBiquad::design(f_lo, f_hi, dt).unwrap();
        let pi = core::f64::consts::PI;
        let (u_lo, u_hi) = ((pi * f_lo * dt).tan(), (pi * f_hi * dt).tan());
        let at_u = |u: f64| biquad_power(&q, u.atan() / (pi * dt), dt);
        assert!((biquad_power(&q, f_lo, dt) - 0.5).abs() < 1.0e-12, "lower edge {}", biquad_power(&q, f_lo, dt));
        assert!((biquad_power(&q, f_hi, dt) - 0.5).abs() < 1.0e-12, "upper edge {}", biquad_power(&q, f_hi, dt));
        assert!((at_u((u_lo * u_hi).sqrt()) - 1.0).abs() < 1.0e-12, "centre {}", at_u((u_lo * u_hi).sqrt()));
        // One section: |H| goes as f below the band and as 1/f above it (in prewarped frequency), 6.02 dB per octave.
        let lower = (at_u(u_lo / 8.0) / at_u(u_lo / 16.0)).sqrt();
        let upper = (at_u(8.0 * u_hi) / at_u(16.0 * u_hi)).sqrt();
        assert!((lower - 2.0).abs() < 0.02 && (upper - 2.0).abs() < 0.02, "skirts x{lower} and x{upper} per octave");
        assert_eq!(biquad_power(&q, 0.0, dt), 0.0, "no DC");
        // Edges out of order, non-positive, or at or above Nyquist: no section.
        for (lo, hi) in [(2500.0, 1000.0), (0.0, 1000.0), (1000.0, 24_000.0), (1000.0, f64::INFINITY), (f64::NAN, 2500.0)] {
            assert!(BandBiquad::design(lo, hi, dt).is_none(), "{lo}-{hi} Hz");
        }
    }

    #[test]
    fn butterworth_in_band_variance_matches_brute_force_and_the_unwarped_half_power_limit() {
        for (f_lo, f_hi, dt) in [(1000.0, 2500.0, 1.0 / 48_000.0), (956.0, 2442.0, 1.0 / 288_000.0), (5000.0, 12_780.0, 1.0 / 48_000.0)] {
            let q = BandBiquad::design(f_lo, f_hi, dt).unwrap();
            let brute = biquad_variance_brute(&q, f_lo, f_hi, dt, 200_000);
            assert!((q.in_band_variance - brute).abs() < 1.0e-9 * brute, "{f_lo}-{f_hi} Hz: {} against {brute}", q.in_band_variance);
            // The whole band holds more: the in-band fraction is above the analog one half, the more so the nearer the band
            // is to Nyquist (0.548, 0.508 and 0.722 here by an independent quadrature).
            let total = biquad_variance_brute(&q, 1.0e-9, 0.5 / dt - 1.0e-9, dt, 2_000_000);
            let fraction = q.in_band_variance / total;
            assert!(fraction > 0.5 && fraction < 0.75, "{f_lo}-{f_hi} Hz: in-band fraction {fraction}");
        }
        // Far below Nyquist the warping vanishes and the analog result holds: b / 2 with b = u_hi - u_lo.
        let dt = 1.0e-7;
        let q = BandBiquad::design(1000.0, 2500.0, dt).unwrap();
        let b = (core::f64::consts::PI * 2500.0 * dt).tan() - (core::f64::consts::PI * 1000.0 * dt).tan();
        assert!((q.in_band_variance / (0.5 * b) - 1.0).abs() < 1.0e-5, "{} against b/2 = {}", q.in_band_variance, 0.5 * b);
    }

    /// The planted negative: the in-band rms of a long stream is one for this filter, and the same estimator puts
    /// the old band (two one-pole low-passes, whole-band normalisation) 5.8 dB below it.
    #[test]
    fn butterworth_band_noise_has_unit_rms_inside_its_edges_where_the_old_band_is_5_8_db_low() {
        let (f_lo, f_hi, dt, n) = (1000.0, 2500.0, 1.0 / 48_000.0, 1usize << 19);
        let mut g = ButterworthBandGaussian::new(41, FLUE_JET_NOISE_KERNEL);
        let mut old = BandedGaussian::new(41, FLUE_JET_NOISE_KERNEL);
        for _ in 0..4096 {
            g.step(f_lo, f_hi, dt);
            old.step(f_lo, f_hi, dt);
        }
        let xs: Vec<f64> = (0..n).map(|_| g.step(f_lo, f_hi, dt)).collect();
        let olds: Vec<f64> = (0..n).map(|_| old.step(f_lo, f_hi, dt)).collect();
        let in_band = in_band_mean_square(&xs, f_lo, f_hi, dt);
        assert!((in_band - 1.0).abs() < 0.06, "in-band mean square {in_band} (rms within 3 %)");
        let old_db = 10.0 * in_band_mean_square(&olds, f_lo, f_hi, dt).log10();
        assert!((old_db + 5.8).abs() < 0.5, "the old band's in-band level {old_db} dB");
        // The whole stream carries the out-of-band skirts too: mean square 1 / (in-band fraction).
        let q = BandBiquad::design(f_lo, f_hi, dt).unwrap();
        let expected = biquad_variance_brute(&q, 1.0e-9, 0.5 / dt - 1.0e-9, dt, 2_000_000) / q.in_band_variance;
        let total = xs.iter().map(|x| x * x).sum::<f64>() / n as f64;
        assert!((total / expected - 1.0).abs() < 0.05, "total mean square {total} against {expected}");
    }

    #[test]
    fn butterworth_band_noise_draws_one_normal_per_step_open_shut_or_out_of_range_and_is_seeded() {
        let (f_lo, f_hi, dt) = (1000.0, 2500.0, 1.0 / 48_000.0);
        let mut open = ButterworthBandGaussian::new(7, FLUE_JET_NOISE_KERNEL);
        let mut shut = ButterworthBandGaussian::new(7, FLUE_JET_NOISE_KERNEL);
        for _ in 0..6 {
            open.step(f_lo, f_hi, dt);
            assert_eq!(shut.skip(), 0.0);
        }
        assert_eq!(open.white.index(), shut.white.index());
        assert!(open.white.index() > 0);
        assert_eq!((shut.x1, shut.x2, shut.y1, shut.y2), (0.0, 0.0, 0.0, 0.0));
        let before = open.white.index();
        let per_step = before / 6;
        assert_eq!(per_step * 6, before, "six open steps must advance the stream evenly");
        assert_eq!(open.step(2500.0, 1000.0, dt), 0.0);
        assert_eq!(open.step(30_000.0, 40_000.0, dt), 0.0, "a lower edge above Nyquist is a closed band");
        assert_eq!(open.white.index(), before + 2 * per_step);
        assert_eq!((open.x1, open.y1), (0.0, 0.0), "a closed band clears the section");
        let run = |seed, kernel| {
            let mut g = ButterworthBandGaussian::new(seed, kernel);
            (0..64).map(|_| g.step(f_lo, f_hi, dt).to_bits()).collect::<Vec<u64>>()
        };
        assert_eq!(run(7, FLUE_JET_NOISE_KERNEL), run(7, FLUE_JET_NOISE_KERNEL));
        assert_ne!(run(7, FLUE_JET_NOISE_KERNEL), run(8, FLUE_JET_NOISE_KERNEL));
        let banded: Vec<u64> = {
            let mut g = BandedGaussian::new(7, FLUE_JET_NOISE_KERNEL);
            (0..64).map(|_| g.step(f_lo, f_hi, dt).to_bits()).collect()
        };
        assert_ne!(run(7, FLUE_JET_NOISE_KERNEL), banded, "same white stream, another filter");
    }

    fn high_pass_power(h: &TruncatedHighPass, f_hz: f64, dt: f64) -> f64 {
        let w = core::f64::consts::TAU * f_hz * dt;
        let (nr, ni) = (h.b0 * (1.0 - w.cos()), h.b0 * w.sin());
        let (dr, di) = (1.0 + h.a1 * w.cos(), -h.a1 * w.sin());
        (nr * nr + ni * ni) / (dr * dr + di * di)
    }

    /// The Nyquist rule. A band below Nyquist is the band-pass; one straddling it keeps `[f_lo, f_N)` at the whole band's
    /// density, power `(f_N - f_lo) / (f_hi - f_lo)` there (planted negative: renormalising the whole band into the
    /// narrower one would put 1 there, 3 dB more at this band); one wholly above it is silent and still draws one normal.
    #[test]
    fn butterworth_band_noise_straddling_nyquist_keeps_the_band_density_below_it_and_is_silent_wholly_above_it() {
        let dt = 1.0 / 48_000.0;
        let mut g = ButterworthBandGaussian::new(23, FLUE_JET_NOISE_KERNEL);
        assert_eq!(g.step_with_coverage(1000.0, 2500.0, dt).1, BandCoverage::Whole);
        assert_eq!(g.step_with_coverage(12_000.0, 36_000.0, dt).1, BandCoverage::Truncated);
        assert_eq!(g.step_with_coverage(1000.0, 24_000.0, dt).1, BandCoverage::Truncated, "an upper edge AT Nyquist");
        assert_eq!(g.step_with_coverage(24_000.0, 36_000.0, dt).1, BandCoverage::Silent, "a lower edge AT Nyquist");
        assert_eq!(g.step_with_coverage(30_000.0, 40_000.0, dt).1, BandCoverage::Silent);
        assert_eq!(g.step_with_coverage(1000.0, f64::INFINITY, dt).1, BandCoverage::Silent);
        // The high-pass: -3 dB at f_lo, 1 at Nyquist, 6 dB per octave below, and its variance on [f_lo, f_N) by brute force.
        let (f_lo, f_hi) = (12_000.0, 36_000.0);
        let h = TruncatedHighPass::design(f_lo, f_hi, dt).unwrap();
        assert!((high_pass_power(&h, f_lo, dt) - 0.5).abs() < 1.0e-12, "{}", high_pass_power(&h, f_lo, dt));
        assert!((high_pass_power(&h, 0.5 / dt, dt) - 1.0).abs() < 1.0e-12);
        let pi = core::f64::consts::PI;
        let at_u = |u: f64| high_pass_power(&h, u.atan() / (pi * dt), dt);
        let u_lo = (pi * f_lo * dt).tan();
        let slope = (at_u(u_lo / 8.0) / at_u(u_lo / 16.0)).sqrt();
        assert!((slope - 2.0).abs() < 0.02, "skirt x{slope} per octave");
        // Both branches of the variance (closed form, and Gauss-Legendre near u_lo = 1 where the closed form is 0/0).
        for lo in [1000.0, 8000.0, 11_800.0, 12_000.0, 20_000.0, 23_000.0] {
            let h = TruncatedHighPass::design(lo, f_hi, dt).unwrap();
            let intervals = 200_000;
            let step = (0.5 / dt - lo) / intervals as f64;
            let mut brute = high_pass_power(&h, lo, dt) + high_pass_power(&h, 0.5 / dt, dt);
            for i in 1..intervals {
                brute += if i % 2 == 1 { 4.0 } else { 2.0 } * high_pass_power(&h, lo + i as f64 * step, dt);
            }
            brute *= step / 3.0 * 2.0 * dt;
            assert!((h.in_band_variance - brute).abs() < 1.0e-9 * brute, "f_lo {lo}: {} against {brute}", h.in_band_variance);
        }
        assert!((h.fraction - 0.5).abs() < 1.0e-15);
        // The stream: mean square on [f_lo, f_N) is the fraction (the Nyquist bin itself left out of the estimate).
        let n = 1usize << 19;
        let mut g = ButterworthBandGaussian::new(23, FLUE_JET_NOISE_KERNEL);
        for _ in 0..4096 {
            g.step(f_lo, f_hi, dt);
        }
        let xs: Vec<f64> = (0..n).map(|_| g.step(f_lo, f_hi, dt)).collect();
        let below = in_band_mean_square(&xs, f_lo, 0.5 / dt - 1.0, dt);
        assert!((below / 0.5 - 1.0).abs() < 0.06, "mean square {below} on [f_lo, f_N) against the density's 0.5");
        // One normal per step on every branch.
        let mut a = ButterworthBandGaussian::new(9, FLUE_JET_NOISE_KERNEL);
        let mut b = ButterworthBandGaussian::new(9, FLUE_JET_NOISE_KERNEL);
        for _ in 0..6 {
            a.step(1000.0, 2500.0, dt);
        }
        let (start, per_step) = (a.white.index(), a.white.index() / 6);
        assert_ne!(a.step(f_lo, f_hi, dt), 0.0);
        assert_eq!(a.step(30_000.0, 40_000.0, dt), 0.0);
        assert_eq!(a.white.index(), start + 2 * per_step);
        for _ in 0..8 {
            b.step(f_lo, f_hi, dt);
        }
        assert_eq!(b.white.index(), 8 * per_step);
    }

    /// The flue pages never reach Nyquist, so their bytes must not move: below Nyquist every sample is, bit for bit, the
    /// section as it stood before the truncation rule (that arithmetic written out here), along a sweep of edges.
    #[test]
    fn butterworth_band_noise_below_nyquist_is_bit_for_bit_the_band_pass_alone() {
        let dt = 1.0 / 288_000.0;
        let mut g = ButterworthBandGaussian::new(13, FLUE_JET_NOISE_KERNEL);
        let mut white = StreamKey { seed: 13, kernel: FLUE_JET_NOISE_KERNEL, tile: 0 }.stream();
        let (mut x1, mut x2, mut y1, mut y2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for i in 0..20_000u32 {
            let speed = 5.0 + 30.0 * f64::from(i % 1000) / 1000.0;
            let per_hz = speed / 1.13e-3;
            let (f_lo, f_hi) = (0.045 * per_hz, 0.115 * per_hz);
            let q = BandBiquad::design(f_lo, f_hi, dt).unwrap();
            let x = white.next_normal();
            let y = q.b0 * (x - x2) - q.a1 * y1 - q.a2 * y2;
            (x2, x1, y2, y1) = (x1, x, y1, y);
            let (got, coverage) = g.step_with_coverage(f_lo, f_hi, dt);
            assert_eq!(coverage, BandCoverage::Whole);
            assert_eq!(got.to_bits(), (y / det::sqrt(q.in_band_variance)).to_bits(), "sample {i}");
        }
    }

    const RHO: f64 = 1.15;
    const MU: f64 = 1.9e-5;

    fn jet_card() -> JetNoiseCard {
        JetNoiseCard { seed: 21, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: 1700.0 }
    }

    /// Gauntlet G3. The series pressure's rms INSIDE the band `St_h` in `[st_lo, st_hi]` is `gain rho U_j^2 / 2` at the jet
    /// speed `Q / (L h)`, so it goes as the square of the speed (the pressure-drop-proportional form). Planted negative: the
    /// band this row drew through before, at the same edges and scale, measures 5.8 dB low on the same estimator.
    #[test]
    fn glottal_jet_noise_in_band_rms_is_gain_times_the_dynamic_pressure_and_goes_as_the_square_of_the_jet_speed() {
        let (l, h, dt, n) = (0.014, 0.8e-3, 1.0 / 192_000.0, 1usize << 19);
        let band = |speed: f64| (0.045 * speed / h, 0.115 * speed / h);
        let in_band_rms = |speed: f64| {
            let mut s = GlottalJetNoise::new(jet_card()).unwrap();
            s.record(speed * l * h, h);
            for _ in 0..4096 {
                s.pressure_pa(l, 1.0e3, RHO, MU, dt);
            }
            let xs: Vec<f64> = (0..n).map(|_| s.pressure_pa(l, 1.0e3, RHO, MU, dt)).collect();
            assert_eq!(s.coverage(), JetNoiseCoverage { whole: (n + 4096) as u64, truncated: 0, silent: 0 });
            let (f_lo, f_hi) = band(speed);
            in_band_mean_square(&xs, f_lo, f_hi, dt).sqrt()
        };
        let u = 40.0;
        let expected = 1.3e-3 * 0.5 * RHO * u * u;
        let measured = in_band_rms(u);
        assert!((measured / expected - 1.0).abs() < 0.03, "in-band rms {measured} against gain rho U^2 / 2 = {expected}");
        let doubled = in_band_rms(2.0 * u);
        assert!((doubled / measured - 4.0).abs() < 0.2, "twice the speed must give four times the rms, got {}", doubled / measured);
        let (f_lo, f_hi) = band(u);
        let mut old = BandedGaussian::new(21, GLOTTAL_JET_NOISE_KERNEL);
        for _ in 0..4096 {
            old.step(f_lo, f_hi, dt);
        }
        let olds: Vec<f64> = (0..n).map(|_| expected * old.step(f_lo, f_hi, dt)).collect();
        let old_db = 20.0 * (in_band_mean_square(&olds, f_lo, f_hi, dt).sqrt() / expected).log10();
        assert!((old_db + 5.8).abs() < 0.5, "the old band's in-band level {old_db} dB re gain rho U^2 / 2");
    }

    /// With `Re_h = U_j h / nu` and `U_j = Q / (L h)` the gap cancels: the gate is a flow threshold `Q = Re_on L nu`.
    #[test]
    fn glottal_jet_noise_gate_is_a_flow_threshold_because_the_gap_cancels() {
        let (l, dt) = (0.014, 1.0 / 192_000.0);
        let q_on = 1700.0 * l * (MU / RHO);
        for h in [0.2e-3, 0.6e-3, 1.5e-3] {
            let mut shut = GlottalJetNoise::new(jet_card()).unwrap();
            shut.record(0.99 * q_on, h);
            assert!((0..64).all(|_| shut.pressure_pa(l, 1.0e3, RHO, MU, dt) == 0.0), "below the flow threshold at gap {h}");
            let mut open = GlottalJetNoise::new(jet_card()).unwrap();
            open.record(1.01 * q_on, h);
            assert!((0..64).any(|_| open.pressure_pa(l, 1.0e3, RHO, MU, dt) != 0.0), "above the flow threshold at gap {h}");
        }
        // No flow, no gap, no driving pressure, or no fold: no jet and no noise.
        let mut none = GlottalJetNoise::new(jet_card()).unwrap();
        none.record(2.0 * q_on, 0.0);
        assert_eq!(none.pressure_pa(l, 1.0e3, RHO, MU, dt), 0.0);
        none.record(2.0 * q_on, 1e-3);
        assert_eq!(none.pressure_pa(l, 0.0, RHO, MU, dt), 0.0, "a zero drive bound means no jet");
    }

    /// The jet cannot be faster than the Bernoulli speed of the whole driving pressure: a tiny gap with a finite flow must
    /// not give an unbounded speed. Ungated here (onset 0) so the clamp alone is read.
    #[test]
    fn glottal_jet_noise_speed_is_bounded_by_the_bernoulli_speed_of_the_drive() {
        let (l, dt, vmax) = (0.014, 1.0 / 192_000.0, 50.0);
        let card = JetNoiseCard { reynolds_onset: 0.0, ..jet_card() };
        let mut s = GlottalJetNoise::new(card).unwrap();
        s.record(5.0e-4, 0.3e-3); // raw speed Q / (L h) = 119 m/s, above the bound
        for _ in 0..4096 {
            s.pressure_pa(l, vmax, RHO, MU, dt);
        }
        let n = 1usize << 18;
        let xs: Vec<f64> = (0..n).map(|_| s.pressure_pa(l, vmax, RHO, MU, dt)).collect();
        // The band is that of the CLAMPED speed: St_h * vmax / h.
        let rms = in_band_mean_square(&xs, 0.045 * vmax / 0.3e-3, 0.115 * vmax / 0.3e-3, dt).sqrt();
        let expected = 1.3e-3 * 0.5 * RHO * vmax * vmax;
        assert!((rms / expected - 1.0).abs() < 0.05, "in-band rms {rms} against the clamped {expected}");
        // A gap so small that the clamped jet is nowhere near turbulent is shut by the gate (Re_h on the CLAMPED speed).
        let mut tiny = GlottalJetNoise::new(jet_card()).unwrap();
        tiny.record(5.0e-4, 1.0e-9);
        assert!((0..64).all(|_| tiny.pressure_pa(l, vmax, RHO, MU, dt) == 0.0));
    }

    #[test]
    fn the_jet_noise_card_without_a_usable_field_is_refused_and_the_stream_advances_one_sample_per_step() {
        let ok = jet_card();
        assert!(ok.validate().is_ok());
        for bad in [
            JetNoiseCard { gain: 0.0, ..ok },
            JetNoiseCard { gain: f64::NAN, ..ok },
            JetNoiseCard { strouhal_lo: 0.0, ..ok },
            JetNoiseCard { strouhal_lo: 0.115, strouhal_hi: 0.045, ..ok },
            JetNoiseCard { strouhal_hi: f64::INFINITY, ..ok },
            JetNoiseCard { reynolds_onset: -1.0, ..ok },
        ] {
            assert!(GlottalJetNoise::new(bad).is_err(), "{bad:?}");
        }
        let (l, dt) = (0.014, 1.0 / 192_000.0);
        let q_on = 1700.0 * l * (MU / RHO);
        let (mut open, mut shut) = (GlottalJetNoise::new(ok).unwrap(), GlottalJetNoise::new(ok).unwrap());
        open.record(2.0 * q_on, 1e-3);
        shut.record(0.1 * q_on, 1e-3);
        for _ in 0..50 {
            open.pressure_pa(l, 1.0e3, RHO, MU, dt);
            shut.pressure_pa(l, 1.0e3, RHO, MU, dt);
        }
        assert_eq!(open.noise.white.index(), shut.noise.white.index(), "one white sample per call, gate open or shut");
        assert!(open.noise.white.index() > 0);
        let run = |seed| {
            let mut s = GlottalJetNoise::new(JetNoiseCard { seed, ..ok }).unwrap();
            s.record(2.0 * q_on, 1e-3);
            (0..64).map(|_| s.pressure_pa(l, 1.0e3, RHO, MU, dt).to_bits()).collect::<Vec<u64>>()
        };
        assert_eq!(run(7), run(7));
        assert_ne!(run(7), run(8));
        assert_ne!(GLOTTAL_JET_NOISE_KERNEL, GLOTTAL_NOISE_KERNEL);
        assert_ne!(GLOTTAL_JET_NOISE_KERNEL, FLUE_JET_NOISE_KERNEL);
    }

    const LIP_DT: f64 = 1.0 / 32_000.0;
    const LIP_WIDTH: f64 = 9.0e-3;

    /// Gauntlet G3, the POSITIVE case planted on an open gap: a 1 mm lip gap passing 9e-4 m^3/s is a 100 m/s jet at
    /// `Re_h = Q / (w nu) = 6000`, well past the 1700 gate; the source is then nonzero and its rms INSIDE the band
    /// (4500-11500 Hz, below the 16 kHz Nyquist) is `gain rho U_j^2 / 2`. Planted negative: the band this row drew through
    /// before measures 2.75 dB low there on the same estimator.
    #[test]
    fn lip_jet_noise_speaks_on_a_planted_open_gap_with_in_band_rms_gain_times_the_dynamic_pressure() {
        let mut s = LipJetNoise::new(jet_card()).unwrap();
        s.record(9.0e-4, 1.0e-3);
        let bound = 1.0e3; // far above the 100 m/s jet: the clamp is not what is read here
        for _ in 0..4096 {
            s.pressure_pa(LIP_WIDTH, bound, RHO, MU, LIP_DT);
        }
        let n = 1usize << 18;
        let xs: Vec<f64> = (0..n).map(|_| s.pressure_pa(LIP_WIDTH, bound, RHO, MU, LIP_DT)).collect();
        assert!(xs.iter().filter(|x| **x != 0.0).count() > n * 99 / 100, "an open gate speaks on (almost) every sample");
        assert_eq!(s.coverage(), JetNoiseCoverage { whole: (n + 4096) as u64, truncated: 0, silent: 0 });
        let (f_lo, f_hi) = (4500.0, 11_500.0);
        let rms = in_band_mean_square(&xs, f_lo, f_hi, LIP_DT).sqrt();
        let expected = 1.3e-3 * 0.5 * RHO * 100.0 * 100.0;
        assert!((rms / expected - 1.0).abs() < 0.03, "in-band rms {rms} against gain rho U^2 / 2 = {expected}");
        let mut old = BandedGaussian::new(21, LIP_JET_NOISE_KERNEL);
        for _ in 0..4096 {
            old.step(f_lo, f_hi, LIP_DT);
        }
        let olds: Vec<f64> = (0..n).map(|_| expected * old.step(f_lo, f_hi, LIP_DT)).collect();
        let old_db = 20.0 * (in_band_mean_square(&olds, f_lo, f_hi, LIP_DT).sqrt() / expected).log10();
        assert!((old_db + 2.75).abs() < 0.5, "the old band's in-band level {old_db} dB re gain rho U^2 / 2");
    }

    /// Gauntlet G0/G3 planted negative (MB 16406): reverse flow, where the mouth pressure is at or below the cup pressure, and
    /// zero flow or a closed lip give no source, with the SAME gap and magnitude that speak when the flow is forward.
    #[test]
    fn lip_jet_noise_is_silent_on_reverse_or_zero_flow_and_a_closed_lip() {
        let speaks = |q: f64, h: f64| {
            let mut s = LipJetNoise::new(jet_card()).unwrap();
            s.record(q, h);
            (0..256).any(|_| s.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT) != 0.0)
        };
        assert!(speaks(9.0e-4, 1.0e-3), "forward flow through an open gap speaks (the negatives below are not vacuous)");
        assert!(!speaks(-9.0e-4, 1.0e-3), "reverse flow (blow at or below the cup pressure) has no jet and no source");
        assert!(!speaks(0.0, 1.0e-3), "no flow, no jet");
        assert!(!speaks(9.0e-4, 0.0), "a closed lip has no jet");
        // A zero bound (no drive) shuts it as well.
        let mut s = LipJetNoise::new(jet_card()).unwrap();
        s.record(9.0e-4, 1.0e-3);
        assert_eq!(s.pressure_pa(LIP_WIDTH, 0.0, RHO, MU, LIP_DT), 0.0);
    }

    /// With `Re_h = U_j h / nu` and `U_j = Q / (w h)` the gap cancels: the gate is the flow threshold `Q = Re_on w nu`, at
    /// every gap (the gate's opening is read from the open-gate count, which does not depend on where the band lies). At
    /// 1.01 times the threshold the jet is `0.0284 / h` m/s: at 1 mm its band (1280-3260 Hz) lies below the 16 kHz Nyquist
    /// and the source speaks; at 0.35 mm (10.4-26.6 kHz) it straddles Nyquist and the part below speaks; at 0.2 mm
    /// (31.9-81.6 kHz) it lies wholly above Nyquist, and the open gate is silent.
    #[test]
    fn lip_jet_noise_gate_is_a_flow_threshold_because_the_gap_cancels() {
        let q_on = 1700.0 * LIP_WIDTH * (MU / RHO);
        for h in [0.2e-3, 0.35e-3, 0.6e-3, 1.0e-3, 1.5e-3] {
            let mut shut = LipJetNoise::new(jet_card()).unwrap();
            shut.record(0.99 * q_on, h);
            assert!((0..64).all(|_| shut.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT) == 0.0), "below the flow threshold at gap {h}");
            assert_eq!(shut.coverage().open_gate(), 0, "the gate stays shut below the flow threshold at gap {h}");
            let mut open = LipJetNoise::new(jet_card()).unwrap();
            open.record(1.01 * q_on, h);
            let xs: Vec<f64> = (0..64).map(|_| open.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT)).collect();
            assert_eq!(open.coverage().open_gate(), 64, "the gate opens above the flow threshold at gap {h}");
            let speaks = xs.iter().any(|x| *x != 0.0);
            match h {
                h if h > 0.3e-3 && h < 0.4e-3 => {
                    assert_eq!(open.coverage().truncated, 64, "the band straddles Nyquist at gap {h}");
                    assert!(speaks, "the part of the band below Nyquist speaks at gap {h}");
                }
                h if h < 0.3e-3 => {
                    assert_eq!(open.coverage().silent, 64, "the band lies wholly above Nyquist at gap {h}");
                    assert!(!speaks, "nothing in the band can be represented at gap {h}");
                }
                _ => {
                    assert_eq!(open.coverage().whole, 64, "the band fits below Nyquist at gap {h}");
                    assert!(speaks, "the source speaks at gap {h}");
                }
            }
        }
    }

    #[test]
    fn lip_jet_noise_card_is_refused_unusable_and_the_stream_is_its_own_one_sample_per_call() {
        let ok = jet_card();
        for bad in [
            JetNoiseCard { gain: 0.0, ..ok },
            JetNoiseCard { strouhal_lo: 0.115, strouhal_hi: 0.045, ..ok },
            JetNoiseCard { reynolds_onset: -1.0, ..ok },
        ] {
            assert!(LipJetNoise::new(bad).is_err(), "{bad:?}");
        }
        let q_on = 1700.0 * LIP_WIDTH * (MU / RHO);
        let (mut open, mut shut) = (LipJetNoise::new(ok).unwrap(), LipJetNoise::new(ok).unwrap());
        open.record(2.0 * q_on, 1e-3);
        shut.record(0.1 * q_on, 1e-3);
        for _ in 0..50 {
            open.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT);
            shut.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT);
        }
        assert_eq!(open.noise.white.index(), shut.noise.white.index(), "one white sample per call, gate open or shut");
        assert!(open.noise.white.index() > 0);
        let run = |seed| {
            let mut s = LipJetNoise::new(JetNoiseCard { seed, ..ok }).unwrap();
            s.record(2.0 * q_on, 1e-3);
            (0..64).map(|_| s.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT).to_bits()).collect::<Vec<u64>>()
        };
        assert_eq!(run(7), run(7));
        assert_ne!(run(7), run(8));
        for other in [GLOTTAL_NOISE_KERNEL, GLOTTAL_JET_NOISE_KERNEL, FLUE_JET_NOISE_KERNEL] {
            assert_ne!(LIP_JET_NOISE_KERNEL, other, "the lip's stream is its own");
        }
    }
}
