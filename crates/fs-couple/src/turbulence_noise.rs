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
//! see [`JetNoiseCard`] for the law, what each part rests on, and the label a page must carry. It uses the band-limited
//! core [`BandedGaussian`] under its own kernel id [`GLOTTAL_JET_NOISE_KERNEL`] and, like the flue's `jet_noise` row, five
//! required fields and no default.

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
/// ([`BandedGaussian`] under this id) never touches the glottal stream or the flue island's onset-seed stream.
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
/// `gain * rho U_j^2 / 2` over the band `St_h = f h / U_j` in `[st_lo, st_hi]`, shut while `Re_h = U_j h / nu` is below
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
/// MEASURED SHORTFALL (SunnySnow, 2026-11-03, rendered on the 14 straight-tone site cards with this constant): the noise the
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

/// The glottal jet-noise source of one island: the card, its band-limited stream and the last step's flow and gap.
#[derive(Clone, Debug)]
pub struct GlottalJetNoise {
    card: JetNoiseCard,
    noise: BandedGaussian,
    last_flow_m3_s: f64,
    last_gap_m: f64,
}

impl GlottalJetNoise {
    /// The source for `card`.
    ///
    /// # Errors
    /// The reason [`JetNoiseCard::validate`] gives.
    pub fn new(card: JetNoiseCard) -> Result<Self, &'static str> {
        card.validate()?;
        Ok(Self {
            noise: BandedGaussian::new(card.seed, GLOTTAL_JET_NOISE_KERNEL),
            card,
            last_flow_m3_s: 0.0,
            last_gap_m: 0.0,
        })
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
        let w = self.noise.step(self.card.strouhal_lo * per_hz, self.card.strouhal_hi * per_hz, dt);
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
/// gate is open, 61 dB under the cup's AC rms. Kept additive on purpose: it repeats the arithmetic of [`GlottalJetNoise`] so that
/// the landed voice path is untouched.
#[derive(Clone, Debug)]
pub struct LipJetNoise {
    card: JetNoiseCard,
    noise: BandedGaussian,
    last_flow_m3_s: f64,
    last_gap_m: f64,
}

impl LipJetNoise {
    /// The source for `card`.
    ///
    /// # Errors
    /// The reason [`JetNoiseCard::validate`] gives.
    pub fn new(card: JetNoiseCard) -> Result<Self, &'static str> {
        card.validate()?;
        Ok(Self {
            noise: BandedGaussian::new(card.seed, LIP_JET_NOISE_KERNEL),
            card,
            last_flow_m3_s: 0.0,
            last_gap_m: 0.0,
        })
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
        let w = self.noise.step(self.card.strouhal_lo * per_hz, self.card.strouhal_hi * per_hz, dt);
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

    const RHO: f64 = 1.15;
    const MU: f64 = 1.9e-5;

    fn jet_card() -> JetNoiseCard {
        JetNoiseCard { seed: 21, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: 1700.0 }
    }

    /// Gauntlet G3. The series pressure's rms is `gain rho U_j^2 / 2` at the jet speed `Q / (L h)`, so it goes as the
    /// square of the speed (the pressure-drop-proportional form).
    #[test]
    fn glottal_jet_noise_rms_is_gain_times_the_dynamic_pressure_and_goes_as_the_square_of_the_jet_speed() {
        let (l, h, dt) = (0.014, 0.8e-3, 1.0 / 192_000.0);
        let rms_at = |speed: f64| {
            let mut s = GlottalJetNoise::new(jet_card()).unwrap();
            s.record(speed * l * h, h);
            for _ in 0..4096 {
                s.pressure_pa(l, 1.0e3, RHO, MU, dt);
            }
            let n = 1usize << 17;
            let ss: f64 = (0..n).map(|_| { let p = s.pressure_pa(l, 1.0e3, RHO, MU, dt); p * p }).sum();
            (ss / n as f64).sqrt()
        };
        let u = 40.0;
        let expected = 1.3e-3 * 0.5 * RHO * u * u;
        let measured = rms_at(u);
        assert!((measured / expected - 1.0).abs() < 0.04, "rms {measured} against gain rho U^2 / 2 = {expected}");
        let doubled = rms_at(2.0 * u);
        assert!((doubled / measured - 4.0).abs() < 0.25, "twice the speed must give four times the rms, got {}", doubled / measured);
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
        let n = 1usize << 16;
        let xs: Vec<f64> = (0..n).map(|_| s.pressure_pa(l, vmax, RHO, MU, dt)).collect();
        let rms = (xs.iter().map(|x| x * x).sum::<f64>() / n as f64).sqrt();
        let expected = 1.3e-3 * 0.5 * RHO * vmax * vmax;
        assert!((rms / expected - 1.0).abs() < 0.05, "rms {rms} against the clamped {expected}");
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
    /// `Re_h = Q / (w nu) = 6000`, well past the 1700 gate; the source is then nonzero and its rms is `gain rho U_j^2 / 2`.
    #[test]
    fn lip_jet_noise_speaks_on_a_planted_open_gap_with_rms_gain_times_the_dynamic_pressure() {
        let mut s = LipJetNoise::new(jet_card()).unwrap();
        s.record(9.0e-4, 1.0e-3);
        let bound = 1.0e3; // far above the 100 m/s jet: the clamp is not what is read here
        for _ in 0..4096 {
            s.pressure_pa(LIP_WIDTH, bound, RHO, MU, LIP_DT);
        }
        let n = 1usize << 17;
        let xs: Vec<f64> = (0..n).map(|_| s.pressure_pa(LIP_WIDTH, bound, RHO, MU, LIP_DT)).collect();
        assert!(xs.iter().filter(|x| **x != 0.0).count() > n * 99 / 100, "an open gate speaks on (almost) every sample");
        let rms = (xs.iter().map(|x| x * x).sum::<f64>() / n as f64).sqrt();
        let expected = 1.3e-3 * 0.5 * RHO * 100.0 * 100.0;
        assert!((rms / expected - 1.0).abs() < 0.04, "rms {rms} against gain rho U^2 / 2 = {expected}");
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

    /// With `Re_h = U_j h / nu` and `U_j = Q / (w h)` the gap cancels: the gate is the flow threshold `Q = Re_on w nu`.
    #[test]
    fn lip_jet_noise_gate_is_a_flow_threshold_because_the_gap_cancels() {
        let q_on = 1700.0 * LIP_WIDTH * (MU / RHO);
        for h in [0.2e-3, 0.6e-3, 1.5e-3] {
            let mut shut = LipJetNoise::new(jet_card()).unwrap();
            shut.record(0.99 * q_on, h);
            assert!((0..64).all(|_| shut.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT) == 0.0), "below the flow threshold at gap {h}");
            let mut open = LipJetNoise::new(jet_card()).unwrap();
            open.record(1.01 * q_on, h);
            assert!((0..64).any(|_| open.pressure_pa(LIP_WIDTH, 1.0e3, RHO, MU, LIP_DT) != 0.0), "above the flow threshold at gap {h}");
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
