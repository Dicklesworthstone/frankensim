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
}
