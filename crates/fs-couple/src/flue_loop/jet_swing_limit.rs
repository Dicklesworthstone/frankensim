//! Opt-in limit on the forced jet's swing at the labium: card row `jet_swing_limit fletcher-1976`,
//! after `jet_spreading`. Without the row nothing here runs and every result is bit-identical.
//!
//! The forced jet's displacement `eta_n = sum_k g_k v_(n-k)` is linear, so nothing stops it
//! growing. The organ B4 card swings 25 mm across a 7 mm cut-up, on a jet 1.4 mm wide.
//!
//! # What is printed: Fletcher's model choice, not a measurement of this jet
//!
//! Fletcher, "Sound production by organ flue pipes", J. Acoust. Soc. Am. 60 (1976) 926-936. The
//! exponential law assumes a displacement small against the jet's thickness. Past that "the growth
//! law probably tends to linear rather than exponential behavior" (p. 927), and his amplification
//! `G_i = cosh(mu_i l)` (eq. 12) holds "with appropriate modification if y_i > W". His choice (p.
//! 933): "The limit for exponential growth was taken equal to the spread width W of the jet at the
//! pipe lip, and the subsequent linear growth to give an envelope with a semiangle of about 15".
//! Of those two numbers he says: "Variation of the last two of these over reasonable ranges has
//! very little effect on the calculated results."
//!
//! His eq. 9-12 give the lip amplitude `y = (v / w) cosh(mu X)` for an exit velocity of amplitude
//! `v` at frequency `w`, so growth reaches `W` at a distance `x_t` that those values fix:
//!
//! ```text
//! x_t = X acosh(W w / v) / acosh(y w / v),   0 when W w / v <= 1
//! E   = W + (X - x_t) tan 15 deg
//! ```
//!
//! `W w / v <= 1` means the jet's bodily displacement `v / w` is already past `W` at the exit. `W`
//! is the station's `b(X)` under the named spreading law (`b = W`), and `X` the cut-up. Neither
//! this row nor the card carries a coefficient.
//!
//! # Measured anchors, neither on this pipe
//!
//! * de la Cuadra and Fabre, ISMA 2004: the amplitude at the transition point is about the channel
//!   height, independent of `Vac/Ub`, at h 1 mm, Re 500, 448 Hz, `Vac/Ub` about 1e-2 and under. The
//!   organ B4 card runs at `u/U_j` 0.287, about 30 times outside that range.
//! * Verge, PhD thesis, TU Eindhoven 1995, ch. 4: linear theory fits while the displacement is under
//!   h/2, for b 0.40 to 0.60 h.
//!
//! # What is ours: the way the limit is imposed in time
//!
//! The print gives an envelope, not a waveform. This module's time form is not printed:
//!
//! 1. A gain. `eta_n = g (known_n + g_0 v_n)`, `g = min(1, E / y)`. A saturation of `eta` itself can't
//!    do what the limit is for: it leaves `eta` crossing zero at the same speed, so the cut switches
//!    just as fast. Open loop on the organ B4 render's `eta`, `E tanh(eta / E)` left the cut's third
//!    to eleventh harmonics within 0.3 dB for E from 1.41 to 3.29 mm, and a gain moved them by 9 to
//!    28 dB.
//! 2. `y` and `v` are the largest `|eta|` and `|v|` over the last two completed half cycles of the
//!    linear displacement, and `w` comes from those half cycles' length. `g` steps only at a sign
//!    change, where it multiplies almost nothing, and holds in between. In a steady state it is
//!    constant, so it adds no line of its own. He works mode by mode, with each mode's amplitude;
//!    this is one gain on the whole `eta`. The gain lags: a half cycle larger than the last two
//!    keeps the old gain until it closes, and a displacement that never changes sign closes no
//!    half cycle and isn't limited at all.
//! 3. `mu` comes from the kernel's own gain `y w / v` read through his cosh. The kernel's response
//!    along the jet is not a cosh.
//!
//! The limit acts on amplitude only, so the jet's convection delay is the linear kernel's. The
//! current-feedback guard of the same-sample solve still reads the linear slope, which bounds the
//! limited one. NOT CLAIMED: a measured limit on any pipe here, the vortex street he sets aside, or
//! any phase change past the transition.

use fs_math::det;

/// The row's one accepted law.
pub const NAME: &str = "fletcher-1976";

/// The semiangle of Fletcher's linear envelope [degrees], J. Acoust. Soc. Am. 60 (1976) p. 933,
/// "about 15".
pub const FLETCHER_1976_ENVELOPE_SEMIANGLE_DEG: f64 = 15.0;

/// The envelope, as the voice provenance and the render sidecar print it.
pub const LAW: &str =
    "E = W + (X - x_t) tan 15 deg; x_t = X acosh(W w/v) / acosh(y w/v), 0 when W w/v <= 1";

/// What the print gives, as the voice provenance and the render sidecar print it.
pub const PRINTED: &str = "Fletcher 1976 pp. 927, 933, his model choice, carried to this jet, not measured on it: \
     exponential growth up to the spread width W at the lip, then linear growth with an envelope of semiangle about 15 deg; \
     he says varying these \"has very little effect\". Measured anchors, neither on this pipe: de la Cuadra and Fabre 2004, \
     amplitude at transition about h at Vac/Ub about 1e-2 and under; Verge 1995 ch. 4, linear to h/2";

/// What the print doesn't give, as the voice provenance and the render sidecar print it.
pub const OURS: &str = "the time form is ours, not printed: one gain g = min(1, E/y) on the whole linear displacement, \
     stepped at its sign changes, with y and v its peaks over the last two half cycles and w from their length; \
     mu from the kernel's gain y w/v read through his cosh";

/// `acosh` for `x >= 1`, on the deterministic log and root.
fn acosh(x: f64) -> f64 {
    det::ln(x + det::sqrt(x * x - 1.0))
}

/// What one half cycle of the linear displacement holds.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct HalfCycle {
    peak_m: f64,
    peak_velocity_m_s: f64,
    samples: u64,
}

/// The limit's state on one voice.
#[derive(Debug, Clone, PartialEq)]
pub struct JetSwingLimit {
    spread_width_m: f64,
    cut_up_m: f64,
    dt_s: f64,
    gain: f64,
    positive: bool,
    running: HalfCycle,
    last: HalfCycle,
}

impl JetSwingLimit {
    /// A limit at unit gain, for spread width `W` at the lip, cut-up `X` and the voice's time step.
    ///
    /// # Errors
    /// Any of them not finite and positive.
    pub fn new(spread_width_m: f64, cut_up_m: f64, dt_s: f64) -> Result<Self, &'static str> {
        if ![spread_width_m, cut_up_m, dt_s]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
        {
            return Err(
                "jet swing limit needs a finite positive spread width, cut-up and time step",
            );
        }
        Ok(Self {
            spread_width_m,
            cut_up_m,
            dt_s,
            gain: 1.0,
            positive: true,
            running: HalfCycle::default(),
            last: HalfCycle::default(),
        })
    }

    /// Move the lip: a jet-length track's spread width `W` at the new station and its cut-up `X`.
    /// The half cycles carry across, and the next gain is set from these at the next sign change.
    ///
    /// # Errors
    /// Either not finite and positive; the limit is then unchanged.
    pub fn set_lip(&mut self, spread_width_m: f64, cut_up_m: f64) -> Result<(), &'static str> {
        if ![spread_width_m, cut_up_m].iter().all(|x| x.is_finite() && *x > 0.0) {
            return Err("jet swing limit needs a finite positive spread width and cut-up");
        }
        self.spread_width_m = spread_width_m;
        self.cut_up_m = cut_up_m;
        Ok(())
    }

    /// The spread width `W` at the lip [m].
    #[must_use]
    pub const fn spread_width_m(&self) -> f64 {
        self.spread_width_m
    }

    /// The cut-up `X` [m].
    #[must_use]
    pub const fn cut_up_m(&self) -> f64 {
        self.cut_up_m
    }

    /// The gain on this sample's linear displacement, in `(0, 1]`.
    #[must_use]
    pub const fn gain(&self) -> f64 {
        self.gain
    }

    /// The transition distance `x_t` [m] and the envelope `E` [m], for linear lip amplitude
    /// `swing_m` above `W`, exit velocity amplitude `velocity_m_s` and frequency `omega_rad_s`.
    /// With no velocity at all, growth never leaves the exit's own displacement: `x_t = X`.
    #[must_use]
    pub fn envelope_m(&self, swing_m: f64, velocity_m_s: f64, omega_rad_s: f64) -> (f64, f64) {
        let at_width = self.spread_width_m * omega_rad_s / velocity_m_s;
        let transition_m = if at_width <= 1.0 {
            0.0
        } else {
            let share = acosh(at_width) / acosh(swing_m * omega_rad_s / velocity_m_s);
            if share.is_finite() {
                self.cut_up_m * share
            } else {
                self.cut_up_m
            }
        };
        let tan = det::tan(FLETCHER_1976_ENVELOPE_SEMIANGLE_DEG.to_radians());
        (
            transition_m,
            self.spread_width_m + (self.cut_up_m - transition_m) * tan,
        )
    }

    /// One accepted sample's linear displacement [m] and exit velocity [m/s]. A sign change closes
    /// the half cycle and sets the next samples' gain.
    pub fn accept(&mut self, linear_eta_m: f64, velocity_m_s: f64) {
        if linear_eta_m != 0.0 && (linear_eta_m > 0.0) != self.positive {
            self.gain = self.cycle_gain();
            self.positive = linear_eta_m > 0.0;
            self.last = self.running;
            self.running = HalfCycle::default();
        }
        self.running.peak_m = self.running.peak_m.max(linear_eta_m.abs());
        self.running.peak_velocity_m_s = self.running.peak_velocity_m_s.max(velocity_m_s.abs());
        self.running.samples += 1;
    }

    /// The gain from the half cycle just closed and the one before it. The first half cycle has
    /// no predecessor and stands for a whole cycle at twice its length.
    fn cycle_gain(&self) -> f64 {
        let swing_m = self.running.peak_m.max(self.last.peak_m);
        if swing_m <= self.spread_width_m {
            return 1.0;
        }
        let velocity = self
            .running
            .peak_velocity_m_s
            .max(self.last.peak_velocity_m_s);
        let samples = if self.last.samples == 0 {
            2 * self.running.samples
        } else {
            self.running.samples + self.last.samples
        };
        let omega = core::f64::consts::TAU / (samples as f64 * self.dt_s);
        let (_, envelope_m) = self.envelope_m(swing_m, velocity, omega);
        (envelope_m / swing_m).min(1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const W: f64 = 1.5e-3;
    const X: f64 = 7e-3;
    const DT: f64 = 1.0 / 720_000.0;
    const OMEGA_500: f64 = core::f64::consts::TAU * 500.0;

    fn limit() -> JetSwingLimit {
        JetSwingLimit::new(W, X, DT).unwrap()
    }

    fn close(a: f64, b: f64) {
        assert!((a / b - 1.0).abs() < 1e-12, "{a} against {b}");
    }

    // Hand values from Python's libm, independent of this code.
    #[test]
    fn the_envelope_is_the_printed_law_at_a_hand_value_and_at_both_ends() {
        // W w/v 2.356, y w/v 15.71: growth reaches W inside the cut-up.
        let (xt, e) = limit().envelope_m(10e-3, 2.0, OMEGA_500);
        close(xt, 3.050_352_518_064_275_6e-3);
        close(e, 2.558_304_853_172_295e-3);
        // W w/v 0.942: v/w is already past W at the exit, so the linear part takes the whole cut-up.
        let (xt, e) = limit().envelope_m(10e-3, 5.0, OMEGA_500);
        assert_eq!(xt, 0.0);
        close(e, 3.375_644_347_017_859e-3);
        // A swing exactly at W: the transition is at the lip and the envelope is W.
        let (xt, e) = limit().envelope_m(W, 0.5, OMEGA_500);
        close(xt, X);
        close(e, W);
    }

    #[test]
    fn no_velocity_at_all_takes_the_width_itself() {
        let (xt, e) = limit().envelope_m(10e-3, 0.0, OMEGA_500);
        assert_eq!(xt, X);
        assert_eq!(e, W);
    }

    #[test]
    fn it_refuses_a_scale_that_is_not_finite_and_positive() {
        for (w, x, dt) in [
            (0.0, X, DT),
            (W, -X, DT),
            (W, X, f64::NAN),
            (f64::INFINITY, X, DT),
        ] {
            assert!(JetSwingLimit::new(w, x, dt).is_err(), "{w} {x} {dt}");
        }
    }

    /// `cycles` cycles of `eta = a sin(w t)` for `eta > 0`, `b sin(w t)` for `eta < 0`, with exit
    /// velocity `eta w / G`. Returns each sample's gain.
    fn drive(limit: &mut JetSwingLimit, a: f64, b: f64, cycles: usize) -> Vec<f64> {
        const GAIN: f64 = 12.0;
        let n = 1440;
        let mut gains = Vec::new();
        for i in 0..cycles * n {
            let phase = core::f64::consts::TAU * ((i % n) as f64 + 0.5) / n as f64;
            let s = phase.sin();
            let eta = if s > 0.0 { a * s } else { b * s };
            gains.push(limit.gain());
            limit.accept(eta, eta * OMEGA_500 / GAIN);
        }
        gains
    }

    #[test]
    fn under_the_width_the_gain_is_exactly_one() {
        let mut l = limit();
        assert!(
            drive(&mut l, 0.99 * W, 0.99 * W, 8)
                .iter()
                .all(|g| *g == 1.0)
        );
    }

    #[test]
    fn the_gain_steps_only_at_a_sign_change_and_a_steady_swing_lands_on_the_envelope() {
        let mut l = limit();
        let gains = drive(&mut l, 10e-3, 10e-3, 8);
        // Half a cycle is 720 samples, and the sign change's own sample still has the old gain: a
        // new gain can start only on the sample after it, index 1 of a half cycle.
        for (i, pair) in gains.windows(2).enumerate() {
            if pair[0] != pair[1] {
                assert_eq!(
                    (i + 1) % 720,
                    1,
                    "gain moved inside a half cycle at {}",
                    i + 1
                );
            }
        }
        // Steady: the same gain on every half cycle, E / y, at this swing's own E. The peak is one
        // half sample off the sine's crest.
        let tail = &gains[gains.len() - 2 * 1440..];
        assert!(tail.iter().all(|g| *g == tail[0]));
        let y = 10e-3 * (core::f64::consts::TAU * 359.5 / 1440.0).sin();
        let (_, e) = l.envelope_m(y, y * OMEGA_500 / 12.0, OMEGA_500);
        assert!(
            (tail[0] * y / e - 1.0).abs() < 1e-12,
            "{} against {}",
            tail[0] * y,
            e
        );
        assert!(tail[0] < 1.0);
    }

    // One gain over the whole cycle: a swing that differs by polarity keeps its own shape.
    #[test]
    fn unequal_half_cycles_share_one_gain_and_the_gain_adds_no_even_line() {
        let mut l = limit();
        let gains = drive(&mut l, 10e-3, 8e-3, 8);
        let tail = &gains[gains.len() - 2 * 1440..];
        assert!(tail.iter().all(|g| *g == tail[0]), "the gain alternates");
        assert!(tail[0] < 1.0, "the limit did not bite: {}", tail[0]);
    }

    // A moved lip is a limit built there: the same envelope, and the next gain from the new lip
    // over the half cycles the old one saw.
    #[test]
    fn a_moved_lip_reads_the_new_width_and_cut_up_and_keeps_its_half_cycles() {
        let (w2, x2) = (1.53e-3, 6.46e-3);
        let mut moved = limit();
        drive(&mut moved, 10e-3, 10e-3, 2);
        let held = moved.clone();
        moved.set_lip(w2, x2).unwrap();
        assert_eq!((moved.spread_width_m(), moved.cut_up_m()), (w2, x2));
        let built = JetSwingLimit::new(w2, x2, DT).unwrap();
        for (y, v) in [(10e-3, 2.0), (10e-3, 5.0), (w2, 0.5), (10e-3, 0.0)] {
            assert_eq!(moved.envelope_m(y, v, OMEGA_500), built.envelope_m(y, v, OMEGA_500));
        }
        assert_eq!(moved.gain(), held.gain(), "the gain holds until the next sign change");
        let mut rebuilt = JetSwingLimit { spread_width_m: w2, cut_up_m: x2, ..held };
        assert_eq!(drive(&mut moved, 10e-3, 10e-3, 2), drive(&mut rebuilt, 10e-3, 10e-3, 2));
        for (w, x) in [(0.0, x2), (w2, f64::NAN), (f64::INFINITY, x2), (w2, -x2)] {
            assert!(moved.clone().set_lip(w, x).is_err(), "{w} {x}");
        }
        let before = moved.clone();
        assert!(moved.set_lip(0.0, x2).is_err());
        assert_eq!(moved, before, "a refused lip leaves the limit unchanged");
    }

    #[test]
    fn exact_zeros_close_no_half_cycle() {
        let mut l = limit();
        drive(&mut l, 10e-3, 10e-3, 2);
        let g = l.gain();
        for _ in 0..1000 {
            l.accept(0.0, 0.0);
        }
        assert_eq!(l.gain(), g);
    }
}
