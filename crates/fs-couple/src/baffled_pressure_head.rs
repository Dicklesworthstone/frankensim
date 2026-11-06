//! A rigid sphere in place of the infinite baffle: the mouth as a cap on a head (opt-in).
//!
//! Geometry, in the mouth-plane frame of [`CircularOutletReceiver`] (origin at the mouth centre, +z out of the mouth): a
//! rigid sphere of radius `a` centred at `(0, 0, -a)`, and the mouth as the spherical cap of the disk's own area `pi b^2`,
//! centred on the pole `(0, 0, 0)`. The cap moves with uniform normal velocity `Q / (pi b^2)`, the rest of the sphere is
//! rigid, and the exterior is the receiver's homogeneous lossless medium. The microphone stands outside the sphere.
//!
//! Law (Morse and Ingard 1968, Theoretical Acoustics, sec. 7.2; time factor `e^{+j w t}`, outgoing `h_n = h_n^(2)`): per
//! unit outward volume flow, at the microphone's distance `r` and polar angle `theta` about the sphere's centre,
//!
//! ```text
//! S(f) = sum_n  -j rho c u_n h_n(k r) / h_n'(k a) P_n(cos theta)
//! u_0 = 1 / (4 pi a^2),   u_n = (P_{n-1}(cos al) - P_{n+1}(cos al)) / (4 pi a^2 (1 - cos al)),   1 - cos al = b^2 / (2 a^2)
//! ```
//!
//! `[Pa s / m^3]`. The Hankel quotients are carried as ratios (`h_n / h_{n-1}` by the upward recurrence, stable for
//! `h^(2)` because its `y_n` part dominates), so no `h_n` itself is formed and nothing overflows at low `k a`. The series
//! stops at `n = k r + 30 + 39 / ln(r / a)`, where a term is below `e^-39` of the leading ones twice over; more than
//! [`MAX_ORDER`] terms is refused.
//!
//! `B(f)` is the baffled disk's pressure per unit flow on the receiver's OWN quadrature, the same points and areas
//! [`BaffledPressure::circular_surface`] gives every circular-outlet image: `B = j w rho sum_s A_s v_s e^{-j k R_s} /
//! (2 pi R_s)`, `v_s = 1 / (pi b^2)`. The head is `G = S / B` applied as a causal FIR to that receiver's pressure. The
//! receiver is causal and linear in the mouth flow, so the composition is the sphere's pressure for the same flow,
//! without a second resampler, derivative or travel kernel. `G` is finite at 0 Hz (`S` and `B` both vanish like `w`
//! there) and causal: for a microphone the mouth sees, the first arrival on the sphere and on the baffle is the same
//! straight path, and everything the sphere adds (its curvature, its creeping waves) arrives later.
//!
//! Band: above the receiver's declared resolved band nothing is a prediction, and `S / B` is singular where the disk's
//! directivity has a null (near `k b sin(theta) = 3.83`). So `G` is blended to 1 by a raised cosine over `[band, 1.25
//! band]`, and above that the baffle's image is kept as it was. A ratio over 4 inside the blend is refused.
//!
//! FIR: `G` on `N = 2^m >= rate / 3` bins (about 3 Hz apart), times `e^{-j w D}` with `D = ceil(rate / 750)` samples
//! (1.33 ms at 48 kHz) so that the short acausal part of the band-limited response falls inside the window; `L = N / 16`
//! taps (21 ms at 48 kHz), the last quarter tapered by a half Hann. `D` is filter latency, reported by
//! [`HeadSphere::latency_samples`], never physical travel. At 48 kHz the FIR matches `G` within 0.005 dB and 0.02 degrees
//! from 50 Hz to the band (tested below).
//!
//! Limits, not hidden: the head is a rigid sphere (no torso, neck, nose, protruding lips or hair); the mouth is a curved
//! cap of the disk's area, not the lips' shape; and the lips' radiation LOAD is not changed here (the tract card's
//! termination still declares it), only the observed pressure.

use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_math::det;

/// Largest series order admitted; a geometry that needs more is refused, never truncated.
pub const MAX_ORDER: usize = 4096;
/// Head FIR latency `ceil(rate / LATENCY_DIVISOR)` output samples (1.33 ms).
const LATENCY_DIVISOR: u32 = 750;
/// Largest admitted `|S / B|` inside the blended band; above it the baffle's directivity null is in the band.
const MAX_RATIO: f64 = 4.0;

/// Minimal complex arithmetic in plain IEEE operations (no fused multiply-add), so every build rounds the same way.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const ONE: Self = Self { re: 1.0, im: 0.0 };
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: Self) -> Self {
        Self::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
    fn scale(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
    fn recip(self) -> Self {
        Self::ONE.div(self)
    }
    fn abs(self) -> f64 {
        det::hypot(self.re, self.im)
    }
    fn finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    /// `e^{-j x}`.
    fn expj_neg(x: f64) -> Self {
        let (s, c) = det::sin_cos(x);
        Self::new(c, -s)
    }
}

/// The sphere's series, prepared once per geometry: Legendre values at the microphone and the cap's coefficients.
#[derive(Clone, Debug)]
struct Series {
    a: f64,
    r: f64,
    /// `P_n(cos theta)` at the microphone.
    p_mic: Vec<f64>,
    /// `u_n`, the cap's normal-velocity coefficients per unit flow [1/m^2].
    u: Vec<f64>,
    medium: RayleighMedium,
}

impl Series {
    fn new(a: f64, b: f64, mic: [f64; 3], medium: RayleighMedium, top_hz: f64, margin: usize) -> Result<Self, String> {
        let rel = [mic[0], mic[1], mic[2] + a];
        let r = det::sqrt(rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]);
        if !(r > a) {
            return Err("the microphone must stand outside the head sphere".into());
        }
        let orders = Self::orders(a, r, top_hz, medium, margin)?;
        let cos_theta = rel[2] / r;
        let one_minus_cos_al = b * b / (2.0 * a * a);
        let cos_al = 1.0 - one_minus_cos_al;
        let p_mic = legendre(cos_theta, orders + 1);
        let p_cap = legendre(cos_al, orders + 2);
        let surface = 4.0 * core::f64::consts::PI * a * a;
        let mut u = Vec::with_capacity(orders + 1);
        u.push(1.0 / surface);
        for n in 1..=orders {
            u.push((p_cap[n - 1] - p_cap[n + 1]) / (surface * one_minus_cos_al));
        }
        Ok(Self { a, r, p_mic, u, medium })
    }

    /// The last order kept at `f_hz`: `k r + margin + 39 / ln(r / a)`.
    fn orders(a: f64, r: f64, f_hz: f64, medium: RayleighMedium, margin: usize) -> Result<usize, String> {
        let k = core::f64::consts::TAU * f_hz / medium.sound_speed;
        let n = k * r + margin as f64 + 39.0 / det::ln(r / a);
        if !n.is_finite() || n > MAX_ORDER as f64 {
            return Err(format!("the head series needs more than {MAX_ORDER} orders here"));
        }
        Ok(n.ceil() as usize)
    }

    /// `S(f)`, the sphere's pressure per unit outward volume flow [Pa s/m^3]; `f_hz > 0`.
    fn pressure(&self, f_hz: f64, margin: usize) -> Result<Cx, String> {
        let k = core::f64::consts::TAU * f_hz / self.medium.sound_speed;
        let (ka, kr) = (k * self.a, k * self.r);
        let last = Self::orders(self.a, self.r, f_hz, self.medium, margin)?.min(self.u.len() - 1);
        // rho_n(x) = h_n(x) / h_{n-1}(x): rho_1 = 1/x + j, rho_{n+1} = (2n+1)/x - 1/rho_n.
        let mut rho_a = Cx::new(1.0 / ka, 1.0);
        let mut rho_r = Cx::new(1.0 / kr, 1.0);
        // q_n = h_n(kr) / h_n(ka), q_0 = (a/r) e^{-jk(r-a)}; d_n = h_n'(ka) / h_n(ka), d_0 = -rho_1(ka).
        let mut q = Cx::expj_neg(k * (self.r - self.a)).scale(self.a / self.r);
        let mut d = Cx::new(-rho_a.re, -rho_a.im);
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut sum = Cx::new(0.0, 0.0);
        for n in 0..=last {
            if n >= 1 {
                if n >= 2 {
                    let m = (2 * n - 1) as f64;
                    rho_a = Cx::new(m / ka, 0.0).sub(rho_a.recip());
                    rho_r = Cx::new(m / kr, 0.0).sub(rho_r.recip());
                }
                q = q.mul(rho_r.div(rho_a));
                d = rho_a.recip().sub(Cx::new((n + 1) as f64 / ka, 0.0));
            }
            // -j rho c u_n (q_n / d_n) P_n(cos theta)
            let t = q.div(d).scale(rho_c * self.u[n] * self.p_mic[n]);
            sum = sum.add(Cx::new(t.im, -t.re));
        }
        if sum.finite() { Ok(sum) } else { Err("the head series is not finite".into()) }
    }
}

/// `P_0..=P_last` at `x` by the three-term recurrence.
fn legendre(x: f64, last: usize) -> Vec<f64> {
    let mut p = Vec::with_capacity(last + 1);
    p.push(1.0);
    if last >= 1 {
        p.push(x);
    }
    for n in 1..last {
        let nf = n as f64;
        p.push(((2.0 * nf + 1.0) * x * p[n] - nf * p[n - 1]) / (nf + 1.0));
    }
    p
}

/// The baffled disk's pressure per unit flow on the receiver's own quadrature, at `f_hz > 0`.
fn baffle(f_hz: f64, surface: &[super::SurfaceSample], mic: [f64; 3], medium: RayleighMedium) -> Cx {
    let k = core::f64::consts::TAU * f_hz / medium.sound_speed;
    let mut sum = Cx::new(0.0, 0.0);
    for s in surface {
        let dx = [mic[0] - s.position_m[0], mic[1] - s.position_m[1], mic[2] - s.position_m[2]];
        let dist = det::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        let w = s.area_m2 * s.mode_shape[0] / (core::f64::consts::TAU * dist);
        sum = sum.add(Cx::expj_neg(k * dist).scale(w));
    }
    // j w rho * sum
    let jw = core::f64::consts::TAU * f_hz * medium.density;
    Cx::new(-sum.im * jw, sum.re * jw)
}

/// The prepared head: its geometry's transfer and the causal FIR that turns the baffle's pressure into the sphere's.
#[derive(Clone, Debug)]
pub struct HeadSphere {
    series: Series,
    surface: Vec<super::SurfaceSample>,
    mic: [f64; 3],
    band_hz: f64,
    rate: f64,
    /// Taps in reverse order, so a window of history oldest-first is one dot product.
    taps_rev: Vec<f64>,
    /// Two copies of the last `L` inputs, so the window is always contiguous.
    history: Vec<f64>,
    cursor: usize,
    latency: usize,
}

impl HeadSphere {
    /// Prepare the head for a circular mouth of radius `mouth_radius_m` observed by `receiver` at `output_rate` Hz.
    ///
    /// # Errors
    /// A non-physical radius, a mouth not smaller than the head, a microphone inside the head, a blend band above
    /// Nyquist, a series beyond [`MAX_ORDER`] or a ratio over 4 in the blend, or any circular-receiver refusal.
    pub fn new(radius_m: f64, mouth_radius_m: f64, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32) -> Result<Self, String> {
        Self::with_margin(radius_m, mouth_radius_m, receiver, medium, output_rate, 30)
    }

    fn with_margin(radius_m: f64, mouth_radius_m: f64, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32, margin: usize) -> Result<Self, String> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err("the head radius must be finite and positive".into());
        }
        if !(mouth_radius_m < radius_m) {
            return Err("the mouth must be smaller than the head".into());
        }
        let rate = f64::from(output_rate);
        let band_hz = receiver.maximum_frequency_hz;
        if !(1.25 * band_hz < 0.5 * rate) {
            return Err("the head's blend band [band, 1.25 band] must lie below Nyquist".into());
        }
        let surface = BaffledPressure::circular_surface(mouth_radius_m, receiver, medium)?;
        let series = Series::new(radius_m, mouth_radius_m, receiver.position_m, medium, 1.25 * band_hz, margin)?;
        let mut head = Self {
            series,
            surface,
            mic: receiver.position_m,
            band_hz,
            rate,
            taps_rev: Vec::new(),
            history: Vec::new(),
            cursor: 0,
            latency: output_rate.div_ceil(LATENCY_DIVISOR) as usize,
        };
        head.design(margin)?;
        Ok(head)
    }

    /// `G(f) = S / B`, blended to 1 over `[band, 1.25 band]` and 1 above; `f_hz >= 0` (0 Hz takes the limit).
    fn ratio(&self, f_hz: f64, df: f64, margin: usize) -> Result<Cx, String> {
        let top = 1.25 * self.band_hz;
        if f_hz >= top {
            return Ok(Cx::ONE);
        }
        // G is finite at 0 Hz, where S and B both vanish like w; a thousandth of a bin stands in for it.
        let f = if f_hz > 0.0 { f_hz } else { 1e-3 * df };
        let g = self.series.pressure(f, margin)?.div(baffle(f, &self.surface, self.mic, self.series.medium));
        if !g.finite() || g.abs() > MAX_RATIO {
            return Err("the head's ratio is singular in its band (the baffle's directivity null is too low)".into());
        }
        if f_hz <= self.band_hz {
            return Ok(g);
        }
        let blend = 0.5 * (1.0 - det::cos(core::f64::consts::PI * (f_hz - self.band_hz) / (top - self.band_hz)));
        Ok(g.scale(1.0 - blend).add(Cx::new(blend, 0.0)))
    }

    /// `G` at a frequency, for reports and tests [`Re`, `Im`].
    ///
    /// # Errors
    /// The series' refusals.
    pub fn transfer(&self, f_hz: f64) -> Result<(f64, f64), String> {
        let g = self.ratio(f_hz, 1.0, 30)?;
        Ok((g.re, g.im))
    }

    fn design(&mut self, margin: usize) -> Result<(), String> {
        let n_fft = ((self.rate / 3.0).ceil() as usize).next_power_of_two();
        let half = n_fft / 2;
        let taps = n_fft / 16;
        let df = self.rate / n_fft as f64;
        let mut bins = Vec::with_capacity(half + 1);
        for k in 0..=half {
            let g = self.ratio(k as f64 * df, df, margin)?;
            // e^{-j 2 pi k D / N}: the latency, so the band-limited response's acausal part is kept.
            let shift = Cx::expj_neg(core::f64::consts::TAU * ((k * self.latency) % n_fft) as f64 / n_fft as f64);
            bins.push(g.mul(shift));
        }
        let twiddle: Vec<(f64, f64)> = (0..n_fft)
            .map(|m| det::sin_cos(core::f64::consts::TAU * m as f64 / n_fft as f64))
            .collect();
        let q = taps / 4;
        let mut fir = Vec::with_capacity(taps);
        for n in 0..taps {
            // g[n] = (1/N) [G_0 + 2 sum_{k=1}^{N/2-1} Re(G_k e^{+j 2 pi k n / N}) + G_{N/2} (-1)^n]
            let mut acc = bins[0].re;
            for (k, bin) in bins.iter().enumerate().take(half).skip(1) {
                let (s, c) = twiddle[(k * n) % n_fft];
                acc += 2.0 * (bin.re * c - bin.im * s);
            }
            acc += if n % 2 == 0 { bins[half].re } else { -bins[half].re };
            let w = if n < taps - q {
                1.0
            } else {
                0.5 * (1.0 + det::cos(core::f64::consts::PI * (n - (taps - q)) as f64 / q as f64))
            };
            fir.push(acc / n_fft as f64 * w);
        }
        fir.reverse();
        self.taps_rev = fir;
        self.history = vec![0.0; 2 * taps];
        self.cursor = 0;
        Ok(())
    }

    /// Head FIR latency in output samples (filter latency, not travel).
    #[must_use]
    pub const fn latency_samples(&self) -> usize {
        self.latency
    }

    /// The FIR's taps in time order.
    #[must_use]
    pub fn taps(&self) -> Vec<f64> {
        self.taps_rev.iter().rev().copied().collect()
    }

    /// One output sample: the baffled receiver's pressure in, the head's pressure out (delayed by the latency).
    pub fn step(&mut self, baffled_pa: f64) -> f64 {
        let l = self.taps_rev.len();
        self.cursor = if self.cursor + 1 == l { 0 } else { self.cursor + 1 };
        self.history[self.cursor] = baffled_pa;
        self.history[self.cursor + l] = baffled_pa;
        super::kernel_tap(&self.taps_rev, &self.history[self.cursor + 1..=self.cursor + l])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const AIR: RayleighMedium = RayleighMedium { density: 1.1336, sound_speed: 356.53 };
    const MIC: [f64; 3] = [0.212_132, 0.0, 0.212_132];

    fn mic(rings: usize, points: usize, band: f64) -> CircularOutletReceiver {
        CircularOutletReceiver { position_m: MIC, radial_rings: rings, angular_points: points, maximum_frequency_hz: band }
    }

    fn g(head: &HeadSphere, f: f64) -> Cx {
        let (re, im) = head.transfer(f).unwrap();
        Cx::new(re, im)
    }

    fn db(x: f64) -> f64 {
        20.0 * x.log10()
    }

    /// A small head is a free-field monopole at its centre: `|G| -> R / (2 r)` (half the baffle's pressure, the
    /// centre's distance), and the remainder, the small sphere's own dipole, shrinks with the radius.
    #[test]
    fn a_small_head_reads_free_field() {
        let dist = |a: f64| -> f64 {
            let head = HeadSphere::new(a, 0.1 * a, mic(4, 16, 8000.0), AIR, 48_000).unwrap();
            let r = (MIC[0] * MIC[0] + (MIC[2] + a) * (MIC[2] + a)).sqrt();
            let rr = (MIC[0] * MIC[0] + MIC[2] * MIC[2]).sqrt();
            let k = core::f64::consts::TAU * 200.0 / AIR.sound_speed;
            let free = Cx::expj_neg(k * (r - rr)).scale(rr / (2.0 * r));
            g(&head, 200.0).sub(free).abs() / free.abs()
        };
        let (mm, tenth) = (dist(1e-3), dist(1e-4));
        assert!(mm < 0.01, "1 mm head off free field by {mm}");
        assert!(tenth < mm / 5.0, "0.1 mm head {tenth} vs 1 mm {mm}");
    }

    /// A large head is the baffle near its pole: `G -> 1` like `1 / a` (scipy on the same quadrature: |G - 1| at 8 kHz
    /// is 3.6e-2, 1.9e-2, 1.0e-2 and 5.2e-3 for 0.25, 0.5, 1 and 2 m).
    #[test]
    fn a_large_head_reads_the_baffle() {
        let off = |a: f64| -> f64 {
            let head = HeadSphere::new(a, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000).unwrap();
            g(&head, 8000.0).sub(Cx::ONE).abs()
        };
        let devs: Vec<f64> = [0.25, 0.5, 1.0, 2.0].into_iter().map(off).collect();
        assert!(devs.windows(2).all(|w| (0.4..0.65).contains(&(w[1] / w[0]))), "{devs:?}");
        assert!(devs[3] < 6e-3, "2 m head at 8 kHz: |G - 1| = {}", devs[3]);
    }

    /// Thirty more orders change nothing a test or a page could read.
    #[test]
    fn the_series_is_converged() {
        let base = HeadSphere::new(0.0875, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000).unwrap();
        let more = HeadSphere::with_margin(0.0875, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000, 60).unwrap();
        for f in [50.0, 700.0, 3000.0, 7900.0] {
            let (a, b) = (base.ratio(f, 1.0, 30).unwrap(), more.ratio(f, 1.0, 60).unwrap());
            assert!(a.sub(b).abs() < 1e-12 * a.abs(), "{f} Hz: {a:?} vs {b:?}");
        }
    }

    /// The sphere's pressure against an independent evaluation (scipy.special spherical_jn/yn, 400 orders) of the
    /// same series at the C3 card's geometry.
    #[test]
    fn the_series_matches_an_independent_evaluation() {
        let s = Series::new(0.0875, 0.013_548, MIC, AIR, 10_000.0, 30).unwrap();
        for (f, re, im) in GOLDEN {
            let p = s.pressure(f, 30).unwrap();
            let want = Cx::new(re, im);
            assert!(p.sub(want).abs() < 1e-9 * want.abs(), "{f} Hz: {p:?} vs {want:?}");
        }
    }

    /// The FIR is the ratio: its response with the latency removed matches `G` from 50 Hz to the band.
    #[test]
    fn the_fir_matches_the_ratio_in_band() {
        let head = HeadSphere::new(0.0875, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000).unwrap();
        assert_eq!(head.latency_samples(), 64);
        let taps = head.taps();
        assert_eq!(taps.len(), 1024);
        let mut f = 50.0;
        while f <= 8000.0 {
            let w = core::f64::consts::TAU * f / 48_000.0;
            let mut h = Cx::new(0.0, 0.0);
            for (n, t) in taps.iter().enumerate() {
                h = h.add(Cx::expj_neg(w * (n as f64 - 64.0)).scale(*t));
            }
            let want = g(&head, f);
            let ratio = h.div(want);
            assert!(db(ratio.abs()).abs() < 0.005, "{f} Hz: {} dB", db(ratio.abs()));
            assert!(ratio.im.atan2(ratio.re).to_degrees().abs() < 0.02, "{f} Hz phase");
            f += 50.0;
        }
    }

    /// The head lowers the lows toward free field and leaves the band's top near the baffle (C3 card geometry).
    #[test]
    fn the_head_tilts_toward_the_baffle_with_frequency() {
        let head = HeadSphere::new(0.0875, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000).unwrap();
        let lo = db(g(&head, 100.0).abs());
        let hi = db(g(&head, 8000.0).abs());
        assert!((-5.6..-4.9).contains(&lo), "{lo}");
        assert!((-0.7..-0.2).contains(&hi), "{hi}");
        assert_eq!(g(&head, 10_000.0), Cx::ONE);
    }

    /// `step` is the FIR, causal: an impulse returns the taps.
    #[test]
    fn step_is_the_causal_fir() {
        let mut head = HeadSphere::new(0.0875, 0.008, mic(4, 16, 4000.0), AIR, 48_000).unwrap();
        let taps = head.taps();
        let out: Vec<f64> = (0..taps.len() + 5).map(|i| head.step(if i == 0 { 1.0 } else { 0.0 })).collect();
        assert_eq!(&out[..taps.len()], &taps[..]);
        assert!(out[taps.len()..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn non_physical_heads_are_refused() {
        let ok = mic(4, 16, 8000.0);
        for (a, b, m, rate, why) in [
            (0.0, 0.01, ok, 48_000, "finite and positive"),
            (f64::NAN, 0.01, ok, 48_000, "finite and positive"),
            (0.01, 0.01, ok, 48_000, "smaller than the head"),
            (0.0875, 0.01, mic(4, 16, 20_000.0), 48_000, "below Nyquist"),
        ] {
            let err = HeadSphere::new(a, b, m, AIR, rate).unwrap_err();
            assert!(err.contains(why), "{err}");
        }
        // Every microphone in front of the mouth plane is outside the sphere; one behind it, inside, is refused.
        let err = Series::new(0.1, 0.01, [0.0, 0.0, -0.05], AIR, 8000.0, 30).unwrap_err();
        assert!(err.contains("outside the head"), "{err}");
    }

    /// `(f, Re S, Im S)` [Hz, Pa s/m^3]: scipy 1.14 spherical_jn/yn, 400 orders, head 8.75 cm, cap of radius 13.548 mm,
    /// microphone (0.212132, 0, 0.212132) m, rho 1.1336, c 356.53.
    const GOLDEN: [(f64, f64, f64); 3] = [
        (300.0, 6.273927074538665e+02, 9.111865394641030e+01),
        (1000.0, -2.569048431744369e+03, 9.166961093269957e+02),
        (3000.0, 9.514481124514730e+01, -9.554082937523743e+03),
    ];
}
