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
//! Limits, not hidden: the head is a rigid sphere (no neck, nose, protruding lips or hair, and no torso unless a
//! [`TorsoSphere`] is given); the mouth is a curved cap of the disk's area, not the lips' shape; and [`HeadSphere`]
//! changes only the observed pressure. The lips' radiation LOAD changes only through the separate, opt-in
//! [`head_mouth_load`] (or [`head_torso_mouth_load`]) below.
//!
//! The torso ([`HeadSphere::with_torso`], [`head_torso_mouth_load`]; the law in the `torso` child module): the snowman
//! of Algazi, Duda and Thompson 2002, a rigid torso sphere coaxially under the head, solved by multipole expansions about
//! both centres. Below `1.25 band` the torso's change is added to `S` (so `G = (S + dS) / B`, blended as before) and to
//! `Z_cap` (weighted 1 up to the band, by the same raised cosine to 0 at `1.25 band`, 0 above), so the microphone and the
//! lips see one geometry over the same band. Without a torso every path is exactly the head's alone.
//!
//! The lips' load ([`head_mouth_load`]): the same cap on the same sphere, seen from the tract. Its self-impedance per
//! unit volume flow (the cap's mean pressure over its flow; Morse and Ingard sec. 7.2, `e^{+j w t}`, `h_n = h_n^(2)`) is
//!
//! ```text
//! Z_cap = (rho c / A) sum_n w_n (-j h_n(k a) / h_n'(k a)),   A = pi b^2 = 2 pi a^2 (1 - cos al)
//! w_n = (P_{n-1}(cos al) - P_{n+1}(cos al))^2 / (2 (2n+1) (1 - cos al)),   P_{-1} = 1,   sum_n w_n = 1
//! ```
//!
//! with the same ratio recurrences as `S`. The whole sphere pulsating (`al = pi`) gives `rho c j k a / (1 + j k a) / A`;
//! as `k a -> 0` the resistance tends to `rho c k^2 / (4 pi)`, HALF the baffled piston's, because the sphere radiates into
//! the whole space; a large sphere at a fixed mouth gives the baffled piston, its difference falling like `1 / a`. The
//! terms fall about as `n^-3` past `n = k a` (the cap's rim is a velocity edge): the tail after `N` terms is a mass
//! reactance of about `k a / (pi al N^2)` in units of `rho c / A`, so the sums to `N` and `2 N` differ by `3/4` of it.
//! Relative to the whole that is about `1 / (3.2 al^2 N^2)` at small `k a` and `0.75 k a / (pi al N^2)` once `Z_cap` is
//! near `rho c / A`, which grows with `k a`. The sum runs to `N = max(2 k a + 20, 1000 a / b, 1000 sqrt(2 k a a / (pi
//! b)))` and to `2 N` (`b / a <= al`), refuses if the two differ by more than `1e-6` relative, and refuses a geometry
//! needing more than [`LOAD_MAX_ORDER`] terms. The lips' load is then
//!
//! ```text
//! Z_L = Z_flanged + (conj(Z_cap) - Z_piston)          (fs-duct's e^{-i w t}; the conjugate carries Z_cap over)
//! ```
//!
//! `Z_flanged` is the declared flanged pipe's own multimodal load across the gas interface
//! (`fs_duct::flanged_pipe_impedance_with_exterior` with `fs_duct::FLANGED_MODES`): its radiation in the exterior gas,
//! its eliminated higher duct modes in the gas inside the duct (the tract's; a voice card's `tract_gas`, or the ambient
//! without one). With one gas it is `fs_duct::compact_radiation_impedance`'s flanged load bit for bit. Everything else in
//! the table is in the exterior gas: the cap, the piston, the torso and the mouth record. `Z_piston` is
//! the closed-form baffled piston of the same radius (`fs_phs::baffled_piston_impedance`). Both added terms have uniform
//! velocity, so their singular near fields cancel and what is added is the smooth change in the exterior: the
//! flanged load's own velocity-profile part stays. As the sphere grows, `Z_L` returns to `Z_flanged`. The table's rows
//! are exactly the bins of the finest line rung (`w_k = 2 pi k / (N_max dt)`, `k = 1..=N_max / 2`), so every rung of a
//! characteristic or mouth-flow line reads its own bin and nothing is interpolated there; the DC limit is 0, as for any
//! exterior load.
//!
//! Openings off the pole ([`HeadSphere::for_opening`], [`head_opening_load`]; the nasal rows' nostrils, N3b). An
//! opening of radius `b_i` at `(x_i, y_i)` on the baffle becomes the cap of the disk's area `pi b_i^2` centred on the ray
//! from the head's centre `(0, 0, -a)` through `(x_i, y_i, 0)`: the cap's centre is `a e_i` from the head's centre,
//! `e_i = (x_i, y_i, a) / |(x_i, y_i, a)|`. That placement is a modelling choice (no print places a nostril on a
//! sphere), as the lips' cap at the pole is: the cap's centre is `|(x_i, y_i, a)| - a` from the opening's place along
//! the ray and `a (1 - e_z)` behind the mouth's plane; for the C3-a card's nostrils at `(+-7.07, 27.63)` mm on an 8.7 cm
//! head, 17.56 degrees up, 4.65 degrees aside and 4.33 mm behind the plane. `S_i` is the series
//! above with `cos theta` measured from `e_i`; `B_i` is the baffled disk on the opening's own quadrature, the points and
//! areas the opening's receiver uses ([`BaffledPressure::resolved_circular_grids`]), moved to its centre; `G_i = S_i /
//! B_i` is designed exactly as `G`, with the same band, blend, latency and length. The openings' pressures are summed
//! after their FIRs: the openings do not load one another, and no cap scatters another's field. With a torso the cap's
//! change takes `P_n^m(cos theta_c)` (the `torso` child module). The cap's own load `Z_cap(b_i)` does not depend on where
//! the cap sits (the sphere is the same in every direction); only the torso's change does. An opening at `(0, 0)` is
//! the lips' cap at the pole, and every path above is taken unchanged. An opening's flanged term takes the gas inside the
//! path that ends at it: the nasal rows fill every edge with the tract's gas, so a nostril's eliminated duct modes are in
//! the tract's gas, as the lips' are.

use super::openings::CircularOpening;
use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_duct::{Duct, Segment, TabulatedLoad, TabulatedMouth};
use fs_material::gas::GasState;
use fs_math::det;

/// The snowman's torso under the head (module doc, "The torso").
#[path = "baffled_pressure_head_torso.rs"]
mod torso;
pub use torso::{TORSO_MAX_ORDER, TorsoSphere};

/// Largest series order admitted; a geometry that needs more is refused, never truncated.
pub const MAX_ORDER: usize = 4096;
/// Largest order of the cap's self-impedance series ([`head_mouth_load`]): its terms fall only as `n^-3`, so it needs
/// `2 max(1000 a / b, 1000 sqrt(2 k a a / (pi b)))` of them at the table's top; a mouth so small against the head, or a
/// table so high, that it needs more is refused, never truncated. It is `2^17` so that a nostril's cap (`b` about 5 mm on
/// an 8.7 cm head) is admitted on the 262144-tap bins at 192 kHz (about 80,500 and 86,100 orders for Dang's subject 2);
/// the convergence rule is unchanged.
pub const LOAD_MAX_ORDER: usize = 1 << 17;
/// Largest admitted change of the cap's self-impedance between `N` and `2 N` terms, relative.
const LOAD_SERIES_TOL: f64 = 1.0e-6;
/// Head FIR latency `ceil(rate / LATENCY_DIVISOR)` output samples (1.33 ms).
const LATENCY_DIVISOR: u32 = 750;
/// Largest admitted `|S / B|` inside the blended band; above it the baffle's directivity null is in the band.
const MAX_RATIO: f64 = 4.0;
/// The lips' cap axis: from the head's centre through the mouth's centre, `+z`.
const POLE: [f64; 3] = [0.0, 0.0, 1.0];

/// The cap axis of an opening centred at `center_m` on the baffle of a head of radius `a` (module doc, "Openings off the
/// pole"): exactly [`POLE`] at the origin.
fn cap_axis(a: f64, center_m: [f64; 2]) -> [f64; 3] {
    if center_m == [0.0, 0.0] {
        return POLE;
    }
    let [x, y] = center_m;
    let n = det::sqrt(x * x + y * y + a * a);
    [x / n, y / n, a / n]
}

/// Minimal complex arithmetic in plain IEEE operations (no fused multiply-add), so every build rounds the same way.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
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
    /// The series for a cap of radius `b` centred on `axis` (a unit vector from the head's centre; [`POLE`] for the lips).
    fn new(a: f64, b: f64, axis: [f64; 3], mic: [f64; 3], medium: RayleighMedium, top_hz: f64, margin: usize)
        -> Result<Self, String> {
        let rel = [mic[0], mic[1], mic[2] + a];
        let r = det::sqrt(rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]);
        if !(r > a) {
            return Err("the microphone must stand outside the head sphere".into());
        }
        let orders = Self::orders(a, r, top_hz, medium, margin)?;
        // The microphone's angle from the cap's centre; the lips' exactly as before.
        let cos_theta = if axis == POLE {
            rel[2] / r
        } else {
            (rel[0] * axis[0] + rel[1] * axis[1] + rel[2] * axis[2]) / r
        };
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
    /// The snowman's torso (module doc, "The torso"); `None` is the head alone.
    torso: Option<torso::Snowman>,
}

impl HeadSphere {
    /// Prepare the head for a circular mouth of radius `mouth_radius_m` observed by `receiver` at `output_rate` Hz.
    ///
    /// # Errors
    /// A non-physical radius, a mouth not smaller than the head, a microphone inside the head, a blend band above
    /// Nyquist, a series beyond [`MAX_ORDER`] or a ratio over 4 in the blend, or any circular-receiver refusal.
    pub fn new(radius_m: f64, mouth_radius_m: f64, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32) -> Result<Self, String> {
        Self::with_margin(radius_m, mouth_radius_m, receiver, medium, output_rate, 30, None)
    }

    /// The head with the snowman's torso under it (module doc, "The torso"): the torso's change joins `S` below
    /// `1.25 band`.
    ///
    /// # Errors
    /// Every refusal of [`HeadSphere::new`], a torso refusal ([`TorsoSphere::validate`]), a microphone inside the
    /// torso, or a band needing more than [`TORSO_MAX_ORDER`] orders.
    pub fn with_torso(radius_m: f64, mouth_radius_m: f64, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32, torso: TorsoSphere) -> Result<Self, String> {
        Self::with_margin(radius_m, mouth_radius_m, receiver, medium, output_rate, 30, Some(torso))
    }

    /// The head for one opening of the baffle (module doc, "Openings off the pole"): its cap on the ray from the head's
    /// centre through `opening`'s centre, `B` on the opening's own quadrature, the torso's change at that cap when
    /// `torso` is given. An opening at `(0, 0)` whose declared grid resolves the band is [`HeadSphere::new`] (or
    /// [`HeadSphere::with_torso`]) bit for bit.
    ///
    /// # Errors
    /// Every refusal of [`HeadSphere::new`] and [`HeadSphere::with_torso`], or an opening array refusal.
    pub fn for_opening(radius_m: f64, opening: CircularOpening, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32, torso: Option<TorsoSphere>) -> Result<Self, String> {
        Self::admit(radius_m, opening.radius_m, receiver, output_rate)?;
        let (rings, points) = BaffledPressure::resolved_circular_grids(&[opening], receiver, medium)?[0];
        let local = CircularOutletReceiver { radial_rings: rings, angular_points: points, ..receiver };
        let mut surface = BaffledPressure::circular_surface(opening.radius_m, local, medium)?;
        for sample in &mut surface {
            sample.position_m[0] += opening.center_m[0];
            sample.position_m[1] += opening.center_m[1];
        }
        let axis = cap_axis(radius_m, opening.center_m);
        Self::on_surface(radius_m, opening.radius_m, axis, surface, receiver, medium, output_rate, 30, torso)
    }

    fn with_margin(radius_m: f64, mouth_radius_m: f64, receiver: CircularOutletReceiver, medium: RayleighMedium,
        output_rate: u32, margin: usize, torso: Option<TorsoSphere>) -> Result<Self, String> {
        // The head's own refusals come before the receiver's, as they always have.
        Self::admit(radius_m, mouth_radius_m, receiver, output_rate)?;
        let surface = BaffledPressure::circular_surface(mouth_radius_m, receiver, medium)?;
        Self::on_surface(radius_m, mouth_radius_m, POLE, surface, receiver, medium, output_rate, margin, torso)
    }

    fn admit(radius_m: f64, cap_radius_m: f64, receiver: CircularOutletReceiver, output_rate: u32) -> Result<(), String> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err("the head radius must be finite and positive".into());
        }
        if !(cap_radius_m < radius_m) {
            return Err("the mouth must be smaller than the head".into());
        }
        if !(1.25 * receiver.maximum_frequency_hz < 0.5 * f64::from(output_rate)) {
            return Err("the head's blend band [band, 1.25 band] must lie below Nyquist".into());
        }
        Ok(())
    }

    /// The head for a cap of radius `cap_radius_m` on `axis`, against the baffled disk sampled at `surface`.
    #[allow(clippy::too_many_arguments)]
    fn on_surface(radius_m: f64, cap_radius_m: f64, axis: [f64; 3], surface: Vec<super::SurfaceSample>,
        receiver: CircularOutletReceiver, medium: RayleighMedium, output_rate: u32, margin: usize,
        torso: Option<TorsoSphere>) -> Result<Self, String> {
        Self::admit(radius_m, cap_radius_m, receiver, output_rate)?;
        let rate = f64::from(output_rate);
        let band_hz = receiver.maximum_frequency_hz;
        let series = Series::new(radius_m, cap_radius_m, axis, receiver.position_m, medium, 1.25 * band_hz, margin)?;
        let torso = match torso {
            Some(t) => Some(torso::Snowman::at_cap(radius_m, cap_radius_m, axis, t, Some(receiver.position_m),
                medium.density, medium.sound_speed, 1.25 * band_hz)?),
            None => None,
        };
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
            torso,
        };
        head.design(margin)?;
        Ok(head)
    }

    /// The frequency a bin is evaluated at: `G` is finite at 0 Hz, where `S` and `B` both vanish like `w`, and a
    /// thousandth of a bin stands in for it.
    fn probe(f_hz: f64, df: f64) -> f64 {
        if f_hz > 0.0 { f_hz } else { 1e-3 * df }
    }

    /// `G(f) = (S + dS) / B` (`dS` the torso's change, when there is one), blended to 1 over `[band, 1.25 band]` and 1
    /// above; `f_hz >= 0` (0 Hz takes the limit).
    fn ratio(&self, f_hz: f64, df: f64, margin: usize) -> Result<Cx, String> {
        let ds = match &self.torso {
            Some(t) if f_hz < 1.25 * self.band_hz => Some(t.change(Self::probe(f_hz, df))?.0),
            _ => None,
        };
        self.ratio_with(f_hz, df, margin, ds)
    }

    /// [`Self::ratio`] with the torso's change `ds` already solved (or none: the head alone, exactly as before).
    fn ratio_with(&self, f_hz: f64, df: f64, margin: usize, ds: Option<Cx>) -> Result<Cx, String> {
        let top = 1.25 * self.band_hz;
        if f_hz >= top {
            return Ok(Cx::ONE);
        }
        let f = Self::probe(f_hz, df);
        let s = self.series.pressure(f, margin)?;
        let s = match ds {
            Some(d) => s.add(d),
            None => s,
        };
        let g = s.div(baffle(f, &self.surface, self.mic, self.series.medium));
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
        // The torso's change on every bin below 1.25 band (a prefix of the bins), each bin solved alone.
        let ds: Vec<Cx> = match &self.torso {
            Some(t) => {
                let top = 1.25 * self.band_hz;
                let n = (0..=half).take_while(|&k| (k as f64 * df) < top).count();
                torso::par_map(n, |k| Ok(t.change(Self::probe(k as f64 * df, df))?.0))?
            }
            None => Vec::new(),
        };
        let mut bins = Vec::with_capacity(half + 1);
        for k in 0..=half {
            let g = self.ratio_with(k as f64 * df, df, margin, ds.get(k).copied())?;
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

/// The cap's self-impedance series (module doc, "The lips' load"), prepared once per geometry.
#[derive(Clone, Debug)]
struct CapLoad {
    a: f64,
    sound_speed: f64,
    /// `rho c / A` [Pa s/m^3].
    scale: f64,
    /// `1000 a / b`: the first sum's length never falls below it (the small-`k a` tail).
    floor: f64,
    /// `2 a / (pi b)`: the first sum runs to at least `1000 sqrt(edge k a)` (the large-`k a` tail).
    edge: f64,
    /// `w_0..=w_last`.
    w: Vec<f64>,
}

impl CapLoad {
    /// The first sum's length `N` at `x = k a` (module doc, "The lips' load"); the second runs to `2 N`.
    fn orders(floor: f64, edge: f64, x: f64) -> f64 {
        (2.0 * x + 20.0).max(floor).max(1000.0 * (edge * x).sqrt()).ceil()
    }

    /// A cap of the disk's area `pi b^2` on a sphere of radius `a`, for angular frequencies up to `omega_top`.
    fn new(a: f64, b: f64, density: f64, sound_speed: f64, omega_top: f64) -> Result<Self, String> {
        let one_minus_cos_al = b * b / (2.0 * a * a);
        let floor = 1000.0 * a / b;
        let edge = 2.0 * a / (core::f64::consts::PI * b);
        let last = 2.0 * Self::orders(floor, edge, omega_top / sound_speed * a);
        if !(last.is_finite() && last <= LOAD_MAX_ORDER as f64) {
            return Err(format!(
                "the head's load series needs more than {LOAD_MAX_ORDER} orders here (mouth radius {b} m on a head of {a} m)"
            ));
        }
        let last = last as usize;
        let p = legendre(1.0 - one_minus_cos_al, last + 1);
        let mut w = Vec::with_capacity(last + 1);
        w.push(0.5 * one_minus_cos_al);
        for n in 1..=last {
            let d = p[n - 1] - p[n + 1];
            w.push(d * d / (2.0 * (2 * n + 1) as f64 * one_minus_cos_al));
        }
        let area = core::f64::consts::PI * b * b;
        Ok(Self { a, sound_speed, scale: density * sound_speed / area, floor, edge, w })
    }

    /// `Z_cap` at `omega > 0` under `e^{+j w t}` [Pa s/m^3]: the sums to `N` and `2 N` terms, refused unless they agree
    /// within [`LOAD_SERIES_TOL`].
    fn impedance(&self, omega: f64) -> Result<Cx, String> {
        let x = omega / self.sound_speed * self.a;
        let n1 = Self::orders(self.floor, self.edge, x);
        if !(x > 0.0 && 2.0 * n1 < self.w.len() as f64) {
            return Err(format!("the head's load series was not prepared for {omega} rad/s"));
        }
        let n1 = n1 as usize;
        let n2 = 2 * n1;
        // rho_n = h_n / h_{n-1}, rho_1 = 1/x + j. Order 0: -j h_0 / h_0' = j h_0 / h_1 = j / rho_1.
        let mut rho = Cx::new(1.0 / x, 1.0);
        let inv = rho.recip();
        let mut sum = Cx::new(-inv.im, inv.re).scale(self.w[0]);
        let mut first = None;
        for n in 1..=n2 {
            let inv = rho.recip(); // h_{n-1} / h_n
            let d = inv.sub(Cx::new((n + 1) as f64 / x, 0.0)); // h_n' / h_n
            let t = d.recip();
            sum = sum.add(Cx::new(t.im, -t.re).scale(self.w[n])); // -j h_n / h_n'
            if n == n1 {
                first = Some(sum);
            }
            rho = Cx::new((2 * n + 1) as f64 / x, 0.0).sub(inv); // rho_{n+1}
        }
        let first = first.ok_or("the head's load series has no first sum")?;
        if !sum.finite() || first.sub(sum).abs() > LOAD_SERIES_TOL * sum.abs() {
            return Err(format!("the head's load series did not converge at {omega} rad/s"));
        }
        Ok(sum.scale(self.scale))
    }
}

/// The lips' radiation load with the head (module doc, "The lips' load"): `Z_flanged + (conj(Z_cap) - Z_piston)`
/// [Pa s/m^3, `e^{-i w t}`] for a cap of the mouth's area on a rigid sphere of radius `radius_m`, on the bins
/// `w_k = 2 pi k / (max_taps dt)`, `k = 1..=max_taps / 2`, of the finest rung of the lines that will read it, bound to
/// the duct's own mouth plane and radius and to `gas`, the exterior, with DC 0. `tract` is the gas inside the duct: only
/// `Z_flanged`'s eliminated duct modes take it (module doc); with `tract` equal to `gas` the table is the one-gas table.
///
/// # Errors
/// A non-physical head, an empty duct or a mouth not smaller than the head, a zero rate or a tap count that is not a
/// power of two of at least 2, a series beyond [`LOAD_MAX_ORDER`] or not converged, a flanged refusal (unequal static
/// pressures, or `ka` above the flanged ceiling in either gas) or a piston refusal, or a row the table refuses (not
/// finite or not passive within 1e-12).
pub fn head_mouth_load(radius_m: f64, duct: &Duct, tract: &GasState, gas: &GasState, sample_rate_hz: u32,
    max_taps: usize) -> Result<TabulatedLoad, String> {
    mouth_load(radius_m, None, duct, tract, gas, sample_rate_hz, max_taps)
}

/// [`head_mouth_load`] with the snowman's torso under the head (module doc, "The torso"): the cap's self-impedance gains
/// the torso's change `dZ`, weighted 1 up to `band_hz` (the microphone's resolved band), by a raised cosine to 0 at
/// `1.25 band_hz`, and 0 above, so the lips see the geometry the microphone does over the same band. Each weighted row is
/// a convex mix of two passive self-impedances, so the table stays passive.
///
/// # Errors
/// Every refusal of [`head_mouth_load`], a torso refusal, a band that is not finite and positive, or a band needing more
/// than [`TORSO_MAX_ORDER`] orders.
#[allow(clippy::too_many_arguments)]
pub fn head_torso_mouth_load(radius_m: f64, torso: TorsoSphere, band_hz: f64, duct: &Duct, tract: &GasState,
    gas: &GasState, sample_rate_hz: u32, max_taps: usize) -> Result<TabulatedLoad, String> {
    mouth_load(radius_m, Some((torso, band_hz)), duct, tract, gas, sample_rate_hz, max_taps)
}

/// The radiation load of one opening of the baffle with the head (module doc, "Openings off the pole"): the
/// [`head_mouth_load`] law for a cap of the opening's area centred on the ray from the head's centre through the
/// opening's centre, with the torso's change at that cap (weighted as [`head_torso_mouth_load`] weights it) when `torso`
/// is given as `(torso, band_hz)`; on the same bins, bound to `axial_position_m` (the opening's plane, from the inlet of
/// the path that reads it), the opening's radius and `gas`, the exterior, with DC 0; `tract` is the gas inside that
/// path, which only `Z_flanged`'s eliminated duct modes take. An opening at `(0, 0)` of a duct's mouth radius, bound to
/// that duct's length, is [`head_mouth_load`] (or [`head_torso_mouth_load`]) bit for bit.
///
/// # Errors
/// Every refusal of [`head_mouth_load`] and [`head_torso_mouth_load`], or a binding the table refuses.
#[allow(clippy::too_many_arguments)]
pub fn head_opening_load(radius_m: f64, opening: CircularOpening, torso: Option<(TorsoSphere, f64)>,
    axial_position_m: f64, tract: &GasState, gas: &GasState, sample_rate_hz: u32, max_taps: usize)
    -> Result<TabulatedLoad, String> {
    if !radius_m.is_finite() || radius_m <= 0.0 {
        return Err("the head radius must be finite and positive".into());
    }
    if !opening.center_m.iter().all(|x| x.is_finite()) {
        return Err("the opening's centre must be finite".into());
    }
    let axis = cap_axis(radius_m, opening.center_m);
    cap_load(radius_m, opening.radius_m, axis, torso, axial_position_m, tract, gas, sample_rate_hz, max_taps)
}

/// The lips' declared flanged load with no head, for `tract_gas` without `head_load` (MB 23833, 24424): fs-duct's
/// flanged load across a stationary gas interface ([`fs_duct::flanged_pipe_impedance_with_exterior`] with
/// [`fs_duct::FLANGED_MODES`]), whose radiation is in `exterior` and whose eliminated higher duct modes are in `tract`,
/// for the duct's mouth radius, on the bins of [`head_mouth_load`], bound to the duct's own mouth plane and radius and
/// to `exterior` with DC 0. Every row is that law at its own bin; with one gas (`tract` equal to `exterior`) it is
/// [`fs_duct::compact_radiation_impedance`]'s flanged load bit for bit.
///
/// # Errors
/// An empty duct, a zero rate or a tap count that is not a power of two of at least 2, a flanged refusal (unequal
/// static pressures, or `ka` above the flanged ceiling in either gas), or a row the table refuses (not finite or not
/// passive within 1e-12).
pub fn flanged_mouth_load(duct: &Duct, tract: &GasState, exterior: &GasState, sample_rate_hz: u32, max_taps: usize)
    -> Result<TabulatedLoad, String> {
    let mouth = duct.segments.last().ok_or("the flanged load needs a duct with segments")?.outlet_radius();
    let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
    if sample_rate_hz == 0 || max_taps < 2 || !max_taps.is_power_of_two() {
        return Err("the flanged load needs a positive rate and a power-of-two tap count of at least 2".into());
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    // `cap_load`'s bins, the lines' own expression, so each rung's query lands on its row bit for bit.
    let omegas: Vec<f64> =
        (1..=max_taps / 2).map(|k| core::f64::consts::TAU * k as f64 / (max_taps as f64 * dt)).collect();
    // Each row is solved alone (the first refusal in row order is returned): the same table on any thread count.
    let rows = torso::par_map(omegas.len(), |i| {
        let omega = omegas[i];
        fs_duct::flanged_pipe_impedance_with_exterior(tract, exterior, mouth, omega, fs_duct::FLANGED_MODES)
            .map(|z| (z.re, z.im))
            .map_err(|e| format!("the flanged load refused at {omega} rad/s: {e:?}"))
    })?;
    let z: Vec<fs_math::c64::C64> = rows.iter().map(|&(re, im)| fs_math::c64::C64::new(re, im)).collect();
    let source = "flanged across the gas interface (fs-duct flanged_pipe_impedance_with_exterior)";
    TabulatedLoad::try_new(omegas, z, source, 1.0e-12)
        .and_then(|table| {
            table.with_mouth(TabulatedMouth {
                axial_position_m: length, radius_m: mouth, gas: *exterior, dc_impedance: 0.0,
            })
        })
        .map_err(|e| format!("the flanged load table refused: {e:?}"))
}

#[allow(clippy::too_many_arguments)]
fn mouth_load(radius_m: f64, torso: Option<(TorsoSphere, f64)>, duct: &Duct, tract: &GasState, gas: &GasState,
    sample_rate_hz: u32, max_taps: usize) -> Result<TabulatedLoad, String> {
    if !radius_m.is_finite() || radius_m <= 0.0 {
        return Err("the head radius must be finite and positive".into());
    }
    let mouth = duct.segments.last().ok_or("the head's load needs a duct with segments")?.outlet_radius();
    let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
    cap_load(radius_m, mouth, POLE, torso, length, tract, gas, sample_rate_hz, max_taps)
}

/// The table of a cap of radius `mouth` on `axis` (module doc, "The lips' load" and "Openings off the pole"): `gas` is
/// the exterior, `tract` the gas inside the path, which only the flanged term's eliminated duct modes take.
#[allow(clippy::too_many_arguments)]
fn cap_load(radius_m: f64, mouth: f64, axis: [f64; 3], torso: Option<(TorsoSphere, f64)>, length: f64, tract: &GasState,
    gas: &GasState, sample_rate_hz: u32, max_taps: usize) -> Result<TabulatedLoad, String> {
    if !(mouth > 0.0 && mouth < radius_m) {
        return Err("the mouth must be smaller than the head".into());
    }
    if sample_rate_hz == 0 || max_taps < 2 || !max_taps.is_power_of_two() {
        return Err("the head's load needs a positive rate and a power-of-two tap count of at least 2".into());
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    // The lines' own bin expression, so each rung's query lands on its row bit for bit.
    let omegas: Vec<f64> =
        (1..=max_taps / 2).map(|k| core::f64::consts::TAU * k as f64 / (max_taps as f64 * dt)).collect();
    let top = omegas[omegas.len() - 1];
    let cap = CapLoad::new(radius_m, mouth, gas.density, gas.sound_speed, top)?;
    // The torso's change on every row below 1.25 band (a prefix of the rows), each row solved alone, and its weight.
    let (dz, weights): (Vec<Cx>, Vec<f64>) = match torso {
        Some((t, band)) => {
            if !(band.is_finite() && band > 0.0) {
                return Err("the torso's band must be finite and positive".into());
            }
            let top_hz = 1.25 * band;
            let snow = torso::Snowman::at_cap(radius_m, mouth, axis, t, None, gas.density, gas.sound_speed, top_hz)?;
            let n = omegas.iter().take_while(|&&w| w / core::f64::consts::TAU < top_hz).count();
            let dz = torso::par_map(n, |i| Ok(snow.change(omegas[i] / core::f64::consts::TAU)?.1))?;
            let weights = omegas[..n]
                .iter()
                .map(|&w| {
                    let f = w / core::f64::consts::TAU;
                    if f <= band { 1.0 } else { 0.5 * (1.0 + det::cos(core::f64::consts::PI * (f - band) / (top_hz - band))) }
                })
                .collect();
            (dz, weights)
        }
        None => (Vec::new(), Vec::new()),
    };
    // Each row is solved alone (the first refusal in row order is returned), so the table is the same on any thread count.
    let rows = torso::par_map(omegas.len(), |i| {
        let omega = omegas[i];
        // Radiation in the exterior, the eliminated duct modes in the tract (module doc); one gas gives the one-gas load.
        let flanged = fs_duct::flanged_pipe_impedance_with_exterior(tract, gas, mouth, omega, fs_duct::FLANGED_MODES)
            .map_err(|e| format!("the flanged load refused at {omega} rad/s: {e:?}"))?;
        let (r_piston, x_piston) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, mouth, omega)
            .map_err(|e| format!("the baffled piston refused at {omega} rad/s: {e:?}"))?;
        let c = cap.impedance(omega)?;
        let c = match (dz.get(i), weights.get(i)) {
            (Some(d), Some(&w)) => c.add(d.scale(w)),
            _ => c,
        };
        // conj(Z_cap) carries the series' e^{+j w t} over to fs-duct's e^{-i w t}.
        Ok(Cx::new(flanged.re + (c.re - r_piston), flanged.im + (-c.im - x_piston)))
    })?;
    let z: Vec<fs_math::c64::C64> = rows.iter().map(|r| fs_math::c64::C64::new(r.re, r.im)).collect();
    let cap_place = if axis == POLE {
        String::new()
    } else {
        format!(" cap on e=({:e},{:e},{:e}) (central projection, a modelling choice)", axis[0], axis[1], axis[2])
    };
    let source = match torso {
        None => format!(
            "head-cap-on-sphere a={radius_m:e} b={mouth:e}{cap_place} (morse-ingard-1968-7.2) + flanged - baffled piston"
        ),
        Some((t, band)) => format!(
            "head-cap-on-sphere a={radius_m:e} b={mouth:e}{cap_place} (morse-ingard-1968-7.2) + torso sphere b={:e} h={:e} \
             below {band:e} Hz (algazi-duda-thompson-2002, multipole) + flanged - baffled piston",
            t.radius_m, t.neck_gap_m
        ),
    };
    TabulatedLoad::try_new(omegas, z, &source, 1.0e-12)
        .and_then(|table| {
            table.with_mouth(TabulatedMouth { axial_position_m: length, radius_m: mouth, gas: *gas, dc_impedance: 0.0 })
        })
        .map_err(|e| format!("the head's load table refused: {e:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_duct::Termination;

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
        let more = HeadSphere::with_margin(0.0875, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000, 60, None).unwrap();
        for f in [50.0, 700.0, 3000.0, 7900.0] {
            let (a, b) = (base.ratio(f, 1.0, 30).unwrap(), more.ratio(f, 1.0, 60).unwrap());
            assert!(a.sub(b).abs() < 1e-12 * a.abs(), "{f} Hz: {a:?} vs {b:?}");
        }
    }

    /// The sphere's pressure against an independent evaluation (scipy.special spherical_jn/yn, 400 orders) of the
    /// same series at the C3 card's geometry.
    #[test]
    fn the_series_matches_an_independent_evaluation() {
        let s = Series::new(0.0875, 0.013_548, POLE, MIC, AIR, 10_000.0, 30).unwrap();
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

    /// The print's KEMAR torso (Algazi, Duda and Thompson 2002).
    const KEMAR: TorsoSphere = TorsoSphere { radius_m: 0.169, neck_gap_m: 0.053 };

    /// With the torso the FIR is still its ratio (the design's precomputed change against the ratio's own, 50 Hz to the
    /// band), the latency is the head's, the ratio is 1 above `1.25 band`, and the torso moves the C3 microphone by at
    /// most 1.847 dB at 50 Hz steps, at 1 kHz (numpy prototype `torso/proto.py`, against the same head alone).
    #[test]
    fn the_torso_joins_the_ratio_and_the_fir_is_still_the_ratio() {
        let head = HeadSphere::new(0.087, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000).unwrap();
        let snow = HeadSphere::with_torso(0.087, 0.013_548, mic(16, 64, 8000.0), AIR, 48_000, KEMAR).unwrap();
        assert_eq!(snow.latency_samples(), head.latency_samples());
        let taps = snow.taps();
        let mut worst = (0.0, 0.0_f64);
        let mut f = 50.0;
        while f <= 8000.0 {
            let w = core::f64::consts::TAU * f / 48_000.0;
            let mut h = Cx::new(0.0, 0.0);
            for (n, t) in taps.iter().enumerate() {
                h = h.add(Cx::expj_neg(w * (n as f64 - 64.0)).scale(*t));
            }
            let want = g(&snow, f);
            let ratio = h.div(want);
            assert!(db(ratio.abs()).abs() < 0.005, "{f} Hz: {} dB", db(ratio.abs()));
            assert!(ratio.im.atan2(ratio.re).to_degrees().abs() < 0.02, "{f} Hz phase");
            let moved = db(want.abs() / g(&head, f).abs());
            if moved.abs() > worst.1.abs() {
                worst = (f, moved);
            }
            f += 50.0;
        }
        assert!(worst.0 == 1000.0 && (worst.1 + 1.847).abs() < 0.005, "{worst:?}");
        assert_eq!(g(&snow, 10_000.0), Cx::ONE);
    }

    /// The lips' load with the torso: each row below `1.25 band` is the head's plus the snowman's `dZ` (conjugated to
    /// fs-duct's time factor) times the raised-cosine weight, every row above is the head's bit for bit, and a torso
    /// 500 m down gives the head's table back.
    #[test]
    fn the_torso_adds_its_change_to_the_lips_load() {
        let gas = air_state();
        let (b, band) = (0.013_548, 4000.0);
        let head = head_mouth_load(0.087, &tube(b), &gas, &gas, 48_000, 4096).unwrap();
        let snow = head_torso_mouth_load(0.087, KEMAR, band, &tube(b), &gas, &gas, 48_000, 4096).unwrap();
        let far = head_torso_mouth_load(0.087, TorsoSphere { radius_m: 0.169, neck_gap_m: 500.0 }, band, &tube(b), &gas,
            &gas, 48_000, 4096).unwrap();
        let s = torso::Snowman::new(0.087, b, KEMAR, None, gas.density, gas.sound_speed, 1.25 * band).unwrap();
        let dt = 1.0 / 48_000.0;
        let mut rows = 0;
        for k in 1..=2048_usize {
            let omega = TAU * k as f64 / (4096.0 * dt);
            let f = omega / TAU;
            let (zh, zs, zf) = (head.z_at(omega).unwrap(), snow.z_at(omega).unwrap(), far.z_at(omega).unwrap());
            if f >= 1.25 * band {
                assert_eq!((zs.re.to_bits(), zs.im.to_bits()), (zh.re.to_bits(), zh.im.to_bits()), "{f} Hz");
                continue;
            }
            rows += 1;
            let w = if f <= band { 1.0 } else { 0.5 * (1.0 + (core::f64::consts::PI * (f - band) / (0.25 * band)).cos()) };
            let dz = s.change(f).unwrap().1.scale(w);
            let scale = zh.abs();
            assert!((zs.re - zh.re - dz.re).abs() < 1e-12 * scale && (zs.im - zh.im + dz.im).abs() < 1e-12 * scale,
                "{f} Hz: {zs:?} - {zh:?} against conj {dz:?}");
            assert!((zf - zh).abs() < 1e-6 * scale, "{f} Hz far torso: {zf:?} vs {zh:?}");
        }
        assert_eq!(rows, 426);
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
        let err = Series::new(0.1, 0.01, POLE, [0.0, 0.0, -0.05], AIR, 8000.0, 30).unwrap_err();
        assert!(err.contains("outside the head"), "{err}");
    }

    /// `(f, Re S, Im S)` [Hz, Pa s/m^3]: scipy 1.14 spherical_jn/yn, 400 orders, head 8.75 cm, cap of radius 13.548 mm,
    /// microphone (0.212132, 0, 0.212132) m, rho 1.1336, c 356.53.
    const GOLDEN: [(f64, f64, f64); 3] = [
        (300.0, 6.273927074538665e+02, 9.111865394641030e+01),
        (1000.0, -2.569048431744369e+03, 9.166961093269957e+02),
        (3000.0, 9.514481124514730e+01, -9.554082937523743e+03),
    ];

    const TAU: f64 = core::f64::consts::TAU;

    fn air_state() -> GasState {
        GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
    }

    fn tube(radius: f64) -> Duct {
        Duct { segments: vec![Segment::Cylinder { radius, length: 0.17 }] }
    }

    /// The cap's self-impedance against an independent Python evaluation of the same series, summed to 51672 orders,
    /// whose terms match mpmath's exact spherical Hankel functions to 6e-16 and which reproduces Aarts and Janssen's
    /// printed cap-on-sphere power curve (NAG-Journaal summary, Fig. 6, on its axial profile) within 0.02 dB median.
    #[test]
    fn the_load_series_matches_an_independent_evaluation() {
        let cap = CapLoad::new(0.0875, 0.013_548, AIR.density, AIR.sound_speed, TAU * 3000.0).unwrap();
        for (f, re, im) in GOLDEN_LOAD {
            let z = cap.impedance(TAU * f).unwrap();
            let want = Cx::new(re, im);
            assert!(z.sub(want).abs() < 1e-6 * want.abs(), "{f} Hz: {z:?} vs {want:?}");
        }
    }

    /// The whole sphere pulsating (a cap of area `4 pi a^2`, `b = 2 a`) is `rho c j k a / (1 + j k a)` over that area.
    #[test]
    fn the_whole_sphere_loads_as_the_pulsating_sphere() {
        let a = 0.0875;
        let cap = CapLoad::new(a, 2.0 * a, AIR.density, AIR.sound_speed, TAU * 10_000.0).unwrap();
        for f in [50.0, 1000.0, 10_000.0] {
            let x = TAU * f / AIR.sound_speed * a;
            let area = 4.0 * core::f64::consts::PI * a * a;
            let want = Cx::new(0.0, x).div(Cx::new(1.0, x)).scale(AIR.density * AIR.sound_speed / area);
            let z = cap.impedance(TAU * f).unwrap();
            assert!(z.sub(want).abs() < 1e-12 * want.abs(), "{f} Hz: {z:?} vs {want:?}");
        }
    }

    /// As `k a -> 0` the cap radiates into the whole space: its resistance tends to `rho c k^2 / (4 pi)`, half the
    /// baffled piston's (Python: 0.999985 at 5 Hz for a 1 cm mouth).
    #[test]
    fn the_low_frequency_load_is_half_the_baffled_pistons() {
        let omega = TAU * 5.0;
        let cap = CapLoad::new(0.0875, 0.01, AIR.density, AIR.sound_speed, omega).unwrap();
        let k = omega / AIR.sound_speed;
        let half = AIR.density * AIR.sound_speed * k * k / (4.0 * core::f64::consts::PI);
        let on_sphere = cap.impedance(omega).unwrap().re / half;
        let (r_piston, _) = fs_phs::baffled_piston_impedance(AIR.density, AIR.sound_speed, 0.01, omega).unwrap();
        assert!((on_sphere - 1.0).abs() < 1e-4, "{on_sphere}");
        assert!((r_piston / half - 2.0).abs() < 1e-6, "{}", r_piston / half);
    }

    /// A large head at a fixed mouth loads it as the baffle: `|conj(Z_cap) / Z_piston - 1|`, the largest over 200 and
    /// 1000 Hz, falls like `1 / a` (Python: 9.81e-2, 5.27e-2 and 2.72e-2 at 0.4, 0.8 and 1.6 m, 5 cm mouth; the same
    /// with 3000 Hz, whose deviation is smaller, but there the 1.6 m head needs 83k orders, past `LOAD_MAX_ORDER`).
    #[test]
    fn a_large_head_loads_the_lips_as_the_baffle() {
        let off = |a: f64| -> f64 {
            let cap = CapLoad::new(a, 0.05, AIR.density, AIR.sound_speed, TAU * 1000.0).unwrap();
            [200.0, 1000.0]
                .into_iter()
                .map(|f| {
                    let c = cap.impedance(TAU * f).unwrap();
                    let (r, x) = fs_phs::baffled_piston_impedance(AIR.density, AIR.sound_speed, 0.05, TAU * f).unwrap();
                    Cx::new(c.re, -c.im).div(Cx::new(r, x)).sub(Cx::ONE).abs()
                })
                .fold(0.0, f64::max)
        };
        let devs: Vec<f64> = [0.4, 0.8, 1.6].into_iter().map(off).collect();
        assert!(devs.windows(2).all(|w| (0.45..0.6).contains(&(w[1] / w[0]))), "{devs:?}");
        assert!(devs[2] < 0.03, "{devs:?}");
    }

    /// At the live table's top (96 kHz, k a 148) on the cards' smallest mouth (8.63 mm), and at k a 26, where the
    /// series once refused under the `1000 a / b` rule alone, the sums to `N` and `2 N` agree; the reactance's tail
    /// there grows with `k a` (module doc).
    #[test]
    fn the_load_series_converges_at_the_table_top() {
        let cap = CapLoad::new(0.0875, 0.008_63, AIR.density, AIR.sound_speed, TAU * 96_000.0).unwrap();
        for f in [16_242.0, 48_000.0, 96_000.0] {
            let z = cap.impedance(TAU * f).unwrap();
            assert!(z.re > 0.0 && (z.abs() / cap.scale - 1.0).abs() < 0.2, "{f} Hz: {z:?}");
        }
    }

    /// A mouth so small against the head that the series would need more than `LOAD_MAX_ORDER` terms is refused.
    #[test]
    fn a_load_series_too_long_is_refused() {
        // 2 * 1000 a / b = 175,000 orders, past 2^17 (a 2 mm mouth's 87,500 was past the old 2^16).
        let err = CapLoad::new(0.0875, 0.001, AIR.density, AIR.sound_speed, TAU * 1000.0).unwrap_err();
        assert!(err.contains("more than"), "{err}");
    }

    /// The C3-a card's left nostril on the baffle (MB 18119) and Dang's subject 2 nostril radii (0.81 and 0.62 cm^2).
    const NOSTRIL: [f64; 2] = [0.007_07, 0.027_63];
    const NOSTRIL_RADII: [f64; 2] = [5.078e-3, 4.442e-3];

    /// Every row of a table on the bins of `taps` at `rate`, as bits (the rows' own expression, so each query hits its row).
    fn bits(table: &TabulatedLoad, rate: u32, taps: usize) -> Vec<(u64, u64)> {
        let dt = 1.0 / f64::from(rate);
        (1..=taps / 2)
            .map(|k| table.z_at(TAU * k as f64 / (taps as f64 * dt)).unwrap())
            .map(|z| (z.re.to_bits(), z.im.to_bits()))
            .collect()
    }

    /// N3b test (2) in the head's own terms: an opening at the mouth's point and radius is the mouth, bit for bit, for
    /// the FIR (head alone and with the torso, NavyTower 18811 point b) and for the load table (alone and with the torso).
    #[test]
    fn an_opening_at_the_pole_is_the_mouth_bit_for_bit() {
        let (a, b, band) = (0.087, 0.013_548, 4000.0);
        let rx = mic(16, 64, band);
        let pole = CircularOpening { center_m: [0.0, 0.0], radius_m: b };
        let taps = |h: &HeadSphere| h.taps().iter().map(|t| t.to_bits()).collect::<Vec<_>>();
        let head = HeadSphere::new(a, b, rx, AIR, 48_000).unwrap();
        let opening = HeadSphere::for_opening(a, pole, rx, AIR, 48_000, None).unwrap();
        assert_eq!(taps(&head), taps(&opening));
        assert_eq!(head.latency_samples(), opening.latency_samples());
        let head = HeadSphere::with_torso(a, b, rx, AIR, 48_000, KEMAR).unwrap();
        let opening = HeadSphere::for_opening(a, pole, rx, AIR, 48_000, Some(KEMAR)).unwrap();
        assert_eq!(taps(&head), taps(&opening));
        let gas = air_state();
        let mouth = head_mouth_load(a, &tube(b), &gas, &gas, 48_000, 4096).unwrap();
        let cap = head_opening_load(a, pole, None, 0.17, &gas, &gas, 48_000, 4096).unwrap();
        assert_eq!(bits(&mouth, 48_000, 4096), bits(&cap, 48_000, 4096));
        assert_eq!((mouth.mouth(), mouth.source_id()), (cap.mouth(), cap.source_id()));
        let mouth = head_torso_mouth_load(a, KEMAR, band, &tube(b), &gas, &gas, 48_000, 4096).unwrap();
        let cap = head_opening_load(a, pole, Some((KEMAR, band)), 0.17, &gas, &gas, 48_000, 4096).unwrap();
        assert_eq!(bits(&mouth, 48_000, 4096), bits(&cap, 48_000, 4096));
        assert_eq!((mouth.mouth(), mouth.source_id()), (cap.mouth(), cap.source_id()));
    }

    /// An opening off the pole (module doc, "Openings off the pole"): its cap sits on the ray through its place, 4.33 mm
    /// behind the mouth's plane for the C3-a nostril; the cap's own load does not depend on where the cap sits, so
    /// without a torso the nostril's table is the pole cap's of the same radius bit for bit, while the torso's change does
    /// depend on it and the table stays passive; and a head far larger than the openings reads each opening's own
    /// baffle (G near 1).
    #[test]
    fn an_opening_off_the_pole_is_a_cap_on_its_ray() {
        let a = 0.087;
        let e = cap_axis(a, NOSTRIL);
        let behind = a * (1.0 - e[2]);
        assert!((behind - 4.33e-3).abs() < 0.01e-3, "{behind}");
        assert!((f64::asin(e[1]).to_degrees() - 17.56).abs() < 0.01, "{e:?}");
        assert!((f64::atan2(e[0], e[2]).to_degrees() - 4.65).abs() < 0.01, "{e:?}");
        let gas = air_state();
        let nostril = CircularOpening { center_m: NOSTRIL, radius_m: NOSTRIL_RADII[0] };
        let pole = CircularOpening { center_m: [0.0, 0.0], ..nostril };
        let alone = head_opening_load(a, nostril, None, 0.2, &gas, &gas, 48_000, 4096).unwrap();
        let alone_pole = head_opening_load(a, pole, None, 0.2, &gas, &gas, 48_000, 4096).unwrap();
        assert_eq!(bits(&alone, 48_000, 4096), bits(&alone_pole, 48_000, 4096));
        let snow = head_opening_load(a, nostril, Some((KEMAR, 4000.0)), 0.2, &gas, &gas, 48_000, 4096).unwrap();
        let snow_pole = head_opening_load(a, pole, Some((KEMAR, 4000.0)), 0.2, &gas, &gas, 48_000, 4096).unwrap();
        let (off, on) = (bits(&snow, 48_000, 4096), bits(&snow_pole, 48_000, 4096));
        assert_ne!(off, on, "the torso's change depends on where the cap sits");
        snow.validate_mouth(0.2, NOSTRIL_RADII[0], &gas).unwrap();
        // As `a_large_head_reads_the_baffle` for the mouth: G -> 1 like 1 / a as the head grows, each doubling scaling
        // |G - 1| by 0.3 to 0.7 (declared before the first run). scipy on the same quadrature, with the cap on its ray
        // (NavyTower's golden method; the mouth's values reproduced as its control): |G - 1| at 8 kHz is 1.21e-1,
        // 5.85e-2, 2.86e-2 and 1.408e-2 for 0.25, 0.5, 1 and 2 m, its cap tilted and recessed by amounts that fall as
        // 1 / a. The bound 1.6e-2 is that value plus 14%, near the mouth test's margin (6e-3 over 5.2e-3, 15%). The
        // first run's bound, 1e-2, set without an estimate, was wrong: that run read 0.01408 (MB 19228).
        let off = |a: f64| -> f64 {
            let head = HeadSphere::for_opening(a, nostril, mic(16, 64, 8000.0), AIR, 48_000, None).unwrap();
            g(&head, 8000.0).sub(Cx::ONE).abs()
        };
        let devs: Vec<f64> = [0.25, 0.5, 1.0, 2.0].into_iter().map(off).collect();
        assert!(devs.windows(2).all(|w| (0.3..0.7).contains(&(w[1] / w[0]))), "{devs:?}");
        assert!(devs[3] < 1.6e-2, "2 m head at 8 kHz: |G - 1| = {}", devs[3]);
    }

    /// The nostrils' caps on the 262144-tap bins at 192 kHz need more orders than the old 2^16 and fit 2^17 (MB 18809).
    #[test]
    fn a_nostrils_cap_fits_the_load_order_cap() {
        let top = core::f64::consts::PI * 192_000.0;
        for b in NOSTRIL_RADII {
            let cap = CapLoad::new(0.087, b, 1.1336, 353.0, top).unwrap();
            assert!(cap.w.len() - 1 > 1 << 16 && cap.w.len() - 1 <= LOAD_MAX_ORDER, "{b}: {}", cap.w.len() - 1);
        }
    }

    /// The timing probe MB 18809 asked for: one nostril's table with the torso on the 262144-tap bins at 192 kHz, as the
    /// nasal card's build makes it; wall time and peak RSS (Linux `VmHWM`) as one JSON line.
    #[test]
    #[ignore = "timing probe for the build host (N3b): run by name with --ignored --nocapture"]
    fn a_nostril_table_build_time_at_the_live_clock() {
        let gas = air_state();
        let nostril = CircularOpening { center_m: NOSTRIL, radius_m: NOSTRIL_RADII[0] };
        let started = std::time::Instant::now();
        let table =
            head_opening_load(0.087, nostril, Some((KEMAR, 8000.0)), 0.2, &gas, &gas, 192_000, 262_144).unwrap();
        let wall = started.elapsed().as_secs_f64();
        let hwm = std::fs::read_to_string("/proc/self/status")
            .ok()
            .and_then(|s| {
                let line = s.lines().find(|l| l.starts_with("VmHWM:"))?;
                Some(line.split_whitespace().nth(1).unwrap_or("").to_string())
            })
            .unwrap_or_else(|| "n/a".into());
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"head-nostril-table-build-time\",\"rows\":{},\"wall_s\":{wall:.2},\
             \"peak_rss_kib\":\"{hwm}\",\"workers\":{}}}",
            bits(&table, 192_000, 262_144).len(),
            std::thread::available_parallelism().map_or(1, |p| p.get()).min(8)
        );
    }

    /// The table: rows on the finest rung's bins, bound to the duct's mouth, passive; at 200 Hz less mass and less
    /// resistance than the flanged load, an end correction of 0.69-0.76 radii against the flanged 0.82 (the stamped
    /// prediction, mail 7212, read 0.71-0.74 for the card mouths); and a very large head gives the flanged load back.
    #[test]
    fn the_lips_load_is_the_flanged_load_plus_the_spheres_change() {
        let gas = air_state();
        let b = 0.013_548;
        let table = head_mouth_load(0.0875, &tube(b), &gas, &gas, 48_000, 4096).unwrap();
        let dt = 1.0 / 48_000.0;
        assert_eq!(table.support(), (TAU * 1.0 / (4096.0 * dt), TAU * 2048.0 / (4096.0 * dt)));
        table.validate_mouth(0.17, b, &gas).unwrap();
        let omega = TAU * 17.0 / (4096.0 * dt);
        let z = table.z_at(omega).unwrap();
        let flanged = fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, b, omega).unwrap();
        assert!(z.im > flanged.im, "less mass on the sphere: {z:?} vs {flanged:?}");
        assert!(z.re > 0.0 && z.re < flanged.re, "less resistance on the sphere: {z:?} vs {flanged:?}");
        let end = -z.im * core::f64::consts::PI * b * b / (gas.density * omega) / b;
        assert!((0.69..0.76).contains(&end), "end correction {end} radii");
        // At 2 kHz, so the 1.6 m head's table (top k a 29) stays within LOAD_MAX_ORDER; bin 120 is 937.5 Hz.
        let back = |a: f64| -> f64 {
            let far = head_mouth_load(a, &tube(0.05), &gas, &gas, 2_000, 256).unwrap();
            let omega = TAU * 120.0 / (256.0 / 2_000.0);
            let flanged = fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, 0.05, omega).unwrap();
            let z = far.z_at(omega).unwrap();
            Cx::new(z.re, z.im).div(Cx::new(flanged.re, flanged.im)).sub(Cx::ONE).abs()
        };
        let (near, far) = (back(0.4), back(1.6));
        assert!(far < 0.03 && far < 0.4 * near, "1.6 m head {far}, 0.4 m head {near}");
    }

    /// The table with no head (MB 24424). With one gas every row is today's flanged load
    /// (`compact_radiation_impedance`) bit for bit. With a tract gas unlike the exterior every row is fs-duct's
    /// two-gas flanged load (radiation in the exterior, the eliminated duct modes in the tract) bit for bit, called
    /// directly, and no row is the one-gas load. Both tables lie on the head table's support and carry the duct's mouth
    /// in the exterior with DC 0. The rate, tap-count, empty-duct and static-pressure refusals hold.
    #[test]
    fn the_flanged_table_is_the_two_gas_flanged_load_on_every_row() {
        let gas = air_state();
        let warm = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 310.15, 101_325.0).unwrap();
        assert_ne!(warm.sound_speed.to_bits(), gas.sound_speed.to_bits(), "the test needs two gases");
        let b = 0.013_548;
        let (rate, taps) = (48_000_u32, 4096_usize);
        let one = flanged_mouth_load(&tube(b), &gas, &gas, rate, taps).unwrap();
        let two = flanged_mouth_load(&tube(b), &warm, &gas, rate, taps).unwrap();
        let head = head_mouth_load(0.0875, &tube(b), &gas, &gas, rate, taps).unwrap();
        for table in [&one, &two] {
            assert_eq!(table.support(), head.support());
            let mouth = table.validate_mouth(0.17, b, &gas).unwrap();
            assert_eq!((mouth.dc_impedance.to_bits(), mouth.gas == gas), (0.0_f64.to_bits(), true));
        }
        assert!(two.validate_mouth(0.17, b, &warm).is_err(), "the table is bound to the exterior, not the tract");
        let dt = 1.0 / f64::from(rate);
        let bits = |z: fs_math::c64::C64| (z.re.to_bits(), z.im.to_bits());
        for k in 1..=taps / 2 {
            let omega = TAU * k as f64 / (taps as f64 * dt);
            let today = fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, b, omega).unwrap();
            assert_eq!(bits(one.z_at(omega).unwrap()), bits(today), "one gas, bin {k}");
            let law = fs_duct::flanged_pipe_impedance_with_exterior(&warm, &gas, b, omega, fs_duct::FLANGED_MODES)
                .unwrap();
            assert_eq!(bits(two.z_at(omega).unwrap()), bits(law), "two gases, bin {k}");
            assert_ne!(bits(law), bits(today), "two gases, bin {k}: the tract's gas moved nothing");
        }
        for (rate, taps) in [(0, 256), (48_000, 300), (48_000, 1)] {
            let err = flanged_mouth_load(&tube(b), &gas, &gas, rate, taps).unwrap_err();
            assert!(err.contains("power-of-two"), "{err}");
        }
        let err = flanged_mouth_load(&Duct { segments: vec![] }, &gas, &gas, 48_000, 256).unwrap_err();
        assert!(err.contains("segments"), "{err}");
        let thin = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
        let err = flanged_mouth_load(&tube(b), &thin, &gas, 48_000, 256).unwrap_err();
        assert!(err.contains("refused"), "{err}");
    }

    /// The head's table with a tract gas (MB 24424, 25423): its flanged term is fs-duct's two-gas load, the radiation in
    /// the exterior and the eliminated duct modes in the tract, and nothing else in the table moves. Without a torso,
    /// for the lips' cap and a nostril's cap off the pole, every row is `cap_load`'s row expression on the same cap and
    /// piston terms: with one gas on today's flanged load (`compact_radiation_impedance`), bit for bit, and with a tract
    /// gas unlike the exterior on `flanged_pipe_impedance_with_exterior` called directly, bit for bit; no row is the
    /// one-gas row. With the torso each row moves by the direct flanged term's change to within 1e-12 of the row (the
    /// row is a floating-point sum, so the change is not exact there) and every row moves. Both tables stay bound to the
    /// exterior, and unequal static pressures are refused by the flanged load.
    #[test]
    fn the_head_tables_flanged_term_carries_the_duct_modes_in_the_tract_gas() {
        let gas = air_state();
        let warm = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 310.15, 101_325.0).unwrap();
        assert_ne!(warm.sound_speed.to_bits(), gas.sound_speed.to_bits(), "the test needs two gases");
        let (a, lips, rate, taps) = (0.0875, 0.013_548, 48_000_u32, 4096_usize);
        let dt = 1.0 / f64::from(rate);
        let omegas: Vec<f64> = (1..=taps / 2).map(|k| TAU * k as f64 / (taps as f64 * dt)).collect();
        let nostril = CircularOpening { center_m: NOSTRIL, radius_m: NOSTRIL_RADII[0] };
        let cases = [
            ("lips", lips, head_mouth_load(a, &tube(lips), &gas, &gas, rate, taps).unwrap(),
                head_mouth_load(a, &tube(lips), &warm, &gas, rate, taps).unwrap()),
            ("nostril", nostril.radius_m, head_opening_load(a, nostril, None, 0.2, &gas, &gas, rate, taps).unwrap(),
                head_opening_load(a, nostril, None, 0.2, &warm, &gas, rate, taps).unwrap()),
        ];
        for (name, b, one, two) in &cases {
            let cap = CapLoad::new(a, *b, gas.density, gas.sound_speed, omegas[omegas.len() - 1]).unwrap();
            let (one, two) = (bits(one, rate, taps), bits(two, rate, taps));
            for (k, &omega) in omegas.iter().enumerate() {
                let c = cap.impedance(omega).unwrap();
                let (r_piston, x_piston) =
                    fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, *b, omega).unwrap();
                // `cap_load`'s row on a given flanged term.
                let row = |f: fs_math::c64::C64| {
                    let z = Cx::new(f.re + (c.re - r_piston), f.im + (-c.im - x_piston));
                    (z.re.to_bits(), z.im.to_bits())
                };
                let today = fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, *b, omega).unwrap();
                let law = fs_duct::flanged_pipe_impedance_with_exterior(&warm, &gas, *b, omega, fs_duct::FLANGED_MODES)
                    .unwrap();
                assert_eq!(one[k], row(today), "{name}, one gas, bin {}", k + 1);
                assert_eq!(two[k], row(law), "{name}, two gases, bin {}", k + 1);
                assert_ne!(two[k], one[k], "{name}, two gases, bin {}: the tract's gas moved nothing", k + 1);
            }
        }
        let band = 4000.0;
        let one = head_torso_mouth_load(a, KEMAR, band, &tube(lips), &gas, &gas, rate, taps).unwrap();
        let two = head_torso_mouth_load(a, KEMAR, band, &tube(lips), &warm, &gas, rate, taps).unwrap();
        for (k, &omega) in omegas.iter().enumerate() {
            let (z1, z2) = (one.z_at(omega).unwrap(), two.z_at(omega).unwrap());
            let today = fs_duct::compact_radiation_impedance(Termination::FlangedOpen, &gas, lips, omega).unwrap();
            let law = fs_duct::flanged_pipe_impedance_with_exterior(&warm, &gas, lips, omega, fs_duct::FLANGED_MODES)
                .unwrap();
            let moved = Cx::new(z2.re - z1.re, z2.im - z1.im);
            let direct = Cx::new(law.re - today.re, law.im - today.im);
            let miss = moved.sub(direct).abs();
            let bound = 1e-12 * Cx::new(z2.re, z2.im).abs();
            assert!(miss <= bound, "torso, bin {}: moved {moved:?}, direct {direct:?}", k + 1);
            assert!(moved.abs() > 0.0, "torso, bin {}: the tract's gas moved nothing", k + 1);
        }
        for table in [&cases[0].3, &two] {
            table.validate_mouth(0.17, lips, &gas).unwrap();
            let tract_bound = table.validate_mouth(0.17, lips, &warm);
            assert!(tract_bound.is_err(), "the table is bound to the exterior, not the tract");
        }
        let thin = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
        let err = head_mouth_load(a, &tube(lips), &thin, &gas, 48_000, 256).unwrap_err();
        assert!(err.contains("flanged load refused"), "{err}");
    }

    #[test]
    fn non_physical_lip_loads_are_refused() {
        let gas = air_state();
        for (a, r, rate, taps, why) in [
            (0.0, 0.01, 48_000, 256, "finite and positive"),
            (f64::NAN, 0.01, 48_000, 256, "finite and positive"),
            (0.01, 0.01, 48_000, 256, "smaller than the head"),
            (0.0875, 0.01, 0, 256, "power-of-two"),
            (0.0875, 0.01, 48_000, 300, "power-of-two"),
            // 2N = 175,000 orders at 48 kHz, over LOAD_MAX_ORDER (2^17 since N3b; 2 mm needs 87,500 and is admitted).
            (0.0875, 0.001, 48_000, 256, "more than"),
        ] {
            let err = head_mouth_load(a, &tube(r), &gas, &gas, rate, taps).unwrap_err();
            assert!(err.contains(why), "{err}");
        }
        let err = head_mouth_load(0.0875, &Duct { segments: vec![] }, &gas, &gas, 48_000, 256).unwrap_err();
        assert!(err.contains("segments"), "{err}");
    }

    /// `(f, Re Z_cap, Im Z_cap)` [Hz, Pa s/m^3, `e^{+j w t}`]: `voice-slow/hl/golden.py` (numpy, 51672 orders; the change
    /// from half as many is 2e-8), head 8.75 cm, cap of radius 13.548 mm, rho 1.1336, c 356.53.
    const GOLDEN_LOAD: [(f64, f64, f64); 3] = [
        (300.0, 8.838839269047968e+02, 3.715138455023602e+04),
        (1000.0, 1.305850074888683e+04, 1.256795001559119e+05),
        (3000.0, 1.384121113452622e+05, 3.485720634727473e+05),
    ];
}
