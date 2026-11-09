//! The snowman: a rigid torso sphere under the head sphere (opt-in, [`TorsoSphere`]), for the microphone and the lips' load.
//!
//! Print: Algazi, Duda and Thompson 2002, "The use of head-and-torso models for improved spatial sound synthesis", AES
//! 113th Convention, Los Angeles, Fig. 2a and p. 3 ("a = 8.7 cm, b = 16.9 cm, and h = 5.3 cm" for the KEMAR mannequin).
//! The head (radius `a`, centre `H`) sits directly above the torso (radius `b`, centre `T`); the two are coaxial and
//! `h` is the neck gap between their surfaces, so `|HT| = D = a + h + b` (the print's own delay estimate
//! `2 (sqrt(D^2 + a^2) - b) / c = 0.89 ms` holds only with that `D`). The print computed the snowman by multipole
//! re-expansion from 0 to 5 kHz; KEMAR is a mannequin, not a singer, and the print says the sphere is not a good
//! approximation of the human torso, particularly at low elevations.
//!
//! Placement, our choices (no print makes them): the mouth cap is centred on the head's horizontal equator, facing front,
//! and the card frame's `+y` is up. In the card frame (origin at the mouth centre, `+z` out of the mouth) `H = (0, 0, -a)`
//! and `T = (0, -D, -a)`. The snowman frame puts its polar axis up (card `+y`) and azimuth 0 at the front (card `+z`):
//! about `H`, `x_s = z + a`, `y_s = x`, `z_s = y`. The cap's centre is at polar angle `pi / 2`, azimuth 0.
//!
//! An opening off the mouth's pole (the nostrils; the head module doc, "Openings off the pole") puts its cap on the ray
//! from `H` through the opening's centre on the baffle, a modelling choice. With that ray's unit vector `e` (card frame),
//! the cap's centre is at `cos theta_c = e_y` and azimuth `phi_c` with `(cos phi_c, sin phi_c) = (e_z, e_x) / hypot(e_z,
//! e_x)`. The torso is coaxial with the polar axis, so the field is even in `phi - phi_c`: the addition theorem below
//! carries `P_n^m(cos theta_c)` in place of `P_n^m(0)` (in `V_n^m` and in `dZ`), and the microphone's azimuth is
//! measured from `phi_c`. The lips' cap passes `e = (0, 0, 1)`: `cos theta_c = 0.0` and `(1.0, 0.0)` exactly, so every
//! product is the one it was.
//!
//! Law (time factor `e^{+j w t}`, outgoing `h_n = h_n^(2)`, as the head's own series). The exterior field is the head's
//! outgoing multipoles about `H` plus the torso's about `T`:
//!
//! ```text
//! p = sum_m cos(m phi) [ sum_n alpha_n^m h_n(k r_H) / h_n'(k a) P_n^m(cos th_H) + sum_n beta_n^m h_n(k r_T) / h_n'(k b) P_n^m(cos th_T) ]
//! ```
//!
//! `P_n^m` is normalised so that `int_{-1}^{1} (P_n^m)^2 dx = 1` (no Condon-Shortley phase; it cancels in every product
//! used here). Both centres lie on the polar axis, so they share `phi`; the cap is centred at `phi = 0`, so the field is
//! even in `phi` and every `m` stands alone. Boundary conditions: `dp / dr_H = -j w rho v` on the head (`v = 1 / A` on the
//! cap per unit volume flow, 0 elsewhere) and `dp / dr_T = 0` on the torso, each projected on its own sphere's `P_n^m` by
//! Gauss-Legendre quadrature in `cos th` (the `m`-th Fourier part is exact by the symmetry). The cap's velocity is the
//! head's own `sum_n u_n P_n(cos gamma)`, and the addition theorem `P_n(cos gamma) = sum_m eps_m 2 / (2n+1) P_n^m(cos th)
//! P_n^m(0) cos(m phi)` (`eps_0 = 1`, `eps_m = 2`) gives its `m`-th part `V_n^m`. The isolated head's own coefficients are
//! `alpha_iso = -j rho c V` (exactly the head's series), so the unknowns are the change `delta = alpha - alpha_iso` and
//! `beta`:
//!
//! ```text
//! k delta + K_HT beta = 0,   k beta + K_TH (alpha_iso + delta) = 0,   so   (k^2 I - K_HT K_TH) delta = K_HT K_TH alpha_iso
//! K_HT[n][n'] = int_{-1}^{1} d/dr_H [ h_n'(k r_T) / h_n''(k b) P_n'^m(cos th_T) ] P_n^m(x) dx   on the head (r_H = a)
//! K_TH[n'][n] = int_{-1}^{1} d/dr_T [ h_n(k r_H) / h_n'(k a) P_n^m(cos th_H) ] P_n'^m(x) dx     on the torso (r_T = b)
//! ```
//!
//! Outputs, both per unit volume flow: the change of the microphone's pressure (added to the head's own `S`), and the
//! change of the cap's mean pressure, which is the change of the cap's self-impedance (added to the head's own `Z_cap`):
//! `dZ = (a^2 / A) sum_m sum_n dc_n^m lam_n P_n^m(0)`, with `dc = delta h_n(k a) / h_n'(k a) + L_HT beta` the change of the
//! head's surface pressure (`L_HT` as `K_HT` with the value in place of the derivative) and `lam_n = 2 pi int_{cos al}^{1}
//! P_n` (Funk-Hecke over the cap).
//!
//! Truncation, at each frequency: head orders to `ceil(k a) + 20`, torso orders to `ceil(k b) + 35`, and on each sphere
//! 15 quadrature nodes more than its top order at the band's top. The torso's field on the head converges like
//! `(a / D)^n`, the head's on the torso like `(b / (D - a))^n` past `k b`. A numpy prototype of the same law (scratch
//! `torso/proto.py`) agrees with twice these margins within 1e-11 of the microphone's pressure and 5e-10 of `dZ` from
//! 50 Hz to 10 kHz; [`tests`] checks it here. Radial functions are carried as ratios (as the head's series), so
//! nothing overflows at small `k`; plain IEEE operations only, so every build rounds the same way, and each frequency is
//! solved alone, so the result does not depend on how many threads share the bins.

use super::Cx;
use fs_math::det;

/// Head orders kept past `k a`.
const HEAD_MARGIN: usize = 20;
/// Torso orders kept past `k b`.
const TORSO_MARGIN: usize = 35;
/// Quadrature nodes on each sphere past its top order.
const NODE_MARGIN: usize = 15;
/// Largest order admitted on either sphere; a geometry or band that needs more is refused, never truncated.
pub const TORSO_MAX_ORDER: usize = 400;
/// Most worker threads for the bins (each bin is solved alone, so the count changes no result).
const MAX_WORKERS: usize = 8;

/// A torso sphere under the head (module doc): its radius and the neck gap between the two surfaces [m].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsoSphere {
    /// Torso radius `b` [m].
    pub radius_m: f64,
    /// Neck gap `h` between the head's and the torso's surfaces [m]; the centres are `a + h + b` apart.
    pub neck_gap_m: f64,
}

impl TorsoSphere {
    /// # Errors
    /// A radius or gap that is not finite and positive (the spheres must not touch: the expansions need `D > a + b`).
    pub fn validate(&self) -> Result<(), String> {
        if !(self.radius_m.is_finite() && self.radius_m > 0.0) {
            return Err("the torso radius must be finite and positive".into());
        }
        if !(self.neck_gap_m.is_finite() && self.neck_gap_m > 0.0) {
            return Err("the neck gap must be finite and positive (the head and torso must not touch)".into());
        }
        Ok(())
    }
}

/// The microphone in the snowman frame (module doc), with its Legendre values per `m`.
#[derive(Clone, Debug)]
struct Mic {
    r_h: f64,
    r_t: f64,
    /// `cos(m phi)`, `m = 0..=nh`.
    cos_m: Vec<f64>,
    /// Per `m`: `P_n^m(cos th_H)`, `n = m..=nh`, and `P_n^m(cos th_T)`, `n = m..=nt`.
    p_h: Vec<Vec<f64>>,
    p_t: Vec<Vec<f64>>,
}

/// The angular tables of one `m`, at the band's top orders (rows `n - m`, columns nodes).
#[derive(Clone, Debug)]
struct Tables {
    /// Head orders at the head's nodes: `P_n^m(x_H)`.
    head_at_head: Vec<Vec<f64>>,
    /// Torso orders at the head's nodes: `P_n^m(cos th_T)` and `d/dth_T`.
    torso_at_head: Vec<Vec<f64>>,
    torso_at_head_d: Vec<Vec<f64>>,
    /// Torso orders at the torso's nodes: `P_n^m(x_T)`.
    torso_at_torso: Vec<Vec<f64>>,
    /// Head orders at the torso's nodes: `P_n^m(cos th_H)` and `d/dth_H`.
    head_at_torso: Vec<Vec<f64>>,
    head_at_torso_d: Vec<Vec<f64>>,
    /// `P_n^m(cos theta_c)`, `n = m..=nh` (the cap's centre; `P_n^m(0)` for the lips).
    at_cap: Vec<f64>,
}

/// Nodes of one sphere: weights, the other centre's distance and the two direction cosines of `d/dr` (module doc).
#[derive(Clone, Debug)]
struct Nodes {
    w: Vec<f64>,
    /// Distance from the OTHER sphere's centre [m].
    r_other: Vec<f64>,
    /// `rhat_own . rhat_other` and `rhat_own . thhat_other`.
    rr: Vec<f64>,
    rt: Vec<f64>,
}

/// The snowman prepared for one geometry and one band (module doc).
#[derive(Clone, Debug)]
pub(super) struct Snowman {
    a: f64,
    b: f64,
    density: f64,
    sound_speed: f64,
    head_margin: usize,
    torso_margin: usize,
    /// Top orders (at the band's top).
    nh: usize,
    nt: usize,
    head: Nodes,
    torso: Nodes,
    tables: Vec<Tables>,
    /// The cap: `u_n` [1/m^2] and `lam_n`, `n = 0..=nh`, and its area [m^2].
    u: Vec<f64>,
    lam: Vec<f64>,
    area: f64,
    mic: Option<Mic>,
}

impl Snowman {
    /// The snowman for a cap of the mouth's area (radius `mouth_m`) on a head of radius `a`, the torso `torso` under it,
    /// the microphone `mic` in the card frame (or none, for the lips' load alone), up to `top_hz`.
    ///
    /// # Errors
    /// A torso refusal, a mouth not smaller than the head, a microphone inside either sphere, a band that is not finite
    /// and positive, or orders past [`TORSO_MAX_ORDER`].
    pub(super) fn new(a: f64, mouth_m: f64, torso: TorsoSphere, mic: Option<[f64; 3]>, density: f64, sound_speed: f64,
        top_hz: f64) -> Result<Self, String> {
        Self::at_cap(a, mouth_m, super::POLE, torso, mic, density, sound_speed, top_hz)
    }

    /// [`Self::new`] for a cap centred on the unit vector `cap_axis` from the head's centre (card frame; module doc):
    /// [`super::POLE`] is the lips' cap of [`Self::new`].
    ///
    /// # Errors
    /// Every refusal of [`Self::new`], or an axis that is not a unit vector facing front (`e_z > 0`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn at_cap(a: f64, cap_m: f64, cap_axis: [f64; 3], torso: TorsoSphere, mic: Option<[f64; 3]>,
        density: f64, sound_speed: f64, top_hz: f64) -> Result<Self, String> {
        Self::with_margins(a, cap_m, cap_axis, torso, mic, density, sound_speed, top_hz,
            [HEAD_MARGIN, TORSO_MARGIN, NODE_MARGIN])
    }

    #[allow(clippy::too_many_arguments)]
    fn with_margins(a: f64, mouth_m: f64, cap_axis: [f64; 3], torso: TorsoSphere, mic: Option<[f64; 3]>, density: f64,
        sound_speed: f64, top_hz: f64, margins: [usize; 3]) -> Result<Self, String> {
        torso.validate()?;
        if !(a.is_finite() && a > 0.0 && mouth_m > 0.0 && mouth_m < a) {
            return Err("the snowman needs a positive head and a mouth smaller than it".into());
        }
        let [ex, ey, ez] = cap_axis;
        let norm = ex * ex + ey * ey + ez * ez;
        if !(ez > 0.0 && (norm - 1.0).abs() <= 1.0e-12) {
            return Err("the snowman's cap axis must be a unit vector facing front".into());
        }
        // The cap's centre (module doc): cos theta_c and the azimuth's (cos, sin), exactly (0, 1, 0) for the lips.
        let cos_cap = ey;
        let h = det::hypot(ez, ex);
        let (cos_az, sin_az) = (ez / h, ex / h);
        if !(top_hz.is_finite() && top_hz > 0.0 && density > 0.0 && sound_speed > 0.0) {
            return Err("the snowman's band, density and sound speed must be positive".into());
        }
        let b = torso.radius_m;
        let d = a + torso.neck_gap_m + b;
        let [head_margin, torso_margin, node_margin] = margins;
        let k_top = core::f64::consts::TAU * top_hz / sound_speed;
        let nh = (k_top * a).ceil() as usize + head_margin;
        let nt = (k_top * b).ceil() as usize + torso_margin;
        if nh.max(nt) > TORSO_MAX_ORDER {
            return Err(format!("the snowman needs more than {TORSO_MAX_ORDER} orders at {top_hz} Hz"));
        }
        let (qh, qt) = (nh + node_margin, nt + node_margin);
        // Head nodes (relative to H, meridian plane): (a sin th, a cos th); relative to T: (a sin th, a cos th + D).
        let (xh, wh) = gauss_legendre(qh);
        let (xt, wt) = gauss_legendre(qt);
        let mut head = Nodes { w: wh, r_other: Vec::new(), rr: Vec::new(), rt: Vec::new() };
        let mut cos_t_at_head = Vec::with_capacity(qh);
        for &x in &xh {
            let s = det::sqrt(1.0 - x * x);
            let (rho, z) = (a * s, a * x + d);
            let r = det::hypot(rho, z);
            let (ct, st) = (z / r, rho / r);
            head.r_other.push(r);
            head.rr.push(s * st + x * ct);
            head.rt.push(s * ct - x * st);
            cos_t_at_head.push(ct);
        }
        // Torso nodes (relative to T): (b sin th, b cos th); relative to H: (b sin th, b cos th - D).
        let mut torso_nodes = Nodes { w: wt, r_other: Vec::new(), rr: Vec::new(), rt: Vec::new() };
        let mut cos_h_at_torso = Vec::with_capacity(qt);
        for &x in &xt {
            let s = det::sqrt(1.0 - x * x);
            let (rho, z) = (b * s, b * x - d);
            let r = det::hypot(rho, z);
            let (ch, sh) = (z / r, rho / r);
            torso_nodes.r_other.push(r);
            torso_nodes.rr.push(s * sh + x * ch);
            torso_nodes.rt.push(s * ch - x * sh);
            cos_h_at_torso.push(ch);
        }
        let columns = |xs: &[f64], m: usize, top: usize| -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
            let mut p = vec![Vec::with_capacity(xs.len()); top + 1 - m];
            let mut dp = vec![Vec::with_capacity(xs.len()); top + 1 - m];
            for &x in xs {
                let (v, dv) = pbar(x, m, top);
                for i in 0..v.len() {
                    p[i].push(v[i]);
                    dp[i].push(dv[i]);
                }
            }
            (p, dp)
        };
        let mut tables = Vec::with_capacity(nh + 1);
        // A torso smaller than the head in `k` units has no orders past `nt`: those `m` get no torso tables (and `change`
        // stops there, the torso's change being zero for them).
        for m in 0..=nh {
            let (head_at_head, _) = columns(&xh, m, nh);
            let (torso_at_head, torso_at_head_d, torso_at_torso) = if m <= nt {
                let (p, d) = columns(&cos_t_at_head, m, nt);
                (p, d, columns(&xt, m, nt).0)
            } else {
                (Vec::new(), Vec::new(), Vec::new())
            };
            let (head_at_torso, head_at_torso_d) = columns(&cos_h_at_torso, m, nh);
            let at_cap = pbar(cos_cap, m, nh).0;
            tables.push(Tables { head_at_head, torso_at_head, torso_at_head_d, torso_at_torso, head_at_torso,
                head_at_torso_d, at_cap });
        }
        // The cap (the head's own u_n), and lam_n = 2 pi int_{cos al}^1 P_n.
        let one_minus_cos_al = mouth_m * mouth_m / (2.0 * a * a);
        let p_cap = legendre(1.0 - one_minus_cos_al, nh + 2);
        let surface = 4.0 * core::f64::consts::PI * a * a;
        let mut u = vec![1.0 / surface];
        let mut lam = vec![core::f64::consts::TAU * one_minus_cos_al];
        for n in 1..=nh {
            let dn = p_cap[n - 1] - p_cap[n + 1];
            u.push(dn / (surface * one_minus_cos_al));
            lam.push(core::f64::consts::TAU * dn / (2 * n + 1) as f64);
        }
        let mic = match mic {
            None => None,
            Some(p) => {
                let (xs, ys, zs) = (p[2] + a, p[0], p[1]);
                let r_h = det::sqrt(xs * xs + ys * ys + zs * zs);
                let r_t = det::sqrt(xs * xs + ys * ys + (zs + d) * (zs + d));
                if !(r_h > a && r_t > b) {
                    return Err("the microphone must stand outside the head and the torso".into());
                }
                let rho = det::hypot(xs, ys);
                // cos(phi - phi_c): the microphone's azimuth from the cap's (xs / rho exactly for the lips).
                let c1 = if rho > 0.0 { (xs * cos_az + ys * sin_az) / rho } else { 1.0 };
                let mut cos_m = vec![1.0, c1];
                for m in 2..=nh {
                    cos_m.push(2.0 * c1 * cos_m[m - 1] - cos_m[m - 2]);
                }
                cos_m.truncate(nh + 1);
                let p_h = (0..=nh).map(|m| pbar(zs / r_h, m, nh).0).collect();
                let p_t = (0..=nh).map(|m| if m <= nt { pbar((zs + d) / r_t, m, nt).0 } else { Vec::new() }).collect();
                Some(Mic { r_h, r_t, cos_m, p_h, p_t })
            }
        };
        Ok(Self { a, b, density, sound_speed, head_margin, torso_margin, nh, nt, head, torso: torso_nodes, tables, u, lam,
            area: core::f64::consts::PI * mouth_m * mouth_m, mic })
    }

    /// The torso's change at `f_hz > 0`, per unit volume flow, `e^{+j w t}`: (the microphone's pressure [Pa s/m^3], zero
    /// without a microphone; the cap's self-impedance [Pa s/m^3]).
    ///
    /// # Errors
    /// A frequency above the prepared band, or a solve that is singular or not finite.
    pub(super) fn change(&self, f_hz: f64) -> Result<(Cx, Cx), String> {
        let k = core::f64::consts::TAU * f_hz / self.sound_speed;
        let nh = (k * self.a).ceil() as usize + self.head_margin;
        let nt = (k * self.b).ceil() as usize + self.torso_margin;
        if !(f_hz > 0.0) || nh > self.nh || nt > self.nt {
            return Err(format!("the snowman was not prepared for {f_hz} Hz"));
        }
        // Radial ratios: torso orders at the head's nodes, head orders at the torso's nodes, the head's surface, the mic.
        let (qh, qt) = (self.head.w.len(), self.torso.w.len());
        let mut ft = vec![vec![Cx::new(0.0, 0.0); qh]; nt + 1];
        let mut gt = ft.clone();
        let (mut f, mut g) = (vec![Cx::new(0.0, 0.0); nt.max(nh) + 1], vec![Cx::new(0.0, 0.0); nt.max(nh) + 1]);
        for q in 0..qh {
            radial(k * self.head.r_other[q], k * self.b, nt, &mut f, &mut g);
            for n in 0..=nt {
                ft[n][q] = f[n];
                gt[n][q] = g[n].scale(k);
            }
        }
        let mut fh = vec![vec![Cx::new(0.0, 0.0); qt]; nh + 1];
        let mut gh = fh.clone();
        for q in 0..qt {
            radial(k * self.torso.r_other[q], k * self.a, nh, &mut f, &mut g);
            for n in 0..=nh {
                fh[n][q] = f[n];
                gh[n][q] = g[n].scale(k);
            }
        }
        radial(k * self.a, k * self.a, nh, &mut f, &mut g);
        let on_head: Vec<Cx> = f[..=nh].to_vec(); // h_n(k a) / h_n'(k a)
        let (mic_h, mic_t) = match &self.mic {
            Some(mic) => {
                radial(k * mic.r_h, k * self.a, nh, &mut f, &mut g);
                let h = f[..=nh].to_vec();
                radial(k * mic.r_t, k * self.b, nt, &mut f, &mut g);
                (h, f[..=nt].to_vec())
            }
            None => (Vec::new(), Vec::new()),
        };
        let rho_c = self.density * self.sound_speed;
        let mut dp = Cx::new(0.0, 0.0);
        let mut dz = Cx::new(0.0, 0.0);
        for m in 0..=nh {
            if nt < m {
                break;
            }
            let tab = &self.tables[m];
            let (n_h, n_t) = (nh - m + 1, nt - m + 1);
            let eps = if m == 0 { 1.0 } else { 2.0 };
            // alpha_iso = -j rho c V, V_n = u_n eps 2/(2n+1) P_n^m(0)
            let a_iso: Vec<Cx> = (0..n_h)
                .map(|i| {
                    let n = m + i;
                    Cx::new(0.0, -rho_c * self.u[n] * eps * 2.0 / (2 * n + 1) as f64 * tab.at_cap[i])
                })
                .collect();
            // The torso's basis at the head's nodes: d/dr_H and value; then projected on the head's orders.
            let mut k_ht = vec![Cx::new(0.0, 0.0); n_h * n_t];
            let mut l_ht = vec![Cx::new(0.0, 0.0); n_h * n_t];
            let mut der = vec![Cx::new(0.0, 0.0); qh];
            let mut val = vec![Cx::new(0.0, 0.0); qh];
            for j in 0..n_t {
                let n = m + j;
                for q in 0..qh {
                    let p = tab.torso_at_head[j][q];
                    let dpt = tab.torso_at_head_d[j][q];
                    der[q] = gt[n][q].scale(self.head.rr[q] * p).add(ft[n][q].scale(self.head.rt[q] * dpt / self.head.r_other[q]));
                    val[q] = ft[n][q].scale(p);
                }
                for i in 0..n_h {
                    let row = &tab.head_at_head[i];
                    let (mut s, mut v) = (Cx::new(0.0, 0.0), Cx::new(0.0, 0.0));
                    for q in 0..qh {
                        let w = self.head.w[q] * row[q];
                        s = s.add(der[q].scale(w));
                        v = v.add(val[q].scale(w));
                    }
                    k_ht[i * n_t + j] = s;
                    l_ht[i * n_t + j] = v;
                }
            }
            // The head's basis at the torso's nodes: d/dr_T, projected on the torso's orders.
            let mut k_th = vec![Cx::new(0.0, 0.0); n_t * n_h];
            let mut der = vec![Cx::new(0.0, 0.0); qt];
            for i in 0..n_h {
                let n = m + i;
                for q in 0..qt {
                    der[q] = gh[n][q]
                        .scale(self.torso.rr[q] * tab.head_at_torso[i][q])
                        .add(fh[n][q].scale(self.torso.rt[q] * tab.head_at_torso_d[i][q] / self.torso.r_other[q]));
                }
                for j in 0..n_t {
                    let row = &tab.torso_at_torso[j];
                    let mut s = Cx::new(0.0, 0.0);
                    for q in 0..qt {
                        s = s.add(der[q].scale(self.torso.w[q] * row[q]));
                    }
                    k_th[j * n_h + i] = s;
                }
            }
            // (k^2 I - K_HT K_TH) delta = K_HT K_TH alpha_iso
            let mut w = vec![Cx::new(0.0, 0.0); n_h * n_h];
            for i in 0..n_h {
                for c in 0..n_h {
                    let mut s = Cx::new(0.0, 0.0);
                    for j in 0..n_t {
                        s = s.add(k_ht[i * n_t + j].mul(k_th[j * n_h + c]));
                    }
                    w[i * n_h + c] = s;
                }
            }
            let mut rhs: Vec<Cx> = (0..n_h)
                .map(|i| (0..n_h).fold(Cx::new(0.0, 0.0), |s, c| s.add(w[i * n_h + c].mul(a_iso[c]))))
                .collect();
            let mut mat: Vec<Cx> = w.iter().map(|x| Cx::new(-x.re, -x.im)).collect();
            for i in 0..n_h {
                mat[i * n_h + i] = mat[i * n_h + i].add(Cx::new(k * k, 0.0));
            }
            solve(&mut mat, &mut rhs, n_h)?;
            let delta = rhs;
            // beta = -(1/k) K_TH (alpha_iso + delta)
            let total: Vec<Cx> = a_iso.iter().zip(&delta).map(|(x, y)| x.add(*y)).collect();
            let beta: Vec<Cx> = (0..n_t)
                .map(|j| {
                    let s = (0..n_h).fold(Cx::new(0.0, 0.0), |s, i| s.add(k_th[j * n_h + i].mul(total[i])));
                    s.scale(-1.0 / k)
                })
                .collect();
            if let Some(mic) = &self.mic {
                let mut s = Cx::new(0.0, 0.0);
                for i in 0..n_h {
                    s = s.add(delta[i].mul(mic_h[m + i]).scale(mic.p_h[m][i]));
                }
                for j in 0..n_t {
                    s = s.add(beta[j].mul(mic_t[m + j]).scale(mic.p_t[m][j]));
                }
                dp = dp.add(s.scale(mic.cos_m[m]));
            }
            // dc = delta h/h' + L_HT beta; dZ += (a^2 / A) sum dc lam P(0)
            let mut s = Cx::new(0.0, 0.0);
            for i in 0..n_h {
                let mut dc = delta[i].mul(on_head[m + i]);
                for j in 0..n_t {
                    dc = dc.add(l_ht[i * n_t + j].mul(beta[j]));
                }
                s = s.add(dc.scale(self.lam[m + i] * tab.at_cap[i]));
            }
            dz = dz.add(s.scale(self.a * self.a / self.area));
        }
        if dp.finite() && dz.finite() { Ok((dp, dz)) } else { Err(format!("the snowman's solve is not finite at {f_hz} Hz")) }
    }
}

/// `f_n = h_n(x) / h_n'(s)` and `g_n = h_n'(x) / h_n'(s)`, `n = 0..=last`, for a field point `x = k r` outside the sphere
/// `s = k R` (`x >= s`), by the head's ratio recurrences: `rho_n = h_n / h_{n-1}` (`rho_1 = 1/x + j`, `rho_{n+1} =
/// (2n+1)/x - 1/rho_n`), `q_n = h_n(x) / h_n(s)` (`q_0 = (s/x) e^{-j(x - s)}`), `d_n = h_n' / h_n = 1/rho_n - (n+1)/x`
/// (`d_0 = -rho_1`).
fn radial(x: f64, s: f64, last: usize, f: &mut [Cx], g: &mut [Cx]) {
    let mut q = Cx::expj_neg(x - s).scale(s / x);
    let mut rx = Cx::new(1.0 / x, 1.0);
    let mut rs = Cx::new(1.0 / s, 1.0);
    let mut dx = Cx::new(-rx.re, -rx.im);
    let mut ds = Cx::new(-rs.re, -rs.im);
    for n in 0..=last {
        if n >= 1 {
            if n >= 2 {
                let c = (2 * n - 1) as f64;
                rx = Cx::new(c / x, 0.0).sub(rx.recip());
                rs = Cx::new(c / s, 0.0).sub(rs.recip());
            }
            q = q.mul(rx.div(rs));
            dx = rx.recip().sub(Cx::new((n + 1) as f64 / x, 0.0));
            ds = rs.recip().sub(Cx::new((n + 1) as f64 / s, 0.0));
        }
        f[n] = q.div(ds);
        g[n] = q.mul(dx).div(ds);
    }
}

/// `P_n^m(x)` and `d P_n^m / d th` for `n = m..=top` (index `n - m`), normalised so that `int_{-1}^{1} (P_n^m)^2 dx = 1`,
/// no Condon-Shortley phase. Recurrences: `P_m^m = sqrt((2m+1)/(2m)) s P_{m-1}^{m-1}` from `P_0^0 = sqrt(1/2)`,
/// `P_{m+1}^m = sqrt(2m+3) x P_m^m`, `P_n^m = A (x P_{n-1}^m - B P_{n-2}^m)` with `A = sqrt((4n^2-1)/(n^2-m^2))`,
/// `B = sqrt(((n-1)^2-m^2)/(4(n-1)^2-1))`; and `d/dth = (n x P_n^m - sqrt((n^2-m^2)(2n+1)/(2n-1)) P_{n-1}^m) / s`
/// (0 where `s = 0`, never read there).
fn pbar(x: f64, m: usize, top: usize) -> (Vec<f64>, Vec<f64>) {
    let s = det::sqrt((1.0 - x * x).max(0.0));
    let mut pmm = det::sqrt(0.5);
    for i in 1..=m {
        pmm *= det::sqrt((2 * i + 1) as f64 / (2 * i) as f64) * s;
    }
    let mut p = Vec::with_capacity(top + 1 - m);
    p.push(pmm);
    if top > m {
        p.push(det::sqrt((2 * m + 3) as f64) * x * pmm);
    }
    for n in m + 2..=top {
        let (nf, mf) = (n as f64, m as f64);
        let aa = det::sqrt((4.0 * nf * nf - 1.0) / (nf * nf - mf * mf));
        let bb = det::sqrt(((nf - 1.0) * (nf - 1.0) - mf * mf) / (4.0 * (nf - 1.0) * (nf - 1.0) - 1.0));
        let v = aa * (x * p[n - m - 1] - bb * p[n - m - 2]);
        p.push(v);
    }
    let mut d = Vec::with_capacity(p.len());
    for n in m..=top {
        let (nf, mf) = (n as f64, m as f64);
        let mut t = nf * x * p[n - m];
        if n > m {
            t -= det::sqrt((nf * nf - mf * mf) * (2.0 * nf + 1.0) / (2.0 * nf - 1.0)) * p[n - m - 1];
        }
        d.push(if s > 0.0 { t / s } else { 0.0 });
    }
    (p, d)
}

/// `P_0..=P_last` at `x` (the head's recurrence).
fn legendre(x: f64, last: usize) -> Vec<f64> {
    let mut p = vec![1.0, x];
    for n in 1..last {
        let nf = n as f64;
        p.push(((2.0 * nf + 1.0) * x * p[n] - nf * p[n - 1]) / (nf + 1.0));
    }
    p.truncate(last + 1);
    p
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`, `n >= 1`, by Newton's method on `P_n` from `cos(pi (i + 3/4) / (n + 1/2))`.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = Vec::with_capacity(n);
    let mut ws = Vec::with_capacity(n);
    for i in 0..n {
        let mut z = det::cos(core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for j in 2..=n {
                let jf = j as f64;
                let p2 = ((2.0 * jf - 1.0) * z * p1 - (jf - 1.0) * p0) / jf;
                p0 = p1;
                p1 = p2;
            }
            if n == 0 {
                break;
            }
            dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let step = p1 / dp;
            z -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        xs.push(z);
        ws.push(2.0 / ((1.0 - z * z) * dp * dp));
    }
    (xs, ws)
}

/// Gaussian elimination with partial pivoting (the first largest pivot) on the `n x n` row-major `a`; `b` becomes the
/// solution.
fn solve(a: &mut [Cx], b: &mut [Cx], n: usize) -> Result<(), String> {
    for col in 0..n {
        let mut piv = col;
        let mut best = a[col * n + col].abs();
        for r in col + 1..n {
            let v = a[r * n + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if !(best > 0.0 && best.is_finite()) {
            return Err("the snowman's system is singular".into());
        }
        if piv != col {
            for c in 0..n {
                a.swap(col * n + c, piv * n + c);
            }
            b.swap(col, piv);
        }
        let inv = a[col * n + col].recip();
        for r in col + 1..n {
            let f = a[r * n + col].mul(inv);
            if f.re == 0.0 && f.im == 0.0 {
                continue;
            }
            for c in col..n {
                a[r * n + c] = a[r * n + c].sub(f.mul(a[col * n + c]));
            }
            b[r] = b[r].sub(f.mul(b[col]));
        }
    }
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in r + 1..n {
            s = s.sub(a[r * n + c].mul(b[c]));
        }
        b[r] = s.div(a[r * n + r]);
    }
    Ok(())
}

/// `f(0..n)` on up to [`MAX_WORKERS`] threads, each on one contiguous block; every value is computed alone, so the result
/// is the same for any thread count. The first refusal in bin order is returned.
pub(super) fn par_map<T, F>(n: usize, f: F) -> Result<Vec<T>, String>
where
    T: Copy + Send + Default,
    F: Fn(usize) -> Result<T, String> + Sync,
{
    let workers = std::thread::available_parallelism().map_or(1, |p| p.get()).clamp(1, MAX_WORKERS);
    let chunk = n.div_ceil(workers).max(1);
    let mut out = vec![T::default(); n];
    let results: Vec<Result<(), String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = out
            .chunks_mut(chunk)
            .enumerate()
            .map(|(c, slot)| {
                let f = &f;
                scope.spawn(move || -> Result<(), String> {
                    for (i, o) in slot.iter_mut().enumerate() {
                        *o = f(c * chunk + i)?;
                    }
                    Ok(())
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err("a snowman worker panicked".into()))).collect()
    });
    for r in results {
        r?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RHO: f64 = 1.1336;
    const C: f64 = 356.53;
    const MIC: [f64; 3] = [0.212_132, 0.0, 0.212_132];
    const MOUTH: f64 = 0.013_548;
    /// The print's KEMAR snowman.
    const KEMAR: TorsoSphere = TorsoSphere { radius_m: 0.169, neck_gap_m: 0.053 };

    fn snowman(torso: TorsoSphere, top: f64) -> Snowman {
        Snowman::new(0.087, MOUTH, torso, Some(MIC), RHO, C, top).unwrap()
    }

    /// The head alone at the microphone (the head's own series, at the snowman's head).
    fn head_alone(f: f64) -> Cx {
        let medium = super::super::RayleighMedium { density: RHO, sound_speed: C };
        super::super::Series::new(0.087, MOUTH, super::super::POLE, MIC, medium, f, 30)
            .unwrap()
            .pressure(f, 30)
            .unwrap()
    }

    #[test]
    fn the_legendre_tables_are_orthonormal_and_their_slopes_are_the_derivatives() {
        let (x, w) = gauss_legendre(80);
        assert!((w.iter().sum::<f64>() - 2.0).abs() < 1e-13);
        for m in [0, 1, 7, 30] {
            let cols: Vec<Vec<f64>> = x.iter().map(|&xi| pbar(xi, m, 60).0).collect();
            for i in 0..=(60 - m) {
                for j in 0..=(60 - m) {
                    let g: f64 = (0..x.len()).map(|q| w[q] * cols[q][i] * cols[q][j]).sum();
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!((g - want).abs() < 1e-11, "m {m} rows {i} {j}: {g}");
                }
            }
            for th in [0.4_f64, 2.1] {
                let e = 1e-6;
                let (p1, _) = pbar((th + e).cos(), m, 60);
                let (p0, _) = pbar((th - e).cos(), m, 60);
                let (_, d) = pbar(th.cos(), m, 60);
                for i in 0..d.len() {
                    let fd = (p1[i] - p0[i]) / (2.0 * e);
                    assert!((fd - d[i]).abs() < 1e-6 * (1.0 + d[i].abs()), "m {m} n {} th {th}: {fd} vs {}", m + i, d[i]);
                }
            }
        }
    }

    /// The snowman against the numpy prototype of the same law (`torso/proto.py`, 60 head and 100 torso orders past
    /// `k a` and `k b`, 60 more nodes), KEMAR geometry, cap of radius 13.548 mm, the head's test microphone and medium.
    #[test]
    fn the_snowman_matches_an_independent_evaluation() {
        let s = snowman(KEMAR, 3000.0);
        for (f, dp_re, dp_im, dz_re, dz_im) in GOLDEN {
            let (dp, dz) = s.change(f).unwrap();
            let (want_p, want_z) = (Cx::new(dp_re, dp_im), Cx::new(dz_re, dz_im));
            assert!(dp.sub(want_p).abs() < 1e-9 * head_alone(f).abs(), "{f} Hz: {dp:?} vs {want_p:?}");
            assert!(dz.sub(want_z).abs() < 1e-7 * want_z.abs(), "{f} Hz: {dz:?} vs {want_z:?}");
        }
    }

    /// Twice the margins change nothing a test or a page could read.
    #[test]
    fn the_snowman_is_converged() {
        let base = snowman(KEMAR, 10_000.0);
        let more =
            Snowman::with_margins(0.087, MOUTH, super::super::POLE, KEMAR, Some(MIC), RHO, C, 10_000.0, [40, 70, 30])
                .unwrap();
        for f in [50.0, 1000.0, 4000.0, 9900.0] {
            let (a, b) = (base.change(f).unwrap(), more.change(f).unwrap());
            let p = head_alone(f).add(b.0).abs();
            assert!(a.0.sub(b.0).abs() < 1e-9 * p, "{f} Hz mic: {:?} vs {:?}", a.0, b.0);
            assert!(a.1.sub(b.1).abs() < 1e-7 * b.1.abs(), "{f} Hz cap: {:?} vs {:?}", a.1, b.1);
        }
    }

    /// A vanishing torso leaves the head: the change falls like `b^3` (Rayleigh) and is gone at 0.1 mm.
    #[test]
    fn a_vanishing_torso_leaves_the_head() {
        let rel = |b: f64| -> f64 {
            let s = snowman(TorsoSphere { radius_m: b, neck_gap_m: 0.053 }, 1000.0);
            s.change(1000.0).unwrap().0.abs() / head_alone(1000.0).abs()
        };
        let (cm, mm, tenth) = (rel(1e-2), rel(1e-3), rel(1e-4));
        assert!((800.0..1250.0).contains(&(cm / mm)) && (800.0..1250.0).contains(&(mm / tenth)), "{cm} {mm} {tenth}");
        assert!(tenth < 1e-9, "{tenth}");
        // A 1 mm torso to 20 kHz needs fewer torso orders (36) than head orders (51): the orders past the torso's are
        // empty, and the change is still the small torso's (numpy, the torso on all 51 orders: 1.08e-5 of the head's
        // pressure at 15 kHz).
        let small = snowman(TorsoSphere { radius_m: 1e-3, neck_gap_m: 0.053 }, 20_000.0);
        assert!(small.nh > small.nt, "{} {}", small.nh, small.nt);
        let (dp, dz) = small.change(15_000.0).unwrap();
        assert!(dp.finite() && dz.finite() && dp.abs() < 1e-4 * head_alone(15_000.0).abs(), "{dp:?} {dz:?}");
    }

    /// A far torso leaves the head: past a few metres the change falls like `1 / D^2` (numpy: 5.23e-6 and 5.26e-8 of the
    /// head's pressure at 1 kHz for gaps of 50 and 500 m), so the snowman returns today's single sphere.
    #[test]
    fn a_far_torso_leaves_the_head() {
        let rel = |h: f64| -> (f64, f64) {
            let s = snowman(TorsoSphere { radius_m: 0.169, neck_gap_m: h }, 1000.0);
            let (dp, dz) = s.change(1000.0).unwrap();
            (dp.abs() / head_alone(1000.0).abs(), dz.abs())
        };
        let (near, far) = (rel(50.0), rel(500.0));
        assert!((90.0..110.0).contains(&(near.0 / far.0)), "{near:?} {far:?}");
        assert!((90.0..110.0).contains(&(near.1 / far.1)), "{near:?} {far:?}");
        assert!(far.0 < 1e-7, "{far:?}");
    }

    /// The print's own delay (Algazi, Duda and Thompson 2002, p. 5): for a source on the ray from the torso's centre
    /// through the ear, the torso reflection trails the main pulse by about `2 (sqrt(D^2 + a^2) - b) / c = 0.89 ms` at
    /// 343 m/s (their multipole solution's peak separation at 75 degrees is "quite close" to it). By reciprocity the
    /// ear's response to that source is a small cap at the ear (on the equator, as the print puts the ears) heard far out
    /// on the ray; the torso's part trails the head's alone by its group delay over 1-4 kHz (numpy: 0.892 ms; read off
    /// the 5 kHz impulse response the two pulses overlap, the main one being 0.8 ms wide). Loose check: within 5 %.
    #[test]
    fn the_torso_reflection_trails_by_the_prints_delay() {
        let (a, c, rho) = (0.087, 343.0, 1.2);
        let d = a + KEMAR.neck_gap_m + KEMAR.radius_m;
        let len = (a * a + d * d).sqrt();
        let (r, ux, uz) = (2.0, a / len, d / len);
        // snowman frame: ear (a, 0, 0) + R u; card frame: x = y_s, y = z_s, z = x_s - a
        let mic = [0.0, r * uz, r * ux];
        let s = Snowman::new(a, 0.003, KEMAR, Some(mic), rho, c, 4000.0).unwrap();
        let mic_r = ((a + r * ux).powi(2) + (r * uz).powi(2)).sqrt();
        let series = super::super::Series::new(a, 0.003, super::super::POLE, mic,
            super::super::RayleighMedium { density: rho, sound_speed: c }, 4000.0, 30).unwrap();
        assert!((mic_r - series.r).abs() < 1e-12);
        let mut freqs = Vec::new();
        let mut phase = Vec::new();
        let mut last = 0.0;
        let mut unwrap = 0.0;
        let mut f = 1000.0;
        while f <= 4000.0 {
            let ratio = s.change(f).unwrap().0.div(series.pressure(f, 30).unwrap());
            let ph = ratio.im.atan2(ratio.re);
            if !freqs.is_empty() {
                while ph + unwrap - last > core::f64::consts::PI {
                    unwrap -= core::f64::consts::TAU;
                }
                while ph + unwrap - last < -core::f64::consts::PI {
                    unwrap += core::f64::consts::TAU;
                }
            }
            last = ph + unwrap;
            freqs.push(core::f64::consts::TAU * f);
            phase.push(last);
            f += 10.0;
        }
        let n = freqs.len() as f64;
        let (mx, my) = (freqs.iter().sum::<f64>() / n, phase.iter().sum::<f64>() / n);
        let sxy: f64 = freqs.iter().zip(&phase).map(|(x, y)| (x - mx) * (y - my)).sum();
        let sxx: f64 = freqs.iter().map(|x| (x - mx) * (x - mx)).sum();
        let delay = -sxy / sxx;
        let printed = 2.0 * (len - KEMAR.radius_m) / c;
        assert!((delay / printed - 1.0).abs() < 0.05, "group delay {delay} s against the print's {printed} s");
    }

    /// The snowman against a boundary-element solve of the same two rigid spheres: fs-bem, Burton-Miller under the dense
    /// second-order `PanelRule::Triangle`, a level-3 head icosphere (1280 flat panels) and a level-4 torso (5120), KEMAR
    /// geometry, the test medium. The cap is the head's panels within 0.3 rad of an icosahedron vertex (a 5-fold symmetric
    /// patch) at unit volume flow; the snowman gets a cap of the same area. The torso's change is read on both sides as
    /// (both spheres) - (head alone), at the C3 microphone and as the cap's mean pressure. fs-bem's time factor is
    /// `e^{-i w t}`, so its pressures are conjugated (the drive is real).
    ///
    /// The reference was chosen by refinement against the 3 % bound declared before any run (the bound did not move).
    /// Relative error of the BEM's torso change against the snowman at each level's own mouth, mic / cap, in %:
    ///
    /// | rule, (head, torso) level | 250 Hz      | 600 Hz      | 1200 Hz     |
    /// |---------------------------|-------------|-------------|-------------|
    /// | Centroid (3,3)            | 1.25 / 1.54 | 2.03 / 3.06 | 3.32 / 5.12 |
    /// | Centroid (3,4)            | 0.78 / 1.04 | 1.12 / 2.07 | 1.56 / 3.39 |
    /// | Centroid (4,3)            | 1.18 / 1.40 | 1.78 / 2.44 | 2.98 / 3.88 |
    /// | Triangle (3,3)            | 0.62 / 0.78 | 1.15 / 1.08 | 1.85 / 2.03 |
    /// | Triangle (3,4), this test | 0.16 / 0.15 | 0.31 / 0.18 | 0.42 / 0.33 |
    /// | Triangle (4,3)            | 0.61 / 0.83 | 1.13 / 1.18 | 1.90 / 2.33 |
    ///
    /// The default `PanelRule::Centroid` is first order, as its own doc says: one point kernel per off-diagonal panel.
    /// Refining the torso cuts the Triangle error 3.7-6.2x per level, within 4 degrees of straight toward the snowman.
    /// Refining the head moves it by at most 0.3 % (the cap patch changes with it, 30 -> 115 panels). (4,4) is past
    /// fs-bem's dense cap. At 2400 Hz, Triangle (3,4) reads 1.68 / 1.57; at 4800 Hz fs-bem refuses that mesh as too
    /// coarse (5.2 panels per wavelength). Above 2400 Hz the torso's change rests on the tests above.
    #[test]
    fn the_snowman_matches_a_boundary_element_solve() {
        against_boundary_elements(0.0, &[250.0, 600.0, 1200.0, 2400.0]);
    }

    /// The same check with the cap tilted up off the equator (N3b: a nostril's cap on the ray through its place on the
    /// baffle, `cos theta_c = 0.3`, about 17.5 degrees up), so `P_n^m(cos theta_c)` carries the cap: the mesh's cap patch
    /// stays at its icosahedron vertex and the snowman's up turns toward it, the torso with it. The bound is the
    /// equator's 3 %, declared before any run.
    #[test]
    fn a_tilted_cap_matches_a_boundary_element_solve() {
        against_boundary_elements(0.3, &[600.0, 1200.0]);
    }

    /// The snowman against fs-bem for a cap at `cos theta_c = cos_cap` (azimuth 0), at `freqs` (the doc above).
    fn against_boundary_elements(cos_cap: f64, freqs: &[f64]) {
        use fs_bem::helmholtz::{Formulation, Medium, PanelRule, solve_radiation_with_rule, with_dense_workers};
        use fs_bem::near_field::{Geometry, Options};
        use fs_bem::panel3d::SpherePanels;
        use fs_math::c64::C64;
        const TOL: f64 = 0.03;
        let (a, b) = (0.087, KEMAR.radius_m);
        let d = a + KEMAR.neck_gap_m + b;
        let unit = |v: [f64; 3]| {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let phi = 0.5 * (1.0 + 5.0_f64.sqrt());
        // Snowman frame in the mesh's coordinates: the cap's centre = an icosahedron vertex of fs-rep-mesh's icosphere,
        // up at cos_cap to it (perpendicular at 0), front = the cap's centre's horizontal part (the cap's centre itself at
        // 0), side = up x front. At cos_cap 0 every vector is the equator test's, bit for bit.
        let centre = unit([phi, 0.0, 1.0]);
        let up0 = unit([-1.0, 0.0, phi]);
        let sin_cap = (1.0 - cos_cap * cos_cap).sqrt();
        let up = [0, 1, 2].map(|i| cos_cap * centre[i] + sin_cap * up0[i]);
        let front = [0, 1, 2].map(|i| (centre[i] - cos_cap * up[i]) / sin_cap);
        let side = [up[1] * front[2] - up[2] * front[1], up[2] * front[0] - up[0] * front[2], up[0] * front[1] - up[1] * front[0]];
        let world = |s: [f64; 3]| -> [f64; 3] {
            [0, 1, 2].map(|i| s[0] * front[i] + s[1] * side[i] + s[2] * up[i])
        };
        // The cap's axis in the card frame (card z = front, x = side, y = up): POLE at cos_cap 0.
        let axis = [0.0, cos_cap, sin_cap];
        let head = SpherePanels::icosphere(a, 3).unwrap();
        let head_tris = head.triangles().unwrap().to_vec();
        let torso_at = world([0.0, 0.0, -d]);
        let mut both_tris = head_tris.clone();
        for t in SpherePanels::icosphere(b, 4).unwrap().triangles().unwrap() {
            both_tris.push((*t).map(|p| [p[0] + torso_at[0], p[1] + torso_at[1], p[2] + torso_at[2]]));
        }
        let both = SpherePanels::from_triangles(both_tris).unwrap();
        let cap: Vec<usize> = (0..head_tris.len())
            .filter(|&i| {
                let c = unit(head.centroids()[i]);
                dot(c, centre).clamp(-1.0, 1.0).acos() < 0.3
            })
            .collect();
        assert!(!cap.is_empty() && cap.len() % 5 == 0, "the cap patch should be 5-fold symmetric: {} panels", cap.len());
        let area: f64 = cap.iter().map(|&i| head.areas()[i]).sum();
        let mut v_head = vec![C64::new(0.0, 0.0); head_tris.len()];
        for &i in &cap {
            v_head[i] = C64::new(1.0 / area, 0.0);
        }
        let mut v_both = v_head.clone();
        v_both.resize(both.areas().len(), C64::new(0.0, 0.0));
        let medium = Medium { density: RHO, sound_speed: C };
        let mic_s = [MIC[2] + a, MIC[0], MIC[1]];
        let mic = world(mic_s);
        let mouth = (area / core::f64::consts::PI).sqrt();
        let ours = Snowman::at_cap(a, mouth, axis, KEMAR, Some(MIC), RHO, C, 2400.0).unwrap();
        for &f in freqs {
            let k = core::f64::consts::TAU * f / C;
            let read = |surface: &SpherePanels, v: &[C64]| -> (C64, C64) {
                let n = surface.areas().len();
                let rule = PanelRule::Triangle { max_pairs: n * n };
                let sol = with_dense_workers(8, || {
                    solve_radiation_with_rule(surface, k, medium, v, Formulation::BurtonMiller, rule)
                })
                .unwrap();
                let at_mic = Geometry::new(surface, &[mic], 1e-3).unwrap().prepare(k, medium, Options::default()).unwrap()
                    .evaluate(&sol).unwrap().pressure[0];
                let mut mean = C64::new(0.0, 0.0);
                for &i in &cap {
                    mean = mean + sol.pressure[i].scale(head.areas()[i] / area);
                }
                (at_mic, mean)
            };
            let (p_head, z_head) = read(&head, &v_head);
            let (p_both, z_both) = read(&both, &v_both);
            let dp_bem = (p_both - p_head).conj();
            let dz_bem = (z_both - z_head).conj();
            let (dp, dz) = ours.change(f).unwrap();
            let e_mic = dp.sub(Cx::new(dp_bem.re, dp_bem.im)).abs() / dp.abs();
            let e_cap = dz.sub(Cx::new(dz_bem.re, dz_bem.im)).abs() / dz.abs();
            let series = super::super::Series::new(a, mouth, axis, MIC,
                super::super::RayleighMedium { density: RHO, sound_speed: C }, f, 30).unwrap().pressure(f, 30).unwrap();
            let e_head = series.sub(Cx::new(p_head.re, -p_head.im)).abs() / series.abs();
            eprintln!("BEM cos_cap {cos_cap} {f} Hz ({} cap panels, mouth {mouth:.6} m): torso change at the mic {e_mic:.3e}, \
                on the cap {e_cap:.3e}; head alone at the mic {e_head:.3e}", cap.len());
            assert!(e_mic < TOL && e_cap < TOL, "{f} Hz: mic {e_mic:.3e} ({dp:?} vs {dp_bem:?}), cap {e_cap:.3e} ({dz:?} vs {dz_bem:?})");
        }
    }

    /// The addition theorem the cap's velocity rests on, off the equator (N3b): `P_n(cos gamma) = sum_m eps_m 2 / (2n+1)
    /// P_n^m(cos th) P_n^m(cos th_c) cos(m (phi - phi_c))` with these tables, for caps at and above the equator.
    #[test]
    fn the_addition_theorem_holds_for_a_cap_off_the_equator() {
        let half_pi = core::f64::consts::FRAC_PI_2;
        for (th, phi, th_c, phi_c) in [(0.7, 0.3, 1.2662, 0.0816), (2.2, -1.1, half_pi, 0.0), (1.9, 2.5, 0.5, -0.4)] {
            let (x, x_c) = (f64::cos(th), f64::cos(th_c));
            let cos_gamma = x * x_c + f64::sin(th) * f64::sin(th_c) * f64::cos(phi - phi_c);
            let top = 40;
            let want = legendre(cos_gamma, top);
            let mut got = vec![0.0; top + 1];
            for m in 0..=top {
                let (p, p_c) = (pbar(x, m, top).0, pbar(x_c, m, top).0);
                let eps = if m == 0 { 1.0 } else { 2.0 };
                for n in m..=top {
                    got[n] += eps * 2.0 / (2 * n + 1) as f64 * p[n - m] * p_c[n - m] * f64::cos(m as f64 * (phi - phi_c));
                }
            }
            for n in 0..=top {
                let (g, w) = (got[n], want[n]);
                assert!((g - w).abs() < 1e-12, "n {n} at {th} {phi} / {th_c} {phi_c}: {g} vs {w}");
            }
        }
    }

    /// The torso is coaxial with the polar axis, so a cap turned in azimuth (on the equator) with the microphone turned
    /// the same way about that axis is the lips' cap: the microphone's change agrees within rounding, and the cap's own
    /// change (which reads no microphone) is the lips' bit for bit. A cap facing back is refused.
    #[test]
    fn a_cap_turned_about_the_torso_is_the_lips_cap_turned() {
        let base = snowman(KEMAR, 4000.0);
        for turn in [0.0816_f64, -0.3, 0.9] {
            let axis = [turn.sin(), 0.0, turn.cos()];
            // The microphone about the head's centre in the snowman frame (x_s = z + a, y_s = x), turned by `turn`.
            let (xs, ys) = (MIC[2] + 0.087, MIC[0]);
            let (xt, yt) = (xs * turn.cos() - ys * turn.sin(), xs * turn.sin() + ys * turn.cos());
            let mic = [yt, MIC[1], xt - 0.087];
            let turned = Snowman::at_cap(0.087, MOUTH, axis, KEMAR, Some(mic), RHO, C, 4000.0).unwrap();
            for f in [300.0, 1500.0, 4000.0] {
                let ((dp0, dz0), (dp, dz)) = (base.change(f).unwrap(), turned.change(f).unwrap());
                assert!(dp.sub(dp0).abs() < 1e-11 * dp0.abs(), "turn {turn} {f} Hz: {dp:?} vs {dp0:?}");
                assert_eq!((dz.re.to_bits(), dz.im.to_bits()), (dz0.re.to_bits(), dz0.im.to_bits()), "turn {turn} {f} Hz");
            }
        }
        let err = Snowman::at_cap(0.087, MOUTH, [0.0, 0.0, -1.0], KEMAR, Some(MIC), RHO, C, 1000.0).unwrap_err();
        assert!(err.contains("facing front"), "{err}");
    }

    #[test]
    fn non_physical_snowmen_are_refused() {
        for (t, why) in [
            (TorsoSphere { radius_m: 0.0, neck_gap_m: 0.05 }, "torso radius"),
            (TorsoSphere { radius_m: f64::NAN, neck_gap_m: 0.05 }, "torso radius"),
            (TorsoSphere { radius_m: 0.169, neck_gap_m: 0.0 }, "must not touch"),
        ] {
            let err = Snowman::new(0.087, MOUTH, t, Some(MIC), RHO, C, 1000.0).unwrap_err();
            assert!(err.contains(why), "{err}");
        }
        // A microphone inside the torso, and a band past the largest order, are refused by name.
        let err = Snowman::new(0.087, MOUTH, KEMAR, Some([0.0, -0.309, 0.0]), RHO, C, 1000.0).unwrap_err();
        assert!(err.contains("outside the head and the torso"), "{err}");
        let err = Snowman::new(0.087, MOUTH, KEMAR, Some(MIC), RHO, C, 200_000.0).unwrap_err();
        assert!(err.contains("more than"), "{err}");
        let s = snowman(KEMAR, 1000.0);
        assert!(s.change(2000.0).unwrap_err().contains("not prepared"));
    }

    #[test]
    fn par_map_is_the_same_for_any_split() {
        let f = |i: usize| -> Result<f64, String> { Ok((i as f64).sqrt()) };
        let v = par_map(1001, f).unwrap();
        assert!(v.iter().enumerate().all(|(i, x)| *x == (i as f64).sqrt()));
        let err = par_map(100, |i| if i == 37 || i == 90 { Err(format!("bin {i}")) } else { Ok(0.0) }).unwrap_err();
        assert_eq!(err, "bin 37");
    }

    /// `(f, Re dp, Im dp, Re dZ, Im dZ)` [Hz, Pa s/m^3, `e^{+j w t}`]: `torso/proto.py` (numpy), KEMAR snowman, head
    /// 8.7 cm, cap of radius 13.548 mm, microphone (0.212132, 0, 0.212132) m, rho 1.1336, c 356.53.
    const GOLDEN: [(f64, f64, f64, f64, f64); 3] = [
        (300.0, 8.027145248952766e+01, 2.187268083485177e+00, 1.145493680745585e+02, 1.065856194033701e+02),
        (1000.0, 4.953475683962679e+02, -1.655421611640986e+02, -9.066309050057920e+02, 1.462277558411516e+02),
        (3000.0, 1.122316775501993e+03, 1.013006861890167e+03, -6.720062911932004e+02, -3.072372226841324e+03),
    ];
}
