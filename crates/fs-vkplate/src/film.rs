//! A head film's loss as a generalized-Maxwell (Prony) law on a plate's
//! bending operator, so that every grid mode carries it, including those above
//! the radiation table's shapes (.370 option C).
//!
//! The law is fitted to the film card's loss-tangent curve alone, never to a
//! recording. A causal loss makes the storage modulus vary with frequency
//! (Kramers-Kronig); the law's storage modulus equals the card's E at
//! [`matsuzawa_anchor_hz`]. Time stepping is fs-material's exact linear-ramp
//! branch recursion ([`GeneralizedMaxwell::step`]), applied to the bending
//! operator's field: with p_j = K q_j,
//! `p_j <- a_j p_j + b_j (K w^n - K w^(n-1))`, `a_j = exp(-dt/tau_j)` and
//! `b_j = (1 - a_j) tau_j / dt`. The bending force is then
//! `-(e_inf K w + sum_j e_j p_j)`. Only the bending operator is viscoelastic:
//! the tension (a pre-stress) and the Berger in-plane stretching stay elastic.
//!
//! Discrete energy. Leapfrog with this recursion satisfies, per branch,
//! `p^n . (w^(n+1) - w^(n-1)) / 2 = Phi^(n+1/2) - Phi^(n-1/2) + D^n`, with
//! `Phi^(n+1/2) = (1 + a)/(4 b) |q^n|_K^2 + q^n . K (w^(n+1) - w^n) / 2` and
//! `D^n = (dt/tau)/4 |q^(n-1) + q^n|_K^2 >= 0`, the midpoint heat. The choice
//! `(1 + a)/(4 b)` is the only coefficient of `|q|^2` that makes the
//! remainder a nonnegative form. That coefficient tends to the physical 1/2
//! as dt/tau -> 0, and the excess is numerical stored energy. The heat this
//! module reports is `D`, which closes `Body::discrete_energy` exactly. The
//! closure check is in `film_tests`.

use crate::{Csr, dot};
use fs_material::visco::{GeneralizedMaxwell, nnls};
use fs_math::det;
use std::f64::consts::PI;

/// The band of Matsuzawa 1959's resonance measurements on the PET film card
/// ("from about 100 to 4000 c/s"), over which he found the film's modulus
/// almost independent of frequency.
pub const MATSUZAWA_BAND_HZ: (f64, f64) = (100.0, 4000.0);

/// Relaxation ladder density: this many candidate branches per decade.
pub const BRANCHES_PER_DECADE: f64 = 3.0;

/// The ladder reaches this many decades beyond each end of the curve's band.
pub const LADDER_MARGIN_DECADES: f64 = 1.0;

/// Log-spaced points of the curve's band in the least-squares fit and in the
/// verification of its certificate.
pub const FIT_POINTS: usize = 400;
pub const VERIFY_POINTS: usize = 2000;

/// Where a film law's storage modulus equals the card's E [Hz]: the log-centre
/// of [`MATSUZAWA_BAND_HZ`].
pub fn matsuzawa_anchor_hz() -> f64 {
    det::sqrt(MATSUZAWA_BAND_HZ.0 * MATSUZAWA_BAND_HZ.1)
}

/// A film's generalized-Maxwell law, with moduli relative to the card's E.
#[derive(Debug, Clone, PartialEq)]
pub struct FilmLaw {
    /// `e_inf`: the relaxed (0 Hz) storage modulus over E. `terms`:
    /// `(E_j / E, tau_j [s])`, the branches the fit kept.
    pub law: GeneralizedMaxwell,
    /// The curve's band [Hz]; the law's tan delta is fitted and verified there.
    pub band_hz: (f64, f64),
    /// The storage modulus is exactly E here [Hz].
    pub anchor_hz: f64,
    /// Candidate branches on the ladder (`law.terms` keeps the nonzero ones).
    pub ladder: usize,
    /// Measured `sup |tan delta_law / tan delta_curve - 1|`, over
    /// `verification_points` log-spaced points of the band.
    pub sup_rel_tan_delta_err: f64,
    pub verification_points: usize,
}

impl FilmLaw {
    pub fn tan_delta(&self, hz: f64) -> f64 {
        self.law.loss_factor(2.0 * PI * hz)
    }

    /// Storage modulus over the card's E at `hz`.
    pub fn storage(&self, hz: f64) -> f64 {
        self.law.modulus(2.0 * PI * hz).0
    }

    /// The instantaneous (infinite-frequency) modulus over E: the stiffness
    /// the explicit step's admission must allow for.
    pub fn glassy(&self) -> f64 {
        self.law.e_inf + self.law.terms.iter().map(|t| t.0).sum::<f64>()
    }

    /// The law's receipt: branches, fit, certificate, anchor and dispersion.
    pub fn receipt(&self) -> String {
        let (lo, hi) = self.band_hz;
        let relax = |tau: f64| 1.0 / (2.0 * PI * tau);
        let fast = self
            .law
            .terms
            .iter()
            .map(|t| relax(t.1))
            .fold(0.0_f64, f64::max);
        let slow = self
            .law
            .terms
            .iter()
            .map(|t| relax(t.1))
            .fold(f64::INFINITY, f64::min);
        format!(
            "film law: generalized Maxwell, {} of {} branches kept (relaxation {slow:.1} Hz to {fast:.0} Hz; the ladder holds {BRANCHES_PER_DECADE} per decade from {LADDER_MARGIN_DECADES} decade below to {LADDER_MARGIN_DECADES} decade above the curve's band); fitted by non-negative least squares to the head card's tan delta curve alone over {lo}-{hi} Hz, never to a recording; sup |tan delta / curve - 1| {:.2} % on {} log-spaced points; storage modulus = the card's E at {:.1} Hz (the log-centre of Matsuzawa 1959's {}-{} c/s, where he found it almost constant); causality then gives E' {:+.1} % at {lo} Hz and {:+.1} % at {hi} Hz, relaxed x{:.4}, glassy x{:.4}",
            self.law.terms.len(),
            self.ladder,
            100.0 * self.sup_rel_tan_delta_err,
            self.verification_points,
            self.anchor_hz,
            MATSUZAWA_BAND_HZ.0,
            MATSUZAWA_BAND_HZ.1,
            100.0 * (self.storage(lo) - 1.0),
            100.0 * (self.storage(hi) - 1.0),
            self.law.e_inf,
            self.glassy(),
        )
    }
}

fn log_point(lo: f64, hi: f64, i: usize, n: usize) -> f64 {
    if i == 0 {
        lo
    } else if i + 1 == n {
        hi
    } else {
        det::exp(det::ln(lo) + (det::ln(hi) - det::ln(lo)) * i as f64 / (n - 1) as f64)
            .clamp(lo, hi)
    }
}

/// Fit a film law to `tan_delta` (the card's curve: `Ok(None)` outside its
/// band) over `band_hz`, with the storage modulus pinned to the card's E at
/// `anchor_hz`. The relative tan delta error weighted by E' is minimized over
/// nonnegative branch moduli on a fixed ladder. This is linear, because
/// `E'' - tan_delta E'` is linear in the branch moduli once `E'(anchor) = 1`
/// fixes `e_inf`. The certificate is measured on a denser grid.
pub fn fit(
    tan_delta: &dyn Fn(f64) -> Result<Option<f64>, String>,
    band_hz: (f64, f64),
    anchor_hz: f64,
) -> Result<FilmLaw, String> {
    let (lo, hi) = band_hz;
    if !(lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo) {
        return Err(format!("film law: band {lo}-{hi} Hz"));
    }
    if !(anchor_hz >= lo && anchor_hz <= hi) {
        return Err(format!(
            "film law: anchor {anchor_hz} Hz lies outside the curve's band {lo}-{hi} Hz"
        ));
    }
    let curve = |hz: f64| -> Result<f64, String> {
        match tan_delta(hz)? {
            Some(t) if t.is_finite() && t > 0.0 && t < 1.0 => Ok(t),
            Some(t) => Err(format!(
                "film law: the curve's tan delta {t} at {hz} Hz is not in (0, 1)"
            )),
            None => Err(format!(
                "film law: the curve has no value at {hz} Hz, inside its own band"
            )),
        }
    };
    let ln10 = det::ln(10.0);
    let (l0, l1) = (
        det::ln(lo) - LADDER_MARGIN_DECADES * ln10,
        det::ln(hi) + LADDER_MARGIN_DECADES * ln10,
    );
    let ladder = ((l1 - l0) / ln10 * BRANCHES_PER_DECADE).round() as usize + 1;
    let tau: Vec<f64> = (0..ladder)
        .map(|j| 1.0 / (2.0 * PI * det::exp(l0 + (l1 - l0) * j as f64 / (ladder - 1) as f64)))
        .collect();
    let response = |omega: f64, tau: f64| {
        let x = omega * tau;
        let den = 1.0 + x * x;
        (x * x / den, x / den)
    };
    let anchored: Vec<f64> = tau
        .iter()
        .map(|&t| response(2.0 * PI * anchor_hz, t).0)
        .collect();
    // Row i: (E''(w_i) - t_i E'(w_i)) / t_i with E' = 1 + sum_j (storage_j - anchored_j) e_j.
    let mut a = vec![0.0; FIT_POINTS * ladder];
    for i in 0..FIT_POINTS {
        let hz = log_point(lo, hi, i, FIT_POINTS);
        let t = curve(hz)?;
        for (j, &tj) in tau.iter().enumerate() {
            let (storage, loss) = response(2.0 * PI * hz, tj);
            a[i * ladder + j] = (loss - t * (storage - anchored[j])) / t;
        }
    }
    // nnls's entry threshold is absolute: give every column unit norm.
    let norms: Vec<f64> = (0..ladder)
        .map(|j| {
            det::sqrt(
                (0..FIT_POINTS)
                    .map(|i| a[i * ladder + j] * a[i * ladder + j])
                    .sum::<f64>(),
            )
        })
        .collect();
    if norms.iter().any(|n| !(n.is_finite() && *n > 0.0)) {
        return Err("film law: a relaxation branch has no response in the curve's band".into());
    }
    for i in 0..FIT_POINTS {
        for j in 0..ladder {
            a[i * ladder + j] /= norms[j];
        }
    }
    let x = nnls(&a, &vec![1.0; FIT_POINTS], FIT_POINTS, ladder);
    let moduli: Vec<f64> = x.iter().zip(&norms).map(|(x, n)| x / n).collect();
    let e_inf = 1.0
        - moduli
            .iter()
            .zip(&anchored)
            .map(|(e, an)| e * an)
            .sum::<f64>();
    if !(e_inf > 0.0 && e_inf.is_finite()) {
        return Err(format!(
            "film law: the fit leaves a relaxed modulus of {e_inf} E"
        ));
    }
    let terms: Vec<(f64, f64)> = moduli
        .iter()
        .zip(&tau)
        .filter(|(e, _)| **e > 0.0)
        .map(|(e, t)| (*e, *t))
        .collect();
    if terms.is_empty() {
        return Err("film law: the fit kept no branch".into());
    }
    let law = GeneralizedMaxwell::new(e_inf, terms).map_err(|e| e.to_string())?;
    let mut sup = 0.0_f64;
    for i in 0..VERIFY_POINTS {
        let hz = log_point(lo, hi, i, VERIFY_POINTS);
        sup = sup.max((law.loss_factor(2.0 * PI * hz) / curve(hz)? - 1.0).abs());
    }
    Ok(FilmLaw {
        law,
        band_hz,
        anchor_hz,
        ladder,
        sup_rel_tan_delta_err: sup,
        verification_points: VERIFY_POINTS,
    })
}

/// One branch's step coefficients at clock `k`: `a = exp(-r)`, `b = (1 - a)/r`
/// (as fs-material's ramp), `r = k / tau`, and the discrete stored-energy
/// coefficient `(1 + a)/(4 b)`.
#[derive(Clone, Copy, Debug)]
struct Branch {
    a: f64,
    b: f64,
    r: f64,
    stored: f64,
}

/// A body's film-law state: the law, its clock's coefficients and the branch
/// fields `p_j = K q_j` over the body's nodes.
#[derive(Clone)]
pub(crate) struct FilmField {
    pub(crate) law: FilmLaw,
    e_inf: f64,
    moduli: Vec<f64>,
    k: f64,
    coef: Vec<Branch>,
    p: Vec<Vec<f64>>,
    /// `q_j' K q_j` per branch.
    s: Vec<f64>,
    /// K w at this step (filled by `linear_forces`) and at the previous one.
    pub(crate) sw: Vec<f64>,
    prev: Vec<f64>,
    delta: Vec<f64>,
    y: Vec<f64>,
    primed: bool,
    pub(crate) heat_j: f64,
    pub(crate) last_j: f64,
}

impl FilmField {
    pub(crate) fn new(law: FilmLaw, n: usize) -> FilmField {
        let moduli: Vec<f64> = law.law.terms.iter().map(|t| t.0).collect();
        let branches = moduli.len();
        FilmField {
            e_inf: law.law.e_inf,
            moduli,
            law,
            k: 0.0,
            coef: vec![],
            p: vec![vec![0.0; n]; branches],
            s: vec![0.0; branches],
            sw: vec![0.0; n],
            prev: vec![0.0; n],
            delta: vec![0.0; n],
            y: vec![0.0; n],
            primed: false,
            heat_j: 0.0,
            last_j: 0.0,
        }
    }

    fn prepare(&mut self, k: f64) {
        if self.k == k && self.coef.len() == self.moduli.len() {
            return;
        }
        self.coef = self
            .law
            .law
            .terms
            .iter()
            .map(|&(_, tau)| {
                let r = k / tau;
                let a = det::exp(-r);
                let b = -det::expm1(-r) / r;
                Branch {
                    a,
                    b,
                    r,
                    stored: (1.0 + a) / (4.0 * b),
                }
            })
            .collect();
        self.k = k;
    }

    /// The branches relax (initial conditions set directly on the body).
    pub(crate) fn reset(&mut self) {
        for p in &mut self.p {
            p.fill(0.0);
        }
        self.s.fill(0.0);
        self.primed = false;
    }

    /// Advance the branches from `wp` to `w` at clock `k`. `f` already holds
    /// `-K w` (the card's E) from `linear_forces`; this adds the law's
    /// difference from it. `prime` fills `K wp` on the first step after a reset.
    pub(crate) fn advance(&mut self, w: &[f64], wp: &[f64], f: &mut [f64], k: f64, stiff: &Csr) {
        self.prepare(k);
        let n = self.sw.len();
        if !self.primed {
            self.prev.fill(0.0);
            stiff.mul_add(&wp[..n], 1.0, &mut self.prev);
            self.primed = true;
        }
        self.step(w, wp, f);
    }

    /// The recursion itself, with `self.sw` = K w and `self.prev` = K wp.
    fn step(&mut self, w: &[f64], wp: &[f64], f: &mut [f64]) {
        let n = self.sw.len();
        for i in 0..n {
            self.delta[i] = w[i] - wp[i];
            self.y[i] = self.sw[i] - self.prev[i];
        }
        let dy = dot(&self.delta, &self.y);
        let relaxed = 1.0 - self.e_inf;
        for (f, s) in f[..n].iter_mut().zip(&self.sw) {
            *f += relaxed * s;
        }
        let mut heat = 0.0;
        for (j, p) in self.p.iter_mut().enumerate() {
            let Branch { a, b, r, .. } = self.coef[j];
            let e = self.moduli[j];
            let pd = dot(p, &self.delta);
            let s0 = self.s[j];
            // (r/4) |q^(n-1) + q^n|_K^2 with q^n = a q^(n-1) + b delta
            let mid = (1.0 + a) * (1.0 + a) * s0 + 2.0 * (1.0 + a) * b * pd + b * b * dy;
            heat += e * 0.25 * r * mid;
            self.s[j] = a * a * s0 + 2.0 * a * b * pd + b * b * dy;
            for ((p, y), f) in p.iter_mut().zip(&self.y).zip(f[..n].iter_mut()) {
                *p = a * *p + b * y;
                *f -= e * *p;
            }
        }
        self.last_j = heat;
        self.heat_j += heat;
        std::mem::swap(&mut self.prev, &mut self.sw);
    }

    /// The branches' discrete stored energy at the half step after the last
    /// one (`w`, `wp` the body's current pair), in joules.
    pub(crate) fn stored_j(&self, w: &[f64], wp: &[f64]) -> f64 {
        if !self.primed {
            return 0.0;
        }
        let delta: Vec<f64> = w.iter().zip(wp).map(|(a, b)| a - b).collect();
        self.p
            .iter()
            .enumerate()
            .map(|(j, p)| self.moduli[j] * (self.coef[j].stored * self.s[j] + 0.5 * dot(p, &delta)))
            .sum()
    }

    /// The branches' physical stored energy `sum_j E_j/2 |q_j|_K^2` [J].
    pub(crate) fn elastic_j(&self) -> f64 {
        self.moduli
            .iter()
            .zip(&self.s)
            .map(|(e, s)| 0.5 * e * s)
            .sum()
    }

    pub(crate) fn e_inf(&self) -> f64 {
        self.e_inf
    }

    #[cfg(test)]
    pub(crate) fn step_with(
        &mut self,
        k: f64,
        sw: &[f64],
        prev: Option<&[f64]>,
        w: &[f64],
        wp: &[f64],
        f: &mut [f64],
    ) {
        self.prepare(k);
        self.sw.copy_from_slice(sw);
        if let Some(prev) = prev {
            self.prev.copy_from_slice(prev);
            self.primed = true;
        }
        self.step(w, wp, f);
    }

    #[cfg(test)]
    pub(crate) fn branch(&self, j: usize) -> &[f64] {
        &self.p[j]
    }
}
