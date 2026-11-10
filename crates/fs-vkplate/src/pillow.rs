//! A drum pillow as a porous layer against the heads, from a measured JCAL
//! card (opt-in; the kick's `--kick-pillow`, MB 22635).
//!
//! The layer's normal surface impedance `Zs(w)` comes from fs-material's
//! JCAL law (`porous`) for a layer of declared thickness, backed either
//! rigidly or by a half-space of air: the bounding pair, since a pillow's
//! far side is neither. Only the resistance `Re Zs` acts. Over a band it is
//! fitted by the passive law
//!
//! `R(w) = R_inf + sum_j c_j / (1 + (w tau_j)^2)`,  `R_inf, c_j >= 0`,
//!
//! on a fixed ladder of relaxation times (non-negative least squares on the
//! relative error, against the card's curve alone, never a recording). Each
//! term is passive: `R_inf` is a dashpot per unit area, and each `c_j` a
//! Maxwell element, a spring `c_j / tau_j` in series with a dashpot `c_j`.
//! A relaxation ladder's resistance cannot rise with frequency, so where the
//! card's curve rises (a rigid-backed layer near its quarter-wave frequency)
//! the fit misses by its stated residual.
//!
//! Dropped, and named on the receipt: the layer's own reactance (its added
//! mass and stiffness; the law's reactance is only the compliance its
//! branches imply, the causal partner of their resistance), the frame's
//! motion (JCAL takes the frame rigid), and any pressure beyond the patch.
//!
//! Time stepping. `R_inf` joins the plate's centred node damping `s0` on the
//! patch, as a muffle dashpot does. The branches use the film module's exact
//! linear-ramp recursion on the diagonal operator `A = h^2` of the patch
//! nodes: `p_j <- a_j p_j + b_j A (w^n - w^(n-1))`, force `-e_j p_j`, with
//! `e_j = c_j / tau_j`, `a_j = exp(-dt/tau_j)` and `b_j = (1 - a_j) tau_j / dt`.
//! Its discrete energy and midpoint heat `e_j (dt/tau_j)/4 |q^(n-1) + q^n|_A^2`
//! close as the film law's do (see `film`).

use fs_material::porous::{PorousAmbient, PorousCard};
use fs_material::visco::nnls;
use fs_math::c64::C64;
use fs_math::det;
use std::f64::consts::PI;

/// The band the law is fitted and certified over [Hz]: down to the kick's
/// lowest line (44.5 Hz in the recording) and up to 5 kHz. Outside the
/// card's own band the JCAL law is extrapolated, and the receipt says so.
pub const FIT_BAND_HZ: (f64, f64) = (40.0, 5000.0);

/// Relaxation ladder density: this many candidate branches per decade.
pub const BRANCHES_PER_DECADE: f64 = 3.0;

/// The ladder reaches this many decades beyond each end of the fit band.
pub const LADDER_MARGIN_DECADES: f64 = 1.0;

/// Log-spaced points of the fit band in the least-squares fit, and in the
/// verification of its residual.
pub const FIT_POINTS: usize = 400;
pub const VERIFY_POINTS: usize = 2000;

/// The frequencies the receipt quotes the card's resistance at [Hz].
pub const QUOTED_HZ: [f64; 4] = [45.0, 250.0, 1000.0, 4000.0];

/// The measured cards a pillow can take.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PillowCard {
    /// FOAM 02's Pinta polyester fleece (a felt-class nonwoven), standing in
    /// for pillow fill: a borrowed material.
    PintaFleeceFoam02,
}

impl PillowCard {
    /// The command-line spellings, `|`-separated.
    pub const ROWS: &'static str = "pinta-fleece-foam02";

    pub fn parse(s: &str) -> Option<PillowCard> {
        match s {
            "pinta-fleece-foam02" => Some(PillowCard::PintaFleeceFoam02),
            _ => None,
        }
    }

    pub fn card(self) -> PorousCard {
        match self {
            PillowCard::PintaFleeceFoam02 => PorousCard::pinta_fleece_foam02(),
        }
    }

    /// The thickness of the sample the card's JCAL row was identified on [m]
    /// (FOAM 02's Pinta row: thickness 50 mm, diameter 98 mm).
    pub fn sample_thickness_m(self) -> f64 {
        match self {
            PillowCard::PintaFleeceFoam02 => 0.050,
        }
    }
}

/// What backs the layer on its far side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backing {
    Rigid,
    /// A half-space of air at the card's ambient (rho0 c0).
    Air,
}

impl Backing {
    /// The command-line spellings, `|`-separated.
    pub const ROWS: &'static str = "rigid|air";

    pub fn parse(s: &str) -> Option<Backing> {
        match s {
            "rigid" => Some(Backing::Rigid),
            "air" => Some(Backing::Air),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Backing::Rigid => "rigid-backed",
            Backing::Air => "air-backed",
        }
    }
}

/// A pillow, as declared: card, thickness, backing, and the patch it covers
/// on each head (centre x, y and radius [m]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pillow {
    pub card: PillowCard,
    pub thickness_m: f64,
    pub backing: Backing,
    pub patch: (f64, f64, f64),
}

/// The layer's normal surface impedance [Pa s/m] at `hz` (e^{+jwt}).
///
/// # Errors
/// The card's own refusals (parameters, frequency, thickness).
pub fn surface_impedance(card: &PorousCard, backing: Backing, thickness_m: f64, hz: f64) -> Result<C64, String> {
    let air = PorousAmbient::air();
    let omega = 2.0 * PI * hz;
    let rigid = card
        .surface_impedance_rigid(air, omega, thickness_m)
        .map_err(|e| format!("pillow: {} at {hz} Hz: {e}", card.name))?;
    match backing {
        Backing::Rigid => Ok(rigid),
        Backing::Air => {
            let (zc, _) = card
                .characteristic(air, omega)
                .map_err(|e| format!("pillow: {} at {hz} Hz: {e}", card.name))?;
            Ok(backed(zc, rigid, C64::from_re(air.rho0 * air.c0)))
        }
    }
}

/// The layer's surface impedance over a backing impedance `zb`, from its
/// rigid-backed one: `Zs_rigid = Zc coth(jkd)`, so `t = tanh(jkd) = Zc / Zs_rigid`
/// and `Zs = Zc (zb + Zc t) / (Zc + zb t)`, which tends to `Zs_rigid` as `zb`
/// grows without bound.
fn backed(zc: C64, rigid: C64, zb: C64) -> C64 {
    let t = zc / rigid;
    zc * (zb + zc * t) / (zc + zb * t)
}

fn log_point(lo: f64, hi: f64, i: usize, n: usize) -> f64 {
    if i == 0 {
        lo
    } else if i + 1 == n {
        hi
    } else {
        det::exp(det::ln(lo) + (det::ln(hi) - det::ln(lo)) * i as f64 / (n - 1) as f64).clamp(lo, hi)
    }
}

/// A pillow's fitted law.
#[derive(Debug, Clone, PartialEq)]
pub struct PillowLaw {
    pub pillow: Pillow,
    pub card_name: &'static str,
    pub provenance: &'static str,
    /// The card's own band [Hz].
    pub card_band_hz: (f64, f64),
    /// Dashpot per unit area [Pa s/m]: the law's resistance at high frequency.
    pub r_inf: f64,
    /// `(c_j [Pa s/m], tau_j [s])`: each kept branch's resistance at low
    /// frequency and its relaxation time.
    pub terms: Vec<(f64, f64)>,
    /// sup |R / Re Zs - 1| on `VERIFY_POINTS` log-spaced points of the fit
    /// band, inside the card's band and outside it.
    pub residual_in: f64,
    pub residual_out: f64,
    /// The card's `Re Zs` at `QUOTED_HZ` [Pa s/m].
    pub quoted: [f64; 4],
    /// The card's `Im Zs` at `QUOTED_HZ` [Pa s/m], the reactance the law drops.
    pub quoted_im: [f64; 4],
    /// The ladder's slowest and fastest rungs [Hz]; a branch kept on either
    /// is a fit the ladder, not the card, limits there.
    pub ladder_hz: (f64, f64),
}

impl PillowLaw {
    /// The law's resistance per unit area [Pa s/m].
    pub fn resistance(&self, hz: f64) -> f64 {
        let omega = 2.0 * PI * hz;
        self.r_inf
            + self
                .terms
                .iter()
                .map(|&(c, tau)| {
                    let x = omega * tau;
                    c / (1.0 + x * x)
                })
                .sum::<f64>()
    }

    /// The law's own reactance per unit area [Pa s/m] (e^{+jwt}): its branches'
    /// springs, `-sum_j c_j w tau_j / (1 + (w tau_j)^2)`, never positive.
    pub fn reactance(&self, hz: f64) -> f64 {
        let omega = 2.0 * PI * hz;
        -self
            .terms
            .iter()
            .map(|&(c, tau)| {
                let x = omega * tau;
                c * x / (1.0 + x * x)
            })
            .sum::<f64>()
    }

    /// The kept branches that sit on the ladder's slowest or fastest rung [Hz].
    pub fn end_rungs_hz(&self) -> Vec<f64> {
        let (slow, fast) = self.ladder_hz;
        self.terms
            .iter()
            .map(|t| 1.0 / (2.0 * PI * t.1))
            .filter(|hz| (hz / slow - 1.0).abs() < 1e-9 || (hz / fast - 1.0).abs() < 1e-9)
            .collect()
    }

    /// The branches as the plate steps them: `(e_j = c_j / tau_j [Pa/m], tau_j [s])`.
    pub fn branches(&self) -> Vec<(f64, f64)> {
        self.terms.iter().map(|&(c, tau)| (c / tau, tau)).collect()
    }

    /// The receipt lines.
    pub fn receipt(&self) -> Vec<String> {
        let p = self.pillow;
        let (lo, hi) = FIT_BAND_HZ;
        let (clo, chi) = self.card_band_hz;
        let ladder = if self.terms.is_empty() {
            "no relaxation branch".to_string()
        } else {
            let slow = self.terms.iter().map(|t| t.1).fold(0.0, f64::max);
            let fast = self.terms.iter().map(|t| t.1).fold(f64::INFINITY, f64::min);
            format!(
                "{} relaxation branches ({:.1} Hz to {:.0} Hz)",
                self.terms.len(),
                1.0 / (2.0 * PI * slow),
                1.0 / (2.0 * PI * fast)
            )
        };
        vec![
            format!(
                "pillow (opt-in, MB 22635): a {:.0} mm {} layer of {} on each head over a {:.0} mm radius centred at ({:.0}, {:.0}) mm, in place of the muffle dashpot; its normal surface resistance Re Zs from the card's JCAL law ({}; the row was identified on a {:.0} mm sample and is applied at the declared thickness): {:.0}, {:.0}, {:.0} and {:.0} Pa s/m at 45, 250, 1000 and 4000 Hz; the card's band is {clo}-{chi} Hz, and outside it the law is extrapolated",
                p.thickness_m * 1e3,
                p.backing.name(),
                self.card_name,
                p.patch.2 * 1e3,
                p.patch.0 * 1e3,
                p.patch.1 * 1e3,
                self.provenance,
                p.card.sample_thickness_m() * 1e3,
                self.quoted[0],
                self.quoted[1],
                self.quoted[2],
                self.quoted[3],
            ),
            format!(
                "pillow law: a dashpot of {:.1} Pa s/m plus {ladder}, fitted by non-negative least squares to the card's Re Zs alone over {lo}-{hi} Hz, never to a recording; sup |R / Re Zs - 1| {:.2} % inside the card's band and {:.2} % outside it, on {VERIFY_POINTS} log-spaced points; a relaxation ladder's resistance cannot rise with frequency",
                self.r_inf,
                100.0 * self.residual_in,
                100.0 * self.residual_out,
            ),
            format!(
                "pillow reactance: the card's Im Zs, dropped, is {:.0}, {:.0}, {:.0} and {:.0} Pa s/m at 45, 250, 1000 and 4000 Hz; the law's own, its branches' springs, is {:.0}, {:.0}, {:.0} and {:.0} Pa s/m there{}",
                self.quoted_im[0],
                self.quoted_im[1],
                self.quoted_im[2],
                self.quoted_im[3],
                self.reactance(QUOTED_HZ[0]),
                self.reactance(QUOTED_HZ[1]),
                self.reactance(QUOTED_HZ[2]),
                self.reactance(QUOTED_HZ[3]),
                match self.end_rungs_hz().as_slice() {
                    [] => String::new(),
                    ends => format!(
                        "; kept branches sit on the ladder's end rungs ({}), where the ladder, not the card, limits the fit",
                        ends.iter().map(|hz| format!("{hz:.1} Hz")).collect::<Vec<_>>().join(", ")
                    ),
                },
            ),
            "pillow, dropped and declared: the layer's own reactance (its added mass and stiffness) is dropped, the law's reactance being only the compliance its branches imply; JCAL takes the frame rigid, so the fill's own motion is not modelled; no pressure acts beyond the patch; that a pillow is there, its patch, thickness and backing are declared, and the fleece stands in for pillow fill, a borrowed material".to_string(),
        ]
    }
}

/// Fit a pillow's law (see the module docs).
///
/// # Errors
/// A non-finite or non-positive thickness or patch radius; the card's
/// refusals; a non-positive resistance on the card's curve; a fit that
/// leaves no resistance.
pub fn fit(pillow: Pillow) -> Result<PillowLaw, String> {
    let (_, _, radius) = pillow.patch;
    if !(pillow.thickness_m.is_finite() && pillow.thickness_m > 0.0) {
        return Err(format!("pillow: thickness {} m", pillow.thickness_m));
    }
    if !(radius.is_finite() && radius > 0.0 && pillow.patch.0.is_finite() && pillow.patch.1.is_finite()) {
        return Err(format!("pillow: patch {:?}", pillow.patch));
    }
    let card = pillow.card.card();
    card.validate().map_err(|e| format!("pillow: {}: {e}", card.name))?;
    let resistance = |hz: f64| -> Result<f64, String> {
        let r = surface_impedance(&card, pillow.backing, pillow.thickness_m, hz)?.re;
        if r.is_finite() && r > 0.0 {
            Ok(r)
        } else {
            Err(format!("pillow: the card's Re Zs is {r} Pa s/m at {hz} Hz"))
        }
    };
    let (lo, hi) = FIT_BAND_HZ;
    let ln10 = det::ln(10.0);
    let (l0, l1) = (det::ln(lo) - LADDER_MARGIN_DECADES * ln10, det::ln(hi) + LADDER_MARGIN_DECADES * ln10);
    let ladder = ((l1 - l0) / ln10 * BRANCHES_PER_DECADE).round() as usize + 1;
    let tau: Vec<f64> = (0..ladder)
        .map(|j| 1.0 / (2.0 * PI * det::exp(l0 + (l1 - l0) * j as f64 / (ladder - 1) as f64)))
        .collect();
    let cols = ladder + 1;
    let basis = |hz: f64, j: usize| -> f64 {
        if j == 0 {
            1.0
        } else {
            let x = 2.0 * PI * hz * tau[j - 1];
            1.0 / (1.0 + x * x)
        }
    };
    // Row i: R_fit(w_i) / R(w_i), against 1: the relative error.
    let mut a = vec![0.0; FIT_POINTS * cols];
    for i in 0..FIT_POINTS {
        let hz = log_point(lo, hi, i, FIT_POINTS);
        let r = resistance(hz)?;
        for j in 0..cols {
            a[i * cols + j] = basis(hz, j) / r;
        }
    }
    // nnls's entry threshold is absolute: give every column unit norm.
    let norms: Vec<f64> = (0..cols)
        .map(|j| det::sqrt((0..FIT_POINTS).map(|i| a[i * cols + j] * a[i * cols + j]).sum::<f64>()))
        .collect();
    if norms.iter().any(|n| !(n.is_finite() && *n > 0.0)) {
        return Err("pillow: a relaxation branch has no response in the fit band".into());
    }
    for i in 0..FIT_POINTS {
        for j in 0..cols {
            a[i * cols + j] /= norms[j];
        }
    }
    let x = nnls(&a, &[1.0; FIT_POINTS], FIT_POINTS, cols);
    let coef: Vec<f64> = x.iter().zip(&norms).map(|(x, n)| x / n).collect();
    let r_inf = coef[0];
    let terms: Vec<(f64, f64)> = coef[1..].iter().zip(&tau).filter(|(c, _)| **c > 0.0).map(|(c, t)| (*c, *t)).collect();
    if !(r_inf >= 0.0 && r_inf.is_finite()) || r_inf + terms.iter().map(|t| t.0).sum::<f64>() <= 0.0 {
        return Err("pillow: the fit leaves no resistance".into());
    }
    let mut quoted = [0.0; 4];
    let mut quoted_im = [0.0; 4];
    for ((q, qi), &hz) in quoted.iter_mut().zip(quoted_im.iter_mut()).zip(&QUOTED_HZ) {
        *q = resistance(hz)?;
        *qi = surface_impedance(&card, pillow.backing, pillow.thickness_m, hz)?.im;
    }
    let mut law = PillowLaw {
        pillow,
        card_name: card.name,
        provenance: card.provenance,
        card_band_hz: card.valid_hz,
        r_inf,
        terms,
        residual_in: 0.0,
        residual_out: 0.0,
        quoted,
        quoted_im,
        ladder_hz: (1.0 / (2.0 * PI * tau[0]), 1.0 / (2.0 * PI * tau[ladder - 1])),
    };
    let (clo, chi) = card.valid_hz;
    for i in 0..VERIFY_POINTS {
        let hz = log_point(lo, hi, i, VERIFY_POINTS);
        let miss = (law.resistance(hz) / resistance(hz)? - 1.0).abs();
        let slot = if hz >= clo && hz <= chi { &mut law.residual_in } else { &mut law.residual_out };
        *slot = slot.max(miss);
    }
    Ok(law)
}

#[derive(Clone, Copy, Debug)]
struct Branch {
    a: f64,
    b: f64,
    r: f64,
    /// `(1 + a)/(4 b)`: the coefficient of `|q|_A^2` in the discrete stored energy.
    stored: f64,
}

/// A pillow's branch state on one head: the patch nodes, their area, and
/// per branch the field `A q_j` on those nodes.
#[derive(Clone, Debug)]
pub(crate) struct PatchBranches {
    nodes: Vec<usize>,
    area: f64,
    moduli: Vec<f64>,
    taus: Vec<f64>,
    k: f64,
    coef: Vec<Branch>,
    p: Vec<Vec<f64>>,
    /// `q_j' A q_j` per branch.
    s: Vec<f64>,
    delta: Vec<f64>,
    pub(crate) heat_j: f64,
}

impl PatchBranches {
    pub(crate) fn new(nodes: Vec<usize>, area: f64, branches: &[(f64, f64)]) -> PatchBranches {
        let m = nodes.len();
        PatchBranches {
            nodes,
            area,
            moduli: branches.iter().map(|b| b.0).collect(),
            taus: branches.iter().map(|b| b.1).collect(),
            k: 0.0,
            coef: vec![],
            p: vec![vec![0.0; m]; branches.len()],
            s: vec![0.0; branches.len()],
            delta: vec![0.0; m],
            heat_j: 0.0,
        }
    }

    pub(crate) fn nodes(&self) -> &[usize] {
        &self.nodes
    }

    fn prepare(&mut self, k: f64) {
        if self.k == k && self.coef.len() == self.taus.len() {
            return;
        }
        self.coef = self
            .taus
            .iter()
            .map(|&tau| {
                let r = k / tau;
                let a = det::exp(-r);
                let b = -det::expm1(-r) / r;
                Branch { a, b, r, stored: (1.0 + a) / (4.0 * b) }
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
    }

    /// Advance the branches from `wp` to `w` at clock `k`, adding their forces to `f`.
    pub(crate) fn advance(&mut self, w: &[f64], wp: &[f64], f: &mut [f64], k: f64) {
        self.prepare(k);
        for (d, &i) in self.delta.iter_mut().zip(&self.nodes) {
            *d = w[i] - wp[i];
        }
        let area = self.area;
        // delta' A delta
        let dy = area * self.delta.iter().map(|d| d * d).sum::<f64>();
        let mut heat = 0.0;
        for (j, p) in self.p.iter_mut().enumerate() {
            let Branch { a, b, r, .. } = self.coef[j];
            let e = self.moduli[j];
            let pd: f64 = p.iter().zip(&self.delta).map(|(p, d)| p * d).sum();
            let s0 = self.s[j];
            // (r/4) |q^(n-1) + q^n|_A^2 with q^n = a q^(n-1) + b delta
            heat += e * 0.25 * r * ((1.0 + a) * (1.0 + a) * s0 + 2.0 * (1.0 + a) * b * pd + b * b * dy);
            self.s[j] = a * a * s0 + 2.0 * a * b * pd + b * b * dy;
            for ((p, d), &i) in p.iter_mut().zip(&self.delta).zip(&self.nodes) {
                *p = a * *p + b * area * d;
                f[i] -= e * *p;
            }
        }
        self.heat_j += heat;
    }

    /// The branches' discrete stored energy at the half step after the last
    /// one (`w`, `wp` the body's current pair) [J].
    pub(crate) fn stored_j(&self, w: &[f64], wp: &[f64]) -> f64 {
        if self.coef.is_empty() {
            return 0.0;
        }
        self.p
            .iter()
            .enumerate()
            .map(|(j, p)| {
                let pd: f64 = p.iter().zip(&self.nodes).map(|(p, &i)| p * (w[i] - wp[i])).sum();
                self.moduli[j] * (self.coef[j].stored * self.s[j] + 0.5 * pd)
            })
            .sum()
    }

    /// The branches' physical stored energy `sum_j e_j/2 |q_j|_A^2` [J].
    pub(crate) fn elastic_j(&self) -> f64 {
        self.moduli.iter().zip(&self.s).map(|(e, s)| 0.5 * e * s).sum()
    }

    /// Subtract the branches' glassy spring `sum_j e_j A` at each patch node
    /// times `w` from `f`: the stiffest they can be within a step (stability bound).
    pub(crate) fn subtract_glassy(&self, w: &[f64], f: &mut [f64]) {
        let glassy = self.area * self.moduli.iter().sum::<f64>();
        for &i in &self.nodes {
            f[i] -= glassy * w[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Body, Edge, Instrument, Nonlinear, PlateSpec};

    fn pillow(thickness_m: f64, backing: Backing) -> Pillow {
        Pillow { card: PillowCard::PintaFleeceFoam02, thickness_m, backing, patch: (0.0, -0.16, 0.0766) }
    }

    #[test]
    fn p0_air_backing_tends_to_the_rigid_layer_as_the_backing_stiffens() {
        let card = PillowCard::PintaFleeceFoam02.card();
        let omega = 2.0 * PI * 300.0;
        let rigid = card.surface_impedance_rigid(PorousAmbient::air(), omega, 0.1).unwrap();
        let (zc, _) = card.characteristic(PorousAmbient::air(), omega).unwrap();
        let far = backed(zc, rigid, C64::from_re(1.0e15));
        assert!((far - rigid).abs() < 1e-9 * rigid.abs(), "{far:?} against {rigid:?}");
        // A backing matched to the layer leaves its characteristic impedance.
        let matched = backed(zc, rigid, zc);
        assert!((matched - zc).abs() < 1e-12 * zc.abs());
    }

    #[test]
    fn p0_the_fit_meets_its_residual_off_its_own_grids() {
        for (d, backing, bound_in) in [
            (0.10, Backing::Rigid, 0.22),
            (0.10, Backing::Air, 0.05),
            (0.127, Backing::Rigid, 0.12),
            (0.127, Backing::Air, 0.04),
        ] {
            let law = fit(pillow(d, backing)).unwrap();
            eprintln!(
                "{d} m {backing:?}: Re Zs {:?}, Im Zs {:?} Pa s/m at {QUOTED_HZ:?} Hz; R_inf {:.2}, terms {:?}; residual {:.4} in, {:.4} out; end rungs {:?} Hz",
                law.quoted, law.quoted_im, law.r_inf, law.terms, law.residual_in, law.residual_out, law.end_rungs_hz()
            );
            for line in law.receipt() {
                eprintln!("  {line}");
            }
            assert!(QUOTED_HZ.iter().all(|&hz| law.reactance(hz) <= 0.0), "a Maxwell ladder is springs");
            assert!(law.receipt().iter().any(|l| l.contains("identified on a 50 mm sample")));
            assert!(law.r_inf >= 0.0 && law.terms.iter().all(|t| t.0 > 0.0 && t.1 > 0.0));
            assert!(law.residual_in <= bound_in, "{d} {backing:?}: {}", law.residual_in);
            let card = PillowCard::PintaFleeceFoam02.card();
            let (clo, chi) = card.valid_hz;
            // 97 points interleaved between both of the fit's grids.
            for i in 0..97 {
                let hz = det::exp(det::ln(41.0) + (det::ln(4950.0) - det::ln(41.0)) * (i as f64 + 0.37) / 97.0);
                let r = surface_impedance(&card, backing, d, hz).unwrap().re;
                let miss = (law.resistance(hz) / r - 1.0).abs();
                let stated = if hz >= clo && hz <= chi { law.residual_in } else { law.residual_out };
                assert!(miss <= stated * 1.02 + 1e-9, "{d} {backing:?} {hz} Hz: {miss} > {stated}");
            }
            assert!(law.quoted.iter().all(|q| q.is_finite() && *q > 0.0));
        }
    }

    #[test]
    fn p0_the_branches_dissipate_the_law_s_resistance() {
        // One node of unit area driven at w(t) = W sin(wt): the branches' force in phase with the
        // velocity, per unit velocity, is the law's resistance less its dashpot.
        let law = fit(pillow(0.10, Backing::Air)).unwrap();
        // Each period a whole number of steps at 192 kHz.
        for hz in [48.0, 250.0, 1000.0] {
            let k = 1.0 / 192_000.0;
            let omega = 2.0 * PI * hz;
            let mut field = PatchBranches::new(vec![0], 1.0, &law.branches());
            let period = (192_000.0 / hz).round() as usize;
            let (mut inphase, mut norm) = (0.0, 0.0);
            let w = |n: usize| 1e-4 * (omega * n as f64 * k).sin();
            for n in 1..=(40 * period) {
                let mut f = [0.0];
                field.advance(&[w(n)], &[w(n - 1)], &mut f, k);
                if n > 20 * period {
                    // the exact velocity at step n, where the recursion's state is
                    let v = 1e-4 * omega * (omega * n as f64 * k).cos();
                    inphase += -f[0] * v;
                    norm += v * v;
                }
            }
            let measured = inphase / norm;
            let expected = law.resistance(hz) - law.r_inf;
            assert!((measured / expected - 1.0).abs() < 2e-3, "{hz} Hz: {measured} against {expected}");
        }
    }

    fn head() -> Body {
        head_with(Body::plate_without_surface_loss)
    }

    fn head_with(build: fn(&str, &PlateSpec) -> Result<Body, String>) -> Body {
        let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
        build(
            "pillow head",
            &PlateSpec {
                radius: 0.2286,
                hole: 0.0,
                grid: 0.02,
                edge: Edge::Clamped,
                thickness: &thickness,
                rest_height: &rest,
                young: 4.9e9,
                poisson: 0.38,
                density: 1390.0,
                tension: 2000.0,
                bending: true,
                nonlinear: Nonlinear::Berger,
                t60: [(100.0, 1.0), (1000.0, 1.0)],
                support: None,
                muffle: None,
            },
        )
        .unwrap()
    }

    #[test]
    fn p0_the_branch_heat_closes_the_discrete_energy_and_the_dashpot_joins_s0() {
        let mut law = fit(pillow(0.10, Backing::Air)).unwrap();
        // The dashpot is s0's, which has no joule account: here the branches alone.
        let dashpot = law.r_inf;
        law.r_inf = 0.0;
        let mut body = head();
        let receipt = body.attach_pillow(&law).unwrap();
        assert!(receipt.contains("pillow branches on"), "{receipt}");
        assert!(body.attach_pillow(&law).is_err(), "a pillow attaches once");
        for (i, p) in body.xy.clone().iter().enumerate() {
            let (x, y) = (p[0] / 0.2286, p[1] / 0.2286);
            body.w[i] = 1e-4 * (1.0 - x * x - y * y).max(0.0) * (1.0 + 0.3 * y);
            body.wp[i] = 0.95 * body.w[i];
        }
        body.sync();
        let mut ins = Instrument::new(48_000.0);
        ins.add(body);
        assert!(ins.stability().0 < 1.0);
        let k = ins.k();
        let initial = ins.bodies[0].discrete_energy(k);
        let mut previous = initial;
        for _ in 0..4800 {
            ins.step(&mut []);
            let b = &ins.bodies[0];
            let energy = b.discrete_energy(k);
            assert!(energy <= previous + 1e-12 * initial, "the branches only remove energy");
            assert!((energy + b.pillow_heat_j() - initial).abs() < 1e-9 * initial, "{energy} + {} against {initial}", b.pillow_heat_j());
            previous = energy;
        }
        assert!(ins.bodies[0].pillow_heat_j() > 1e-2 * initial);
        // A clone steps exactly as the original (the branches are state).
        let mut copy = ins.clone();
        for _ in 0..64 {
            ins.step(&mut []);
            copy.step(&mut []);
        }
        assert_eq!(ins.bodies[0].w, copy.bodies[0].w);
        assert_eq!(ins.bodies[0].pillow_heat_j(), copy.bodies[0].pillow_heat_j());
        // The dashpot joins s0 on the patch nodes only, as a muffle of the same value does.
        law.r_inf = dashpot;
        let (mut with, plain) = (head(), head());
        let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
        let muffled = Body::plate_without_surface_loss(
            "muffled head",
            &PlateSpec {
                radius: 0.2286,
                hole: 0.0,
                grid: 0.02,
                edge: Edge::Clamped,
                thickness: &thickness,
                rest_height: &rest,
                young: 4.9e9,
                poisson: 0.38,
                density: 1390.0,
                tension: 2000.0,
                bending: true,
                nonlinear: Nonlinear::Berger,
                t60: [(100.0, 1.0), (1000.0, 1.0)],
                support: None,
                muffle: Some((law.pillow.patch.0, law.pillow.patch.1, law.pillow.patch.2, dashpot)),
            },
        )
        .unwrap();
        with.attach_pillow(&law).unwrap();
        assert_eq!(with.node_s0(), muffled.node_s0());
        assert_ne!(with.node_s0(), plain.node_s0());
    }

    #[test]
    fn p0_a_pillow_refuses_after_stepping_and_without_patch_nodes() {
        let law = fit(pillow(0.10, Backing::Rigid)).unwrap();
        let mut ins = Instrument::new(48_000.0);
        ins.add(head());
        ins.step(&mut []);
        assert!(ins.bodies[0].attach_pillow(&law).unwrap_err().contains("before the body steps"));
        let mut off = law.clone();
        off.pillow.patch = (5.0, 5.0, 0.01);
        assert!(head().attach_pillow(&off).unwrap_err().contains("covers no node"));
        assert!(head_with(Body::plate).attach_pillow(&law).unwrap_err().contains("fitted T60"));
        let mut wire = Body::string("wire", 8, 0.01, 1e-3, 10.0, 1.0, 0.0, 0.0);
        assert!(wire.attach_pillow(&law).is_err(), "a string carries no pillow");
        assert!(fit(Pillow { thickness_m: 0.0, ..pillow(0.1, Backing::Air) }).is_err());
        assert!(fit(Pillow { patch: (0.0, 0.0, -1.0), ..pillow(0.1, Backing::Air) }).is_err());
    }
}
