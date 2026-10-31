//! Geometric extension of the sounding strings (`Config::string_stretching`).
//!
//! The same per-key file and physics as the reference's `--string-stretching`
//! (`grand_piano/string_stretching.rs`, `longitudinal_string.rs`): fs-nlmodal's
//! `LongitudinalString` owns the moderate-slope strain `e = u_x + y_x^2/2`, its
//! storage and its gradient. A `stretch` row retains no axial coordinate (the
//! averaged Kirchhoff--Carrier tension); a `waves` row retains the axial sine
//! modes of its band, their damping, and an axial bridge port `u(L) = eta b` with
//! the authored `LIVE_BRIDGE_AXIAL_RATIO` (the reference now derives a per-mode row
//! instead; see that constant). Slopes come from the segment's own coordinates and the
//! bridge, `y = sum phi_k (z_k - beta_k b) + x/L b`, with the reference's slope
//! scales (the Kirchhoff--Carrier channel's diagonal).
//!
//! The reference resolves the extension force implicitly, a discrete-gradient
//! fixed point iterated with its contact solve every tick (about 80 us per 48 kHz
//! step for one A3 in `examples/validate.rs`). Here the linear motion of strings
//! and board keeps its exact transitions and the extension potential acts as one
//! velocity kick per step, evaluated where the step left the coordinates. The
//! axial modes are linear oscillators forced by the strings; with the strings'
//! start and end of the step both known, they take the exact response to that
//! force varying linearly over the step. Against the reference both converge to
//! the same motion at second order, and at 48 kHz the two differ by about as much
//! as the reference differs from itself at 192 kHz. The method's numerical work
//! is booked in `Ledger::coupling_j`.

use crate::refmods::linear::string_stretching::{Material, WaveMaterial};
use fs_math::det;
use fs_nlmodal::{KcStringParams, kirchhoff_carrier_string, longitudinal::LongitudinalString};

/// AUTHORED: the live engine's axial bridge coupling, `u(L) = LIVE_BRIDGE_AXIAL_RATIO * b` with `b` the bridge's
/// vertical displacement. It is the former single scalar of the offline renderer (`grand_piano`), whose comment
/// read ~1.8 degrees of downbearing plus a 28 mm bridge-cap lever arm and which never recorded a source. The
/// offline renderer no longer uses it: `grand_piano` derives a per-mode axial participation row from its board's
/// retained rotations (`BoardGeometry::bridge_attachments`, bead pexfp.44). This engine keeps the scalar, and
/// only here, until its board carries slope rows (`gen_board`'s `.bin`; bead frankensim-live-piano-axial-row-3eh0l);
/// until then the live and offline renders of a `waves` key couple the axial motion to the board differently.
pub const LIVE_BRIDGE_AXIAL_RATIO: f64 = 0.06;

/// Young's modulus [Pa] and density [kg/m^3] of piano wire (high-carbon steel), as the
/// reference's graded A3 file (`axial-a3.fspx`) uses them.
pub const STEEL_YOUNG_PA: f64 = 200.0e9;
pub const STEEL_DENSITY_KG_M3: f64 = 7850.0;

/// A complete `--string-stretching` file for `courses`: `waves` rows for plain steel
/// strings and `linear` rows for the rest. A course is plain steel when its own linear
/// density matches steel of the diameter its flexural rigidity implies
/// (`EI = E pi d^4/64`) within `tolerance`; then `EA = E pi d^2/4`. A wound string's
/// density is several times its core's, and its axial rigidity is not in the scale,
/// so it stays linear. The other columns are the graded A3 file's: slope bound 0.2,
/// the axial band, a mode cap of 8, axial damping ratio `zeta` (an assumption there
/// too) and a strain bound of 0.01.
pub fn plain_steel_waves_spec(courses: &[crate::refmods::geometry::Course], band_hz: f64, zeta: f64, tolerance: f64) -> String {
    use std::fmt::Write;
    let mut text = String::from("frankensim-piano-string-stretching-v1\n");
    for c in courses {
        let d4 = 64.0 * c.flexural_rigidity_nm2 / (std::f64::consts::PI * STEEL_YOUNG_PA);
        let area = 0.25 * std::f64::consts::PI * det::sqrt(d4);
        let implied = STEEL_DENSITY_KG_M3 * area;
        if (implied / c.linear_density_kg_m - 1.0).abs() <= tolerance {
            let _ = writeln!(text, "waves,{},{:e},0.2,{band_hz:e},8,{zeta:e},0.01", c.midi, STEEL_YOUNG_PA * area);
        } else {
            let _ = writeln!(text, "linear,{}", c.midi);
        }
    }
    text
}

/// The extension state of one sounding string segment.
#[derive(Clone)]
pub(crate) struct Stretch {
    model: LongitudinalString,
    length: f64,
    ea: f64,
    rest_tension: f64,
    material: Material,
    maximum_strain: f64,
    /// Axial bridge port (zero for `stretch`).
    port: f64,
    /// Slope coefficient per unit relative coordinate of each retained mode.
    scales: Vec<f64>,
    /// Axial coordinate per unit axial bridge displacement (static correction).
    betas: Vec<f64>,
    omega2: Vec<f64>,
    /// Axial damping coefficient 2 zeta w [1/s] and the exact free transition per step.
    damping: Vec<f64>,
    flow: Vec<[f64; 4]>,
    u: Vec<f64>,
    w: Vec<f64>,
    /// Axial force of the strings and bridge where the last step left them.
    force: Vec<f64>,
    /// Extension potential less the axial oscillators' own, where the last kick left it [J].
    potential: f64,
    slopes: Vec<f64>,
    slope_gradient: Vec<f64>,
    relative: Vec<f64>,
    harmonics: Vec<f64>,
    strain: Vec<f64>,
    weights: Vec<f64>,
}

/// What one step's kick did.
pub(crate) struct Kick {
    /// d(potential)/d(bridge displacement) [N]: the board takes minus this.
    pub reaction: f64,
    /// Numerical work of the method on this segment [J].
    pub work: f64,
    /// Axial modal loss of the step [J].
    pub loss: f64,
    /// The strain left the supplied slope or strain domain, or tension went nonpositive.
    pub out_of_domain: bool,
}

impl Stretch {
    /// `n` retained transverse modes of a segment of `length`, `density` and installed
    /// `tension`, at the mechanical `rate`. Refuses as the reference does (band above the
    /// 0.45 rate guard, a mode cap that cannot cover the band); a damping ratio of one or
    /// more is refused because the axial flow is the underdamped exact transition.
    pub fn new(material: Material, waves: Option<WaveMaterial>, length: f64, density: f64, tension: f64,
        n: usize, rate: u32) -> Result<Self, String> {
        let params = KcStringParams { length, tension, lin_density: density, ea: material.axial_rigidity_n };
        let count = n.max(1);
        let owner = kirchhoff_carrier_string(&params, count).map_err(|e| e.to_string())?;
        if owner.channels.len() != 1 || owner.channels[0].coupling.len() != count * count {
            return Err("string owner did not return its single complete stress channel".into());
        }
        let coupling = &owner.channels[0].coupling;
        let scales: Vec<f64> = (0..n).map(|i| det::sqrt(coupling[i * count + i])).collect();
        if scales.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return Err("derived string extension channel is unrepresentable".into());
        }
        let mut indices = Vec::new();
        let (mut port, mut damping, mut maximum_strain) = (0.0, 0.0, f64::INFINITY);
        if let Some(w) = waves {
            if w.bandwidth_hz > 0.45 * f64::from(rate) {
                return Err("longitudinal band exceeds the mechanics clock's 0.45*rate guard".into());
            }
            if !(0.0..1.0).contains(&w.damping_ratio) {
                return Err("live longitudinal waves need an axial damping ratio in [0, 1)".into());
            }
            let fundamental = 0.5 * det::sqrt(material.axial_rigidity_n / density) / length;
            for mode in 1..=w.mode_cap + 1 {
                if mode as f64 * fundamental > w.bandwidth_hz {
                    break;
                }
                if mode > w.mode_cap {
                    return Err("longitudinal mode cap cannot cover the supplied axial band".into());
                }
                indices.push(mode);
            }
            port = LIVE_BRIDGE_AXIAL_RATIO;
            damping = w.damping_ratio;
            maximum_strain = w.maximum_strain;
        }
        let model = LongitudinalString::new(params, n, 1, &indices).map_err(|e| e.to_string())?;
        let dt = 1.0 / f64::from(rate);
        let omega2: Vec<f64> = model.axial_omegas().iter().map(|w| w * w).collect();
        let coefficients: Vec<f64> = model.axial_omegas().iter().map(|&w| 2.0 * damping * w).collect();
        let flow = omega2.iter().zip(&coefficients).map(|(&w2, &c)| crate::kernel::damped_transition(w2, c, dt)).collect();
        let betas = model.axial_modes().iter().zip(model.axial_slope_scales())
            .map(|(&mode, &scale)| 2.0 * if mode % 2 == 0 { -1.0 } else { 1.0 } / (length * scale)).collect();
        let axial = indices.len();
        Ok(Self {
            slopes: vec![0.0; n + 1], slope_gradient: vec![0.0; n + 1],
            relative: vec![0.0; axial], harmonics: vec![0.0; axial], strain: vec![0.0; axial], weights: vec![0.0; axial],
            u: vec![0.0; axial], w: vec![0.0; axial], force: vec![0.0; axial],
            model, length, ea: material.axial_rigidity_n, rest_tension: tension, material, maximum_strain,
            port, scales, betas, omega2, damping: coefficients, flow, potential: 0.0,
        })
    }

    /// Retained axial coordinates.
    pub fn axial_modes(&self) -> usize {
        self.u.len()
    }

    fn fill_slopes(&mut self, q: &[f64], beta: &[f64], b: f64) {
        self.slopes[0] = b / self.length;
        for (((slot, &scale), &z), &bk) in self.slopes[1..].iter_mut().zip(&self.scales).zip(q).zip(beta) {
            *slot = scale * (z - bk * b);
        }
    }

    fn fill_relative(&mut self, b: f64) {
        let axial_bridge = self.port * b;
        for ((slot, &u), &bj) in self.relative.iter_mut().zip(&self.u).zip(&self.betas) {
            *slot = u - bj * axial_bridge;
        }
    }

    /// Total extension potential (with the axial stiffness), mean strain (with the
    /// bridge port's) and nonuniform strain bound at the filled slopes and axial state,
    /// from the harmonics of `geometric_strain` (`mean` is their geometric mean).
    fn storage(&mut self, mean: f64, b: f64) -> (f64, f64, f64) {
        let (base, spread) = self.model.storage_from_strain(mean, &self.harmonics, &self.relative, &mut self.strain)
            .unwrap_or((f64::NAN, f64::NAN));
        let eps = self.port * b / self.length;
        (base + 0.5 * self.ea * self.length * eps * (2.0 * mean + eps), mean + eps, spread)
    }

    /// Stored extension and axial kinetic energy at coordinates `q` (this segment's
    /// modes), mode constants `beta` and bridge displacement `b` [J]. Cold path.
    pub fn stored(&self, q: &[f64], beta: &[f64], b: f64) -> f64 {
        let mut copy = self.clone();
        copy.fill_slopes(q, beta, b);
        copy.fill_relative(b);
        let mean = copy.model.geometric_strain(&copy.slopes, &mut copy.harmonics).unwrap_or(f64::NAN);
        copy.storage(mean, b).0 + 0.5 * self.w.iter().map(|x| x * x).sum::<f64>()
    }

    /// One step of the axial modes, then the kick of the extension potential on the segment's
    /// velocities `v`, from the coordinates the step left (`q`, `b`). One pass over the slopes
    /// gives the strain harmonics; everything else follows from them.
    ///
    /// The axial modes are linear oscillators forced by the strings and the bridge:
    /// `g_j = -(dE/du_j - w_j^2 u_j)` does not depend on `u`. With the strings' start and end
    /// of the step both known, they take the exact response to that force varying linearly
    /// over the step (a first-order hold): an impulse would sample oscillators of up to
    /// `w dt` ~ 1.3 at 48 kHz too coarsely.
    pub fn step(&mut self, q: &[f64], beta: &[f64], b: f64, dt: f64, v: &mut [f64]) -> Kick {
        self.fill_slopes(q, beta, b);
        let mean = self.model.geometric_strain(&self.slopes, &mut self.harmonics).unwrap_or(f64::NAN);
        let c = self.model.strain_energy_coefficient();
        let (mut loss, mut axial_change) = (0.0, 0.0);
        let axial_bridge = self.port * b;
        for j in 0..self.u.len() {
            let (w2, damping, p) = (self.omega2[j], self.damping[j], self.flow[j]);
            let d = self.model.axial_slope_scales()[j];
            let g1 = w2 * self.betas[j] * axial_bridge - c * d * self.harmonics[j];
            let g0 = std::mem::replace(&mut self.force[j], g1);
            let (u0, w0) = (self.u[j], self.w[j]);
            // Particular response to g0 + (g1 - g0) t/dt: u_p = a + rate t.
            let rate = (g1 - g0) / (dt * w2);
            let a = (g0 - damping * rate) / w2;
            let (du, dw) = (u0 - a, w0 - rate);
            let (u1, w1) = (a + rate * dt + p[0] * du + p[1] * dw, rate + p[2] * du + p[3] * dw);
            self.u[j] = u1;
            self.w[j] = w1;
            loss += 0.5 * damping * dt * (w0 * w0 + w1 * w1);
            axial_change += 0.5 * (w1 * w1 - w0 * w0 + w2 * (u1 * u1 - u0 * u0));
        }
        self.fill_relative(b);
        let (total, mean_total, spread) = self.storage(mean, b);
        let axial_own: f64 = self.u.iter().zip(&self.omega2).map(|(u, w2)| 0.5 * w2 * u * u).sum();
        let potential = total - axial_own;
        let slope_bound = self.slopes[0].abs() + self.slopes[1..].iter().map(|a| a.abs()).sum::<f64>();
        let out_of_domain = !(slope_bound <= self.material.maximum_slope)
            || !(mean_total.abs() + spread <= self.maximum_strain)
            || self.rest_tension + self.ea * (mean_total - spread) <= 0.0;
        for (weight, &s) in self.weights.iter_mut().zip(&self.strain) {
            *weight = c * s;
        }
        // The bridge port's mean extension adds EA L eps to the mean's weight
        // (`longitudinal_string.rs`, `evaluate`): 2 c (s_0 + eps) in all.
        let eps = self.port * b / self.length;
        let sg = &mut self.slope_gradient;
        if self.model.geometric_strain_gradient(&self.slopes, 2.0 * c * (mean + eps), &self.weights, sg).is_err() {
            sg.fill(0.0);
        }
        let mut reaction = sg[0] / self.length + self.port * self.ea * mean_total;
        let mut kinetic = 0.0;
        for (((vk, &g), &scale), &bk) in v.iter_mut().zip(&sg[1..]).zip(&self.scales).zip(beta) {
            let force = g * scale;
            reaction -= bk * force;
            let old = *vk;
            *vk = old - dt * force;
            kinetic += 0.5 * (*vk * *vk - old * old);
        }
        for ((&weight, &d), &bj) in self.weights.iter().zip(self.model.axial_slope_scales()).zip(&self.betas) {
            reaction -= self.port * bj * weight * d;
        }
        // Exact dynamics would trade these against each other; what is left is the method's.
        let work = kinetic + axial_change + loss + potential - self.potential;
        self.potential = potential;
        Kick { reaction, work, loss, out_of_domain }
    }
}
