//! Geometric extension of retained physical strings (`Config::string_stretching`).
//!
//! The same per-key file and physics as the reference's `--string-stretching`
//! (`grand_piano/string_stretching.rs`, `longitudinal_string.rs`): fs-nlmodal's
//! `LongitudinalString` owns the moderate-slope strain `e = u_x + y_x^2/2`, its
//! storage and its gradient. A `stretch` row retains no axial coordinate (the
//! averaged Kirchhoff--Carrier tension); a `waves` row retains the axial sine
//! modes of its band, their damping, and an axial bridge port whose displacement
//! comes from the reference's geometry-derived axial participation row. Slopes
//! come from the segment's own coordinates and the
//! bridge, `y = sum phi_k (z_k - beta_k b) + x/L b`, with the reference's slope
//! scales (the Kirchhoff--Carrier channel's diagonal).
//!
//! The reference resolves the extension force implicitly, a discrete-gradient
//! fixed point iterated with its contact solve every tick. Here the linear motion of strings
//! and board keeps its exact transitions and the extension potential acts as one
//! velocity kick per step, evaluated where the step left the coordinates. The
//! axial modes are linear oscillators forced by the strings; with the strings'
//! start and end of the step both known, they take the exact response to that
//! force varying linearly over the step. The shared `fs-modal` owner also
//! integrates their physical viscous heat along that same trajectory. Earlier
//! timestep comparisons used the former scalar bridge port; that convergence
//! measurement must be repeated with the geometry-derived rows. The two
//! integration methods do not promise bitwise parity. The method's numerical
//! work is booked in `Ledger::coupling_j`.

use crate::refmods::linear::string_stretching::{Material, WaveMaterial};
use fs_math::det;
use fs_modal::time::{LinearForceTransition, ModalAcousticState};
use fs_nlmodal::{KcStringParams, kirchhoff_carrier_string, longitudinal::LongitudinalString};

/// Exact transition of y'' + c y' + w2 y = 0 over `dt` for (y, y'), underdamped.
pub(crate) fn damped_transition(w2: f64, c: f64, dt: f64) -> [f64; 4] {
    let w = det::sqrt(w2);
    let zeta = c / (2.0 * w);
    let wd = w * det::sqrt(1.0 - zeta * zeta);
    let (s, co) = det::sin_cos(wd * dt);
    let e = det::exp(-0.5 * c * dt);
    let h = 0.5 * c / wd;
    [e * (co + h * s), e * s / wd, -e * w2 * s / wd, e * (co - h * s)]
}

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

/// The extension state of one physical string segment. Active and sympathetic
/// partials feed the same retained transverse prefix and the same axial history.
#[derive(Clone)]
pub(crate) struct Stretch {
    model: LongitudinalString,
    length: f64,
    ea: f64,
    rest_tension: f64,
    material: Material,
    maximum_strain: f64,
    /// Only the `waves` model has an axial bridge port.
    axial_port: bool,
    /// Slope coefficient per unit relative coordinate of each retained mode.
    scales: Vec<f64>,
    /// Axial coordinate per unit axial bridge displacement (static correction).
    betas: Vec<f64>,
    omega2: Vec<f64>,
    /// Shared exact linearly forced transitions, including independent viscous heat.
    flow: Vec<LinearForceTransition>,
    step_s: f64,
    u: Vec<f64>,
    w: Vec<f64>,
    /// Axial force of the strings and bridge where the last step left them.
    force: Vec<f64>,
    /// Rejected axial trials must not advance any accepted coordinate or force history.
    next_u: Vec<f64>,
    next_w: Vec<f64>,
    next_force: Vec<f64>,
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
    /// d(potential)/d(axial termination displacement) [N]. Its own modal row
    /// returns the force to the board; it is not a vertical-force multiplier.
    pub axial_reaction: f64,
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
    /// 0.45 rate guard, a mode cap that cannot cover the band). The file owner's
    /// complete damping domain [0, 1], including critical damping, is retained.
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
        let (mut damping, mut maximum_strain) = (0.0, f64::INFINITY);
        let axial_port = waves.is_some();
        if let Some(w) = waves {
            if w.bandwidth_hz > 0.45 * f64::from(rate) {
                return Err("longitudinal band exceeds the mechanics clock's 0.45*rate guard".into());
            }
            if !(0.0..=1.0).contains(&w.damping_ratio) {
                return Err("live longitudinal waves need an axial damping ratio in [0, 1]".into());
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
            damping = w.damping_ratio;
            maximum_strain = w.maximum_strain;
        }
        let model = LongitudinalString::new(params, n, 1, &indices).map_err(|e| e.to_string())?;
        let dt = 1.0 / f64::from(rate);
        let omega2: Vec<f64> = model.axial_omegas().iter().map(|w| w * w).collect();
        let flow = model.axial_omegas().iter().map(|&omega|
            LinearForceTransition::try_new(omega, damping, dt).map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        let betas = model.axial_modes().iter().zip(model.axial_slope_scales())
            .map(|(&mode, &scale)| 2.0 * if mode % 2 == 0 { -1.0 } else { 1.0 } / (length * scale)).collect();
        let axial = indices.len();
        Ok(Self {
            slopes: vec![0.0; n + 1], slope_gradient: vec![0.0; n + 1],
            relative: vec![0.0; axial], harmonics: vec![0.0; axial], strain: vec![0.0; axial], weights: vec![0.0; axial],
            u: vec![0.0; axial], w: vec![0.0; axial], force: vec![0.0; axial],
            next_u: vec![0.0; axial], next_w: vec![0.0; axial], next_force: vec![0.0; axial],
            model, length, ea: material.axial_rigidity_n, rest_tension: tension, material, maximum_strain,
            axial_port, scales, betas, omega2, flow, step_s: dt, potential: 0.0,
        })
    }

    /// Retained axial coordinates.
    pub fn axial_modes(&self) -> usize {
        self.u.len()
    }

    /// Number of physical transverse partials driving this extension state.
    pub fn transverse_modes(&self) -> usize {
        self.scales.len()
    }

    /// Cold projection onto a changed transverse prefix of the SAME string.
    /// The axial sine coordinates have not changed basis: retain their physical
    /// displacement and velocity exactly. The new prefix's current strain is the
    /// next FOH drive's starting value, and its extension potential is the split's
    /// new energy baseline. The caller accounts separately for projection energy.
    ///
    /// When the retained union is unchanged, the kernel keeps the original owner
    /// without calling this method, preserving its drive history bit for bit.
    /// This is initialization of a replacement trial; an error never changes the
    /// previous accepted owner.
    pub fn reproject_from(&mut self, previous: Option<&Self>, q: &[f64], beta: &[f64],
        b: f64, axial_bridge: f64) -> Result<(), String> {
        if q.len() != self.scales.len() || beta.len() != q.len()
            || !b.is_finite() || !axial_bridge.is_finite()
            || q.iter().chain(beta).any(|x| !x.is_finite()) {
            return Err("string extension projection needs its complete finite transverse prefix".into());
        }
        if let Some(old) = previous {
            if self.length != old.length || self.ea != old.ea
                || self.rest_tension != old.rest_tension || self.step_s != old.step_s
                || self.axial_port != old.axial_port || self.omega2 != old.omega2
                || self.betas != old.betas || self.model.axial_modes() != old.model.axial_modes() {
                return Err("string extension projection cannot change its physical axial coordinates".into());
            }
            self.u.copy_from_slice(&old.u);
            self.w.copy_from_slice(&old.w);
        }
        let axial_bridge = if self.axial_port { axial_bridge } else { 0.0 };
        self.fill_slopes(q, beta, b);
        let mean = self.model.geometric_strain(&self.slopes, &mut self.harmonics)
            .map_err(|e| e.to_string())?;
        self.fill_relative(axial_bridge);
        let total = self.storage(mean, axial_bridge).0;
        let axial_own: f64 = self.u.iter().zip(&self.omega2).map(|(u, w2)| 0.5 * w2 * u * u).sum();
        let potential = total - axial_own;
        let c = self.model.strain_energy_coefficient();
        for j in 0..self.force.len() {
            self.force[j] = self.omega2[j] * self.betas[j] * axial_bridge
                - c * self.model.axial_slope_scales()[j] * self.harmonics[j];
        }
        if !total.is_finite() || !potential.is_finite()
            || self.force.iter().chain(&self.u).chain(&self.w).any(|x| !x.is_finite()) {
            return Err("string extension projection is not representable".into());
        }
        self.potential = potential;
        Ok(())
    }

    fn fill_slopes(&mut self, q: &[f64], beta: &[f64], b: f64) {
        self.slopes[0] = b / self.length;
        for (((slot, &scale), &z), &bk) in self.slopes[1..].iter_mut().zip(&self.scales).zip(q).zip(beta) {
            *slot = scale * (z - bk * b);
        }
    }

    fn fill_relative(&mut self, axial_bridge: f64) {
        for ((slot, &u), &bj) in self.relative.iter_mut().zip(&self.u).zip(&self.betas) {
            *slot = u - bj * axial_bridge;
        }
    }

    /// Total extension potential (with the axial stiffness), mean strain (with the
    /// bridge port's) and nonuniform strain bound at the filled slopes and axial state,
    /// from the harmonics of `geometric_strain` (`mean` is their geometric mean).
    fn storage(&mut self, mean: f64, axial_bridge: f64) -> (f64, f64, f64) {
        let (base, spread) = self.model.storage_from_strain(mean, &self.harmonics, &self.relative, &mut self.strain)
            .unwrap_or((f64::NAN, f64::NAN));
        let eps = axial_bridge / self.length;
        (base + 0.5 * self.ea * self.length * eps * (2.0 * mean + eps), mean + eps, spread)
    }

    /// Stored extension and axial kinetic energy at coordinates `q` (this segment's
    /// modes), mode constants `beta`, transverse bridge displacement `b` and
    /// axial termination displacement `axial_bridge` [J]. Cold path.
    pub fn stored(&self, q: &[f64], beta: &[f64], b: f64, axial_bridge: f64) -> f64 {
        let axial_bridge = if self.axial_port { axial_bridge } else { 0.0 };
        let mut copy = self.clone();
        copy.fill_slopes(q, beta, b);
        copy.fill_relative(axial_bridge);
        let mean = copy.model.geometric_strain(&copy.slopes, &mut copy.harmonics).unwrap_or(f64::NAN);
        copy.storage(mean, axial_bridge).0 + 0.5 * self.w.iter().map(|x| x * x).sum::<f64>()
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
    /// On an axial-transition refusal, no axial coordinate, force history,
    /// extension potential or supplied transverse velocity is published.
    pub fn step(&mut self, q: &[f64], beta: &[f64], b: f64, axial_bridge: f64, dt: f64, v: &mut [f64]) -> Result<Kick, String> {
        if dt != self.step_s {
            return Err("longitudinal step must use its prepared mechanics clock".into());
        }
        let axial_bridge = if self.axial_port { axial_bridge } else { 0.0 };
        self.fill_slopes(q, beta, b);
        let mean = self.model.geometric_strain(&self.slopes, &mut self.harmonics).unwrap_or(f64::NAN);
        let c = self.model.strain_energy_coefficient();
        let (mut loss, mut axial_change) = (0.0, 0.0);
        for j in 0..self.u.len() {
            let w2 = self.omega2[j];
            let d = self.model.axial_slope_scales()[j];
            let g1 = w2 * self.betas[j] * axial_bridge - c * d * self.harmonics[j];
            let (u0, w0) = (self.u[j], self.w[j]);
            let frame = self.flow[j].advance_with_energy(ModalAcousticState {
                displacement_m_sqrt_kg: u0,
                velocity_m_sqrt_kg_per_s: w0,
            }, self.force[j], g1).map_err(|e| e.to_string())?;
            let (u1, w1) = (frame.state.displacement_m_sqrt_kg, frame.state.velocity_m_sqrt_kg_per_s);
            self.next_u[j] = u1;
            self.next_w[j] = w1;
            self.next_force[j] = g1;
            loss += frame.viscous_dissipation_j;
            axial_change += 0.5 * (w1 * w1 - w0 * w0 + w2 * (u1 * u1 - u0 * u0));
        }
        std::mem::swap(&mut self.u, &mut self.next_u);
        std::mem::swap(&mut self.w, &mut self.next_w);
        std::mem::swap(&mut self.force, &mut self.next_force);
        self.fill_relative(axial_bridge);
        let (total, mean_total, spread) = self.storage(mean, axial_bridge);
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
        let eps = axial_bridge / self.length;
        let sg = &mut self.slope_gradient;
        if self.model.geometric_strain_gradient(&self.slopes, 2.0 * c * (mean + eps), &self.weights, sg).is_err() {
            sg.fill(0.0);
        }
        let mut reaction = sg[0] / self.length;
        let mut axial_reaction = self.ea * mean_total;
        let mut kinetic = 0.0;
        for (((vk, &g), &scale), &bk) in v.iter_mut().zip(&sg[1..]).zip(&self.scales).zip(beta) {
            let force = g * scale;
            reaction -= bk * force;
            let old = *vk;
            *vk = old - dt * force;
            kinetic += 0.5 * (*vk * *vk - old * old);
        }
        for ((&weight, &d), &bj) in self.weights.iter().zip(self.model.axial_slope_scales()).zip(&self.betas) {
            axial_reaction -= bj * weight * d;
        }
        // Exact dynamics would trade these against each other; what is left is the method's.
        let work = kinetic + axial_change + loss + potential - self.potential;
        self.potential = potential;
        Ok(Kick { reaction, axial_reaction: if self.axial_port { axial_reaction } else { 0.0 }, work, loss, out_of_domain })
    }
}

#[cfg(test)]
mod axial_transition_tests {
    use super::*;

    fn fixture(fundamental: f64, damping: f64, rate: u32, axial_modes: usize) -> Stretch {
        let density = 0.006;
        Stretch::new(
            Material { axial_rigidity_n: density * (2.0 * fundamental).powi(2), maximum_slope: 0.2 },
            Some(WaveMaterial {
                bandwidth_hz: (axial_modes as f64 + 0.01) * fundamental,
                mode_cap: axial_modes,
                damping_ratio: damping,
                maximum_strain: 0.01,
            }),
            1.0, density, 500.0, 0, rate,
        ).unwrap()
    }

    // Independent closed form, including the coalesced critical poles.
    fn free_motion(omega: f64, zeta: f64, q0: f64, v0: f64, t: f64) -> (f64, f64) {
        let a = zeta * omega;
        let decay = (-a * t).exp();
        if zeta == 1.0 {
            let slope = v0 + omega * q0;
            (decay * (q0 + slope * t), decay * (v0 - omega * slope * t))
        } else {
            let wd = omega * (1.0 - zeta * zeta).sqrt();
            let (s, c) = (wd * t).sin_cos();
            (decay * (q0 * c + (v0 + a * q0) * s / wd),
             decay * (v0 * c - (a * v0 + omega * omega * q0) * s / wd))
        }
    }

    /// G1: high-frequency axial motion and physical heat, independently
    /// integrated from the continuous viscous power. Endpoint trapezoids do
    /// not pass this test, including at the small material-loss limit.
    #[test]
    fn free_axial_ringdown_and_heat_follow_the_continuous_oscillator() {
        let rate = 48_000;
        let dt = 1.0 / f64::from(rate);
        for zeta in [0.0, 1e-10, 0.01, 0.8, 1.0] {
            let mut model = fixture(20_000.0, zeta, rate, 1);
            let omega = model.model.axial_omegas()[0];
            let (q0, v0) = (2e-8, 0.01);
            model.u[0] = q0;
            model.w[0] = v0;
            let kick = model.step(&[], &[], 0.0, 0.0, dt, &mut []).unwrap();
            let (q, v) = free_motion(omega, zeta, q0, v0, dt);
            let scale = (v0 * v0 + (omega * q0).powi(2)).sqrt();
            assert!(((model.u[0] - q) * omega).abs() < 2e-12 * scale);
            assert!((model.w[0] - v).abs() < 2e-12 * scale);

            // Composite Simpson integrates only the independent analytic v(t).
            let intervals = 4096;
            let integral = (0..=intervals).map(|i| {
                let weight = if i == 0 || i == intervals { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
                let velocity = free_motion(omega, zeta, q0, v0, dt * i as f64 / intervals as f64).1;
                weight * 2.0 * zeta * omega * velocity * velocity
            }).sum::<f64>() * dt / (3.0 * intervals as f64);
            if zeta == 0.0 {
                assert_eq!(kick.loss, 0.0);
            } else {
                assert!(kick.loss > 0.0);
                assert!((kick.loss - integral).abs() < 3e-11 * integral,
                    "zeta={zeta}: heat {}, independent viscous integral {integral}", kick.loss);
            }
            let e0 = 0.5 * scale * scale;
            assert!(kick.work.abs() < 5e-13 * e0,
                "free axial motion must not acquire coupling work: {}", kick.work);
        }
    }

    /// G1: the actual axial bridge force builds the lossless resonant secular
    /// response. Only interpolation of the sine remains, with second-order
    /// convergence as the existing mechanics clock is refined.
    #[test]
    fn axial_bridge_resonance_converges_without_a_numerical_frequency_shift() {
        let mut errors = Vec::new();
        for rate in [32_000, 64_000, 128_000] {
            let mut model = fixture(2000.0, 0.0, rate, 1);
            let omega = model.model.axial_omegas()[0];
            let amplitude = 1e-6;
            let steps = 8 * rate / 2000;
            let dt = 1.0 / f64::from(rate);
            for i in 1..=steps {
                let b = amplitude * (omega * f64::from(i) * dt).sin();
                let kick = model.step(&[], &[], 0.0, b, dt, &mut []).unwrap();
                assert_eq!(kick.loss, 0.0);
                assert!(!kick.out_of_domain);
            }
            let t = f64::from(steps) * dt;
            let scale = model.betas[0] * amplitude * 0.5;
            let q = scale * ((omega * t).sin() - omega * t * (omega * t).cos());
            let v = scale * omega * omega * t * (omega * t).sin();
            let error = ((model.u[0] - q).powi(2) + ((model.w[0] - v) / omega).powi(2)).sqrt() / q.abs();
            errors.push(error);
        }
        assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9, "{errors:?}");
        assert!(errors[2] < 0.001, "{errors:?}");
    }

    /// G0/G4: a later axial refusal cannot consume earlier coordinates or the
    /// endpoint-force history; the next admitted step exactly matches a retry.
    #[test]
    fn axial_transition_refusal_preserves_state_and_force_history() {
        let mut model = fixture(2000.0, 0.01, 48_000, 2);
        let dt = 1.0 / 48_000.0;
        model.step(&[], &[], 0.0, 1e-6, dt, &mut []).unwrap();
        let mut reference = model.clone();
        let beta = model.betas[1];
        // Fault the second mode only, after the first has completed its trial.
        model.betas[1] = f64::MAX;
        assert!(model.step(&[], &[], 0.0, 1.0, dt, &mut []).is_err());
        assert_eq!(model.u, reference.u);
        assert_eq!(model.w, reference.w);
        assert_eq!(model.force, reference.force);
        assert_eq!(model.potential.to_bits(), reference.potential.to_bits());
        model.betas[1] = beta;
        assert!(model.step(&[], &[], 0.0, 2e-6, 2.0 * dt, &mut []).is_err());
        let got = model.step(&[], &[], 0.0, 2e-6, dt, &mut []).unwrap();
        let expected = reference.step(&[], &[], 0.0, 2e-6, dt, &mut []).unwrap();
        assert_eq!(model.u, reference.u);
        assert_eq!(model.w, reference.w);
        assert_eq!(model.force, reference.force);
        assert_eq!(got.loss.to_bits(), expected.loss.to_bits());
        assert_eq!(got.work.to_bits(), expected.work.to_bits());
    }

    /// G0/G3: a transverse reduction does not reset an axial resonator. Its
    /// cross term can INCREASE extension energy, so projection is not heat.
    #[test]
    fn changing_the_transverse_prefix_retains_axial_motion_and_rebases_the_drive() {
        let make = |n| Stretch::new(
            Material { axial_rigidity_n: 96_000.0, maximum_slope: 0.2 },
            Some(WaveMaterial { bandwidth_hz: 8_100.0, mode_cap: 4,
                damping_ratio: 0.01, maximum_strain: 0.01 }),
            1.0, 0.006, 500.0, n, 48_000,
        ).unwrap();
        let mut full = make(2);
        assert_eq!(full.axial_modes(), 4);
        let slope = 0.01_f64;
        let geometric = 0.25 * slope * slope;
        let d4 = full.model.axial_slope_scales()[3];
        full.u[3] = -3.0 * geometric / d4;
        full.w[0] = 0.02;
        let q = [0.0, slope / full.scales[1]];
        let before = full.stored(&q, &[0.0; 2], 0.0, 0.0);
        let c = 0.5 * full.ea * full.length;
        // The second transverse sine has a constant strain and a fourth
        // cosine strain. Its latter term partly cancels the axial compression.
        let expected_before = 3.0 * c * geometric * geometric + 0.5 * 0.02_f64.powi(2);
        assert!((before - expected_before).abs() < 2e-14 * expected_before);

        let mut reduced = make(1);
        reduced.reproject_from(Some(&full), &[0.0], &[0.0], 0.0, 0.0).unwrap();
        assert_eq!(reduced.u, full.u);
        assert_eq!(reduced.w, full.w);
        let after = reduced.stored(&[0.0], &[0.0], 0.0, 0.0);
        let expected_after = 4.5 * c * geometric * geometric + 0.5 * 0.02_f64.powi(2);
        assert!(after > before, "removing a cancelling strain component can add extension storage");
        assert!((after - expected_after).abs() < 2e-14 * expected_after);
        assert!(reduced.force.iter().all(|&f| f == 0.0));
        assert!(reduced.potential.abs() < 2e-14 * after);

        // All remaining axial motion is now a freely decaying oscillator.
        // An obsolete pre-projection force or potential would introduce work.
        let omega = reduced.model.axial_omegas()[3];
        let (want_q, want_v) = free_motion(omega, 0.01, reduced.u[3], reduced.w[3], reduced.step_s);
        let kick = reduced.step(&[0.0], &[0.0], 0.0, 0.0, reduced.step_s, &mut [0.0]).unwrap();
        assert!((reduced.u[3] - want_q).abs() < 2e-12 * full.u[3].abs());
        assert!((reduced.w[3] - want_v).abs() < 2e-12 * (omega * full.u[3]).abs());
        assert!(kick.work.abs() < 2e-12 * after);
        assert!(kick.loss > 0.0);

        let mut expanded = make(5);
        expanded.reproject_from(Some(&reduced), &[0.0; 5], &[0.0; 5], 0.0, 0.0).unwrap();
        assert_eq!(expanded.u, reduced.u);
        assert_eq!(expanded.w, reduced.w);
        let kept = reduced.stored(&[0.0], &[0.0], 0.0, 0.0);
        assert_eq!(expanded.stored(&[0.0; 5], &[0.0; 5], 0.0, 0.0).to_bits(), kept.to_bits());
    }
}
