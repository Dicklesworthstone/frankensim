//! Coupled quasi-steady orifice loss at the first `N <= 3` open tone holes (card row `open_hole_orifice N`).
//!
//! [`RegisterHoleShuntLine`] throttles the flow of the FIRST open hole only and leaves every other hole linear. This
//! module puts the same sharp-edged separation law, `Delta p_i = K_i |U_i| U_i` with `K_i = rho / (2 C_vc^2 A_i^2)`,
//! on the first `N` open holes at once. The holes are coupled through the bore, so they are solved together.
//!
//! The linear network behind the ports is the TMM chain with every port's orifice terminals opened. Seen from
//! those terminals it is the admittance matrix `Y` of `fs_duct::exterior::input_impedance_holes_and_ports`
//! (`U = u_0 - Y Delta p`, `U_i` the flow INTO hole `i`, `Delta p_i` the drop across its orifice). With the numerical
//! WDF port references `R = diag(Z_0,i)` the port waves are `p = a + b`, `R U_net = b - a`, and the
//! scattering matrix is `S = (I + R Y)^-1 (I - R Y)`; the ordering `R Y` (not `Y R`) matters as soon as the
//! references differ. Per sample, with the direct taps `D` and the histories `h` of the `S_ij` kernels
//! (`b = D a + h`):
//!
//! ```text
//!   Delta p = (I + D) a + h,           U = u_0 + R^-1 (b - a),
//!   U = u_th - A Delta p,              A = R^-1 (I - D)(I + D)^-1,
//!   u_th = u_0 + 2 R^-1 (I + D)^-1 h,  Delta p_i = K_i |U_i| U_i.
//! ```
//!
//! `U + A Phi(U) = u_th` has one solution because `A` is positive in the `R`-weighted inner product whenever the
//! weighted direct matrix has norm below one (the realization enforces it). It is found by a damped Newton
//! iteration on the `N`-vector (`N <= 3`). `N = 1` reduces to the closed form of [`RegisterHoleShuntLine`]; the
//! voice keeps that struct for `N = 1` so a card without the row and one with `open_hole_orifice 1` use
//! the same scalar solve. Both realizations retain phase when bounding reflected-wave gain.
//!
//! The returned reed perturbation is `sum_i s12_i * Delta p_i`. Radiated hole sound still comes from the linear
//! network, as for the first-hole shunt.
//!
//! With `open_hole_interaction` on a continuous gas field, that network is the COUPLED one: `Y` is
//! `fs_duct::interaction::port_admittance` (the coupled chimney admittance seen from the junctions less its transfer through the whole
//! bore, with the matched inlet and every opening stamped), the same network the reed's base reflectance, the hole flows and `s12` come
//! from, so the nonlinear drops perturb the network they belong to. Without the row `Y` is the independent network of the chimneys. At a
//! frequency where the coupled network is singular (a half-wave run) the independent `Y` is kept there.

use super::{OutletError, RegisterHoleShuntLine, invalid};
use super::super::RegionalDuct;
use crate::driving_point::{DrivingPointError, reflectance_rolloff, realize_reflectance_ir_with_rolloff};
use fs_duct::exterior::MAX_ORIFICE_PORTS;
use fs_duct::{LossModel, Termination, WallModel};
use fs_fft::{C64 as FftC64, Fft};
use fs_math::{c64::C64, det};
use fs_vfit::scattering::FirScatteringPort;

const NP: usize = MAX_ORIFICE_PORTS;
const _: () = assert!(NP == 16, "the `1..=16` refusal text in `realize_orifice_ports_with` names the cap: update it with the cap");
type Mat = [[f64; NP]; NP];
type CMat = [[C64; NP]; NP];

/// The voice's hole shunt: the unchanged first-hole line, or the coupled orifice ports.
#[derive(Clone, Debug)]
pub enum RegisterHoleShunt {
    /// The first open hole alone (`N = 1`, the default).
    Single(RegisterHoleShuntLine),
    /// The first `2..=MAX_ORIFICE_PORTS` open holes together.
    Ports(RegisterHoleOrificePorts),
}

impl RegisterHoleShunt {
    /// Incoming reed pressure perturbation from the orifice drops accepted so far.
    #[must_use]
    pub fn s12_history(&self) -> f64 {
        match self {
            Self::Single(line) => line.s12_history(),
            Self::Ports(ports) => ports.s12_history(),
        }
    }

    /// Prime sample `-1`; see [`RegisterHoleShuntLine::prime`].
    ///
    /// # Errors
    /// Nonfinite input or discrete filter failure.
    pub fn prime(&mut self, p_plus_prev: f64) -> Result<(), DrivingPointError> {
        match self {
            Self::Single(line) => line.prime(p_plus_prev),
            Self::Ports(ports) => ports.prime(p_plus_prev),
        }
    }

    /// Advance one sample; returns the flow through the first hole.
    ///
    /// # Errors
    /// Nonfinite state or an orifice solve that does not converge.
    pub fn step(&mut self, p_plus: f64) -> Result<f64, DrivingPointError> {
        match self {
            Self::Single(line) => line.step(p_plus),
            Self::Ports(ports) => ports.step(p_plus).map(|u| u[0]),
        }
    }
}

/// Realized kernels of the coupled orifice ports (plain impulse responses, before they become filters).
#[derive(Clone, Debug)]
pub(crate) struct RealizedPorts {
    pub(crate) ports: usize,
    /// `s12[i]`: reed-side incident perturbation per unit drop at port `i`.
    pub(crate) s12: Vec<Vec<f64>>,
    /// `s21[i]`: flow into port `i` per unit incident reed wave.
    pub(crate) s21: Vec<Vec<f64>>,
    /// `s22[i * ports + j]`: reflected wave at port `i` per unit incident wave at port `j`.
    pub(crate) s22: Vec<Vec<f64>>,
    pub(crate) r_ref: Vec<f64>,
    pub(crate) k_quad: Vec<f64>,
    /// Factor the passivity treatment applied to the `S22` tails (`1.0` = none needed).
    pub(crate) tail_damp: f64,
}

/// Time-domain WDF multiport: `N` coupled quasi-steady orifices at the first `N` open holes.
#[derive(Clone, Debug)]
pub struct RegisterHoleOrificePorts {
    ports: usize,
    s12: Vec<FirScatteringPort>,
    s21: Vec<FirScatteringPort>,
    s22: Vec<FirScatteringPort>,
    r_ref: Vec<f64>,
    k_quad: Vec<f64>,
}

impl RegisterHoleOrificePorts {
    pub(crate) fn from_realized(kernels: &RealizedPorts, dt: f64) -> Result<Self, OutletError> {
        let make = |ir: &Vec<f64>| {
            FirScatteringPort::new(dt, ir.clone())
                .map_err(|_| OutletError::from(invalid("register-hole multiport scattering port failed")))
        };
        let s12 = kernels.s12.iter().map(make).collect::<Result<Vec<_>, _>>()?;
        let s21 = kernels.s21.iter().map(make).collect::<Result<Vec<_>, _>>()?;
        let s22 = kernels.s22.iter().map(make).collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            ports: kernels.ports,
            s12,
            s21,
            s22,
            r_ref: kernels.r_ref.clone(),
            k_quad: kernels.k_quad.clone(),
        })
    }

    /// Number of coupled holes.
    #[must_use]
    pub const fn ports(&self) -> usize {
        self.ports
    }

    /// Incoming reed pressure perturbation from the accepted orifice drops, summed over the ports.
    #[must_use]
    pub fn s12_history(&self) -> f64 {
        let mut total = 0.0;
        for line in &self.s12 {
            total += line.response().history;
        }
        total
    }

    /// Copy with every `K_i = 0` (the linear, vortex-free limit).
    #[must_use]
    pub fn without_vortex_loss(&self) -> Self {
        let mut copy = self.clone();
        for k in &mut copy.k_quad {
            *k = 0.0;
        }
        copy
    }

    /// Prime sample `-1` with the outgoing reed wave and zero port waves.
    ///
    /// # Errors
    /// Nonfinite input or discrete filter failure.
    pub fn prime(&mut self, p_plus_prev: f64) -> Result<(), DrivingPointError> {
        for line in &mut self.s21 {
            line.push(p_plus_prev).map_err(|_| invalid("register-hole multiport s21 priming failed"))?;
        }
        for line in &mut self.s22 {
            line.push(0.0).map_err(|_| invalid("register-hole multiport s22 priming failed"))?;
        }
        for line in &mut self.s12 {
            line.push(0.0).map_err(|_| invalid("register-hole multiport s12 priming failed"))?;
        }
        Ok(())
    }

    /// Solve the coupled memoryless orifices at this sample and advance every kernel with the accepted
    /// `p_plus`, port waves `a` and drops `Delta p`. Returns the flows into the ports.
    ///
    /// # Errors
    /// Nonfinite state, a direct matrix outside the passive domain, or a solve that does not converge.
    pub fn step(&mut self, p_plus: f64) -> Result<[f64; NP], DrivingPointError> {
        self.step_law(p_plus, None)
    }

    /// As [`Self::step`], with the orifices replaced by the linear resistors `Delta p_i = rho_i U_i`
    /// (a test seam: the frequency-domain network predicts that case exactly).
    #[cfg(test)]
    pub(crate) fn step_linear(&mut self, p_plus: f64, rho: &[f64]) -> Result<[f64; NP], DrivingPointError> {
        self.step_law(p_plus, Some(rho))
    }

    fn step_law(&mut self, p_plus: f64, linear: Option<&[f64]>) -> Result<[f64; NP], DrivingPointError> {
        let n = self.ports;
        let mut d = [[0.0_f64; NP]; NP];
        let mut h = [0.0_f64; NP];
        let mut u0 = [0.0_f64; NP];
        for i in 0..n {
            u0[i] = self.s21[i].response().history;
            for j in 0..n {
                let response = self.s22[i * n + j].response();
                d[i][j] = response.direct;
                h[i] += response.history;
            }
        }
        if d.iter().flatten().chain(h.iter()).chain(u0.iter()).any(|x| !x.is_finite()) {
            return Err(invalid("register-hole multiport left the finite domain"));
        }
        let mut i_plus_d = d;
        for i in 0..n {
            i_plus_d[i][i] += 1.0;
        }
        let Some(inverse) = invert_real(&i_plus_d, n) else {
            return Err(invalid("register-hole multiport I + D is singular"));
        };
        // u_th = u_0 + 2 R^-1 (I + D)^-1 h ;  A = R^-1 (2 (I + D)^-1 - I).
        let mut u_th = [0.0_f64; NP];
        let mut a_mat = [[0.0_f64; NP]; NP];
        for i in 0..n {
            let mut g = 0.0;
            for j in 0..n {
                g += inverse[i][j] * h[j];
                a_mat[i][j] = (2.0 * inverse[i][j] - if i == j { 1.0 } else { 0.0 }) / self.r_ref[i];
            }
            u_th[i] = u0[i] + 2.0 * g / self.r_ref[i];
        }
        let mut dp = [0.0_f64; NP];
        let flows = if let Some(rho) = linear {
            let mut jac = [[0.0_f64; NP]; NP];
            for i in 0..n {
                for j in 0..n {
                    jac[i][j] = a_mat[i][j] * rho[j] + if i == j { 1.0 } else { 0.0 };
                }
            }
            let Some(jac_inv) = invert_real(&jac, n) else {
                return Err(invalid("register-hole multiport linear solve is singular"));
            };
            let mut flows = [0.0_f64; NP];
            for i in 0..n {
                for j in 0..n {
                    flows[i] += jac_inv[i][j] * u_th[j];
                }
                dp[i] = 0.0;
            }
            for i in 0..n {
                dp[i] = rho[i] * flows[i];
            }
            flows
        } else {
            let flows = solve_orifices(&a_mat, &self.k_quad, &u_th, n)?;
            for i in 0..n {
                dp[i] = self.k_quad[i] * flows[i].abs() * flows[i];
            }
            flows
        };
        // a = (I + D)^-1 (Delta p - h)
        let mut a_wave = [0.0_f64; NP];
        for i in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                sum += inverse[i][j] * (dp[j] - h[j]);
            }
            a_wave[i] = sum;
        }
        if a_wave.iter().chain(dp.iter()).any(|x| !x.is_finite()) {
            return Err(invalid("register-hole multiport reflected wave is not finite"));
        }
        for line in &mut self.s21 {
            line.push(p_plus).map_err(|_| invalid("register-hole multiport s21 step failed"))?;
        }
        for i in 0..n {
            for j in 0..n {
                self.s22[i * n + j]
                    .push(a_wave[j])
                    .map_err(|_| invalid("register-hole multiport s22 step failed"))?;
            }
        }
        for (line, &drop) in self.s12.iter_mut().zip(dp.iter()) {
            line.push(drop).map_err(|_| invalid("register-hole multiport s12 step failed"))?;
        }
        Ok(flows)
    }
}

/// Solve `U + A Phi(U) = u_th`, `Phi_i(U) = K_i |U_i| U_i`, by damped Newton from `U = u_th`.
///
/// Starting from `u_th` (not from the last sample) makes the result a pure function of this sample's inputs.
pub(crate) fn solve_orifices(a: &Mat, k: &[f64], u_th: &[f64; NP], n: usize) -> Result<[f64; NP], DrivingPointError> {
    let scale = u_th[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 || k[..n].iter().all(|&kk| kk == 0.0) {
        return Ok(*u_th);
    }
    let residual = |u: &[f64; NP]| -> [f64; NP] {
        let mut f = [0.0_f64; NP];
        for i in 0..n {
            let mut sum = u[i] - u_th[i];
            for j in 0..n {
                sum += a[i][j] * k[j] * u[j].abs() * u[j];
            }
            f[i] = sum;
        }
        f
    };
    let merit = |f: &[f64; NP]| f[..n].iter().map(|x| x * x).sum::<f64>();
    let mut u = *u_th;
    let mut f = residual(&u);
    for _ in 0..80 {
        let norm = f[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        if norm <= 1.0e-14 * scale {
            return Ok(u);
        }
        let mut jac = [[0.0_f64; NP]; NP];
        for i in 0..n {
            for j in 0..n {
                jac[i][j] = a[i][j] * 2.0 * k[j] * u[j].abs() + if i == j { 1.0 } else { 0.0 };
            }
        }
        let Some(jac_inv) = invert_real(&jac, n) else {
            return Err(invalid("register-hole multiport orifice Jacobian is singular"));
        };
        let mut step = [0.0_f64; NP];
        for i in 0..n {
            for j in 0..n {
                step[i] += jac_inv[i][j] * f[j];
            }
        }
        let m0 = merit(&f);
        let mut t = 1.0_f64;
        let mut accepted = false;
        for _ in 0..40 {
            let mut trial = u;
            for i in 0..n {
                trial[i] = u[i] - t * step[i];
            }
            let f_trial = residual(&trial);
            if f_trial[..n].iter().all(|x| x.is_finite()) && merit(&f_trial) <= (1.0 - 1.0e-4 * t) * m0 {
                u = trial;
                f = f_trial;
                accepted = true;
                break;
            }
            t *= 0.5;
        }
        if !accepted {
            // Roundoff floor: no step lowers the residual any further.
            let norm = f[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            if norm <= 1.0e-10 * scale {
                return Ok(u);
            }
            return Err(invalid("register-hole multiport orifice solve stalled"));
        }
    }
    let norm = f[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if norm <= 1.0e-10 * scale {
        Ok(u)
    } else {
        Err(invalid("register-hole multiport orifice solve did not converge"))
    }
}

fn invert_real(a: &Mat, n: usize) -> Option<Mat> {
    let mut m = *a;
    let mut inv = [[0.0_f64; NP]; NP];
    for (i, row) in inv.iter_mut().enumerate().take(n) {
        row[i] = 1.0;
    }
    for col in 0..n {
        let mut pivot = col;
        for r in col + 1..n {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        let p = m[pivot][col];
        if !p.is_finite() || p == 0.0 {
            return None;
        }
        m.swap(col, pivot);
        inv.swap(col, pivot);
        let scale = 1.0 / p;
        for c in 0..n {
            m[col][c] *= scale;
            inv[col][c] *= scale;
        }
        for r in 0..n {
            if r != col {
                let factor = m[r][col];
                if factor != 0.0 {
                    for c in 0..n {
                        m[r][c] -= factor * m[col][c];
                        inv[r][c] -= factor * inv[col][c];
                    }
                }
            }
        }
    }
    Some(inv)
}

fn invert_complex(a: &CMat, n: usize) -> Option<CMat> {
    let mut m = *a;
    let mut inv = [[C64::ZERO; NP]; NP];
    for (i, row) in inv.iter_mut().enumerate().take(n) {
        row[i] = C64::ONE;
    }
    for col in 0..n {
        let mut pivot = col;
        for r in col + 1..n {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        let p = m[pivot][col];
        let magnitude = p.abs();
        if !magnitude.is_finite() || magnitude == 0.0 {
            return None;
        }
        m.swap(col, pivot);
        inv.swap(col, pivot);
        let scale = p.recip();
        for c in 0..n {
            m[col][c] = m[col][c] * scale;
            inv[col][c] = inv[col][c] * scale;
        }
        for r in 0..n {
            if r != col {
                let factor = m[r][col];
                if factor.abs() != 0.0 {
                    for c in 0..n {
                        let sub_m = factor * m[col][c];
                        let sub_i = factor * inv[col][c];
                        m[r][c] = m[r][c] - sub_m;
                        inv[r][c] = inv[r][c] - sub_i;
                    }
                }
            }
        }
    }
    Some(inv)
}

/// `S = (I + R Y)^-1 (I - R Y) = 2 (I + R Y)^-1 - I` for the row-major port admittance `y`.
pub(crate) fn scattering_from_admittance(y: &[C64], r_ref: &[f64], n: usize) -> Option<CMat> {
    let mut a = [[C64::ZERO; NP]; NP];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = y[i * n + j].scale(r_ref[i]);
            if i == j {
                a[i][j] = a[i][j] + C64::ONE;
            }
        }
    }
    let inv = invert_complex(&a, n)?;
    let mut s = [[C64::ZERO; NP]; NP];
    for i in 0..n {
        for j in 0..n {
            s[i][j] = inv[i][j].scale(2.0);
            if i == j {
                s[i][j] = s[i][j] - C64::ONE;
            }
        }
    }
    Some(s)
}

/// Largest singular value of `T~ = R^-1/2 T R^1/2` (entries `T_ij sqrt(R_j / R_i)`): the power-wave norm, which
/// is the one that bounds a passive port set when the references differ.
pub(crate) fn weighted_spectral_norm(t: &CMat, r_ref: &[f64], n: usize) -> f64 {
    let mut w = [[C64::ZERO; NP]; NP];
    for i in 0..n {
        for j in 0..n {
            w[i][j] = t[i][j].scale(det::sqrt(r_ref[j] / r_ref[i]));
        }
    }
    if n == 1 {
        return w[0][0].abs();
    }
    // G = W^H W
    let mut g = [[C64::ZERO; NP]; NP];
    for j in 0..n {
        for k in 0..n {
            let mut re = 0.0;
            let mut im = 0.0;
            for i in 0..n {
                let c = C64::new(w[i][j].re, -w[i][j].im) * w[i][k];
                re += c.re;
                im += c.im;
            }
            g[j][k] = C64::new(re, im);
        }
    }
    let mut best = 0.0_f64;
    for start in 0..2 {
        let mut v = [C64::ZERO; NP];
        for (i, slot) in v.iter_mut().enumerate().take(n) {
            let x = i as f64;
            *slot = C64::from_re(if start == 0 { 1.0 / (1.0 + x) } else { 0.3 + 0.9 * x - 1.3 * (x * x) });
        }
        let mut lambda = 0.0_f64;
        for _ in 0..60 {
            let mut norm_sq = 0.0;
            for slot in v.iter().take(n) {
                norm_sq += slot.re * slot.re + slot.im * slot.im;
            }
            let inv_norm = 1.0 / det::sqrt(norm_sq);
            for slot in v.iter_mut().take(n) {
                *slot = slot.scale(inv_norm);
            }
            let mut gv = [C64::ZERO; NP];
            for j in 0..n {
                let mut re = 0.0;
                let mut im = 0.0;
                for k in 0..n {
                    let c = g[j][k] * v[k];
                    re += c.re;
                    im += c.im;
                }
                gv[j] = C64::new(re, im);
            }
            let mut gv_sq = 0.0;
            for slot in gv.iter().take(n) {
                gv_sq += slot.re * slot.re + slot.im * slot.im;
            }
            let next = det::sqrt(gv_sq);
            if next == 0.0 {
                lambda = 0.0;
                break;
            }
            lambda = next;
            v = gv;
        }
        best = best.max(lambda);
    }
    det::sqrt(best)
}

impl RegionalDuct {
    /// The shunt the voice uses: the unchanged first-hole line, or (for `open_hole_orifice N >= 2` with at
    /// least two open holes) the coupled orifice ports of the first `min(N, open holes)` holes.
    ///
    /// `single` is what [`Self::register_hole_shunt_line`] (or the baffled observer build) already produced.
    ///
    /// # Errors
    /// The multiport realization refuses (a segment-gas duct with several regions, vented hole without an
    /// orifice law, nonconvergent passivity treatment, grid budget).
    pub fn select_register_hole_shunt(
        &self,
        single: Option<RegisterHoleShuntLine>,
        sample_rate_hz: u32,
        fft_samples: usize,
        wall: WallModel<'_>,
    ) -> Result<Option<RegisterHoleShunt>, OutletError> {
        let wanted = self.open_hole_orifice_ports();
        let ports = wanted.min(self.open_hole_sites().len());
        if ports < 2 {
            return Ok(single.map(RegisterHoleShunt::Single));
        }
        let kernels = self.realize_orifice_ports(sample_rate_hz, fft_samples, wall, ports)?;
        let line = RegisterHoleOrificePorts::from_realized(&kernels, 1.0 / f64::from(sample_rate_hz))?;
        Ok(Some(RegisterHoleShunt::Ports(line)))
    }

    /// Frequency sweep, realization and passivity treatment of the `ports`-hole network (any `1..=MAX_ORIFICE_PORTS`).
    pub(crate) fn realize_orifice_ports(
        &self,
        sample_rate_hz: u32,
        n: usize,
        wall: WallModel<'_>,
        ports: usize,
    ) -> Result<RealizedPorts, OutletError> {
        self.realize_orifice_ports_with(sample_rate_hz, n, wall, ports, PassivityRule::Exact)
    }

    /// As [`Self::realize_orifice_ports`] with the passivity rule named (the tests compare the
    /// phase-aware rule with the former triangle bound).
    pub(crate) fn realize_orifice_ports_with(
        &self,
        sample_rate_hz: u32,
        n: usize,
        wall: WallModel<'_>,
        ports: usize,
        rule: PassivityRule,
    ) -> Result<RealizedPorts, OutletError> {
        let sites = self.open_hole_sites();
        let m = ports;
        if m == 0 || m > NP || m > sites.len() || sample_rate_hz == 0 || !(256..=super::MAX_OUTLET_FFT_SAMPLES).contains(&n)
            || !n.is_power_of_two()
        {
            return Err(invalid("orifice ports need 1..=16 open holes and a power-of-two grid in 256..=65536").into());
        }
        let zc = self.inlet_impedance();
        if !zc.is_finite() || zc <= 0.0 {
            return Err(invalid("orifice ports need a positive characteristic impedance").into());
        }
        let rate = f64::from(sample_rate_hz);
        let half = n / 2;
        let mut s12_pos = vec![vec![C64::ZERO; half + 1]; m];
        let mut s22_pos = vec![vec![C64::ZERO; half + 1]; m * m];
        let mut r_ref = vec![0.0_f64; m];
        let mut k_quad = vec![0.0_f64; m];
        for k in 1..=half {
            let theta = core::f64::consts::TAU * k as f64 / n as f64;
            let (response, flows, y) = self
                .response_with_wall_holes_and_ports(theta * rate, LossModel::Bessel, Termination::FlangedOpen, wall, m)
                .map_err(|e| OutletError::Transfer(DrivingPointError::Duct(e)))?;
            if flows.len() != sites.len() || y.len() != m * m {
                return Err(invalid("orifice-port transfer count mismatch across frequencies").into());
            }
            if k == 1 {
                for i in 0..m {
                    r_ref[i] = flows[i].shunt_port_reference_pa_s_m3;
                    k_quad[i] = flows[i].orifice_quadratic_coeff_pa_s2_m6;
                    if !(r_ref[i].is_finite() && r_ref[i] > 0.0) || !(k_quad[i].is_finite() && k_quad[i] >= 0.0) {
                        return Err(invalid("an open hole in the orifice ports has no valid orifice law").into());
                    }
                }
            }
            let denominator = C64::new(response.impedance.re + zc, response.impedance.im);
            let two_over_denom = denominator.recip().scale(2.0);
            let zc_over_denom = two_over_denom.scale(0.5 * zc);
            let Some(s) = scattering_from_admittance(&y, &r_ref, m) else {
                return Err(invalid("orifice-port scattering matrix is singular").into());
            };
            for i in 0..m {
                s12_pos[i][k] = flows[i].u_chimney_in_over_u_in * zc_over_denom;
                for j in 0..m {
                    s22_pos[i * m + j][k] = s[i][j];
                }
            }
        }
        // DC bins from a near-DC evaluation, as the first-hole line does.
        if let Ok((r0, flows0, y0)) = self.response_with_wall_holes_and_ports(
            1.0e-3, LossModel::Bessel, Termination::FlangedOpen, wall, m,
        ) {
            let denom0 = C64::new(r0.impedance.re + zc, r0.impedance.im);
            let zc_over_denom0 = denom0.recip().scale(zc);
            for i in 0..m {
                let s12_0 = (flows0[i].u_chimney_in_over_u_in * zc_over_denom0).re;
                if s12_0.is_finite() {
                    s12_pos[i][0] = C64::new(s12_0, 0.0);
                }
            }
            if let Some(s0) = scattering_from_admittance(&y0, &r_ref, m) {
                for i in 0..m {
                    for j in 0..m {
                        if s0[i][j].re.is_finite() {
                            s22_pos[i * m + j][0] = C64::new(s0[i][j].re.clamp(-1.0, 1.0), 0.0);
                        }
                    }
                }
            }
        }
        let rolloff = reflectance_rolloff(n);
        let mut s12 = Vec::with_capacity(m);
        for (i, pos) in s12_pos.iter().enumerate() {
            let mut ir = realize_reflectance_ir_with_rolloff(pos, &rolloff, n, false);
            let reg_delay = ((sites[i].bore_one_way_time_s * rate).floor().max(1.0) as usize).min(ir.len());
            for slot in &mut ir[..reg_delay] {
                *slot = 0.0;
            }
            s12.push(ir);
        }
        let scale_21 = 2.0 / zc;
        let s21: Vec<Vec<f64>> = s12.iter().map(|ir| ir.iter().map(|&h| h * scale_21).collect()).collect();
        let mut s22: Vec<Vec<f64>> = s22_pos
            .iter()
            .map(|pos| {
                let mut ir = realize_reflectance_ir_with_rolloff(pos, &rolloff, n, false);
                for slot in &mut ir[half + 1..] {
                    *slot = 0.0;
                }
                ir
            })
            .collect();
        let tail_damp = enforce_passive_matrix(&mut s22, &r_ref, m, n, rule);
        Ok(RealizedPorts { ports: m, s12, s21, s22, r_ref, k_quad, tail_damp })
    }
}

/// How the realized `S_ij` kernels are bounded to be passive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PassivityRule {
    /// Clamp the direct matrix to norm `0.999`, then scale the whole tail by the largest factor for which the
    /// realized matrix `D + lambda T(w)` has power-wave norm at most `1 - 1e-3` at every bin of the 2n grid.
    Exact,
    /// The former first-hole rule, retained as a comparison oracle: scale the tail so that
    /// `|D| + max_w |T(w)| <= 1 - 1e-3`. That sum ignores the phase between `D` and `T`, so it damps more than
    /// passivity needs even for a single port. It is not the production admission rule.
    Triangle,
}

/// Direct taps and tail of the `S_ij` kernels bounded in the power-wave norm by `rule`. Returns the tail factor
/// applied (`1.0` when none was needed).
fn enforce_passive_matrix(ir: &mut [Vec<f64>], r_ref: &[f64], m: usize, n: usize, rule: PassivityRule) -> f64 {
    if m == 1 {
        // A single port's reference cancels from its power-wave norm. Avoid
        // allocating one maximum-sized port matrix per frequency for this case.
        return enforce_passive_single(&mut ir[0], n, rule);
    }
    let half = n / 2;
    let mut d = [[C64::ZERO; NP]; NP];
    for i in 0..m {
        for j in 0..m {
            let direct = ir[i * m + j][0].clamp(-0.999, 0.999);
            ir[i * m + j][0] = direct;
            d[i][j] = C64::from_re(direct);
        }
    }
    let mut norm_d = weighted_spectral_norm(&d, r_ref, m);
    if norm_d > 0.999 {
        let scale = 0.999 / norm_d;
        for entry in ir.iter_mut() {
            entry[0] *= scale;
        }
        norm_d = 0.999;
    }
    let fft = Fft::new(2 * n);
    let mut scratch = vec![FftC64::new(0.0, 0.0); 2 * n];
    let mut spectra: Vec<Vec<FftC64>> = Vec::with_capacity(m * m);
    for entry in ir.iter() {
        let mut buf = vec![FftC64::new(0.0, 0.0); 2 * n];
        for (dst, &src) in buf[1..=half].iter_mut().zip(&entry[1..=half]) {
            *dst = FftC64::new(src, 0.0);
        }
        fft.forward(&mut buf, &mut scratch);
        spectra.push(buf);
    }
    // The tail's matrix at every bin of the oversampled grid, kept so the exact rule can re-evaluate the norm
    // of `D + lambda T` for several tail factors.
    let mut tails: Vec<CMat> = Vec::with_capacity(2 * n);
    for bin in 0..2 * n {
        let mut t = [[C64::ZERO; NP]; NP];
        for i in 0..m {
            for j in 0..m {
                let c = spectra[i * m + j][bin];
                t[i][j] = C64::new(c.re, c.im);
            }
        }
        tails.push(t);
    }
    let limit = 1.0 - 1.0e-3;
    let factor = match rule {
        PassivityRule::Triangle => {
            let max_tail = tails.iter().fold(0.0_f64, |top, t| top.max(weighted_spectral_norm(t, r_ref, m)));
            let max_allowed = limit - norm_d;
            if max_tail > max_allowed && max_tail > 0.0 { max_allowed / max_tail } else { 1.0 }
        }
        PassivityRule::Exact => {
            // The norm of the realized matrix itself, D + lambda T(w), is convex in lambda and equals |D| <= limit at
            // lambda = 0, so the passive factors form an interval [0, lambda_max]: bisect for its end.
            let mut direct = [[C64::ZERO; NP]; NP];
            for i in 0..m {
                for j in 0..m {
                    direct[i][j] = C64::from_re(ir[i * m + j][0]);
                }
            }
            let worst = |lambda: f64| -> f64 {
                tails.iter().fold(0.0_f64, |top, t| {
                    let mut s = [[C64::ZERO; NP]; NP];
                    for i in 0..m {
                        for j in 0..m {
                            s[i][j] = direct[i][j] + t[i][j].scale(lambda);
                        }
                    }
                    top.max(weighted_spectral_norm(&s, r_ref, m))
                })
            };
            passive_tail_factor(worst, limit)
        }
    };
    if factor < 1.0 {
        for entry in ir.iter_mut() {
            for slot in &mut entry[1..=half] {
                *slot *= factor;
            }
        }
    }
    factor
}

/// Scalar specialization of the same sampled reflected-wave norm admission.
/// The existing direct-tap margin, 2n frequency grid and tail-factor resolution
/// are unchanged. This is a sampled realization check, not a continuous-band
/// passivity certificate or a new physical loss law.
pub(super) fn enforce_passive_single(ir: &mut [f64], n: usize, rule: PassivityRule) -> f64 {
    let limit = 1.0 - 1.0e-3;
    let direct = ir[0].clamp(-limit, limit);
    ir[0] = direct;
    let fft = Fft::new(2 * n);
    let mut tail = vec![FftC64::new(0.0, 0.0); 2 * n];
    let mut scratch = vec![FftC64::new(0.0, 0.0); 2 * n];
    for (dst, &src) in tail[1..=n / 2].iter_mut().zip(&ir[1..=n / 2]) {
        *dst = FftC64::new(src, 0.0);
    }
    fft.forward(&mut tail, &mut scratch);
    let factor = match rule {
        PassivityRule::Triangle => {
            let largest = tail.iter().map(|t| det::hypot(t.re, t.im)).fold(0.0_f64, f64::max);
            let allowed = limit - direct.abs();
            if largest > allowed && largest > 0.0 { allowed / largest } else { 1.0 }
        }
        PassivityRule::Exact => passive_tail_factor(
            |lambda| {
                tail.iter()
                    .map(|t| det::hypot(direct + lambda * t.re, lambda * t.im))
                    .fold(0.0_f64, f64::max)
            },
            limit,
        ),
    };
    if factor < 1.0 {
        for tap in &mut ir[1..=n / 2] {
            *tap *= factor;
        }
    }
    factor
}

// The admitted direct operator already lies in the norm ball. Convexity of
// ||D + lambda T|| makes the feasible factors an interval starting at zero.
// Keep the original 2^-14 bisection resolution for scalar and matrix paths.
fn passive_tail_factor(worst: impl Fn(f64) -> f64, limit: f64) -> f64 {
    if worst(1.0) <= limit {
        return 1.0;
    }
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..14 {
        let mid = 0.5 * (lo + hi);
        if worst(mid) <= limit { lo = mid } else { hi = mid }
    }
    lo
}

#[cfg(test)]
mod tests;
