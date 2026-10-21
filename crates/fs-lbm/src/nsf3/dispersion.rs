//! Exact linearized Navier-Stokes-Fourier modes of a uniform ideal gas (bead pexfp.134, M3.2-W3b test support, shipped so that
//! M3.3-W1 can reuse it).
//!
//! A perturbation `exp(i k x + s t)` of density, velocity and temperature about a gas at rest obeys
//! `d rho'/dt = -rho0 i k u'`, `d u'/dt = -(R/rho0) i k (T0 rho' + rho0 T') - k^2 nu_L u'` with
//! `nu_L = (4 mu / 3 + mu_b) / rho0`, `d T'/dt = -(gamma - 1) T0 i k u' - chi k^2 T'` with
//! `chi = (gamma - 1) kappa / (rho0 R)`; the two transverse velocity components decouple and decay at `mu k^2 / rho0`.
//! The longitudinal 3 x 3 system has the characteristic polynomial (real coefficients)
//! `s^3 + k^2 (nu_L + chi) s^2 + (c^2 k^2 + nu_L chi k^4) s + R T0 chi k^4`, `c^2 = gamma R T0`.
//! Its roots are the acoustic pair, the entropy mode (`s about -kappa k^2 / (rho0 c_p)`) and, with the two shear modes, the
//! five modes of the bead. Convention: `e^{s t}` growth rates; `-Re s` is a decay rate and `Im s / k` a phase speed.
//!
//! No tuning: every coefficient is a gas receipt or a declared input; `GasPoint` carries SI numbers only.

use super::{NsfError, NsfLaw, THETA0};

/// A complex number (only what the modes need).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

/// A uniform ideal gas at rest with its transport coefficients, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasPoint {
    /// Density [kg/m^3].
    pub rho0: f64,
    /// Temperature [K].
    pub t0: f64,
    /// Specific gas constant [J/(kg K)].
    pub gas_constant: f64,
    /// Heat-capacity ratio.
    pub gamma: f64,
    /// Shear viscosity [Pa s].
    pub mu: f64,
    /// Bulk viscosity [Pa s].
    pub mu_bulk: f64,
    /// Thermal conductivity [W/(m K)].
    pub kappa: f64,
}

impl GasPoint {
    /// Adiabatic sound speed `sqrt(gamma R T0)` [m/s].
    #[must_use]
    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.gas_constant * self.t0).sqrt()
    }
    /// Isobaric specific heat `gamma R / (gamma - 1)` [J/(kg K)].
    #[must_use]
    pub fn specific_heat_cp(&self) -> f64 {
        self.gamma * self.gas_constant / (self.gamma - 1.0)
    }
    /// The point of a running grid: reference density and temperature, the receipt transport coefficients at the reference
    /// temperature converted from lattice units, and the declared bulk-to-shear ratio.
    #[must_use]
    pub fn from_law(law: &NsfLaw) -> GasPoint {
        let unit = law.rho_ref_kg_m3 * law.dx_m * law.dx_m / law.dt_s;
        let mu = law.mu_lattice(THETA0) * unit;
        GasPoint {
            rho0: law.rho_ref_kg_m3,
            t0: law.temperature_ref_k,
            gas_constant: law.gas_constant,
            gamma: law.gamma,
            mu,
            mu_bulk: law.bulk_to_shear_ratio * mu,
            kappa: law.kappa_lattice(THETA0) * law.gas_constant * unit,
        }
    }
}

fn valid(p: &GasPoint, k: f64) -> bool {
    [p.rho0, p.t0, p.gas_constant, p.gamma - 1.0, k].iter().all(|v| v.is_finite() && *v > 0.0)
        && p.mu.is_finite() && p.mu >= 0.0 && p.mu_bulk.is_finite() && p.mu_bulk >= 0.0 && p.kappa.is_finite() && p.kappa >= 0.0
}

fn cubic(s: f64, b: [f64; 3]) -> f64 {
    ((s + b[0]) * s + b[1]) * s + b[2]
}

/// The five modes `[acoustic +, acoustic -, entropy, shear, shear]` of wavenumber `k` [1/m] as `e^{s t}` rates.
///
/// # Errors
/// [`NsfError::Input`] for a non-positive or non-finite state or wavenumber.
pub fn linearized_modes(p: &GasPoint, k: f64) -> Result<[C64; 5], NsfError> {
    if !valid(p, k) {
        return Err(NsfError::Input("linearized modes need a positive finite gas state, non-negative transport and wavenumber"));
    }
    let nu_l = (4.0 * p.mu / 3.0 + p.mu_bulk) / p.rho0;
    let chi = (p.gamma - 1.0) * p.kappa / (p.rho0 * p.gas_constant);
    let c2 = p.gamma * p.gas_constant * p.t0;
    let q2 = k * k;
    let b = [q2 * (nu_l + chi), c2 * q2 + nu_l * chi * q2 * q2, p.gas_constant * p.t0 * chi * q2 * q2];
    // a real root in (-bound, 0]: the polynomial is -inf at -bound and b2 >= 0 at 0 (all coefficients non-negative)
    let bound = 1.0 + b[0].max(b[1]).max(b[2]);
    let (mut lo, mut hi) = (-bound.max(2.0 * (c2 * q2).sqrt() + 2.0 * b[0]), 0.0);
    while cubic(lo, b) > 0.0 {
        lo *= 2.0;
    }
    // bisection toward the root nearest zero when several exist is not needed: take the one bisection finds, then polish
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if cubic(mid, b) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let mut r = 0.5 * (lo + hi);
    for _ in 0..8 {
        let d = (3.0 * r + 2.0 * b[0]) * r + b[1];
        if d != 0.0 {
            r -= cubic(r, b) / d;
        }
    }
    // deflate: s^2 + (b0 + r) s + (b1 + r (b0 + r)) = 0
    let (qa, qb) = (b[0] + r, b[1] + r * (b[0] + r));
    let disc = qa * qa - 4.0 * qb;
    let (m1, m2) = if disc < 0.0 {
        let w = 0.5 * (-disc).sqrt();
        (C64 { re: -0.5 * qa, im: w }, C64 { re: -0.5 * qa, im: -w })
    } else {
        // real pair (overdamped): stable quadratic formula
        let sq = disc.sqrt();
        let big = -0.5 * (qa + qa.signum() * sq);
        let (x1, x2) = (big, if big != 0.0 { qb / big } else { 0.0 });
        (C64 { re: x1, im: 0.0 }, C64 { re: x2, im: 0.0 })
    };
    let shear = C64 { re: -p.mu * q2 / p.rho0, im: 0.0 };
    // the entropy mode is the real root of smallest magnitude when the pair is complex; otherwise order by magnitude
    if m1.im != 0.0 {
        Ok([m1, m2, C64 { re: r, im: 0.0 }, shear, shear])
    } else {
        let mut real = [r, m1.re, m2.re];
        real.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
        Ok([C64 { re: real[1], im: 0.0 }, C64 { re: real[2], im: 0.0 }, C64 { re: real[0], im: 0.0 }, shear, shear])
    }
}

/// The five modes of a running grid's gas at wavenumber `k_m` [1/m].
///
/// # Errors
/// As [`linearized_modes`].
pub fn linearized_nsf_modes(law: &NsfLaw, k_m: f64) -> Result<[C64; 5], NsfError> {
    linearized_modes(&GasPoint::from_law(law), k_m)
}
