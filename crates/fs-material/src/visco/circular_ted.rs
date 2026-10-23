//! Zener thermoelastic damping of a solid circular section in bending.
//!
//! Bending strain `-kappa r sin(theta)` drives a temperature field through the
//! thermal expansion coupling. With an adiabatic surface, the field expands in
//! the disk's Neumann eigenfunctions `J1(xi_n r/R) sin(theta)`, where
//! `J1'(xi_n) = 0`. Each one relaxes independently with
//! `tau_n = rho c_p R^2 / (k xi_n^2)` (Zener 1938; for wires, Heptonstall et al.,
//! LIGO P1400229 eq. 6: f_c = 2.16 k/(rho C d^2) = 4 xi_1^2/(2 pi) k/(rho C d^2)).
//! Projecting the bending strain onto the modes gives the weights
//! `w_n = <r, J1>^2 / (<r, r> <J1, J1>) = 8 / (xi_n^2 (xi_n^2 - 1))`, using
//! `J2(xi_n) = J1(xi_n)/xi_n` at a zero of `J1'`. They sum to one, so
//!
//! `tan(delta)(omega) = Delta sum_n w_n omega tau_n / (1 + (omega tau_n)^2)`,
//! `Delta = E alpha^2 T / (rho c_p)`.
//!
//! Only `xi_n` is computed numerically. `J1'(x) = (1/2 pi) int_0^2pi sin(t)
//! sin(t - x sin(t)) dt` is periodic and analytic in t, so the trapezoid rule
//! converges exponentially; each zero is bisected in `[(n - 1/2) pi, (n + 1/2) pi]`.
//! The series keeps [`CIRCULAR_TED_MODES`] modes; the omitted weight is below 1e-6.
//! Homogeneous isotropic sections only: no layered heat flow, surface heat
//! exchange or anisotropic expansion.

use super::ThermoelasticZener;
use fs_math::det;
use std::sync::OnceLock;

/// Number of retained thermal modes of the circular-section series.
pub const CIRCULAR_TED_MODES: usize = 32;

const QUADRATURE_POINTS: usize = 512;

/// `J1'(x)` by the periodic trapezoid rule on its integral representation.
fn bessel_j1_derivative(x: f64) -> f64 {
    let step = core::f64::consts::TAU / QUADRATURE_POINTS as f64;
    let sum: f64 = (0..QUADRATURE_POINTS).map(|j| {
        let t = j as f64 * step;
        det::sin(t) * det::sin(t - x * det::sin(t))
    }).sum();
    sum / QUADRATURE_POINTS as f64
}

/// `(xi_n, w_n)` for n = 1..=CIRCULAR_TED_MODES: the zeros of `J1'` and the
/// bending-strain weights `8 / (xi_n^2 (xi_n^2 - 1))`, computed once.
#[must_use]
pub fn circular_thermal_modes() -> &'static [(f64, f64); CIRCULAR_TED_MODES] {
    static MODES: OnceLock<[(f64, f64); CIRCULAR_TED_MODES]> = OnceLock::new();
    MODES.get_or_init(|| {
        let mut modes = [(0.0, 0.0); CIRCULAR_TED_MODES];
        for (index, mode) in modes.iter_mut().enumerate() {
            let n = (index + 1) as f64;
            let pi = core::f64::consts::PI;
            let (mut lo, mut hi) = ((n - 0.5) * pi, (n + 0.5) * pi);
            let lo_sign = bessel_j1_derivative(lo) < 0.0;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if (bessel_j1_derivative(mid) < 0.0) == lo_sign { lo = mid; } else { hi = mid; }
            }
            let xi = 0.5 * (lo + hi);
            let xi2 = xi * xi;
            *mode = (xi, 8.0 / (xi2 * (xi2 - 1.0)));
        }
        modes
    })
}

/// `sum_n w_n omega tau_n / (1 + (omega tau_n)^2)` with `tau_n = base_time / xi_n^2`:
/// the loss tangent per unit relaxation strength. Zero for nonpositive,
/// nonfinite frequency or time.
#[must_use]
pub fn circular_relaxation_shape(omega: f64, base_time: f64) -> f64 {
    if !(omega.is_finite() && omega > 0.0 && base_time.is_finite() && base_time > 0.0) {
        return 0.0;
    }
    circular_thermal_modes().iter().map(|&(xi, weight)| {
        let wt = omega * base_time / (xi * xi);
        weight * wt / (1.0 + wt * wt)
    }).sum()
}

impl ThermoelasticZener {
    /// Thermal time scale `rho c_p R^2 / k` of a solid circular section of
    /// radius `radius`; mode n relaxes in this time divided by `xi_n^2`.
    #[must_use]
    pub fn circular_base_time(&self, radius: f64) -> f64 {
        self.rho * self.cp * radius * radius / self.conductivity
    }

    /// Bending loss tangent of a solid circular section at angular frequency
    /// `omega`, summed over the retained thermal modes.
    #[must_use]
    pub fn circular_loss_factor(&self, omega: f64, radius: f64) -> f64 {
        self.relaxation_strength() * circular_relaxation_shape(omega, self.circular_base_time(radius))
    }
}
