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
//! The scalar loss evaluates the complete diffusion response with a bounded
//! Bessel continued fraction and scaled high-frequency asymptotics; it does
//! not truncate this spectrum. The explicit [`circular_thermal_modes`] table
//! remains available with [`CIRCULAR_TED_MODES`] roots (omitted weight below
//! 2e-6). Its `J1'(x)` integral uses exponentially convergent periodic
//! trapezoid quadrature and root bisection. A small omitted total weight does
//! not bound relative high-frequency loss, so that table is not the scalar
//! loss oracle. The constitutive law retains weak thermoelastic coupling.
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

/// Complete insulated-cylinder diffusion loss per unit relaxation strength.
///
/// With `u = omega * base_time`, this evaluates the infinite thermal-mode
/// sum rather than truncating [`circular_thermal_modes`]. It is odd in signed
/// frequency, zero at DC, and returns NaN for nonfinite frequency or a
/// nonpositive/nonfinite thermal time. The physical loss is nonnegative at
/// positive frequency. This remains first-order weak thermoelastic coupling.
///
/// Norris and Photiadis, Q. J. Mech. Appl. Math. 58 (2005), Appendix C:
/// `F = 1 - 4 J2(z)/(z^2 J1'(z))`, `z^2 = -i u`, and loss = Im(F).
/// A positive Bessel continued fraction avoids cancellation at low frequency;
/// scaled modified-Bessel asymptotics avoid exponentially large intermediates
/// at high frequency. The scalar path neither allocates nor initializes roots.
#[must_use]
pub fn circular_relaxation_shape(omega: f64, base_time: f64) -> f64 {
    if !(omega.is_finite() && base_time.is_finite() && base_time > 0.0) {
        return f64::NAN;
    }
    if omega == 0.0 {
        return omega;
    }
    let frequency = omega.abs();
    let u = frequency * base_time;
    let shape = if u <= 1.0e-8 {
        // Exact first spectral moment sum(w_n / xi_n^2) = 7/24.
        // The omitted cubic term is below roundoff in this branch.
        (7.0 / 24.0) * u
    } else if u <= 1024.0 {
        circular_fraction_shape(u)
    } else if u.is_infinite() {
        // Preserve representable high-frequency loss even if omega * tau
        // overflows. Here both positive factors are necessarily greater than 1.
        (4.0 / frequency) / base_time
    } else {
        circular_asymptotic_shape(u)
    };
    shape.copysign(omega)
}

/// D = 4 + i*u/(6 + i*u/(8 + ...)); the response is Im(-4/(D+i*u)).
/// Backward evaluation uses positive real denominators and no subtraction of
/// nearly equal storage terms. At most 48 levels are needed on this interval.
fn circular_fraction_shape(u: f64) -> f64 {
    let depth = 16 + det::sqrt(u).ceil() as usize;
    let mut re = (2 * depth) as f64;
    let mut im = 0.0;
    for order in (2..depth).rev() {
        let scale = u / (re * re + im * im);
        (re, im) = (2.0 * order as f64 + scale * im, scale * re);
    }
    im += u;
    (4.0 / (re * re + im * im)) * im
}

/// For w^2 = i*u, Im(F) = (4/u) * (1 - Re(I1(w)/(w*I1'(w)))).
/// The common exp(w)/sqrt(2*pi*w) scale cancels. Twelve inverse-w terms at
/// |w| >= 32 resolve the quotient to binary64 accuracy; at higher frequencies
/// the expansion gets more accurate and the exponentially small term recedes.
fn circular_asymptotic_shape(u: f64) -> f64 {
    // 1/w = scale * (1-i). Keep the two scalar components explicit.
    let scale = 1.0 / (core::f64::consts::SQRT_2 * det::sqrt(u));
    let series = |order_squared: f64| {
        let (mut re, mut im) = (1.0, 0.0);
        let (mut term_re, mut term_im) = (1.0, 0.0);
        for k in 1..=12 {
            let odd = f64::from(2 * k - 1);
            let coefficient = (odd * odd - 4.0 * order_squared) / (8.0 * f64::from(k));
            let next_re = scale * (term_re + term_im) * coefficient;
            term_im = scale * (term_im - term_re) * coefficient;
            term_re = next_re;
            re += term_re;
            im += term_im;
        }
        (re, im)
    };
    let (i0_re, i0_im) = series(0.0);
    let (i1_re, i1_im) = series(1.0);
    let numerator_re = scale * (i1_re + i1_im);
    let numerator_im = scale * (i1_im - i1_re);
    // I1'(w) = I0(w) - I1(w)/w, with the common scale suppressed.
    let denominator_re = i0_re - numerator_re;
    let denominator_im = i0_im - numerator_im;
    let quotient_re = (numerator_re * denominator_re + numerator_im * denominator_im)
        / (denominator_re * denominator_re + denominator_im * denominator_im);
    (4.0 / u) * (1.0 - quotient_re)
}

impl ThermoelasticZener {
    /// Thermal time scale `rho c_p R^2 / k` of a solid circular section of
    /// radius `radius`; mode n relaxes in this time divided by `xi_n^2`.
    #[must_use]
    pub fn circular_base_time(&self, radius: f64) -> f64 {
        self.rho * self.cp * radius * radius / self.conductivity
    }

    /// Bending loss tangent of a solid circular section at angular frequency
    /// `omega`, including the complete insulated-cylinder thermal spectrum.
    #[must_use]
    pub fn circular_loss_factor(&self, omega: f64, radius: f64) -> f64 {
        self.relaxation_strength() * circular_relaxation_shape(omega, self.circular_base_time(radius))
    }
}

#[cfg(test)]
mod tests {
    use super::{circular_relaxation_shape, circular_thermal_modes};

    fn relative_near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual / expected - 1.0).abs() <= tolerance,
            "{actual:e} != {expected:e}, relative tolerance {tolerance:e}"
        );
    }

    #[test]
    fn complete_diffusion_matches_independent_bessel_reference() {
        // Independently evaluated with 90-decimal-digit defining Taylor sums
        // for J0 and J2, J1'=(J0-J2)/2, z^2=-i*u. No production continued
        // fraction, asymptotic coefficients, or truncated thermal sum is used.
        let cases = [
            (1.0e-6, 2.916_666_666_666_413e-7),
            (0.01, 2.916_641_319_665_006_3e-3),
            (0.1, 2.914_134_148_164_576_8e-2),
            (1.0, 2.683_485_156_018_991_7e-1),
            (3.39, 4.950_063_002_095_472_3e-1),
            (10.0, 3.036_931_813_990_054e-1),
            (100.0, 3.716_528_311_509_312e-2),
            (1024.0, 3.819_920_311_619_411_6e-3),
            (1024.000_000_01, 3.819_920_311_582_529_4e-3),
            (4096.0, 9.657_725_405_216_063e-4),
            (10_000.0, 3.971_715_349_931_839e-4),
        ];
        for (u, expected) in cases {
            relative_near(circular_relaxation_shape(u, 1.0), expected, 3.0e-14);
            relative_near(circular_relaxation_shape(-u, 1.0), -expected, 3.0e-14);
            relative_near(circular_relaxation_shape(u / 8.0, 8.0), expected, 3.0e-14);
        }
    }

    #[test]
    fn complete_spectrum_retains_both_moments_and_high_frequency_tail() {
        let low = 1.0e-10;
        relative_near(
            circular_relaxation_shape(low, 1.0),
            (7.0 / 24.0) * low,
            2.0e-15,
        );
        let high = 1.0e20;
        let complete = circular_relaxation_shape(high, 1.0);
        relative_near(complete, 4.0 / high, 2.0e-10);
        let truncated: f64 = circular_thermal_modes()
            .iter()
            .map(|&(xi, weight)| {
                let wt = high / (xi * xi);
                weight * wt / (1.0 + wt * wt)
            })
            .sum();
        // Tiny omitted static weight still causes over half a percent loss
        // error in the tail: the high-frequency moment weights each xi^2.
        assert!((complete - truncated) / complete > 0.005);
    }

    #[test]
    fn finite_extreme_scales_preserve_weak_loss_and_invalid_inputs_refuse() {
        relative_near(
            circular_relaxation_shape(1.0e200, 1.0),
            4.0e-200,
            3.0e-14,
        );
        relative_near(
            circular_relaxation_shape(1.0e200, 1.0e110),
            4.0e-310,
            3.0e-14,
        );
        relative_near(
            circular_relaxation_shape(1.0e-300, 1.0),
            (7.0 / 24.0) * 1.0e-300,
            3.0e-14,
        );
        assert_eq!(circular_relaxation_shape(0.0, 1.0), 0.0);
        for (omega, time) in [
            (f64::NAN, 1.0),
            (f64::INFINITY, 1.0),
            (1.0, f64::NAN),
            (1.0, f64::INFINITY),
            (1.0, 0.0),
            (1.0, -1.0),
        ] {
            assert!(circular_relaxation_shape(omega, time).is_nan());
        }
    }
}
