//! Acoustic root of homogeneous linear Navier--Stokes--Fourier propagation.
//!
//! With `exp(i k x - i omega t)`, `nu = (4 mu/3 + zeta)/rho`,
//! `chi = kappa/(rho Cp)`, `gamma = Cp/Cv`, `x = omega nu/c²`,
//! `y = omega chi/c²` and `z = (k c/omega)²`, elimination of density,
//! temperature and longitudinal velocity gives
//!
//! `(gamma*x*y + i*y) z² + [1 - i*(x + gamma*y)] z - 1 = 0`.
//!
//! This is the spatially attenuated, real-frequency dispersion problem, not
//! the temporally attenuated real-wavenumber problem. See Zaidni & Benjelloun,
//! arXiv:2107.08886, sections 2--3 and 5, <https://arxiv.org/abs/2107.08886>.
//! Their Laplace/spatial convention is converted here to the crate convention.
//! The acoustic root continues from z=1 as both transport coefficients vanish;
//! the other root is the thermal mode. All operations are allocation-free.

use fs_math::c64::C64;

pub(super) fn acoustic_wavenumber(
    omega: f64,
    sound_speed: f64,
    viscous_time: f64,
    thermal_time: f64,
    gamma: f64,
) -> Result<C64, &'static str> {
    if !omega.is_finite()
        || omega < 0.0
        || !sound_speed.is_finite()
        || sound_speed <= 0.0
        || !viscous_time.is_finite()
        || viscous_time < 0.0
        || !thermal_time.is_finite()
        || thermal_time < 0.0
        || !gamma.is_finite()
        || gamma < 1.0
    {
        return Err(
            "bulk dispersion needs finite passive transport, Cp/Cv >= 1 and nonnegative frequency",
        );
    }
    if omega == 0.0 {
        return Ok(C64::ZERO);
    }
    let x = omega * viscous_time;
    let y = omega * thermal_time;
    let thermal_coupling = (gamma - 1.0) * y;
    if !x.is_finite() || !y.is_finite() || !thermal_coupling.is_finite() {
        return Err("bulk dispersion transport scales are unrepresentable");
    }
    // z^-1 is a complex squared phase speed, relative to c². For gamma=1
    // thermal diffusion decouples EXACTLY, including at a thermal-root crossing.
    let inverse_z = if thermal_coupling == 0.0 {
        C64::new(1.0, -x)
    } else {
        let difference = x - gamma * y;
        // Expand b²+4a analytically: the large cross terms cancel before
        // evaluation. No subtraction of nearly equal quadratic roots occurs.
        let discriminant = C64::new(
            1.0 - difference * difference,
            -2.0 * (x + (gamma - 2.0) * y),
        );
        if !discriminant.re.is_finite() || !discriminant.im.is_finite() {
            return Err("bulk dispersion discriminant is unrepresentable");
        }
        let scale = 1.0 + difference * difference + discriminant.im.abs();
        // Along the negative real cut the acoustic continuation cannot be
        // distinguished from the thermal root. Roundoff in the cancelling
        // imaginary terms must not select opposite roots at adjacent omega.
        let imaginary_scale = x.abs() + ((gamma - 2.0) * y).abs();
        let unresolved_cut = discriminant.re < 0.0
            && discriminant.im.abs() <= 128.0 * f64::EPSILON * imaginary_scale;
        if discriminant.abs() <= 64.0 * f64::EPSILON * scale || unresolved_cut {
            return Err("acoustic and thermal dispersion branches are unresolved");
        }
        let root = discriminant.sqrt();
        // (b + sqrt(D))/2 = 1-i*x - 2i*(gamma-1)*y /
        // [sqrt(D) + 1-i*(x-gamma*y)]. This rationalized form retains tiny
        // thermal coupling when Cp/Cv is close to one. The denominator's
        // real part is >= 1, so it does not subtract nearly equal quantities.
        C64::new(1.0, -x)
            + C64::new(0.0, -2.0 * thermal_coupling) / (root + C64::new(1.0, -difference))
    };
    let wave = inverse_z.sqrt().recip().scale(omega / sound_speed);
    if !wave.re.is_finite()
        || wave.re <= 0.0
        || !wave.im.is_finite()
        || wave.im < 0.0
        || ((viscous_time > 0.0 || (gamma > 1.0 && thermal_time > 0.0)) && wave.im == 0.0)
    {
        return Err("bulk acoustic wavenumber is unrepresentable or nonpassive");
    }
    Ok(wave)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn relative(actual: C64, expected: C64) -> f64 {
        (actual - expected).abs() / expected.abs()
    }

    #[test]
    fn g1_independent_continuity_momentum_energy_equations() {
        // Normalized rho=c=omega=1; real constant coefficients. Eliminate
        // temperature in the original energy equation, then compare its
        // thermodynamic pressure with the original viscous momentum equation.
        // This checks the governing equations, not the quadratic solver again.
        for gamma in [1.000_000_01, 1.2, 1.4, 5.0 / 3.0, 2.5] {
            for viscous in [1e-7, 0.001, 0.1, 1.0, 10.0] {
                for thermal in [1e-7, 0.003, 0.2, 2.0] {
                    let k = acoustic_wavenumber(1.0, 1.0, viscous, thermal, gamma).unwrap();
                    let density = k; // continuity, with velocity amplitude 1
                    let temperature_pressure = k.scale((gamma - 1.0) / gamma)
                        / (C64::ONE + C64::new(0.0, 1.0) * k * k.scale(gamma * thermal));
                    let pressure = density.scale(1.0 / gamma) + temperature_pressure;
                    let momentum_pressure = k.recip() + C64::new(0.0, 1.0) * k.scale(viscous);
                    assert!(
                        relative(pressure, momentum_pressure) < 2e-12,
                        "gamma={gamma}, nu={viscous}, chi={thermal}, k={k:?}"
                    );
                    assert!(k.re > 0.0 && k.im > 0.0);
                }
            }
        }
    }

    #[test]
    fn g1_zero_thermal_coupling_is_exact_viscous_solution() {
        for x in [1e-12, 0.01, 1.0, 100.0] {
            // Independent real-arithmetic inverse square root of 1-i*x.
            let norm = 1.0_f64.hypot(x);
            let re = ((norm + 1.0) / (2.0 * norm * norm)).sqrt();
            let im = x / (2.0 * norm * norm * re);
            let expected = C64::new(re, im);
            for thermal in [0.0, 1e-6, 1.0, 100.0] {
                let wave = acoustic_wavenumber(1.0, 1.0, x, thermal, 1.0).unwrap();
                assert!(relative(wave, expected) < 8e-16);
            }
        }
        assert_eq!(
            acoustic_wavenumber(12.0, 3.0, 0.0, 0.0, 1.4).unwrap(),
            C64::from_re(4.0)
        );
        assert_eq!(
            acoustic_wavenumber(0.0, 3.0, 0.2, 0.3, 1.4).unwrap(),
            C64::ZERO
        );
    }

    #[test]
    fn g1_stokes_kirchhoff_is_second_order_weak_loss_limit() {
        let mut errors = Vec::new();
        for omega in [0.04, 0.02, 0.01] {
            let wave = acoustic_wavenumber(omega, 2.0, 0.7, 0.4, 1.4).unwrap();
            let classical = omega * omega / 4.0 * (0.7 + 0.4 * 0.4);
            errors.push((wave.im / classical - 1.0).abs());
            assert!(wave.re < omega / 2.0);
        }
        for pair in errors.windows(2) {
            assert!((pair[0] / pair[1] - 4.0).abs() < 0.02, "{errors:?}");
        }
    }

    #[test]
    fn g3_dimensionless_similarity_and_weak_thermal_loss_are_retained() {
        let base = acoustic_wavenumber(100.0, 4.0, 1e-3, 2e-3, 1.4).unwrap();
        let scaled = acoustic_wavenumber(700.0, 12.0, 1e-3 / 7.0, 2e-3 / 7.0, 1.4).unwrap();
        assert!(relative(scaled, base.scale(7.0 / 3.0)) < 5e-16);
        let gamma = 1.0 + 2.0 * f64::EPSILON;
        let tiny = acoustic_wavenumber(1.0, 1.0, 0.0, 1e-8, gamma).unwrap();
        let expected = 0.5 * (gamma - 1.0) * 1e-8;
        assert!((tiny.im / expected - 1.0).abs() < 2e-15, "{tiny:?}");
    }

    #[test]
    fn g0_bad_inputs_unrepresentable_loss_and_ambiguous_branches_refuse() {
        for omega in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(acoustic_wavenumber(omega, 1.0, 0.1, 0.2, 1.4).is_err());
        }
        assert!(acoustic_wavenumber(1.0, 1.0, 0.1, 0.2, 0.9).is_err());
        assert!(acoustic_wavenumber(f64::MAX, 1.0, 1.0, 1.0, 1.4).is_err());
        assert!(acoustic_wavenumber(f64::from_bits(1), 1.0, 0.1, 0.2, 1.4).is_err());
        // gamma=1.5, x=0.5*y: D is real and its roots coalesce at y=1.
        assert!(acoustic_wavenumber(1.0, 1.0, 0.5, 1.0, 1.5).is_err());
        assert!(acoustic_wavenumber(1.0, 1.0, 1.0, 2.0, 1.5).is_err());
        // Continuation is resolved on either side of that material ratio.
        for x in [0.999, 1.001] {
            let wave = acoustic_wavenumber(1.0, 1.0, x, 2.0, 1.5).unwrap();
            assert!(wave.re > 0.0 && wave.im > 0.0);
        }
    }

    #[test]
    fn g0_branch_cut_roundoff_cannot_switch_the_selected_mode() {
        let gamma = 1.4;
        let thermal_time = 0.314_159_265_35;
        let viscous_time = (2.0 - gamma) * thermal_time;
        // Cancellation in Im(D) has opposite signs at 5.1 and 6.2 in binary64.
        // Both are on the unresolved cut, not distinct admissible modes.
        for omega in [5.0, 5.1, 5.2, 6.1, 6.2, 6.3] {
            assert!(acoustic_wavenumber(omega, 1.0, viscous_time, thermal_time, gamma).is_err());
        }
    }
}
