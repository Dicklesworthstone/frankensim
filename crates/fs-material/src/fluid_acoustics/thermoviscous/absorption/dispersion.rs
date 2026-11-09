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

/// The same spatial acoustic branch with an ideal gas's original complex
/// caloric response. `cv_over_r = cv(omega)/R_specific`; viscosity and
/// conductivity times are respectively `4 mu/(3 rho c0^2)` and
/// `kappa/(rho R_specific c0^2)`. No molecular loss enters either time.
///
/// Eliminating temperature from `-i omega rho cv T' + kappa k^2 T' =
/// -i k p u`, together with continuity, viscous momentum and `p = rho R T`,
/// gives `A z^2 + B z - 1 = 0`, for `z=(k c0/omega)^2`, where
/// `g=(1+R/cv)/gamma0`, `x=omega*viscous_time`, `y=omega*kappa/(rho cv c0^2)`,
/// `A=(x+i/gamma0)y` and `B=g-i(x+y)`. This uses the SAME relaxing heat
/// capacity in both compression and heat conduction; an equilibrium-loss
/// multiplier on the molecular wavenumber does not satisfy these equations.
pub(crate) fn relaxing_acoustic_wavenumber(
    omega: f64,
    sound_speed: f64,
    viscous_time: f64,
    conductivity_time: f64,
    equilibrium_gamma: f64,
    cv_over_r: C64,
) -> Result<C64, &'static str> {
    if !omega.is_finite()
        || omega < 0.0
        || !sound_speed.is_finite()
        || sound_speed <= 0.0
        || !viscous_time.is_finite()
        || viscous_time < 0.0
        || !conductivity_time.is_finite()
        || conductivity_time < 0.0
        || !equilibrium_gamma.is_finite()
        || equilibrium_gamma <= 1.0
        || !cv_over_r.re.is_finite()
        || cv_over_r.re <= 0.0
        || !cv_over_r.im.is_finite()
        || cv_over_r.im < 0.0
    {
        return Err("relaxing dispersion needs finite passive transport and caloric response");
    }
    if omega == 0.0 {
        return Ok(C64::ZERO);
    }
    let inverse_cv = cv_over_r.recip();
    let r = inverse_cv.scale(1.0 / equilibrium_gamma);
    let g = C64::from_re(1.0 / equilibrium_gamma) + r;
    if cv_over_r.im == 0.0 {
        // A real caloric response has the existing classical branch policy,
        // with its actual (possibly frozen) sound speed and Cp. Preserve the
        // full strong-loss domain established by that owner.
        return acoustic_wavenumber(
            omega,
            sound_speed * fs_math::det::sqrt(g.re),
            viscous_time / g.re,
            conductivity_time / (cv_over_r.re + 1.0) / g.re,
            1.0 + inverse_cv.re,
        );
    }
    let x = omega * viscous_time;
    let y = inverse_cv.scale(omega * conductivity_time);
    if !x.is_finite() || !y.re.is_finite() || !y.im.is_finite() {
        return Err("relaxing dispersion transport scales are unrepresentable");
    }
    let viscous = g + C64::new(0.0, -x);
    let inverse_z = if conductivity_time == 0.0 {
        viscous
    } else {
        // D = u^2 - 4 i y r, u = g-i(x-y). Continue its square root from
        // g as BOTH transport coefficients increase at this fixed frequency.
        let u = viscous + C64::new(-y.im, y.re);
        let yr = y * r;
        let correction = C64::new(4.0 * yr.im, -4.0 * yr.re);
        let discriminant = u * u + correction;
        let scale = u.norm_sq() + correction.abs();
        if !discriminant.re.is_finite() || !discriminant.im.is_finite() || !scale.is_finite() {
            return Err("relaxing dispersion discriminant is unrepresentable");
        }
        let product = g * (C64::from_re(x) + y);
        let d1 = C64::new(2.0 * product.im, -2.0 * product.re)
            + C64::new(
                -4.0 * y.im / equilibrium_gamma,
                4.0 * y.re / equilibrium_gamma,
            );
        let difference = C64::from_re(x) - y;
        let d2 = -(difference * difference);
        let root = continued_discriminant_sqrt(g, d1, d2, discriminant)?;
        let denominator = root + u;
        if denominator.abs() >= 0.5 * (root.abs() + u.abs()) {
            // Retain arbitrarily weak thermal coupling without subtracting
            // nearby roots. A continued root need not be on the principal
            // sheet, so the old positive-real denominator guarantee is gone.
            viscous + correction.scale(0.5) / denominator
        } else {
            // w=z^-1 satisfies w^2-Bw-A=0. After a sheet change use whichever
            // equivalent formula adds the separated roots instead of
            // subtracting them; the product of the two w roots is -A.
            let b = viscous + C64::new(y.im, -y.re);
            let sum = b + root;
            let difference = b - root;
            if sum.abs() >= difference.abs() {
                sum.scale(0.5)
            } else {
                let a = (C64::from_re(x) + C64::new(0.0, 1.0 / equilibrium_gamma)) * y;
                -a / difference.scale(0.5)
            }
        }
    };
    let wave = inverse_z.sqrt().recip().scale(omega / sound_speed);
    if !wave.re.is_finite()
        || wave.re <= 0.0
        || !wave.im.is_finite()
        || wave.im < 0.0
        || ((viscous_time > 0.0 || conductivity_time > 0.0 || cv_over_r.im > 0.0) && wave.im == 0.0)
    {
        return Err("relaxing acoustic wavenumber is unrepresentable or nonpassive");
    }
    Ok(wave)
}

/// Continue sqrt(D(t)) from `start_root`, for the quadratic discriminant
/// D(t)=start_root^2 + linear*t + quadratic*t^2, 0<=t<=1. This is transport
/// continuation at ONE frequency, not a global causal-analyticity claim.
///
/// On each interval [a,b], the quadratic lies in the convex hull of its three
/// exact Bezier controls D(a), D(a)+(b-a)D'(a)/2, D(b). Enclose that hull in
/// a disk centered at D(a), with radius at most |D(a)|/2. The whole interval
/// then avoids zero and changes argument by at most pi/6; the square root
/// with positive projection on the preceding root is uniquely continued.
/// This tests an enclosure of the interval, not just sampled endpoints.
///
/// The coefficient-scaled roundoff margin also rejects near coalescence.
/// Storage is constant and at most 256 attempted intervals are evaluated;
/// unresolved or excessively conditioned paths refuse rather than selecting
/// a mode from passivity alone. No previous frequency evaluation is needed.
fn continued_discriminant_sqrt(
    start_root: C64,
    linear: C64,
    quadratic: C64,
    end_value: C64,
) -> Result<C64, &'static str> {
    let constant = start_root * start_root;
    let scale = constant.abs() + linear.abs() + quadratic.abs();
    if !scale.is_finite() || scale == 0.0 {
        return Err("relaxing dispersion continuation is unrepresentable");
    }
    let margin = 512.0 * f64::EPSILON * scale;
    // The common weak/moderate-loss case needs no subdivision: an analytic
    // minimum of Re(D) establishes the principal sheet on the whole path.
    let mut minimum_real = constant.re.min(end_value.re);
    if quadratic.re > 0.0 {
        let vertex = -0.5 * (linear.re / quadratic.re);
        if vertex > 0.0 && vertex < 1.0 {
            minimum_real =
                minimum_real.min(constant.re + vertex * (linear.re + vertex * quadratic.re));
        }
    }
    if minimum_real > margin {
        return Ok(end_value.sqrt());
    }
    let mut at: f64 = 0.0;
    let mut step: f64 = 1.0;
    let mut value = constant;
    let mut root = start_root;
    for _ in 0..256 {
        let next = (at + step).min(1.0);
        let width = next - at;
        if width == 0.0 || value.abs() <= margin {
            return Err("relaxing acoustic and thermal branches are unresolved");
        }
        let next_value = if next == 1.0 {
            end_value
        } else {
            constant + (linear + quadratic.scale(next)).scale(next)
        };
        let control = value + (linear + quadratic.scale(2.0 * at)).scale(0.5 * width);
        let radius = (control - value).abs().max((next_value - value).abs());
        let magnitude = value.abs();
        if radius <= 0.5 * magnitude && magnitude - radius > margin {
            let candidate = next_value.sqrt();
            // The ratio stays near a unit complex value on an accepted step;
            // it avoids overflowing an unscaled complex dot product.
            root = if (candidate / root).re >= 0.0 {
                candidate
            } else {
                -candidate
            };
            if next == 1.0 {
                return Ok(root);
            }
            at = next;
            value = next_value;
            step = (2.0 * width).min(1.0 - at);
        } else {
            step = 0.5 * width;
        }
    }
    Err("relaxing acoustic continuation exceeds its bounded resolution")
}

pub(super) fn acoustic_wavenumber(
    omega: f64,
    sound_speed: f64,
    viscous_time: f64,
    thermal_time: f64,
    gamma: f64,
) -> Result<C64, &'static str> {
    acoustic_response(omega, sound_speed, viscous_time, thermal_time, gamma).map(|wave| wave.0)
}

/// Acoustic k and the traction/velocity and pressure/velocity impedances,
/// normalized by rho*c. Keep the caloric part of the squared phase speed
/// separate: recovering it by adding i*x after forming 1-i*x loses pressure
/// in the strongly viscous limit, even while k remains accurately resolved.
/// The fourth value retains q+i*gamma*y, q=(omega/(k*c))², BEFORE
/// cancellation at coincident viscous/thermal scales. Only a thermal-field
/// consumer divides by it; ordinary k/pressure/traction callers do no extra solve.
pub(super) fn acoustic_response(
    omega: f64,
    sound_speed: f64,
    viscous_time: f64,
    thermal_time: f64,
    gamma: f64,
) -> Result<(C64, C64, C64, C64), &'static str> {
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
        return Ok((C64::ZERO, C64::ONE, C64::ONE, C64::ONE));
    }
    let x = omega * viscous_time;
    let y = omega * thermal_time;
    let thermal_coupling = (gamma - 1.0) * y;
    if !x.is_finite() || !y.is_finite() || !thermal_coupling.is_finite() {
        return Err("bulk dispersion transport scales are unrepresentable");
    }
    // z^-1 is a complex squared phase speed, relative to c². For gamma=1
    // thermal diffusion decouples EXACTLY, including at a thermal-root crossing.
    let (inverse_z, caloric, thermal_denominator) = if thermal_coupling == 0.0 {
        (
            C64::new(1.0, -x),
            C64::ONE,
            C64::new(1.0, -(x - gamma * y)),
        )
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
        let denominator = root + C64::new(1.0, -difference);
        let correction = C64::new(0.0, -2.0 * thermal_coupling) / denominator;
        (
            C64::new(1.0, -x) + correction,
            C64::ONE + correction,
            denominator.scale(0.5),
        )
    };
    let traction = inverse_z.sqrt();
    let inverse_traction = traction.recip();
    let wave = inverse_traction.scale(omega / sound_speed);
    if !wave.re.is_finite()
        || wave.re <= 0.0
        || !wave.im.is_finite()
        || wave.im < 0.0
        || ((viscous_time > 0.0 || (gamma > 1.0 && thermal_time > 0.0)) && wave.im == 0.0)
    {
        return Err("bulk acoustic wavenumber is unrepresentable or nonpassive");
    }
    Ok((wave, traction, caloric * inverse_traction, thermal_denominator))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn relative(actual: C64, expected: C64) -> f64 {
        (actual - expected).abs() / expected.abs()
    }

    #[test]
    fn g1_relaxing_wave_closes_original_mass_momentum_and_energy() {
        // rho=R=c0=1, T=p=1/gamma0, velocity amplitude 1. A single
        // independent reservoir with frozen Cv=2 and equilibrium Cv=5.
        let gamma = 1.2;
        for omega in [0.001, 0.1, 1.0, 10.0, 100.0] {
            let cv = C64::new(
                2.0 + 3.0 / (1.0 + omega * omega),
                3.0 * omega / (1.0 + omega * omega),
            );
            for viscosity in [1e-9, 1e-5, 0.002] {
                for conductivity in [1e-9, 3e-5, 0.007] {
                    let k = relaxing_acoustic_wavenumber(
                        omega,
                        1.0,
                        viscosity,
                        conductivity,
                        gamma,
                        cv,
                    )
                    .unwrap();
                    let density = k.scale(omega.recip());
                    let temperature = k.scale(gamma.recip())
                        / (cv.scale(omega) + C64::new(0.0, conductivity) * k * k);
                    let pressure = density.scale(gamma.recip()) + temperature;
                    let momentum_pressure = k.recip().scale(omega) + C64::new(0.0, viscosity) * k;
                    assert!(
                        relative(pressure, momentum_pressure) < 3e-12,
                        "omega={omega}, nu={viscosity}, kappa={conductivity}, k={k:?}"
                    );
                    assert!(k.re > 0.0 && k.im > 0.0);
                }
            }
        }
    }

    #[test]
    fn g3_equilibrium_multiplier_fails_the_coupled_equations() {
        let cv = C64::new(3.5, 1.5);
        let gamma: f64 = 1.2;
        let molecular = (C64::ONE + cv.recip()).scale(gamma.recip()).sqrt().recip();
        let delta = 0.2 + 2.0 * (gamma - 1.0) * (gamma - 1.0) / gamma;
        let old = molecular * C64::new(1.0, -delta).sqrt().recip();
        let pressure = old.scale(gamma.recip())
            + old.scale(gamma.recip()) / (cv + C64::new(0.0, 2.0) * old * old);
        let momentum_pressure = old.recip() + C64::new(0.0, 0.2) * old;
        assert!(relative(pressure, momentum_pressure) > 0.05);
        let coupled = relaxing_acoustic_wavenumber(1.0, 1.0, 0.2, 2.0, gamma, cv).unwrap();
        assert!(relative(coupled, old) > 0.02);
    }

    #[test]
    fn g1_zero_reservoir_limit_is_full_classical_dispersion() {
        for gamma in [1.1, 1.4, 5.0 / 3.0] {
            for omega in [1e-8, 0.1, 1.0, 10.0] {
                let cv = C64::from_re(1.0 / (gamma - 1.0));
                for thermal in [0.0, 0.1, 0.7] {
                    let coupled = relaxing_acoustic_wavenumber(
                        omega,
                        2.0,
                        0.2,
                        thermal * gamma / (gamma - 1.0),
                        gamma,
                        cv,
                    )
                    .unwrap();
                    let classical = acoustic_wavenumber(omega, 2.0, 0.2, thermal, gamma).unwrap();
                    assert!(relative(coupled, classical) < 3e-14);
                    if classical.im != 0.0 {
                        assert!((coupled.im / classical.im - 1.0).abs() < 3e-14);
                    }
                }
            }
        }
    }

    #[test]
    fn g3_relaxing_dispersion_limits_and_dimensional_similarity() {
        let cv = C64::new(3.5, 1.5);
        let base = relaxing_acoustic_wavenumber(1.0, 2.0, 0.02, 0.07, 1.2, cv).unwrap();
        let scaled = relaxing_acoustic_wavenumber(7.0, 6.0, 0.02 / 7.0, 0.01, 1.2, cv).unwrap();
        assert!(relative(scaled, base.scale(7.0 / 3.0)) < 1e-15);
        let lossless_transport = relaxing_acoustic_wavenumber(1.0, 2.0, 0.0, 0.0, 1.2, cv).unwrap();
        let molecular = (C64::ONE + cv.recip())
            .scale(1.0 / 1.2)
            .sqrt()
            .recip()
            .scale(0.5);
        assert!(relative(lossless_transport, molecular) < 3e-16);
        assert_eq!(
            relaxing_acoustic_wavenumber(0.0, 2.0, 0.2, 0.7, 1.2, cv).unwrap(),
            C64::ZERO
        );
    }

    #[test]
    fn g0_relaxing_dispersion_refuses_invalid_or_unresolved_modes() {
        let cv = C64::new(3.5, 1.5);
        for omega in [-1.0, f64::NAN, f64::INFINITY, f64::MAX, f64::from_bits(1)] {
            assert!(relaxing_acoustic_wavenumber(omega, 2.0, 0.2, 0.7, 1.2, cv).is_err());
        }
        for cv in [
            C64::ZERO,
            C64::new(1.0, -1.0),
            C64::new(f64::NAN, 1.0),
            C64::new(1.0, f64::INFINITY),
        ] {
            assert!(relaxing_acoustic_wavenumber(1.0, 2.0, 0.2, 0.7, 1.2, cv).is_err());
        }
        // The same acoustic/thermal coalescence as the nonrelaxing owner.
        assert!(relaxing_acoustic_wavenumber(1.0, 1.0, 0.5, 3.0, 1.5, C64::from_re(2.0)).is_err());
        assert!(relaxing_acoustic_wavenumber(1.0, 1.0, 1.0, 6.0, 1.5, C64::from_re(2.0)).is_err());
    }

    #[test]
    fn g1_molecular_thermal_cut_crossing_continues_the_acoustic_mode() {
        // Positive Cv_inf=2R and reservoir C=0.25R. The transport ratio is
        // Eucken-consistent. A principal square root without continuation
        // jumps by about 25% between these neighboring frequencies, although
        // both roots have positive attenuation. The continued root is smooth.
        let mut waves = Vec::new();
        let mut principal_waves = Vec::new();
        for omega in [49.425, 49.430, 100.0, 1000.0] {
            let wt = 0.03 * omega;
            let cv = C64::new(2.0 + 0.25 / (1.0 + wt * wt), 0.25 * wt / (1.0 + wt * wt));
            let gamma = 13.0 / 9.0;
            let k = relaxing_acoustic_wavenumber(omega, 1.0, 0.16, 0.54, gamma, cv).unwrap();
            // Independent reference: follow the two roots of the ORIGINAL
            // quadratic in inverse_z through 8192 small transport steps.
            // The step-size assertion establishes a uniquely nearest root;
            // this does not use the production discriminant enclosure.
            let g = (C64::ONE + cv.recip()).scale(gamma.recip());
            let mut inverse_z = g;
            for step in 1..=8192 {
                let fraction = f64::from(step) / 8192.0;
                let x = omega * 0.16 * fraction;
                let y = cv.recip().scale(omega * 0.54 * fraction);
                let a = (C64::from_re(x) + C64::new(0.0, gamma.recip())) * y;
                let b = g + C64::new(y.im, -x - y.re);
                let root = (b * b + a.scale(4.0)).sqrt();
                let first = (b + root).scale(0.5);
                let second = (b - root).scale(0.5);
                let next = if (first - inverse_z).abs() < (second - inverse_z).abs() {
                    first
                } else {
                    second
                };
                assert!((next - inverse_z).abs() < 0.05 * (first - second).abs());
                inverse_z = next;
                if step == 8192 {
                    principal_waves.push(first.sqrt().recip().scale(omega));
                }
            }
            let expected = inverse_z.sqrt().recip().scale(omega);
            assert!(
                relative(k, expected) < 3e-14,
                "omega={omega}: {k:?} vs {expected:?}"
            );
            waves.push(k);
        }
        assert!(relative(waves[1], waves[0]) < 1e-4);
        assert!(relative(principal_waves[1], principal_waves[0]) > 0.2);
    }

    #[test]
    fn g1_strong_molecular_continuation_closes_physics_and_scales() {
        let gamma: f64 = 13.0 / 9.0;
        for omega in [49.425, 49.430, 100.0, 1000.0, 1e6] {
            let wt = 0.03 * omega;
            let cv = C64::new(2.0 + 0.25 / (1.0 + wt * wt), 0.25 * wt / (1.0 + wt * wt));
            let k = relaxing_acoustic_wavenumber(omega, 1.0, 0.16, 0.54, gamma, cv).unwrap();
            let g = (C64::ONE + cv.recip()).scale(gamma.recip());
            let x = omega * 0.16;
            let y = cv.recip().scale(omega * 0.54);
            let a = (C64::from_re(x) + C64::new(0.0, gamma.recip())) * y;
            let b = g + C64::new(y.im, -x - y.re);
            let z = k.scale(omega.recip()) * k.scale(omega.recip());
            let quadratic = a * z * z;
            let linear = b * z;
            assert!(
                (quadratic + linear - C64::ONE).abs()
                    < 3e-14 * (quadratic.abs() + linear.abs() + 1.0)
            );
            // Original continuity, energy, EOS and momentum equations with
            // rho=R=c0=1, T=p=1/gamma and unit peak longitudinal velocity.
            let density = k.scale(omega.recip());
            let momentum_pressure = k.recip().scale(omega) + C64::new(0.0, 0.16) * k;
            let temperature = momentum_pressure - density.scale(gamma.recip());
            let caloric = cv.scale(omega) * temperature;
            let conduction = C64::new(0.0, 0.54) * k * k * temperature;
            let compression = k.scale(gamma.recip());
            // At extreme transport the caloric and conduction terms nearly
            // cancel. Check the original energy residual, rather than
            // amplifying rounding by dividing by that small difference.
            assert!(
                (caloric + conduction - compression).abs()
                    < 3e-14 * (caloric.abs() + conduction.abs() + compression.abs()),
                "omega={omega}",
            );
            if omega <= 1000.0 {
                let energy_temperature =
                    k.scale(gamma.recip()) / (cv.scale(omega) + C64::new(0.0, 0.54) * k * k);
                let energy_pressure = density.scale(gamma.recip()) + energy_temperature;
                assert!(relative(energy_pressure, momentum_pressure) < 3e-13);
            }
            assert!(k.re > 0.0 && k.im > 0.0);
            let scaled =
                relaxing_acoustic_wavenumber(7.0 * omega, 3.0, 0.16 / 7.0, 0.54 / 7.0, gamma, cv)
                    .unwrap();
            assert!(relative(scaled, k.scale(7.0 / 3.0)) < 3e-14);
        }
    }

    #[test]
    fn g3_quadratic_continuation_tracks_exact_sheet_and_scale() {
        // D(t)=(g+t*q)^2 has the known continuous square root g+t*q.
        // It crosses the negative-real cut, including exactly at t=1.
        let g = C64::new(1.0, -0.125);
        for q in [-1.0, -1.2, -3.0] {
            let q = C64::from_re(q);
            for scale in [1e-100, 1.0, 1e100] {
                let root_scale = fs_math::det::sqrt(scale);
                let expected = (g + q).scale(root_scale);
                let actual = continued_discriminant_sqrt(
                    g.scale(root_scale),
                    (g * q).scale(2.0 * scale),
                    (q * q).scale(scale),
                    ((g + q) * (g + q)).scale(scale),
                )
                .unwrap();
                assert!(relative(actual, expected) < 3e-14, "q={q:?}, scale={scale}");
            }
        }
    }

    #[test]
    fn g0_quadratic_continuation_refuses_coalescence_and_roundoff_neighbours() {
        // D(t)=(1-2t)^2-i*offset. A nonzero endpoint does not establish
        // separated modes in the interior of the transport path.
        for offset in [0.0, 1e-15, 1e-13] {
            let endpoint = C64::new(1.0, -offset);
            assert!(
                continued_discriminant_sqrt(
                    endpoint.sqrt(),
                    C64::from_re(-4.0),
                    C64::from_re(4.0),
                    endpoint,
                )
                .is_err()
            );
        }
        let endpoint = C64::new(1.0, -1e-4);
        let root = continued_discriminant_sqrt(
            endpoint.sqrt(),
            C64::from_re(-4.0),
            C64::from_re(4.0),
            endpoint,
        )
        .unwrap();
        assert!(relative(root, endpoint.sqrt()) < 1e-14);
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
