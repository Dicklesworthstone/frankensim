//! Positive local loss and stable analytic depth integration.
use fs_math::{c64::C64, det};

pub(super) fn electric_heat_density(electric: &[C64; 3], heat_scale: f64) -> f64 {
    // Scale before squaring: finite physical heat can have an unrepresentable
    // unweighted |E|² when the material loss coefficient is very small.
    electric
        .iter()
        .map(|value| value.scale(heat_scale).norm_sq())
        .sum()
}

pub(super) fn mode_absorbed_power(
    entering_w: f64,
    negative_kz_im: f64,
    start_m: f64,
    end_m: f64,
) -> f64 {
    // The owner validates a passive mode, finite start and ordered end. +inf
    // is a legal end. Depth before two avoids overflowing the coefficient
    // alone when a very strongly absorbing material has a very thin interval.
    let depth = end_m - start_m;
    if entering_w == 0.0 || depth == 0.0 || negative_kz_im == 0.0 {
        return 0.0;
    }
    let start_exponent = (negative_kz_im * start_m) * 2.0;
    let at_start = det::exp(start_exponent);
    let fraction = -det::expm1((negative_kz_im * depth) * 2.0);
    if at_start.is_normal() && fraction.is_normal() {
        return entering_w * at_start * fraction;
    }
    // A subnormal/vanished transmission or opacity can still deposit finite
    // watts when multiplied by a large entering power. Combine logarithms
    // before exponentiation; in the subnormal-opacity limit expm1(x)=x to
    // far better than binary64 precision. Reconstruct its product in log
    // space rather than taking the log of an already rounded subnormal.
    let log_fraction = if fraction.is_normal() {
        det::ln(fraction)
    } else {
        core::f64::consts::LN_2 + det::ln(-negative_kz_im) + det::ln(depth)
    };
    det::exp(det::ln(entering_w) + start_exponent + log_fraction)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn relative(a: f64, b: f64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-300),
            "{a:e} != {b:e}"
        );
    }
    #[test]
    fn local_material_heat_matches_oblique_maxwell_flux_for_te_and_tm() {
        let omega: f64 = 2.0;
        let mu: f64 = 1.3;
        for epsilon in [
            C64::new(2.1, -0.3),
            C64::new(-5.0, -2.0),
            C64::new(0.2, -1e-8),
        ] {
            for kx in [0.0, 0.7, 2.0] {
                let kz = (epsilon.scale(omega * omega * mu) - C64::from_re(kx * kx)).sqrt();
                let amplitude = C64::new(0.8, -0.3);
                let te = [C64::ZERO, amplitude, C64::ZERO];
                let tm = [amplitude, C64::ZERO, -(C64::from_re(kx) / kz) * amplitude];
                for (electric, admittance) in [
                    (te, kz.scale(1.0 / (omega * mu))),
                    (tm, epsilon.scale(omega) / kz),
                ] {
                    let entering = 0.5 * amplitude.norm_sq() * admittance.re;
                    let heat =
                        electric_heat_density(&electric, det::sqrt(-0.5 * omega * epsilon.im));
                    relative(heat, -2.0 * kz.im * entering, 2e-13);
                    let end = 0.4;
                    let intervals = 2048;
                    let mut volume = 0.0;
                    for j in 0..=intervals {
                        let depth = end * j as f64 / intervals as f64;
                        let weight = if j == 0 || j == intervals {
                            1.0
                        } else if j % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                        volume += weight * heat * det::exp(2.0 * kz.im * depth);
                    }
                    volume *= end / (3.0 * intervals as f64);
                    relative(
                        mode_absorbed_power(entering, kz.im, 0.0, end),
                        volume,
                        2e-11,
                    );
                }
            }
        }
    }
    #[test]
    fn weak_absorption_survives_endpoint_roundoff_and_partitions_additively() {
        let absorbed = mode_absorbed_power(0.5, -5e-15, 0.0, 1e-6);
        relative(absorbed, 5e-21, 2e-15);
        assert_eq!(0.5 - 0.5 * det::exp(-1e-20), 0.0);
        for attenuation in [-1e-20, -0.2, -1e3] {
            let first = mode_absorbed_power(0.7, attenuation, 0.0, 0.1);
            let rest = mode_absorbed_power(0.7, attenuation, 0.1, f64::INFINITY);
            relative(first + rest, 0.7, 3e-15);
            assert!(first >= 0.0 && rest >= 0.0);
        }
    }
    #[test]
    fn large_incident_power_preserves_representable_deep_and_thin_layer_heat() {
        // Independent 90-digit decimal exp(-800) and exp(-720), scaled by
        // 1e200. The first transmission is zero in binary64 before scaling;
        // the second is subnormal and loses significant relative precision.
        relative(
            mode_absorbed_power(1e200, -1.0, 400.0, f64::INFINITY),
            3.667_874_584_177_687e-148,
            3e-13,
        );
        relative(
            mode_absorbed_power(1e200, -1.0, 360.0, f64::INFINITY),
            2.032_230_802_424_293_2e-113,
            3e-13,
        );
        // Opacities 2e-400 and 2e-320 need multiplication by entering watts
        // before rounding; the expm1 correction is negligible at these scales.
        relative(
            mode_absorbed_power(1e200, -1e-200, 0.0, 1e-200),
            2e-200,
            3e-13,
        );
        relative(
            mode_absorbed_power(1e300, -1e-160, 0.0, 1e-160),
            2e-20,
            3e-13,
        );
    }
    #[test]
    fn opaque_and_large_coefficient_limits_remain_finite() {
        assert_eq!(mode_absorbed_power(0.7, -1.0, 1e300, f64::INFINITY), 0.0);
        assert_eq!(mode_absorbed_power(0.7, -1.0, 0.0, f64::INFINITY), 0.7);
        relative(
            mode_absorbed_power(0.7, -1e308, 0.0, 1e-308),
            0.7 * (1.0 - (-2.0f64).exp()),
            3e-15,
        );
    }
    #[test]
    fn weak_local_loss_does_not_overflow_an_unweighted_electric_norm() {
        let electric = [C64::new(1e200, -1e200), C64::from_re(2e200), C64::ZERO];
        assert!(electric[0].norm_sq().is_infinite());
        relative(electric_heat_density(&electric, 1e-150), 6e100, 3e-15);
        assert_eq!(electric_heat_density(&electric, 0.0), 0.0);
    }
}
