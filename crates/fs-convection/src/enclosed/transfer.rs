//! Original horizontal arithmetic and inclined Hollands air-gap transfer.
//!
//! Inclined equation authority: EnergyPlus 8.3 Engineering Reference, Solar
//! Collectors / Heat Transfer Coefficients, reproducing Hollands et al. (1976).
//! <https://bigladdersoftware.com/epx/docs/8-3/engineering-reference/solar-collectors.html#heat-transfer-coefficients>
//! The parent source card admits physical domains before this O(1) kernel.

use fs_math::det;

const RADIANS_PER_DEGREE: f64 = core::f64::consts::PI / 180.0;

// Original horizontal implementation moved without changing its arithmetic.
pub(crate) fn horizontal_nusselt(signed_rayleigh: f64) -> f64 {
    if signed_rayleigh <= 1708.0 {
        return 1.0;
    }
    let onset = 1.44 * ((signed_rayleigh - 1708.0) / signed_rayleigh);
    let high_rayleigh = if signed_rayleigh > 5830.0 {
        det::pow(signed_rayleigh / 5830.0, 1.0 / 3.0) - 1.0
    } else {
        0.0
    };
    1.0 + onset + high_rayleigh
}

pub(crate) fn horizontal_partial(signed_rayleigh: f64, axis: &str) -> Option<f64> {
    if axis == "Re" {
        return None;
    }
    if axis != "Ra" {
        return Some(0.0);
    }
    if signed_rayleigh == 1708.0 || signed_rayleigh == 5830.0 {
        return None;
    }
    if signed_rayleigh < 1708.0 {
        return Some(0.0);
    }
    let onset = (1.44 * 1708.0 / signed_rayleigh) / signed_rayleigh;
    let high_rayleigh = if signed_rayleigh > 5830.0 {
        det::pow(signed_rayleigh / 5830.0, 1.0 / 3.0) / (3.0 * signed_rayleigh)
    } else {
        0.0
    };
    Some(onset + high_rayleigh)
}

pub(crate) fn inclined_nusselt(rayleigh: f64, degrees: f64) -> f64 {
    if degrees == 0.0 {
        return horizontal_nusselt(rayleigh);
    }
    let angle = degrees * RADIANS_PER_DEGREE;
    let projected = rayleigh * det::cos(angle);
    if projected <= 1708.0 {
        return 1.0;
    }
    let tilted_onset = 1708.0 * det::pow(det::sin(1.8 * angle), 1.6);
    // The unclipped sine factor is nonnegative whenever the source's
    // positive-part onset is active: sin(1.8 theta)^1.6 <= 1 in our domain.
    // Subtract before dividing to retain the small onset excess.
    let onset =
        1.44 * ((projected - 1708.0) / projected) * ((projected - tilted_onset) / projected);
    let high_rayleigh = if projected > 5830.0 {
        det::pow(projected / 5830.0, 1.0 / 3.0) - 1.0
    } else {
        0.0
    };
    1.0 + onset + high_rayleigh
}

/// Absolute partials in full-gravity Ra and angle DEGREES, without division
/// by either input. Source kinks refuse; domain endpoints are handled upstream.
pub(crate) fn inclined_partial(rayleigh: f64, degrees: f64, axis: &str) -> Option<f64> {
    if axis == "Re" {
        return None;
    }
    if axis != "Ra" && axis != "inclination_degrees" {
        return Some(0.0);
    }
    if degrees == 0.0 {
        return if axis == "Ra" {
            horizontal_partial(rayleigh, axis)
        } else {
            Some(0.0)
        };
    }
    let angle = degrees * RADIANS_PER_DEGREE;
    let cosine = det::cos(angle);
    let projected = rayleigh * cosine;
    if projected == 1708.0 || projected == 5830.0 {
        return None;
    }
    if projected < 1708.0 {
        return Some(0.0);
    }
    let sine = det::sin(1.8 * angle);
    let tilted_onset = 1708.0 * det::pow(sine, 1.6);
    let first = (projected - 1708.0) / projected;
    let second = (projected - tilted_onset) / projected;
    let mut by_projected =
        (1.44 / projected) * ((1708.0 / projected) * second + first * (tilted_onset / projected));
    if projected > 5830.0 {
        by_projected += det::pow(projected / 5830.0, 1.0 / 3.0) / (3.0 * projected);
    }
    if axis == "Ra" {
        return Some(by_projected * cosine);
    }
    let projected_change = -rayleigh * det::sin(angle) * RADIANS_PER_DEGREE;
    let onset_change =
        1708.0 * 1.6 * det::pow(sine, 0.6) * det::cos(1.8 * angle) * (1.8 * RADIANS_PER_DEGREE);
    Some(by_projected * projected_change - 1.44 * first * (onset_change / projected))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual / expected - 1.0).abs() < 3e-13,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    // Independent 80-digit Python Decimal evaluation of the printed
    // EnergyPlus equation, with Taylor-series sin/cos and Decimal ln/exp
    // powers. Derivative references use central differences with step 1e-20
    // in that arithmetic, not fs-math or the analytic implementation below.
    const REFERENCE: [(f64, f64, f64, f64, f64); 7] = [
        (
            4000.0,
            30.0,
            1.473_578_270_121_595_1,
            0.000_179_256_471_743_889_9,
            -0.016_589_669_847_445_275,
        ),
        (
            4000.0,
            45.0,
            1.232_724_625_844_332_7,
            0.000_173_118_374_462_818_06,
            -0.014_774_493_269_449_2,
        ),
        (
            4000.0,
            60.0,
            1.044_547_313_539_069_3,
            0.000_106_565_989_701_474_48,
            -0.010_179_821_359_751_9,
        ),
        (
            10_000.0,
            30.0,
            2.134_580_810_339_825_6,
            0.000_078_685_519_424_621_3,
            -0.013_860_556_662_474_592,
        ),
        (
            10_000.0,
            45.0,
            1.899_983_267_836_469_4,
            0.000_087_957_625_570_388_23,
            -0.017_410_556_438_046_24,
        ),
        (
            10_000.0,
            60.0,
            1.649_213_630_778_42,
            0.000_063_571_617_279_677_49,
            -0.014_336_273_000_233_435,
        ),
        (
            100_000.0,
            60.0,
            3.393_883_225_052_502_6,
            0.000_007_737_903_125_204_401,
            -0.022_675_592_927_924_86,
        ),
    ];

    #[test]
    fn inclined_source_matches_independent_values_and_both_input_partials() {
        for (rayleigh, degrees, expected, by_rayleigh, by_degree) in REFERENCE {
            near(inclined_nusselt(rayleigh, degrees), expected);
            near(
                inclined_partial(rayleigh, degrees, "Ra").unwrap(),
                by_rayleigh,
            );
            near(
                inclined_partial(rayleigh, degrees, "inclination_degrees").unwrap(),
                by_degree,
            );
        }
    }

    #[test]
    fn zero_tilt_preserves_original_horizontal_arithmetic_and_conduction() {
        for (rayleigh, expected) in [
            (0.0, 1.0),
            (1708.0, 1.0),
            (2000.0, 1.210_24),
            (3416.0, 1.72),
            (5830.0, 2.018_126_929_674_099_5),
            (10_000.0, 2.391_093_013_195_09),
            (100_000.0, 3.994_360_101_957_176),
        ] {
            near(inclined_nusselt(rayleigh, 0.0), expected);
            assert_eq!(
                inclined_nusselt(rayleigh, 0.0).to_bits(),
                horizontal_nusselt(rayleigh).to_bits()
            );
            assert_eq!(
                inclined_partial(rayleigh, 0.0, "Ra"),
                horizontal_partial(rayleigh, "Ra")
            );
        }
        for degrees in [0.0, 30.0, 45.0, 60.0] {
            for rayleigh in [0.0, 1e-300, 1000.0, 1700.0] {
                assert_eq!(inclined_nusselt(rayleigh, degrees), 1.0);
                assert_eq!(inclined_partial(rayleigh, degrees, "Ra"), Some(0.0));
                assert_eq!(
                    inclined_partial(rayleigh, degrees, "inclination_degrees"),
                    Some(0.0)
                );
            }
        }
        assert_eq!(horizontal_nusselt(-f64::MAX), 1.0);
        assert_eq!(horizontal_partial(-10_000.0, "Ra"), Some(0.0));
    }

    #[test]
    fn angle_correction_cannot_be_replaced_with_gravity_projection_alone() {
        let rayleigh = 4000.0;
        let degrees = 45.0;
        let projected_only = horizontal_nusselt(rayleigh * det::cos(degrees * RADIANS_PER_DEGREE));
        assert!(projected_only - inclined_nusselt(rayleigh, degrees) > 0.3);
        // Full-gravity Ra can be supercritical while its normal component
        // remains subcritical. No unprojected onset or double-counted Fourier
        // term may add convection to this layer.
        assert_eq!(inclined_nusselt(3000.0, 60.0), 1.0);
        assert!(horizontal_nusselt(3000.0) > 1.0);
    }

    #[test]
    fn projected_branch_switches_are_continuous_and_refuse_smooth_partials() {
        let degrees = 30.0;
        let cosine = det::cos(degrees * RADIANS_PER_DEGREE);
        for kink in [1708.0, 5830.0] {
            let rayleigh = kink / cosine;
            assert_eq!(rayleigh * cosine, kink);
            assert_eq!(inclined_partial(rayleigh, degrees, "Ra"), None);
            assert_eq!(
                inclined_partial(rayleigh, degrees, "inclination_degrees"),
                None
            );
            assert_eq!(inclined_partial(rayleigh, degrees, "Pr"), Some(0.0));
            let at = inclined_nusselt(rayleigh, degrees);
            near(inclined_nusselt(rayleigh.next_down(), degrees), at);
            near(inclined_nusselt(rayleigh.next_up(), degrees), at);
        }
        for degrees in [0.0, 30.0, 45.0, 60.0] {
            let mut previous = 1.0;
            for rayleigh in [
                0.0, 1708.0, 2000.0, 3416.0, 4000.0, 5830.0, 10_000.0, 100_000.0,
            ] {
                let value = inclined_nusselt(rayleigh, degrees);
                assert!(value >= previous && value.is_finite());
                previous = value;
            }
        }
    }
}
