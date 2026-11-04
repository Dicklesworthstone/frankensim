//! Shared arithmetic of original mean flat-plate relations and their strips.
//! Source domains and dimensional/evidence conversion stay with the parent.

use fs_math::det;

const LAMINAR: f64 = 0.664;
const TURBULENT: f64 = 0.037;
const LEADING_EDGE: f64 = 871.0;
const TURBULENT_POWER: f64 = 0.8;

// Preserve the original mean-card arithmetic, including operation ordering.
pub(crate) fn flat_plate_laminar_mean(reynolds: f64, prandtl: f64) -> f64 {
    LAMINAR * det::sqrt(reynolds) * det::pow(prandtl, 1.0 / 3.0)
}

pub(crate) fn flat_plate_turbulent_mean(reynolds: f64, prandtl: f64) -> f64 {
    (TURBULENT * det::pow(reynolds, TURBULENT_POWER) - LEADING_EDGE) * det::pow(prandtl, 1.0 / 3.0)
}

pub(crate) fn laminar_section(downstream_mean: f64, start_fraction: f64) -> f64 {
    downstream_mean / (1.0 + det::sqrt(start_fraction))
}

pub(crate) fn turbulent_section(downstream_mean: f64, prandtl: f64, start_fraction: f64) -> f64 {
    if start_fraction == 0.0 {
        return downstream_mean;
    }
    // Reuse the original mean instead of evaluating Re^0.8 twice. The source
    // domain bounds its leading-edge subtraction away from catastrophic loss.
    let uncorrected = downstream_mean + LEADING_EDGE * det::pow(prandtl, 1.0 / 3.0);
    // The input is the represented fraction itself, so det::ln retains its
    // exact difference from one. expm1 then avoids 1 - pow(f, 0.8) cancellation.
    let quotient = -det::expm1(TURBULENT_POWER * det::ln(start_fraction)) / (1.0 - start_fraction);
    uncorrected * quotient
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-13 * expected.abs(),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_original_means_and_laminar_local_response() {
        // Independent exactly soluble Pr=8, Re=10^4 laminar reference.
        let mean = flat_plate_laminar_mean(10_000.0, 8.0);
        close(mean, 132.8);
        assert_eq!(laminar_section(mean, 0.0), mean);
        close(laminar_section(mean, 0.25), 132.8 / 1.5);
        close(laminar_section(mean, 1.0_f64.next_down()), 66.4);
        // The original turbulent relation contains its laminar entrance.
        for (re, pr) in [(500_000.0_f64, 0.6_f64), (2_000_000.0, 8.0), (1.0e7, 60.0)] {
            close(
                flat_plate_turbulent_mean(re, pr),
                (0.037 * re.powf(0.8) - 871.0) * pr.cbrt(),
            );
        }
    }

    #[test]
    fn g1_g3_physical_subdivision_recovers_total_heat_without_resetting_the_entrance() {
        let positions = [0.0, 0.125, 0.2, 0.36, 0.5, 0.85, 1.0];
        for turbulent in [false, true] {
            let re = if turbulent { 8.0e6 } else { 4.0e5 };
            let pr = 8.0;
            let mean = |re| {
                if turbulent {
                    flat_plate_turbulent_mean(re, pr)
                } else {
                    flat_plate_laminar_mean(re, pr)
                }
            };
            let section = |mean, fraction| {
                if turbulent {
                    turbulent_section(mean, pr, fraction)
                } else {
                    laminar_section(mean, fraction)
                }
            };
            let mut integrated = 0.0;
            let mut prior_h = f64::INFINITY;
            for pair in positions.windows(2) {
                let [a, b] = [pair[0], pair[1]];
                let strip_nu = section(mean(re * b), a / b);
                integrated += strip_nu * ((b - a) / b);
                // For laminar strips, downstream convection weakens with the
                // growing boundary layer. This is a spatially varying h.
                if !turbulent {
                    assert!(strip_nu / b < prior_h);
                    prior_h = strip_nu / b;
                }
            }
            close(integrated, mean(re));
        }
        let reset_heat: f64 = positions
            .windows(2)
            .map(|pair| flat_plate_laminar_mean(4.0e5 * (pair[1] - pair[0]), 8.0))
            .sum();
        assert!(reset_heat > 2.0 * flat_plate_laminar_mean(4.0e5, 8.0));
    }

    #[test]
    fn g1_turbulent_short_strips_retain_the_differential_limit_and_first_corrections() {
        let re = 2.0e6_f64;
        let pr = 0.71_f64;
        let mean = flat_plate_turbulent_mean(re, pr);
        assert_eq!(turbulent_section(mean, pr, 0.0), mean);
        let power_amplitude = 0.037 * re.powf(0.8) * pr.cbrt();
        for bits in [8, 16, 24, 40, 52, 53] {
            let width = 2.0_f64.powi(-bits);
            // Independent binomial integral: p + p(1-p)w/2 +
            // p(1-p)(2-p)w^2/6 + ... . Here w <= 1/256.
            let mut term = 0.8;
            let mut integral = term;
            for n in 1..12 {
                term *= (f64::from(n) - 0.8) * width / f64::from(n + 1);
                integral += term;
            }
            close(
                turbulent_section(mean, pr, 1.0 - width),
                power_amplitude * integral,
            );
        }
        let fraction = 1.0_f64.next_down();
        let naive = (1.0 - fraction.powf(0.8)) / (1.0 - fraction);
        assert!(
            (naive - 0.8).abs() > 0.1,
            "fixture must expose subtractive cancellation"
        );
        close(
            turbulent_section(mean, pr, fraction),
            0.0296 * re.powf(0.8) * pr.cbrt(),
        );
    }

    #[test]
    fn g1_turbulent_cumulative_difference_and_streamwise_decay() {
        let re = 8.0e6_f64;
        let pr = 2.0_f64;
        let mean = flat_plate_turbulent_mean(re, pr);
        for f in [0.125_f64, 0.25, 0.5, 0.9] {
            let independent = 0.037 * pr.cbrt() * (re.powf(0.8) - (re * f).powf(0.8)) / (1.0 - f);
            close(turbulent_section(mean, pr, f), independent);
        }
        let first = turbulent_section(flat_plate_turbulent_mean(re * 0.5, pr), pr, 0.5) / 0.5;
        let later = turbulent_section(mean, pr, 0.75);
        assert!(first > later);
    }
}
