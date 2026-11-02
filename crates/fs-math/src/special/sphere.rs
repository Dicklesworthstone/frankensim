//! Spherical moments for interference integrals over the complete unit sphere.
//!
//! The Legendre integral and series are DLMF 10.54.2 and 10.53.1:
//! https://dlmf.nist.gov/10.54.E2 and https://dlmf.nist.gov/10.53.E1.
//! Small remainders are represented separately from their nonzero limits so
//! a weak radiating residual is not erased before the source contraction.

use crate::{dd::Dd, det};

/// Spherical Bessel values and cancellation-preserving remainders.
///
/// For a separation vector with magnitude x and unit direction e, the
/// normalized angular moments of exp(-i x n.e) are
/// M0 = j0, M1 = -i j1 e, and M2 = (j1/x) I - j2 e e^T.
/// At zero separation M0 = 1, M1 = 0, and M2 = I/3.
///
/// Keep the two remainder fields separate through a source contraction;
/// adding their constant limits and subtracting them again can erase
/// subwavelength radiation even in double-double arithmetic.
#[derive(Debug, Clone, Copy)]
pub struct SphericalBesselMoments {
    /// j0(x) - 1, even, with its continuous value zero at x = 0.
    pub j0_minus_one: Dd,
    /// j1(x), odd.
    pub j1: Dd,
    /// j1(x)/x - 1/3, even, with its continuous value zero at x = 0.
    pub j1_over_x_minus_one_third: Dd,
    /// j2(x), even, with its continuous value zero at x = 0.
    pub j2: Dd,
}

/// Evaluate the scalar, vector, and tensor spherical moments.
///
/// Accepts finite signed x; nonfinite inputs return None. For |x| <= 1,
/// eighteen fixed series updates accumulate the small remainders directly
/// in double-double arithmetic. Beyond that range the explicit spherical
/// formulas use one deterministic sine/cosine evaluation and reciprocal
/// powers, avoiding overflow from x*x. No allocation is performed.
///
/// The series truncation at the switch is below 1e-40; this does not bound
/// floating-point rounding. Representable range and underflow still apply,
/// and no uniform relative accuracy at nulls or interval enclosure is
/// claimed. The existing strict-core golden is unchanged.
#[must_use]
pub fn spherical_bessel_moments(x: f64) -> Option<SphericalBesselMoments> {
    if !x.is_finite() {
        return None;
    }

    let third = Dd::ONE / Dd::from_f64(3.0);
    let result = if x.abs() <= 1.0 {
        let xd = Dd::from_f64(x);
        let z = -(xd * xd);
        let mut scalar_term = Dd::ONE;
        let mut ratio_term = third;
        let mut quadrupole_term = -z / Dd::from_f64(15.0);
        let mut scalar_remainder = Dd::ZERO;
        let mut ratio_remainder = Dd::ZERO;
        let mut quadrupole = quadrupole_term;

        for n in 1..=18_u32 {
            let n = f64::from(n);
            scalar_term = scalar_term * z / Dd::from_f64(2.0 * n * (2.0 * n + 1.0));
            ratio_term = ratio_term * z / Dd::from_f64(2.0 * n * (2.0 * n + 3.0));
            quadrupole_term =
                quadrupole_term * z / Dd::from_f64(2.0 * n * (2.0 * n + 5.0));
            scalar_remainder = scalar_remainder + scalar_term;
            ratio_remainder = ratio_remainder + ratio_term;
            quadrupole = quadrupole + quadrupole_term;
        }

        SphericalBesselMoments {
            j0_minus_one: scalar_remainder,
            j1: xd * (third + ratio_remainder),
            j1_over_x_minus_one_third: ratio_remainder,
            j2: quadrupole,
        }
    } else {
        let (sine, cosine) = det::sin_cos(x);
        let inverse = Dd::ONE / Dd::from_f64(x);
        let sinc = Dd::from_f64(sine) * inverse;
        let j1 = (sinc - Dd::from_f64(cosine)) * inverse;
        let ratio = j1 * inverse;

        SphericalBesselMoments {
            j0_minus_one: sinc - Dd::ONE,
            j1,
            j1_over_x_minus_one_third: ratio - third,
            j2: Dd::from_f64(3.0) * ratio - sinc,
        }
    };

    [
        result.j0_minus_one,
        result.j1,
        result.j1_over_x_minus_one_third,
        result.j2,
    ]
    .iter()
    .all(|value| value.hi.is_finite() && value.lo.is_finite())
    .then_some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn values(x: f64) -> [f64; 4] {
        let m = spherical_bessel_moments(x).expect("finite spherical moments");
        [
            m.j0_minus_one.to_f64(),
            m.j1.to_f64(),
            m.j1_over_x_minus_one_third.to_f64(),
            m.j2.to_f64(),
        ]
    }

    fn close(actual: f64, expected: f64, absolute: f64) {
        assert!(
            (actual - expected).abs() <= absolute,
            "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}"
        );
    }

    #[test]
    fn finite_domain_zero_and_parity() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(spherical_bessel_moments(x).is_none());
        }
        for x in [0.0, -0.0] {
            for value in values(x) {
                assert_eq!(value, 0.0);
            }
        }
        for x in [1e-100, 1e-8, 0.5, 1.0, 1.01, 7.0, 1e100] {
            let positive = values(x);
            let negative = values(-x);
            for i in [0, 2, 3] {
                assert_eq!(positive[i], negative[i]);
            }
            assert_eq!(positive[1], -negative[1]);
        }
    }

    #[test]
    fn tiny_remainders_survive_below_the_constant_baselines() {
        for x in [1.0 / 1_073_741_824.0, 1e-20, 1e-100] {
            let [scalar, dipole, ratio, quadrupole] = values(x);
            let squared = x * x;
            for (actual, expected) in [
                (scalar, -squared / 6.0),
                (dipole, x / 3.0),
                (ratio, -squared / 30.0),
                (quadrupole, squared / 15.0),
            ] {
                assert_ne!(actual, 0.0, "representable remainder lost at x={x:e}");
                close(actual, expected, expected.abs() * 2e-14);
            }
        }
    }

    // Independent real angular integration, using platform trig only in
    // this oracle. It does not call any spherical Bessel implementation.
    fn sphere_average_polar(f: impl Fn(f64) -> f64) -> f64 {
        const STEPS: u32 = 20_000;
        let h = 2.0 / f64::from(STEPS);
        let mut sum = f(-1.0) + f(1.0);
        for i in 1..STEPS {
            let mu = -1.0 + f64::from(i) * h;
            let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
            sum += weight * f(mu);
        }
        sum * h / 6.0
    }

    #[test]
    fn scalar_vector_tensor_match_independent_sphere_quadrature() {
        let third = Dd::ONE / Dd::from_f64(3.0);
        for x in [0.2, 1.0, 1.0001, 3.0, 10.0, 25.0] {
            let m = spherical_bessel_moments(x).unwrap();
            let scalar = sphere_average_polar(|mu| (x * mu).cos());
            let dipole = sphere_average_polar(|mu| mu * (x * mu).sin());
            let axial = sphere_average_polar(|mu| mu * mu * (x * mu).cos());
            close((Dd::ONE + m.j0_minus_one).to_f64(), scalar, 2e-11);
            close(m.j1.to_f64(), dipole, 2e-11);
            close(
                (third + m.j1_over_x_minus_one_third - m.j2).to_f64(),
                axial,
                2e-11,
            );
        }
    }

    #[test]
    fn series_switch_is_continuous() {
        let center = values(1.0);
        for x in [1.0_f64.next_down(), 1.0_f64.next_up()] {
            let near = values(x);
            for i in 0..4 {
                close(near[i], center[i], 3e-15);
            }
        }
    }

    #[test]
    fn finite_huge_arguments_do_not_square_the_argument() {
        for x in [1e100, -1e100, 1e300, -1e300, f64::MAX, -f64::MAX] {
            for value in values(x) {
                assert!(value.is_finite(), "nonfinite moment at x={x:e}");
            }
        }
    }
}
