//! Peak physical strain of a complex harmonic ellipse, in original units.

use fs_math::{c64::C64, dd::Dd, det};

// A scaled two-column singular value avoids both overflowing squared strains
// and the overly conservative Frobenius norm of the complex phasor itself.
pub(super) fn peak_norm<const N: usize>(strain: &[C64; N]) -> Option<f64> {
    if strain
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return None;
    }
    let scale = strain
        .iter()
        .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
    if scale == 0.0 {
        return Some(0.0);
    }
    let (mut a, mut b, mut c) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for z in strain {
        let re = Dd::from_f64(z.re / scale);
        let im = Dd::from_f64(z.im / scale);
        a = a + re * re;
        b = b + im * im;
        c = c + re * im;
    }
    let trace = (a + b).to_f64();
    let gap = det::hypot((a - b).to_f64(), 2.0 * c.to_f64());
    let peak = scale * det::sqrt(0.5 * (trace + gap));
    peak.is_finite().then_some(peak)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_strain_does_not_sum_independent_component_maxima() {
        assert_eq!(
            peak_norm(&[C64::from_re(0.75), C64::new(0.0, 0.75)]),
            Some(0.75)
        );
        assert_eq!(
            peak_norm(&[C64::new(0.75, 0.75)]),
            Some(0.75 * det::sqrt(2.0))
        );
    }

    #[test]
    fn general_ellipse_matches_analytic_gram_eigenvalue() {
        // A=25, B=9, C=12; neither quadrature phase attains the maximum.
        let expected = det::sqrt(0.5 * (34.0 + det::sqrt(832.0)));
        let actual = peak_norm(&[C64::new(3.0, 0.0), C64::new(4.0, 3.0)]).unwrap();
        assert!((actual - expected).abs() < 2e-15 * expected);
        let rotation = C64::new(0.6, 0.8);
        let rotated = [C64::new(3.0, 0.0) * rotation, C64::new(4.0, 3.0) * rotation];
        assert!((peak_norm(&rotated).unwrap() - actual).abs() < 2e-15 * actual);
    }

    #[test]
    fn retains_finite_extreme_and_subnormal_physical_scales() {
        for scale in [f64::from_bits(1), 1e-250, 1e250, f64::MAX] {
            assert_eq!(
                peak_norm(&[C64::from_re(scale), C64::new(0.0, scale)]),
                Some(scale)
            );
        }
        assert_eq!(peak_norm(&[C64::new(f64::MAX, f64::MAX)]), None);
    }

    #[test]
    fn zero_and_nonfinite_inputs_have_explicit_outcomes() {
        assert_eq!(peak_norm(&[C64::ZERO; 6]), Some(0.0));
        assert_eq!(peak_norm(&[C64::new(f64::NAN, 0.0)]), None);
        assert_eq!(peak_norm(&[C64::new(0.0, f64::INFINITY)]), None);
    }
}
