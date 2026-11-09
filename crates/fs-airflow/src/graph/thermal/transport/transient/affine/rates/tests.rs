use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_swept_rates_match_independent_signed_flux_quadrature() {
    for [initial, final_flow] in [
        [2.0, -2.0],
        [-3.0, 1.0],
        [1.0, -3.0],
        [0.0, 4.0],
        [4.0, 0.0],
        [-3.0, -1.0],
    ] {
        let mean = 0.5 * initial + 0.5 * final_flow;
        let actual = AffineVolumeRates::try_new(initial, mean, final_flow).unwrap();
        // Composite trapezoids on the actual affine trajectory. The selected
        // crossings coincide with grid points, so both piecewise-linear
        // positive/negative flux integrals are exact apart from roundoff.
        let mut integrals = [0.0; 2];
        let n = 4096;
        for i in 0..=n {
            let fraction = i as f64 / n as f64;
            let q = (1.0 - fraction) * initial + fraction * final_flow;
            let weight = if i == 0 || i == n { 0.5 } else { 1.0 };
            integrals[0] += weight * q.max(0.0) / n as f64;
            integrals[1] += weight * (-q).max(0.0) / n as f64;
        }
        near(actual.forward_m3_s, integrals[0], 2e-14);
        near(actual.reverse_m3_s, integrals[1], 2e-14);
        near(actual.forward_m3_s - actual.reverse_m3_s, mean, 2e-15);
    }
    let reversal = AffineVolumeRates::try_new(2.0, 0.0, -2.0).unwrap();
    assert_eq!(reversal.forward_m3_s, 0.5);
    assert_eq!(reversal.reverse_m3_s, 0.5);
}

#[test]
fn g3_direction_relabeling_time_partition_and_no_reversal_parity() {
    for (a, b) in [(7.0, -2.0), (-9.0, 4.0), (2.0, 7.0), (0.0, 0.0)] {
        let mean = 0.5 * a + 0.5 * b;
        let whole = AffineVolumeRates::try_new(a, mean, b).unwrap();
        let reverse = AffineVolumeRates::try_new(-a, -mean, -b).unwrap();
        assert_eq!(whole.forward_m3_s, reverse.reverse_m3_s);
        assert_eq!(whole.reverse_m3_s, reverse.forward_m3_s);
        for fraction in [0.125, 0.25, 0.5, 0.875] {
            let split = (1.0 - fraction) * a + fraction * b;
            let first = AffineVolumeRates::try_new(a, 0.5 * a + 0.5 * split, split).unwrap();
            let second = AffineVolumeRates::try_new(split, 0.5 * split + 0.5 * b, b).unwrap();
            near(
                whole.forward_m3_s,
                fraction * first.forward_m3_s + (1.0 - fraction) * second.forward_m3_s,
                2e-15,
            );
            near(
                whole.reverse_m3_s,
                fraction * first.reverse_m3_s + (1.0 - fraction) * second.reverse_m3_s,
                2e-15,
            );
        }
    }
    let mean = f64::from_bits(1.0f64.to_bits() + 1);
    let unchanged = AffineVolumeRates::try_new(0.5, mean, 1.5).unwrap();
    assert_eq!(unchanged.forward_m3_s.to_bits(), mean.to_bits());
    assert_eq!(unchanged.reverse_m3_s, 0.0);
}

#[test]
fn g0_scaled_quadrature_preserves_extreme_representable_swept_rates() {
    for scale in [1e-290, 1.0, 1e290, f64::MAX] {
        let actual = AffineVolumeRates::try_new(scale, 0.0, -scale).unwrap();
        near(actual.forward_m3_s / scale, 0.25, 4e-16);
        near(actual.reverse_m3_s / scale, 0.25, 4e-16);
    }
    let steady = AffineVolumeRates::try_new(f64::MAX, f64::MAX, f64::MAX).unwrap();
    assert_eq!(steady.forward_m3_s, f64::MAX);
    assert_eq!(steady.reverse_m3_s, 0.0);
}

#[test]
fn g0_unsupported_mean_nonfinite_and_unrepresentable_lobes_refuse() {
    for [a, mean, b] in [
        [1.0, 0.0, 1.0],
        [1.0, 1.0, -1.0],
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0],
        [1e-300, -0.5e300, -1e300],
        [f64::from_bits(1), 0.0, 0.0],
    ] {
        assert!(AffineVolumeRates::try_new(a, mean, b).is_err());
    }
}
