//! G0/G1/G3 range, sign, ray-invariance and manufactured Newton regressions.
use super::*;
use super::super::{TrustRegionState, TrustRegionProgress, StopRule, StopReason};

#[test]
fn binary_fallback_scaling_roundtrips_every_exponent_bin_and_signed_zero() {
    let mask = (1u64 << 52) - 1;
    for exponent in 0u64..=2046 {
        for fraction in [0, 1, mask >> 1, mask] {
            for sign in [0, 1u64 << 63] {
                let x = f64::from_bits(sign | (exponent << 52) | fraction);
                let (mantissa, power) = split(x);
                assert_eq!(scale(mantissa, power).to_bits(), x.to_bits());
            }
        }
    }
    assert_eq!(scale(0.0, 4096).to_bits(), 0.0f64.to_bits());
    assert_eq!(scale(-0.0, 4096).to_bits(), (-0.0f64).to_bits());
}

#[test]
fn overflowing_curvature_keeps_a_representable_positive_newton_ratio() {
    let (negative, alpha) = curvature(2.0, &[1.0, 1.0], &[1e308, 1e308]);
    assert!(!negative);
    assert!((alpha / 1e-308 - 1.0).abs() < 1e-14);
    // Individually overflowing terms cancel; the remaining curvature is positive.
    let (negative, alpha) = curvature(9.0, &[2.0, 2.0, 1.0], &[1e308, -1e308, 1e308]);
    assert!(!negative);
    assert!((alpha / 9e-308 - 1.0).abs() < 1e-14);
}

#[test]
fn underflowing_curvature_is_not_mislabelled_as_negative() {
    let (negative, alpha) = curvature(1e-300, &[1e-200], &[1e-200]);
    assert!(!negative);
    assert!((alpha / 1e100 - 1.0).abs() < 1e-14);
    assert!(curvature(1e-300, &[1e-200], &[-1e-200]).0);
    assert!(curvature(1.0, &[1.0], &[0.0]).0);
}

#[test]
fn ordinary_curvature_preserves_its_operation_tree() {
    let d = [0.25, -0.125, 0.0625];
    let h = [2.0, 0.75, -0.5];
    let expected = 3.0 / d.iter().zip(&h).map(|(a, b)| a*b).sum::<f64>();
    let (negative, actual) = curvature(3.0, &d, &h);
    assert!(!negative);
    assert_eq!(actual.to_bits(), expected.to_bits());
}

#[test]
fn boundary_point_survives_extreme_radius_and_ray_units() {
    for radius in [1e-250, 1.0, 1e250] {
        for direction_scale in [f64::from_bits(16), 1e-200, 1.0, 1e200, 1e300] {
            let mut p = [0.0, 0.0];
            let d = [3.0*direction_scale, 4.0*direction_scale];
            assert!(to_boundary(&mut p, &d, radius));
            assert!((p[0] / radius - 0.6).abs() < 4e-15, "{radius} {direction_scale}: {p:?}");
            assert!((p[1] / radius - 0.8).abs() < 4e-15);
        }
    }
}

#[test]
fn boundary_fallback_preserves_the_off_origin_ray() {
    for radius in [1e-250, 1.0, 1e250] {
        for direction_scale in [1e-200, 1e200] {
            let mut p = [0.25*radius, 0.5*radius];
            assert!(to_boundary(&mut p, &[direction_scale, 0.0], radius));
            assert!((p[0]/radius - 0.75f64.sqrt()).abs() < 4e-15);
            assert_eq!(p[1], 0.5*radius);
        }
    }
}

#[test]
fn near_boundary_outward_step_retains_the_representable_increment() {
    let mut p = [0.6, 0.7999999999999999];
    assert!(to_boundary(&mut p, &[1.0, 0.0], 1.0));
    // sqrt(1-y^2), independently rounded from the exact input float using
    // 120-digit decimal arithmetic. No re-use of this implementation's root.
    assert_eq!(p[0].to_bits(), 0.6000000000000001f64.to_bits());
    assert_eq!(p[1], 0.7999999999999999);

    // The old pp reduction discards this tiny x^2 contribution, inflating the
    // horizontal crossing by about 0.2%. The exact-input decimal oracle gives
    // 2^-26 as the correctly rounded x coordinate.
    let mut tangent = [1e-9, 0.9999999999999999];
    assert!(to_boundary(&mut tangent, &[1.0, 0.0], 1.0));
    assert_eq!(tangent[0].to_bits(), (1.0 / 67_108_864.0f64).to_bits());
    assert_eq!(tangent[1], 0.9999999999999999);
}

#[test]
fn invalid_boundary_request_cannot_mutate_the_retained_step() {
    for (p, d, radius) in [
        ([0.5, 0.0], [0.0, 0.0], 1.0),
        ([2.0, 0.0], [1.0, 0.0], 1.0),
        ([0.5, 0.0], [f64::INFINITY, 0.0], 1.0),
        ([0.5, 0.0], [1.0, 0.0], f64::INFINITY),
    ] {
        let mut result = p;
        assert!(!to_boundary(&mut result, &d, radius));
        assert_eq!(result, p);
    }
}

#[test]
fn giant_finite_stiffness_solves_and_restarts_under_a_hessian_ceiling() {
    let mut fg = |x: &[f64]| {
        (x.iter().map(|v| v + (0.5e308 * v) * v).sum(),
         x.iter().map(|v| 1.0 + 1e308 * v).collect())
    };
    let mut hv = |_: &[f64], v: &[f64]| v.iter().map(|v| 1e308 * v).collect();
    let rule = StopRule::GradNorm(1e-14);
    let initial = TrustRegionState::new(&[0.0, 0.0], &mut fg);
    let mut direct = initial.clone();
    let result = direct.run(&mut fg, &mut hv, &rule, 3);
    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert!(result.solution.f < 0.0 && result.solution.f.is_finite());
    for x in result.solution.x { assert!((x / -1e-308 - 1.0).abs() < 1e-14); }
    let mut resumed = initial;
    let result = resumed.run_bounded(&mut fg, &mut hv, &rule, 3, 1);
    assert_eq!(result.progress, TrustRegionProgress::HessianBudget);
    assert_eq!(result.solution.x, vec![0.0; 2]);
    assert_eq!(resumed.radius(), 1.0);
    let result = resumed.run_bounded(&mut fg, &mut hv, &rule, 3, 3);
    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert_eq!(resumed.report().x, direct.report().x);
    assert_eq!(resumed.history(), direct.history());
    assert_eq!(resumed.report().hv_evals, direct.report().hv_evals + 1);
}
