//! Focused G0/G1/G3 tests of sampled geometry and its real quadrature consumer.
use super::*;
use crate::HexCell;
use crate::quad3::{QuadratureControl3, QuadratureOptions3, cut_cell_rules3};

fn field(shape: [usize; 3], origin: [f64; 3], spacing: [f64; 3],
    f: impl Fn([f64; 3]) -> f64) -> SampledSdf3 {
    let mut values = Vec::new();
    for z in 0..shape[2] {
        for y in 0..shape[1] {
            for x in 0..shape[0] {
                let index = [x, y, z];
                values.push(f(std::array::from_fn(|a| origin[a] + index[a] as f64 * spacing[a])));
            }
        }
    }
    SampledSdf3::with_control(origin, spacing, shape, values,
        &mut || ControlFlow::Continue(())).unwrap()
}
fn contains(range: Interval, value: f64) {
    assert!(range.lo() <= value && value <= range.hi(), "{value} outside {range:?}");
}

#[test]
fn g0_x_fastest_knots_and_translated_anisotropic_affine_field() {
    let f = |p: [f64; 3]| 2.0 * p[0] - 3.0 * p[1] + 4.0 * p[2] - 0.375;
    let sampled = field([4, 3, 5], [10.0, -2.0, 0.0], [0.25, 0.5, 0.125], f);
    assert_eq!(sampled.bounds(), ([10.0, -2.0, 0.0], [10.75, -1.0, 0.5]));
    for z in 0..5 {
        for y in 0..3 {
            for x in 0..4 {
                let p = [10.0 + x as f64 * 0.25, -2.0 + y as f64 * 0.5, z as f64 * 0.125];
                assert_eq!(sampled.value(p), f(p));
            }
        }
    }
    let lo = [10.1, -1.9, 0.05];
    let hi = [10.6, -1.1, 0.4];
    for (axis, derivative) in [(HeightAxis::X, 2.0), (HeightAxis::Y, -3.0), (HeightAxis::Z, 4.0)] {
        let range = sampled.derivative_enclose(lo, hi, axis);
        contains(range, derivative);
        assert!(range.hi() - range.lo() < 1e-10);
    }
    let p = [10.375, -1.25, 0.3125];
    contains(sampled.enclose(p, p), f(p));
    assert!((sampled.value(p) - f(p)).abs() < 1e-12);
}

#[test]
fn g3_multiaffine_ranges_and_derivatives_cover_interior_values() {
    let f = |p: [f64; 3]| p[0] * p[1] - 2.0 * p[1] * p[2]
        + 3.0 * p[0] * p[2] + p[0] * p[1] * p[2] - 0.25;
    let sampled = field([5; 3], [0.0; 3], [0.25; 3], f);
    let lo = [0.125, 0.2, 0.3];
    let hi = [0.875, 0.7, 0.8];
    let ranges = [sampled.enclose(lo, hi),
        sampled.derivative_enclose(lo, hi, HeightAxis::X),
        sampled.derivative_enclose(lo, hi, HeightAxis::Y),
        sampled.derivative_enclose(lo, hi, HeightAxis::Z)];
    for i in 0..1000 {
        let t = [i % 97, (i * 13) % 97, (i * 37) % 97];
        let p: [f64; 3] = std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * t[a] as f64 / 96.0);
        contains(ranges[0], f(p));
        contains(ranges[1], p[1] + 3.0 * p[2] + p[1] * p[2]);
        contains(ranges[2], p[0] - 2.0 * p[2] + p[0] * p[2]);
        contains(ranges[3], -2.0 * p[1] + 3.0 * p[0] + p[0] * p[1]);
        assert!((sampled.value(p) - f(p)).abs() < 1e-12);
    }
}


#[test]
fn g0_surface_aligned_subbox_tightens_roundoff_without_hiding_a_real_sliver() {
    let sampled = field([3; 3], [0.0; 3], [0.5; 3], |p| {
        let roof = 0.625 + if p[0] == 0.5 && p[1] == 0.5 { 0.125 } else { 0.0 };
        p[2] - roof
    });
    let lo = [0.25, 0.25, 0.75];
    let hi = [0.5, 0.5, 1.0];
    let range = sampled.enclose(lo, hi);
    assert_eq!(range.lo(), 0.0, "surface contact must not invent negative material: {range:?}");
    assert!(range.hi() > 0.0);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let rules = cut_cell_rules3(&sampled, HexCell::try_new(lo, hi).unwrap(), &mut control).unwrap();
    assert!(rules.bulk().is_empty());
    assert_eq!(rules.cut_boxes(), 0);

    let sliver = field([2; 3], [0.0; 3], [1.0; 3], |p| p[2] - 0.750_000_1);
    let range = sliver.enclose([0.0, 0.0, 0.75], [1.0, 1.0, 1.0]);
    assert!(range.lo() < 0.0, "exact sign tightening must retain a genuine negative sliver: {range:?}");
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let rules = cut_cell_rules3(&sliver,
        HexCell::try_new([0.0, 0.0, 0.75], [1.0, 1.0, 1.0]).unwrap(), &mut control).unwrap();
    assert!(rules.volume() > 0.0);
}

#[test]
fn g0_seam_queries_include_both_one_sided_slopes() {
    let sampled = field([3, 2, 2], [0.0; 3], [1.0; 3], |p| {
        if p[0] <= 1.0 { 2.0 - 2.0 * p[0] } else { 3.0 * (p[0] - 1.0) }
    });
    let p = [1.0, 0.5, 0.5];
    let seam = sampled.derivative_enclose(p, p, HeightAxis::X);
    contains(seam, -2.0);
    contains(seam, 3.0);
    let right = sampled.derivative_enclose([1.25, 0.2, 0.2], [1.5, 0.8, 0.8], HeightAxis::X);
    contains(right, 3.0);
    assert!(right.lo() > 2.99);
    contains(sampled.enclose(p, p), 0.0);
}

#[test]
fn g0_interior_samples_cannot_be_hidden_by_outer_corner_signs() {
    let sampled = field([3, 2, 2], [0.0; 3], [1.0; 3], |p| (p[0] - 1.0).abs() * 2.0 - 1.0);
    let bounds = sampled.enclose([0.0; 3], [2.0, 1.0, 1.0]);
    contains(bounds, -1.0);
    contains(bounds, 1.0);
    assert_eq!(sampled.value([1.0, 0.5, 0.5]), -1.0);
}

#[test]
fn g0_large_box_uses_sound_global_value_and_slope_bounds() {
    let sampled = field([10; 3], [0.0; 3], [1.0; 3], |p| p[0] + 2.0 * p[1] + 3.0 * p[2]);
    let range = sampled.enclose([0.0; 3], [9.0; 3]);
    contains(range, 0.0);
    contains(range, 54.0);
    for (axis, slope) in [(HeightAxis::X, 1.0), (HeightAxis::Y, 2.0), (HeightAxis::Z, 3.0)] {
        let range = sampled.derivative_enclose([0.0; 3], [9.0; 3], axis);
        contains(range, slope);
        assert!(range.lo() > 0.0);
    }
}

#[test]
fn g0_no_extrapolation_or_silent_sample_repair() {
    let sampled = field([2; 3], [0.0; 3], [1.0; 3], |p| p[2] - 0.5);
    assert!(sampled.value([-0.01, 0.5, 0.5]).is_nan());
    assert!(sampled.value([f64::NAN, 0.5, 0.5]).is_nan());
    for (lo, hi) in [([-0.01, 0.0, 0.0], [1.0; 3]), ([0.0; 3], [1.01; 3]), ([1.0; 3], [0.0; 3])] {
        let range = sampled.enclose(lo, hi);
        assert_eq!(range.lo(), f64::NEG_INFINITY);
        assert_eq!(range.hi(), f64::INFINITY);
    }
    assert_eq!(sampled.value([1.0; 3]), 0.5);
    let make = |origin, spacing, shape, values| SampledSdf3::with_control(origin, spacing, shape, values,
        &mut || ControlFlow::Continue(()));
    assert!(make([0.0; 3], [1.0; 3], [2; 3], vec![0.0; 7]).is_err());
    assert!(make([0.0; 3], [1.0; 3], [usize::MAX; 3], vec![]).is_err());
    assert!(make([0.0; 3], [0.0; 3], [2; 3], vec![0.0; 8]).is_err());
    assert!(make([0.0; 3], [-1.0; 3], [2; 3], vec![0.0; 8]).is_err());
    assert!(make([1.0; 3], [1e-100; 3], [2; 3], vec![0.0; 8]).is_err());
    assert!(make([0.0; 3], [1.0; 3], [2; 3], vec![f64::NAN; 8]).is_err());
    assert!(make([0.0; 3], [1.0; 3], [2; 3], vec![1e101; 8]).is_err());
}

#[test]
fn g4_construction_checks_cancellation_during_the_sample_sweep() {
    let mut calls = 0;
    let result = SampledSdf3::with_control([0.0; 3], [1.0; 3], [9; 3], vec![0.0; 729], &mut || {
        calls += 1;
        if calls == 7 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    });
    assert!(matches!(result, Err(SampledSdf3Error::Cancelled)));
    assert_eq!(calls, 7);
}

#[test]
fn g1_real_cut_quadrature_integrates_sampled_plane_volume() {
    let sampled = field([5; 3], [0.0; 3], [0.25; 3], |p| p[2] - 0.375);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let rules = cut_cell_rules3(&sampled, HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), &mut control).unwrap();
    assert!((rules.volume() - 0.375).abs() < 1e-8, "volume={}", rules.volume());
    contains(rules.volume_bounds(), 0.375);
    assert!(rules.bulk().iter().all(|(_, weight)| *weight > 0.0));
}
