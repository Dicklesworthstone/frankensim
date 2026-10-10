use super::*;
use crate::{HexCell, VerticalLineRoot, isolate_certified_height_root};

fn contains(range: Interval, value: f64) {
    assert!(range.lo() <= value && value <= range.hi(), "{value} outside [{}, {}]", range.lo(), range.hi());
}

#[test]
fn quadrics_enclose_values_and_analytic_partials_on_translated_anisotropic_boxes() {
    let mut builder = CsgBuilder3::new();
    let center = [0.25, -0.4, 0.6];
    let radii = [0.5, 0.8, 0.3];
    let root = builder.ellipsoid(center, radii).unwrap();
    let domain = builder.finish(root).unwrap();
    let lo = [-0.2, -0.6, 0.1];
    let hi = [0.8, 0.7, 0.9];
    let value_range = domain.enclose(lo, hi);
    let axes = [HeightAxis::X, HeightAxis::Y, HeightAxis::Z];
    for i in 0..=6 {
        for j in 0..=6 {
            for k in 0..=6 {
                let sample = [i, j, k];
                let p = std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * sample[a] as f64 / 6.0);
                contains(value_range, domain.value(p));
                contains(domain.enclose(p, p), domain.value(p));
                for a in 0..3 {
                    let derivative = 0.3 * ((p[a] - center[a]) / radii[a]) / radii[a];
                    contains(domain.derivative_enclose(lo, hi, axes[a]), derivative);
                }
            }
        }
    }
    assert!(domain.value(center) < 0.0);
    assert!(domain.value([2.0, 0.0, 0.0]) > 0.0);
}

#[test]
fn all_cylinder_axes_have_zero_axial_slope_and_the_declared_cross_section() {
    for (axis, index) in [(HeightAxis::X, 0), (HeightAxis::Y, 1), (HeightAxis::Z, 2)] {
        let mut builder = CsgBuilder3::new();
        let node = builder.cylinder(axis, [0.0; 3], 0.25).unwrap();
        let domain = builder.finish(node).unwrap();
        let mut p = [0.0; 3];
        p[index] = 100.0;
        assert_eq!(domain.value(p), -0.125);
        let slope = domain.derivative_enclose([-1.0; 3], [1.0; 3], axis);
        assert_eq!(slope.lo(), 0.0);
        assert_eq!(slope.hi(), 0.0);
        p[(index + 1) % 3] = 0.5;
        assert!(domain.value(p) > 0.0);
    }
}

#[test]
fn bored_plate_reaches_the_existing_certified_height_line_consumer() {
    let mut builder = CsgBuilder3::new();
    let plate = builder.half_space([0.0, 0.0, 1.0], 0.63).unwrap();
    let hole = builder.cylinder(HeightAxis::Y, [0.55, 0.0, 0.28], 0.12).unwrap();
    let root = builder.combine(CsgOp3::Difference, plate, hole, 0.0).unwrap();
    let domain = builder.finish(root).unwrap();
    assert!(domain.value([0.55, 0.2, 0.28]) > 0.0, "the bore removes material");
    assert!(domain.value([0.1, 0.2, 0.28]) < 0.0);
    assert!(domain.value([0.1, 0.2, 0.9]) > 0.0);
    let cell = HexCell::try_new([0.6, 0.1, 0.29], [0.72, 0.3, 0.31]).unwrap();
    let result = isolate_certified_height_root(&domain, cell, [0.2, 0.3], 8).unwrap();
    let VerticalLineRoot::CertifiedRoot { height_axis, enclosure } = result else {
        panic!("a separated cylindrical side must be a genuine certified crossing: {result:?}");
    };
    assert_eq!(height_axis, HeightAxis::X);
    contains(enclosure, 0.55 + (0.12_f64.powi(2) - (0.3_f64 - 0.28).powi(2)).sqrt());
}

#[test]
fn smooth_boolean_derivative_uses_blend_weights_not_a_selected_child() {
    let mut builder = CsgBuilder3::new();
    let x = builder.half_space([1.0, 0.0, 0.0], 0.0).unwrap();
    let y = builder.half_space([0.0, 1.0, 0.0], 0.0).unwrap();
    let blended = builder.combine(CsgOp3::Union, x, y, 0.2).unwrap();
    let hard = builder.combine(CsgOp3::Union, x, y, 0.0).unwrap();
    let hard_domain = CsgDomain3 { nodes: builder.nodes.clone(), root: hard.0 };
    let domain = builder.finish(blended).unwrap();
    let p = [0.05, 0.1, 0.0];
    contains(domain.derivative_enclose(p, p, HeightAxis::X), 0.625);
    contains(domain.derivative_enclose(p, p, HeightAxis::Y), 0.375);
    assert!(hard_domain.value([0.01, 0.01, 0.0]) > 0.0);
    assert!(domain.value([0.01, 0.01, 0.0]) < 0.0, "blend is an explicit geometry change");
    for op in [CsgOp3::Union, CsgOp3::Intersection, CsgOp3::Difference] {
        let mut b = CsgBuilder3::new();
        let a = b.sphere([0.2, 0.0, 0.0], 0.5).unwrap();
        let c = b.box_region([0.0; 3], [0.4, 0.3, 0.2]).unwrap();
        let root = b.combine(op, a, c, 0.1).unwrap();
        let d = b.finish(root).unwrap();
        let lo = [-0.3, -0.2, -0.1];
        let hi = [0.6, 0.4, 0.5];
        let range = d.enclose(lo, hi);
        for i in 0..=20 {
            let p = std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * i as f64 / 20.0);
            contains(range, d.value(p));
            contains(d.enclose(p, p), d.value(p));
        }
    }
}

#[test]
fn hard_crease_does_not_invent_a_monotone_height_axis() {
    let mut b = CsgBuilder3::new();
    let a = b.half_space([1.0, 0.0, 0.0], 0.0).unwrap();
    let c = b.half_space([-1.0, 0.0, 0.0], 0.0).unwrap();
    let root = b.combine(CsgOp3::Intersection, a, c, 0.0).unwrap();
    let d = b.finish(root).unwrap();
    let cell = HexCell::try_new([-0.2; 3], [0.2; 3]).unwrap();
    let slope = d.derivative_enclose(cell.lo(), cell.hi(), HeightAxis::X);
    contains(slope, -1.0);
    contains(slope, 1.0);
    assert_eq!(isolate_certified_height_root(&d, cell, [0.0, 0.0], 8).unwrap(), VerticalLineRoot::Ambiguous);
}

#[test]
fn malformed_construction_and_work_overflow_refuse_without_losing_the_recipe() {
    let mut b = CsgBuilder3::new();
    assert!(b.sphere([f64::NAN, 0.0, 0.0], 1.0).is_err());
    assert!(b.ellipsoid([0.0; 3], [1.0, -1.0, 1.0]).is_err());
    assert!(b.half_space([0.0; 3], 1.0).is_err());
    let mut root = b.sphere([0.0; 3], 1.0).unwrap();
    assert!(b.combine(CsgOp3::Union, root, root, f64::NAN).is_err());
    for _ in 1..MAX_CSG3_NODES { root = b.combine(CsgOp3::Union, root, root, 0.0).unwrap(); }
    assert!(b.combine(CsgOp3::Union, root, root, 0.0).is_err());
    let d = b.finish(root).unwrap();
    assert_eq!(d.node_count(), MAX_CSG3_NODES);
    assert_eq!(d.value([0.0; 3]), -0.5);
    assert!(d.value([f64::INFINITY; 3]).is_nan());
    let invalid = d.enclose([1.0; 3], [-1.0; 3]);
    assert_eq!(invalid.lo(), f64::NEG_INFINITY);
    assert_eq!(invalid.hi(), f64::INFINITY);
}
