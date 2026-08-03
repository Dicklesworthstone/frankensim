//! G0/G3 contact-support checks for axisymmetric line/arc charts.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Vec3;
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricSupportAuthority, AxisymmetricSupportError, MeridianPoint,
    MeridianSegment, SquatDiscEdgeTreatment,
};

fn point(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}

fn line(start: MeridianPoint, end: MeridianPoint) -> MeridianSegment {
    MeridianSegment::Line { start, end }
}

fn cylinder(radius: f64, half_height: f64, axial_offset: f64) -> AxisymmetricChart {
    AxisymmetricChart::try_new(vec![
        line(
            point(0.0, axial_offset - half_height),
            point(radius, axial_offset - half_height),
        ),
        line(
            point(radius, axial_offset - half_height),
            point(radius, axial_offset + half_height),
        ),
        line(
            point(radius, axial_offset + half_height),
            point(0.0, axial_offset + half_height),
        ),
        line(
            point(0.0, axial_offset + half_height),
            point(0.0, axial_offset - half_height),
        ),
    ])
    .expect("closed CCW cylinder")
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x53555050,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g0_sharp_cylinder_support_selects_the_closed_form_rim_point() {
    let chart = cylinder(2.0, 1.0, 0.0);
    let support = with_cx(false, |cx| {
        chart
            .minimum_support_point(Vec3::new(3.0e200, 4.0e200, -5.0e200), cx)
            .expect("tipped direction has a unique sharp rim")
    });

    assert_close(support.point.x, -1.2);
    assert_close(support.point.y, -1.6);
    assert_close(support.point.z, 1.0);
    assert_close(support.support_value, -3.0 / 2.0_f64.sqrt());
    assert_eq!(support.source_feature, 1);
    assert_eq!(support.authority, AxisymmetricSupportAuthority::Estimate);
}

#[test]
fn g0_true_fillet_uses_an_arc_interior_extremum_not_a_polygon_corner() {
    let chart = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
    )
    .expect("true circular fillets");
    let support = with_cx(false, |cx| {
        chart
            .minimum_support_point(Vec3::new(1.0, 0.0, -1.0), cx)
            .expect("tipped direction selects upper fillet interior")
    });
    let root_half = 0.5_f64.sqrt() * 0.5;

    assert_close(support.point.x, -(1.5 + root_half));
    assert_close(support.point.y, 0.0);
    assert_close(support.point.z, 0.5 + root_half);
    assert_eq!(support.source_feature, 3);
}

#[test]
fn g3_axial_translation_and_z_rotation_transform_support_equivariantly() {
    let base = cylinder(2.0, 1.0, 0.0);
    let shifted = cylinder(2.0, 1.0, 7.0);
    let direction = Vec3::new(3.0, 4.0, -5.0);
    let angle: f64 = 0.713;
    let rotated_direction = Vec3::new(
        direction.x * angle.cos() - direction.y * angle.sin(),
        direction.x * angle.sin() + direction.y * angle.cos(),
        direction.z,
    );
    let (base_support, shifted_support, rotated_support) = with_cx(false, |cx| {
        (
            base.minimum_support_point(direction, cx)
                .expect("base support"),
            shifted
                .minimum_support_point(direction, cx)
                .expect("translated support"),
            base.minimum_support_point(rotated_direction, cx)
                .expect("rotated support"),
        )
    });

    assert_close(shifted_support.point.x, base_support.point.x);
    assert_close(shifted_support.point.y, base_support.point.y);
    assert_close(shifted_support.point.z, base_support.point.z + 7.0);
    assert_close(
        shifted_support.support_value,
        base_support.support_value + 7.0 * direction.z / 50.0_f64.sqrt(),
    );
    assert_close(
        rotated_support.point.x,
        base_support.point.x * angle.cos() - base_support.point.y * angle.sin(),
    );
    assert_close(
        rotated_support.point.y,
        base_support.point.x * angle.sin() + base_support.point.y * angle.cos(),
    );
    assert_close(rotated_support.point.z, base_support.point.z);
    assert_close(rotated_support.support_value, base_support.support_value);
}

#[test]
fn g0_sharp_cylinder_radial_and_cap_normal_support_are_non_unique() {
    let chart = cylinder(2.0, 1.0, 0.0);
    with_cx(false, |cx| {
        assert!(matches!(
            chart.minimum_support_point(Vec3::new(1.0, 0.0, 0.0), cx),
            Err(AxisymmetricSupportError::NonUniqueFeatureSupport { .. })
        ));
        assert!(matches!(
            chart.minimum_support_point(Vec3::new(0.0, 0.0, 1.0), cx),
            Err(AxisymmetricSupportError::NonUniqueFeatureSupport { .. })
        ));
    });
}

#[test]
fn g0_hostile_and_cancelled_support_requests_refuse() {
    let chart = cylinder(2.0, 1.0, 0.0);
    with_cx(false, |cx| {
        assert!(matches!(
            chart.minimum_support_point(Vec3::new(0.0, 0.0, 0.0), cx),
            Err(AxisymmetricSupportError::ZeroDirection)
        ));
        assert!(matches!(
            chart.minimum_support_point(Vec3::new(f64::NAN, 0.0, 1.0), cx),
            Err(AxisymmetricSupportError::NonFiniteDirection { .. })
        ));
    });
    with_cx(true, |cx| {
        assert!(matches!(
            chart.minimum_support_point(Vec3::new(1.0, 0.0, -1.0), cx),
            Err(AxisymmetricSupportError::Cancelled)
        ));
    });
}
