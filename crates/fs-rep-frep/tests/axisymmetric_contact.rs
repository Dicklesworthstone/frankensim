//! G0/G3 checks for the validated axisymmetric line/arc chart.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{BettiBounds, Chart, Point3, TraceStepClaim};
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricError, MAX_AXISYMMETRIC_SEGMENTS, MeridianPoint,
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
    .expect("closed CCW cylinder meridian")
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xA815,
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

#[test]
fn g0_cylinder_axis_inside_outside_and_non_authoritative_distance_estimate() {
    let chart = cylinder(2.0, 1.0, 0.0);
    with_cx(|cx| {
        let inside = chart.eval(Point3::new(1.0, 0.0, 0.0), cx);
        let side = chart.eval(Point3::new(3.0, 0.0, 0.0), cx);
        let cap = chart.eval(Point3::new(0.0, 0.0, 2.0), cx);
        let axis = chart.eval(Point3::new(0.0, 0.0, 0.0), cx);
        let seam = chart.eval(Point3::new(3.0, 0.0, 2.0), cx);
        assert_eq!(inside.signed_distance, -1.0);
        assert_eq!(side.signed_distance, 1.0);
        assert_eq!(cap.signed_distance, 1.0);
        assert_eq!(axis.signed_distance, -1.0);
        assert!(axis.gradient.is_none(), "axis closest points form a circle");
        assert!(
            seam.gradient.is_none(),
            "equal line-feature minima are set-valued"
        );
        assert_eq!(side.error.kind, fs_evidence::NumericalKind::Estimate);
        assert_eq!(side.error.lo, side.signed_distance);
        assert_eq!(side.error.hi, side.signed_distance);
    });
    assert_eq!(chart.trace_step_claim(), TraceStepClaim::NoClaim);
    assert_eq!(chart.construction_certificate().surfaced_feature_count, 3);
    assert!(chart.verify_construction().is_ok());
}

#[test]
fn g0_cancelled_search_refuses_without_promoting_a_partial_feature_scan() {
    let chart = cylinder(2.0, 1.0, 0.0);
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xA815,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let sample = chart.eval(Point3::new(3.0, 0.0, 0.0), &cx);
        assert!(sample.signed_distance.is_nan());
        assert!(sample.lipschitz.is_none());
        assert_eq!(sample.error.kind, fs_evidence::NumericalKind::NoClaim);
    });
}

#[test]
fn g0_squat_disc_sharp_and_zero_fillet_match_closed_form_samples() {
    let sharp = AxisymmetricChart::squat_disc(2.0, 2.0, SquatDiscEdgeTreatment::Sharp)
        .expect("positive sharp disc");
    let zero = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.0 },
    )
    .expect("zero fillet canonicalizes to sharp");
    assert_eq!(sharp.segments(), zero.segments());
    with_cx(|cx| {
        assert_eq!(
            sharp.eval(Point3::new(0.0, 0.0, 0.0), cx).signed_distance,
            -1.0
        );
        assert_eq!(
            sharp.eval(Point3::new(3.0, 0.0, 0.0), cx).signed_distance,
            1.0
        );
        assert_eq!(
            sharp.eval(Point3::new(0.0, 0.0, 2.0), cx).signed_distance,
            1.0
        );
        assert_eq!(
            sharp.eval(Point3::new(3.0, 0.0, 2.0), cx).signed_distance,
            2.0_f64.sqrt()
        );
    });
}

#[test]
fn g0_squat_disc_true_fillets_handle_maximum_radius_without_degenerate_lines() {
    let filleted = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 1.0 },
    )
    .expect("edge radius equals thickness/2");
    assert_eq!(filleted.construction_certificate().input_feature_count, 5);
    assert_eq!(
        filleted.construction_certificate().surfaced_feature_count,
        4
    );
    with_cx(|cx| {
        let cylindrical = filleted.eval(Point3::new(2.25, 0.0, 0.0), cx);
        assert_eq!(cylindrical.signed_distance, 0.25);
        let corner = filleted.eval(Point3::new(2.5, 0.0, 1.0), cx);
        assert!((corner.signed_distance - (3.25_f64.sqrt() - 1.0)).abs() < 1e-12);
    });

    let all_fillet = AxisymmetricChart::squat_disc(
        1.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 1.0 },
    )
    .expect("edge radius equals both closed bounds");
    assert_eq!(all_fillet.construction_certificate().input_feature_count, 3);
}

#[test]
fn g0_squat_disc_rejects_invalid_dimensions_and_edge_radii() {
    for (radius, thickness) in [(0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, f64::NAN)] {
        assert!(matches!(
            AxisymmetricChart::squat_disc(radius, thickness, SquatDiscEdgeTreatment::Sharp),
            Err(AxisymmetricError::NonPositiveDimension { .. })
        ));
    }
    for radius in [-0.1, 1.01, f64::NAN] {
        assert!(matches!(
            AxisymmetricChart::squat_disc(
                2.0,
                2.0,
                SquatDiscEdgeTreatment::CircularFillet { radius },
            ),
            Err(AxisymmetricError::InvalidEdgeRadius { .. })
        ));
    }
}

#[test]
fn g0_true_circular_fillets_and_sloped_chamfers_are_not_polygonized() {
    let filleted = AxisymmetricChart::try_new(vec![
        line(point(0.0, -1.0), point(1.0, -1.0)),
        MeridianSegment::Arc {
            start: point(1.0, -1.0),
            end: point(1.5, -0.5),
            center: point(1.0, -0.5),
            clockwise: false,
        },
        line(point(1.5, -0.5), point(1.5, 0.5)),
        MeridianSegment::Arc {
            start: point(1.5, 0.5),
            end: point(1.0, 1.0),
            center: point(1.0, 0.5),
            clockwise: false,
        },
        line(point(1.0, 1.0), point(0.0, 1.0)),
        line(point(0.0, 1.0), point(0.0, -1.0)),
    ])
    .expect("two quarter-circle fillets bound a simple profile");
    let chamfered = AxisymmetricChart::try_new(vec![
        line(point(0.0, -1.0), point(1.0, -1.0)),
        line(point(1.0, -1.0), point(1.4, -0.6)),
        line(point(1.4, -0.6), point(1.4, 0.6)),
        line(point(1.4, 0.6), point(1.0, 1.0)),
        line(point(1.0, 1.0), point(0.0, 1.0)),
        line(point(0.0, 1.0), point(0.0, -1.0)),
    ])
    .expect("sloped line faces form admissible chamfers");
    with_cx(|cx| {
        let rounded = filleted.eval(Point3::new(1.5, 0.0, -1.0), cx);
        assert!(
            (rounded.signed_distance - (0.5_f64.sqrt() - 0.5)).abs() < 1e-12,
            "true quarter-circle distance"
        );
        let chamfer = chamfered.eval(Point3::new(1.4, 0.0, -1.0), cx);
        assert!(chamfer.signed_distance > 0.25 && chamfer.signed_distance < 0.3);
    });
}

#[test]
fn g0_bore_is_supported_and_axis_closures_are_not_false_boundary_features() {
    let bore = AxisymmetricChart::try_new(vec![
        line(point(1.0, -1.0), point(2.0, -1.0)),
        line(point(2.0, -1.0), point(2.0, 1.0)),
        line(point(2.0, 1.0), point(1.0, 1.0)),
        line(point(1.0, 1.0), point(1.0, -1.0)),
    ])
    .expect("annular meridian is a valid bore");
    with_cx(|cx| {
        assert_eq!(
            bore.eval(Point3::new(0.0, 0.0, 0.0), cx).signed_distance,
            1.0
        );
        assert_eq!(
            bore.eval(Point3::new(1.5, 0.0, 0.0), cx).signed_distance,
            -0.5
        );
    });
    assert_eq!(bore.topology_hint(), BettiBounds::unknown());
}

#[test]
fn g0_hostile_profiles_refuse_before_a_chart_is_exposed() {
    let open = AxisymmetricChart::try_new(vec![
        line(point(0.0, 0.0), point(1.0, 0.0)),
        line(point(1.0, 0.0), point(1.0, 1.0)),
        line(point(1.0, 1.0), point(0.0, 1.1)),
    ]);
    assert!(matches!(open, Err(AxisymmetricError::OpenLoop { .. })));

    let crossing = AxisymmetricChart::try_new(vec![
        line(point(0.0, 0.0), point(2.0, 2.0)),
        line(point(2.0, 2.0), point(0.0, 2.0)),
        line(point(0.0, 2.0), point(2.0, 0.0)),
        line(point(2.0, 0.0), point(0.0, 0.0)),
    ]);
    assert!(matches!(
        crossing,
        Err(AxisymmetricError::SelfIntersection { .. })
            | Err(AxisymmetricError::NonPositiveOrientation)
    ));

    let reversed_cylinder = AxisymmetricChart::try_new(vec![
        line(point(0.0, -1.0), point(0.0, 1.0)),
        line(point(0.0, 1.0), point(1.0, 1.0)),
        line(point(1.0, 1.0), point(1.0, -1.0)),
        line(point(1.0, -1.0), point(0.0, -1.0)),
    ]);
    assert!(matches!(
        reversed_cylinder,
        Err(AxisymmetricError::NonPositiveOrientation)
    ));

    let too_many = AxisymmetricChart::try_new(vec![
        line(point(0.0, 0.0), point(1.0, 0.0));
        MAX_AXISYMMETRIC_SEGMENTS + 1
    ]);
    assert!(matches!(
        too_many,
        Err(AxisymmetricError::SegmentCount { .. })
    ));
}

#[test]
fn g0_interior_arc_tangent_to_axis_refuses_singular_pinch_topology() {
    let singular_pinch = AxisymmetricChart::try_new(vec![
        MeridianSegment::Arc {
            start: point(1.0, 1.0),
            end: point(1.0, -1.0),
            center: point(1.0, 0.0),
            clockwise: false,
        },
        line(point(1.0, -1.0), point(3.0, -1.0)),
        line(point(3.0, -1.0), point(3.0, 1.0)),
        line(point(3.0, 1.0), point(1.0, 1.0)),
    ]);
    assert!(matches!(
        singular_pinch,
        Err(AxisymmetricError::AxisTangentArc { index: 0 })
    ));
}

#[test]
fn g0_filleted_join_heights_keep_the_meridian_inside_rule_half_open() {
    let chart = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
    )
    .expect("filleted disc");
    with_cx(|cx| {
        for axial in [-0.5, 0.5] {
            let sample = chart.eval(Point3::new(1.75, 0.0, axial), cx);
            assert!(
                sample.signed_distance < 0.0,
                "meridian point at fillet/vertical join z={axial} must stay inside"
            );
        }
    });
}

#[test]
fn g3_axial_translation_unit_rescaling_and_z_rotation_preserve_distance_laws() {
    let base = cylinder(2.0, 1.0, 0.0);
    let translated = cylinder(2.0, 1.0, 7.0);
    let scaled = cylinder(20.0, 10.0, 0.0);
    with_cx(|cx| {
        let p = Point3::new(2.4, -0.6, 0.25);
        let base_value = base.eval(p, cx).signed_distance;
        let angle: f64 = 0.713;
        let rotated = Point3::new(
            p.x * angle.cos() - p.y * angle.sin(),
            p.x * angle.sin() + p.y * angle.cos(),
            p.z,
        );
        assert_eq!(
            base_value.to_bits(),
            base.eval(rotated, cx).signed_distance.to_bits()
        );
        let shifted = Point3::new(p.x, p.y, p.z + 7.0);
        assert_eq!(
            base_value.to_bits(),
            translated.eval(shifted, cx).signed_distance.to_bits()
        );
        let ten_x = Point3::new(p.x * 10.0, p.y * 10.0, p.z * 10.0);
        assert!((scaled.eval(ten_x, cx).signed_distance - 10.0 * base_value).abs() < 2e-12);
    });
}
