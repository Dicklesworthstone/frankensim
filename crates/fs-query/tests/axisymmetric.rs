//! Comprehensive test battery for axisymmetric query adapters (bead frankensim-b8bxd.2).
//!
//! Verifies:
//! - [`AxisymmetricSupportMap`]: global analytic maximization over line/arc meridian,
//!   deterministic zero-direction handling, certified support slack, contained ball inradius;
//! - Refusal on non-convex profiles (annular bore / central hole, concave arcs);
//! - Contact inflation: monotone outward expansion of support slack and brackets;
//! - [`axisymmetric_normal`]: unique Smooth/Axis normals and fail-closed boundary admission;
//! - [`axisymmetric_curvature`]: meridional and azimuthal curvatures (Meusnier theorem),
//!   mean/Gaussian magnitude diagnostics, and explicit Estimate authority;
//! - [`axisymmetric_reach`]: explicit Estimate-only reach-related feature scales;
//! - [`AxisymmetricGapOracle`]: pointwise gap, separation upper bounds, overlap witnesses;
//! - Convex separation between [`AxisymmetricSupportMap`] and [`ConvexSphere`] / [`ConvexBox`].

use asupersync::types::Budget;
use fs_evidence::{Certified, Evidence, NumericalKind, ProvenanceHash};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Point3, Vec3};
use fs_query::{
    AxisymmetricGapOracle, AxisymmetricSupportMap, ContactInflation, ConvexSphere,
    ConvexSupportMap, NormalClassification, QueryError, axisymmetric_curvature,
    axisymmetric_normal, axisymmetric_reach, convex_separation,
};
use fs_rep_frep::axisymmetric::{
    AxisymmetricChart, MeridianPoint, MeridianSegment, SquatDiscEdgeTreatment,
};

fn meridian_line(start: (f64, f64), end: (f64, f64)) -> MeridianSegment {
    MeridianSegment::Line {
        start: MeridianPoint::new(start.0, start.1),
        end: MeridianPoint::new(end.0, end.1),
    }
}

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xB8B_2001,
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

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate_cx(&CancelGate::new(), f)
}

fn conversion_receipt(radius: f64) -> Certified<f64> {
    Evidence::enclosed(
        radius,
        0.0,
        radius,
        ProvenanceHash::of_bytes(b"fs-query/axisymmetric/test"),
    )
    .certified()
    .expect("valid certified evidence")
}

#[test]
fn ax_001_convex_support_sharp_disc_exact_extrema() {
    let outer_radius = 1.0;
    let thickness = 0.4;
    let disc =
        AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
            .expect("sharp disc");

    let map = AxisymmetricSupportMap::try_new(disc).expect("support map");
    assert_eq!(map.max_radius(), 1.0);
    assert_eq!(map.axial_bounds(), (-0.2, 0.2));
    assert!(map.support_slack() > 0.0);
    assert_eq!(map.interior_point(), Point3::new(0.0, 0.0, 0.0));

    // Support in +z direction -> should hit upper flat face z = 0.2
    let s_top = map.support_point(Vec3::new(0.0, 0.0, 1.0));
    assert!((s_top.z - 0.2).abs() < 1e-12);

    // Support in -z direction -> should hit lower flat face z = -0.2
    let s_bot = map.support_point(Vec3::new(0.0, 0.0, -1.0));
    assert!((s_bot.z - (-0.2)).abs() < 1e-12);

    // Support in +x direction -> should hit rim x = 1.0, y = 0.0
    let s_x = map.support_point(Vec3::new(1.0, 0.0, 0.0));
    assert!((s_x.x - 1.0).abs() < 1e-12);
    assert!(s_x.y.abs() < 1e-12);

    // G3: support selection is homogeneous even when a naive squared norm
    // would overflow or underflow.
    assert_eq!(map.support_point(Vec3::new(f64::MAX, 0.0, 0.0)), s_x);
    assert_eq!(
        map.support_point(Vec3::new(f64::from_bits(1), 0.0, 0.0)),
        s_x
    );

    // Any nonzero radial component selects the rim of the top face. Treating
    // a tiny component as axial would return a point one metre from the true
    // supporting set while declaring only rounding-scale support slack.
    let almost_axial = map.support_point(Vec3::new(f64::from_bits(1), 0.0, 1.0));
    assert_eq!(almost_axial.x, 1.0);
    assert_eq!(almost_axial.z, 0.2);

    // Support in diagonal (+x, +z)
    let s_diag = map.support_point(Vec3::new(1.0, 0.0, 1.0));
    assert!((s_diag.x - 1.0).abs() < 1e-12);
    assert!((s_diag.z - 0.2).abs() < 1e-12);

    // Deterministic zero-direction behavior
    let s_zero = map.support_point(Vec3::new(0.0, 0.0, 0.0));
    assert!(s_zero.x.is_finite() && s_zero.y.is_finite() && s_zero.z.is_finite());
}

#[test]
fn ax_002_convex_support_filleted_disc_smooth_transitions() {
    let outer_radius = 1.0;
    let thickness = 0.4;
    let fillet = 0.1;
    let disc = AxisymmetricChart::squat_disc(
        outer_radius,
        thickness,
        SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
    )
    .expect("filleted disc");

    let map = AxisymmetricSupportMap::try_new(disc).expect("support map");

    // Diagonal direction (1.0, 0.0, 1.0) normalized is (1/sqrt(2), 0, 1/sqrt(2))
    // Center of fillet is (0.9, 0.1)
    // Support on circular fillet should be (0.9 + 0.1/sqrt(2), 0, 0.1 + 0.1/sqrt(2))
    let s_diag = map.support_point(Vec3::new(1.0, 0.0, 1.0));
    let expected_x = 0.9 + 0.1 / core::f64::consts::SQRT_2;
    let expected_z = 0.1 + 0.1 / core::f64::consts::SQRT_2;
    assert!((s_diag.x - expected_x).abs() < 1e-10);
    assert!((s_diag.z - expected_z).abs() < 1e-10);
}

#[test]
fn ax_003_refuses_nonconvex_annulus_and_concave_profiles() {
    let outer_radius = 1.0;
    let inner_radius = 0.3;
    let thickness = 0.4;
    let fillet = 0.05;

    // Annular disc has a central bore -> non-convex
    let annulus = AxisymmetricChart::annular_disc_outer_fillets(
        outer_radius,
        inner_radius,
        thickness,
        fillet,
    )
    .expect("annular disc");

    let result = AxisymmetricSupportMap::try_new(annulus);
    assert!(matches!(result, Err(QueryError::ConvexInvalidShape { .. })));

    // This is a valid simple CCW chart, but its inward notch is a reflex
    // line-line join. A support map for the chart itself must refuse instead
    // of silently representing its convex hull.
    let concave = AxisymmetricChart::try_new(vec![
        meridian_line((0.0, -1.0), (1.0, -1.0)),
        meridian_line((1.0, -1.0), (0.5, 0.0)),
        meridian_line((0.5, 0.0), (1.0, 1.0)),
        meridian_line((1.0, 1.0), (0.0, 1.0)),
        meridian_line((0.0, 1.0), (0.0, -1.0)),
    ])
    .expect("valid concave chart");
    assert!(matches!(
        AxisymmetricSupportMap::try_new(concave),
        Err(QueryError::ConvexInvalidShape { .. })
    ));
}

#[test]
fn ax_004_contact_inflation_monotonically_widens_support() {
    let outer_radius = 1.0;
    let thickness = 0.4;
    let disc =
        AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
            .expect("sharp disc");

    let base_map = AxisymmetricSupportMap::try_new(disc.clone()).expect("base support map");
    let cert_1 = conversion_receipt(0.05);
    let cert_2 = conversion_receipt(0.10);
    let inflation_1 = ContactInflation::from_conversion(&cert_1).expect("valid inflation 1");
    let inflation_2 = ContactInflation::from_conversion(&cert_2).expect("valid inflation 2");

    let map_1 = AxisymmetricSupportMap::try_new_with_inflation(disc.clone(), inflation_1)
        .expect("inflated map 1");
    let map_2 =
        AxisymmetricSupportMap::try_new_with_inflation(disc, inflation_2).expect("inflated map 2");

    assert!(map_1.support_slack() > base_map.support_slack());
    assert!(map_2.support_slack() > map_1.support_slack());

    let dir = Vec3::new(1.0, 0.0, 0.0);
    let p_base = base_map.support_point(dir);
    let p_1 = map_1.support_point(dir);
    let p_2 = map_2.support_point(dir);

    assert!((p_1.x - (p_base.x + 0.05)).abs() < 1e-12);
    assert!((p_2.x - (p_base.x + 0.10)).abs() < 1e-12);
}

#[test]
fn ax_005_contained_ball_inradius_evaluation() {
    let outer_radius = 1.0;
    let thickness = 0.4;
    let disc =
        AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
            .expect("sharp disc");

    let map = AxisymmetricSupportMap::try_new(disc).expect("support map");

    // Center (0, 0, 0): distance to top/bottom is 0.2, distance to outer wall is 1.0
    // Inradius should be 0.2
    let inrad = map.contained_ball_radius(Point3::new(0.0, 0.0, 0.0));
    assert!(inrad.is_some());
    assert!((inrad.unwrap() - 0.2).abs() < 1e-10);

    // Outside point (0, 0, 0.5) -> None
    assert!(
        map.contained_ball_radius(Point3::new(0.0, 0.0, 0.5))
            .is_none()
    );

    let rounded_disc = AxisymmetricChart::squat_disc(
        outer_radius,
        thickness,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.1 },
    )
    .expect("filleted disc");
    let rounded_map = AxisymmetricSupportMap::try_new(rounded_disc).expect("filleted support map");
    assert!(
        rounded_map
            .contained_ball_radius(Point3::new(0.0, 0.0, 0.0))
            .is_some()
    );
    assert!(
        rounded_map
            .contained_ball_radius(Point3::new(0.99, 0.0, 0.19))
            .is_none()
    );

    // A cone exposes why bounding-cylinder membership is insufficient: this
    // probe is within r/z extents but lies well outside the sloped face.
    let cone = AxisymmetricChart::try_new(vec![
        meridian_line((0.0, -1.0), (1.0, -1.0)),
        meridian_line((1.0, -1.0), (0.0, 1.0)),
        meridian_line((0.0, 1.0), (0.0, -1.0)),
    ])
    .expect("valid convex cone");
    let cone_map = AxisymmetricSupportMap::try_new(cone).expect("cone support map");
    assert!(
        cone_map
            .contained_ball_radius(Point3::new(0.9, 0.0, 0.9))
            .is_none()
    );
    assert!(
        cone_map
            .contained_ball_radius(Point3::new(0.1, 0.0, 0.0))
            .is_some()
    );
}

#[test]
fn ax_006_surface_normal_and_curvature_adapters() {
    with_cx(|cx| {
        let outer_radius = 1.0;
        let thickness = 0.4;
        let fillet = 0.1;
        let disc = AxisymmetricChart::squat_disc(
            outer_radius,
            thickness,
            SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
        )
        .expect("filleted disc");

        // Top flat cap: z = 0.2, r = 0.5
        let n_top = axisymmetric_normal(&disc, Point3::new(0.5, 0.0, 0.2), cx).expect("normal");
        assert_eq!(n_top.classification, NormalClassification::Smooth);
        assert!((n_top.normal.z - 1.0).abs() < 1e-10);

        // Cylindrical outer wall: r = 1.0, z = 0.0
        let n_wall = axisymmetric_normal(&disc, Point3::new(1.0, 0.0, 0.0), cx).expect("normal");
        assert_eq!(n_wall.classification, NormalClassification::Smooth);
        assert!((n_wall.normal.x - 1.0).abs() < 1e-10);

        // Axis point: (0, 0, 0.2)
        let n_axis = axisymmetric_normal(&disc, Point3::new(0.0, 0.0, 0.2), cx).expect("normal");
        assert_eq!(n_axis.classification, NormalClassification::Axis);

        // Curvature on flat top: meridional = 0, azimuthal = 0 (since n_r = 0)
        let curv_top =
            axisymmetric_curvature(&disc, Point3::new(0.5, 0.0, 0.2), cx).expect("curvature");
        assert_eq!(curv_top.meridional_curvature, 0.0);
        assert!(curv_top.azimuthal_curvature.abs() < 1e-12);
        assert!(curv_top.mean_curvature.abs() < 1e-12);
        assert_eq!(curv_top.gaussian_curvature, 0.0);
        assert_eq!(curv_top.authority, NumericalKind::Estimate);
        assert!(curv_top.uncertainty_m_inverse.is_finite());

        // Curvature on cylindrical wall (r = 1.0, z = 0.0):
        // meridional = 0, azimuthal = 1.0 / 1.0 = 1.0
        let curv_wall =
            axisymmetric_curvature(&disc, Point3::new(1.0, 0.0, 0.0), cx).expect("curvature");
        assert_eq!(curv_wall.meridional_curvature, 0.0);
        assert!((curv_wall.azimuthal_curvature - 1.0).abs() < 1e-10);
        assert!((curv_wall.mean_curvature - 0.5).abs() < 1e-10);
        assert_eq!(curv_wall.gaussian_curvature, 0.0);
        assert_eq!(curv_wall.authority, NumericalKind::Estimate);
    });
}

#[test]
fn ax_007_reach_scale_is_explicitly_estimate_only() {
    let outer_radius = 1.0;
    let thickness = 0.4;
    let fillet = 0.1;
    let disc = AxisymmetricChart::squat_disc(
        outer_radius,
        thickness,
        SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
    )
    .expect("filleted disc");

    let reach = axisymmetric_reach(&disc).expect("reach");
    assert!((reach.min_meridional_radius_estimate - 0.1).abs() < 1e-10);
    assert!(reach.global_reach_estimate > 0.0);
    assert_eq!(reach.authority, NumericalKind::Estimate);
}

#[test]
fn ax_007a_differential_queries_refuse_points_without_unique_surface_geometry() {
    with_cx(|cx| {
        let disc = AxisymmetricChart::squat_disc(
            1.0,
            0.4,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.1 },
        )
        .expect("filleted disc");

        // Interior points on or near the revolution axis are not surface
        // points. The old r<1e-12 shortcut fabricated a +/-z normal here.
        assert!(matches!(
            axisymmetric_normal(&disc, Point3::new(0.0, 0.0, 0.0), cx),
            Err(QueryError::NotOnBoundary { .. })
        ));
        assert!(matches!(
            axisymmetric_normal(&disc, Point3::new(f64::from_bits(1), 0.0, 0.0), cx),
            Err(QueryError::NotOnBoundary { .. })
        ));

        // This point lies on the upper fillet's full supporting circle but
        // outside the retained quarter-arc. Bounded-arc projection must not
        // turn that circle point into a surface point.
        assert!(matches!(
            axisymmetric_normal(&disc, Point3::new(0.9, 0.0, 0.0), cx),
            Err(QueryError::NotOnBoundary { .. })
        ));

        // The cap/fillet endpoint has agreeing incident normals, so the normal
        // is unique. Curvature still refuses because the local second
        // derivative changes across the feature boundary.
        let tangent_join = Point3::new(0.9, 0.0, 0.2);
        let tangent_normal =
            axisymmetric_normal(&disc, tangent_join, cx).expect("unique tangent normal");
        assert!((tangent_normal.normal.z - 1.0).abs() < 1e-12);
        assert!(matches!(
            axisymmetric_curvature(&disc, tangent_join, cx),
            Err(QueryError::InvalidPointArithmetic { .. })
        ));

        // A sharp disc rim has distinct incident normals and must refuse.
        let sharp = AxisymmetricChart::squat_disc(1.0, 0.4, SquatDiscEdgeTreatment::Sharp)
            .expect("sharp disc");
        let sharp_join = Point3::new(1.0, 0.0, 0.2);
        assert!(matches!(
            axisymmetric_normal(&sharp, sharp_join, cx),
            Err(QueryError::InvalidPointArithmetic { .. })
        ));
        assert!(matches!(
            axisymmetric_curvature(&sharp, sharp_join, cx),
            Err(QueryError::InvalidPointArithmetic { .. })
        ));
    });
}

#[test]
fn ax_007b_surface_normal_observes_preexisting_cancellation() {
    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        let disc = AxisymmetricChart::squat_disc(1.0, 0.4, SquatDiscEdgeTreatment::Sharp)
            .expect("sharp disc");
        assert_eq!(
            axisymmetric_normal(&disc, Point3::new(0.5, 0.0, 0.2), cx),
            Err(QueryError::Cancelled)
        );
    });
}

#[test]
fn ax_008_pointwise_gap_oracle_with_inflation() {
    with_cx(|cx| {
        let outer_radius = 1.0;
        let thickness = 0.4;
        let disc =
            AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
                .expect("sharp disc");

        let sphere = SphereChart {
            center: Point3::new(3.0, 0.0, 0.0),
            radius: 0.5,
        };
        let oracle = AxisymmetricGapOracle::new(&disc, &sphere).expect("gap oracle");

        // Probe midpoint between shapes: (2.0, 0.0, 0.0)
        let sample = oracle
            .gap_at(Point3::new(2.0, 0.0, 0.0), cx)
            .expect("gap sample");
        assert!(sample.separation_upper.is_some());
        assert!(sample.sum_lo > 0.0);
        assert!(sample.overlap_inradius.is_none());
    });
}

#[test]
fn ax_009_convex_separation_against_sphere() {
    with_cx(|cx| {
        let outer_radius = 1.0;
        let thickness = 0.4;
        let disc =
            AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
                .expect("sharp disc");

        let map = AxisymmetricSupportMap::try_new(disc).expect("support map");
        let sphere = ConvexSphere::new(Point3::new(3.0, 0.0, 0.0), 0.5).expect("sphere");

        // Distance between disc (x in [-1, 1]) and sphere (center x = 3, r = 0.5 -> x in [2.5, 3.5])
        // Exact distance is 2.5 - 1.0 = 1.5
        let sep = convex_separation(&map, &sphere, 256, cx).expect("separation");
        assert!(sep.lo <= 1.5 + 1e-6);
        assert!(sep.hi >= 1.5 - 1e-6);
        assert!((sep.hi - 1.5).abs() < 1e-4);
    });
}
