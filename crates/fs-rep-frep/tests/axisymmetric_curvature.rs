//! G0 local differential-geometry checks for selected axisymmetric support.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricCurvatureAuthority, AxisymmetricCurvatureError,
    SquatDiscEdgeTreatment,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4355_5256,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-12,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g0_true_fillet_support_reports_local_line_arc_curvatures_as_estimates() {
    let chart = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
    )
    .expect("filleted squat disc");
    let estimate = with_cx(|cx| {
        chart
            .principal_curvatures_at_support_direction(Vec3::new(1.0, 0.0, -1.0), cx)
            .expect("interior arc curvature")
    });
    assert_eq!(estimate.source_feature, 3);
    assert_eq!(estimate.authority, AxisymmetricCurvatureAuthority::Estimate);
    assert_close(estimate.meridional_m_inverse, 2.0);
    assert_close(
        estimate.azimuthal_m_inverse,
        0.5_f64.sqrt() / (1.5 + 0.5_f64.sqrt() * 0.5),
    );
    assert!(estimate.uncertainty_m_inverse.is_finite());
    assert!(estimate.uncertainty_m_inverse > 0.0);
}

#[test]
fn g0_axis_feature_boundaries_and_wrong_feature_points_refuse() {
    let chart = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
    )
    .expect("filleted squat disc");
    with_cx(|cx| {
        assert!(matches!(
            chart.principal_curvatures_at_feature_point(3, Point3::new(2.0, 0.0, 0.5), cx),
            Err(AxisymmetricCurvatureError::NonsmoothFeatureBoundary { .. })
        ));
        assert!(matches!(
            chart.principal_curvatures_at_feature_point(0, Point3::new(0.0, 0.0, -1.0), cx),
            Err(AxisymmetricCurvatureError::AxisPoint { .. })
        ));
        assert!(matches!(
            chart.principal_curvatures_at_feature_point(3, Point3::new(1.0, 0.0, 0.0), cx),
            Err(AxisymmetricCurvatureError::PointNotOnSelectedFeature { .. })
        ));
    });
}

#[test]
fn g3_directional_entrypoint_cannot_promote_a_forged_on_feature_point_to_support() {
    let chart = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
    )
    .expect("filleted squat disc");
    with_cx(|cx| {
        // The lower fillet interior is a legitimate local point, but for this
        // direction the global support search selects the upper fillet.  The
        // selected-support API takes only direction, so it cannot be fed a
        // forged feature/point record.
        let local_only = chart
            .principal_curvatures_at_feature_point(
                1,
                Point3::new(1.5 + 0.5_f64.sqrt() * 0.5, 0.0, -0.5 - 0.5_f64.sqrt() * 0.5),
                cx,
            )
            .expect("lower fillet is locally smooth");
        assert_eq!(local_only.source_feature, 1);
        let selected = chart
            .principal_curvatures_at_support_direction(Vec3::new(1.0, 0.0, -1.0), cx)
            .expect("recomputed selected support");
        assert_eq!(selected.source_feature, 3);
    });
}
