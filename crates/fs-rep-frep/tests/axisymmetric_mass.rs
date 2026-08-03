//! G0/G3 mass-property checks for exact axisymmetric line/arc solids.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricMassError, MeridianPoint, MeridianSegment,
    SquatDiscEdgeTreatment,
};

fn point(radius: f64, axial: f64) -> MeridianPoint {
    MeridianPoint::new(radius, axial)
}

fn line(start: MeridianPoint, end: MeridianPoint) -> MeridianSegment {
    MeridianSegment::Line { start, end }
}

fn cylinder(radius: f64, thickness: f64, axial_offset: f64) -> AxisymmetricChart {
    let half = 0.5 * thickness;
    AxisymmetricChart::try_new(vec![
        line(
            point(0.0, axial_offset - half),
            point(radius, axial_offset - half),
        ),
        line(
            point(radius, axial_offset - half),
            point(radius, axial_offset + half),
        ),
        line(
            point(radius, axial_offset + half),
            point(0.0, axial_offset + half),
        ),
        line(
            point(0.0, axial_offset + half),
            point(0.0, axial_offset - half),
        ),
    ])
    .expect("closed CCW cylinder meridian")
}

fn refined_cylinder(radius: f64, thickness: f64) -> AxisymmetricChart {
    let half = 0.5 * thickness;
    AxisymmetricChart::try_new(vec![
        line(point(0.0, -half), point(0.5 * radius, -half)),
        line(point(0.5 * radius, -half), point(radius, -half)),
        line(point(radius, -half), point(radius, 0.0)),
        line(point(radius, 0.0), point(radius, half)),
        line(point(radius, half), point(0.0, half)),
        line(point(0.0, half), point(0.0, -half)),
    ])
    .expect("feature-refined cylinder is the same closed CCW profile")
}

fn sphere(radius: f64, axial_offset: f64) -> AxisymmetricChart {
    AxisymmetricChart::try_new(vec![
        MeridianSegment::Arc {
            start: point(0.0, axial_offset - radius),
            end: point(radius, axial_offset),
            center: point(0.0, axial_offset),
            clockwise: false,
        },
        MeridianSegment::Arc {
            start: point(radius, axial_offset),
            end: point(0.0, axial_offset + radius),
            center: point(0.0, axial_offset),
            clockwise: false,
        },
        line(
            point(0.0, axial_offset + radius),
            point(0.0, axial_offset - radius),
        ),
    ])
    .expect("two exact semicircle arcs and an axis closure bound a sphere")
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
                seed: 0x4D415353,
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
    let tolerance = 2e-10 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.17e}"
    );
}

#[test]
fn g0_cylinder_mass_center_and_principal_inertia_match_closed_forms() {
    let radius = 2.0;
    let thickness = 4.0;
    let density = 3.0;
    let properties = with_cx(false, |cx| {
        cylinder(radius, thickness, 0.0)
            .mass_properties(density, cx)
            .expect("valid cylinder mass")
    });
    let volume = core::f64::consts::PI * radius * radius * thickness;
    let mass = density * volume;

    assert_close(properties.volume, volume);
    assert_close(properties.mass, mass);
    assert_eq!(properties.center_of_mass, Point3::new(0.0, 0.0, 0.0));
    assert_close(
        properties.principal_inertia.axial,
        0.5 * mass * radius * radius,
    );
    assert_close(
        properties.principal_inertia.transverse,
        mass * (3.0 * radius * radius + thickness * thickness) / 12.0,
    );
    assert_eq!(properties.principal_inertia, properties.origin_inertia);
    assert!(properties.roundoff_diagnostics.volume_term_scale >= properties.volume.abs());
    assert!(
        properties
            .roundoff_diagnostics
            .centroidal_transverse_term_scale
            >= properties.principal_inertia.transverse.abs()
    );
    for diagnostic in [
        properties.roundoff_diagnostics.volume_term_scale,
        properties
            .roundoff_diagnostics
            .axial_first_moment_term_scale,
        properties
            .roundoff_diagnostics
            .centroidal_transverse_term_scale,
        properties.roundoff_diagnostics.axial_inertia_term_scale,
    ] {
        assert!(
            diagnostic.is_finite(),
            "published diagnostic must be finite"
        );
    }
}

#[test]
fn g0_exact_arc_sphere_matches_volume_and_isotropic_inertia() {
    let radius = 1.75;
    let density = 2.5;
    let properties = with_cx(false, |cx| {
        sphere(radius, 0.0)
            .mass_properties(density, cx)
            .expect("exact circular meridian sphere mass")
    });
    let volume = 4.0 * core::f64::consts::PI * radius.powi(3) / 3.0;
    let mass = density * volume;
    let inertia = 0.4 * mass * radius * radius;

    assert_close(properties.volume, volume);
    assert_close(properties.mass, mass);
    assert_close(properties.center_of_mass.z, 0.0);
    assert_close(properties.principal_inertia.transverse, inertia);
    assert_close(properties.principal_inertia.axial, inertia);
}

#[test]
fn g0_nominal_millimetre_filleted_disc_constructs_and_feeds_support_and_mass() {
    let chart = AxisymmetricChart::squat_disc(
        0.038,
        0.006,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    )
    .expect("nominal 38 mm by 6 mm disc with a 1 mm fillet");
    with_cx(false, |cx| {
        let support = chart
            .minimum_support_point(Vec3::new(1.0, 0.0, 1.0), cx)
            .expect("tilted nominal disc support is unique");
        assert!(support.point.x.is_finite());
        assert!(support.point.z.is_finite());

        let properties = chart
            .mass_properties(7_800.0, cx)
            .expect("nominal disc mass properties");
        assert!(properties.volume.is_finite() && properties.volume > 0.0);
        assert!(properties.mass.is_finite() && properties.mass > 0.0);
        assert!(properties.principal_inertia.transverse > 0.0);
        assert!(properties.principal_inertia.axial > 0.0);
    });
}

#[test]
fn g3_fillet_removes_material_monotonically_without_polygonizing_the_rim() {
    let sharp =
        AxisymmetricChart::squat_disc(2.0, 2.0, SquatDiscEdgeTreatment::Sharp).expect("sharp disc");
    let small = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.25 },
    )
    .expect("true small fillet");
    let large = AxisymmetricChart::squat_disc(
        2.0,
        2.0,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.75 },
    )
    .expect("true large fillet");
    let (sharp_mass, small_mass, large_mass) = with_cx(false, |cx| {
        (
            sharp.mass_properties(1.0, cx).expect("sharp mass"),
            small.mass_properties(1.0, cx).expect("small fillet mass"),
            large.mass_properties(1.0, cx).expect("large fillet mass"),
        )
    });

    assert!(sharp_mass.volume > small_mass.volume);
    assert!(small_mass.volume > large_mass.volume);
    assert_eq!(small_mass.center_of_mass, Point3::new(0.0, 0.0, 0.0));
    assert_eq!(large_mass.center_of_mass, Point3::new(0.0, 0.0, 0.0));
}

#[test]
fn g3_scale_and_axial_translation_obey_dimensional_and_parallel_axis_laws() {
    let density = 4.0;
    let base = with_cx(false, |cx| {
        cylinder(1.25, 0.8, 0.0)
            .mass_properties(density, cx)
            .expect("base mass")
    });
    let scale = 3.0;
    let scaled = with_cx(false, |cx| {
        cylinder(scale * 1.25, scale * 0.8, 0.0)
            .mass_properties(density, cx)
            .expect("scaled mass")
    });
    assert_close(scaled.volume, scale.powi(3) * base.volume);
    assert_close(scaled.mass, scale.powi(3) * base.mass);
    assert_close(
        scaled.principal_inertia.transverse,
        scale.powi(5) * base.principal_inertia.transverse,
    );
    assert_close(
        scaled.principal_inertia.axial,
        scale.powi(5) * base.principal_inertia.axial,
    );

    let axial_shift = 2.75;
    let translated = with_cx(false, |cx| {
        cylinder(1.25, 0.8, axial_shift)
            .mass_properties(density, cx)
            .expect("translated mass")
    });
    assert_close(translated.center_of_mass.z, axial_shift);
    assert_close(
        translated.principal_inertia.transverse,
        base.principal_inertia.transverse,
    );
    assert_close(
        translated.principal_inertia.axial,
        base.principal_inertia.axial,
    );
    assert_close(
        translated.origin_inertia.transverse,
        translated.principal_inertia.transverse + translated.mass * axial_shift * axial_shift,
    );
    assert_close(
        translated.origin_inertia.axial,
        translated.principal_inertia.axial,
    );
}

#[test]
fn g3_huge_axial_translation_keeps_centroidal_inertia_out_of_parallel_axis_cancellation() {
    let density = 4.0;
    let base = with_cx(false, |cx| {
        cylinder(1.25, 4.0, 0.0)
            .mass_properties(density, cx)
            .expect("base mass")
    });
    let axial_shift = 1.0e12;
    let translated = with_cx(false, |cx| {
        cylinder(1.25, 4.0, axial_shift)
            .mass_properties(density, cx)
            .expect("huge translated mass")
    });

    assert_close(translated.center_of_mass.z, axial_shift);
    assert_close(
        translated.principal_inertia.transverse,
        base.principal_inertia.transverse,
    );
    assert_close(
        translated.principal_inertia.axial,
        base.principal_inertia.axial,
    );
    assert_close(
        translated.origin_inertia.transverse,
        translated.principal_inertia.transverse + translated.mass * axial_shift * axial_shift,
    );
}

#[test]
fn g3_feature_refinement_and_repeated_evaluation_are_deterministic() {
    let unrefined = cylinder(1.8, 1.2, 0.0);
    let refined = refined_cylinder(1.8, 1.2);
    let (first, second, refinement) = with_cx(false, |cx| {
        let first = unrefined
            .mass_properties(1.5, cx)
            .expect("first exact evaluation");
        let second = unrefined
            .mass_properties(1.5, cx)
            .expect("second exact evaluation");
        let refinement = refined
            .mass_properties(1.5, cx)
            .expect("refined profile mass");
        (first, second, refinement)
    });
    assert_eq!(first.volume.to_bits(), second.volume.to_bits());
    assert_eq!(
        first.principal_inertia.transverse.to_bits(),
        second.principal_inertia.transverse.to_bits()
    );
    assert_eq!(
        first.principal_inertia.axial.to_bits(),
        second.principal_inertia.axial.to_bits()
    );
    assert_close(refinement.volume, first.volume);
    assert_close(
        refinement.principal_inertia.transverse,
        first.principal_inertia.transverse,
    );
    assert_close(
        refinement.principal_inertia.axial,
        first.principal_inertia.axial,
    );
}

#[test]
fn g0_invalid_density_and_cancelled_work_refuse_without_partial_properties() {
    let chart = cylinder(1.0, 1.0, 0.0);
    with_cx(false, |cx| {
        assert!(matches!(
            chart.mass_properties(0.0, cx),
            Err(AxisymmetricMassError::InvalidDensity { .. })
        ));
        assert!(matches!(
            chart.mass_properties(f64::NAN, cx),
            Err(AxisymmetricMassError::InvalidDensity { .. })
        ));
    });
    with_cx(true, |cx| {
        assert!(matches!(
            chart.mass_properties(1.0, cx),
            Err(AxisymmetricMassError::Cancelled)
        ));
    });
}

#[test]
fn g0_positive_density_with_representational_mass_underflow_refuses() {
    let chart = cylinder(1.0e-8, 1.0e-8, 0.0);
    with_cx(false, |cx| {
        assert!(matches!(
            chart.mass_properties(f64::MIN_POSITIVE, cx),
            Err(AxisymmetricMassError::NonPositiveMass { mass: 0.0 })
        ));
    });
}

#[test]
fn g0_positive_mass_with_underflowed_principal_inertia_refuses() {
    // This ordinary cylinder has volume above one, so multiplying by the least
    // positive subnormal leaves a positive representable mass. Its axial
    // inertia coefficient is below one half, which underflows the principal
    // moment without relying on an extreme meridian that construction should
    // reject before mechanics begins.
    let chart = cylinder(0.5, 2.0, 0.0);
    let density = f64::from_bits(1);
    with_cx(false, |cx| {
        assert!(matches!(
            chart.mass_properties(density, cx),
            Err(AxisymmetricMassError::NonPositivePrincipalInertia { .. })
        ));
    });
}
