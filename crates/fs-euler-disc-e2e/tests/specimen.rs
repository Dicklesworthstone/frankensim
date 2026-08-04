//! G0/G3 checks for one-source-of-truth Euler-disc specimen resolution.

use fs_euler_disc_e2e::specimen::{DiscProfileError, DiscProfileSpec};
use fs_exec::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_frep::SquatDiscEdgeTreatment;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_5052,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 3.0e-10 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.17e}"
    );
}

#[test]
fn g0_solid_cylinder_resolution_matches_closed_form_mass_and_inertia() {
    let radius = 0.038;
    let thickness = 0.006;
    let density = 7_800.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: radius,
            thickness_m: thickness,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(density, cx)
        .expect("admitted sharp cylinder")
    });
    let volume = core::f64::consts::PI * radius.powi(2) * thickness;
    let mass = density * volume;
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(resolved.mass_properties.center_of_mass.z, 0.0);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.5 * mass * radius.powi(2),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (3.0 * radius.powi(2) + thickness.powi(2)) / 12.0,
    );
}

#[test]
fn g0_annular_cylinder_resolution_uses_bore_for_mass_not_only_inertia() {
    let outer = 0.038;
    let inner = 0.021;
    let thickness = 0.006;
    let density = 2_680.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::AnnularCylinder {
            outer_radius_m: outer,
            inner_radius_m: inner,
            thickness_m: thickness,
        }
        .resolve(density, cx)
        .expect("admitted annular cylinder")
    });
    let volume = core::f64::consts::PI * (outer.powi(2) - inner.powi(2)) * thickness;
    let mass = density * volume;
    assert!(
        !resolved.chart.construction_certificate().touches_axis,
        "an annular profile must retain its bore rather than gain an axis closure"
    );
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.5 * mass * (outer.powi(2) + inner.powi(2)),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (3.0 * (outer.powi(2) + inner.powi(2)) + thickness.powi(2)) / 12.0,
    );
}

#[test]
fn g0_outer_filleted_annulus_uses_one_chart_for_equal_mass_geometry_and_inertia() {
    let outer = 0.038;
    let inner = 0.021;
    let thickness = 0.006;
    let fillet = 0.001;
    let solid_density = 7_800.0;
    let solid = DiscProfileSpec::SolidCylinder {
        outer_radius_m: outer,
        thickness_m: thickness,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
    };
    let annulus = DiscProfileSpec::OuterFilletedAnnularCylinder {
        outer_radius_m: outer,
        inner_radius_m: inner,
        thickness_m: thickness,
        outer_fillet_radius_m: fillet,
    };
    let (solid, unit_density_annulus) = with_cx(|cx| {
        (
            solid
                .resolve(solid_density, cx)
                .expect("filleted solid baseline"),
            annulus
                .resolve(1.0, cx)
                .expect("unit-density filleted annulus"),
        )
    });
    let equal_mass_density =
        solid.mass_properties.mass / unit_density_annulus.mass_properties.volume;
    let ring = with_cx(|cx| {
        annulus
            .resolve(equal_mass_density, cx)
            .expect("equal-mass filleted annulus")
    });

    assert!(
        !ring.chart.construction_certificate().touches_axis,
        "the equal-mass control must retain a real through bore"
    );
    assert!(
        ring.chart
            .segments()
            .iter()
            .filter(|segment| matches!(segment, fs_rep_frep::MeridianSegment::Arc { .. }))
            .count()
            == 2,
        "the ring's outer contact geometry must be circular rather than sharp"
    );
    assert!(
        equal_mass_density > solid_density,
        "equal mass is achieved by the declared density control, not an inertia multiplier"
    );
    assert_close(ring.mass_properties.mass, solid.mass_properties.mass);
    assert_close(ring.mass_properties.center_of_mass.z, 0.0);
    assert!(
        ring.mass_properties.principal_inertia.axial
            > solid.mass_properties.principal_inertia.axial
    );
    assert!(
        ring.mass_properties
            .principal_inertia
            .transverse
            .is_finite()
    );
}

#[test]
fn g0_true_symmetric_bicone_matches_closed_forms_and_is_centered() {
    let radius = 0.038;
    let thickness = 0.006;
    let density = 2_680.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::SymmetricTapered {
            outer_radius_m: radius,
            face_radius_m: 0.0,
            thickness_m: thickness,
        }
        .resolve(density, cx)
        .expect("admitted symmetric bicone")
    });
    // The profile is two equal right cones, each height thickness / 2.
    let volume = core::f64::consts::PI * radius.powi(2) * thickness / 3.0;
    let mass = density * volume;
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(resolved.mass_properties.center_of_mass.z, 0.0);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.3 * mass * radius.powi(2),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (6.0 * radius.powi(2) + thickness.powi(2)) / 40.0,
    );
    assert!(
        resolved.chart.construction_certificate().touches_axis,
        "the zero-face-radius bicone must retain its two axis tips"
    );
}

#[test]
fn g0_full_height_chamfers_omit_the_collapsed_cylindrical_band() {
    let resolved = with_cx(|cx| {
        DiscProfileSpec::ChamferedCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            chamfer_radial_m: 0.001,
            chamfer_axial_m: 0.003,
        }
        .resolve(2_680.0, cx)
        .expect("full-height chamfers form a valid double-frustum rather than a zero line")
    });
    let certificate = resolved.chart.construction_certificate();
    assert_eq!(certificate.input_feature_count, 5);
    assert_eq!(certificate.surfaced_feature_count, 4);
    assert!(resolved.mass_properties.mass > 0.0);
}

#[test]
fn g0_zero_fillet_is_the_same_resolved_sharp_profile() {
    let sharp = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(7_800.0, cx)
        .expect("sharp")
    });
    let zero_fillet = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.0 },
        }
        .resolve(7_800.0, cx)
        .expect("zero fillet")
    });
    assert_eq!(sharp.identity, zero_fillet.identity);
    assert_eq!(sharp.chart.segments(), zero_fillet.chart.segments());
    assert_eq!(
        sharp.mass_properties.mass.to_bits(),
        zero_fillet.mass_properties.mass.to_bits()
    );
    assert_eq!(
        sharp.mass_properties.principal_inertia,
        zero_fillet.mass_properties.principal_inertia
    );
}

#[test]
fn g3_uniform_scaling_preserves_exact_dimensional_mass_laws_for_chamfers() {
    let base = DiscProfileSpec::ChamferedCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        chamfer_radial_m: 0.001,
        chamfer_axial_m: 0.0015,
    };
    let scale = 2.5;
    let scaled = DiscProfileSpec::ChamferedCylinder {
        outer_radius_m: 0.038 * scale,
        thickness_m: 0.006 * scale,
        chamfer_radial_m: 0.001 * scale,
        chamfer_axial_m: 0.0015 * scale,
    };
    let (base, scaled) = with_cx(|cx| {
        (
            base.resolve(2_680.0, cx).expect("base chamfer"),
            scaled.resolve(2_680.0, cx).expect("scaled chamfer"),
        )
    });
    assert_close(
        scaled.mass_properties.volume,
        scale.powi(3) * base.mass_properties.volume,
    );
    assert_close(
        scaled.mass_properties.mass,
        scale.powi(3) * base.mass_properties.mass,
    );
    assert_close(
        scaled.mass_properties.principal_inertia.axial,
        scale.powi(5) * base.mass_properties.principal_inertia.axial,
    );
    assert_close(
        scaled.mass_properties.principal_inertia.transverse,
        scale.powi(5) * base.mass_properties.principal_inertia.transverse,
    );
}

#[test]
fn g0_invalid_profile_relationships_refuse_before_mass_properties_publish() {
    with_cx(|cx| {
        assert!(matches!(
            DiscProfileSpec::AnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.038,
                thickness_m: 0.006,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::SymmetricTapered {
                outer_radius_m: 0.038,
                face_radius_m: 0.038,
                thickness_m: 0.006,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::ChamferedCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                chamfer_radial_m: 0.001,
                chamfer_axial_m: 0.0031,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::OuterFilletedAnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.021,
                thickness_m: 0.006,
                outer_fillet_radius_m: 0.0,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidParameter {
                field: "outer_fillet_radius_m",
                ..
            })
        ));
        assert!(matches!(
            DiscProfileSpec::OuterFilletedAnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.021,
                thickness_m: 0.006,
                outer_fillet_radius_m: 0.0031,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::SolidCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                edge_treatment: SquatDiscEdgeTreatment::Sharp,
            }
            .resolve(0.0, cx),
            Err(DiscProfileError::InvalidParameter {
                field: "density_kg_per_m3",
                ..
            })
        ));
    });
}
