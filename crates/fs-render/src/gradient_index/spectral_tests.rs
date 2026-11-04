//! Source-band refraction uses the original phase, group and absorption laws.

use super::*;
use crate::dielectric::{
    BeerLambertAbsorption, DielectricEvent, DielectricIorSample, GlassProvenance, TabulatedIor,
};
use crate::gradient_index::bundle::event::{
    BundleEventOptions, propagate_seeded_bundle_to_boundary,
};
use crate::gradient_index::bundle::{
    BundleOptions, IndexHessian, ProjectionOptions, QuadraticIndexField, propagate_bundle,
};
use crate::gradient_index::event::{
    BoundaryOracle, BoundaryQuery, BoundarySide, BoundaryTermination, BoundaryWitness,
    GradientEventOptions, propagate_to_boundary,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame, PolarizedBeam};

fn table(index: impl Fn(usize) -> f64, absorption: BeerLambertAbsorption) -> DielectricGlass {
    let samples = core::array::from_fn(|i| {
        DielectricIorSample::try_new_in_band(800.0 + 100.0 * i as f64, index(i)).unwrap()
    });
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(samples, 100.0).unwrap(),
        absorption,
        GlassProvenance::Custom,
    )
}

fn infrared() -> DielectricGlass {
    // Exact linear phase interpolation: n = 2 - lambda_nm/3200, n_g = 2.
    // The absorption source covers only part of the phase-index source.
    table(
        |i| 1.75 - i as f64 / 32.0,
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
            (900.0, 0.25),
            (1500.0, 0.75),
        ])
        .unwrap(),
    )
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn axis() -> Vec3 {
    Vec3::new(1.0, 0.0, 0.0)
}
fn bounds() -> Aabb {
    Aabb::new(Point3::new(-0.3, -0.3, -0.1), Point3::new(0.3, 0.3, 2.0))
}
fn tight() -> GradientIndexOptions {
    GradientIndexOptions {
        position_tolerance_m: 1e-12,
        direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-12,
        optical_depth_tolerance: 1e-12,
        relative_tolerance: 1e-12,
        max_step_m: 0.04,
        ..GradientIndexOptions::default()
    }
}
fn bundle_options() -> BundleOptions {
    BundleOptions {
        integration: tight(),
        position_derivative_tolerance_m: 1e-12,
        direction_derivative_tolerance: 1e-12,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_infrared_native_path_retains_phase_group_and_absorption() {
    let material = infrared();
    let field = HomogeneousIndexField::new(&material);
    for wavelength in [900.0, 1050.0, 1200.0, 1500.0] {
        let path = propagate(&field, ray(), axis(), wavelength, 0.75, tight()).unwrap();
        let n = 2.0 - wavelength / 3200.0;
        let alpha = 0.25 + (wavelength - 900.0) / 1200.0;
        assert_eq!(
            path.ray,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.75),
                ..ray()
            }
        );
        near(path.optical_path_length_m, 0.75 * n, 2e-14);
        near(path.group_optical_path_length_m.unwrap(), 1.5, 2e-14);
        near(path.power_transmittance(), (-0.75 * alpha).exp(), 2e-14);
        assert_eq!(path.field_evaluations, 1);
    }
}

#[test]
fn g1_infrared_affine_refraction_matches_closed_form() {
    let material = infrared();
    let gradient = 0.3;
    let field = AffineIndexField::try_new(
        &material,
        ray().origin,
        Vec3::new(gradient, 0.0, 0.0),
        bounds(),
    )
    .unwrap();
    let distance = 0.8;
    for wavelength in [950.0, 1450.0] {
        let path = propagate(&field, ray(), axis(), wavelength, distance, tight()).unwrap();
        let n0 = 2.0 - wavelength / 3200.0;
        let n_end = (n0 * n0 + (gradient * distance).powi(2)).sqrt();
        let z = (n0 / gradient) * (gradient * distance / n0).asinh();
        let optical_length = 0.5 * (distance * n_end + n0 * z);
        near(path.ray.origin.x, (n_end - n0) / gradient, 2e-11);
        near(path.ray.origin.z, z, 2e-11);
        near(path.final_sample.phase_index() * path.ray.dir.z, n0, 2e-11);
        near(path.optical_path_length_m, optical_length, 2e-11);
        // The achromatic spatial offset changes phase and group equally.
        near(
            path.group_optical_path_length_m.unwrap(),
            optical_length + (2.0 - n0) * distance,
            2e-11,
        );
    }
}

#[test]
fn g1_infrared_bundle_retains_chromatic_focusing() {
    let material = infrared();
    let curvature = -1.5;
    let field = QuadraticIndexField::try_new(
        &material,
        ray().origin,
        Vec3::new(0.0, 0.0, 0.0),
        IndexHessian::try_new([curvature, curvature, 0.0], [0.0; 3]).unwrap(),
        bounds(),
    )
    .unwrap();
    let distance = 0.8;
    for wavelength in [950.0, 1450.0] {
        let path = propagate_bundle(
            &field,
            ray(),
            axis(),
            wavelength,
            distance,
            bundle_options(),
        )
        .unwrap();
        // On the symmetry axis X'' = (curvature/n_phase) X; X(0)=0, X'(0)=1.
        let frequency = (-curvature / (2.0 - wavelength / 3200.0)).sqrt();
        let footprint = (frequency * distance).sin() / frequency;
        let projection = path
            .project_at_endpoint(ray().dir, ProjectionOptions::default())
            .unwrap();
        near(
            projection.area_per_solid_angle_m2(),
            footprint * footprint,
            2e-11,
        );
        near(
            path.propagation().group_optical_path_length_m.unwrap(),
            2.0 * distance,
            2e-11,
        );
    }
}

#[test]
fn g0_refractive_fields_keep_the_common_source_band_and_interior_extrema() {
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (900.0, 0.25),
        (1500.0, 0.75),
    ])
    .unwrap();
    // Unsupported phase knots would cross the index limits after the spatial
    // offset. Only the original COMMON support is a propagation promise.
    let clipped = table(
        |i| {
            if i == 0 {
                1.0
            } else if i == 8 {
                3.0
            } else {
                1.5
            }
        },
        absorption,
    );
    let affine =
        AffineIndexField::try_new(&clipped, ray().origin, Vec3::new(0.2, 0.0, 0.0), bounds())
            .unwrap();
    let quadratic = QuadraticIndexField::try_new(
        &clipped,
        ray().origin,
        Vec3::new(0.0, 0.0, 0.0),
        IndexHessian::try_new([-1.0, -1.0, 0.0], [0.0; 3]).unwrap(),
        bounds(),
    )
    .unwrap();
    for wavelength in [0.0, -1.0, f64::NAN, f64::INFINITY, 850.0, 1550.0] {
        let expected = Err(GradientIndexError::Dielectric(
            DielectricError::InvalidWavelength,
        ));
        assert_eq!(
            propagate(&affine, ray(), axis(), wavelength, 0.1, tight()),
            expected
        );
        assert_eq!(
            propagate(&quadratic, ray(), axis(), wavelength, 0.1, tight()),
            expected
        );
    }
    let interior_peak = table(|i| if i == 4 { 2.99 } else { 1.5 }, absorption);
    assert!(
        AffineIndexField::try_new(
            &interior_peak,
            ray().origin,
            Vec3::new(0.2, 0.0, 0.0),
            bounds(),
        )
        .is_err()
    );
    assert!(
        QuadraticIndexField::try_new(
            &interior_peak,
            ray().origin,
            Vec3::new(0.0, 0.0, 0.0),
            IndexHessian::try_new([1.0, 1.0, 0.0], [0.0; 3]).unwrap(),
            bounds(),
        )
        .is_err()
    );

    // A shared endpoint is valid for a monochromatic ray; separated sources
    // must not be extrapolated to fabricate an overlap.
    for (lower, allowed) in [(1600.0, true), (1700.0, false)] {
        let material = table(
            |_| 1.5,
            BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
                (lower, 0.0),
                (1800.0, 0.0),
            ])
            .unwrap(),
        );
        let field =
            AffineIndexField::try_new(&material, ray().origin, Vec3::new(0.1, 0.0, 0.0), bounds());
        assert_eq!(field.is_ok(), allowed);
        if let Ok(field) = field {
            assert!(propagate(&field, ray(), axis(), 1600.0, 0.1, tight()).is_ok());
        }
    }
    // Existing analytic glass retains its declared visible support.
    let visible = DielectricGlass::new(
        crate::dielectric::CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    assert!(
        propagate(
            &HomogeneousIndexField::new(&visible),
            ray(),
            axis(),
            1200.0,
            0.1,
            tight()
        )
        .is_err()
    );
}

struct Plane(f64);
impl BoundaryOracle for Plane {
    fn query(&mut self, p: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let gap = p.z - self.0;
        Ok(BoundaryQuery {
            clearance_m: gap.abs().next_down().max(0.0),
            witness: Some(BoundaryWitness {
                boundary_id: 0,
                feature_id: 0,
                surface_point: Point3::new(p.x, p.y, self.0),
                normal: ray().dir,
                distance_upper_m: gap.abs().next_up(),
                side: if gap < 0.0 {
                    BoundarySide::Negative
                } else {
                    BoundarySide::Positive
                },
                departure_radius_m: Some(1.0),
            }),
        })
    }
}

#[test]
fn g1_infrared_boundary_and_seeded_bundle_preserve_the_native_path() {
    let material = infrared();
    let field = HomogeneousIndexField::new(&material);
    let options = GradientEventOptions {
        integration: tight(),
        ..GradientEventOptions::default()
    };
    let scalar = propagate_to_boundary(
        &field,
        &mut Plane(0.75),
        ray(),
        axis(),
        1200.0,
        1.0,
        None,
        options,
    )
    .unwrap();
    assert!(matches!(
        scalar.termination,
        BoundaryTermination::Surface(_)
    ));
    near(
        scalar.propagation.ray.origin.z,
        0.75,
        2.0 * options.surface_tolerance_m,
    );
    let seed = propagate_bundle(&field, ray(), axis(), 1200.0, 0.25, bundle_options()).unwrap();
    let result = propagate_seeded_bundle_to_boundary(
        &field,
        &mut Plane(0.75),
        &seed.continuation(),
        0.75,
        None,
        BundleEventOptions {
            boundary: options,
            position_derivative_tolerance_m: 1e-12,
            direction_derivative_tolerance: 1e-12,
        },
    )
    .unwrap();
    assert!(matches!(
        result.termination,
        BoundaryTermination::Surface(_)
    ));
    let path = result.bundle.propagation();
    near(path.ray.origin.z, 0.75, 2.0 * options.surface_tolerance_m);
    near(
        path.group_optical_path_length_m.unwrap(),
        2.0 * path.distance_m,
        2e-14,
    );
    near(path.optical_path_length_m, 1.625 * path.distance_m, 2e-14);
    let projected = result
        .bundle
        .project_at_endpoint(ray().dir, ProjectionOptions::default())
        .unwrap();
    near(
        projected.area_per_solid_angle_m2(),
        path.distance_m * path.distance_m,
        2e-14,
    );
}

#[test]
fn g1_infrared_polarized_beam_uses_phase_refraction_and_group_delay() {
    let ambient = table(|_| 1.0, BeerLambertAbsorption::CLEAR);
    let glass = infrared();
    let field = HomogeneousIndexField::new(&glass);
    let medium = BeamMedium::try_gradient_at(&field, ray().origin, 1200.0).unwrap();
    let beam = PolarizedBeam::try_new(
        ray().origin,
        PolarizationFrame::try_new(Vec3::new(0.6, 0.0, 0.8), Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        1200.0,
        ambient,
        Stokes::unpolarized(1.0),
    )
    .unwrap();
    let transmitted = beam
        .sample_smooth_interface(ray().dir, medium, 0.99)
        .unwrap();
    let reflected = beam
        .sample_smooth_interface(ray().dir, medium, 0.0)
        .unwrap();
    assert_eq!(transmitted.event, DielectricEvent::Transmission);
    near(transmitted.beam.ray().dir.x, 0.6 / 1.625, 2e-14);
    near(
        transmitted.pdf * transmitted.beam.power().i + reflected.pdf * reflected.beam.power().i,
        1.0,
        2e-14,
    );
    let advanced = transmitted
        .beam
        .advance_gradient(&field, 0.5, tight())
        .unwrap();
    near(
        advanced.propagation.optical_path_length_m,
        0.5 * 1.625,
        2e-14,
    );
    near(
        advanced.propagation.group_optical_path_length_m.unwrap(),
        1.0,
        2e-14,
    );
    near(
        advanced.beam.power().i,
        transmitted.beam.power().i * (-0.25_f64).exp(),
        2e-14,
    );
    assert_eq!(advanced.beam.ray().dir, transmitted.beam.ray().dir);
}
