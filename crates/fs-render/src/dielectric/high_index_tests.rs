use super::*;
use crate::charts::Ray;
use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorSource};
use crate::gradient_index::bundle::{
    BundleError, BundleOptions, IndexHessian, ProjectionOptions, QuadraticIndexField,
    propagate_bundle,
};
use crate::gradient_index::{AffineIndexField, GradientIndexOptions, propagate};
use crate::polarization::Stokes;
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use fs_blake3::ContentHash;
use fs_geom::{Aabb, Point3};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn constant_index(index: f64) -> SellmeierIor {
    SellmeierIor::try_new_in_band([index * index - 1.0, 0.0, 0.0], [0.0; 3], (1200.0, 1800.0))
        .unwrap()
}

#[test]
fn explicit_spectral_laws_admit_high_indices_and_preserve_visible_constructor_domains() {
    let model = constant_index(4.0);
    assert_eq!(model.eval(1550.0), Ok(4.0));
    assert_eq!(model.group_index(1550.0), Ok(4.0));
    assert_eq!(model.eval(1199.0), Err(DielectricError::InvalidWavelength));
    let samples = core::array::from_fn(|i| {
        DielectricIorSample::try_new_in_band(1200.0 + 75.0 * i as f64, 4.0).unwrap()
    });
    let table = TabulatedIor::try_new_in_band(samples, 75.0).unwrap();
    assert_eq!(table.eval(1550.0), Ok(4.0));
    assert_eq!(table.group_index(1550.0), Ok(4.0));
    assert_eq!(
        CauchyIor::try_constant(4.0),
        Err(DielectricError::InvalidIor)
    );
    assert_eq!(
        SellmeierIor::try_new([15.0, 0.0, 0.0], [0.0; 3]),
        Err(DielectricError::InvalidIor)
    );
    assert_eq!(
        DielectricIorSample::try_new(550.0, 4.0),
        Err(DielectricError::InvalidIor)
    );
    let visible_knots = core::array::from_fn(|i| {
        DielectricIorSample::try_new_in_band(380.0 + 50.0 * i as f64, 4.0).unwrap()
    });
    assert_eq!(
        TabulatedIor::try_new(visible_knots),
        Err(DielectricError::InvalidIor)
    );
    for index in [0.999, 1.0e154, f64::INFINITY, f64::NAN] {
        assert_eq!(
            DielectricIorSample::try_new_in_band(1550.0, index),
            Err(DielectricError::InvalidIor)
        );
    }
}

#[test]
fn high_index_interfaces_keep_snell_power_reciprocity_brewster_and_critical_behavior() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    for index in [3.5_f64, 4.0, 9.0] {
        let ct = (1.0 - (0.6 / index).powi(2)).sqrt();
        let rs = (0.8 - index * ct) / (0.8 + index * ct);
        let rp = (index * 0.8 - ct) / (index * 0.8 + ct);
        let expected = 0.5 * (rs * rs + rp * rp);
        let scalar = fresnel_dielectric(0.8, 1.0, index).unwrap();
        near(scalar.reflectance, expected, 3.0e-15);
        near(scalar.reflectance + scalar.transmittance, 1.0, 3.0e-15);
        let reversed = fresnel_dielectric(ct, index, 1.0).unwrap();
        near(reversed.reflectance, expected, 3.0e-14);
        near(reversed.transmitted_cosine, 0.8, 3.0e-14);
        let ray = sample_smooth_dielectric(normal, wo, 1.0, index, 0.99).unwrap();
        assert_eq!(ray.event, DielectricEvent::Transmission);
        near(ray.direction.x, -0.6 / index, 1.0e-15);
        near(ray.direction.z, -ct, 1.0e-15);
        near(ray.radiance_weight, 1.0 / (index * index), 1.0e-15);
        let polarized = fresnel_dielectric_polarized(0.8, 1.0, index).unwrap();
        near(
            polarized.reflection.apply(Stokes::unpolarized(1.0)).i,
            expected,
            3.0e-15,
        );
        let brewster =
            fresnel_dielectric_polarized(1.0 / (1.0 + index * index).sqrt(), 1.0, index).unwrap();
        let p = Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap();
        near(brewster.reflection.apply(p).i, 0.0, 3.0e-15);
        let critical = (1.0 - 1.0 / (index * index)).sqrt();
        assert!(
            fresnel_dielectric(critical - 2.0e-15, index, 1.0)
                .unwrap()
                .total_internal_reflection
        );
        assert!(
            !fresnel_dielectric(critical + 2.0e-15, index, 1.0)
                .unwrap()
                .total_internal_reflection
        );
    }
    let normal_power = fresnel_dielectric(1.0, 1.0, 4.0).unwrap();
    near(normal_power.reflectance, 9.0 / 25.0, 3.0e-16);
    near(normal_power.transmittance, 16.0 / 25.0, 3.0e-16);
}

#[test]
fn high_index_rough_transmission_retains_radiance_reciprocity() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    for (u, v) in [(0.03, 0.2), (0.3, 0.7), (0.85, 0.4)] {
        let sample = sample_rough_dielectric(normal, wo, 1.0, 4.0, 0.2, u, v, 0.99)
            .unwrap()
            .unwrap();
        assert_eq!(sample.event, DielectricEvent::Transmission);
        let forward =
            evaluate_rough_dielectric(normal, wo, sample.direction, 1.0, 4.0, 0.2).unwrap();
        let reverse =
            evaluate_rough_dielectric(normal.scale(-1.0), sample.direction, wo, 4.0, 1.0, 0.2)
                .unwrap();
        // Radiance BTDF reciprocity: f_forward * n_t² = f_reverse * n_i².
        near(
            forward.value * 16.0,
            reverse.value,
            2.0e-11 * reverse.value.max(1.0),
        );
        assert!(forward.pdf.is_finite() && forward.pdf > 0.0);
        assert!(sample.radiance_weight.is_finite() && sample.radiance_weight >= 0.0);
    }
}

#[test]
fn high_index_coherent_films_keep_quarter_wave_matching_and_reverse_power() {
    let source =
        ConductorSource::try_new(ContentHash([84; 32]), ConductorDataStatus::Representative)
            .unwrap();
    for (film_index, exit_index) in [(2.0_f64, 4.0_f64), (4.0, 2.0)] {
        let films = [CoherentDielectricCoating::try_new(
            constant_index(film_index),
            1550.0e-9 / (4.0 * film_index),
            source,
        )
        .unwrap()];
        let stack = DielectricCoatingStack::try_new_in_band(&films, (1200.0, 1800.0)).unwrap();
        // Independent quarter-wave admittance formula at normal incidence.
        let expected = ((exit_index - film_index * film_index)
            / (exit_index + film_index * film_index))
            .powi(2);
        for (reverse, incident, exit) in [(false, 1.0, exit_index), (true, exit_index, 1.0)] {
            let balance = stack
                .layer_power_balance(1550.0, reverse, 1.0, incident, exit)
                .unwrap();
            near(balance.fresnel.scalar.reflectance, expected, 2.0e-14);
            near(
                balance.fresnel.scalar.transmittance,
                1.0 - expected,
                2.0e-14,
            );
            assert_eq!(balance.absorbed, [0.0; 2]);
            assert_eq!(balance.layer_absorbed(), &[[0.0; 2]]);
        }
    }
}

#[test]
fn high_index_gradient_field_matches_closed_form_bending_phase_group_and_beer() {
    let glass = DielectricGlass::new(
        constant_index(4.0),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    );
    let origin = Point3::new(0.0, 0.0, 0.0);
    let field = AffineIndexField::try_new(
        &glass,
        origin,
        Vec3::new(0.0, 0.4, 0.0),
        Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0)),
    )
    .unwrap();
    let length = 0.7;
    let path = propagate(
        &field,
        Ray {
            origin,
            dir: Vec3::new(1.0, 0.0, 0.0),
        },
        Vec3::new(0.0, 0.0, 1.0),
        1550.0,
        length,
        GradientIndexOptions {
            max_step_m: 0.03,
            position_tolerance_m: 1.0e-12,
            direction_tolerance: 1.0e-12,
            optical_path_tolerance_m: 1.0e-12,
            optical_depth_tolerance: 1.0e-12,
            relative_tolerance: 1.0e-12,
            ..GradientIndexOptions::default()
        },
    )
    .unwrap();
    // n(y)=4+.4y: conserved n*dx/ds=4, and d(n*dy/ds)/ds=.4.
    let index = (16.0_f64 + (0.4 * length).powi(2)).sqrt();
    let integral = (0.4 * length / 4.0).asinh();
    near(path.ray.origin.x, 4.0 / 0.4 * integral, 2.0e-11);
    near(path.ray.origin.y, (index - 4.0) / 0.4, 2.0e-11);
    near(path.ray.dir.x, 4.0 / index, 2.0e-11);
    near(path.ray.dir.y, 0.4 * length / index, 2.0e-11);
    let optical_length = 0.5 * (length * index + 16.0 / 0.4 * integral);
    near(path.optical_path_length_m, optical_length, 2.0e-11);
    near(
        path.group_optical_path_length_m.unwrap(),
        optical_length,
        2.0e-11,
    );
    near(path.power_transmittance(), (-0.2 * length).exp(), 2.0e-12);
}

#[test]
fn admitted_extreme_contrast_retains_normal_transmission_and_finite_tir_phase() {
    let index = 1.0e100;
    let normal = fresnel_dielectric(1.0, index, 1.0).unwrap();
    assert!(!normal.total_internal_reflection);
    assert_eq!(normal.transmitted_cosine, 1.0);
    near(normal.transmittance / (4.0 / index), 1.0, 3.0e-15);
    assert_eq!(
        refract(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            index,
            1.0,
            1.0
        )
        .unwrap(),
        Vec3::new(0.0, 0.0, -1.0)
    );
    let polarized = fresnel_dielectric_polarized(0.5, index, 1.0).unwrap();
    assert!(polarized.scalar.total_internal_reflection);
    assert!(
        polarized
            .reflection
            .0
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );
    let result = polarized
        .reflection
        .apply(Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap());
    near(result.i, 1.0, 1.0e-15);
    near(result.degree(), 1.0, 1.0e-15);
}

#[test]
fn high_index_quadratic_lenses_retain_exact_focusing_and_physical_path_integrals() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let axis = Vec3::new(0.0, 0.0, 1.0);
    let distance = 0.6_f64;
    for index in [4.0, 9.0] {
        let glass = DielectricGlass::new(
            constant_index(index),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        for sign in [-1.0, 1.0] {
            let curvature = sign * 4.0 * index;
            let field = QuadraticIndexField::try_new(
                &glass,
                origin,
                Vec3::new(0.0, 0.0, 0.0),
                IndexHessian::try_new([curvature, curvature, 0.0], [0.0; 3]).unwrap(),
                Aabb::new(Point3::new(-0.1, -0.1, -0.1), Point3::new(0.1, 0.1, 1.0)),
            )
            .unwrap();
            let ray = Ray { origin, dir: axis };
            let options = BundleOptions {
                integration: GradientIndexOptions {
                    position_tolerance_m: 1.0e-12,
                    direction_tolerance: 1.0e-12,
                    optical_path_tolerance_m: 1.0e-12,
                    optical_depth_tolerance: 1.0e-12,
                    relative_tolerance: 1.0e-12,
                    ..GradientIndexOptions::default()
                },
                position_derivative_tolerance_m: 1.0e-12,
                direction_derivative_tolerance: 1.0e-12,
            };
            let bundle = propagate_bundle(
                &field,
                ray,
                Vec3::new(1.0, 0.0, 0.0),
                1550.0,
                distance,
                options,
            )
            .unwrap();
            // On axis n is constant, but the exact variations obey X''=+/-4X.
            // This compares the full differential integrator with the analytic
            // ray-tube solution, independently of finitely separated test rays.
            let (position, direction) = if sign < 0.0 {
                (0.5 * (2.0 * distance).sin(), (2.0 * distance).cos())
            } else {
                (0.5 * (2.0 * distance).sinh(), (2.0 * distance).cosh())
            };
            for (column, basis) in bundle
                .angular_derivatives()
                .into_iter()
                .zip(bundle.launch_basis())
            {
                near(column.position_m_per_rad.dot(basis), position, 3.0e-11);
                near(column.direction_per_rad.dot(basis), direction, 3.0e-11);
            }
            let path = bundle.propagation();
            near(path.ray.origin.z, distance, 2.0e-14);
            near(path.optical_path_length_m, index * distance, 3.0e-14);
            near(
                path.group_optical_path_length_m.unwrap(),
                index * distance,
                3.0e-14,
            );
            near(path.optical_depth, 0.2 * distance, 2.0e-15);
            let projection = bundle
                .project_at_endpoint(axis, ProjectionOptions::default())
                .unwrap();
            near(
                projection.area_per_solid_angle_m2(),
                position * position,
                4.0e-11,
            );
            near(
                projection.point_source_irradiance(2.0).unwrap(),
                2.0 * (-0.2 * distance).exp() / (position * position),
                3.0e-10,
            );
            assert!(
                propagate_bundle(
                    &field,
                    ray,
                    Vec3::new(1.0, 0.0, 0.0),
                    1100.0,
                    distance,
                    options,
                )
                .is_err()
            );
        }
    }
}

#[test]
fn quadratic_index_enclosure_still_refuses_nonphysical_or_unrepresentable_regions() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let bounds = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
    for (index, curvature) in [(4.0, -8.0), (1.0e153, 1.0e155)] {
        let glass = DielectricGlass::new(
            constant_index(index),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        );
        assert!(matches!(
            QuadraticIndexField::try_new(
                &glass,
                origin,
                Vec3::new(0.0, 0.0, 0.0),
                IndexHessian::try_new([curvature, 0.0, 0.0], [0.0; 3]).unwrap(),
                bounds,
            ),
            Err(BundleError::InvalidDifferential)
        ));
    }
}
