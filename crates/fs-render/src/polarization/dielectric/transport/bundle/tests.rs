use super::*;
use crate::charts::Ray;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, GlassProvenance,
};
use crate::gradient_index::bundle::interface::{BundleInterfaceError, InterfaceOptions};
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexField, HomogeneousIndexField};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Chart, NormalDifferential, Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane_normal(normal: Vec3) -> NormalDifferential {
    NormalDifferential {
        normal,
        derivative_per_m: [[0.0; 3]; 3],
    }
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn reference_refraction(direction: Vec3, outward: Vec3, ratio: f64) -> Vec3 {
    let normal = if direction.dot(outward) < 0.0 {
        outward
    } else {
        outward.scale(-1.0)
    };
    let cosine = -direction.dot(normal);
    let transmitted = (1.0 - ratio * ratio * (1.0 - cosine * cosine)).sqrt();
    unit(add(
        direction.scale(ratio),
        normal.scale(ratio * cosine - transmitted),
    ))
}

#[derive(Clone, Copy)]
struct BallReference {
    endpoint: Ray,
    before_length: f64,
    inside_length: f64,
    after_length: f64,
    outside_length: f64,
    transmission_s: f64,
    transmission_p: f64,
}

/// Independent quadratic sphere intersections and vector Snell calculation.
/// The production event walk, normal Jacobian and bundle solver are not used.
fn reference_ball(ray: Ray, index: f64, receiver_z: f64) -> BallReference {
    let origin = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z);
    let b = origin.dot(ray.dir);
    let before = -b - (b * b - origin.dot(origin) + 1.0).sqrt();
    assert!(before > 0.0);
    let entry = ray.at(before);
    let entry_normal = unit(Vec3::new(entry.x, entry.y, entry.z));
    let inside_direction = reference_refraction(ray.dir, entry_normal, 1.0 / index);
    let inside_length = -2.0 * inside_direction.dot(Vec3::new(entry.x, entry.y, entry.z));
    let exit = entry.offset(inside_direction.scale(inside_length));
    let exit_normal = unit(Vec3::new(exit.x, exit.y, exit.z));
    let direction = reference_refraction(inside_direction, exit_normal, index);
    let after = (receiver_z - exit.z) / direction.z;
    let ci = -ray.dir.dot(entry_normal);
    let ct = -inside_direction.dot(entry_normal);
    let rs = ((ci - index * ct) / (ci + index * ct)).powi(2);
    let rp = ((index * ci - ct) / (index * ci + ct)).powi(2);
    BallReference {
        endpoint: Ray {
            origin: exit.offset(direction.scale(after)),
            dir: direction,
        },
        before_length: before,
        inside_length,
        after_length: after,
        outside_length: before + after,
        transmission_s: 1.0 - rs,
        transmission_p: 1.0 - rp,
    }
}

fn launch_basis(direction: Vec3) -> [Vec3; 2] {
    let first = unit(add(Vec3::new(1.0, 0.0, 0.0), direction.scale(-direction.x)));
    [first, cross(direction, first)]
}

fn reference_ball_columns(ray: Ray, index: f64, receiver_z: f64) -> [Vec3; 2] {
    let h = 1e-6;
    launch_basis(ray.dir).map(|axis| {
        let perturbed = |angle| Ray {
            origin: ray.origin,
            dir: unit(add(ray.dir, axis.scale(angle))),
        };
        let plus = reference_ball(perturbed(h), index, receiver_z)
            .endpoint
            .origin;
        let minus = reference_ball(perturbed(-h), index, receiver_z)
            .endpoint
            .origin;
        plus.delta_from(minus).scale(0.5 / h)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tolerance {tolerance:.3e})"
    );
}

fn through_ball(
    cx: &Cx<'_>,
    ray: Ray,
    material: &DielectricGlass,
    wavelength: f64,
    receiver_z: f64,
) -> PolarizedRayBundle {
    let reference = reference_ball(ray, material.ior().eval(wavelength).unwrap(), receiver_z);
    let air = glass(1.0, 0.0);
    let outside = HomogeneousIndexField::new(&air);
    let inside = HomogeneousIndexField::new(material);
    let sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let options = BundleOptions::default();
    let first = propagate_bundle(
        &outside,
        ray,
        launch_basis(ray.dir)[0],
        wavelength,
        reference.before_length,
        options,
    )
    .unwrap();
    let mut path = PolarizedRayBundle::unpolarized(first).unwrap();
    let entry = path.propagation().propagation().ray.origin;
    path = path
        .smooth_interface(
            sphere.trace_normal_differential(entry, cx).unwrap(),
            inside.sample(entry, wavelength).unwrap(),
            DielectricEvent::Transmission,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&inside, reference.inside_length, options)
        .unwrap();
    let exit = path.propagation().propagation().ray.origin;
    path.smooth_interface(
        sphere.trace_normal_differential(exit, cx).unwrap(),
        outside.sample(exit, wavelength).unwrap(),
        DielectricEvent::Transmission,
        InterfaceOptions::default(),
    )
    .unwrap()
    .continue_in(&outside, reference.after_length, options)
    .unwrap()
}

#[test]
#[cfg(feature = "photometry")]
fn g1_dispersive_thick_ball_lens_matches_abcd_and_calibrated_spectral_irradiance() {
    use crate::instances::RigidTransform;
    use crate::photometry::source::{
        AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, PhotometricSource,
        SPECTRAL_BINS, SourceSpectrum, photopic_weight_nm, wavelength_nm,
    };
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, -3.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let lens = DielectricGlass::new(
            CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.15).unwrap(),
            GlassProvenance::Custom,
        );
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -3.0]).unwrap(),
            1200.0,
            0.75,
            AngularDistribution::isotropic(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let receiver_normal = Vec3::new(0.0, 0.0, -1.0);
        let mut actual_integrals = [0.0; 2];
        let mut expected_integrals = [0.0; 2];
        let mut spectral_endpoints = [0.0; 2];
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let path = through_ball(cx, ray, &lens, wavelength, 2.5);
            let geometric = path.propagation();
            assert_eq!(geometric.interface_events(), 2);
            assert_eq!(geometric.initial_ray(), ray);
            close(geometric.propagation().ray.origin.z, 2.5, 2e-12);
            // In (x,n*theta), free propagation is [[1,d/n],[0,1]] and a
            // spherical interface is [[1,0],[-(n_out-n_in)/R,1]]. With
            // front/back R=+1/-1, 2 m before, 2 m glass and 1.5 m after,
            // the complete source-angle column is B=15/n-9.5 in both axes.
            let b = 15.0 / index - 9.5;
            let projection = geometric
                .project_at_endpoint(receiver_normal, ProjectionOptions::default())
                .unwrap();
            close(projection.area_per_solid_angle_m2(), b * b, 2e-10);
            let columns = projection.angular_derivatives();
            close(columns[0].position_m_per_rad.x, b, 2e-10);
            close(columns[1].position_m_per_rad.y, b, 2e-10);
            close(geometric.propagation().distance_m, 5.5, 3e-12);
            close(
                geometric.propagation().optical_path_length_m,
                3.5 + 2.0 * index,
                3e-12,
            );
            let group_index = 1.45 + 0.036 / (wavelength * 1e-3).powi(2);
            close(
                geometric.propagation().group_optical_path_length_m.unwrap(),
                3.5 + 2.0 * group_index,
                3e-12,
            );
            close(geometric.propagation().optical_depth, 0.3, 3e-13);
            let transmission = 4.0 * index / (1.0 + index).powi(2);
            close(path.interface_log_power(), 2.0 * transmission.ln(), 2e-12);
            let intensity = source
                .spectral_intensity_w_sr_nm(ray.dir, wavelength)
                .unwrap();
            let actual = path
                .point_source_irradiance(receiver_normal, ProjectionOptions::default(), intensity)
                .unwrap();
            let expected = intensity * transmission.powi(2) * (-0.3_f64).exp() / (b * b);
            close(actual, expected, 3e-9 * expected);
            actual_integrals[0] += actual * BIN_WIDTH_NM;
            actual_integrals[1] += actual * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
            expected_integrals[0] += expected * BIN_WIDTH_NM;
            expected_integrals[1] += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
            if bin == 0 {
                spectral_endpoints[0] = actual;
            }
            if bin == SPECTRAL_BINS - 1 {
                spectral_endpoints[1] = actual;
            }
        }
        for (actual, expected) in actual_integrals.into_iter().zip(expected_integrals) {
            close(actual, expected, 3e-9 * expected);
        }
        assert!(
            spectral_endpoints[0] > 4.0 * spectral_endpoints[1],
            "chromatic focusing must reach physical irradiance"
        );
    });
}

#[test]
fn g1_off_axis_ball_lens_derivatives_match_independent_sphere_snell_perturbations() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let lens = glass(1.5, 0.17);
        let expected = reference_ball(ray, 1.5, 2.5);
        let expected_columns = reference_ball_columns(ray, 1.5, 2.5);
        let path = through_ball(cx, ray, &lens, 550.0, 2.5);
        let projection = path
            .propagation()
            .project_at_endpoint(Vec3::new(0.0, 0.0, -1.0), ProjectionOptions::default())
            .unwrap();
        let endpoint = path.propagation().propagation().ray;
        close(endpoint.origin.x, expected.endpoint.origin.x, 3e-12);
        close(endpoint.origin.y, expected.endpoint.origin.y, 3e-12);
        close(endpoint.origin.z, expected.endpoint.origin.z, 3e-12);
        for (actual, expected) in projection
            .angular_derivatives()
            .into_iter()
            .zip(expected_columns)
        {
            close(actual.position_m_per_rad.x, expected.x, 2e-7);
            close(actual.position_m_per_rad.y, expected.y, 2e-7);
            close(actual.position_m_per_rad.z, expected.z, 2e-7);
        }
        let area = cross(expected_columns[0], expected_columns[1]).z.abs();
        close(projection.area_per_solid_angle_m2(), area, 2e-7);
        assert!(
            (area - 0.25).abs() > 0.1,
            "off-axis aberration must survive both normal Jacobians"
        );
        close(
            path.propagation().propagation().distance_m,
            expected.inside_length + expected.outside_length,
            3e-12,
        );
        // The first oblique face polarizes an initially unpolarized source.
        // True transmitted power at the second face is mean(Ts²,Tp²),
        // not the square of the unpolarized mean transmission.
        let throughput = 0.5 * (expected.transmission_s.powi(2) + expected.transmission_p.powi(2));
        close(path.interface_log_power().exp(), throughput, 3e-12);
        let irradiance = path
            .point_source_irradiance(Vec3::new(0.0, 0.0, -1.0), ProjectionOptions::default(), 1.0)
            .unwrap();
        close(
            irradiance,
            throughput * (-0.17 * expected.inside_length).exp() / area,
            3e-6,
        );
    });
}

#[test]
fn g1_plane_parallel_slab_retains_refractive_bundle_scale_and_true_transmitted_power() {
    let air = glass(1.0, 0.0);
    let slab = glass(1.5, 0.2);
    let outside = HomogeneousIndexField::new(&air);
    let inside = HomogeneousIndexField::new(&slab);
    let options = BundleOptions::default();
    let wavelength = 550.0;
    let launch = Ray {
        origin: Point3::new(0.0, 0.0, -0.5),
        dir: Vec3::new(0.0, 0.0, 1.0),
    };
    let first = propagate_bundle(
        &outside,
        launch,
        Vec3::new(1.0, 0.0, 0.0),
        wavelength,
        0.5,
        options,
    )
    .unwrap();
    let path = PolarizedRayBundle::unpolarized(first)
        .unwrap()
        .smooth_interface(
            plane_normal(Vec3::new(0.0, 0.0, -1.0)),
            inside
                .sample(Point3::new(0.0, 0.0, 0.0), wavelength)
                .unwrap(),
            DielectricEvent::Transmission,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&inside, 0.4, options)
        .unwrap();
    let exit = path.propagation().propagation().ray.origin;
    let path = path
        .smooth_interface(
            plane_normal(Vec3::new(0.0, 0.0, 1.0)),
            outside.sample(exit, wavelength).unwrap(),
            DielectricEvent::Transmission,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&outside, 0.7, options)
        .unwrap();
    let normal = Vec3::new(0.0, 0.0, -1.0);
    let area = (0.5_f64 + 0.4 / 1.5 + 0.7).powi(2);
    close(
        path.propagation()
            .project_at_endpoint(normal, ProjectionOptions::default())
            .unwrap()
            .area_per_solid_angle_m2(),
        area,
        2e-12,
    );
    close(path.propagation().propagation().ray.origin.z, 1.1, 2e-13);
    close(
        path.point_source_irradiance(normal, ProjectionOptions::default(), 1.0)
            .unwrap(),
        0.96_f64.powi(2) * (-0.08_f64).exp() / area,
        3e-12,
    );
}

#[test]
fn g1_tir_reflection_retains_folded_bundle_area_and_refuses_transmission() {
    let dense = glass(1.5, 0.1);
    let air = glass(1.0, 0.0);
    let field = HomogeneousIndexField::new(&dense);
    let target = HomogeneousIndexField::new(&air);
    let launch = Ray {
        origin: Point3::new(0.0, 0.0, -1.0),
        dir: Vec3::new(0.8, 0.0, 0.6),
    };
    let first = propagate_bundle(
        &field,
        launch,
        launch_basis(launch.dir)[0],
        550.0,
        1.0 / 0.6,
        BundleOptions::default(),
    )
    .unwrap();
    let path = PolarizedRayBundle::unpolarized(first).unwrap();
    let endpoint = path.propagation().propagation().ray.origin;
    let sample = target.sample(endpoint, 550.0).unwrap();
    let normal = plane_normal(Vec3::new(0.0, 0.0, 1.0));
    assert!(matches!(
        PolarizedRayBundle::unpolarized(first)
            .unwrap()
            .smooth_interface(
                normal,
                sample,
                DielectricEvent::Transmission,
                InterfaceOptions::default()
            ),
        Err(PolarizedBundleError::Interface(
            BundleInterfaceError::TotalInternalReflection
        ))
    ));
    let reflected = path
        .smooth_interface(
            normal,
            sample,
            DielectricEvent::Reflection,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&field, 1.0 / 0.6, BundleOptions::default())
        .unwrap();
    let receiver_normal = Vec3::new(0.0, 0.0, 1.0);
    let length = 2.0 / 0.6;
    let area = length * length / 0.6;
    close(reflected.propagation().propagation().ray.dir.x, 0.8, 2e-13);
    close(reflected.propagation().propagation().ray.dir.z, -0.6, 2e-13);
    close(
        reflected
            .propagation()
            .project_at_endpoint(receiver_normal, ProjectionOptions::default())
            .unwrap()
            .area_per_solid_angle_m2(),
        area,
        3e-11,
    );
    close(reflected.interface_log_power(), 0.0, 2e-13);
    close(
        reflected
            .point_source_irradiance(receiver_normal, ProjectionOptions::default(), 1.0)
            .unwrap(),
        (-0.1 * length).exp() / area,
        3e-12,
    );
}

#[test]
fn g0_near_critical_transmission_refuses_the_singular_angular_derivative() {
    let dense = glass(1.5, 0.0);
    let air = glass(1.0, 0.0);
    let field = HomogeneousIndexField::new(&dense);
    let target = HomogeneousIndexField::new(&air);
    // A supported transmitted central ray, but cos(theta_t) ≈1.4e-7 is
    // below the explicit 1e-6 differential admission limit.
    let sine = (1.0_f64 / 1.5) * (1.0 - 1e-14);
    let cosine = (1.0 - sine * sine).sqrt();
    let launch = Ray {
        origin: Point3::new(0.0, 0.0, -1.0),
        dir: Vec3::new(sine, 0.0, cosine),
    };
    let first = propagate_bundle(
        &field,
        launch,
        launch_basis(launch.dir)[0],
        550.0,
        1.0 / cosine,
        BundleOptions::default(),
    )
    .unwrap();
    let path = PolarizedRayBundle::unpolarized(first).unwrap();
    let sample = target
        .sample(path.propagation().propagation().ray.origin, 550.0)
        .unwrap();
    assert!(matches!(
        path.smooth_interface(
            plane_normal(Vec3::new(0.0, 0.0, 1.0)),
            sample,
            DielectricEvent::Transmission,
            InterfaceOptions::default()
        ),
        Err(PolarizedBundleError::Interface(
            BundleInterfaceError::CriticalTransmission
        ))
    ));
    // Reflection has no transmitted-cosine denominator and remains regular.
    let reflected = path
        .smooth_interface(
            plane_normal(Vec3::new(0.0, 0.0, 1.0)),
            sample,
            DielectricEvent::Reflection,
            InterfaceOptions::default(),
        )
        .unwrap();
    close(
        reflected.propagation().propagation().ray.dir.z,
        -cosine,
        2e-13,
    );
    assert!(reflected.interface_log_power().is_finite());
    assert!(reflected.interface_log_power() < 0.0);
}

#[test]
fn g1_quarter_wave_coating_reflects_between_equal_index_bulk_media() {
    use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorSource};
    let material = glass(1.0, 0.0);
    let field = HomogeneousIndexField::new(&material);
    let launch = Ray {
        origin: Point3::new(0.0, 0.0, -1.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    };
    let wavelength = 550.0;
    let first = propagate_bundle(
        &field,
        launch,
        Vec3::new(1.0, 0.0, 0.0),
        wavelength,
        1.0,
        BundleOptions::default(),
    )
    .unwrap();
    let path =
        PolarizedRayBundle::try_new(first, Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap()).unwrap();
    let target = field
        .sample(Point3::new(0.0, 0.0, 0.0), wavelength)
        .unwrap();
    let surface = plane_normal(Vec3::new(0.0, 0.0, -1.0));
    assert!(matches!(
        path.smooth_interface(
            surface,
            target,
            DielectricEvent::Reflection,
            InterfaceOptions::default()
        ),
        Err(PolarizedBundleError::Interface(
            BundleInterfaceError::UnavailableBranch
        ))
    ));
    let source = ConductorSource::try_new(
        fs_blake3::hash_domain("interface-bundle-test", b"explicit-quarter-wave-film"),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let film_index = 1.5;
    let coating = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(film_index).unwrap(),
        wavelength * 1e-9 / (4.0 * film_index),
        source,
    )
    .unwrap();
    let reflected = path
        .smooth_coated_interface(
            surface,
            target,
            &coating,
            DielectricEvent::Reflection,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&field, 1.0, BundleOptions::default())
        .unwrap();
    // At quarter wave between identical unit-index bulk media, the coherent
    // field amplitude is (1-n_f²)/(1+n_f²). Bare bulk Fresnel is exactly zero.
    let reflectance = ((1.0 - film_index * film_index) / (1.0 + film_index * film_index)).powi(2);
    close(reflected.interface_log_power().exp(), reflectance, 2e-12);
    close(reflected.interface_state().i, 1.0, 2e-13);
    close(reflected.interface_state().v, -1.0, 2e-12);
    close(reflected.propagation().propagation().distance_m, 2.0, 2e-13);
    close(
        reflected.propagation().propagation().optical_path_length_m,
        2.0,
        2e-13,
    );
    close(
        reflected
            .point_source_irradiance(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default(), 1.0)
            .unwrap(),
        reflectance / 4.0,
        2e-12,
    );
}

#[test]
fn g1_p_polarized_brewster_branch_is_dark_without_invalid_normalized_stokes() {
    let air = glass(1.0, 0.0);
    let lens = glass(1.5, 0.0);
    let field = HomogeneousIndexField::new(&air);
    let target = HomogeneousIndexField::new(&lens);
    let cosine = 1.0 / 3.25_f64.sqrt();
    let launch = Ray {
        origin: Point3::new(0.0, 0.0, -1.0),
        dir: Vec3::new(1.5 * cosine, 0.0, cosine),
    };
    // +Y is the s direction, hence Q=-I is purely p polarized.
    let first = propagate_bundle(
        &field,
        launch,
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0 / cosine,
        BundleOptions::default(),
    )
    .unwrap();
    let path =
        PolarizedRayBundle::try_new(first, Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap()).unwrap();
    let endpoint = path.propagation().propagation().ray.origin;
    let dark = path
        .smooth_interface(
            plane_normal(Vec3::new(0.0, 0.0, -1.0)),
            target.sample(endpoint, 550.0).unwrap(),
            DielectricEvent::Reflection,
            InterfaceOptions::default(),
        )
        .unwrap()
        .continue_in(&field, 1.0 / cosine, BundleOptions::default())
        .unwrap();
    assert_eq!(dark.interface_log_power(), f64::NEG_INFINITY);
    let state = dark.interface_state();
    assert!(Stokes::try_new(state.i, state.q, state.u, state.v).is_some());
    assert_eq!(state.i, 1.0);
    assert_eq!(
        dark.point_source_irradiance(
            Vec3::new(0.0, 0.0, 1.0),
            ProjectionOptions::default(),
            1e200,
        )
        .unwrap(),
        0.0
    );
    assert_eq!(path.interface_log_power(), 0.0);
    assert_eq!(
        path.interface_state(),
        Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap()
    );
}
