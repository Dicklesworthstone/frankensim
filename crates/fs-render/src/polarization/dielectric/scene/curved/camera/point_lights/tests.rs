use super::super::{CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter};
use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::AffineIndexField;
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::bundle::{IndexHessian, QuadraticIndexField};
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum,
};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::spectral::SpectralReflectance;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Aabb, Point3, Vec3};

mod rough;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, z: f64, center_x: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("curved-point-fixture-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center_x - radius, -radius, z],
                [center_x + radius, -radius, z],
                [center_x, radius, z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn diffuse(object_id: u64, reflectance: usize) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::DepolarizingDiffuse {
            medium: 0,
            reflectance,
        },
    }
}

fn opaque(object_id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Opaque,
    }
}

fn emitter(object_id: u64, radiance: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        object_id,
        RadianceSpectrum::try_constant(radiance).unwrap(),
        EmissionSides::Both,
    )
}

fn source(position: Point3) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [position.x, position.y, position.z])
            .unwrap(),
        100.0,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
    )
    .unwrap()
}

// Isotropic calibration derived from retained lumens and the normalized source
// spectrum, independently of the connector's direction-dependent intensity API.
fn intensity(source: &PhotometricSource) -> f64 {
    source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[0] / (4.0 * core::f64::consts::PI)
}

fn trace_settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 6.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}

fn connection_options() -> BundleConnectionOptions {
    BundleConnectionOptions {
        scene: CurvedBundleSettings {
            trace: trace_settings(),
            position_derivative_tolerance_m: 1e-11,
            direction_derivative_tolerance: 1e-12,
            ..CurvedBundleSettings::default()
        },
        ..BundleConnectionOptions::default()
    }
}

fn ray(z: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, z),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-6 * expected.abs().max(1e-10),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn branch(
    source: &PhotometricSource,
    receiver_object_id: u64,
    initial_direction: Vec3,
) -> PointLightBranch<'_> {
    PointLightBranch {
        source,
        receiver_object_id,
        initial_medium: 0,
        initial_direction,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        policy: BundleBranchPolicy::Sequence(&[]),
        position_tolerance_m: 2e-6,
    }
}

fn settings() -> PointLightCameraSettings {
    PointLightCameraSettings {
        surfaces: trace_settings(),
        connection: connection_options(),
        max_connections: 16,
        max_trace_evaluations: 256,
    }
}

#[test]
fn g1_native_camera_fixture_preserves_calibrated_cosine_beer_and_surface_emission() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 2.0), plane(9, -1.0, 0.0, 2.0)])
                .unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9, 2.0)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let lamp = source(Point3::new(0.6, 0.0, 0.0));
        // The initial ray reaches the same finite diffuse face away from the
        // receiver. The actual native Newton connector must solve the target.
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let ordinary = camera
            .sample_ray_with_samples(cx, ray(0.0), 550.0, 0, trace_settings(), |_| [0.0; 3])
            .unwrap();
        let sample = lit
            .sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| [0.0; 3])
            .unwrap();
        let distance = 1.36_f64.sqrt();
        let expected_fixture =
            0.6 / core::f64::consts::PI * intensity(&lamp) * (-0.2 * (distance + 1.0)).exp()
                / distance.powi(3);
        near(sample.fixture_radiance_w_m2_sr_nm(), expected_fixture);
        near(ordinary.radiance_w_m2_sr_nm(), 2.0 * 0.6 * (-0.6_f64).exp());
        near(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            ordinary.radiance_w_m2_sr_nm() + expected_fixture,
        );
        assert_eq!(sample.radiance_sample().path(), ordinary.path());
        assert_eq!(sample.connections(), 1);
        assert!(sample.trace_evaluations() > 1);
    });
}

#[test]
fn g1_native_camera_fixture_focusing_uses_the_analytic_hessian_bundle_area() {
    with_cx(|cx, _| {
        let distance = 1.2_f64;
        let geometry =
            InstanceScene::try_new(vec![plane(1, distance, 0.0, 0.2), plane(9, -0.4, 0.0, 0.2)])
                .unwrap();
        let material = glass(1.5, 0.07);
        let field = QuadraticIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            IndexHessian::try_new([-1.5, -1.5, 0.0], [0.0; 3]).unwrap(),
            Aabb::new(
                Point3::new(-0.25, -0.25, -0.5),
                Point3::new(0.25, 0.25, 1.5),
            ),
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: ContentHash([71; 32]),
        }];
        let reflectances = [SpectralReflectance::constant(0.7).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, 0.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let sample = lit
            .sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| [0.0; 3])
            .unwrap();
        // On the axis grad(n)=0, but X''=-X and J=sin(D)^2. Treating the
        // segment as homogeneous would incorrectly replace this by D^2.
        let expected = 0.7 * intensity(&lamp) * (-2.0 * 0.07 * distance).exp()
            / (core::f64::consts::PI * distance.sin().powi(2));
        near(sample.fixture_radiance_w_m2_sr_nm(), expected);
        near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
        assert_eq!(
            sample.radiance_sample().path().termination(),
            PathTermination::Opaque
        );
        assert_eq!(sample.connections(), 1);
    });
}

#[test]
fn g1_native_camera_fixture_uses_one_local_endpoint_index_ratio() {
    with_cx(|cx, _| {
        let distance = 1.2_f64;
        let n_source = 1.2_f64;
        let gradient = 0.15;
        let geometry =
            InstanceScene::try_new(vec![plane(1, distance, 0.0, 0.2), plane(9, -0.4, 0.0, 0.2)])
                .unwrap();
        let material = glass(n_source, 0.05);
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, gradient),
            Aabb::new(
                Point3::new(-0.25, -0.25, -0.5),
                Point3::new(0.25, 0.25, 1.5),
            ),
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: ContentHash([72; 32]),
        }];
        let reflectances = [SpectralReflectance::constant(0.4).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, 0.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let sample = lit
            .sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| [0.0; 3])
            .unwrap();
        let n_receiver = n_source + gradient * distance;
        // Snell's invariant n sin(theta) gives dX/dtheta0 = integral n0/n(z) dz.
        let spread = n_source / gradient * (n_receiver / n_source).ln();
        let expected = 0.4 * intensity(&lamp) * (-2.0 * 0.05 * distance).exp()
            / (core::f64::consts::PI * spread * spread)
            * (n_source / n_receiver).powi(2);
        near(sample.fixture_radiance_w_m2_sr_nm(), expected);
        near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
        assert_eq!(sample.connections(), 1);
    });
}

fn bounce_samples(event: u32) -> [f64; 3] {
    match event {
        0 => [0.0; 3],
        1 => [0.5, 0.25, 0.0],
        _ => panic!("the terminal emitter must not draw a continuation"),
    }
}

#[test]
fn g1_every_diffuse_encounter_accumulates_fixture_light_without_changing_the_walk() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 1.0, 0.0, 0.2),
            plane(2, 0.0, 0.0, 0.2),
            plane(9, 2.0, 2.0, 0.5),
        ])
        .unwrap();
        let material = glass(1.25, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [
            SpectralReflectance::constant(0.5).unwrap(),
            SpectralReflectance::constant(0.25).unwrap(),
        ];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), diffuse(2, 1), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(1, 2.0), emitter(2, 4.0), emitter(9, 3.0)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let lamp = source(Point3::new(0.5, 0.0, 0.5));
        let branches = [
            branch(&lamp, 1, unit(Vec3::new(-1.0, 0.0, 1.0))),
            branch(&lamp, 2, unit(Vec3::new(-1.0, 0.0, -1.0))),
        ];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let ordinary = camera
            .sample_ray_with_samples(cx, ray(0.5), 550.0, 0, trace_settings(), bounce_samples)
            .unwrap();
        let sample = lit
            .sample_ray_with_samples(cx, ray(0.5), 550.0, 0, settings(), bounce_samples)
            .unwrap();
        let expected_fixture =
            (0.5 + 0.5 * 0.25) * intensity(&lamp) * 2.0_f64.sqrt() / core::f64::consts::PI;
        near(
            ordinary.radiance_w_m2_sr_nm(),
            2.0 + 0.5 * 4.0 + 0.5 * 0.25 * 3.0,
        );
        near(sample.fixture_radiance_w_m2_sr_nm(), expected_fixture);
        near(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            ordinary.radiance_w_m2_sr_nm() + expected_fixture,
        );
        assert_eq!(sample.radiance_sample().path(), ordinary.path());
        assert_eq!(sample.radiance_sample().path().surface_events(), 3);
        assert_eq!(sample.connections(), 2);
        assert_eq!(sample.trace_evaluations(), 2);
        // Each individual solve fits a one-trace budget. The second solve must
        // still refuse when the COMPLETE camera sample has only one trace left.
        let limited = PointLightCameraSettings {
            max_trace_evaluations: 1,
            ..settings()
        };
        assert!(matches!(
            lit.sample_ray_with_samples(cx, ray(0.5), 550.0, 0, limited, bounce_samples),
            Err(PointLightCameraError::TraceEvaluationLimit)
        ));
    });
}

#[test]
fn g2_black_diffuse_skips_connection_work_but_a_failed_positive_connection_is_not_black() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 1.0)]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let lamp = source(Point3::new(0.0, 0.0, 0.0));
        // This is an admissible launch frame but misses the finite receiver.
        let branches = [branch(&lamp, 1, Vec3::new(1.0, 0.0, 0.0))];
        for rho in [0.0, 0.5] {
            let reflectances = [SpectralReflectance::constant(rho).unwrap()];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![diffuse(1, 0)],
                100_000,
            )
            .unwrap();
            let emitters = [emitter(1, 2.0)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
            let result = lit.sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| {
                panic!("black or failed source connection must not draw a diffuse continuation")
            });
            if rho == 0.0 {
                let sample = result.unwrap();
                assert_eq!(sample.fixture_radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(sample.connections(), 0);
                assert_eq!(sample.trace_evaluations(), 0);
                assert_eq!(sample.radiance_sample().radiance_w_m2_sr_nm(), 2.0);
            } else {
                assert!(matches!(
                    result,
                    Err(PointLightCameraError::Connection { .. })
                ));
            }
        }
    });
}

#[test]
fn g0_original_source_binding_limits_and_cancellation_refuse_before_partial_camera_output() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 1.0)]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.5).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0)],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, 0.0));
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let duplicate = [branch(&lamp, 1, direction), branch(&lamp, 1, direction)];
        assert!(matches!(
            CurvedPointLightCameraScene::try_new(cx, &camera, &duplicate),
            Err(PointLightCameraError::DuplicateBinding {
                first: 0,
                second: 1
            })
        ));
        let foreign_receiver = [branch(&lamp, 99, direction)];
        assert!(matches!(
            CurvedPointLightCameraScene::try_new(cx, &camera, &foreign_receiver),
            Err(PointLightCameraError::InvalidBinding { index: 0 })
        ));
        let branches = [branch(&lamp, 1, direction)];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let invalid = PointLightCameraSettings {
            max_connections: 0,
            ..settings()
        };
        assert!(matches!(
            lit.sample_ray_with_samples(cx, ray(0.0), 550.0, 0, invalid, |_| panic!(
                "invalid settings drew a continuation"
            )),
            Err(PointLightCameraError::InvalidSettings)
        ));
        gate.request();
        assert!(
            lit.sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| panic!(
                "cancelled camera drew a continuation"
            ))
            .is_err()
        );
    });
}

#[test]
fn g2_positive_fixture_irradiance_underflow_refuses_but_exact_directional_darkness_continues() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 1.0), plane(9, 0.9, 0.0, 1.0)]).unwrap();
        // The camera travels only 0.05+0.1 m: its prefix and continuation
        // remain representable. The source leg is about 2 m, with tau>800.
        let material = glass(1.0, 400.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.4).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let ordinary = camera
            .sample_ray_with_samples(cx, ray(0.95), 550.0, 0, trace_settings(), |_| [0.0; 3])
            .unwrap();
        assert!(ordinary.path().beam().power().i > 0.0);
        let direction = unit(Vec3::new(-2.0, 0.0, 0.05));
        for (dark, distribution) in [
            (false, AngularDistribution::isotropic()),
            (true, AngularDistribution::uniform_cone(0.1).unwrap()),
        ] {
            let lamp = PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [2.0, 0.0, 0.95]).unwrap(),
                100.0,
                1.0,
                distribution,
                SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
            )
            .unwrap();
            assert_eq!(
                lamp.spectral_intensity_w_sr_nm(direction, 550.0).unwrap() == 0.0,
                dark
            );
            let branches = [branch(&lamp, 1, direction)];
            let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
            let mut draws = 0;
            let result = lit.sample_ray_with_samples(cx, ray(0.95), 550.0, 0, settings(), |_| {
                draws += 1;
                [0.0; 3]
            });
            if dark {
                let sample = result.unwrap();
                assert_eq!(sample.fixture_radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(*sample.radiance_sample(), ordinary);
                assert_eq!(sample.connections(), 1);
                assert_eq!(draws, 1);
            } else {
                assert!(matches!(result, Err(PointLightCameraError::NumericalRange)));
                assert_eq!(draws, 0);
            }
        }
    });
}
