use super::super::{EmissionSides, RadianceSpectrum, SurfaceEmitter};
use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedDensityField, CurvedRayleighMedium,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

mod thermal;
mod spatial_thermal;
mod microsurface;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
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

fn plane(axis: usize, position: f64) -> GeometryInstance {
    let vertices = if axis == 0 {
        vec![
            [position, -10.0, -10.0],
            [position, 10.0, -10.0],
            [position, 0.0, 10.0],
        ]
    } else {
        vec![
            [-10.0, -10.0, position],
            [10.0, -10.0, position],
            [0.0, 10.0, position],
        ]
    };
    GeometryInstance::try_new(
        9,
        hash_domain("org.frankensim.test.rayleigh-camera", &[axis as u8]),
        SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2]])),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn binding() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Opaque,
    }]
}
fn emitters() -> [SurfaceEmitter; 1] {
    [SurfaceEmitter::new(
        9,
        RadianceSpectrum::try_constant(4.0).unwrap(),
        EmissionSides::Both,
    )]
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn settings() -> CurvedRayleighTraceSettings {
    CurvedRayleighTraceSettings {
        surfaces: CurvedTraceSettings {
            max_surface_events: 8,
            max_distance_m: 4.0,
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
        },
        max_volume_events: 8,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

struct AffineIndex;
impl GradientIndexField for AffineIndex {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + 0.2 * p.x;
        GradientIndexSample::try_new(n, Vec3::new(0.2, 0.0, 0.0), 0.3, Some(n))
    }
}

#[test]
fn g1_zero_scattering_keeps_every_surface_camera_bit_and_consumes_no_flight() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([71; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let zero = MassScatteringSpectrum::constant(0.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let inactive = CurvedRayleighMedium::try_new(&zero, 2.0, &phase).unwrap();
        for law in [None, Some(inactive)] {
            let volumes = [law];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let participating =
                CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            for wavelength in [400.0, 550.0, 700.0] {
                let ordinary = camera
                    .sample_ray(cx, ray(), wavelength, 0, settings().surfaces, 31, 42)
                    .unwrap();
                let molecular = participating
                    .sample_ray(cx, ray(), wavelength, 0, settings(), 31, 42)
                    .unwrap();
                assert_eq!(ordinary, *molecular.radiance_sample());
                assert_eq!(molecular.volume_events(), 0);
                assert_eq!(molecular.free_flight_draws(), 0);
                let custom = participating
                    .sample_ray_with_tracking_samples(
                        cx,
                        ray(),
                        wavelength,
                        0,
                        settings(),
                        |_| panic!("opaque source consumes no surface variate"),
                        |_| panic!("zero scattering consumes no flight variate"),
                    )
                    .unwrap();
                assert_eq!(molecular, custom);
            }
        }
    });
}

#[test]
fn g1_uncollided_grin_radiance_has_scattering_survival_only_in_the_flight_probability() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([73; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let arc = 6.0 * (1.0_f64 / 6.0).sinh();
        let n_endpoint = 1.2 * (1.0_f64 / 6.0).cosh();
        let sample = participating
            .sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("no dielectric"),
                |_| [1.0 - (-0.7 * 2.0 * arc).exp(), 0.5, 0.0],
            )
            .unwrap();
        let conditional = 4.0 * (1.2 / n_endpoint).powi(2) * (-0.3 * arc).exp();
        close(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            conditional,
            2e-6,
        );
        // The analytic probability of this uncollided contribution supplies
        // exp(-sigma_s*s). Including it in packet weight again would square it.
        close(
            sample.radiance_sample().radiance_w_m2_sr_nm() * (-0.7 * arc).exp(),
            4.0 * (1.2 / n_endpoint).powi(2) * (-1.0 * arc).exp(),
            2e-6,
        );
        assert_eq!(sample.volume_events(), 0);
        assert_eq!(sample.free_flight_draws(), 1);
    });
}

#[test]
fn g1_molecular_scattering_reveals_a_finite_source_missed_by_the_camera_ray() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(0, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        assert_eq!(
            camera
                .sample_ray(cx, ray(), 550.0, 0, settings().surfaces, 31, 42)
                .unwrap()
                .radiance_w_m2_sr_nm(),
            0.0
        );
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let participating =
                CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let isotropic = 3.0 * delta / (2.0 + delta);
            let right_angle = isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON;
            let sample = participating
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("no dielectric"),
                    |flight| match flight {
                        0 => [1.0 - (-0.35_f64).exp(), right_angle, 0.0],
                        1 => [1.0 - 1e-12, 0.5, 0.0],
                        _ => panic!("source precedes next collision"),
                    },
                )
                .unwrap();
            close(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                4.0 * (-0.3_f64).exp(),
                2e-6,
            );
            close(sample.radiance_sample().path().traveled_m(), 1.5, 2e-7);
            let power = sample.radiance_sample().path().beam().power();
            close(power.q / power.i, (1.0 - delta) / (1.0 + delta), 2e-12);
            assert_eq!(sample.volume_events(), 1);
            assert_eq!(sample.free_flight_draws(), 2);
        }
    });
}

#[test]
fn g1_second_collision_retains_polarization_in_the_scalar_camera_response() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let participating =
                CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let isotropic = 3.0 * delta / (2.0 + delta);
            // Remapping the mixture CDF can put nominal cosine zero on either
            // side of the angular frame's equator. Approach 90° from positive
            // cosine so azimuth pi consistently returns the second ray to +z.
            let right_angle = isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON;
            let sample = participating
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("no dielectric"),
                    |flight| match flight {
                        0 => [1.0 - (-0.25_f64).exp(), right_angle, 0.0],
                        1 => [1.0 - (-0.5_f64).exp(), right_angle, 0.5],
                        2 => [1.0 - 1e-12, 0.5, 0.0],
                        _ => panic!("source precedes next collision"),
                    },
                )
                .unwrap();
            let polarization = (1.0 - delta) / (1.0 + delta);
            close(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                4.0 * (1.0 + polarization * polarization) * (-0.3_f64).exp(),
                4e-6,
            );
            close(sample.radiance_sample().path().traveled_m(), 1.5, 2e-7);
            assert_eq!(sample.volume_events(), 2);
            assert_eq!(sample.free_flight_draws(), 3);
        }
    });
}

struct VacuumDensity;
impl CurvedDensityField for VacuumDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> {
        Ok(0.0)
    }
}

#[test]
fn g1_null_tracking_preserves_radiance_but_does_not_invent_environment_escape() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.2, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &VacuumDensity,
                1.0,
                ContentHash([79; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let sample = participating
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("no dielectric"),
                |flight| match flight {
                    0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.5],
                    1 => [0.99, 0.5, 0.0, 0.5],
                    _ => panic!("source precedes next candidate"),
                },
            )
            .unwrap();
        close(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            4.0 * (-0.2_f64).exp(),
            2e-6,
        );
        assert_eq!(sample.volume_events(), 0);
        assert_eq!(sample.free_flight_draws(), 2);
        let outward = Ray {
            origin: Point3::new(0.0, 0.0, 2.0),
            ..ray()
        };
        let mut limited = settings();
        limited.surfaces.max_distance_m = 0.5;
        assert!(matches!(
            participating.sample_ray_with_tracking_samples(
                cx,
                outward,
                550.0,
                0,
                limited,
                |_| panic!("no dielectric"),
                |_| [0.99, 0.5, 0.0, 0.5],
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
    });
}

#[test]
fn g0_camera_and_molecular_paths_require_the_identical_original_surface_scene() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let first =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let second =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &first, &sources, None).unwrap();
        let volumes = [None];
        let different = CurvedRayleighScene::try_new(cx, &second, &volumes).unwrap();
        assert!(matches!(
            CurvedRayleighCameraScene::try_new(cx, &camera, &different),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::InvalidBindings
            )))
        ));
    });
}

#[test]
fn g0_volume_distance_and_cancellation_limits_never_become_black_samples() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(2, 1.0)]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, binding(), 100_000).unwrap();
        let sources = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let mut limited = settings();
        limited.max_volume_events = 1;
        assert!(matches!(
            participating.sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                limited,
                |_| panic!("no dielectric"),
                |_| [1.0 - (-0.2_f64).exp(), 0.5, 0.0],
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::VolumeLimit
            ))
        ));
        limited = settings();
        limited.surfaces.max_distance_m = 0.1;
        assert!(matches!(
            participating.sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                limited,
                |_| panic!("no dielectric"),
                |_| [0.99, 0.5, 0.0],
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
        gate.request();
        for result in [
            participating.sample_ray(cx, ray(), 550.0, 0, settings(), 31, 42),
            participating.sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("cancelled surface sample"),
                |_| panic!("cancelled flight sample"),
            ),
            participating.sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("cancelled surface sample"),
                |_| panic!("cancelled flight sample"),
            ),
        ] {
            assert!(matches!(
                result,
                Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                    SceneError::Cancelled
                )))
            ));
        }
    });
}
