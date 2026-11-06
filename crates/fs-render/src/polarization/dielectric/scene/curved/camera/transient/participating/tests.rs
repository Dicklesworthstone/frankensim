mod exposure;

use super::super::{
    LIGHT_M_S, SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError,
};
use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedDensityField, CurvedRayleighMedium,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::sync::Mutex;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 53,
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

fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![
        GeometryInstance::try_new(
            9,
            ContentHash([71; 32]),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-10.0, -10.0, 1.0], [10.0, -10.0, 1.0], [0.0, 10.0, 1.0]],
                vec![[0, 1, 2]],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    ])
    .unwrap()
}

fn bindings() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Opaque,
    }]
}

fn sources() -> [SurfaceEmitter; 1] {
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
                    optical_path_tolerance_m: 1e-11,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-8,
                ..GradientEventOptions::default()
            },
        },
        max_volume_events: 8,
    }
}

struct History<F> {
    evaluate: F,
    queries: Mutex<Vec<SurfaceEmissionQuery>>,
}
impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            queries: Mutex::new(Vec::new()),
        }
    }
    fn binding(&self) -> SurfaceEmissionHistoryBinding<'_>
    where
        F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync,
    {
        SurfaceEmissionHistoryBinding {
            object_id: 9,
            source_identity: ContentHash([73; 32]),
            history: self,
        }
    }
}
impl<F> SurfaceEmissionHistory for History<F>
where
    F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync,
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.evaluate)(query)
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol {tolerance:e}"
    );
}

fn two_turns(flight: u32) -> [f64; 3] {
    match flight {
        0 => [1.0 - (-0.25_f64).exp(), 0.5 + f64::EPSILON, 0.0],
        1 => [1.0 - (-0.5_f64).exp(), 0.5 + f64::EPSILON, 0.5],
        2 => [1.0 - 1e-12, 0.5, 0.0],
        _ => panic!("source precedes the next collision"),
    }
}

#[test]
fn g1_two_molecular_turns_keep_polarization_and_the_dispersive_retarded_clock() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let unity = History::new(|_| Ok(1.0));
        let static_camera =
            TransientCurvedCameraScene::try_new(cx, &camera, &[unity.binding()]).unwrap();
        let static_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &static_camera, &scattering).unwrap();
        let epoch = Dd::from_f64(1_073_741_824.0);
        let pulse = History::new(|query: SurfaceEmissionQuery| {
            Ok((1.0 - ((query.emission_time_s - epoch).to_f64() / 1e-10).abs()).max(0.0))
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[pulse.binding()]).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let time = epoch + Dd::from_f64(1.5 * (1.4 + 0.06 / 0.45_f64.powi(2)) / LIGHT_M_S);
        for wavelength in [450.0, 650.0] {
            let expected = steady
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    wavelength,
                    0,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            let static_sample = static_scene
                .sample_ray_with_samples_at_time(
                    cx,
                    ray(),
                    wavelength,
                    0,
                    time,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            assert_eq!(static_sample.radiance_sample(), expected.radiance_sample());
            assert_eq!(static_sample.volume_events(), expected.volume_events());
            assert_eq!(
                static_sample.free_flight_draws(),
                expected.free_flight_draws()
            );
            let sampled = scene
                .sample_ray_with_samples_at_time(
                    cx,
                    ray(),
                    wavelength,
                    0,
                    time,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            // Two right-angle turns: 0.25 m forward + 0.5 m sideways + 0.75 m
            // forward. The 1.5 m path exceeds the straight source distance.
            near(sampled.radiance_sample().path().traveled_m(), 1.5, 2e-7);
            assert_eq!(sampled.volume_events(), 2);
            assert_eq!(sampled.free_flight_draws(), 3);
            let emitted = if wavelength == 450.0 {
                8.0 * (-0.3_f64).exp()
            } else {
                0.0
            };
            near(
                sampled.radiance_sample().radiance_w_m2_sr_nm(),
                emitted,
                3e-5,
            );
        }
        let queries = pulse.queries.lock().unwrap();
        assert_eq!(queries.len(), 2);
        for (query, wavelength) in queries.iter().zip([450.0_f64, 650.0]) {
            near(
                query.group_delay_s.to_f64(),
                1.5 * (1.4 + 0.06 / (wavelength * 1e-3).powi(2)) / LIGHT_M_S,
                3e-16,
            );
            near(query.position_m.x, 0.5, 2e-7);
            near(query.position_m.z, 1.0, 2e-7);
        }
    });
}

struct AffineIndex;
impl GradientIndexField for AffineIndex {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + 0.2 * p.x;
        GradientIndexSample::try_new(n, Vec3::new(0.2, 0.0, 0.0), 0.3, Some(n))
    }
}
struct VacuumDensity;
impl CurvedDensityField for VacuumDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> {
        Ok(0.0)
    }
}

#[test]
fn g1_null_collision_preserves_the_original_curved_group_integral_and_source_query() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([74; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Ok(2.0));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history.binding()]).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &sigma,
                &VacuumDensity,
                1.0,
                ContentHash([79; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let sample = scene
            .sample_ray_with_tracking_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                Dd::from_f64(1e-8),
                settings(),
                |_| panic!("opaque"),
                |flight| match flight {
                    0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.5],
                    1 => [0.99, 0.5, 0.0, 0.5],
                    _ => panic!("source precedes the next candidate"),
                },
            )
            .unwrap();
        let arc = 6.0 * (1.0_f64 / 6.0).sinh();
        let n_end = 1.2 * (1.0_f64 / 6.0).cosh();
        // n=1.2 cosh(z/6), ds/dz=cosh(z/6), so integral(n ds)
        // from z=0 to z=1 is 0.6 * (1 + 3 sinh(1/3)).
        let group_m = 0.6 * (1.0 + 3.0 * (1.0_f64 / 3.0).sinh());
        near(sample.radiance_sample().path().traveled_m(), arc, 2e-7);
        near(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            8.0 * (1.2 / n_end).powi(2) * (-0.3 * arc).exp(),
            2e-6,
        );
        assert_eq!(sample.volume_events(), 0);
        assert_eq!(sample.free_flight_draws(), 2);
        let queries = history.queries.lock().unwrap();
        assert_eq!(queries.len(), 1);
        near(
            queries[0].group_delay_s.to_f64(),
            group_m / LIGHT_M_S,
            3e-16,
        );
    });
}

#[test]
fn g1_seeded_zero_scattering_preserves_every_transient_surface_bit_and_draw_identity() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([74; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Ok(0.75));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history.binding()]).unwrap();
        let volumes = [None];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let time = Dd::from_f64(1e-8);
        let expected = transient
            .sample_ray_at_time(cx, ray(), 550.0, 0, time, settings().surfaces, 31, 42)
            .unwrap();
        let seeded = scene
            .sample_ray_at_time(cx, ray(), 550.0, 0, time, settings(), 31, 42)
            .unwrap();
        assert_eq!(seeded.transient_sample(), &expected);
        assert_eq!(seeded.volume_events(), 0);
        assert_eq!(seeded.free_flight_draws(), 0);
        let custom = scene
            .sample_ray_with_tracking_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                time,
                settings(),
                |_| panic!("opaque"),
                |_| panic!("zero scattering consumes no flight variate"),
            )
            .unwrap();
        assert_eq!(seeded, custom);
    });
}

#[test]
fn g0_scene_identity_history_failure_and_volume_limits_refuse_partial_observations() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([74; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let other =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history.binding()]).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let wrong = CurvedRayleighScene::try_new(cx, &other, &volumes).unwrap();
        assert!(matches!(
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &wrong),
            Err(TransientCameraError::Camera(CameraRadianceError::Scene(
                CurvedSceneError::Scene(SceneError::InvalidBindings)
            )))
        ));
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let time = Dd::from_f64(1e-8);
        assert!(matches!(
            scene.sample_ray_with_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                time,
                settings(),
                |_| panic!("opaque"),
                |_| [0.99, 0.5, 0.0]
            ),
            Err(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::OutsideTimeDomain
            })
        ));
        let mut limited = settings();
        limited.max_volume_events = 1;
        assert!(matches!(
            scene.sample_ray_with_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                time,
                limited,
                |_| panic!("opaque"),
                |_| [0.01, 0.5, 0.0]
            ),
            Err(TransientCameraError::Camera(
                CameraRadianceError::Unresolved(PathTermination::VolumeLimit)
            ))
        ));
        assert_eq!(history.queries.lock().unwrap().len(), 1);
        gate.request();
        assert!(matches!(
            scene.sample_ray_at_time(cx, ray(), 550.0, 0, time, settings(), 31, 42),
            Err(TransientCameraError::Camera(CameraRadianceError::Scene(
                CurvedSceneError::Scene(SceneError::Cancelled)
            )))
        ));
    });
}
