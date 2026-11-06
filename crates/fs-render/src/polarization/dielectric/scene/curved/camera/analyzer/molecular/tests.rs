use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::raster::CurvedCameraRaster;
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedDensityField, CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

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
fn glass(n: f64, alpha: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(alpha).unwrap(),
        GlassProvenance::Custom,
    )
}
fn geometry(id: u64, vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("molecular-analyzer", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices, triangles)),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn plane(id: u64, x: f64) -> GeometryInstance {
    geometry(
        id,
        vec![[x, -10.0, -10.0], [x, 10.0, -10.0], [x, 0.0, 10.0]],
        vec![[0, 1, 2]],
    )
}
fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}
fn source() -> [SurfaceEmitter; 1] {
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
                    optical_depth_tolerance: 1e-12,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-7,
                ..GradientEventOptions::default()
            },
        },
        max_volume_events: 8,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
fn right_angle(delta: f64) -> f64 {
    let isotropic = 3.0 * delta / (2.0 + delta);
    isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON
}
fn collision(flight: u32, delta: f64) -> [f64; 3] {
    if flight == 0 {
        [1.0 - (-0.35_f64).exp(), right_angle(delta), 0.0]
    } else {
        [1.0 - 1e-12, 0.5, 0.0]
    }
}

#[test]
fn native_right_angle_molecular_polarization_obeys_analyzer_angle_and_frame_rotation() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = source();
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let axis = Vec3::new(0.0, 1.0, 0.0);
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &transport).unwrap();
            let g = (1.0 - delta) / (1.0 + delta);
            // A right-angle molecular collision has M01/M00=g; its native
            // proposal already divides by M00. Remaining 1.5 m travel owns Beer.
            for angle in [
                0.0,
                core::f64::consts::FRAC_PI_8,
                core::f64::consts::FRAC_PI_4,
            ] {
                let analyzer = CameraAnalyzer::try_linear(axis, angle).unwrap();
                let sample = scene
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("opaque source has no surface draw"),
                        |flight| collision(flight, delta),
                    )
                    .unwrap();
                let expected = 4.0 * (-0.3_f64).exp() * 0.5 * (1.0 + g * (2.0 * angle).cos());
                close(
                    sample.radiance_sample().radiance_w_m2_sr_nm(),
                    expected,
                    3e-6,
                );
                close(sample.radiance_sample().path().traveled_m(), 1.5, 2e-7);
                assert_eq!(sample.volume_events(), 1);
                assert_eq!(sample.free_flight_draws(), 2);
                let reversed = CameraAnalyzer::try_linear(axis.scale(-1.0), angle).unwrap();
                let other = scene
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        reversed,
                        settings(),
                        |_| panic!("opaque source"),
                        |flight| collision(flight, delta),
                    )
                    .unwrap();
                close(
                    other.radiance_sample().radiance_w_m2_sr_nm(),
                    expected,
                    3e-6,
                );
            }
            // Incoming travel is -z: rotating the reference y toward +x by
            // 45 degrees and then analyzing another 45 degrees selects -Q.
            // A nonzero depolarization leaves a finite, independently known signal.
            if delta > 0.0 {
                let rotated = Vec3::new(1.0, 1.0, 0.0).scale(core::f64::consts::FRAC_1_SQRT_2);
                let analyzer =
                    CameraAnalyzer::try_linear(rotated, core::f64::consts::FRAC_PI_4).unwrap();
                let sample = scene
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("opaque source"),
                        |flight| collision(flight, delta),
                    )
                    .unwrap();
                close(
                    sample.radiance_sample().radiance_w_m2_sr_nm(),
                    4.0 * (-0.3_f64).exp() * 0.5 * (1.0 - g),
                    3e-6,
                );
            }
            let mut limited = settings();
            limited.max_volume_events = 1;
            assert!(matches!(
                scene.sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    CameraAnalyzer::try_circular(axis, true).unwrap(),
                    limited,
                    |_| panic!("opaque source"),
                    |flight| collision(flight, delta)
                ),
                Err(CameraRadianceError::Unresolved(
                    PathTermination::VolumeLimit
                ))
            ));
        }
    });
}

struct AffineIndex;
impl GradientIndexField for AffineIndex {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + 0.2 * p.x;
        GradientIndexSample::try_new(n, Vec3::new(0.2, 0.0, 0.0), 0.3, Some(n + 0.1))
    }
}
struct HalfDensity;
impl CurvedDensityField for HalfDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> {
        Ok(0.5)
    }
}

#[test]
fn curved_rayleigh_refraction_preserves_physical_beer_index_pdf_and_null_tracking() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0), plane(9, 2.0)]).unwrap();
        let field = AffineIndex;
        let exit = glass(1.5, 0.1);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([107; 32]),
            },
            SceneMedium::Homogeneous(&exit),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = source();
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let delta = 0.03;
        let g = (1.0 - delta) / (1.0 + delta);
        let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let law = CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap();
        let volumes = [Some(law), Some(law)];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &transport).unwrap();
        let axis = Vec3::new(0.0, 1.0, 0.0);
        // At arc s=.5 before the collision, n²=n0²+(g_index*s)². A 90-degree
        // turn conserves subsequent n*dz/ds=-.1 up to and beyond the x=1 face.
        let ci = (1.0_f64 - (0.1_f64 / 1.4).powi(2)).sqrt();
        let ct = (1.0_f64 - (0.1_f64 / 1.5).powi(2)).sqrt();
        let rs = ((1.4 * ci - 1.5 * ct) / (1.4 * ci + 1.5 * ct)).powi(2);
        let rp = ((1.5 * ci - 1.4 * ct) / (1.5 * ci + 1.4 * ct)).powi(2);
        let diattenuation = (rp - rs) / (2.0 - rs - rp);
        let grin_arc = 0.5 + ((1.4_f64.powi(2) - 0.1_f64.powi(2)).sqrt() - 1.2) / 0.2;
        let exit_arc = 1.0 / ct;
        let radiance = 4.0 * (1.2_f64 / 1.5).powi(2) * (-0.3 * grin_arc - 0.1 * exit_arc).exp();
        for angle in [
            0.0,
            core::f64::consts::FRAC_PI_4,
            core::f64::consts::FRAC_PI_2,
        ] {
            let analyzer = CameraAnalyzer::try_linear(axis, angle).unwrap();
            let sample = scene
                .sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| [0.5, 0.5, 0.99],
                    |flight| collision(flight, delta),
                )
                .unwrap();
            let q = 0.5 * (2.0 * angle).cos();
            let expected = radiance * ((0.5 + g * q) + diattenuation * g * (0.5 + q));
            close(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                expected,
                4e-6,
            );
            close(
                sample.radiance_sample().path().traveled_m(),
                grin_arc + exit_arc,
                3e-6,
            );
            close(
                sample.radiance_sample().path().beam().ray().dir.z,
                -0.1 / 1.5,
                3e-7,
            );
            assert_eq!(sample.radiance_sample().path().medium_index(), 1);
            assert_eq!(sample.volume_events(), 1);
            assert_eq!(sample.free_flight_draws(), 3);
            // Same actual density, doubled majorant: the first .25 m candidate
            // rejects, the second accepts. Null tracking must preserve Stokes.
            let tracking_spectrum = MassScatteringSpectrum::constant(1.4).unwrap();
            let tracked_law = CurvedRayleighMedium::try_new_heterogeneous(
                &tracking_spectrum,
                &HalfDensity,
                1.0,
                ContentHash([108; 32]),
                &phase,
            )
            .unwrap();
            let tracked_volumes = [Some(tracked_law), Some(tracked_law)];
            let tracked_transport =
                CurvedRayleighScene::try_new(cx, &native, &tracked_volumes).unwrap();
            let tracked =
                CurvedRayleighCameraScene::try_new(cx, &surface, &tracked_transport).unwrap();
            let tracked_sample = tracked
                .sample_ray_analyzed_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| [0.5, 0.5, 0.99],
                    |flight| match flight {
                        0 => [1.0 - (-0.35_f64).exp(), right_angle(delta), 0.0, 0.75],
                        1 => [1.0 - (-0.35_f64).exp(), right_angle(delta), 0.0, 0.25],
                        _ => [1.0 - 1e-12, 0.5, 0.0, 0.25],
                    },
                )
                .unwrap();
            close(
                tracked_sample.radiance_sample().radiance_w_m2_sr_nm(),
                expected,
                4e-6,
            );
            assert_eq!(tracked_sample.volume_events(), 1);
            assert_eq!(tracked_sample.free_flight_draws(), 4);
        }
        let identity = CameraAnalyzer::try_new(axis, Stokes::unpolarized(1.0)).unwrap();
        let original = scene
            .sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| [0.5, 0.5, 0.99],
                |flight| collision(flight, delta),
            )
            .unwrap();
        assert_eq!(
            original,
            scene
                .sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    identity,
                    settings(),
                    |_| [0.5, 0.5, 0.99],
                    |flight| collision(flight, delta)
                )
                .unwrap()
        );
    });
}

#[test]
fn real_xyz_and_qe_molecular_pixels_retain_paths_counters_and_complementary_channels() {
    with_cx(|cx, _| {
        // Closed finite emission catches every scattered direction without an
        // invented homogeneous escape through a positive scattering medium.
        let cube = geometry(
            9,
            vec![
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ],
            vec![
                [0, 2, 1],
                [0, 3, 2],
                [4, 5, 6],
                [4, 6, 7],
                [0, 1, 5],
                [0, 5, 4],
                [3, 7, 6],
                [3, 6, 2],
                [0, 4, 7],
                [0, 7, 3],
                [1, 2, 6],
                [1, 6, 5],
            ],
        );
        let geometry = InstanceScene::try_new(vec![cube]).unwrap();
        let glass = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = source();
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::rayleigh_depolarized(0.03).unwrap();
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &transport).unwrap();
        let camera = PhysicalCamera::try_look_at(
            ray().origin,
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.15).unwrap(),
            2.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 3, 2).unwrap();
        let mut limits = settings();
        limits.surfaces.max_distance_m = 32.0;
        limits.max_volume_events = 64;
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let identity = CameraAnalyzer::try_new(axis, Stokes::unpolarized(1.0)).unwrap();
        let analyzers = [0.4, -0.4].map(|q| {
            CameraAnalyzer::try_new(
                axis,
                Stokes {
                    i: 0.5,
                    q,
                    u: 0.0,
                    v: 0.0,
                },
            )
            .unwrap()
        });
        let original = scene
            .sample_ray(cx, ray(), 550.0, 0, limits, 73, 42)
            .unwrap();
        assert_eq!(
            original,
            scene
                .sample_ray_analyzed(cx, ray(), 550.0, 0, identity, limits, 73, 42)
                .unwrap()
        );
        // This identity has a first native flight below .25 m in lane zero,
        // strictly inside the cube for every admitted aperture point.
        let scalar = raster
            .sample_pixel_rayleigh(cx, &scene, 0, limits, 73, 4, 0)
            .unwrap();
        let clear = raster
            .sample_pixel_rayleigh_analyzed(cx, &scene, 0, identity, limits, 73, 4, 0)
            .unwrap();
        assert_eq!(scalar, clear);
        assert!(scalar.volume_events().iter().any(|&count| count > 0));
        let channels = analyzers.map(|a| {
            raster
                .sample_pixel_rayleigh_analyzed(cx, &scene, 0, a, limits, 73, 4, 0)
                .unwrap()
        });
        for channel in &channels {
            assert_eq!(channel.volume_events(), scalar.volume_events());
            assert_eq!(channel.free_flight_draws(), scalar.free_flight_draws());
            assert_eq!(
                channel.pixel_sample().camera_rays(),
                scalar.pixel_sample().camera_rays()
            );
            assert_eq!(
                channel.pixel_sample().wavelengths_nm(),
                scalar.pixel_sample().wavelengths_nm()
            );
        }
        for axis in 0..3 {
            close(
                channels[0].pixel_sample().xyz()[axis] + channels[1].pixel_sample().xyz()[axis],
                scalar.pixel_sample().xyz()[axis],
                3e-11,
            );
        }
        for lane in 0..4 {
            let a = channels[0].pixel_sample().spectral_samples()[lane];
            let b = channels[1].pixel_sample().spectral_samples()[lane];
            let original = scalar.pixel_sample().spectral_samples()[lane];
            close(
                a.radiance_w_m2_sr_nm() + b.radiance_w_m2_sr_nm(),
                original.radiance_w_m2_sr_nm(),
                3e-11,
            );
            assert_eq!(a.path().beam().ray(), original.path().beam().ray());
            assert_eq!(b.path().beam().ray(), original.path().beam().ray());
        }
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (380.0, [0.0, 1.0, 0.0]),
            (780.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let scalar = raster
            .sample_pixel_sensor_rayleigh(cx, &scene, &sensor, 0, limits, 73, 4, 0)
            .unwrap();
        assert_eq!(
            scalar,
            raster
                .sample_pixel_sensor_rayleigh_analyzed(
                    cx, &scene, &sensor, 0, identity, limits, 73, 4, 0
                )
                .unwrap()
        );
        let channels = analyzers.map(|a| {
            raster
                .sample_pixel_sensor_rayleigh_analyzed(cx, &scene, &sensor, 0, a, limits, 73, 4, 0)
                .unwrap()
        });
        for channel in &channels {
            assert_eq!(channel.volume_events(), scalar.volume_events());
            assert_eq!(channel.free_flight_draws(), scalar.free_flight_draws());
            assert_eq!(
                channel.pixel_sample().detector_wavelengths(),
                scalar.pixel_sample().detector_wavelengths()
            );
        }
        for channel in 0..3 {
            let value = scalar.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel];
            close(
                channels[0].pixel_sample().qe_weighted_radiance_w_m2_sr()[channel]
                    + channels[1].pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                value,
                3e-11 * value.abs().max(1.0),
            );
            let value = scalar.pixel_sample().photoelectron_rate_per_m2_sr()[channel];
            close(
                channels[0].pixel_sample().photoelectron_rate_per_m2_sr()[channel]
                    + channels[1].pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                value,
                3e-11 * value.abs().max(1.0),
            );
        }
        assert_eq!(
            channels[0].pixel_sample().photoelectron_rate_per_m2_sr()[2],
            0.0
        );
    });
}
