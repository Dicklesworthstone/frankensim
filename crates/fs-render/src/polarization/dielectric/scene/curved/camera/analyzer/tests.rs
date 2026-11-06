use super::super::raster::CurvedCameraRaster;
use super::super::{EmissionSides, RadianceEndpoint, RadianceSpectrum, SurfaceEmitter};
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
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, SceneError, SurfaceBinding,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
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
    }
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn add(p: Point3, d: Vec3, distance: f64) -> Point3 {
    Point3::new(
        p.x + distance * d.x,
        p.y + distance * d.y,
        p.z + distance * d.z,
    )
}

fn patch(id: u64, center: Point3, normal: Vec3, radius: f64) -> GeometryInstance {
    let reference = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = cross(reference, normal);
    let u = u.scale(1.0 / u.norm());
    let v = cross(normal, u);
    let vertices = [(-1.0, -1.0), (1.0, -1.0), (0.0, 2.0)]
        .map(|(a, b)| add(add(center, u, a * radius), v, b * radius))
        .map(|p| [p.x, p.y, p.z])
        .to_vec();
    GeometryInstance::try_new(
        id,
        hash_domain("camera-analyzer-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2]])),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn glass(n: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}

fn source(id: u64, value: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        id,
        RadianceSpectrum::try_constant(value).unwrap(),
        EmissionSides::Both,
    )
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn analyzer_admission_bounds_every_physical_stokes_response() {
    let axis = Vec3::new(0.0, 1.0, 0.0);
    for row in [
        Stokes {
            i: 0.5,
            q: 0.3,
            u: 0.4,
            v: 0.0,
        },
        Stokes {
            i: 0.8,
            q: 0.0,
            u: 0.0,
            v: 0.1,
        },
        Stokes::unpolarized(1.0),
    ] {
        let analyzer = CameraAnalyzer::try_new(axis, row).unwrap();
        assert_eq!(analyzer.row(), row);
        assert_eq!(analyzer.reference_axis(), axis);
        let degree = (row.q * row.q + row.u * row.u + row.v * row.v).sqrt();
        assert!(row.i - degree >= 0.0 && row.i + degree <= 1.0);
    }
    for row in [
        Stokes::unpolarized(0.0),
        Stokes::unpolarized(-0.1),
        Stokes::unpolarized(1.01),
        Stokes {
            i: 0.5,
            q: 0.51,
            u: 0.0,
            v: 0.0,
        },
        Stokes {
            i: 0.75,
            q: 0.0,
            u: 0.0,
            v: 0.3,
        },
        Stokes {
            i: 0.5,
            q: 0.0,
            u: f64::NAN,
            v: 0.0,
        },
    ] {
        assert_eq!(
            CameraAnalyzer::try_new(axis, row),
            Err(BeamError::InvalidStokes)
        );
    }
    for axis in [Vec3::new(0.0, 0.0, 0.0), Vec3::new(f64::INFINITY, 0.0, 0.0)] {
        assert_eq!(
            CameraAnalyzer::try_circular(axis, true),
            Err(BeamError::InvalidFrame)
        );
    }
    assert!(CameraAnalyzer::try_linear(axis, f64::NAN).is_err());
    assert!(CameraAnalyzer::try_linear(axis, f64::MAX).is_err());
    for angle in [0.0, 0.13, 0.71, 1.1, 2.8, -1.7] {
        let row = CameraAnalyzer::try_linear(axis, angle).unwrap().row();
        near(row.i, 0.5, 0.0);
        assert!(det::hypot(row.q, row.u) <= 0.5);
        // Malus law independently evaluates a unit fully linear input.
        let input_angle = 0.37_f64;
        let power = row.i + row.q * (2.0 * input_angle).cos() + row.u * (2.0 * input_angle).sin();
        near(power, (input_angle - angle).cos().powi(2), 4e-15);
    }
}

#[test]
fn native_brewster_and_tir_path_resolves_linear_and_circular_analyzers_reciprocally() {
    with_cx(|cx, _| {
        // Forward source -> Brewster reflection -> rotated 45-degree TIR -> detector.
        // n=3/2 to 1: Brewster Rs=25/169,Rp=0. The subsequent TIR has
        // rs=4/5-3i/5, rp=7/25-24i/25; rotating planes 45 degrees produces
        // S=(r/2,0,-2r/5,-3r/10), r=25/169, at the actual detector.
        let d0 = Vec3::new(0.0, 0.0, 1.0);
        let n1 = Vec3::new(2.0, 0.0, -3.0).scale(1.0 / 13.0_f64.sqrt());
        let d1 = Vec3::new(12.0, 0.0, -5.0).scale(1.0 / 13.0);
        let s1 = Vec3::new(0.0, 1.0, 0.0);
        let v1 = cross(d1, s1);
        let s2 = Vec3::new(s1.x + v1.x, s1.y + v1.y, s1.z + v1.z).scale(1.0 / 2.0_f64.sqrt());
        let side = cross(s2, d1);
        let n2 =
            Vec3::new(-d1.x + side.x, -d1.y + side.y, -d1.z + side.z).scale(1.0 / 2.0_f64.sqrt());
        let d2 = Vec3::new(
            d1.x - 2.0 * d1.dot(n2) * n2.x,
            d1.y - 2.0 * d1.dot(n2) * n2.y,
            d1.z - 2.0 * d1.dot(n2) * n2.z,
        );
        let p1 = Point3::new(0.0, 0.0, 0.0);
        let p2 = add(p1, d1, 0.4);
        let detector = add(p2, d2, 0.4);
        let geometry = InstanceScene::try_new(vec![
            patch(1, p1, n1, 0.08),
            patch(2, p2, n2, 0.08),
            patch(9, add(p1, d0, -0.4), d0, 0.08),
        ])
        .unwrap();
        let bulk = glass(1.5);
        let air = glass(1.0);
        let media = [
            SceneMedium::Homogeneous(&bulk),
            SceneMedium::Homogeneous(&air),
        ];
        let mut bindings: Vec<_> = (1..=2)
            .map(|object_id| SurfaceBinding {
                object_id,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 1,
                    surface: DielectricSurface::SMOOTH,
                },
            })
            .collect();
        bindings.push(opaque(9));
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings, 100_000).unwrap();
        let emitters = [source(9, 1.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let ray = Ray {
            origin: detector,
            dir: d2.scale(-1.0),
        };
        let scalar = camera
            .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| [0.5, 0.5, 0.0])
            .unwrap();
        assert_eq!(
            scalar.endpoint(),
            RadianceEndpoint::Surface { object_id: 9 }
        );
        assert_eq!(scalar.path().surface_events(), 3);
        near(scalar.radiance_w_m2_sr_nm(), 1.0, 2e-11);
        // The original proposal divides by the Brewster scalar reflection
        // probability r/2. These conditional values include it exactly once.
        for (analyzer, expected) in [
            (CameraAnalyzer::try_circular(s2, true).unwrap(), 0.2),
            (CameraAnalyzer::try_circular(s2, false).unwrap(), 0.8),
            (
                CameraAnalyzer::try_linear(s2, core::f64::consts::FRAC_PI_4).unwrap(),
                0.1,
            ),
            (
                CameraAnalyzer::try_linear(s2, -core::f64::consts::FRAC_PI_4).unwrap(),
                0.9,
            ),
        ] {
            let sample = camera
                .sample_ray_analyzed_with_samples(cx, ray, 550.0, 0, analyzer, settings(), |_| {
                    [0.5, 0.5, 0.0]
                })
                .unwrap();
            near(sample.radiance_w_m2_sr_nm(), expected, 2e-10);
            near(sample.path().traveled_m(), 1.2, 2e-6);
            let reversed = CameraAnalyzer::try_new(s2.scale(-1.0), analyzer.row()).unwrap();
            let repeated = camera
                .sample_ray_analyzed_with_samples(cx, ray, 550.0, 0, reversed, settings(), |_| {
                    [0.5, 0.5, 0.0]
                })
                .unwrap();
            near(repeated.radiance_w_m2_sr_nm(), expected, 2e-10);
        }
    });
}

struct AffineField;
impl GradientIndexField for AffineField {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        GradientIndexSample::try_new(
            1.2 + 0.2 * p.x,
            Vec3::new(0.2, 0.0, 0.0),
            0.3,
            Some(1.4 + 0.2 * p.x),
        )
    }
}

#[test]
fn native_grin_analyzers_preserve_bending_endpoint_index_beer_and_failure_semantics() {
    with_cx(|cx, gate| {
        let n0 = 1.2_f64;
        let gradient = 0.2_f64;
        let z = 1.0_f64;
        let end_index = n0 * (gradient * z / n0).cosh();
        let distance = n0 / gradient * (gradient * z / n0).sinh();
        let x = (end_index - n0) / gradient;
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(x, 0.0, z),
            Vec3::new(0.0, 0.0, 1.0),
            0.01,
        )])
        .unwrap();
        let field = AffineField;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([92; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = [source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let analyzer = CameraAnalyzer::try_linear(axis, 0.3).unwrap();
        let scalar = camera
            .sample_ray(cx, ray, 550.0, 0, settings(), 31, 1)
            .unwrap();
        let analyzed = camera
            .sample_ray_analyzed(cx, ray, 550.0, 0, analyzer, settings(), 31, 1)
            .unwrap();
        let expected = 2.5 * (n0 / end_index).powi(2) * (-0.3 * distance).exp();
        near(analyzed.radiance_w_m2_sr_nm(), expected, 2e-6);
        assert!(analyzed.path().beam().ray().origin.x > 0.01);
        assert_eq!(analyzed.path().beam().ray(), scalar.path().beam().ray());
        assert_eq!(analyzed.path().traveled_m(), scalar.path().traveled_m());
        assert_eq!(
            analyzed.path().optical_path_length_m(),
            scalar.path().optical_path_length_m()
        );
        assert_eq!(
            analyzed.path().group_optical_path_length_m(),
            scalar.path().group_optical_path_length_m()
        );
        assert_eq!(
            analyzed.path().optical_depth(),
            scalar.path().optical_depth()
        );
        assert_eq!(
            analyzed.path().medium_snapshot_identity(),
            Some(ContentHash([92; 32]))
        );
        let identity = CameraAnalyzer::try_new(ray.dir, Stokes::unpolarized(1.0)).unwrap();
        assert_eq!(
            scalar,
            camera
                .sample_ray_analyzed(cx, ray, 550.0, 0, identity, settings(), 31, 1)
                .unwrap()
        );
        let parallel = CameraAnalyzer::try_linear(ray.dir, 0.0).unwrap();
        assert!(matches!(
            camera.sample_ray_analyzed(cx, ray, 550.0, 0, parallel, settings(), 31, 1),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Beam(BeamError::InvalidFrame)
            )))
        ));
        let mut limited = settings();
        limited.max_distance_m = 0.5;
        assert!(matches!(
            camera.sample_ray_analyzed(cx, ray, 550.0, 0, analyzer, limited, 31, 1),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
        gate.request();
        assert!(matches!(
            camera.sample_ray_analyzed(cx, ray, 550.0, 0, analyzer, settings(), 31, 1),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}

#[test]
fn real_spectral_raster_and_qe_sensor_share_grin_paths_and_analyzed_radiometry() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 0.0, 1.0),
            50.0,
        )])
        .unwrap();
        let field = AffineField;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([93; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = [source(9, 2.5)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.15).unwrap(),
            2.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 3, 2).unwrap();
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let linear = CameraAnalyzer::try_linear(axis, 0.7).unwrap();
        let identity = CameraAnalyzer::try_new(axis, Stokes::unpolarized(1.0)).unwrap();
        let scalar = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 9)
            .unwrap();
        let analyzed = raster
            .sample_pixel_analyzed(cx, &scene, 0, linear, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(
            scalar,
            raster
                .sample_pixel_analyzed(cx, &scene, 0, identity, settings(), 73, 4, 9)
                .unwrap()
        );
        assert_eq!(scalar.camera_rays(), analyzed.camera_rays());
        assert_eq!(scalar.wavelengths_nm(), analyzed.wavelengths_nm());
        assert_eq!(scalar.diffraction_weights(), analyzed.diffraction_weights());
        for (a, b) in analyzed.xyz().into_iter().zip(scalar.xyz()) {
            near(a, 0.5 * b, 2e-12);
        }
        for (a, b) in analyzed
            .spectral_samples()
            .iter()
            .zip(scalar.spectral_samples())
        {
            near(
                a.radiance_w_m2_sr_nm(),
                0.5 * b.radiance_w_m2_sr_nm(),
                2e-12,
            );
            assert_eq!(a.path().beam().ray(), b.path().beam().ray());
            assert_eq!(a.path().optical_depth(), b.path().optical_depth());
        }
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (380.0, [0.0, 1.0, 0.0]),
            (780.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let scalar = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 4, 9)
            .unwrap();
        let analyzed = raster
            .sample_pixel_sensor_analyzed(cx, &scene, &sensor, 0, linear, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(
            scalar,
            raster
                .sample_pixel_sensor_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    identity,
                    settings(),
                    73,
                    4,
                    9
                )
                .unwrap()
        );
        assert_eq!(scalar.camera_rays(), analyzed.camera_rays());
        assert_eq!(
            scalar.detector_wavelengths(),
            analyzed.detector_wavelengths()
        );
        for (a, b) in analyzed
            .qe_weighted_radiance_w_m2_sr()
            .iter()
            .copied()
            .zip(scalar.qe_weighted_radiance_w_m2_sr().iter().copied())
        {
            near(a, 0.5 * b, 2e-11);
        }
        for (a, b) in analyzed
            .photoelectron_rate_per_m2_sr()
            .iter()
            .copied()
            .zip(scalar.photoelectron_rate_per_m2_sr().iter().copied())
        {
            near(a, 0.5 * b, 2e-12 * b.abs().max(1.0));
        }
        assert_eq!(analyzed.photoelectron_rate_per_m2_sr()[2], 0.0);
    });
}
