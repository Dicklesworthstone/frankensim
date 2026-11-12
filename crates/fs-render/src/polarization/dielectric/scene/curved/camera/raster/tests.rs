use super::*;
use crate::camera::{Aperture, CameraProjection, FocusedDiffraction};
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

pub(super) fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

pub(super) fn camera(radius: f64) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_half_tangent(0.15).unwrap(),
        2.0,
        Aperture::try_circular(radius).unwrap(),
    )
    .unwrap()
}

pub(super) fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![
        GeometryInstance::try_new(
            7,
            hash_domain("curved-camera-raster-test", b"source-plane"),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-50.0, -50.0, 2.0], [50.0, -50.0, 2.0], [0.0, 50.0, 2.0]],
                vec![[0, 1, 2]],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    ])
    .unwrap()
}

/// geometry() closed below: object 8 is an open-top box, x and y in [-2, 2] and z in
/// [-1, 2], whose top edges lie in the plane. The camera, the fixture point sources and
/// every receiver point the tests' cameras see lie inside it; bind it Interaction::Opaque.
/// A trace that samples a volume flight is never certified to escape (curved.rs: Escaped
/// requires no volume flight), so in the open scene a diffuse continuation from the plane
/// into a participating medium runs on to the distance limit, which the camera refuses as
/// Unresolved(DistanceLimit). Inside the box it ends on the box. Every straight segment
/// inside the box is at most its diagonal, sqrt(41) m.
#[cfg(feature = "volumes")]
pub(super) fn enclosed_geometry() -> InstanceScene {
    let mut instances = geometry().instances().to_vec();
    instances.push(
        GeometryInstance::try_new(
            8,
            hash_domain("curved-camera-raster-test", b"black-enclosure"),
            SharedGeometry::mesh(TriMesh::new(
                vec![
                    [-2.0, -2.0, -1.0],
                    [2.0, -2.0, -1.0],
                    [2.0, 2.0, -1.0],
                    [-2.0, 2.0, -1.0],
                    [-2.0, -2.0, 2.0],
                    [2.0, -2.0, 2.0],
                    [2.0, 2.0, 2.0],
                    [-2.0, 2.0, 2.0],
                ],
                vec![
                    [0, 2, 1],
                    [0, 3, 2],
                    [0, 1, 5],
                    [0, 5, 4],
                    [1, 2, 6],
                    [1, 6, 5],
                    [2, 3, 7],
                    [2, 7, 6],
                    [3, 0, 4],
                    [3, 4, 7],
                ],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    );
    InstanceScene::try_new(instances).unwrap()
}

pub(super) fn bindings() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 7,
        interaction: Interaction::Opaque,
    }]
}

pub(super) fn settings() -> CurvedTraceSettings {
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

/// A trace-distance budget that never ends a homogeneous trace inside `geometry`.
/// The box is the axis-aligned box of every mesh vertex and the camera `eye`. Each
/// straight segment between events lies in that box, so it is at most the box
/// diagonal. Escape is certified only once the forward half-line misses the box,
/// widened by the contact tolerance (instances/proximity.rs straight_ray_escapes).
/// The last segment reaches that widened box's boundary within diagonal +
/// 2 sqrt(3) tolerance. A trace has at most `trace.max_surface_events + 1`
/// segments, so it needs at most that many diagonals, each lengthened by
/// 4 tolerances. A diffuse bounce may send its continuation in any direction;
/// this bound holds for every one.
pub(super) fn distance_budget_m(
    geometry: &InstanceScene,
    eye: Point3,
    trace: CurvedTraceSettings,
) -> f64 {
    let mut lower = [eye.x, eye.y, eye.z];
    let mut upper = lower;
    for instance in geometry.instances() {
        let SharedGeometry::Mesh(mesh) = instance.geometry() else {
            panic!(
                "the distance budget reads mesh vertices; object {} is a chart",
                instance.object_id()
            );
        };
        for vertex in &mesh.vertices {
            let world = instance
                .transform()
                .transform_point(Point3::new(vertex[0], vertex[1], vertex[2]));
            for (axis, value) in [world.x, world.y, world.z].into_iter().enumerate() {
                lower[axis] = lower[axis].min(value);
                upper[axis] = upper[axis].max(value);
            }
        }
    }
    let diagonal = (0..3)
        .map(|axis| (upper[axis] - lower[axis]).powi(2))
        .sum::<f64>()
        .sqrt();
    f64::from(trace.max_surface_events + 1) * (diagonal + 4.0 * trace.event.surface_tolerance_m)
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn geometric_aperture_and_spectral_color_preserve_absolute_sample_identity() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let vacuum = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(2.5).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&camera, 3, 2).unwrap();
        let sample = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 9)
            .unwrap();
        let other = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 2)
            .unwrap();
        let repeated = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(sample, repeated);
        assert_ne!(sample.camera_ray(), other.camera_ray());
        let origin = sample.camera_ray().origin;
        assert!(origin.x * origin.x + origin.y * origin.y <= 0.01_f64.powi(2));
        assert_eq!(origin.z, 0.0);
        let mut expected = [0.0; 3];
        for (wavelength, lane) in sample
            .wavelengths_nm()
            .iter()
            .zip(sample.spectral_samples())
        {
            assert_eq!(lane.radiance_w_m2_sr_nm(), 2.5);
            assert_eq!(lane.path().termination(), PathTermination::Opaque);
            for (sum, observer) in expected.iter_mut().zip([
                cie_x(*wavelength),
                cie_y(*wavelength),
                cie_z(*wavelength),
            ]) {
                *sum += 2.5 * observer * 100.0 / y_integral();
            }
        }
        for (value, expected) in sample.xyz().into_iter().zip(expected) {
            near(value, expected, 2e-14);
        }
    });
}

struct DispersiveVerticalIndex;

impl GradientIndexField for DispersiveVerticalIndex {
    fn sample(
        &self,
        p: Point3,
        wavelength: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        let c = (550.0 / wavelength).powi(2);
        let n = 1.28 + 0.03 * p.z + (0.02 + 0.01 * p.z) * c;
        GradientIndexSample::try_new(
            n,
            Vec3::new(0.0, 0.0, 0.03 + 0.01 * c),
            0.12,
            Some(n + 2.0 * (0.02 + 0.01 * p.z) * c),
        )
    }
}

#[test]
fn every_wavelength_follows_its_own_curved_ray_and_beer_length() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = DispersiveVerticalIndex;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([17; 32]),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let sample = raster
            .sample_pixel(cx, &scene, 0, settings(), 79, 0, 6)
            .unwrap();
        let initial = sample.camera_ray();
        let transverse_squared = initial.dir.x * initial.dir.x + initial.dir.y * initial.dir.y;
        assert!(transverse_squared > 0.0);
        for (&wavelength, lane) in sample
            .wavelengths_nm()
            .iter()
            .zip(sample.spectral_samples())
        {
            let c = (550.0 / wavelength).powi(2);
            let n0 = 1.28 + 0.02 * c;
            let n1 = 1.34 + 0.04 * c;
            // Independent planar GRIN invariant p = n sin(theta). Integrating
            // ds/dz = n/sqrt(n²-p²) gives the physical arc length below; the
            // rationalized expression avoids subtracting neighboring roots.
            let p_squared = n0 * n0 * transverse_squared;
            let length =
                2.0 * (n1 + n0) / ((n1 * n1 - p_squared).sqrt() + (n0 * n0 - p_squared).sqrt());
            let expected = 3.0 * (n0 / n1).powi(2) * (-0.12 * length).exp();
            near(lane.radiance_w_m2_sr_nm(), expected, 2e-7);
            near(lane.path().traveled_m(), length, 2e-7);
            near(
                lane.path().beam().ray().dir.x,
                initial.dir.x * n0 / n1,
                2e-8,
            );
            near(
                lane.path().beam().ray().dir.y,
                initial.dir.y * n0 / n1,
                2e-8,
            );
        }
        assert_ne!(
            sample.spectral_samples()[0].path().beam().ray().dir,
            sample.spectral_samples()[1].path().beam().ray().dir,
        );
    });
}

#[test]
fn invalid_raster_diffraction_and_unresolved_paths_do_not_become_black_pixels() {
    let camera = camera(0.01);
    assert!(matches!(
        CurvedCameraRaster::try_new(&camera, 0, 1),
        Err(CurvedRasterError::InvalidRaster)
    ));
    assert!(matches!(
        CurvedCameraRaster::try_new(&camera, u32::MAX, 2),
        Err(CurvedRasterError::InvalidRaster)
    ));
    let diffraction = camera
        .clone()
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.1).unwrap())
        .unwrap();
    assert!(matches!(
        CurvedCameraRaster::try_new(&diffraction, 1, 1),
        Err(CurvedRasterError::Camera(CameraError::InvalidDiffraction))
    ));
    with_cx(|cx, gate| {
        let geometry = geometry();
        let field = DispersiveVerticalIndex;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([17; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        assert!(matches!(
            raster.sample_pixel(cx, &scene, 0, settings(), 3, 1, 0),
            Err(CurvedRasterError::InvalidRaster)
        ));
        let mut limited = settings();
        limited.max_distance_m = 0.5;
        assert!(matches!(
            raster.sample_pixel(cx, &scene, 0, limited, 3, 0, 0),
            Err(CurvedRasterError::Transport(
                CameraRadianceError::Unresolved(PathTermination::DistanceLimit)
            ))
        ));
        gate.request();
        assert!(matches!(
            raster.sample_pixel(cx, &scene, 0, settings(), 3, 0, 0),
            Err(CurvedRasterError::Camera(CameraError::Cancelled))
        ));
    });
}
