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
