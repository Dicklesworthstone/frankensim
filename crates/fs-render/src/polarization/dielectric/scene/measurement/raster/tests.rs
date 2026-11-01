use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, SharedGeometry};
use crate::photometry::source::{SPECTRAL_BINS, wavelength_nm};
use crate::polarization::dielectric::scene::curved::{CurvedTraceSettings, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn air() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, transform: RigidTransform) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("native-raster-plane", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [-2.0, -2.0, 0.0],
                [2.0, -2.0, 0.0],
                [2.0, 2.0, 0.0],
                [-2.0, 2.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        transform,
    )
    .unwrap()
}

fn binding(id: u64, transform: RigidTransform) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Detector {
            reference_axis: transform.transform_vector(Vec3::new(1.0, 0.0, 0.0)),
        },
    }
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 4,
        max_distance_m: 8.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions::default(),
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}

fn trace(
    scene: &CurvedPolarizedScene<'_>,
    cx: &Cx<'_>,
    transform: RigidTransform,
    x: f64,
    y: f64,
) -> CurvedPathResult {
    let origin = transform.transform_point(Point3::new(x, y, -1.0));
    let frame = PolarizationFrame::try_new(
        transform.transform_vector(Vec3::new(0.0, 0.0, 1.0)),
        transform.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
    )
    .unwrap();
    let beam = scene
        .launch_beam(
            origin,
            frame,
            wavelength_nm(0),
            0,
            Stokes::try_new(2.0, 2.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
    scene.trace(cx, beam, 0, settings(), 0, 0).unwrap()
}

#[test]
fn g0_raster_half_open_edges_and_physical_pixel_area_are_explicit() {
    let layout = DetectorRaster::try_new(RigidTransform::identity(), 4.0, 2.0, 4, 2).unwrap();
    assert_eq!((layout.cols(), layout.rows(), layout.len()), (4, 2, 8));
    assert_eq!(layout.pixel_area_m2(), 1.0);
    let map = |x, y| {
        layout
            .classify(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0), 0.0)
            .unwrap()
    };
    assert_eq!(map(-2.0, 1.0), Some(0));
    assert_eq!(map(-1.0, 1.0), Some(1));
    assert_eq!(map(0.0, 0.0), Some(6));
    assert_eq!(map(1.5, -0.5), Some(7));
    assert_eq!(map(2.0, 0.0), None);
    assert_eq!(map(0.0, -1.0), None);
    assert_eq!(map(-2.01, 0.0), None);
    assert_eq!(map(0.0, 1.01), None);
}

#[test]
fn g1_native_contacts_preserve_total_power_crop_and_detector_stokes_frame() {
    with_cx(|cx| {
        let transform = RigidTransform::identity();
        let geometry = InstanceScene::try_new(vec![plane(9, transform)]).unwrap();
        let air = air();
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![binding(9, transform)],
            100_000,
        )
        .unwrap();
        let layout = DetectorRaster::try_new(transform, 2.0, 2.0, 2, 2).unwrap();
        let raster = scene.detector_raster(9, layout).unwrap();
        let mut sum = 0.0;
        let mut cropped = 0.0;
        for (x, y, pixel) in [
            (-0.4, 0.3, Some(0)),
            (0.4, 0.3, Some(1)),
            (-0.4, -0.3, Some(2)),
            (0.4, -0.3, Some(3)),
            (1.5, 0.3, None),
        ] {
            let path = trace(&scene, cx, transform, x, y);
            let contribution = raster.contribution(0, &path).unwrap().unwrap();
            assert_eq!(contribution.pixel, pixel);
            assert!((contribution.stokes_w_nm.i - 2.0).abs() < 1e-12);
            // Source +Q used Y; detector +Q is X. This must share meter rotation.
            assert!((contribution.stokes_w_nm.q + 2.0).abs() < 1e-12);
            let mut meter = scene.power_meter(9).unwrap();
            meter.record_curved(0, &path).unwrap();
            assert_eq!(
                contribution.stokes_w_nm,
                meter.spectral_mean_w_nm(0).unwrap()
            );
            assert_eq!(
                contribution.contact_resolution_m,
                path.max_contact_resolution_m()
            );
            if pixel.is_some() {
                sum += contribution.stokes_w_nm.i;
            } else {
                cropped += contribution.stokes_w_nm.i;
            }
        }
        assert_eq!(sum, 8.0);
        assert_eq!(cropped, 2.0);
    });
}

#[test]
fn g3_native_translated_rotated_sensor_keeps_pixel_and_binding_identity() {
    with_cx(|cx| {
        let q = 0.5_f64.sqrt();
        let transform = RigidTransform::try_new([q, 0.0, 0.0, q], [12.0, -5.0, 7.0]).unwrap();
        let geometry = InstanceScene::try_new(vec![plane(9, transform)]).unwrap();
        let air = air();
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![binding(9, transform)],
            100_000,
        )
        .unwrap();
        let layout = DetectorRaster::try_new(transform, 2.0, 2.0, 2, 2).unwrap();
        let raster = scene.detector_raster(9, layout).unwrap();
        let path = trace(&scene, cx, transform, 0.4, 0.3);
        assert_eq!(
            raster.contribution(0, &path).unwrap().unwrap().pixel,
            Some(1)
        );
        assert_eq!(
            raster.frame_identity(),
            geometry.instances()[0].frame_identity()
        );
        let moved = RigidTransform::try_new([q, 0.0, 0.0, q], [12.0, -6.0, 7.0]).unwrap();
        let moved_geometry = InstanceScene::try_new(vec![plane(9, moved)]).unwrap();
        let moved_scene = CurvedPolarizedScene::try_new(
            cx,
            &moved_geometry,
            &media,
            vec![binding(9, moved)],
            100_000,
        )
        .unwrap();
        let moved_path = trace(&moved_scene, cx, moved, 0.4, 0.3);
        assert_eq!(
            raster.contribution(0, &moved_path),
            Err(RasterError::Meter(MeterError::DetectorGeometryMismatch))
        );
    });
}

#[test]
fn g4_invalid_layout_wrong_plane_wavelength_and_unresolved_path_refuse() {
    let transform = RigidTransform::identity();
    assert_eq!(
        DetectorRaster::try_new(transform, 0.0, 1.0, 2, 2),
        Err(RasterError::InvalidLayout)
    );
    assert_eq!(
        DetectorRaster::try_new(transform, 1.0, 1.0, u32::MAX, u32::MAX),
        Err(RasterError::PixelBudgetExceeded)
    );
    assert_eq!(
        DetectorRaster::try_new(transform, f64::MAX, f64::MAX, 1, 1),
        Err(RasterError::NumericalRange)
    );
    assert_eq!(
        DetectorRaster::try_new(transform, f64::from_bits(1), 1e308, 1, 1),
        Err(RasterError::NumericalRange)
    );
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(9, transform)]).unwrap();
        let air = air();
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![binding(9, transform)],
            100_000,
        )
        .unwrap();
        let layout = DetectorRaster::try_new(transform, 2.0, 2.0, 2, 2).unwrap();
        let raster = scene.detector_raster(9, layout).unwrap();
        let path = trace(&scene, cx, transform, 0.2, 0.3);
        assert_eq!(
            raster.contribution(SPECTRAL_BINS, &path),
            Err(RasterError::Meter(MeterError::InvalidBin))
        );
        assert_eq!(
            raster.contribution(1, &path),
            Err(RasterError::Meter(MeterError::WavelengthMismatch))
        );
        let displaced = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.01]).unwrap();
        let wrong = scene
            .detector_raster(
                9,
                DetectorRaster::try_new(displaced, 2.0, 2.0, 2, 2).unwrap(),
            )
            .unwrap();
        assert_eq!(
            wrong.contribution(0, &path),
            Err(RasterError::PlaneMismatch)
        );
        let beam = scene
            .launch_beam(
                Point3::new(0.2, 0.3, -1.0),
                PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0))
                    .unwrap(),
                wavelength_nm(0),
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let mut limits = settings();
        limits.max_distance_m = 0.1;
        let limited = scene.trace(cx, beam, 0, limits, 0, 0).unwrap();
        assert_eq!(
            raster.contribution(0, &limited),
            Err(RasterError::UnresolvedPath(PathTermination::DistanceLimit))
        );
        assert_eq!(
            scene.detector_raster(42, layout),
            Err(RasterError::Meter(MeterError::NotADetector))
        );

        // A tangent plane at a curved detector does not turn its native
        // geometry into a planar pixel sensor, even at a normal-incidence hit.
        let sphere = GeometryInstance::try_new(
            9,
            hash_domain("native-raster-sphere", b"unit"),
            SharedGeometry::chart(fs_geom::fixtures::SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
            transform,
        )
        .unwrap();
        let curved_geometry = InstanceScene::try_new(vec![sphere]).unwrap();
        let curved_scene = CurvedPolarizedScene::try_new(
            cx,
            &curved_geometry,
            &media,
            vec![binding(9, transform)],
            100_000,
        )
        .unwrap();
        let tangent = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap();
        let curved_raster = curved_scene
            .detector_raster(9, DetectorRaster::try_new(tangent, 2.0, 2.0, 2, 2).unwrap())
            .unwrap();
        let curved_path = trace(&curved_scene, cx, tangent, 0.0, 0.0);
        assert_eq!(
            curved_raster.contribution(0, &curved_path),
            Err(RasterError::NonPlanarContact)
        );
    });
}
