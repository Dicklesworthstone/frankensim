//! G0/G3/G5 battery for aligned cinematic AOV accumulation and EXR artifacts.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::aov::{
    CinematicAovConfig, CinematicAovError, CinematicAovLimits, CinematicAovProfile,
    CinematicAovProvenance, validity,
};
use fs_render::camera::{AnimatedCamera, Aperture, CutSide, PhysicalCamera};
use fs_render::charts::TriMesh;
use fs_render::instances::{RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    Camera, DirectStrategy, Film, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
    render_cinematic, render_cinematic_range_with_aovs, render_cinematic_with_aovs,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate(&gate, f)
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xa0_51,
                kernel_id: 6,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn hash(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.cinematic-aov", label.as_bytes())
}

fn quad() -> TriMesh {
    TriMesh::new(
        vec![
            [-0.75, -0.75, 0.0],
            [0.75, -0.75, 0.0],
            [0.75, 0.75, 0.0],
            [-0.75, 0.75, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn scene(animated: bool) -> Scene {
    let geometry = SharedGeometry::mesh(quad());
    let shape = if animated {
        let start = TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap();
        let end = TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap();
        Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                101,
                hash("geometry"),
                geometry,
                RigidTransformTrajectory::try_new(vec![start, end]).unwrap(),
            )
            .unwrap(),
        )
    } else {
        Shape::Mesh(quad())
    };
    let reflectance = lift_rgb([0.7, 0.4, 0.2]);
    let emission = (reflectance, 2.0);
    Scene {
        primitives: vec![
            Primitive {
                shape,
                material: Material::Lambertian { reflectance },
                emission: Some(emission),
            },
            Primitive {
                shape: Shape::Mesh(TriMesh::new(
                    vec![
                        [9.0, -0.5, 0.0],
                        [10.0, -0.5, 0.0],
                        [10.0, 0.5, 0.0],
                        [9.0, 0.5, 0.0],
                    ],
                    vec![[0, 1, 2], [0, 2, 3]],
                )),
                material: Material::Lambertian { reflectance },
                emission: Some(emission),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(9.0, -0.5, 0.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 1,
            emission,
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.25,
        },
    }
}

fn camera(scene: &Scene) -> AnimatedCamera {
    let physical = PhysicalCamera::try_legacy_compatible(
        scene.camera.eye,
        scene.camera.forward,
        scene.camera.up,
        scene.camera.half_tan,
        2.0,
        Aperture::try_circular(0.0).unwrap(),
    )
    .unwrap();
    AnimatedCamera::try_static(701, 0.0, 1.0, physical).unwrap()
}

fn settings(spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0xa0_51_0001,
    }
}

fn shutter() -> ShutterInterval {
    ShutterInterval::resolve(
        0.5,
        0.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 1_024 },
        ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
    )
    .unwrap()
}

fn config(profile: CinematicAovProfile) -> CinematicAovConfig {
    CinematicAovConfig::new(
        profile,
        CinematicAovProvenance::try_new(
            12,
            0.5,
            0.0,
            1.0,
            hash("trajectory"),
            hash("scene"),
            hash("composition"),
        )
        .unwrap(),
        CinematicAovLimits::default(),
    )
}

fn assert_film_bits_eq(left: &Film, right: &Film) {
    assert_eq!((left.width, left.height), (right.width, right.height));
    assert_eq!(left.spp_done, right.spp_done);
    assert_eq!(left.time_mode, right.time_mode);
    for (left, right) in left.xyz.iter().zip(&right.xyz) {
        assert_eq!(left.map(f64::to_bits), right.map(f64::to_bits));
    }
}

fn channel<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a fs_img::Channel {
    decoded
        .channels
        .iter()
        .find(|channel| channel.name == name)
        .unwrap_or_else(|| panic!("missing EXR channel {name}"))
}

#[test]
fn g5_aov_capture_preserves_cinematic_beauty_bits() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(3);
    let (legacy, with_aovs) = with_cx(|cx| {
        (
            render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter()).unwrap(),
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter(),
                config(CinematicAovProfile::FinalDiagnostic),
            )
            .unwrap(),
        )
    });
    assert_film_bits_eq(&legacy, with_aovs.beauty());
    assert_eq!(with_aovs.sample_count(0), Some(3));
    assert_eq!(with_aovs.primary_count(0), Some(3));
}

#[test]
fn g0_final_exr_round_trip_has_exact_channels_palettes_and_units() {
    let scene = scene(true);
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(2),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let first = film.to_exr().unwrap();
    let second = film.to_exr().unwrap();
    assert_eq!(first, second, "EXR encoding must be byte deterministic");
    let decoded = fs_img::read_exr(&first).unwrap();
    assert_eq!(decoded.channels.len(), 30);
    assert!(
        decoded
            .channels
            .windows(2)
            .all(|pair| pair[0].name < pair[1].name)
    );
    assert!(decoded.channels.iter().all(|plane| {
        plane.ty == fs_img::PixelType::Float && plane.data.iter().all(|value| value.is_finite())
    }));
    assert!((channel(&decoded, "depth.Z").data[0] - 2.0).abs() < 1.0e-6);
    assert_eq!(channel(&decoded, "id.object").data, [1.0]);
    assert_eq!(channel(&decoded, "id.material").data, [1.0]);
    assert_eq!(channel(&decoded, "samples").data, [2.0]);
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_eq!(
        mask & (validity::PRIMARY
            | validity::ALBEDO
            | validity::PREVIOUS_MOTION
            | validity::OBJECT_ID
            | validity::MATERIAL_ID
            | validity::CONTRIBUTION_SPLIT),
        validity::PRIMARY
            | validity::ALBEDO
            | validity::PREVIOUS_MOTION
            | validity::OBJECT_ID
            | validity::MATERIAL_ID
            | validity::CONTRIBUTION_SPLIT
    );
    assert!(channel(&decoded, "motion.prev.X").data[0] < 0.0);
    let attribute = |name: &str| {
        decoded
            .attributes
            .iter()
            .find(|attribute| attribute.name == name)
            .unwrap_or_else(|| panic!("missing EXR attribute {name}"))
    };
    assert_eq!(attribute("frankensim.aov.authority").value, b"raw-estimate");
    assert_eq!(
        attribute("frankensim.aov.objectPalette").value,
        b"0=unavailable;1=101"
    );
    let material_palette =
        std::str::from_utf8(&attribute("frankensim.aov.materialPalette").value).unwrap();
    assert!(material_palette.starts_with("0=unavailable;1="));
}

#[test]
fn g3_progressive_ranges_equal_straight_through_and_mismatch_rolls_back() {
    let mut scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(4);
    let straight = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let mut progressive = fs_render::aov::CinematicAovFilm::try_new(
        1,
        1,
        config(CinematicAovProfile::FinalDiagnostic),
    )
    .unwrap();
    with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            0,
            2,
            shutter(),
        )
        .unwrap();
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            2,
            4,
            shutter(),
        )
        .unwrap();
    });
    assert_film_bits_eq(straight.beauty(), progressive.beauty());
    assert_eq!(straight.to_exr().unwrap(), progressive.to_exr().unwrap());

    let before = progressive.to_exr().unwrap();
    scene.primitives[0].material = Material::Lambertian {
        reflectance: lift_rgb([0.1, 0.8, 0.3]),
    };
    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            4,
            5,
            shutter(),
        )
    });
    assert_eq!(error, Err(CinematicAovError::ProgressiveBindingMismatch));
    assert_eq!(progressive.to_exr().unwrap(), before);
}

#[test]
fn g4_pre_cancelled_render_leaves_every_public_plane_uncommitted() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(2);
    let mut film = fs_render::aov::CinematicAovFilm::try_new(
        1,
        1,
        config(CinematicAovProfile::FinalDiagnostic),
    )
    .unwrap();
    let gate = CancelGate::new();
    gate.request();
    let error = with_gate(&gate, |cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut film,
            0,
            2,
            shutter(),
        )
    });
    assert!(matches!(error, Err(CinematicAovError::Tracer(_))));
    assert_eq!(film.beauty().spp_done, 0);
    assert_eq!(film.sample_count(0), Some(0));
    assert_eq!(film.primary_count(0), Some(0));
    assert!(film.palette().is_none());
    assert_eq!(film.to_exr(), Err(CinematicAovError::UnboundFilm));
}

#[test]
fn g0_direct_mesh_object_identity_and_motion_are_unavailable_not_forged() {
    let scene = scene(false);
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(1),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert_eq!(channel(&decoded, "id.object").data, [0.0]);
    assert_eq!(channel(&decoded, "motion.prev.X").data, [0.0]);
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_eq!(mask & validity::OBJECT_ID, 0);
    assert_eq!(mask & validity::PREVIOUS_MOTION, 0);
    assert_ne!(mask & validity::PRIMARY, 0);
    assert_ne!(mask & validity::MATERIAL_ID, 0);
}
