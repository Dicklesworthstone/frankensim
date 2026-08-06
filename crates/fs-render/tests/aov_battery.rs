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
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::instances::{RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    AdaptiveDecision, AdaptiveSamplingConfig, Camera, DirectStrategy, Film, Material, Primitive,
    RectLight, Sampler, Scene, Settings, Shape, render_cinematic,
    render_cinematic_adaptive_with_aovs, render_cinematic_range_with_aovs,
    render_cinematic_with_aovs, trace_cinematic_pixel_sample,
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

fn camera_at(scene: &Scene, eye: Point3) -> AnimatedCamera {
    let physical = PhysicalCamera::try_legacy_compatible(
        eye,
        scene.camera.forward,
        scene.camera.up,
        scene.camera.half_tan,
        2.0,
        Aperture::try_circular(0.0).unwrap(),
    )
    .unwrap();
    AnimatedCamera::try_static(701, 0.0, 1.0, physical).unwrap()
}

fn camera(scene: &Scene) -> AnimatedCamera {
    camera_at(scene, scene.camera.eye)
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

fn config_with_times(
    profile: CinematicAovProfile,
    previous_frame_time_s: f64,
    frame_time_s: f64,
    next_frame_time_s: f64,
) -> CinematicAovConfig {
    CinematicAovConfig::new(
        profile,
        CinematicAovProvenance::try_new(
            12,
            frame_time_s,
            previous_frame_time_s,
            next_frame_time_s,
            hash("trajectory"),
            hash("scene"),
            hash("composition"),
        )
        .unwrap(),
        CinematicAovLimits::default(),
    )
}

fn config(profile: CinematicAovProfile) -> CinematicAovConfig {
    config_with_times(profile, 0.0, 0.5, 1.0)
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

fn attribute<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a fs_img::ExrAttribute {
    decoded
        .attributes
        .iter()
        .find(|attribute| attribute.name == name)
        .unwrap_or_else(|| panic!("missing EXR attribute {name}"))
}

fn closed_slab_mesh(z_front: f64, z_back: f64) -> TriMesh {
    let half = 4.0;
    TriMesh::new(
        vec![
            [-half, -half, z_back],
            [half, -half, z_back],
            [half, half, z_back],
            [-half, half, z_back],
            [-half, -half, z_front],
            [half, -half, z_front],
            [half, half, z_front],
            [-half, half, z_front],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    )
}

fn dielectric_scene() -> Scene {
    let white = lift_rgb([1.0, 1.0, 1.0]);
    let emission = (white, 3.0);
    let glass = Material::Dielectric {
        glass: DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    };
    let emitter = TriMesh::new(
        vec![
            [-2.0, -2.0, -1.0],
            [2.0, -2.0, -1.0],
            [2.0, 2.0, -1.0],
            [-2.0, 2.0, -1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    );
    Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Mesh(closed_slab_mesh(0.0, -0.2)),
                material: glass,
                emission: None,
            },
            Primitive {
                shape: Shape::Mesh(emitter),
                material: Material::Lambertian { reflectance: white },
                emission: Some(emission),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(-2.0, -2.0, -1.0),
            edge_u: Vec3::new(4.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 4.0, 0.0),
            prim: 1,
            emission,
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.1,
        },
    }
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
    assert_eq!(
        attribute(&decoded, "frankensim.aov.authority").value,
        b"raw-estimate"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.aov.objectPalette").value,
        b"0=unavailable;1=101"
    );
    let material_palette =
        std::str::from_utf8(&attribute(&decoded, "frankensim.aov.materialPalette").value).unwrap();
    assert!(material_palette.starts_with("0=unavailable;1="));
}

#[test]
fn g0_every_profile_has_its_exact_frozen_channel_schema() {
    let scene = scene(true);
    let camera = camera(&scene);
    let cases: &[(CinematicAovProfile, &[&str])] = &[
        (CinematicAovProfile::BeautyOnly, &["B", "G", "R"]),
        (
            CinematicAovProfile::DailyCore,
            &[
                "B",
                "G",
                "R",
                "albedo.B",
                "albedo.G",
                "albedo.R",
                "depth.Z",
                "motion.prev.X",
                "motion.prev.Y",
                "normal.X",
                "normal.Y",
                "normal.Z",
                "primary.coverage",
                "variance.Y",
            ],
        ),
        (
            CinematicAovProfile::FinalDiagnostic,
            &[
                "B",
                "G",
                "R",
                "albedo.B",
                "albedo.G",
                "albedo.R",
                "depth.Z",
                "diagnostic.validity",
                "direct.B",
                "direct.G",
                "direct.R",
                "emission.B",
                "emission.G",
                "emission.R",
                "id.material",
                "id.object",
                "indirect.B",
                "indirect.G",
                "indirect.R",
                "motion.prev.X",
                "motion.prev.Y",
                "normal.X",
                "normal.Y",
                "normal.Z",
                "normal_geom.X",
                "normal_geom.Y",
                "normal_geom.Z",
                "primary.coverage",
                "samples",
                "variance.Y",
            ],
        ),
    ];

    for &(profile, expected_names) in cases {
        let film = with_cx(|cx| {
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings(1),
                shutter(),
                config(profile),
            )
            .unwrap()
        });
        let first = film.to_exr().unwrap();
        assert_eq!(
            first,
            film.to_exr().unwrap(),
            "{profile:?} EXR bytes drifted"
        );
        let decoded = fs_img::read_exr(&first).unwrap();
        let actual_names = decoded
            .channels
            .iter()
            .map(|channel| channel.name.as_str())
            .collect::<Vec<_>>();
        assert_eq!(actual_names, expected_names, "{profile:?} schema drifted");
        assert!(
            decoded
                .channels
                .iter()
                .all(|channel| channel.ty == fs_img::PixelType::Float)
        );
    }
}

#[test]
fn g0_true_background_exports_finite_zero_surface_aovs_and_zero_ids() {
    let mut background = scene(false);
    background.primitives.clear();
    background.lights.clear();
    let camera = camera(&background);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &background,
            &camera,
            CutSide::After,
            cx,
            &settings(3),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    assert_eq!(film.sample_count(0), Some(3));
    assert_eq!(film.primary_count(0), Some(0));

    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    for plane in &decoded.channels {
        assert!(plane.data.iter().all(|value| value.is_finite()));
        match plane.name.as_str() {
            "samples" => assert_eq!(plane.data, [3.0]),
            "diagnostic.validity" => {
                assert_eq!(plane.data, [validity::CONTRIBUTION_SPLIT as f32]);
            }
            _ => assert_eq!(plane.data, [0.0], "background {} was nonzero", plane.name),
        }
    }
    assert_eq!(channel(&decoded, "primary.coverage").data, [0.0]);
    assert_eq!(channel(&decoded, "id.object").data, [0.0]);
    assert_eq!(channel(&decoded, "id.material").data, [0.0]);
    assert_eq!(
        attribute(&decoded, "frankensim.aov.objectPalette").value,
        b"0=unavailable"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.aov.materialPalette").value,
        b"0=unavailable"
    );
}

#[test]
fn g0_camera_visible_emission_split_reconstructs_each_beauty_sample_and_exr() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(4);
    with_cx(|cx| {
        for sample_index in 0..settings.spp {
            let sample = trace_cinematic_pixel_sample(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter(),
                0,
                sample_index,
            )
            .unwrap();
            assert!(sample.primary.is_some());
            assert!(sample.emission_xyz.iter().any(|value| *value > 0.0));
            assert_eq!(sample.direct_xyz, [0.0; 3]);
            assert_eq!(sample.indirect_xyz, [0.0; 3]);
            for component in 0..3 {
                let reconstructed = sample.direct_xyz[component]
                    + sample.indirect_xyz[component]
                    + sample.emission_xyz[component];
                let tolerance = 16.0
                    * f64::EPSILON
                    * sample.xyz[component]
                        .abs()
                        .max(reconstructed.abs())
                        .max(1.0);
                assert!(
                    (sample.xyz[component] - reconstructed).abs() <= tolerance,
                    "sample {sample_index} XYZ component {component} failed reconstruction"
                );
            }
        }
    });

    let film = with_cx(|cx| {
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
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    for suffix in ["R", "G", "B"] {
        let beauty = channel(&decoded, suffix).data[0];
        let reconstructed = channel(&decoded, &format!("direct.{suffix}")).data[0]
            + channel(&decoded, &format!("indirect.{suffix}")).data[0]
            + channel(&decoded, &format!("emission.{suffix}")).data[0];
        let tolerance = 4.0 * f32::EPSILON * beauty.abs().max(reconstructed.abs()).max(1.0);
        assert!((beauty - reconstructed).abs() <= tolerance);
        assert_eq!(channel(&decoded, &format!("direct.{suffix}")).data, [0.0]);
        assert_eq!(channel(&decoded, &format!("indirect.{suffix}")).data, [0.0]);
    }
}

#[test]
fn g3_adaptive_aovs_preserve_the_exact_uniform_prefix_and_export_terminal_counts() {
    let scene = scene(true);
    let camera = camera(&scene);
    let maximum = settings(4);
    let policy = AdaptiveSamplingConfig::try_new(2, 1, 1.0e30, 0.0, 0.0).unwrap();
    let adaptive = with_cx(|cx| {
        render_cinematic_adaptive_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &maximum,
            policy,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    assert_eq!(adaptive.beauty().sample_counts(), [2]);
    assert_eq!(
        adaptive.beauty().decisions(),
        [AdaptiveDecision::ErrorThreshold]
    );

    let prefix_settings = settings(2);
    let uniform_prefix = with_cx(|cx| {
        render_cinematic(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &prefix_settings,
            shutter(),
        )
        .unwrap()
    });
    for (adaptive, uniform) in adaptive.beauty().xyz_sums()[0]
        .iter()
        .zip(uniform_prefix.xyz[0])
    {
        assert_eq!(adaptive.to_bits(), uniform.to_bits());
    }

    let first = adaptive.to_exr().unwrap();
    assert_eq!(first, adaptive.to_exr().unwrap());
    let decoded = fs_img::read_exr(&first).unwrap();
    assert_eq!(channel(&decoded, "samples").data, [2.0]);
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    assert_eq!(
        attribute(&decoded, "frankensim.render.sampleMode").value,
        b"adaptive"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.render.spp").value,
        b"per-pixel-channel"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.render.sppCeiling").value,
        b"4"
    );
    let adaptive_policy =
        std::str::from_utf8(&attribute(&decoded, "frankensim.render.adaptivePolicy").value)
            .unwrap();
    assert!(adaptive_policy.starts_with("version=1;minimum=2;batch=1;"));
}

#[test]
fn g0_dielectric_first_hit_retains_surface_diagnostics_without_forging_albedo() {
    let scene = dielectric_scene();
    let camera = camera(&scene);
    let mut settings = settings(2);
    settings.max_depth = 3;
    let film = with_cx(|cx| {
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
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    assert!((channel(&decoded, "depth.Z").data[0] - 1.0).abs() <= 1.0e-6);
    for name in ["albedo.R", "albedo.G", "albedo.B"] {
        assert_eq!(channel(&decoded, name).data, [0.0]);
    }
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_ne!(mask & validity::PRIMARY, 0);
    assert_ne!(mask & validity::MATERIAL_ID, 0);
    assert_eq!(mask & validity::ALBEDO, 0);
    assert!(channel(&decoded, "id.material").data[0] > 0.0);

    for prefix in ["normal", "normal_geom"] {
        let squared_norm = ["X", "Y", "Z"]
            .into_iter()
            .map(|suffix| {
                let value = channel(&decoded, &format!("{prefix}.{suffix}")).data[0];
                value * value
            })
            .sum::<f32>();
        assert!((squared_norm - 1.0).abs() <= 8.0 * f32::EPSILON);
    }
    assert!(
        ["indirect.R", "indirect.G", "indirect.B"]
            .into_iter()
            .any(|name| channel(&decoded, name).data[0] > 0.0),
        "the transmitted emitter should be classified as a multi-bounce contribution"
    );
    for suffix in ["R", "G", "B"] {
        let beauty = channel(&decoded, suffix).data[0];
        let reconstructed = channel(&decoded, &format!("direct.{suffix}")).data[0]
            + channel(&decoded, &format!("indirect.{suffix}")).data[0]
            + channel(&decoded, &format!("emission.{suffix}")).data[0];
        let tolerance = 4.0 * f32::EPSILON * beauty.abs().max(reconstructed.abs()).max(1.0);
        assert!((beauty - reconstructed).abs() <= tolerance);
    }
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
    let changed_camera = camera_at(&scene, Point3::new(0.125, 0.0, 2.0));
    let camera_error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &changed_camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            4,
            5,
            shutter(),
        )
    });
    assert_eq!(
        camera_error,
        Err(CinematicAovError::ProgressiveBindingMismatch)
    );
    assert_eq!(progressive.to_exr().unwrap(), before);

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
fn g4_reference_time_refusal_precedes_any_aov_publication() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(2);
    let bad_config = config_with_times(CinematicAovProfile::FinalDiagnostic, 0.0, 0.75, 1.0);
    let mut film = fs_render::aov::CinematicAovFilm::try_new(1, 1, bad_config).unwrap();
    let error = with_cx(|cx| {
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
    assert_eq!(
        error,
        Err(CinematicAovError::ReferenceTimesDoNotCoverShutter)
    );
    assert_eq!(film.beauty().spp_done, 0);
    assert_eq!(film.sample_count(0), Some(0));
    assert_eq!(film.primary_count(0), Some(0));
    assert!(film.palette().is_none());
    assert_eq!(film.to_exr(), Err(CinematicAovError::UnboundFilm));
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
