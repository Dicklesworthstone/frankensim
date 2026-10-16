//! G0/G3/G4/G5 battery for progressive cinematic-AOV checkpoints.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::aov::{
    CINEMATIC_AOV_CHECKPOINT_CONTENT_DOMAIN, CinematicAovCheckpointError,
    CinematicAovCheckpointExpectation, CinematicAovCheckpointWriteError, CinematicAovConfig,
    CinematicAovError, CinematicAovFilm, CinematicAovLimits, CinematicAovProfile,
    CinematicAovProvenance,
};
use fs_render::camera::{AnimatedCamera, Aperture, CutSide, PhysicalCamera};
use fs_render::charts::TriMesh;
use fs_render::instances::{RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    Camera, DirectStrategy, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
    TracerError, render_cinematic_range_with_aovs, render_cinematic_with_aovs,
};
use std::convert::Infallible;

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xa0_51,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new(), f)
}

fn hash(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.cinematic-aov-checkpoint",
        label.as_bytes(),
    )
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

fn scene() -> Scene {
    let trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap(),
        TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap(),
    ])
    .unwrap();
    let secondary_trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [4.0, 0.0, 0.0]).unwrap(),
            [0.0, 0.0, 0.0],
        )
        .unwrap(),
        TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [4.0, 0.0, 0.0]).unwrap(),
            [0.0, 0.0, 0.0],
        )
        .unwrap(),
    ])
    .unwrap();
    let reflectance = lift_rgb([0.7, 0.4, 0.2]);
    let emission = (reflectance, 2.0);
    Scene {
        primitives: vec![
            Primitive {
                shape: Shape::AnimatedInstance(
                    AnimatedGeometryInstance::try_new(
                        101,
                        hash("geometry"),
                        SharedGeometry::mesh(quad()),
                        trajectory,
                    )
                    .unwrap(),
                ),
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
            Primitive {
                shape: Shape::AnimatedInstance(
                    AnimatedGeometryInstance::try_new(
                        202,
                        hash("secondary-geometry"),
                        SharedGeometry::mesh(quad()),
                        secondary_trajectory,
                    )
                    .unwrap(),
                ),
                material: Material::Lambertian {
                    reflectance: lift_rgb([0.1, 0.8, 0.3]),
                },
                emission: None,
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

fn settings() -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp: 4,
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
    config_with_limits(profile, CinematicAovLimits::default())
}

fn config_with_limits(
    profile: CinematicAovProfile,
    limits: CinematicAovLimits,
) -> CinematicAovConfig {
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
        limits,
    )
}

fn encode(
    film: &CinematicAovFilm,
    cx: &Cx<'_>,
) -> (Vec<u8>, fs_render::aov::CinematicAovCheckpointReceipt) {
    let mut bytes = Vec::new();
    let receipt = film
        .write_checkpoint::<Infallible>(u64::MAX, cx, |chunk| {
            bytes.extend_from_slice(chunk);
            Ok(())
        })
        .unwrap();
    (bytes, receipt)
}

fn unbound_expectation(profile: CinematicAovProfile) -> CinematicAovCheckpointExpectation {
    CinematicAovCheckpointExpectation::unbound(config(profile))
}

fn bound_expectation(
    profile: CinematicAovProfile,
    committed_samples_per_pixel: u32,
) -> CinematicAovCheckpointExpectation {
    CinematicAovCheckpointExpectation::bound(
        config(profile),
        settings(),
        shutter(),
        701,
        CutSide::After,
        committed_samples_per_pixel,
    )
    .unwrap()
}

fn reseal(bytes: &mut [u8]) {
    let payload_len = bytes.len() - 32;
    let seal = hash_domain(
        CINEMATIC_AOV_CHECKPOINT_CONTENT_DOMAIN,
        &bytes[..payload_len],
    );
    bytes[payload_len..].copy_from_slice(seal.as_bytes());
}

fn schema_v1_first_pixel_offset(bytes: &[u8]) -> usize {
    // FSRAOVC1 fixed header ends with the binding tag at byte 275. The bound
    // extension stores its palette counts at 365 and 369, followed by the two
    // packed palettes. Keep this test parser independent from the production
    // decoder so a bad offset cannot make both accept the same malformed byte.
    assert_eq!(bytes[275], 1);
    let object_count = u32::from_le_bytes(bytes[365..369].try_into().unwrap()) as usize;
    let material_count = u32::from_le_bytes(bytes[369..373].try_into().unwrap()) as usize;
    373 + object_count * 8 + material_count * 32
}

#[test]
fn g5_checkpoint_restore_resume_is_bitwise_equal_to_straight_through() {
    let scene = scene();
    let camera = camera(&scene);
    let settings = settings();
    let reference = with_cx(|cx| {
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

    let mut partial = CinematicAovFilm::try_new(
        settings.width,
        settings.height,
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
            &mut partial,
            0,
            2,
            shutter(),
        )
        .unwrap();
    });
    let palette = partial.palette().unwrap();
    assert_eq!(palette.object_ids(), [101, 202]);
    assert_eq!(palette.material_identities().len(), 2);

    let (bytes, receipt) = with_cx(|cx| encode(&partial, cx));
    assert_eq!(receipt.byte_len(), bytes.len() as u64);
    assert_eq!(receipt.byte_len(), partial.checkpoint_byte_len().unwrap());
    assert_eq!(receipt.samples_per_pixel(), 2);
    let (mut restored, restored_receipt) = with_cx(|cx| {
        CinematicAovFilm::restore_checkpoint(
            bound_expectation(CinematicAovProfile::FinalDiagnostic, 2),
            &bytes,
            bytes.len() as u64,
            cx,
        )
        .unwrap()
    });
    assert_eq!(restored_receipt, receipt);
    assert_eq!(restored, partial);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            bound_expectation(CinematicAovProfile::FinalDiagnostic, 1),
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "committed sample prefix"
        })
    );
    let wrong_cut_side = CinematicAovCheckpointExpectation::bound(
        config(CinematicAovProfile::FinalDiagnostic),
        settings,
        shutter(),
        701,
        CutSide::Before,
        2,
    )
    .unwrap();
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            wrong_cut_side,
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "cut-side convention"
        })
    );
    let (reencoded, second_receipt) = with_cx(|cx| encode(&restored, cx));
    assert_eq!(reencoded, bytes);
    assert_eq!(second_receipt, receipt);

    let cut_side_error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::Before,
            cx,
            &settings,
            &mut restored,
            2,
            3,
            shutter(),
        )
    });
    assert!(matches!(
        cut_side_error,
        Err(fs_render::aov::CinematicAovError::ProgressiveBindingMismatch)
    ));
    assert_eq!(
        with_cx(|cx| encode(&restored, cx)).0,
        bytes,
        "cut-side refusal must leave the restored prefix unchanged"
    );

    with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut restored,
            2,
            4,
            shutter(),
        )
        .unwrap();
    });
    assert_eq!(restored, reference);
    assert_eq!(restored.to_exr().unwrap(), reference.to_exr().unwrap());
}

#[test]
fn g0_all_profiles_round_trip_empty_and_bound_prefixes() {
    for profile in [
        CinematicAovProfile::BeautyOnly,
        CinematicAovProfile::DailyCore,
        CinematicAovProfile::FinalDiagnostic,
    ] {
        let empty = CinematicAovFilm::try_new(1, 1, config(profile)).unwrap();
        let (bytes, receipt) = with_cx(|cx| encode(&empty, cx));
        let (restored, restored_receipt) = with_cx(|cx| {
            CinematicAovFilm::restore_checkpoint(
                unbound_expectation(profile),
                &bytes,
                bytes.len() as u64,
                cx,
            )
            .unwrap()
        });
        assert_eq!(restored, empty, "{profile:?} empty checkpoint drifted");
        assert_eq!(restored_receipt, receipt);

        let scene = scene();
        let camera = camera(&scene);
        let mut bound = CinematicAovFilm::try_new(1, 1, config(profile)).unwrap();
        with_cx(|cx| {
            render_cinematic_range_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings(),
                &mut bound,
                0,
                1,
                shutter(),
            )
            .unwrap();
        });
        let (bytes, receipt) = with_cx(|cx| encode(&bound, cx));
        let (restored, restored_receipt) = with_cx(|cx| {
            CinematicAovFilm::restore_checkpoint(
                bound_expectation(profile, 1),
                &bytes,
                bytes.len() as u64,
                cx,
            )
            .unwrap()
        });
        assert_eq!(restored, bound, "{profile:?} bound checkpoint drifted");
        assert_eq!(restored_receipt, receipt);
    }

    let custom_limits =
        CinematicAovLimits::try_new(2, 1_000_000, 2_000_000, 3_000_000, 4_000, 5_000_000, 99)
            .unwrap();
    let custom_config = config_with_limits(CinematicAovProfile::FinalDiagnostic, custom_limits);
    let film = CinematicAovFilm::try_new(1, 1, custom_config).unwrap();
    let (bytes, _) = with_cx(|cx| encode(&film, cx));
    let (restored, _) = with_cx(|cx| {
        CinematicAovFilm::restore_checkpoint(
            CinematicAovCheckpointExpectation::unbound(custom_config),
            &bytes,
            bytes.len() as u64,
            cx,
        )
        .unwrap()
    });
    assert_eq!(restored.config().limits(), custom_limits);
}

#[test]
fn g0_uniform_renderer_refuses_prefixes_the_checkpoint_contract_cannot_represent() {
    let scene = scene();
    let camera = camera(&scene);
    let mut film =
        CinematicAovFilm::try_new(1, 1, config(CinematicAovProfile::FinalDiagnostic)).unwrap();
    let before = with_cx(|cx| encode(&film, cx)).0;

    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(),
            &mut film,
            0,
            5,
            shutter(),
        )
    });
    assert_eq!(
        error,
        Err(CinematicAovError::Tracer(TracerError::InvalidInput))
    );
    assert_eq!(with_cx(|cx| encode(&film, cx)).0, before);

    let mut inexact_settings = settings();
    inexact_settings.spp = (1 << 24) + 1;
    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &inexact_settings,
            &mut film,
            0,
            inexact_settings.spp,
            shutter(),
        )
    });
    assert_eq!(
        error,
        Err(CinematicAovError::InexactSampleCount {
            samples: (1 << 24) + 1
        })
    );
    assert_eq!(with_cx(|cx| encode(&film, cx)).0, before);

    let mut zero_depth_settings = settings();
    zero_depth_settings.max_depth = 0;
    let mut zero_depth_film =
        CinematicAovFilm::try_new(1, 1, config(CinematicAovProfile::FinalDiagnostic)).unwrap();
    with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &zero_depth_settings,
            &mut zero_depth_film,
            0,
            1,
            shutter(),
        )
        .unwrap();
    });
    let (bytes, _) = with_cx(|cx| encode(&zero_depth_film, cx));
    let expected = CinematicAovCheckpointExpectation::bound(
        config(CinematicAovProfile::FinalDiagnostic),
        zero_depth_settings,
        shutter(),
        701,
        CutSide::After,
        1,
    )
    .unwrap();
    let (restored, _) = with_cx(|cx| {
        CinematicAovFilm::restore_checkpoint(expected, &bytes, bytes.len() as u64, cx).unwrap()
    });
    assert_eq!(restored, zero_depth_film);

    let defaults = CinematicAovLimits::default();
    let palette_limits = CinematicAovLimits::try_new(
        1,
        defaults.max_retained_bytes(),
        defaults.max_export_plane_bytes(),
        defaults.max_export_metadata_bytes(),
        defaults.max_exr_encoder_scratch_bytes(),
        defaults.max_encoded_exr_bytes(),
        1,
    )
    .unwrap();
    let mut palette_limited = CinematicAovFilm::try_new(
        1,
        1,
        config_with_limits(CinematicAovProfile::FinalDiagnostic, palette_limits),
    )
    .unwrap();
    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(),
            &mut palette_limited,
            0,
            1,
            shutter(),
        )
    });
    assert_eq!(
        error,
        Err(CinematicAovError::PaletteLimit {
            kind: "object",
            requested: 2,
            limit: 1
        })
    );
    assert_eq!(palette_limited.beauty().spp_done, 0);
    assert!(palette_limited.palette().is_none());
}

#[test]
fn g4_checkpoint_budget_cancellation_truncation_and_corruption_fail_closed() {
    let scene = scene();
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let required = film.checkpoint_byte_len().unwrap();
    let mut called = false;
    let error = with_cx(|cx| {
        film.write_checkpoint::<Infallible>(required - 1, cx, |_| {
            called = true;
            Ok(())
        })
    });
    assert!(matches!(
        error,
        Err(CinematicAovCheckpointWriteError::Checkpoint(
            CinematicAovCheckpointError::ByteLimit { .. }
        ))
    ));
    assert!(!called, "budget refusal emitted a partial artifact");

    let gate = CancelGate::new();
    gate.request();
    let mut called = false;
    let error = with_gate(&gate, |cx| {
        film.write_checkpoint::<Infallible>(required, cx, |_| {
            called = true;
            Ok(())
        })
    });
    assert!(matches!(
        error,
        Err(CinematicAovCheckpointWriteError::Checkpoint(
            CinematicAovCheckpointError::Cancelled
        ))
    ));
    assert!(!called, "pre-cancellation emitted a partial artifact");

    let gate = CancelGate::new();
    let mut chunks = 0_u32;
    let error = with_gate(&gate, |cx| {
        film.write_checkpoint::<Infallible>(required, cx, |_| {
            chunks += 1;
            if chunks == 2 {
                gate.request();
            }
            Ok(())
        })
    });
    assert!(matches!(
        error,
        Err(CinematicAovCheckpointWriteError::Checkpoint(
            CinematicAovCheckpointError::Cancelled
        ))
    ));
    assert_eq!(chunks, 2, "test must cancel during final seal emission");

    let sink_error = with_cx(|cx| {
        film.write_checkpoint(required, cx, |_| Err::<(), _>("simulated sink refusal"))
    });
    assert!(matches!(
        sink_error,
        Err(CinematicAovCheckpointWriteError::Sink(
            "simulated sink refusal"
        ))
    ));

    let (bytes, _) = with_cx(|cx| encode(&film, cx));
    let expected = bound_expectation(CinematicAovProfile::FinalDiagnostic, 4);
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        with_gate(&gate, |cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::Cancelled)
    );
    for prefix in 0..bytes.len() {
        assert!(
            with_cx(|cx| CinematicAovFilm::restore_checkpoint(
                expected,
                &bytes[..prefix],
                bytes.len() as u64,
                cx
            ))
            .is_err(),
            "truncated prefix {prefix} was accepted"
        );
    }
    let mut corrupted = bytes.clone();
    corrupted[16] ^= 0x80;
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &corrupted,
            corrupted.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::SealMismatch)
    );
    let mut semantically_wrong = bytes.clone();
    semantically_wrong[11] ^= 0x01;
    reseal(&mut semantically_wrong);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &semantically_wrong,
            semantically_wrong.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::SemanticsMismatch)
    );
    let mut wrong_chart_semantics = bytes.clone();
    wrong_chart_semantics[47] ^= 0x01;
    reseal(&mut wrong_chart_semantics);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &wrong_chart_semantics,
            wrong_chart_semantics.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::SemanticsMismatch)
    );

    let pixel = schema_v1_first_pixel_offset(&bytes);
    let common = pixel + 24;
    let mut impossible_albedo_count = bytes.clone();
    impossible_albedo_count[common + 96..common + 100].copy_from_slice(&0_u32.to_le_bytes());
    reseal(&mut impossible_albedo_count);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &impossible_albedo_count,
            impossible_albedo_count.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::InvalidState {
            field: "common pixel"
        })
    );

    let mut impossible_authored_normal = bytes.clone();
    impossible_authored_normal[common + 100..common + 104].copy_from_slice(&1_u32.to_le_bytes());
    reseal(&mut impossible_authored_normal);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &impossible_authored_normal,
            impossible_authored_normal.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::InvalidState {
            field: "common pixel"
        })
    );

    let mut impossible_zero_hit_depth = bytes.clone();
    impossible_zero_hit_depth[common + 48..common + 56]
        .copy_from_slice(&0.0_f64.to_bits().to_le_bytes());
    reseal(&mut impossible_zero_hit_depth);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &impossible_zero_hit_depth,
            impossible_zero_hit_depth.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::InvalidState {
            field: "common pixel"
        })
    );

    let final_diagnostic = common + 108;
    let mut impossible_categorical_rank = bytes.clone();
    impossible_categorical_rank[final_diagnostic + 97..final_diagnostic + 105]
        .copy_from_slice(&1.0_f64.to_bits().to_le_bytes());
    reseal(&mut impossible_categorical_rank);
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &impossible_categorical_rank,
            impossible_categorical_rank.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::InvalidState {
            field: "categorical primary"
        })
    );
    assert!(matches!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            expected,
            &bytes,
            bytes.len() as u64 - 1,
            cx
        )),
        Err(CinematicAovCheckpointError::ByteLimit { .. })
    ));

    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            unbound_expectation(CinematicAovProfile::FinalDiagnostic),
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "bound render state"
        })
    );
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            bound_expectation(CinematicAovProfile::DailyCore, 4),
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "AOV configuration"
        })
    );
    let mut wrong_settings = settings();
    wrong_settings.seed ^= 1;
    let wrong_stream = CinematicAovCheckpointExpectation::bound(
        config(CinematicAovProfile::FinalDiagnostic),
        wrong_settings,
        shutter(),
        701,
        CutSide::After,
        4,
    )
    .unwrap();
    assert_eq!(
        with_cx(|cx| CinematicAovFilm::restore_checkpoint(
            wrong_stream,
            &bytes,
            bytes.len() as u64,
            cx
        )),
        Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "render settings"
        })
    );
}
