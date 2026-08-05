//! G3/E2E coverage for Design-Ledger-backed Euler render checkpoints.

use std::num::NonZeroU32;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::ChannelOwnership;
use fs_euler_disc_e2e::render_checkpoint::{
    EULER_RENDER_CHECKPOINT_ARTIFACT_KIND, EulerRenderCheckpointError,
    EulerRenderCheckpointProvenance, euler_render_checkpoint_frame_identity,
    restore_adaptive_render_checkpoint, restore_uniform_render_checkpoint,
    store_adaptive_render_checkpoint, store_uniform_render_checkpoint,
    try_adaptive_render_checkpoint_binding, try_uniform_render_checkpoint_binding,
};
use fs_euler_disc_e2e::render_scene_bridge::{
    EulerCinematicScene, EulerFrameRequest, EulerSceneConfig, EulerTessellationConfig,
    euler_scene_smoke_settings,
};
use fs_euler_disc_e2e::specimen::{DiscProfileSpec, ResolvedDiscProfile};
use fs_euler_disc_e2e::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact,
    ExposureEventPolicy, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
    RenderContactBranch, RenderMassProperties, RenderSampleDisposition, RenderTrajectory,
    RenderTrajectoryAuthority, RenderTrajectoryCodecBudget, RenderTrajectoryMetadata,
    RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_ledger::Ledger;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_render::camera::{AnimatedCamera, Aperture, CameraProjection, CutSide, PhysicalCamera};
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{
    AdaptiveFilm, AdaptiveSamplingConfig, Film, RenderCheckpointError, RenderCheckpointWriteError,
    RenderExecutionConfig, RenderProgress,
};
use fs_rep_frep::SquatDiscEdgeTreatment;

const END_TIME_S: f64 = 0.02;
const CHECKPOINT_MAX_BYTES: u64 = 4 << 20;
static LEDGER_NONCE: AtomicU64 = AtomicU64::new(0);

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4348_4543_4b50_4f49,
                kernel_id: 0x4555_4c45_525f_4c36,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn with_cancelled_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    gate.request();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4348_4543_4b50_4f49,
                kernel_id: 0x4555_4c45_525f_4c36,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn identity(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.euler-render-checkpoint.v1",
        label.as_bytes(),
    )
}

fn ledger_path(label: &str) -> PathBuf {
    let nonce = LEDGER_NONCE.fetch_add(1, Ordering::Relaxed);
    let epoch_nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "frankensim-euler-render-checkpoint-{label}-{}-{epoch_nanos}-{nonce}.db",
        std::process::id()
    ))
}

fn specimen(cx: &Cx<'_>) -> ResolvedDiscProfile {
    DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    }
    .resolve(7_800.0, cx)
    .expect("real circular-filleted steel specimen")
}

fn render_mass(specimen: &ResolvedDiscProfile) -> MassProperties {
    MassProperties::new(
        specimen.mass_properties.mass,
        Vec3::ZERO,
        Vec3::new(
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.transverse,
            specimen.mass_properties.principal_inertia.axial,
        ),
    )
    .expect("resolved specimen mass properties")
}

fn state(mass: MassProperties) -> RigidBodyState {
    let orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 1.0).expect("tilted pose");
    RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, 0.045), orientation).expect("finite pose"),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 4.0 * mass.principal_inertia_body().z),
    )
    .expect("finite rigid-body state")
}

fn sample(
    interval_start_time_s: f64,
    time_s: f64,
    specimen: &ResolvedDiscProfile,
    mass: MassProperties,
    disposition: RenderSampleDisposition,
    cx: &Cx<'_>,
) -> RenderTrajectorySampleInput {
    let state = state(mass);
    let orientation = state.pose().orientation();
    let contact = fs_euler_disc_e2e::profile_contact_geometry(
        &specimen.chart,
        specimen.mass_properties,
        state.pose(),
        cx,
    )
    .expect("exact open-state support geometry");
    RenderTrajectorySampleInput {
        interval_start_time_s,
        time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: state.pose().position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: RenderContactBranch::Open,
        contact_geometry: None,
        signed_gap_m: contact.contact.gap_m,
        interval_contact_active: false,
        interval_normal_force_n: 0.0,
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 1.0,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois::from_state(state, mass, 0.0).expect("finite Euler QoIs"),
        disposition,
        terminal_event: None,
    }
}

fn artifact(specimen: &ResolvedDiscProfile, cx: &Cx<'_>) -> EulerRenderTrajectoryArtifact {
    let mass = render_mass(specimen);
    let first = sample(
        0.0,
        0.0,
        specimen,
        mass,
        RenderSampleDisposition::Continue,
        cx,
    );
    let last = sample(
        0.0,
        END_TIME_S,
        specimen,
        mass,
        RenderSampleDisposition::HorizonCensored,
        cx,
    );
    let identities = specimen.content_identities();
    let trajectory = RenderTrajectory::try_new(
        RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: identities.chart,
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: mass,
            },
            initial_state: state(mass),
            initial_base_mode: first.base_mode.expect("fixture base state"),
            base_model_identity: identity("base"),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: identity("model"),
            channel_availability: RenderChannelAvailability::NONE_AVAILABLE,
            configuration_identity: identity("trajectory-configuration"),
            configuration_fingerprint: 0x4348_4543_4b50_4f49,
            timestep_s: END_TIME_S,
            producer_version: "render-checkpoint-test-v1".into(),
            applicability: "deterministic visualization checkpoint fixture only".into(),
            no_claims: vec!["render replay does not validate mechanics".into()],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        },
        vec![first, last],
    )
    .expect("valid render trajectory");
    EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("campaign"),
        trajectory,
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .expect("canonical trajectory artifact")
}

fn camera_with_projection(projection: CameraProjection) -> AnimatedCamera {
    let eye = Point3::new(0.24, -0.30, 0.18);
    let target = Point3::new(0.0, 0.0, 0.025);
    let physical = PhysicalCamera::try_look_at(
        eye,
        target,
        GeomVec3::new(0.0, 0.0, 1.0),
        projection,
        target.delta_from(eye).norm(),
        Aperture::try_circular(0.0).expect("pinhole"),
    )
    .expect("scene camera");
    AnimatedCamera::try_static(0x4555_4c45_525f_434b, 0.0, END_TIME_S, physical)
        .expect("static shot")
}

fn camera() -> AnimatedCamera {
    camera_with_projection(CameraProjection::try_half_tangent(0.48).expect("projection"))
}

fn config() -> EulerSceneConfig {
    let mut config = EulerSceneConfig::reference(camera());
    config.tessellation = EulerTessellationConfig {
        azimuthal_segments: 16,
        arc_subdivisions_per_arc: 4,
    };
    config
}

fn request(cut_side: CutSide) -> EulerFrameRequest {
    EulerFrameRequest {
        frame_time_s: 0.0,
        exposure_duration_s: 0.0,
        convention: ShutterConvention::Centered,
        distribution: ShutterDistribution::UniformCounterV1,
        event_policy: ExposureEventPolicy::Refuse,
        cut_side,
    }
}

fn execution(run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(2, 2, 1, 16 << 20, RunId(run))
        .expect("bounded deterministic render execution")
}

fn provenance() -> EulerRenderCheckpointProvenance {
    EulerRenderCheckpointProvenance::try_root(
        identity("producer-build"),
        identity("producer-claim"),
    )
    .expect("explicit root checkpoint provenance")
}

fn assert_film_bits_eq(actual: &Film, expected: &Film, context: &str) {
    assert_eq!(
        (
            actual.width,
            actual.height,
            actual.spp_done,
            actual.time_mode
        ),
        (
            expected.width,
            expected.height,
            expected.spp_done,
            expected.time_mode
        ),
        "{context}: film metadata"
    );
    assert_eq!(actual.xyz.len(), expected.xyz.len(), "{context}: pixels");
    for (pixel, (actual, expected)) in actual.xyz.iter().zip(&expected.xyz).enumerate() {
        for channel in 0..3 {
            assert_eq!(
                actual[channel].to_bits(),
                expected[channel].to_bits(),
                "{context}: pixel={pixel} channel={channel}"
            );
        }
    }
}

fn assert_adaptive_bits_eq(actual: &AdaptiveFilm, expected: &AdaptiveFilm, context: &str) {
    assert_eq!(
        (
            actual.width(),
            actual.height(),
            actual.maximum_samples(),
            actual.policy(),
            actual.sampler(),
            actual.stream_seed(),
            actual.semantics_version(),
            actual.time_mode(),
        ),
        (
            expected.width(),
            expected.height(),
            expected.maximum_samples(),
            expected.policy(),
            expected.sampler(),
            expected.stream_seed(),
            expected.semantics_version(),
            expected.time_mode(),
        ),
        "{context}: adaptive identity"
    );
    assert_eq!(
        actual.sample_counts(),
        expected.sample_counts(),
        "{context}"
    );
    assert_eq!(actual.decisions(), expected.decisions(), "{context}");
    for pixel in 0..actual.xyz_sums().len() {
        for (label, actual, expected) in [
            ("sum", actual.xyz_sums()[pixel], expected.xyz_sums()[pixel]),
            (
                "mean",
                actual.running_means_xyz()[pixel],
                expected.running_means_xyz()[pixel],
            ),
            ("m2", actual.m2_xyz()[pixel], expected.m2_xyz()[pixel]),
        ] {
            for channel in 0..3 {
                assert_eq!(
                    actual[channel].to_bits(),
                    expected[channel].to_bits(),
                    "{context}: {label} pixel={pixel} channel={channel}"
                );
            }
        }
    }
}

fn assert_strict_partial_progress(progress: RenderProgress, context: &str) {
    assert!(
        progress.committed_tile_rows > 0,
        "{context}: checkpoint must contain real rendered progress"
    );
    assert!(
        progress.committed_tile_rows < progress.total_tile_rows,
        "{context}: checkpoint must remain strictly partial"
    );
    assert_eq!(progress.attempts, 1, "{context}: one bounded attempt");
}

#[test]
fn g3_ledger_reopen_restores_uniform_and_adaptive_jobs_exactly() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let prepared = scene
            .prepare_frame(request(CutSide::After))
            .expect("resolved Euler frame");
        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0001);
        let uniform_pending = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("uniform pending render");
        let uniform_binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared,
            0,
            &settings,
            &execution,
            &uniform_pending,
            provenance(),
            cx,
        )
        .expect("uniform checkpoint binding");
        let adaptive =
            AdaptiveSamplingConfig::try_new(2, 2, 0.0, 0.0, 0.0).expect("adaptive policy");
        let adaptive_pending = scene
            .begin_segment_adaptive_render(&prepared, 0, settings, adaptive, execution.clone(), cx)
            .expect("adaptive pending render");
        let adaptive_binding = try_adaptive_render_checkpoint_binding(
            &scene,
            &prepared,
            0,
            &settings,
            &execution,
            adaptive,
            &adaptive_pending,
            provenance(),
            cx,
        )
        .expect("adaptive checkpoint binding");
        assert_ne!(uniform_binding, adaptive_binding);

        let path = ledger_path("reopen");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let (uniform_artifact, uniform_content, adaptive_artifact, adaptive_content) = {
            let mut ledger = Ledger::open(path).expect("create checkpoint ledger");
            let uniform_yield = uniform_pending
                .advance_to_safe_point(cx, NonZeroU32::MIN)
                .expect("advance uniform render to a strict safe point");
            let uniform_progress = uniform_yield.progress();
            assert_strict_partial_progress(uniform_progress, "uniform yield");
            let uniform_pending = uniform_yield.into_pending();
            let uniform = store_uniform_render_checkpoint(
                &mut ledger,
                &uniform_pending,
                uniform_binding,
                CHECKPOINT_MAX_BYTES,
                cx,
            )
            .expect("store uniform checkpoint");
            assert_eq!(uniform.checkpoint().progress(), uniform_progress);
            let adaptive_yield = adaptive_pending
                .advance_to_safe_point(cx, NonZeroU32::MIN)
                .expect("advance adaptive render to a strict safe point");
            let adaptive_progress = adaptive_yield.progress();
            assert_strict_partial_progress(adaptive_progress, "adaptive yield");
            let adaptive_pending = adaptive_yield.into_pending();
            let adaptive = store_adaptive_render_checkpoint(
                &mut ledger,
                &adaptive_pending,
                adaptive_binding,
                CHECKPOINT_MAX_BYTES,
                cx,
            )
            .expect("store adaptive checkpoint");
            assert_eq!(adaptive.checkpoint().progress(), adaptive_progress);
            (
                uniform.artifact().hash,
                uniform.checkpoint().content_hash(),
                adaptive.artifact().hash,
                adaptive.checkpoint().content_hash(),
            )
        };

        let ledger = Ledger::open(path).expect("reopen checkpoint ledger");
        let uniform_seed = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("fresh uniform seed");
        let (uniform_restored, uniform_stored) = restore_uniform_render_checkpoint(
            &ledger,
            uniform_artifact,
            uniform_seed,
            uniform_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore uniform after ledger reopen");
        let uniform_receipt = uniform_stored.checkpoint();
        assert_eq!(uniform_receipt.content_hash(), uniform_content);
        assert_strict_partial_progress(uniform_receipt.progress(), "restored uniform receipt");
        let restored_uniform = uniform_restored
            .resume(cx)
            .expect("finish restored uniform");
        let reference_uniform = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("uniform reference pending")
            .resume(cx)
            .expect("finish uninterrupted uniform");
        assert_film_bits_eq(
            &restored_uniform.film,
            &reference_uniform.film,
            "ledger-restored uniform",
        );

        let adaptive_seed = scene
            .begin_segment_adaptive_render(&prepared, 0, settings, adaptive, execution.clone(), cx)
            .expect("fresh adaptive seed");
        let (adaptive_restored, adaptive_stored) = restore_adaptive_render_checkpoint(
            &ledger,
            adaptive_artifact,
            adaptive_seed,
            adaptive_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore adaptive after ledger reopen");
        let adaptive_receipt = adaptive_stored.checkpoint();
        assert_eq!(adaptive_receipt.content_hash(), adaptive_content);
        assert_strict_partial_progress(adaptive_receipt.progress(), "restored adaptive receipt");
        let restored_adaptive = adaptive_restored
            .resume(cx)
            .expect("finish restored adaptive");
        let reference_adaptive = scene
            .begin_segment_adaptive_render(&prepared, 0, settings, adaptive, execution, cx)
            .expect("adaptive reference pending")
            .resume(cx)
            .expect("finish uninterrupted adaptive");
        assert_adaptive_bits_eq(
            &restored_adaptive.film,
            &reference_adaptive.film,
            "ledger-restored adaptive",
        );
    });
}

#[test]
fn g0_cross_wired_pending_binding_and_wrong_ledger_kind_refuse() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let prepared = scene
            .prepare_frame(request(CutSide::After))
            .expect("resolved Euler frame");
        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0010);
        let pending = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("uniform checkpoint job");
        let binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared,
            0,
            &settings,
            &execution,
            &pending,
            provenance(),
            cx,
        )
        .expect("binding from admitted uniform job");

        let mut changed_settings = settings;
        changed_settings.seed ^= 1;
        let changed_pending = scene
            .begin_segment_render(&prepared, 0, changed_settings, execution.clone(), cx)
            .expect("different admitted render job");
        assert!(matches!(
            try_uniform_render_checkpoint_binding(
                &scene,
                &prepared,
                0,
                &settings,
                &execution,
                &changed_pending,
                provenance(),
                cx,
            ),
            Err(EulerRenderCheckpointError::PendingJobMismatch(_))
        ));

        let path = ledger_path("cross-wire-kind");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create refusal ledger");
        assert!(matches!(
            store_uniform_render_checkpoint(
                &mut ledger,
                &changed_pending,
                binding,
                CHECKPOINT_MAX_BYTES,
                cx,
            ),
            Err(EulerRenderCheckpointError::PendingJobMismatch(_))
        ));
        assert_eq!(
            ledger.table_count("artifacts").expect("artifact count"),
            0,
            "cross-wired pending/binding must publish nothing"
        );

        let wrong_kind = ledger
            .put_artifact(
                "generic-render-bytes",
                b"stored under the wrong artifact contract",
                None,
            )
            .expect("store deliberately wrong-kind artifact");
        let fresh = scene
            .begin_segment_render(&prepared, 0, settings, execution, cx)
            .expect("fresh correct restore target");
        match restore_uniform_render_checkpoint(
            &ledger,
            wrong_kind.hash,
            fresh,
            binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        ) {
            Err(EulerRenderCheckpointError::ArtifactKindMismatch { expected, actual }) => {
                assert_eq!(expected, EULER_RENDER_CHECKPOINT_ARTIFACT_KIND);
                assert_eq!(actual, "generic-render-bytes");
            }
            other => panic!("wrong-kind checkpoint restore produced {other:?}"),
        }
    });
}

#[test]
fn g3_root_successor_reopen_restores_successor_and_preserves_root() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let prepared = scene
            .prepare_frame(request(CutSide::After))
            .expect("resolved Euler frame");
        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0011);
        let pending = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("root render job")
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("advance root generation to a strict safe point")
            .into_pending();
        assert_strict_partial_progress(pending.progress(), "root pending");
        let root_binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared,
            0,
            &settings,
            &execution,
            &pending,
            provenance(),
            cx,
        )
        .expect("root binding from admitted job");

        let path = ledger_path("root-successor");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create generation ledger");
        let root = store_uniform_render_checkpoint(
            &mut ledger,
            &pending,
            root_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store root generation");
        let root_artifact = root.artifact().hash;
        let root_content = root.checkpoint().content_hash();
        let root_bytes = ledger
            .get_artifact_bounded(&root_artifact, CHECKPOINT_MAX_BYTES)
            .expect("read stored root")
            .expect("stored root exists");
        drop(ledger);

        let mut ledger = Ledger::open(path).expect("reopen ledger before deriving successor");
        assert_eq!(
            ledger
                .get_artifact_bounded(&root_artifact, CHECKPOINT_MAX_BYTES)
                .expect("read immutable root after first reopen")
                .expect("stored root remains readable"),
            root_bytes
        );
        let root_seed = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("fresh root restore target");
        let (restored_root, reopened_root) = restore_uniform_render_checkpoint(
            &ledger,
            root_artifact,
            root_seed,
            root_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore root after ledger reopen");
        assert_eq!(reopened_root.artifact().hash, root_artifact);
        assert_eq!(reopened_root.checkpoint().content_hash(), root_content);
        let successor_provenance = EulerRenderCheckpointProvenance::try_successor(
            identity("successor-build"),
            identity("successor-claim"),
            reopened_root,
        )
        .expect("derive successor from the reopened, restored root");
        let successor_pending = restored_root
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("advance successor generation")
            .into_pending();
        assert_eq!(
            successor_pending.progress().completed_tiles,
            successor_pending.progress().total_tiles,
            "second row quota completes this compact frame"
        );
        let successor_binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared,
            0,
            &settings,
            &execution,
            &successor_pending,
            successor_provenance,
            cx,
        )
        .expect("successor binding retains the root lineage");
        assert_eq!(successor_binding.renderer().generation(), 1);
        assert_eq!(
            successor_binding.renderer().predecessor_checkpoint(),
            Some(root_content)
        );
        let successor = store_uniform_render_checkpoint(
            &mut ledger,
            &successor_pending,
            successor_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store successor generation");
        let successor_artifact = successor.artifact().hash;
        let successor_content = successor.checkpoint().content_hash();
        assert_ne!(root_artifact, successor_artifact);
        drop(ledger);

        let ledger = Ledger::open(path).expect("reopen ledger with both generations");
        assert_eq!(
            ledger
                .get_artifact_bounded(&root_artifact, CHECKPOINT_MAX_BYTES)
                .expect("read immutable root after reopen")
                .expect("root remains readable after successor publication"),
            root_bytes
        );
        let fresh = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("fresh successor restore target");
        let (restored, restored_successor) = restore_uniform_render_checkpoint(
            &ledger,
            successor_artifact,
            fresh,
            successor_binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore successor after ledger reopen");
        assert_eq!(restored_successor.artifact().hash, successor_artifact);
        let receipt = restored_successor.checkpoint();
        assert_eq!(receipt.content_hash(), successor_content);
        assert_eq!(receipt.binding().generation(), 1);
        assert_eq!(
            receipt.binding().predecessor_checkpoint(),
            Some(root_content)
        );
        let restored = restored
            .resume(cx)
            .expect("publish restored complete successor without retracing");
        let reference = scene
            .begin_segment_render(&prepared, 0, settings, execution, cx)
            .expect("uninterrupted reference")
            .resume(cx)
            .expect("finish uninterrupted reference");
        assert_film_bits_eq(
            &restored.film,
            &reference.film,
            "reopened successor generation",
        );
        assert_eq!(
            ledger
                .get_artifact_bounded(&root_artifact, CHECKPOINT_MAX_BYTES)
                .expect("read root after successor restore")
                .expect("successor restore must not replace the root"),
            root_bytes
        );
    });
}

#[test]
fn g3_identity_mismatch_and_failed_writer_preserve_prior_artifact() {
    with_cx(|cx| {
        assert_eq!(
            EULER_RENDER_CHECKPOINT_ARTIFACT_KIND,
            "euler-render-checkpoint-v1"
        );
        assert!(matches!(
            EulerRenderCheckpointProvenance::try_root(ContentHash([0; 32]), identity("claim")),
            Err(EulerRenderCheckpointError::InvalidProvenance(_))
        ));
        assert!(matches!(
            EulerRenderCheckpointProvenance::try_root(identity("build"), ContentHash([0; 32])),
            Err(EulerRenderCheckpointError::InvalidProvenance(_))
        ));
        let root = provenance();
        assert_eq!(root.generation(), 0);
        assert_eq!(root.predecessor(), None);

        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let mut focal_sensor_config = config();
        focal_sensor_config.camera = camera_with_projection(
            CameraProjection::try_focal_sensor(1.0, 0.96).expect("equivalent focal projection"),
        );
        let focal_sensor_scene =
            EulerCinematicScene::try_build(&artifact, &specimen, focal_sensor_config, cx)
                .expect("same framing through a different declared projection");
        assert_ne!(
            scene.source_configuration_identity(),
            focal_sensor_scene.source_configuration_identity(),
            "complete source configuration identity includes projection parameterization"
        );
        let prepared_after = scene
            .prepare_frame(request(CutSide::After))
            .expect("after-cut frame");
        let prepared_before = scene
            .prepare_frame(request(CutSide::Before))
            .expect("before-cut frame");
        assert_ne!(
            euler_render_checkpoint_frame_identity(&prepared_after, 0).expect("after identity"),
            euler_render_checkpoint_frame_identity(&prepared_before, 0).expect("before identity"),
            "cut ownership is part of canonical frame identity"
        );

        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0002);
        let pending = scene
            .begin_segment_render(&prepared_after, 0, settings, execution.clone(), cx)
            .expect("uniform pending render");
        let binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared_after,
            0,
            &settings,
            &execution,
            &pending,
            provenance(),
            cx,
        )
        .expect("uniform checkpoint binding");
        assert_eq!(binding.renderer().generation(), 0);
        assert_eq!(binding.renderer().predecessor_checkpoint(), None);
        let path = ledger_path("rollback");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create rollback ledger");
        let stored = store_uniform_render_checkpoint(
            &mut ledger,
            &pending,
            binding,
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store prior checkpoint");
        let successor = EulerRenderCheckpointProvenance::try_successor(
            identity("producer-build"),
            identity("producer-claim"),
            stored,
        )
        .expect("derive successor from stored receipt");
        assert_eq!(successor.generation(), 1);
        assert_eq!(
            successor.predecessor(),
            Some(stored.checkpoint().content_hash())
        );
        let prior_hash = stored.artifact().hash;
        let prior_bytes = ledger
            .get_artifact_bounded(&prior_hash, CHECKPOINT_MAX_BYTES)
            .expect("read prior artifact")
            .expect("prior artifact exists");

        with_cancelled_cx(|cancelled_cx| {
            assert!(matches!(
                store_uniform_render_checkpoint(
                    &mut ledger,
                    &pending,
                    binding,
                    CHECKPOINT_MAX_BYTES,
                    cancelled_cx,
                ),
                Err(EulerRenderCheckpointError::Renderer(
                    RenderCheckpointError::Cancelled
                ))
            ));
        });
        assert_eq!(
            ledger.table_count("artifacts").expect("artifact count"),
            1,
            "cancelled store must not publish a second artifact"
        );

        let mut writer = ledger
            .artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)
            .expect("begin deliberately failing writer");
        let mut callback_count = 0_u32;
        let failure =
            pending.write_checkpoint(binding.renderer(), CHECKPOINT_MAX_BYTES, cx, |chunk| {
                callback_count += 1;
                writer.write(chunk).expect("stage first checkpoint chunk");
                Err::<(), _>("injected sink failure")
            });
        assert!(matches!(
            failure,
            Err(RenderCheckpointWriteError::Sink("injected sink failure"))
        ));
        assert_eq!(callback_count, 1, "injected failure must run in the sink");
        drop(writer);
        drop(ledger);

        let ledger = Ledger::open(path).expect("reopen after failed writer rollback");
        assert_eq!(
            ledger
                .get_artifact_bounded(&prior_hash, CHECKPOINT_MAX_BYTES)
                .expect("read preserved prior artifact")
                .expect("prior artifact survived rollback"),
            prior_bytes
        );

        with_cancelled_cx(|cancelled_cx| {
            let cancelled_seed = scene
                .begin_segment_render(&prepared_after, 0, settings, execution.clone(), cx)
                .expect("fresh pending job for cancelled restore");
            assert!(matches!(
                restore_uniform_render_checkpoint(
                    &ledger,
                    prior_hash,
                    cancelled_seed,
                    binding,
                    CHECKPOINT_MAX_BYTES,
                    cancelled_cx,
                ),
                Err(EulerRenderCheckpointError::Renderer(
                    RenderCheckpointError::Cancelled
                ))
            ));
        });

        let mut changed_settings = settings;
        changed_settings.seed ^= 1;
        let changed_pending = scene
            .begin_segment_render(&prepared_after, 0, changed_settings, execution.clone(), cx)
            .expect("changed pending job");
        let changed_binding = try_uniform_render_checkpoint_binding(
            &scene,
            &prepared_after,
            0,
            &changed_settings,
            &execution,
            &changed_pending,
            provenance(),
            cx,
        )
        .expect("changed job binding");
        assert_ne!(binding, changed_binding);
        assert!(matches!(
            restore_uniform_render_checkpoint(
                &ledger,
                prior_hash,
                changed_pending,
                changed_binding,
                CHECKPOINT_MAX_BYTES,
                cx,
            ),
            Err(EulerRenderCheckpointError::Renderer(_))
        ));
    });
}
