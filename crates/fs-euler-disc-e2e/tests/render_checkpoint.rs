//! G3/E2E coverage for Design-Ledger-backed Euler render checkpoints.

use std::num::NonZeroU32;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{ChannelOwnership, ContactTransitionKind};
use fs_euler_disc_e2e::render_checkpoint::{
    EULER_RENDER_CHECKPOINT_ARTIFACT_KIND, EulerRenderCheckpointError,
    EulerRenderCheckpointExpectation, EulerRenderCheckpointProducer, begin_adaptive_checkpoint_job,
    begin_uniform_checkpoint_job, euler_render_checkpoint_frame_identity,
    restore_adaptive_render_checkpoint, restore_uniform_render_checkpoint,
    store_adaptive_render_checkpoint, store_uniform_render_checkpoint,
};
use fs_euler_disc_e2e::render_scene_bridge::{
    EulerCinematicScene, EulerFrameRequest, EulerSceneConfig, EulerSceneError,
    EulerTessellationConfig, euler_scene_smoke_settings,
};
#[cfg(feature = "render-sharding-ledger")]
use fs_euler_disc_e2e::render_sharding::{
    EULER_RENDER_SHARD_RESULT_ARTIFACT_KIND, EulerRenderFrameInput, EulerRenderShardArtifactRef,
    EulerRenderShardLimits, EulerRenderShardingError, EulerUniformRenderPlan,
    execute_uniform_render_shard, merge_uniform_render_segment_artifacts,
    store_uniform_render_shard_artifact_bytes,
};
use fs_euler_disc_e2e::specimen::{DiscProfileSpec, ResolvedDiscProfile};
use fs_euler_disc_e2e::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact,
    ExposureEventPolicy, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
    RenderContactBranch, RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderSampleDisposition, RenderSupportFeature, RenderTrajectory, RenderTrajectoryAuthority,
    RenderTrajectoryCodecBudget, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
    RenderUnitSystem, RenderWorldFrame,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_ledger::Ledger;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_render::camera::{AnimatedCamera, Aperture, CameraProjection, CutSide, PhysicalCamera};
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{
    AdaptiveFilm, AdaptiveSamplingConfig, Film, RenderCheckpointBinding, RenderCheckpointError,
    RenderExecutionConfig, RenderExecutionError, RenderProgress, TracerError,
};
#[cfg(feature = "render-sharding-ledger")]
use fs_render::tracer::{
    RenderShardError, RenderShardMergeLimits, UniformRenderShardResult, merge_uniform_shards,
};
use fs_rep_frep::SquatDiscEdgeTreatment;

const END_TIME_S: f64 = 0.02;
const EVENT_TIME_S: f64 = 0.01;
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

fn reimpact_sample(
    specimen: &ResolvedDiscProfile,
    mass: MassProperties,
    cx: &Cx<'_>,
) -> RenderTrajectorySampleInput {
    let provisional_state = state(mass);
    let provisional = fs_euler_disc_e2e::profile_contact_geometry(
        &specimen.chart,
        specimen.mass_properties,
        provisional_state.pose(),
        cx,
    )
    .expect("provisional event-fixture support geometry");
    let position = provisional_state.pose().position_world();
    let corrected_pose = Pose::new(
        Vec3::new(
            position.x,
            position.y,
            position.z - provisional.contact.gap_m,
        ),
        provisional_state.pose().orientation(),
    )
    .expect("contact-corrected event-fixture pose");
    let corrected_state = RigidBodyState::new(
        corrected_pose,
        provisional_state.linear_momentum_world(),
        provisional_state.angular_momentum_body(),
    )
    .expect("contact-corrected event-fixture state");
    let contact = fs_euler_disc_e2e::profile_contact_geometry(
        &specimen.chart,
        specimen.mass_properties,
        corrected_pose,
        cx,
    )
    .expect("exact event-fixture support geometry");
    let orientation = corrected_pose.orientation();
    RenderTrajectorySampleInput {
        interval_start_time_s: 0.0,
        time_s: END_TIME_S,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: corrected_pose.position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: corrected_state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: corrected_state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: RenderContactBranch::Closed,
        contact_geometry: Some(RenderContactGeometry {
            point_world_m: contact.contact.point_world_m,
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            support_feature: RenderSupportFeature::ProfileFeature(contact.support_source_feature),
        }),
        signed_gap_m: contact.contact.gap_m,
        interval_contact_active: true,
        interval_normal_force_n: 1.0,
        contact_transitions: vec![RenderContactTransition {
            kind: ContactTransitionKind::Reimpact,
            time_s: EVENT_TIME_S,
            bracket_start_s: EVENT_TIME_S - 0.001,
            bracket_end_s: EVENT_TIME_S + 0.001,
        }],
        base_mode: Some(RenderBaseModeState {
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 1.0,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois::from_state(corrected_state, mass, 0.0)
            .expect("finite event-fixture Euler QoIs"),
        disposition: RenderSampleDisposition::HorizonCensored,
        terminal_event: None,
    }
}

fn artifact(specimen: &ResolvedDiscProfile, cx: &Cx<'_>) -> EulerRenderTrajectoryArtifact {
    artifact_with_event(specimen, false, cx)
}

fn artifact_with_event(
    specimen: &ResolvedDiscProfile,
    with_event: bool,
    cx: &Cx<'_>,
) -> EulerRenderTrajectoryArtifact {
    let mass = render_mass(specimen);
    let first = sample(
        0.0,
        0.0,
        specimen,
        mass,
        RenderSampleDisposition::Continue,
        cx,
    );
    let last = if with_event {
        reimpact_sample(specimen, mass, cx)
    } else {
        sample(
            0.0,
            END_TIME_S,
            specimen,
            mass,
            RenderSampleDisposition::HorizonCensored,
            cx,
        )
    };
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
        identity(if with_event {
            "event-campaign"
        } else {
            "campaign"
        }),
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

fn producer() -> EulerRenderCheckpointProducer {
    EulerRenderCheckpointProducer::try_new(identity("producer-build"), identity("producer-claim"))
        .expect("explicit checkpoint producer")
}

fn successor_producer() -> EulerRenderCheckpointProducer {
    EulerRenderCheckpointProducer::try_new(identity("successor-build"), identity("successor-claim"))
        .expect("explicit successor producer")
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
        let mut uniform_job =
            begin_uniform_checkpoint_job(&scene, &prepared, 0, settings, execution.clone(), cx)
                .expect("uniform durable render job");
        let adaptive =
            AdaptiveSamplingConfig::try_new(2, 2, 0.0, 0.0, 0.0).expect("adaptive policy");
        let mut adaptive_job = begin_adaptive_checkpoint_job(
            &scene,
            &prepared,
            0,
            settings,
            adaptive,
            execution.clone(),
            cx,
        )
        .expect("adaptive durable render job");
        assert_eq!(uniform_job.head(), None);
        assert_eq!(adaptive_job.head(), None);

        let path = ledger_path("reopen");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let (uniform_artifact, uniform_content, adaptive_artifact, adaptive_content) = {
            let mut ledger = Ledger::open(path).expect("create checkpoint ledger");
            let uniform_yield = uniform_job
                .advance_to_safe_point(cx, NonZeroU32::MIN)
                .expect("advance uniform render to a strict safe point");
            let uniform_progress = uniform_yield.progress();
            assert_strict_partial_progress(uniform_progress, "uniform yield");
            uniform_job = uniform_yield.into_job();
            let uniform = store_uniform_render_checkpoint(
                &mut ledger,
                &mut uniform_job,
                producer(),
                CHECKPOINT_MAX_BYTES,
                cx,
            )
            .expect("store uniform checkpoint");
            assert_eq!(uniform.checkpoint().progress(), uniform_progress);
            assert_eq!(uniform_job.head(), Some(uniform));
            let adaptive_yield = adaptive_job
                .advance_to_safe_point(cx, NonZeroU32::MIN)
                .expect("advance adaptive render to a strict safe point");
            let adaptive_progress = adaptive_yield.progress();
            assert_strict_partial_progress(adaptive_progress, "adaptive yield");
            adaptive_job = adaptive_yield.into_job();
            let adaptive = store_adaptive_render_checkpoint(
                &mut ledger,
                &mut adaptive_job,
                producer(),
                CHECKPOINT_MAX_BYTES,
                cx,
            )
            .expect("store adaptive checkpoint");
            assert_eq!(adaptive.checkpoint().progress(), adaptive_progress);
            assert_eq!(adaptive_job.head(), Some(adaptive));
            (
                uniform.artifact().hash,
                uniform.checkpoint().content_hash(),
                adaptive.artifact().hash,
                adaptive.checkpoint().content_hash(),
            )
        };

        let ledger = Ledger::open(path).expect("reopen checkpoint ledger");
        let uniform_restored = restore_uniform_render_checkpoint(
            &ledger,
            uniform_artifact,
            &scene,
            &prepared,
            0,
            settings,
            execution.clone(),
            EulerRenderCheckpointExpectation::root(producer()),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore uniform after ledger reopen");
        let uniform_receipt = uniform_restored
            .head()
            .expect("restored uniform carries its durable head")
            .checkpoint();
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

        let adaptive_restored = restore_adaptive_render_checkpoint(
            &ledger,
            adaptive_artifact,
            &scene,
            &prepared,
            0,
            settings,
            adaptive,
            execution.clone(),
            EulerRenderCheckpointExpectation::root(producer()),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore adaptive after ledger reopen");
        let adaptive_receipt = adaptive_restored
            .head()
            .expect("restored adaptive carries its durable head")
            .checkpoint();
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
fn g0_same_shutter_different_composition_frame_and_wrong_ledger_kind_refuse() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact_with_event(&specimen, true, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("event-bearing renderable Euler scene");
        let event_composition = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s: EVENT_TIME_S,
                exposure_duration_s: END_TIME_S,
                convention: ShutterConvention::Centered,
                distribution: ShutterDistribution::UniformCounterV1,
                event_policy: ExposureEventPolicy::Subdivide,
                cut_side: CutSide::After,
            })
            .expect("event-delimited composition frame");
        assert_eq!(event_composition.segments().len(), 2);
        let singleton = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s: 0.0,
                exposure_duration_s: EVENT_TIME_S,
                convention: ShutterConvention::FrontLoaded,
                distribution: ShutterDistribution::UniformCounterV1,
                event_policy: ExposureEventPolicy::Refuse,
                cut_side: CutSide::After,
            })
            .expect("same shutter as a standalone full-weight frame");
        assert_eq!(singleton.segments().len(), 1);
        assert_eq!(
            event_composition.segments()[0].shutter(),
            singleton.segments()[0].shutter(),
            "the renderer job sees the exact same shutter"
        );
        assert_eq!(event_composition.segments()[0].duration_weight(), 0.5);
        assert_eq!(singleton.segments()[0].duration_weight(), 1.0);
        assert_ne!(
            euler_render_checkpoint_frame_identity(&event_composition, 0)
                .expect("event-composition frame identity"),
            euler_render_checkpoint_frame_identity(&singleton, 0)
                .expect("standalone frame identity"),
            "composition weight must remain part of Euler frame identity"
        );

        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0010);
        let mut job = begin_uniform_checkpoint_job(
            &scene,
            &event_composition,
            0,
            settings,
            execution.clone(),
            cx,
        )
        .expect("event-composition durable job");
        let path = ledger_path("frame-collision-kind");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create refusal ledger");
        let root = store_uniform_render_checkpoint(
            &mut ledger,
            &mut job,
            producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store exact event-composition root");
        let root_bytes = ledger
            .get_artifact_bounded(&root.artifact().hash, CHECKPOINT_MAX_BYTES)
            .expect("read exact-frame root")
            .expect("exact-frame root exists");

        assert!(matches!(
            restore_uniform_render_checkpoint(
                &ledger,
                root.artifact().hash,
                &scene,
                &singleton,
                0,
                settings,
                execution.clone(),
                EulerRenderCheckpointExpectation::root(producer()),
                CHECKPOINT_MAX_BYTES,
                cx,
            ),
            Err(EulerRenderCheckpointError::Renderer(_))
        ));
        assert_eq!(
            ledger
                .get_artifact_bounded(&root.artifact().hash, CHECKPOINT_MAX_BYTES)
                .expect("re-read root after frame refusal")
                .expect("frame refusal preserves root"),
            root_bytes
        );

        let wrong_kind = ledger
            .put_artifact(
                "generic-render-bytes",
                b"stored under the wrong artifact contract",
                None,
            )
            .expect("store deliberately wrong-kind artifact");
        match restore_uniform_render_checkpoint(
            &ledger,
            wrong_kind.hash,
            &scene,
            &event_composition,
            0,
            settings,
            execution,
            EulerRenderCheckpointExpectation::root(producer()),
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
fn g0_fresh_same_job_cannot_mint_successor_but_continued_sealed_job_can() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let prepared = scene
            .prepare_frame(request(CutSide::After))
            .expect("resolved Euler frame");
        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0012);
        let mut continued =
            begin_uniform_checkpoint_job(&scene, &prepared, 0, settings, execution.clone(), cx)
                .expect("initial sealed job");
        assert_eq!(continued.head(), None);

        let path = ledger_path("sealed-ancestry");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create ancestry ledger");
        let root = store_uniform_render_checkpoint(
            &mut ledger,
            &mut continued,
            producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("fresh job publishes only a root");
        assert_eq!(root.checkpoint().binding().generation(), 0);
        assert_eq!(root.checkpoint().binding().predecessor_checkpoint(), None);
        assert_eq!(continued.head(), Some(root));

        let mut fresh_same_job =
            begin_uniform_checkpoint_job(&scene, &prepared, 0, settings, execution, cx)
                .expect("fresh zero-progress job with identical inputs");
        assert_eq!(
            fresh_same_job.head(),
            None,
            "a same-job admission carries no predecessor authority"
        );
        let fresh_root = store_uniform_render_checkpoint(
            &mut ledger,
            &mut fresh_same_job,
            successor_producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("fresh same-job state remains a root even with a successor producer");
        assert_eq!(fresh_root.checkpoint().binding().generation(), 0);
        assert_eq!(
            fresh_root.checkpoint().binding().predecessor_checkpoint(),
            None
        );

        let yielded = continued
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("continue the exact sealed state that owns the root");
        assert_strict_partial_progress(yielded.progress(), "continued successor state");
        continued = yielded.into_job();
        assert_eq!(continued.head(), Some(root));
        let successor = store_uniform_render_checkpoint(
            &mut ledger,
            &mut continued,
            successor_producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("continued sealed state publishes the real successor");
        assert_eq!(successor.checkpoint().binding().generation(), 1);
        assert_eq!(
            successor.checkpoint().binding().predecessor_checkpoint(),
            Some(root.checkpoint().content_hash())
        );
        assert_eq!(continued.head(), Some(successor));
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
        let root_yield =
            begin_uniform_checkpoint_job(&scene, &prepared, 0, settings, execution.clone(), cx)
                .expect("root durable render job")
                .advance_to_safe_point(cx, NonZeroU32::MIN)
                .expect("advance root generation to a strict safe point");
        assert_strict_partial_progress(root_yield.progress(), "root yield");
        let mut root_job = root_yield.into_job();
        assert_eq!(root_job.head(), None);

        let path = ledger_path("root-successor");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create generation ledger");
        let root = store_uniform_render_checkpoint(
            &mut ledger,
            &mut root_job,
            producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store root generation");
        assert_eq!(root_job.head(), Some(root));
        assert_eq!(root.checkpoint().binding().generation(), 0);
        assert_eq!(root.checkpoint().binding().predecessor_checkpoint(), None);
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
        let mut restored_root = restore_uniform_render_checkpoint(
            &ledger,
            root_artifact,
            &scene,
            &prepared,
            0,
            settings,
            execution.clone(),
            EulerRenderCheckpointExpectation::root(producer()),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("restore root after ledger reopen");
        let reopened_root = restored_root
            .head()
            .expect("restored root carries predecessor authority");
        assert_eq!(reopened_root.artifact().hash, root_artifact);
        assert_eq!(reopened_root.checkpoint().content_hash(), root_content);
        let successor_yield = restored_root
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("advance restored state before successor publication");
        assert_eq!(
            successor_yield.progress().completed_tiles,
            successor_yield.progress().total_tiles,
            "second row quota completes this compact frame"
        );
        restored_root = successor_yield.into_job();
        let successor = store_uniform_render_checkpoint(
            &mut ledger,
            &mut restored_root,
            successor_producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store successor generation");
        assert_eq!(restored_root.head(), Some(successor));
        assert_eq!(successor.checkpoint().binding().generation(), 1);
        assert_eq!(
            successor.checkpoint().binding().predecessor_checkpoint(),
            Some(root_content)
        );
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
        let restore_pair_bytes = root
            .artifact()
            .len
            .checked_add(successor.artifact().len)
            .expect("checkpoint-pair byte count");
        assert!(matches!(
            restore_uniform_render_checkpoint(
                &ledger,
                successor_artifact,
                &scene,
                &prepared,
                0,
                settings,
                execution.clone(),
                EulerRenderCheckpointExpectation::successor(successor_producer(), root),
                restore_pair_bytes - 1,
                cx,
            ),
            Err(EulerRenderCheckpointError::Renderer(
                RenderCheckpointError::ByteLimitExceeded { required, limit }
            )) if required == restore_pair_bytes && limit + 1 == required
        ));
        let restored = restore_uniform_render_checkpoint(
            &ledger,
            successor_artifact,
            &scene,
            &prepared,
            0,
            settings,
            execution.clone(),
            EulerRenderCheckpointExpectation::successor(successor_producer(), root),
            restore_pair_bytes,
            cx,
        )
        .expect("restore successor after ledger reopen");
        let restored_successor = restored
            .head()
            .expect("restored successor carries its durable head");
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
fn g0_canonical_raw_regression_cannot_mint_successor_authority() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let prepared = scene
            .prepare_frame(request(CutSide::After))
            .expect("resolved Euler frame");
        let settings = euler_scene_smoke_settings(4, 3);
        let execution = execution(0x4348_4543_4b50_0013);

        let root_yield =
            begin_uniform_checkpoint_job(&scene, &prepared, 0, settings, execution.clone(), cx)
                .expect("sealed root job")
                .advance_to_safe_point(
                    cx,
                    NonZeroU32::new(2).expect("nonzero complete-tile row quota"),
                )
                .expect("advance root beyond the forged successor fixture");
        assert_eq!(
            root_yield.progress().completed_tiles,
            root_yield.progress().total_tiles,
            "root fixture must own more committed state than the forgery"
        );
        let mut root_job = root_yield.into_job();

        let path = ledger_path("forged-successor-regression");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create regression ledger");
        let root = store_uniform_render_checkpoint(
            &mut ledger,
            &mut root_job,
            producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store fully committed root");

        // This lower-level checkpoint is individually canonical and uses only
        // public safe APIs, but it regresses the named predecessor from two
        // committed rows per incomplete tile to one. Metadata alone must not
        // let it enter the sealed Euler lineage.
        let regressed = scene
            .begin_segment_render(&prepared, 0, settings, execution.clone(), cx)
            .expect("fresh raw renderer job")
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("advance raw forgery to less progress with the same attempt count")
            .into_pending();
        assert_eq!(
            regressed.progress().attempts,
            root.checkpoint().progress().attempts,
            "fixture must reach the row-state check rather than fail on attempts"
        );
        assert!(
            regressed.progress().committed_tile_rows
                < root.checkpoint().progress().committed_tile_rows
        );
        let prior = root.checkpoint().binding();
        let next_producer = successor_producer();
        let forged_binding = RenderCheckpointBinding::try_new(
            prior.source_artifact_identity(),
            prior.source_configuration_identity(),
            prior.scene_identity(),
            prior.frame_identity(),
            prior.render_job_identity(),
            next_producer.producer_build(),
            next_producer.producer_claim(),
            1,
            Some(root.checkpoint().content_hash()),
        )
        .expect("syntactically valid forged successor binding");
        let mut writer = ledger
            .artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)
            .expect("stage canonical raw forgery");
        let forged_checkpoint = regressed
            .write_checkpoint(forged_binding, CHECKPOINT_MAX_BYTES, cx, |chunk| {
                writer.write(chunk)
            })
            .expect("serialize canonical but regressed checkpoint");
        let forged_artifact = writer.finish(None).expect("publish adversarial fixture");
        assert_eq!(forged_artifact.len, forged_checkpoint.byte_len());

        assert!(matches!(
            restore_uniform_render_checkpoint(
                &ledger,
                forged_artifact.hash,
                &scene,
                &prepared,
                0,
                settings,
                execution,
                EulerRenderCheckpointExpectation::successor(successor_producer(), root),
                CHECKPOINT_MAX_BYTES,
                cx,
            ),
            Err(EulerRenderCheckpointError::Renderer(
                RenderCheckpointError::SuccessorRegression {
                    field: "tile_next_row"
                }
            ))
        ));
        assert_eq!(
            ledger
                .get_artifact_bounded(&root.artifact().hash, CHECKPOINT_MAX_BYTES)
                .expect("read root after refusal")
                .expect("root survives adversarial successor refusal")
                .len() as u64,
            root.artifact().len
        );
    });
}

#[test]
fn g3_identity_mismatch_and_cancelled_store_preserve_prior_artifact() {
    with_cx(|cx| {
        assert_eq!(
            EULER_RENDER_CHECKPOINT_ARTIFACT_KIND,
            "euler-render-checkpoint-v1"
        );
        assert!(matches!(
            EulerRenderCheckpointProducer::try_new(ContentHash([0; 32]), identity("claim")),
            Err(EulerRenderCheckpointError::InvalidProvenance(_))
        ));
        assert!(matches!(
            EulerRenderCheckpointProducer::try_new(identity("build"), ContentHash([0; 32])),
            Err(EulerRenderCheckpointError::InvalidProvenance(_))
        ));
        assert_eq!(producer().producer_build(), identity("producer-build"));
        assert_eq!(producer().producer_claim(), identity("producer-claim"));

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
        let mut job = begin_uniform_checkpoint_job(
            &scene,
            &prepared_after,
            0,
            settings,
            execution.clone(),
            cx,
        )
        .expect("uniform durable job");
        assert_eq!(job.head(), None);
        let path = ledger_path("rollback");
        let path = path.to_str().expect("UTF-8 temporary ledger path");
        let mut ledger = Ledger::open(path).expect("create rollback ledger");
        let stored = store_uniform_render_checkpoint(
            &mut ledger,
            &mut job,
            producer(),
            CHECKPOINT_MAX_BYTES,
            cx,
        )
        .expect("store prior checkpoint");
        assert_eq!(stored.checkpoint().binding().generation(), 0);
        assert_eq!(stored.checkpoint().binding().predecessor_checkpoint(), None);
        assert_eq!(job.head(), Some(stored));
        let prior_hash = stored.artifact().hash;
        let prior_bytes = ledger
            .get_artifact_bounded(&prior_hash, CHECKPOINT_MAX_BYTES)
            .expect("read prior artifact")
            .expect("prior artifact exists");

        with_cancelled_cx(|cancelled_cx| {
            assert!(matches!(
                store_uniform_render_checkpoint(
                    &mut ledger,
                    &mut job,
                    successor_producer(),
                    CHECKPOINT_MAX_BYTES,
                    cancelled_cx,
                ),
                Err(EulerRenderCheckpointError::Renderer(
                    RenderCheckpointError::Cancelled
                ))
            ));
        });
        assert_eq!(
            job.head(),
            Some(stored),
            "cancelled successor store must not advance the private head"
        );
        assert_eq!(
            ledger.table_count("artifacts").expect("artifact count"),
            1,
            "cancelled store must not publish a second artifact"
        );
        drop(ledger);

        let ledger = Ledger::open(path).expect("reopen after cancelled store rollback");
        assert_eq!(
            ledger
                .get_artifact_bounded(&prior_hash, CHECKPOINT_MAX_BYTES)
                .expect("read preserved prior artifact")
                .expect("prior artifact survived rollback"),
            prior_bytes
        );

        with_cancelled_cx(|cancelled_cx| {
            assert!(matches!(
                restore_uniform_render_checkpoint(
                    &ledger,
                    prior_hash,
                    &scene,
                    &prepared_after,
                    0,
                    settings,
                    execution.clone(),
                    EulerRenderCheckpointExpectation::root(producer()),
                    CHECKPOINT_MAX_BYTES,
                    cancelled_cx,
                ),
                Err(EulerRenderCheckpointError::Scene(
                    EulerSceneError::RenderExecution(RenderExecutionError::Tracer(
                        TracerError::Cancelled
                    ))
                ))
            ));
        });

        let mut changed_settings = settings;
        changed_settings.seed ^= 1;
        assert!(matches!(
            restore_uniform_render_checkpoint(
                &ledger,
                prior_hash,
                &scene,
                &prepared_after,
                0,
                changed_settings,
                execution,
                EulerRenderCheckpointExpectation::root(producer()),
                CHECKPOINT_MAX_BYTES,
                cx,
            ),
            Err(EulerRenderCheckpointError::Renderer(_))
        ));
    });
}

#[cfg(feature = "render-sharding-ledger")]
#[test]
fn g0_euler_render_plan_is_canonical_bounded_and_exactly_covers_one_and_many_frames() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact_with_event(&specimen, true, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("event-bearing renderable Euler scene");
        let early = scene
            .prepare_frame(request(CutSide::After))
            .expect("early singleton frame");
        let event = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s: EVENT_TIME_S,
                exposure_duration_s: END_TIME_S,
                convention: ShutterConvention::Centered,
                distribution: ShutterDistribution::UniformCounterV1,
                event_policy: ExposureEventPolicy::Subdivide,
                cut_side: CutSide::After,
            })
            .expect("event-delimited frame");
        let late = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s: END_TIME_S,
                ..request(CutSide::After)
            })
            .expect("late singleton frame");
        assert_eq!(event.segments().len(), 2, "fixture must cross one event");

        let mut settings = euler_scene_smoke_settings(5, 3);
        settings.spp = 5;
        let limits = EulerRenderShardLimits::try_new(8, 128, 1 << 20, 1_000, 4 << 20, 64 << 20)
            .expect("bounded plan limits");
        let sequence = identity("canonical-shard-sequence");
        let permuted_inputs = [
            EulerRenderFrameInput::new(30, &late),
            EulerRenderFrameInput::new(10, &early),
            EulerRenderFrameInput::new(20, &event),
        ];
        let plan = EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &permuted_inputs,
            settings,
            3,
            2,
            3,
            2,
            1,
            limits,
            cx,
        )
        .expect("canonical irregular render plan");
        let canonical_inputs = [
            EulerRenderFrameInput::new(10, &early),
            EulerRenderFrameInput::new(20, &event),
            EulerRenderFrameInput::new(30, &late),
        ];
        let replay = EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &canonical_inputs,
            settings,
            3,
            2,
            3,
            2,
            1,
            limits,
            cx,
        )
        .expect("canonical-order replay plan");
        assert_eq!(plan, replay, "input order changed canonical plan bytes");
        assert_eq!(plan.plan_identity(), replay.plan_identity());
        let plan_bytes = plan
            .encode_canonical(plan.summary().encoded_plan_bytes, cx)
            .expect("encode canonical plan at its exact byte cap");
        assert_eq!(plan_bytes.len() as u64, plan.summary().encoded_plan_bytes);
        let decoded = EulerUniformRenderPlan::decode_canonical(
            &plan_bytes,
            plan_bytes.len() as u64,
            plan.plan_identity(),
            cx,
        )
        .expect("artifact-only plan replay");
        assert_eq!(decoded, plan, "decoded plan changed canonical semantics");
        assert!(matches!(
            plan.encode_canonical(plan_bytes.len() as u64 - 1, cx),
            Err(EulerRenderShardingError::PlanByteLimit { .. })
        ));
        assert!(matches!(
            EulerUniformRenderPlan::decode_canonical(
                &plan_bytes,
                plan_bytes.len() as u64 - 1,
                plan.plan_identity(),
                cx,
            ),
            Err(EulerRenderShardingError::PlanByteLimit { .. })
        ));
        for prefix in [0, 1, 7, 8, plan_bytes.len() - 1] {
            assert!(
                EulerUniformRenderPlan::decode_canonical(
                    &plan_bytes[..prefix],
                    plan_bytes.len() as u64,
                    plan.plan_identity(),
                    cx,
                )
                .is_err(),
                "accepted truncated plan prefix {prefix}/{}",
                plan_bytes.len()
            );
        }
        let mut trailing = plan_bytes.clone();
        trailing.push(0);
        assert!(
            EulerUniformRenderPlan::decode_canonical(
                &trailing,
                trailing.len() as u64,
                plan.plan_identity(),
                cx,
            )
            .is_err(),
            "accepted a trailing plan byte"
        );
        let mut corrupt = plan_bytes.clone();
        corrupt[100] ^= 1;
        assert!(
            EulerUniformRenderPlan::decode_canonical(
                &corrupt,
                corrupt.len() as u64,
                plan.plan_identity(),
                cx,
            )
            .is_err(),
            "accepted a corrupted plan identity field"
        );
        assert_eq!(
            plan.frames()
                .iter()
                .map(|frame| frame.frame_ordinal())
                .collect::<Vec<_>>(),
            vec![10, 20, 30]
        );
        assert_eq!(
            plan.summary(),
            fs_euler_disc_e2e::render_sharding::EulerRenderPlanSummary {
                frame_count: 3,
                segment_count: 4,
                shard_count: 24,
                total_paths: 300,
                encoded_plan_bytes: plan.summary().encoded_plan_bytes,
            }
        );
        assert!(plan.summary().encoded_plan_bytes > 0);

        let one_frame = EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &[EulerRenderFrameInput::new(7, &early)],
            settings,
            3,
            2,
            3,
            2,
            1,
            limits,
            cx,
        )
        .expect("one-frame render plan");
        assert_eq!(
            (
                one_frame.summary().frame_count,
                one_frame.summary().segment_count,
                one_frame.summary().shard_count,
                one_frame.summary().total_paths,
            ),
            (1, 1, 6, 75)
        );

        let event_segments = plan
            .segments()
            .iter()
            .filter(|segment| segment.frame_ordinal() == 20)
            .collect::<Vec<_>>();
        assert_eq!(event_segments.len(), 2);
        assert_ne!(
            event_segments[0].frame_identity(),
            event_segments[1].frame_identity(),
            "event-delimited segments require distinct frame identities"
        );
        let event_zero_index = plan
            .segments()
            .iter()
            .position(|segment| segment.frame_ordinal() == 20 && segment.segment_index() == 0)
            .expect("event segment zero");
        assert_eq!(
            plan.finishing_neighbors(event_zero_index)
                .expect("event finishing neighbors")
                .map(|segment| (segment.frame_ordinal(), segment.segment_index()))
                .collect::<Vec<_>>(),
            vec![(10, 0), (30, 0)],
            "finishing dependencies include adjacent frames, not sibling segments"
        );
        assert_eq!(
            plan.finishing_neighbors(0)
                .expect("early finishing neighbors")
                .map(|segment| (segment.frame_ordinal(), segment.segment_index()))
                .collect::<Vec<_>>(),
            vec![(20, 0), (20, 1)],
            "every event-delimited neighbor segment must be explicit"
        );

        let layout = plan.tile_layout().expect("plan tile layout");
        assert_eq!(layout.tile_count(), 4);
        assert_eq!(
            layout.bounds(3).expect("irregular final tile"),
            fs_render::tracer::RenderTileBounds {
                x: 3,
                y: 2,
                width: 2,
                height: 1,
            }
        );
        assert!(
            plan.shards()
                .iter()
                .any(|shard| shard.tile_start() == 3 && shard.tile_end() == 4),
            "irregular final tile block is absent"
        );
        assert!(
            plan.shards()
                .iter()
                .any(|shard| shard.sample_start() == 4 && shard.sample_end() == 5),
            "short final sample block is absent"
        );
        for segment in plan.segments() {
            let mut visits = vec![0_u8; layout.tile_count() as usize * settings.spp as usize];
            let first = segment.first_shard() as usize;
            let end = first + segment.shard_count() as usize;
            for shard in &plan.shards()[first..end] {
                assert_eq!(shard.frame_ordinal(), segment.frame_ordinal());
                assert_eq!(shard.segment_index(), segment.segment_index());
                assert_eq!(shard.frame_identity(), segment.frame_identity());
                for tile in shard.tile_start()..shard.tile_end() {
                    for sample in shard.sample_start()..shard.sample_end() {
                        let cell = tile as usize * settings.spp as usize + sample as usize;
                        visits[cell] = visits[cell]
                            .checked_add(1)
                            .expect("coverage count must remain bounded");
                    }
                }
            }
            assert!(
                visits.iter().all(|visits| *visits == 1),
                "segment frame={} part={} has a shard gap or overlap: {visits:?}",
                segment.frame_ordinal(),
                segment.segment_index(),
            );
        }

        assert!(matches!(
            EulerUniformRenderPlan::try_new(
                &scene,
                sequence,
                &[
                    EulerRenderFrameInput::new(10, &early),
                    EulerRenderFrameInput::new(10, &late),
                ],
                settings,
                3,
                2,
                3,
                2,
                1,
                limits,
                cx,
            ),
            Err(EulerRenderShardingError::DuplicateFrameOrdinal(10))
        ));
        for duplicate_inputs in [
            [
                EulerRenderFrameInput::new(30, &late),
                EulerRenderFrameInput::new(10, &early),
                EulerRenderFrameInput::new(30, &late),
                EulerRenderFrameInput::new(10, &early),
            ],
            [
                EulerRenderFrameInput::new(10, &early),
                EulerRenderFrameInput::new(30, &late),
                EulerRenderFrameInput::new(10, &early),
                EulerRenderFrameInput::new(30, &late),
            ],
        ] {
            assert!(matches!(
                EulerUniformRenderPlan::try_new(
                    &scene,
                    sequence,
                    &duplicate_inputs,
                    settings,
                    3,
                    2,
                    3,
                    2,
                    1,
                    limits,
                    cx,
                ),
                Err(EulerRenderShardingError::DuplicateFrameOrdinal(10))
            ));
        }

        let exact_plan_bytes = plan.summary().encoded_plan_bytes;
        let exact_limits =
            EulerRenderShardLimits::try_new(3, 24, exact_plan_bytes, 26, 4 << 20, 64 << 20)
                .expect("exact plan limits");
        EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &canonical_inputs,
            settings,
            3,
            2,
            3,
            2,
            1,
            exact_limits,
            cx,
        )
        .expect("exact frame, shard, plan-byte, and path caps must pass");
        for (label, limited, expected) in [
            (
                "frames",
                EulerRenderShardLimits::try_new(2, 24, exact_plan_bytes, 26, 4 << 20, 64 << 20)
                    .expect("frame-limited policy"),
                "frame",
            ),
            (
                "shards",
                EulerRenderShardLimits::try_new(3, 23, exact_plan_bytes, 26, 4 << 20, 64 << 20)
                    .expect("shard-limited policy"),
                "shard",
            ),
            (
                "plan bytes",
                EulerRenderShardLimits::try_new(3, 24, exact_plan_bytes - 1, 26, 4 << 20, 64 << 20)
                    .expect("byte-limited policy"),
                "plan byte",
            ),
            (
                "paths",
                EulerRenderShardLimits::try_new(3, 24, exact_plan_bytes, 25, 4 << 20, 64 << 20)
                    .expect("path-limited policy"),
                "path",
            ),
        ] {
            let refusal = EulerUniformRenderPlan::try_new(
                &scene,
                sequence,
                &canonical_inputs,
                settings,
                3,
                2,
                3,
                2,
                1,
                limited,
                cx,
            )
            .expect_err(label);
            assert!(
                matches!(
                    (&refusal, expected),
                    (EulerRenderShardingError::FrameLimit { .. }, "frame")
                        | (EulerRenderShardingError::ShardLimit { .. }, "shard")
                        | (EulerRenderShardingError::PlanByteLimit { .. }, "plan byte")
                        | (EulerRenderShardingError::PathLimit { .. }, "path")
                ),
                "{label} one-short cap produced {refusal:?}"
            );
        }
    });
}

#[cfg(feature = "render-sharding-ledger")]
#[test]
fn g3_independent_workers_exchange_strict_bytes_and_reopened_ledger_replays_exact_films() {
    with_cx(|cx| {
        let specimen = specimen(cx);
        let artifact = artifact(&specimen, cx);
        let scene = EulerCinematicScene::try_build(&artifact, &specimen, config(), cx)
            .expect("renderable Euler scene");
        let first = scene
            .prepare_frame(request(CutSide::After))
            .expect("first singleton frame");
        let last = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s: END_TIME_S,
                ..request(CutSide::After)
            })
            .expect("last singleton frame");
        let inputs = [
            EulerRenderFrameInput::new(100, &first),
            EulerRenderFrameInput::new(200, &last),
        ];
        let mut settings = euler_scene_smoke_settings(4, 3);
        settings.spp = 4;
        let sequence = identity("two-worker-byte-shard-sequence");

        // Probe only to derive exact, implementation-independent caps. The
        // final plan below has its own identity because limits are authority.
        let probe_limits =
            EulerRenderShardLimits::try_new(4, 32, 1 << 20, 1 << 20, 4 << 20, 64 << 20)
                .expect("probe limits");
        let probe_plan = EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &inputs,
            settings,
            2,
            2,
            2,
            2,
            1,
            probe_limits,
            cx,
        )
        .expect("probe shard plan");
        assert_eq!(probe_plan.summary().shard_count, 8);
        let probe_results = with_cx(|worker_cx| {
            (0..probe_plan.shards().len())
                .map(|index| {
                    execute_uniform_render_shard(&probe_plan, &scene, &inputs, index, worker_cx)
                        .expect("probe shard execution")
                })
                .collect::<Vec<_>>()
        });
        let exact_paths = probe_plan
            .shards()
            .iter()
            .map(|shard| shard.path_count())
            .max()
            .expect("nonempty probe plan");
        let exact_result_bytes = probe_results
            .iter()
            .map(UniformRenderShardResult::encoded_result_bytes)
            .max()
            .expect("nonempty probe results");
        let segment_input_bytes = probe_plan
            .segments()
            .iter()
            .map(|segment| {
                let start = segment.first_shard() as usize;
                let end = start + segment.shard_count() as usize;
                probe_results[start..end]
                    .iter()
                    .map(UniformRenderShardResult::encoded_result_bytes)
                    .sum::<u64>()
            })
            .collect::<Vec<_>>();
        let exact_aggregate_bytes = *segment_input_bytes
            .iter()
            .max()
            .expect("nonempty segment byte totals");
        let exact_limits = EulerRenderShardLimits::try_new(
            2,
            8,
            probe_plan.summary().encoded_plan_bytes,
            exact_paths,
            exact_result_bytes,
            exact_aggregate_bytes,
        )
        .expect("exact worker and coordinator limits");
        let plan = EulerUniformRenderPlan::try_new(
            &scene,
            sequence,
            &inputs,
            settings,
            2,
            2,
            2,
            2,
            1,
            exact_limits,
            cx,
        )
        .expect("exact-cap two-worker plan");
        assert_eq!(plan.summary().shard_count, 8);
        assert_eq!(
            plan.shards().iter().map(|shard| shard.path_count()).max(),
            Some(exact_limits.max_paths_per_shard())
        );

        assert!(matches!(
            with_cancelled_cx(|cancelled_cx| execute_uniform_render_shard(
                &plan,
                &scene,
                &inputs,
                0,
                cancelled_cx,
            )),
            Err(EulerRenderShardingError::Cancelled)
        ));

        // Each worker has its own cancellation gate and arena pool. Their
        // only coordinator-facing product is a canonical byte vector.
        let even_payloads = with_cx(|worker_cx| {
            (0..plan.shards().len())
                .step_by(2)
                .map(|index| {
                    let result =
                        execute_uniform_render_shard(&plan, &scene, &inputs, index, worker_cx)
                            .expect("even-worker shard execution");
                    let identity = result.result_identity();
                    let bytes = result
                        .encode_canonical(exact_limits.max_result_bytes_per_shard(), worker_cx)
                        .expect("even-worker canonical bytes");
                    (index, identity, bytes)
                })
                .collect::<Vec<_>>()
        });
        let odd_payloads = with_cx(|worker_cx| {
            (1..plan.shards().len())
                .step_by(2)
                .map(|index| {
                    let result =
                        execute_uniform_render_shard(&plan, &scene, &inputs, index, worker_cx)
                            .expect("odd-worker shard execution");
                    let identity = result.result_identity();
                    let bytes = result
                        .encode_canonical(exact_limits.max_result_bytes_per_shard(), worker_cx)
                        .expect("odd-worker canonical bytes");
                    (index, identity, bytes)
                })
                .collect::<Vec<_>>()
        });
        let mut payload_slots = (0..plan.shards().len())
            .map(|_| None)
            .collect::<Vec<Option<(ContentHash, Vec<u8>)>>>();
        for (index, result_identity, bytes) in even_payloads.into_iter().chain(odd_payloads) {
            assert!(
                payload_slots[index]
                    .replace((result_identity, bytes))
                    .is_none(),
                "worker overlap at shard {index}"
            );
        }
        let payloads = payload_slots
            .into_iter()
            .enumerate()
            .map(|(index, payload)| {
                payload.unwrap_or_else(|| panic!("missing worker shard {index}"))
            })
            .collect::<Vec<_>>();

        // A separate single-worker execution is the same frozen-plan oracle;
        // it is not reconstructed from the distributed worker objects.
        let single_worker_results = with_cx(|worker_cx| {
            (0..plan.shards().len())
                .map(|index| {
                    execute_uniform_render_shard(&plan, &scene, &inputs, index, worker_cx)
                        .expect("single-worker same-plan execution")
                })
                .collect::<Vec<_>>()
        });
        for (index, ((worker_identity, worker_bytes), reference)) in
            payloads.iter().zip(&single_worker_results).enumerate()
        {
            assert_eq!(
                *worker_identity,
                reference.result_identity(),
                "worker-context result identity drift at shard {index}"
            );
            assert_eq!(
                worker_bytes.len() as u64,
                reference.encoded_result_bytes(),
                "worker-context canonical length drift at shard {index}"
            );
        }
        assert_eq!(
            single_worker_results
                .iter()
                .map(UniformRenderShardResult::encoded_result_bytes)
                .max(),
            Some(exact_limits.max_result_bytes_per_shard())
        );

        let output_bytes = u64::from(settings.width) * u64::from(settings.height) * 24;
        let direct_films = plan
            .segments()
            .iter()
            .map(|segment| {
                let start = segment.first_shard() as usize;
                let end = start + segment.shard_count() as usize;
                let specs = single_worker_results[start..end]
                    .iter()
                    .map(|result| *result.spec())
                    .collect::<Vec<_>>();
                merge_uniform_shards(
                    &specs,
                    &single_worker_results[start..end],
                    RenderShardMergeLimits::try_new(exact_aggregate_bytes, output_bytes)
                        .expect("exact direct-merge limits"),
                    cx,
                )
                .expect("single-worker same-plan direct merge")
            })
            .collect::<Vec<_>>();

        let path = ledger_path("render-shard-workers");
        let path = path.to_str().expect("UTF-8 shard ledger path").to_owned();
        let (artifact_refs, foreign_kind_hash, corrupt_hash) = {
            let ledger = Ledger::open(&path).expect("create shard coordinator ledger");

            let mut corrupted_worker_bytes = payloads[0].1.clone();
            let corrupt_index = corrupted_worker_bytes.len() / 2;
            corrupted_worker_bytes[corrupt_index] ^= 1;
            assert!(matches!(
                store_uniform_render_shard_artifact_bytes(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &corrupted_worker_bytes,
                    cx,
                ),
                Err(EulerRenderShardingError::Renderer(
                    RenderShardError::Integrity | RenderShardError::NonCanonical(_)
                ))
            ));
            assert!(
                store_uniform_render_shard_artifact_bytes(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &payloads[1].1,
                    cx,
                )
                .is_err(),
                "coordinator accepted another shard's bytes at index zero"
            );
            assert!(matches!(
                with_cancelled_cx(|cancelled_cx| store_uniform_render_shard_artifact_bytes(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &payloads[0].1,
                    cancelled_cx,
                )),
                Err(EulerRenderShardingError::Cancelled)
            ));

            let mut refs = (0..plan.shards().len())
                .map(|_| None)
                .collect::<Vec<Option<EulerRenderShardArtifactRef>>>();
            for index in (0..plan.shards().len()).rev() {
                let receipt = store_uniform_render_shard_artifact_bytes(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    index,
                    &payloads[index].1,
                    cx,
                )
                .unwrap_or_else(|error| panic!("store worker shard {index}: {error:?}"));
                assert!(!receipt.deduped, "first store deduped shard {index}");
                assert_eq!(receipt.result_identity, payloads[index].0);
                assert_eq!(receipt.len, payloads[index].1.len() as u64);
                refs[index] = Some(receipt.artifact);
            }
            let duplicate = store_uniform_render_shard_artifact_bytes(
                &ledger,
                &plan,
                &scene,
                &inputs,
                0,
                &payloads[0].1,
                cx,
            )
            .expect("idempotent duplicate publication");
            assert!(duplicate.deduped);
            assert_eq!(duplicate.artifact, refs[0].expect("stored shard zero"));

            let foreign_kind_hash = ledger
                .put_artifact("foreign-render-shard-kind", b"not a render shard", None)
                .expect("store foreign-kind fixture")
                .hash;
            let corrupt_hash = ledger
                .put_artifact(
                    EULER_RENDER_SHARD_RESULT_ARTIFACT_KIND,
                    b"truncated render shard",
                    None,
                )
                .expect("store corrupt shard-kind fixture")
                .hash;
            (
                refs.into_iter()
                    .enumerate()
                    .map(|(index, artifact)| {
                        artifact.unwrap_or_else(|| panic!("missing artifact ref {index}"))
                    })
                    .collect::<Vec<_>>(),
                foreign_kind_hash,
                corrupt_hash,
            )
        };

        let first_artifact_films = {
            let ledger = Ledger::open(&path).expect("first coordinator reopen");
            let first_segment = plan.segments()[0];
            let first_start = first_segment.first_shard() as usize;
            let first_end = first_start + first_segment.shard_count() as usize;
            let first_refs = &artifact_refs[first_start..first_end];
            assert!(matches!(
                merge_uniform_render_segment_artifacts(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &first_refs[..first_refs.len() - 1],
                    cx,
                ),
                Err(EulerRenderShardingError::MissingArtifact(_))
            ));

            let mut foreign_refs = first_refs.to_vec();
            foreign_refs[0] = EulerRenderShardArtifactRef::new(
                plan.shards()[first_start].logical_shard_identity(),
                foreign_kind_hash,
            );
            assert!(matches!(
                merge_uniform_render_segment_artifacts(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &foreign_refs,
                    cx,
                ),
                Err(EulerRenderShardingError::ForeignArtifactKind { .. })
            ));
            let mut corrupt_refs = first_refs.to_vec();
            corrupt_refs[0] = EulerRenderShardArtifactRef::new(
                plan.shards()[first_start].logical_shard_identity(),
                corrupt_hash,
            );
            assert!(matches!(
                merge_uniform_render_segment_artifacts(
                    &ledger,
                    &plan,
                    &scene,
                    &inputs,
                    0,
                    &corrupt_refs,
                    cx,
                ),
                Err(EulerRenderShardingError::Renderer(
                    RenderShardError::Truncated
                ))
            ));

            plan.segments()
                .iter()
                .enumerate()
                .map(|(segment_index, segment)| {
                    let start = segment.first_shard() as usize;
                    let end = start + segment.shard_count() as usize;
                    let mut arrival = artifact_refs[start..end].to_vec();
                    arrival.reverse();
                    arrival.push(arrival[0]);
                    merge_uniform_render_segment_artifacts(
                        &ledger,
                        &plan,
                        &scene,
                        &inputs,
                        segment_index,
                        &arrival,
                        cx,
                    )
                    .unwrap_or_else(|error| {
                        panic!("merge stored segment {segment_index}: {error:?}")
                    })
                })
                .collect::<Vec<_>>()
        };
        for (segment_index, (artifact_film, direct_film)) in
            first_artifact_films.iter().zip(&direct_films).enumerate()
        {
            assert_film_bits_eq(
                artifact_film,
                direct_film,
                &format!("artifact merge vs one-worker segment {segment_index}"),
            );
        }

        let replay_films = {
            let ledger = Ledger::open(&path).expect("second coordinator reopen");
            plan.segments()
                .iter()
                .enumerate()
                .map(|(segment_index, segment)| {
                    let start = segment.first_shard() as usize;
                    let end = start + segment.shard_count() as usize;
                    merge_uniform_render_segment_artifacts(
                        &ledger,
                        &plan,
                        &scene,
                        &inputs,
                        segment_index,
                        &artifact_refs[start..end],
                        cx,
                    )
                    .unwrap_or_else(|error| {
                        panic!("artifact-only replay segment {segment_index}: {error:?}")
                    })
                })
                .collect::<Vec<_>>()
        };
        for (segment_index, (first_merge, replay)) in
            first_artifact_films.iter().zip(&replay_films).enumerate()
        {
            assert_film_bits_eq(
                first_merge,
                replay,
                &format!("ledger-reopen replay segment {segment_index}"),
            );
        }

        assert!(matches!(
            EulerUniformRenderPlan::try_new(
                &scene,
                sequence,
                &inputs,
                settings,
                2,
                2,
                2,
                2,
                1,
                EulerRenderShardLimits::try_new(
                    2,
                    8,
                    probe_plan.summary().encoded_plan_bytes,
                    exact_paths,
                    exact_result_bytes - 1,
                    exact_aggregate_bytes,
                )
                .expect("one-short result policy"),
                cx,
            ),
            Err(EulerRenderShardingError::Renderer(
                RenderShardError::ResultByteLimit { .. }
            ))
        ));
        assert!(matches!(
            EulerUniformRenderPlan::try_new(
                &scene,
                sequence,
                &inputs,
                settings,
                2,
                2,
                2,
                2,
                1,
                EulerRenderShardLimits::try_new(
                    2,
                    8,
                    probe_plan.summary().encoded_plan_bytes,
                    exact_paths,
                    exact_result_bytes,
                    exact_aggregate_bytes - 1,
                )
                .expect("one-short aggregate policy"),
                cx,
            ),
            Err(EulerRenderShardingError::AggregateResultByteLimit { .. })
        ));
    });
}
