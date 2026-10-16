//! G0/G3/G4 acceptance tests for cinematic CLI discovery and static admission.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_cli::{exit, run, run_cinematic_with_gate};
use fs_euler_disc_e2e::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact,
    RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability, RenderContactBranch,
    RenderMassProperties, RenderSampleDisposition, RenderTrajectory, RenderTrajectoryAuthority,
    RenderTrajectoryCodecBudget, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
    RenderUnitSystem, RenderWorldFrame,
};
use fs_evidence::{
    cinematic_budget::{
        CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION, CinematicQualityProfile, CinematicQualityTier,
    },
    cinematic_config::{CinematicAssetBinding, CinematicAssetInterpretation},
    cinematic_config_codec::CINEMATIC_CONFIG_DOCUMENT_SCHEMA,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const MATERIAL: &[u8] = b"cinematic CLI steel spectrum v1";
const LIGHT: &[u8] = b"cinematic CLI area-light spectrum v1";
const ENVIRONMENT: &[u8] = b"cinematic CLI environment spectrum v1";
const ABUNDANT_RESOURCES: &[&str] = &[
    "--memory-bytes",
    "137438953472",
    "--free-storage-bytes",
    "2199023255552",
    "--wall-time-s",
    "31536000",
    "--workers",
    "256",
    "--paths-per-second",
    "10000000",
];

static SCRATCH_ORDINAL: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
struct Fixture {
    directory: PathBuf,
    config: PathBuf,
    artifact_root: PathBuf,
    material: PathBuf,
}

fn id(label: &str) -> ContentHash {
    hash_domain("org.frankensim.tests.cli-cinematic.v1", label.as_bytes())
}

fn component(label: &str) -> String {
    format!("1:{}", id(label).to_hex())
}

fn asset_id(bytes: &[u8], interpretation: CinematicAssetInterpretation) -> ContentHash {
    CinematicAssetBinding::from_bytes(bytes, interpretation, 1, "relocatable".to_owned())
        .expect("fixture asset")
        .content_identity()
}

fn scratch(tag: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock after epoch")
        .as_nanos();
    let ordinal = SCRATCH_ORDINAL.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "fs-cli-cinematic-{tag}-{}-{nonce}-{ordinal}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).expect("scratch directory");
    path
}

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn write_tiny_trajectory(path: &Path) -> ContentHash {
    with_cx(|cx| {
        let mass = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0))
            .expect("mass properties");
        let orientation =
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.35).expect("tilt");
        let pose = Pose::new(Vec3::new(0.0, 0.0, 0.02), orientation).expect("pose");
        let state = RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO).expect("state");
        let base = RenderBaseModeState {
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        };
        let sample = RenderTrajectorySampleInput {
            interval_start_time_s: 0.0,
            time_s: 0.0,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            center_of_mass_world_m: state.pose().position_world(),
            orientation_body_to_world: orientation.components(),
            linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
            angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
            symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
            contact_branch: RenderContactBranch::Open,
            contact_geometry: None,
            signed_gap_m: 1.0e-3,
            interval_contact_active: false,
            interval_normal_force_n: 0.0,
            contact_transitions: Vec::new(),
            base_mode: Some(base),
            channels: Default::default(),
            mechanical_energy_j: 0.0,
            energy_defect_j: 0.0,
            qois: DerivedEulerQois::from_state(state, mass, 0.0).expect("Euler QoIs"),
            disposition: RenderSampleDisposition::HorizonCensored,
            terminal_event: None,
        };
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: id("trajectory-profile"),
            specimen_chart_identity: id("trajectory-chart"),
            mass_properties: RenderMassProperties {
                identity: id("trajectory-mass"),
                properties: mass,
            },
            initial_state: state,
            initial_base_mode: base,
            base_model_identity: id("trajectory-base"),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: id("trajectory-model"),
            channel_availability: RenderChannelAvailability::NONE_AVAILABLE,
            configuration_identity: id("trajectory-configuration"),
            configuration_fingerprint: 1,
            timestep_s: 1.0 / 240.0,
            producer_version: "fs-cli-cinematic-test-v1".to_owned(),
            applicability: "tiny static-admission fixture only".to_owned(),
            no_claims: vec!["does not validate physical fidelity".to_owned()],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        let trajectory = RenderTrajectory::try_new(metadata, vec![sample]).expect("trajectory");
        let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
            id("trajectory-campaign"),
            trajectory,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .expect("trajectory artifact");
        let expected = artifact.receipt();
        let mut file = File::create(path).expect("trajectory file");
        let written = artifact
            .write_to(&mut file, RenderTrajectoryCodecBudget::DEFAULT, cx)
            .expect("trajectory bytes");
        assert_eq!(written, expected);
        expected.artifact_identity()
    })
}

fn profile_name(tier: CinematicQualityTier) -> &'static str {
    match tier {
        CinematicQualityTier::StoryboardSmoke => "storyboard-smoke",
        CinematicQualityTier::Daily1080p => "daily-1080p",
        CinematicQualityTier::Qualification4kFrame => "qualification-4k-frame",
        CinematicQualityTier::Final4k => "final-4k",
    }
}

fn write_fixture(tag: &str, tier: CinematicQualityTier, mux: bool) -> Fixture {
    write_fixture_with_trajectory(
        tag,
        tier,
        mux,
        u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION),
        id("trajectory"),
    )
}

fn write_fixture_with_trajectory(
    tag: &str,
    tier: CinematicQualityTier,
    mux: bool,
    trajectory_version: u32,
    trajectory_identity: ContentHash,
) -> Fixture {
    let directory = scratch(tag);
    let assets = directory.join("assets");
    std::fs::create_dir_all(&assets).expect("asset directory");
    let material = assets.join("private-material-token.spectrum");
    std::fs::write(&material, MATERIAL).expect("material fixture");
    std::fs::write(assets.join("light.spectrum"), LIGHT).expect("light fixture");
    std::fs::write(assets.join("environment.spectrum"), ENVIRONMENT).expect("environment fixture");

    let profile = CinematicQualityProfile::canonical(tier).expect("profile");
    let profile_ref = format!(
        "{}:{}",
        CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION,
        profile.identity().to_hex()
    );
    let material_id = asset_id(MATERIAL, CinematicAssetInterpretation::SpectralReflectance);
    let light_id = asset_id(LIGHT, CinematicAssetInterpretation::SpectralEmission);
    let environment_id = asset_id(ENVIRONMENT, CinematicAssetInterpretation::SpectralEmission);
    let (capabilities, mux_request) = if mux {
        (
            "render,audio,quarantined-mux",
            format!("av1-opus-matroska:1:{}", id("mux-adapter").to_hex()),
        )
    } else {
        ("render,audio", "none".to_owned())
    };
    let source = format!(
        "schema={CINEMATIC_CONFIG_DOCUMENT_SCHEMA}\n\
         quality_profile={}\n\
         units=si-m-kg-s-rad\n\
         seed=7\n\
         capabilities={capabilities}\n\
         render_budget_profile={profile_ref}\n\
         audio_budget_profile={profile_ref}\n\
         trajectory={}\n\
         timeline={}\n\
         camera={}\n\
         scene_geometry={}\n\
         instance_mapping={}\n\
         renderer={}\n\
         image_pipeline={}\n\
         audio_excitation={}\n\
         sound_model={}\n\
         microphone={}\n\
         room={}\n\
         material_asset=spectral-reflectance:1:{material_id}:assets/private-material-token.spectrum\n\
         light_asset=spectral-emission:1:{light_id}:assets/light.spectrum\n\
         environment_asset=spectral-emission:1:{environment_id}:assets/environment.spectrum\n\
         artifact_namespace=euler/{tag}\n\
         artifact_root=outputs/{tag}\n\
         mux={mux_request}\n",
        profile_name(tier),
        format!("{trajectory_version}:{}", trajectory_identity.to_hex()),
        component("timeline"),
        component("camera"),
        component("scene-geometry"),
        component("instance-mapping"),
        component("renderer"),
        component("image-pipeline"),
        component("audio-excitation"),
        component("sound-model"),
        component("microphone"),
        component("room"),
    );
    let config = directory.join("reference.fscine");
    std::fs::write(&config, source).expect("config fixture");
    Fixture {
        artifact_root: directory.join("outputs").join(tag),
        directory,
        config,
        material,
    }
}

fn attach_tiny_trajectory(fixture: &Fixture) -> (PathBuf, ContentHash) {
    let trajectory = fixture.directory.join("tiny.fseultrj");
    let trajectory_identity = write_tiny_trajectory(&trajectory);
    let source = std::fs::read_to_string(&fixture.config).expect("config source");
    let source = source.replace(
        &format!(
            "trajectory={}:{}",
            EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            id("trajectory").to_hex()
        ),
        &format!(
            "trajectory={}:{}",
            EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            trajectory_identity.to_hex()
        ),
    );
    std::fs::write(&fixture.config, source).expect("trajectory-bound config");
    (trajectory, trajectory_identity)
}

fn cinematic_args(mode: &str, config: &Path, dry_run: bool, resources: bool) -> Vec<String> {
    let mut arguments = vec![
        "--json".to_owned(),
        "cinematic".to_owned(),
        mode.to_owned(),
        config.to_string_lossy().into_owned(),
        "--run-reduced".to_owned(),
    ];
    if dry_run {
        arguments.push("--dry-run".to_owned());
    }
    if resources {
        arguments.extend(ABUNDANT_RESOURCES.iter().map(|value| (*value).to_owned()));
    }
    arguments
}

fn cinematic_existing_args(
    mode: &str,
    config: &Path,
    trajectory: &Path,
    dry_run: bool,
    resources: bool,
) -> Vec<String> {
    let mut arguments = vec![
        "--json".to_owned(),
        "cinematic".to_owned(),
        mode.to_owned(),
        config.to_string_lossy().into_owned(),
        "--trajectory".to_owned(),
        trajectory.to_string_lossy().into_owned(),
    ];
    if dry_run {
        arguments.push("--dry-run".to_owned());
    }
    if resources {
        arguments.extend(ABUNDANT_RESOURCES.iter().map(|value| (*value).to_owned()));
    }
    arguments
}

#[test]
fn help_and_strict_trajectory_grammar_are_discoverable() {
    let help = run(["cinematic".to_owned(), "help".to_owned()]);
    assert_eq!(help.exit_code, exit::SUCCESS);
    assert!(help.stdout.contains("representative-4k-frame"));
    assert!(help.stdout.contains(CINEMATIC_CONFIG_DOCUMENT_SCHEMA));
    assert!(help.stdout.contains("verify/mux"));
    assert!(help.stdout.contains("--trajectory"));

    let fixture = write_fixture("grammar", CinematicQualityTier::StoryboardSmoke, false);
    let missing_source = run([
        "cinematic".to_owned(),
        "inspect".to_owned(),
        fixture.config.to_string_lossy().into_owned(),
    ]);
    assert_eq!(missing_source.exit_code, exit::USAGE);
    assert!(missing_source.stderr.contains("cinematic-cli-usage"));

    let missing_verify_source = run([
        "cinematic".to_owned(),
        "verify".to_owned(),
        fixture.config.to_string_lossy().into_owned(),
    ]);
    assert_eq!(missing_verify_source.exit_code, exit::USAGE);
    assert!(
        missing_verify_source
            .stderr
            .contains("require `--trajectory <artifact>`")
    );

    let reduced_verify_source = run([
        "cinematic".to_owned(),
        "verify".to_owned(),
        fixture.config.to_string_lossy().into_owned(),
        "--run-reduced".to_owned(),
    ]);
    assert_eq!(reduced_verify_source.exit_code, exit::USAGE);
    assert!(
        reduced_verify_source
            .stderr
            .contains("`--run-reduced` is not a valid source")
    );

    let conflicting_sources = run([
        "cinematic".to_owned(),
        "inspect".to_owned(),
        fixture.config.to_string_lossy().into_owned(),
        "--run-reduced".to_owned(),
        "--trajectory".to_owned(),
        "ignored.trajectory".to_owned(),
    ]);
    assert_eq!(conflicting_sources.exit_code, exit::USAGE);

    let unknown = run([
        "--json".to_owned(),
        "cinematic".to_owned(),
        "secret-mode".to_owned(),
    ]);
    assert_eq!(unknown.exit_code, exit::USAGE);
    assert!(unknown.stdout.contains("\"mode\":\"unknown\""));
    assert!(unknown.stderr.contains("\"mode\":\"unknown\""));
}

#[test]
fn inspect_is_deterministic_and_never_creates_the_artifact_root() {
    let fixture = write_fixture(
        "inspect-no-write",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    assert!(!fixture.artifact_root.exists());
    let arguments = cinematic_args("inspect", &fixture.config, false, false);
    let first = run(arguments.clone());
    let second = run(arguments);
    assert_eq!(first, second);
    assert_eq!(first.exit_code, exit::SUCCESS);
    assert!(first.stdout.contains("\"status\":\"inspected\""));
    assert!(first.stdout.contains("\"would_write\":false"));
    assert!(first.stdout.contains("requested-reduced-campaign"));
    assert!(first.stderr.is_empty());
    assert!(!fixture.artifact_root.exists());
    assert!(
        !first
            .stdout
            .contains(&fixture.directory.to_string_lossy()[..])
    );
}

#[test]
fn inspect_decodes_and_identity_binds_a_real_trajectory_artifact() {
    let fixture = write_fixture(
        "inspect-existing-trajectory",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    let (trajectory, trajectory_identity) = attach_tiny_trajectory(&fixture);
    let arguments = cinematic_existing_args("inspect", &fixture.config, &trajectory, false, false);
    let first = run(arguments.clone());
    let second = run(arguments.clone());
    assert_eq!(first, second);
    assert_eq!(first.exit_code, exit::SUCCESS, "{}", first.stderr);
    assert!(first.stdout.contains("\"status\":\"inspected\""));
    assert!(first.stdout.contains("\"source\":\"verified-artifact\""));
    assert!(first.stdout.contains("\"verified\":true"));
    assert!(first.stdout.contains("\"sample_count\":1"));
    assert!(first.stdout.contains("\"transition_count\":0"));
    assert!(first.stdout.contains("\"chunk_count\":1"));
    assert!(
        first
            .stdout
            .contains("\"resource_admission\":\"not-requested\"")
    );
    assert!(first.stderr.is_empty());
    assert!(!fixture.artifact_root.exists());
    assert!(
        !first
            .stdout
            .contains(&fixture.directory.to_string_lossy()[..])
    );

    let source = std::fs::read_to_string(&fixture.config).expect("bound config source");
    let stale = source.replace(
        &trajectory_identity.to_hex(),
        &id("wrong-trajectory-artifact").to_hex(),
    );
    assert_ne!(
        stale, source,
        "fixture must replace the configured identity"
    );
    std::fs::write(&fixture.config, stale).expect("stale trajectory reference");
    let refused = run(arguments);
    assert_eq!(refused.exit_code, exit::REFUSED);
    assert!(
        refused
            .stderr
            .contains("cinematic-trajectory-identity-mismatch")
    );
    assert!(!refused.stderr.contains(&trajectory.to_string_lossy()[..]));
    assert!(!fixture.artifact_root.exists());
}

#[test]
fn reduced_campaign_still_refuses_an_unsupported_trajectory_schema_reference() {
    let fixture = write_fixture_with_trajectory(
        "trajectory-version",
        CinematicQualityTier::StoryboardSmoke,
        false,
        u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION) + 1,
        id("trajectory"),
    );
    let output = run(cinematic_args("inspect", &fixture.config, false, false));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(
        output
            .stderr
            .contains("cinematic-trajectory-version-mismatch")
    );
    assert!(!fixture.artifact_root.exists());
}

#[test]
fn all_four_launch_modes_bind_the_exact_named_profile() {
    let cases = [
        (
            "storyboard",
            CinematicQualityTier::StoryboardSmoke,
            "storyboard-smoke",
            "960",
            "540",
        ),
        (
            "daily",
            CinematicQualityTier::Daily1080p,
            "daily-1080p",
            "1920",
            "1080",
        ),
        (
            "representative-4k-frame",
            CinematicQualityTier::Qualification4kFrame,
            "qualification-4k-frame",
            "3840",
            "2160",
        ),
        (
            "final",
            CinematicQualityTier::Final4k,
            "final-4k",
            "3840",
            "2160",
        ),
    ];
    for (mode, tier, expected_profile, width, height) in cases {
        let fixture = write_fixture(mode, tier, false);
        let output = run(cinematic_args(mode, &fixture.config, true, true));
        assert_eq!(output.exit_code, exit::SUCCESS, "{mode}: {}", output.stderr);
        assert!(output.stdout.contains("\"status\":\"planned\""));
        assert!(
            output
                .stdout
                .contains(&format!("\"quality_profile\":\"{expected_profile}\""))
        );
        assert!(output.stdout.contains(&format!("\"width_pixels\":{width}")));
        assert!(
            output
                .stdout
                .contains(&format!("\"height_pixels\":{height}"))
        );
        assert!(!fixture.artifact_root.exists());
    }
}

#[test]
fn final_mode_refuses_preview_profile_instead_of_downgrading() {
    let fixture = write_fixture("final-guard", CinematicQualityTier::StoryboardSmoke, false);
    let output = run(cinematic_args("final", &fixture.config, true, true));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stderr.contains("cinematic-profile-mode-conflict"));
    assert!(output.stderr.contains("config.quality_profile"));
    assert!(!fixture.artifact_root.exists());
}

#[test]
fn malformed_config_asset_substitution_and_capability_omission_refuse_cleanly() {
    let unsupported = write_fixture(
        "unsupported-schema",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    let source = std::fs::read_to_string(&unsupported.config).expect("config source");
    std::fs::write(
        &unsupported.config,
        source.replace(CINEMATIC_CONFIG_DOCUMENT_SCHEMA, "unknown.schema.v999"),
    )
    .expect("mutated schema");
    let output = run(cinematic_args("inspect", &unsupported.config, false, false));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(
        output
            .stderr
            .contains("cinematic-document-unsupported-schema")
    );

    let stale = write_fixture(
        "stale-secret-asset",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    std::fs::write(&stale.material, b"substituted secret bytes").expect("stale asset");
    let output = run(cinematic_args("inspect", &stale.config, false, false));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(
        output
            .stderr
            .contains("cinematic-document-asset-identity-mismatch")
    );
    assert!(!output.stderr.contains("private-material-token"));
    assert!(!output.stderr.contains("substituted secret bytes"));

    let missing_capability = write_fixture(
        "missing-capability",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    let source = std::fs::read_to_string(&missing_capability.config).expect("config source");
    std::fs::write(
        &missing_capability.config,
        source.replace("capabilities=render,audio", "capabilities=render"),
    )
    .expect("mutated capabilities");
    let output = run(cinematic_args(
        "inspect",
        &missing_capability.config,
        false,
        false,
    ));
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stderr.contains("config.capabilities"));
}

#[test]
fn resource_deficits_report_exact_units_bounds_and_ranked_repairs() {
    let fixture = write_fixture(
        "resource-deficit",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    let mut arguments = vec![
        "--json".to_owned(),
        "cinematic".to_owned(),
        "storyboard".to_owned(),
        fixture.config.to_string_lossy().into_owned(),
        "--run-reduced".to_owned(),
        "--dry-run".to_owned(),
    ];
    arguments.extend(
        [
            "--memory-bytes",
            "0",
            "--free-storage-bytes",
            "0",
            "--wall-time-s",
            "0",
            "--workers",
            "0",
            "--paths-per-second",
            "10000000",
        ]
        .into_iter()
        .map(str::to_owned),
    );
    let output = run(arguments);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stderr.contains("cinematic-budget-host-memory"));
    assert!(output.stderr.contains("\"unit\":\"bytes\""));
    assert!(output.stderr.contains("\"required\":"));
    assert!(output.stderr.contains("\"available\":0"));
    assert!(output.stderr.contains("\"ranked_fixes\":["));
    assert!(output.stderr.contains("increase-host-memory"));
    assert!(output.stderr.contains("increase-free-storage"));
    assert!(!fixture.artifact_root.exists());
}

#[test]
fn pre_cancelled_admission_returns_130_without_reading_or_writing() {
    let directory = scratch("cancelled");
    let config = directory.join("private-config-token.fscine");
    let artifact_root = directory.join("must-not-exist");
    let gate = CancelGate::new_clock_free();
    gate.request();
    let output = run_cinematic_with_gate(
        [
            "inspect".to_owned(),
            config.to_string_lossy().into_owned(),
            "--run-reduced".to_owned(),
        ],
        true,
        &gate,
    );
    assert_eq!(output.exit_code, exit::CANCELLED);
    assert!(output.stdout.contains("\"status\":\"cancelled\""));
    assert!(output.stderr.contains("cinematic-cancelled"));
    assert!(!output.stderr.contains("private-config-token"));
    assert!(!artifact_root.exists());
}

#[test]
fn deferred_execution_modes_fail_closed_with_their_owner_beads() {
    let final_fixture = write_fixture("deferred-final", CinematicQualityTier::Final4k, false);
    let final_output = run(cinematic_args("final", &final_fixture.config, false, true));
    assert_eq!(final_output.exit_code, exit::UNAVAILABLE);
    assert!(final_output.stdout.contains("frankensim-h7xu5.8.2"));
    assert!(final_output.stderr.contains("cinematic-stage-unavailable"));

    let preview = write_fixture(
        "deferred-other",
        CinematicQualityTier::StoryboardSmoke,
        false,
    );
    let storyboard = run(cinematic_args("storyboard", &preview.config, false, true));
    assert_eq!(storyboard.exit_code, exit::UNAVAILABLE);
    assert!(storyboard.stdout.contains("frankensim-h7xu5.8.3"));

    let (preview_trajectory, _) = attach_tiny_trajectory(&preview);
    let verify = run(cinematic_existing_args(
        "verify",
        &preview.config,
        &preview_trajectory,
        false,
        false,
    ));
    assert_eq!(verify.exit_code, exit::UNAVAILABLE);
    assert!(verify.stdout.contains("frankensim-h7xu5.8.4"));

    let mux_fixture = write_fixture("deferred-mux", CinematicQualityTier::StoryboardSmoke, true);
    let (mux_trajectory, _) = attach_tiny_trajectory(&mux_fixture);
    let mux = run(cinematic_existing_args(
        "mux",
        &mux_fixture.config,
        &mux_trajectory,
        false,
        false,
    ));
    assert_eq!(mux.exit_code, exit::UNAVAILABLE);
    assert!(mux.stdout.contains("frankensim-h7xu5.8.5"));

    let resume = run(cinematic_args("resume", &preview.config, false, true));
    assert_eq!(resume.exit_code, exit::UNAVAILABLE);
    assert!(resume.stdout.contains("frankensim-h7xu5.8.2"));
}

#[test]
fn real_binary_dispatches_cinematic_help() {
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cinematic", "help"])
        .output()
        .expect("frankensim binary executes");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("UTF-8 stdout");
    assert!(stdout.contains("frankensim.cinematic.cli-result.v1"));
    assert!(stdout.contains(CINEMATIC_CONFIG_DOCUMENT_SCHEMA));
    assert!(output.stderr.is_empty());
}
