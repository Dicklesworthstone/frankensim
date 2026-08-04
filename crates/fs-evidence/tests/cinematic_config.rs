//! G0/G3 configuration identity, invalidation, and fail-closed tests.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::cinematic_config::{
    CINEMATIC_CONFIG_SCHEMA_VERSION, CinematicArtifactRoot, CinematicAssetBinding,
    CinematicAssetInterpretation, CinematicCapabilities, CinematicComponentRef,
    CinematicComponentRole, CinematicConfig, CinematicConfigError, CinematicConfigInput,
    CinematicConfigUnits, CinematicMuxCodec, CinematicMuxRequest,
};

fn id(label: &str) -> ContentHash {
    hash_domain("org.frankensim.tests.cinematic-config.v1", label.as_bytes())
}

fn component(role: CinematicComponentRole, label: &str) -> CinematicComponentRef {
    CinematicComponentRef::try_new(role, id(label), 1).expect("valid component")
}

fn asset(label: &str, interpretation: CinematicAssetInterpretation) -> CinematicAssetBinding {
    CinematicAssetBinding::from_bytes(label.as_bytes(), interpretation, 1, format!("/a/{label}"))
        .expect("valid asset")
}

fn input() -> CinematicConfigInput {
    CinematicConfigInput {
        schema_version: CINEMATIC_CONFIG_SCHEMA_VERSION,
        units: Some(CinematicConfigUnits::SiMetersKilogramsSecondsRadians),
        seed: Some(42),
        capabilities: Some(
            CinematicCapabilities::try_new(
                CinematicCapabilities::RENDER
                    | CinematicCapabilities::AUDIO
                    | CinematicCapabilities::QUARANTINED_MUX,
            )
            .expect("capabilities"),
        ),
        render_budget_profile: Some(component(
            CinematicComponentRole::RenderBudgetProfile,
            "render-budget",
        )),
        audio_budget_profile: Some(component(
            CinematicComponentRole::AudioBudgetProfile,
            "audio-budget",
        )),
        trajectory: component(CinematicComponentRole::Trajectory, "trajectory"),
        timeline: component(CinematicComponentRole::Timeline, "timeline"),
        camera: component(CinematicComponentRole::Camera, "camera"),
        scene_geometry: component(CinematicComponentRole::SceneGeometry, "scene"),
        instance_mapping: component(CinematicComponentRole::InstanceMapping, "instances"),
        renderer: component(CinematicComponentRole::Renderer, "renderer"),
        image_pipeline: component(CinematicComponentRole::ImagePipeline, "image"),
        audio_excitation: component(CinematicComponentRole::AudioExcitation, "excitation"),
        sound_model: component(CinematicComponentRole::SoundModel, "sound"),
        microphone: component(CinematicComponentRole::Microphone, "microphone"),
        room: component(CinematicComponentRole::Room, "room"),
        material_assets: vec![
            asset("steel", CinematicAssetInterpretation::SpectralReflectance),
            asset(
                "tungsten",
                CinematicAssetInterpretation::SpectralReflectance,
            ),
        ],
        light_assets: vec![asset(
            "softbox",
            CinematicAssetInterpretation::SpectralEmission,
        )],
        environment_asset: asset(
            "environment",
            CinematicAssetInterpretation::SpectralEmission,
        ),
        artifact_root: CinematicArtifactRoot::try_new(
            "euler/reference-v1".to_owned(),
            "/renders/a".to_owned(),
        )
        .expect("root"),
        mux_request: CinematicMuxRequest::QuarantinedAdapter {
            adapter_identity: id("mux-adapter"),
            adapter_version: 1,
            codec: CinematicMuxCodec::Av1OpusMatroska,
        },
    }
}

fn config(input: CinematicConfigInput) -> CinematicConfig {
    CinematicConfig::try_new(input).expect("configuration must admit")
}

#[test]
fn equivalent_asset_order_and_relocated_paths_have_identical_composition() {
    let first = config(input());
    let mut relocated = input();
    relocated.material_assets.reverse();
    relocated.material_assets[0] = relocated.material_assets[0]
        .with_locator_hint("/different/machine/asset".to_owned())
        .expect("relocated asset");
    relocated.artifact_root = CinematicArtifactRoot::try_new(
        "euler/reference-v1".to_owned(),
        "/different/output/root".to_owned(),
    )
    .expect("relocated root");
    let second = config(relocated);
    assert_eq!(first.composition_identity(), second.composition_identity());
    assert_eq!(first.image_identity(), second.image_identity());
    assert_eq!(first.canonical_bytes(), second.canonical_bytes());
}

#[test]
fn one_field_invalidation_matrix_is_minimal_and_never_under_invalidates() {
    let base = config(input());

    let mut changed = input();
    changed.camera = component(CinematicComponentRole::Camera, "camera-2");
    let camera = config(changed);
    assert_eq!(base.trajectory_identity(), camera.trajectory_identity());
    assert_ne!(base.image_identity(), camera.image_identity());
    assert_eq!(base.audio_identity(), camera.audio_identity());
    assert_ne!(base.mux_identity(), camera.mux_identity());

    let mut changed = input();
    changed.material_assets[0] =
        asset("steel-2", CinematicAssetInterpretation::SpectralReflectance);
    let material = config(changed);
    assert_eq!(base.trajectory_identity(), material.trajectory_identity());
    assert_ne!(base.image_identity(), material.image_identity());
    assert_eq!(base.audio_identity(), material.audio_identity());

    let mut changed = input();
    changed.room = component(CinematicComponentRole::Room, "room-2");
    let room = config(changed);
    assert_eq!(base.trajectory_identity(), room.trajectory_identity());
    assert_eq!(base.image_identity(), room.image_identity());
    assert_ne!(base.audio_identity(), room.audio_identity());

    let mut changed = input();
    changed.trajectory = component(CinematicComponentRole::Trajectory, "trajectory-2");
    let trajectory = config(changed);
    assert_ne!(base.trajectory_identity(), trajectory.trajectory_identity());
    assert_ne!(base.image_identity(), trajectory.image_identity());
    assert_ne!(base.audio_identity(), trajectory.audio_identity());
    assert_ne!(
        base.composition_identity(),
        trajectory.composition_identity()
    );
}

#[test]
fn artifact_root_and_mux_changes_only_invalidate_their_downstream_products() {
    let base = config(input());
    let mut changed = input();
    changed.artifact_root = CinematicArtifactRoot::try_new(
        "euler/alternate-delivery".to_owned(),
        "/renders/a".to_owned(),
    )
    .expect("root");
    let root = config(changed);
    assert_eq!(base.trajectory_identity(), root.trajectory_identity());
    assert_eq!(base.image_identity(), root.image_identity());
    assert_eq!(base.audio_identity(), root.audio_identity());
    assert_eq!(base.mux_identity(), root.mux_identity());
    assert_ne!(base.composition_identity(), root.composition_identity());

    let mut changed = input();
    changed.mux_request = CinematicMuxRequest::QuarantinedAdapter {
        adapter_identity: id("mux-adapter-2"),
        adapter_version: 1,
        codec: CinematicMuxCodec::H265PcmQuickTime,
    };
    let mux = config(changed);
    assert_eq!(base.image_identity(), mux.image_identity());
    assert_eq!(base.audio_identity(), mux.audio_identity());
    assert_ne!(base.mux_identity(), mux.mux_identity());
    assert_ne!(base.composition_identity(), mux.composition_identity());
}

#[test]
fn five_explicits_unknown_versions_and_role_cross_wiring_refuse() {
    let mut candidate = input();
    candidate.units = None;
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingUnits)
    );
    let mut candidate = input();
    candidate.seed = None;
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingSeed)
    );
    let mut candidate = input();
    candidate.capabilities = None;
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingCapabilities)
    );
    let mut candidate = input();
    candidate.render_budget_profile = None;
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingBudget)
    );
    let mut candidate = input();
    candidate.audio_budget_profile = None;
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingBudget)
    );
    let mut candidate = input();
    candidate.schema_version += 1;
    assert!(matches!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::UnsupportedSchemaVersion(_))
    ));
    let mut candidate = input();
    candidate.camera = component(CinematicComponentRole::Room, "wrong-role");
    assert!(matches!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::ComponentRoleMismatch { .. })
    ));
}

#[test]
fn zero_versions_unknown_capabilities_and_zero_identities_refuse_at_the_leaf() {
    assert_eq!(
        CinematicComponentRef::try_new(CinematicComponentRole::Camera, id("camera"), 0),
        Err(CinematicConfigError::InvalidComponentVersion(
            CinematicComponentRole::Camera
        ))
    );
    assert!(matches!(
        CinematicCapabilities::try_new(1 << 31),
        Err(CinematicConfigError::InvalidCapabilities(_))
    ));
    assert_eq!(
        CinematicComponentRef::try_new(CinematicComponentRole::Camera, ContentHash([0; 32]), 1,),
        Err(CinematicConfigError::MissingContentIdentity)
    );
    assert_eq!(
        CinematicAssetBinding::from_bytes(
            b"asset",
            CinematicAssetInterpretation::LinearSceneTexture,
            0,
            "/asset".to_owned(),
        ),
        Err(CinematicConfigError::InvalidAssetVersion)
    );

    let mut candidate = input();
    candidate.mux_request = CinematicMuxRequest::QuarantinedAdapter {
        adapter_identity: id("mux"),
        adapter_version: 0,
        codec: CinematicMuxCodec::Av1OpusMatroska,
    };
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::InvalidMuxAdapterVersion)
    );
}

#[test]
fn unused_mux_capability_does_not_invalidate_image_audio_or_composition() {
    let mut with_mux_capability = input();
    with_mux_capability.mux_request = CinematicMuxRequest::None;
    let first = config(with_mux_capability);

    let mut without_mux_capability = input();
    without_mux_capability.mux_request = CinematicMuxRequest::None;
    without_mux_capability.capabilities = Some(
        CinematicCapabilities::try_new(
            CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO,
        )
        .expect("base capabilities"),
    );
    let second = config(without_mux_capability);
    assert_eq!(first.canonical_bytes(), second.canonical_bytes());
    assert_eq!(first.composition_identity(), second.composition_identity());
}

#[test]
fn composition_identity_known_answer_vector_is_stable() {
    let admitted = config(input());
    assert_eq!(
        admitted.composition_identity().to_string(),
        "b78a64017c61a41ebedad4235030ae68868ccc76d92a8c126bf004eda9eaef34"
    );
}

#[test]
fn stale_asset_bytes_duplicate_bindings_and_mux_without_capability_refuse() {
    let bound = CinematicAssetBinding::from_bytes(
        b"asset-v1",
        CinematicAssetInterpretation::LinearSceneTexture,
        1,
        "/asset".to_owned(),
    )
    .expect("binding");
    assert!(bound.verify_bytes(b"asset-v1").is_ok());
    assert_eq!(
        bound.verify_bytes(b"asset-v2"),
        Err(CinematicConfigError::AssetContentMismatch)
    );

    let mut candidate = input();
    candidate
        .material_assets
        .push(candidate.material_assets[0].clone());
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::DuplicateAsset)
    );

    let mut candidate = input();
    candidate.material_assets.clear();
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingMaterialAssets)
    );

    let mut candidate = input();
    candidate.light_assets.clear();
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingLightAssets)
    );

    let mut candidate = input();
    candidate.material_assets = (0..1025)
        .map(|index| {
            asset(
                &format!("material-{index}"),
                CinematicAssetInterpretation::SpectralReflectance,
            )
        })
        .collect();
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::TooManyAssets)
    );

    let mut candidate = input();
    candidate.capabilities = Some(
        CinematicCapabilities::try_new(
            CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO,
        )
        .expect("base capabilities"),
    );
    assert_eq!(
        CinematicConfig::try_new(candidate),
        Err(CinematicConfigError::MissingMuxCapability)
    );
}

#[test]
fn shared_explicits_invalidate_image_and_audio_but_partitioned_budgets_do_not_crosswire() {
    let base = config(input());
    for mutate in [
        |candidate: &mut CinematicConfigInput| candidate.seed = Some(43),
        |candidate: &mut CinematicConfigInput| {
            candidate.timeline = component(CinematicComponentRole::Timeline, "timeline-2");
        },
    ] {
        let mut candidate = input();
        mutate(&mut candidate);
        let changed = config(candidate);
        assert_eq!(base.trajectory_identity(), changed.trajectory_identity());
        assert_ne!(base.image_identity(), changed.image_identity());
        assert_ne!(base.audio_identity(), changed.audio_identity());
    }

    let mut candidate = input();
    candidate.render_budget_profile = Some(component(
        CinematicComponentRole::RenderBudgetProfile,
        "render-budget-2",
    ));
    let render_budget = config(candidate);
    assert_ne!(base.image_identity(), render_budget.image_identity());
    assert_eq!(base.audio_identity(), render_budget.audio_identity());

    let mut candidate = input();
    candidate.audio_budget_profile = Some(component(
        CinematicComponentRole::AudioBudgetProfile,
        "audio-budget-2",
    ));
    let audio_budget = config(candidate);
    assert_eq!(base.image_identity(), audio_budget.image_identity());
    assert_ne!(base.audio_identity(), audio_budget.audio_identity());
}
